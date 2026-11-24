//
// Created by gillard on 2026-07-03.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_nispSCA.cxx
//
// Crescendo test suite for nisp::SCA: (a) construction/initialize basics,
// (b) a single small-detector MACC exposure with statistical noise-budget
// checks, (c) persistence, (d) 1/f noise, (e) all three combined as a smoke
// test, (f) a full 2048x2048 exposure timed for a CPU-cost sanity check,
// (g) per-exposure parameters (nisp::readout, expose(flux, readout)),
// (h) full well and non-linearity in every readout mode,
// (i) photon transfer curve: gain and read noise recovered from CDS and DPU exposures.
//
// Most tests use the default MACC mode. The readout modes (selected with
// "exposureMode") are validated on a full detector in gTest_nispReadoutModes.cxx;
// here the CDS noise is checked by ExposeCdsTwoReadsOfReadNoiseAndKtcCancels and by (i).
//

#include <gtest/gtest.h>
#include <nisp/nispSCA.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <valarray>
#include <vector>

using namespace nisp;
using namespace phot2pix;

// ── Helpers ──────────────────────────────────────────────────────────────────

// Write a test output to ./build/<name> for inspection. CFITSIO refuses to create
// a file that already exists (status 105), so the file left by a previous run is
// purged first; the directory is created if missing.
static void writeTestOutput(const std::shared_ptr<DSL::FITScube>& img, const std::string& name)
{
    const std::filesystem::path dir  = "./build";
    const std::filesystem::path file = dir / name;

    std::filesystem::create_directories(dir);
    std::filesystem::remove(file);   // no-op when absent
    img->Write(file.string());
}

static std::valarray<double> makeUniformFlux(size_t w, size_t h, double value)
{
    return std::valarray<double>(value, w * h);
}

// The matrix starts at detector pixel (0, 0) by default, inside the 4-pixel reference
// border (QE = 0). Tests that need every pixel active place it on the active area.
static nlohmann::json activeOrigin()
{
    return nlohmann::json::array({4, 4});
}

// Write a W × H double map whose pixel (x, y) holds value(x, y) to test_data/<name> and
// return its path. CFITSIO does not overwrite: a file left by a previous run is removed first.
template<typename F>
static std::string writeTestMap(const std::string& name, size_t W, size_t H, F value)
{
    const std::filesystem::path dir  = "test_data";
    const std::filesystem::path file = dir / name;

    std::filesystem::create_directories(dir);
    std::filesystem::remove(file);

    DSL::FITSimg<double> map({W, H});
    for (size_t y = 0; y < H; ++y)
        for (size_t x = 0; x < W; ++x)
            map.SetPixelValue(value(x, y), y * W + x);
    map.Write(file.string());
    return file.string();
}

// Zero intrinsic pixel noise/response contribution (dark=0, readNoise=0,
// gain=1, QE=1, KTC=0, baseline=0) so only the effect under test shows up;
// the matrix is on the active area, so no pixel is a reference pixel.
static PixelConfig makeQuietPixelCfg()
{
    return PixelConfig{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 1.0},
        {"KTCnoise", 0.0}, {"baseline", 0.0}, {"detectorOrigin", activeOrigin()}
    };
}

static double populationVariance(const std::valarray<uint16_t>& d)
{
    double m = 0.0;
    for (auto v : d) m += static_cast<double>(v);
    m /= static_cast<double>(d.size());

    double var = 0.0;
    for (auto v : d) { const double x = static_cast<double>(v) - m; var += x * x; }
    return var / static_cast<double>(d.size());
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (a) Construction / initialize basics
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCA, DefaultConstruction2048)
{
    nisp::SCA sca;
    EXPECT_EQ(sca.getWidth(),  2048u);
    EXPECT_EQ(sca.getHeight(), 2048u);
    EXPECT_EQ(sca.size(),      2048u * 2048u);
}

TEST(SCA, SizedConstruction)
{
    nisp::SCA sca(16, 8);
    EXPECT_EQ(sca.getWidth(),  16u);
    EXPECT_EQ(sca.getHeight(), 8u);
    EXPECT_EQ(sca.size(),      128u);
}

TEST(SCA, InitializeAppliesBasePixelMaps)
{
    nisp::SCA sca(4, 4);
    PixelConfig cfg{
        {"gain", 2.5}, {"dark", 0.01}, {"readNoise", 3.0},
        {"QE", 0.8}, {"KTCnoise", 1.0}, {"baseline", 50.0},
        {"detectorOrigin", activeOrigin()}
    };
    sca.initialize(cfg);

    for (size_t i = 0; i < sca.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(sca.at(i).getGain(),              2.5);
        EXPECT_DOUBLE_EQ(sca.at(i).getDarkCurrent(),       0.01);
        EXPECT_DOUBLE_EQ(sca.at(i).getReadNoise(),         3.0);
        EXPECT_DOUBLE_EQ(sca.at(i).getQuantumEfficiency(), 0.8);
        EXPECT_DOUBLE_EQ(sca.at(i).getKTCnoise(),          1.0);
        EXPECT_DOUBLE_EQ(sca.at(i).getBaseline(),          50.0);
    }
}

TEST(SCA, InitializeAcceptsDpuAliasKeys)
{
    // dpuScale/dpuOffset are only applied by _toUint16() (windowed expose()
    // path); expose_macc's ramp path converts via FITScube::cast<uint16_t>(),
    // which does NOT apply them (they're written as descriptive HDU metadata
    // only in that path) — so verify offset/scale via the windowed call, and
    // frame_time via the MACC path's FRTIME header key.
    nisp::SCA sca(2, 2);
    PixelConfig cfg{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 0.0},
        {"KTCnoise", 0.0}, {"baseline", 0.0},
        {"dpu_offset", 100.0}, {"dpu_scale", 2.0}, {"frame_time", 0.5}
    };
    sca.initialize(cfg);

    auto flux = makeUniformFlux(2, 2, 0.0);
    auto windowed = sca.expose(flux, 1.0, 0, 0, 2, 2);

    writeTestOutput(windowed, "test_expose_macc_windows.fits");

    ASSERT_NE(windowed, nullptr);
    const auto& data = *(windowed->GetData<uint16_t>());
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 100u); // round(0*2.0 + 100.0 + 0.5) truncated -> 100

    auto ramp = sca.expose(flux, nisp::macc{2, 1, 0});
    writeTestOutput(ramp, "test_expose_macc_ramp.fits");
    ASSERT_NE(ramp, nullptr);
    EXPECT_NEAR(ramp->HDU().GetDoubleValueForKey("FRTIME"), 0.5, 1e-9);
}

TEST(SCA, AtAccessorsReadOnly)
{
    nisp::SCA sca(4, 3);
    sca.initialize(PixelConfig{{"gain", 3.0}});

    for (size_t row = 0; row < 3; ++row)
        for (size_t col = 0; col < 4; ++col)
            EXPECT_DOUBLE_EQ(sca.at(col, row).getGain(), sca.at(row * 4 + col).getGain());
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (b) Single MACC exposure, small detector
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCA, ExposeMaccReturnsExpectedLayerCount)
{
    nisp::SCA sca(4, 4);
    sca.initialize(PixelConfig{{"gain", 1.0}});
    auto flux = makeUniformFlux(4, 4, 10.0);

    auto result = sca.expose(flux, nisp::macc{5, 1, 0});

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(3), 5u);

    writeTestOutput(result, "test_expose_macc01.fits");
}

TEST(SCA, ExposeMaccLayerMeansGrowLinearlyStatistical)
{
    // NOTE: pixelresponse()'s shot-noise sigma is sqrt(signal),
    // where signal = flux*QE*t + dark*t. With flux=20, QE=1 that signal is
    // itself nonzero, so photon shot noise is present regardless of
    // readNoise/KTC/baseline being zeroed -- there is no way to get an
    // exactly deterministic nonzero ramp through the public API. Average
    // over many pixels instead, so shot noise averages down and the
    // cumulative-sum ramp-building logic can be checked statistically.
    constexpr double F = 20.0;
    constexpr double T = 1.0;
    constexpr size_t W = 64, H = 64; // 4096 pixels for spatial averaging
    nisp::SCA sca(W, H);
    sca.initialize(PixelConfig{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 1.0},
        {"KTCnoise", 0.0}, {"baseline", 0.0}, {"frame_time", T},
        {"detectorOrigin", activeOrigin()}
    });
    auto flux = makeUniformFlux(W, H, F);

    constexpr size_t nG = 4;
    auto result = sca.expose(flux, nisp::macc{nG, 1, 0});
    writeTestOutput(result, "test_expose_macc02.fits");

    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->Size(3), nG);

    const auto& data = *(result->GetData<uint16_t>());
    const size_t planeSize = result->Size(1) * result->Size(2);
    const double n = static_cast<double>(planeSize);

    constexpr double k = 5.0;
    for (size_t g = 0; g < nG; ++g)
    {
        // gain=1 -> Var[cumulative after (g+1) independent frame draws] = (g+1)*F*T
        const double expectedMean = static_cast<double>(g + 1) * F * T;
        const double expectedVar  = static_cast<double>(g + 1) * F * T;

        double sum = 0.0;
        for (size_t p = 0; p < planeSize; ++p)
            sum += static_cast<double>(data[g * planeSize + p]);
        const double sampleMean = sum / n;

        EXPECT_NEAR(sampleMean, expectedMean, k * std::sqrt(expectedVar / n) + 0.5)
            << "group " << g;
    }
}

TEST(SCA, ExposeMaccStatisticalNoiseBudget)
{
    constexpr double dark = 1.0, readNoise = 2.0, gain = 1.5, QE = 0.9, flux = 20.0, T = 1.0;
    nisp::SCA sca(1, 1);
    sca.initialize(PixelConfig{
        {"gain", gain}, {"dark", dark}, {"readNoise", readNoise}, {"QE", QE},
        {"KTCnoise", 0.0}, {"baseline", 0.0}, {"frame_time", T},
        {"detectorOrigin", activeOrigin()}
    });
    auto fluxArr = makeUniformFlux(1, 1, flux);

    constexpr double signal_e = (flux * QE + dark) * T;
    const double mu  = signal_e / gain;
    const double var = (signal_e + readNoise * readNoise) / (gain * gain);

    constexpr size_t nSample = 3000;
    const double n = static_cast<double>(nSample);

    double sigMean = 0.0, sigVar = 0.0;
    for (size_t i = 0; i < nSample; ++i)
    {
        auto result = sca.expose(fluxArr, nisp::macc{1, 1, 0});
        const double x = static_cast<double>((*(result->GetData<uint16_t>()))[0]);
        sigMean += x;
        sigVar  += x * x;
    }
    sigMean /= n;
    sigVar   = sigVar / n - sigMean * sigMean;

    // k=6 (vs k=5 elsewhere) plus a flat slack term: absorbs the +/-0.5 uint16
    // rounding quantization on top of the underlying Gaussian noise budget.
    constexpr double k = 6.0;
    EXPECT_NEAR(sigMean, mu,  k * std::sqrt(var / n) + 0.5);
    EXPECT_NEAR(sigVar,  var, k * std::sqrt(2.0) * var / std::sqrt(n) + 0.05 * var);
}

TEST(SCA, ExposeMaccBaselineDrawnOnceNotPerLayer)
{
    nisp::SCA sca(1, 1);
    sca.initialize(PixelConfig{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 1.0},
        {"KTCnoise", 5.0}, {"baseline", 100.0}, {"frame_time", 1.0}
    });
    auto flux = makeUniformFlux(1, 1, 0.0);

    auto result = sca.expose(flux, nisp::macc{4, 1, 0});
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->Size(3), 4u);

    const auto& data = *(result->GetData<uint16_t>());
    // 1x1 detector -> one pixel per layer; baseline/KTC drawn once on layer 0,
    // then forced off (cfg "baseline"="off") for every subsequent layer.
    EXPECT_EQ(data[1], data[0]);
    EXPECT_EQ(data[2], data[0]);
    EXPECT_EQ(data[3], data[0]);
}

TEST(SCA, ExposeMaccReadNoiseNotAccumulated)
{
    // Zero flux: every read is baseline + its own read noise. Read noise is a
    // read property, not charge, so every layer of the ramp must show RN²
    // (not k*RN² for the k-th read), and two consecutive reads differ by 2*RN².
    // Baseline 1000 keeps every value far from the uint16 clamp at 0.
    constexpr double RN = 13.0, base = 1000.0;
    constexpr size_t W = 256, H = 256, nG = 4, nR = 16, nD = 4;
    nisp::SCA sca(W, H);
    sca.initialize(PixelConfig{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", RN}, {"QE", 1.0},
        {"KTCnoise", 0.0}, {"baseline", base}, {"frame_time", 1.0}
    });
    auto flux = makeUniformFlux(W, H, 0.0);

    auto result = sca.expose(flux, nisp::macc{nG, nR, nD});
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->Size(3), nG * nR);

    const auto&      data  = *(result->GetData<uint16_t>());
    constexpr size_t plane = W * H;   // constant expression: usable in the lambda without capture
    const double     n     = static_cast<double>(plane);

    // Mean and population variance of plane a, or of (plane a − plane b) when b >= 0.
    auto stats = [&data, n](size_t a, long b) -> std::pair<double, double>
    {
        double sum = 0.0, sum2 = 0.0;
        for (size_t p = 0; p < plane; ++p)
        {
            double x = static_cast<double>(data[a * plane + p]);
            if (b >= 0)
                x -= static_cast<double>(data[static_cast<size_t>(b) * plane + p]);
            sum  += x;
            sum2 += x * x;
        }
        const double mean = sum / n;
        return {mean, sum2 / n - mean * mean};
    };

    constexpr double k       = 5.0;
    const double     varRead = RN * RN + 1.0 / 12.0;   // read noise + uint16 rounding
    for (size_t L = 0; L < nG * nR; ++L)
    {
        const auto [mean, var] = stats(L, -1);
        EXPECT_NEAR(mean, base,    k * std::sqrt(varRead / n) + 0.01)                  << "layer " << L;
        EXPECT_NEAR(var,  varRead, k * std::sqrt(2.0 / n) * varRead + 0.01 * varRead) << "layer " << L << " (read noise must not accumulate)";

        if (L == 0)
            continue;
        const auto [dMean, dVar] = stats(L, static_cast<long>(L - 1));
        EXPECT_NEAR(dMean, 0.0,         k * std::sqrt(2.0 * varRead / n) + 0.01)                    << "read " << L << " - read " << L - 1;
        EXPECT_NEAR(dVar,  2.0 * varRead, k * std::sqrt(2.0 / n) * 2.0 * varRead + 0.02 * varRead) << "read " << L << " - read " << L - 1;
    }
}

TEST(SCA, ExposeCdsTwoReadsOfReadNoiseAndKtcCancels)
{
    // Zero flux, CDS over 10 frames: CDS = read N − read 1. The reset level
    // (baseline + kTC) is common to both reads and cancels; each read keeps
    // its own read noise, so the CDS variance is 2*RN² (not RN², not + kTC²).
    // The DPU offset keeps the zero-mean CDS away from the uint16 clamp at 0.
    constexpr double RN = 13.0, KTC = 20.0, base = 1000.0, offset = 500.0;
    constexpr size_t W = 256, H = 256;
    nisp::SCA sca(W, H);
    sca.initialize(PixelConfig{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", RN}, {"QE", 1.0},
        {"KTCnoise", KTC}, {"baseline", base}, {"frame_time", 1.0},
        {"offset", offset}, {"exposureMode", std::string("CDS")}
    });
    auto flux = makeUniformFlux(W, H, 0.0);

    auto result = sca.expose(flux, 10.0);
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->GetDimension(), 2u);

    const auto&  data = *(result->GetData<uint16_t>());
    const double n    = static_cast<double>(data.size());
    double sum = 0.0;
    for (auto v : data) sum += static_cast<double>(v);
    const double mean = sum / n;
    const double var  = populationVariance(data);

    constexpr double k      = 5.0;
    const double     varCds = 2.0 * RN * RN + 1.0 / 12.0;
    EXPECT_NEAR(mean, offset, k * std::sqrt(varCds / n) + 0.01) << "baseline must cancel";
    EXPECT_NEAR(var,  varCds, k * std::sqrt(2.0 / n) * varCds + 0.01 * varCds) << "CDS variance: 2*RN², kTC must cancel";

    EXPECT_EQ  (result->HDU().GetValueForKey("MODE"), "CDS");
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("INTTIME"), 9.0,  1e-9);   // read 1 → read 10
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("EXPTIME"), 10.0, 1e-9);   // reset → read 10
}

TEST(SCA, ExposeMaccHDUKeys)
{
    nisp::SCA sca(2, 2);
    sca.initialize(PixelConfig{{"gain", 1.0}});
    auto flux = makeUniformFlux(2, 2, 5.0);

    auto result = sca.expose(flux, nisp::macc{3, 2, 1});
    ASSERT_NE(result, nullptr);

    EXPECT_EQ  (result->HDU().GetValueForKey("MODE"), "MACC");
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("NG"),   3.0,  1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("NR"),   2.0,  1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("ND"),   1.0,  1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("GAIN"), 1.0,  1e-9); // mean pixel gain ("gain": 1.0)
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("RN"),   0.0,  1e-9); // mean pixel read noise (pixsim default 0)
    EXPECT_FALSE(result->HDU().Exists("PERSIST"));
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (c) Persistence
//
//  NOTE: persistence::scheduler stores its PixelConfig by reference, so every
//  cfg passed to initialize() when "persistence" is present must be a named
//  local that outlives the SCA object — never a temporary.
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCA, PersistenceIncreasesSecondExposureSignal)
{
    // The persistence contribution itself is deterministic (it's computed
    // from the flux*QE*t charge formula, not from the noisy pixelresponse
    // draw), but the base pixel signal here is nonzero (flux=100, QE=1), so
    // each single-pixel readout still carries photon shot noise on top of
    // it. A single-pixel comparison can therefore be dominated by shot-noise
    // scatter rather than the (much smaller) persistence shift -- average
    // over many pixels so the deterministic persistence shift dominates.
    const size_t W = 64, H = 64;

    PixelConfig cfgOn = makeQuietPixelCfg();
    nlohmann::json traps = nlohmann::json::array();
    traps.push_back({{"f", 0.02}, {"t", 50.0}});
    cfgOn.set("persistence", traps);

    nisp::SCA scaOn(W, H);
    scaOn.initialize(cfgOn); // cfgOn outlives scaOn (declared first, destroyed last)

    PixelConfig cfgOff = makeQuietPixelCfg(); // no "persistence" key: control
    nisp::SCA scaOff(W, H);
    scaOff.initialize(cfgOff);

    auto flux = makeUniformFlux(W, H, 100.0);

    auto onFirst  = scaOn.expose(flux, nisp::macc{1, 1, 0});
    auto offFirst = scaOff.expose(flux, nisp::macc{1, 1, 0});
    ASSERT_NE(onFirst,  nullptr);
    ASSERT_NE(offFirst, nullptr);

    auto onSecond  = scaOn.expose(flux, nisp::macc{1, 1, 0});
    auto offSecond = scaOff.expose(flux, nisp::macc{1, 1, 0});
    ASSERT_NE(onSecond,  nullptr);
    ASSERT_NE(offSecond, nullptr);

    auto meanOf = [](const std::shared_ptr<DSL::FITScube>& cube) -> double
    {
        const auto& data = *(cube->GetData<uint16_t>());
        double sum = 0.0;
        for (auto v : data) sum += static_cast<double>(v);
        return sum / static_cast<double>(data.size());
    };

    EXPECT_GT(meanOf(onSecond), meanOf(offSecond));
}

TEST(SCA, PersistenceHDUKeyPresent)
{
    PixelConfig cfgOn = makeQuietPixelCfg();
    nlohmann::json traps = nlohmann::json::array();
    traps.push_back({{"f", 0.02}, {"t", 50.0}});
    cfgOn.set("persistence", traps);

    nisp::SCA scaOn(2, 2);
    scaOn.initialize(cfgOn);

    PixelConfig cfgOff = makeQuietPixelCfg();
    nisp::SCA scaOff(2, 2);
    scaOff.initialize(cfgOff);

    auto flux = makeUniformFlux(2, 2, 10.0);
    auto onResult  = scaOn.expose(flux, nisp::macc{1, 1, 0});
    auto offResult = scaOff.expose(flux, nisp::macc{1, 1, 0});

    EXPECT_EQ(onResult->HDU().GetValueForKey("PERSIST"), "ON");
    EXPECT_FALSE(offResult->HDU().Exists("PERSIST"));
}

TEST(SCA, PersistenceCannotBeDisabledOnceEnabled)
{
    // See PersistenceIncreasesSecondExposureSignal: average over many pixels
    // so the deterministic persistence shift isn't swamped by per-pixel
    // photon shot noise (base signal is nonzero here: flux=50, QE=1).
    const size_t W = 64, H = 64;

    PixelConfig cfgOn = makeQuietPixelCfg();
    nlohmann::json traps = nlohmann::json::array();
    traps.push_back({{"f", 0.02}, {"t", 50.0}});
    cfgOn.set("persistence", traps);

    nisp::SCA sca(W, H);
    sca.initialize(cfgOn); // enables persistence; scheduler stores cfgOn by reference

    PixelConfig cfgNoPersist = makeQuietPixelCfg(); // no "persistence" key
    sca.initialize(cfgNoPersist);                   // re-initialize WITHOUT the key

    PixelConfig cfgControl = makeQuietPixelCfg();
    nisp::SCA scaControl(W, H);
    scaControl.initialize(cfgControl);

    auto flux = makeUniformFlux(W, H, 50.0);

    sca.expose(flux, nisp::macc{1, 1, 0});        // first exposure: registers charge
    scaControl.expose(flux, nisp::macc{1, 1, 0}); // control's first exposure

    auto second        = sca.expose(flux, nisp::macc{1, 1, 0});
    auto controlSecond = scaControl.expose(flux, nisp::macc{1, 1, 0});

    auto meanOf = [](const std::shared_ptr<DSL::FITScube>& cube) -> double
    {
        const auto& data = *(cube->GetData<uint16_t>());
        double sum = 0.0;
        for (auto v : data) sum += static_cast<double>(v);
        return sum / static_cast<double>(data.size());
    };

    EXPECT_GT(meanOf(second), meanOf(controlSecond));
}

TEST(SCA, PersistenceReleaseAndDelayBetweenExposures)
{
    // Exposure 1 stores E = 1e5 e- per pixel (flux 1e5 ph/s, QE 1, one 1 s frame).
    // Exposure 2 is dark and noiseless (gain 1), so its single read is the charge
    // released during its frame, which starts dT after exposure 1:
    //   P = E·f·τ·e^(−dT/τ)·(1 − e^(−t_f/τ))
    // ≈ 99.5 ADU with dT = 0 and ≈ 36.6 ADU with dT = 100 s. (Before the fixes: the
    // interval [−t_f, 0] gave ≈ 100.5 ADU, and dT was ignored.)
    const double E = 1e5, f = 1e-3, tau = 100.0, tf = 1.0;
    const size_t W = 4, H = 4;

    auto secondRead = [&](double delay)
    {
        PixelConfig cfg = makeQuietPixelCfg();   // declared first: outlives the SCA
        cfg.set("frame_time", tf);
        nlohmann::json p;
        p["traps"] = nlohmann::json::array({ {{"f", f}, {"t", tau}} });
        p["dT"]    = delay;
        cfg.set("persistence", p);

        nisp::SCA sca(W, H);
        sca.initialize(cfg);
        sca.expose(makeUniformFlux(W, H, E / tf), nisp::macc{1, 1, 0});
        auto second = sca.expose(makeUniformFlux(W, H, 0.0), nisp::macc{1, 1, 0});

        const auto& data = *(second->GetData<uint16_t>());
        double mean = 0.0;
        for (auto v : data) mean += static_cast<double>(v);
        return mean / static_cast<double>(data.size());
    };

    for (double delay : {0.0, 100.0})
    {
        const double expected = E * f * tau * std::exp(-delay / tau) * (1.0 - std::exp(-tf / tau));
        EXPECT_NEAR(secondRead(delay), expected, 0.51) << "dT = " << delay << " s";   // uint16 rounding
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (d) 1/f noise
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCA, OneOverFNoiseAddsVarianceVsDisabled)
{
    const size_t W = 128, H = 32;

    // baseline 1000 ADU keeps the reads away from 0: with baseline 0, a channel
    // whose 1/f level happens to be negative would be clamped to 0 by the uint16
    // conversion and hide its noise.
    PixelConfig cfgOn{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 0.0},
        {"KTCnoise", 0.0}, {"baseline", 1000.0}, {"oneOverF", 15.0}
    };
    nisp::SCA scaOn(W, H);
    scaOn.initialize(cfgOn);

    PixelConfig cfgOff{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 0.0},
        {"KTCnoise", 0.0}, {"baseline", 1000.0}
    };
    nisp::SCA scaOff(W, H);
    scaOff.initialize(cfgOff);

    auto flux = makeUniformFlux(W, H, 0.0);

    // QE=0, dark=0, readNoise=0, KTC=0 -> intrinsic pixelresponse is exactly the
    // baseline everywhere; all frame structure comes from the noise map.
    auto onResult  = scaOn.expose(flux, nisp::macc{1, 1, 0});
    auto offResult = scaOff.expose(flux, nisp::macc{1, 1, 0});

    EXPECT_GT(populationVariance(*(onResult->GetData<uint16_t>())), 1.0);
    EXPECT_NEAR(populationVariance(*(offResult->GetData<uint16_t>())), 0.0, 1e-9);
}

TEST(SCA, OneOverFNoisePerChannelColumnMeanVarianceElevated)
{
    const size_t W = 128, H = 32; // at (0, 0): the 2 whole physical channels 0 and 1

    PixelConfig cfg{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 0.0},
        {"KTCnoise", 0.0}, {"baseline", 1000.0}, {"oneOverF", 15.0}   // baseline: no clamping at 0
    };
    nisp::SCA sca(W, H);
    sca.initialize(cfg);

    auto flux = makeUniformFlux(W, H, 0.0);
    auto result = sca.expose(flux, nisp::macc{1, 1, 0});
    const auto& data = *(result->GetData<uint16_t>());

    std::vector<double> colMeans(W, 0.0);
    for (size_t c = 0; c < W; ++c)
    {
        double sum = 0.0;
        for (size_t r = 0; r < H; ++r)
            sum += static_cast<double>(data[r * W + c]);
        colMeans[c] = sum / static_cast<double>(H);
    }

    double m = 0.0;
    for (double v : colMeans) m += v;
    m /= static_cast<double>(colMeans.size());

    double var = 0.0;
    for (double v : colMeans) var += (v - m) * (v - m);
    var /= static_cast<double>(colMeans.size());

    // Banding: column means retain structure after averaging over rows
    // (unlike white noise, which would average toward 0 variance).
    EXPECT_GT(var, 1.0);
}

TEST(SCA, OneOverFNoiseHDUAndPersistenceIndependentFlags)
{
    // Object form with a non-default colPerChannel (4 divides the 2048-wide
    // detector: 512 channels); the 4 × 4 matrix is on its active area.
    PixelConfig cfg = makeQuietPixelCfg();
    nlohmann::json oneOverF;
    oneOverF["amplitude"]     = 5.0;
    oneOverF["colPerChannel"] = 4;
    cfg.set("oneOverF", oneOverF); // no "persistence" key

    nisp::SCA sca(4, 4);
    sca.initialize(cfg);

    auto flux = makeUniformFlux(4, 4, 10.0);
    auto result = sca.expose(flux, nisp::macc{1, 1, 0});

    EXPECT_FALSE(result->HDU().Exists("PERSIST"));
}

// ── Readout noise: units, bias voltage noise, pedestal drift ─────────────────

TEST(SCA, ReadoutNoiseIsConvertedToAduWithPixelGain)
{
    // Constant pedestal of 400 e- (drift = a0 only), gain 4 e-/ADU, every other
    // effect off: each read is 400 / 4 = 100 ADU.
    const size_t W = 128, H = 8;
    PixelConfig cfg{
        {"gain", 4.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 0.0},
        {"KTCnoise", 0.0}, {"baseline", 0.0}
    };
    nlohmann::json ped;
    ped["drift"]["poly"] = nlohmann::json::array({400.0});
    cfg.set("pedestal_drift", ped);

    nisp::SCA sca(W, H);
    sca.initialize(cfg);

    auto result = sca.expose(makeUniformFlux(W, H, 0.0), nisp::macc{2, 1, 0});
    for (auto v : *(result->GetData<uint16_t>()))
        ASSERT_EQ(v, 100u);
}

TEST(SCA, BiasVoltageNoiseIsIdenticalInEveryChannel)
{
    // Bias voltage noise only (baseline 1000 ADU keeps the reads away from 0): pixel s
    // of a row reads the same value in every channel, unlike the per-channel 1/f noise.
    const size_t W = 128, H = 32;   // 2 channels of 64 columns
    PixelConfig cfg{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 0.0},
        {"KTCnoise", 0.0}, {"baseline", 1000.0}, {"bias_voltage_noise", 15.0}
    };
    nisp::SCA sca(W, H);
    sca.initialize(cfg);

    auto result = sca.expose(makeUniformFlux(W, H, 0.0), nisp::macc{2, 1, 0});
    const auto& data = *(result->GetData<uint16_t>());

    EXPECT_GT(populationVariance(data), 1.0);   // the noise is there
    for (size_t layer = 0; layer < 2; ++layer)
        for (size_t y = 0; y < H; ++y)
            for (size_t s = 0; s < 64; ++s)
                ASSERT_EQ(data[layer * W * H + y * W + s], data[layer * W * H + y * W + 64 + s])
                    << "layer " << layer << ", row " << y << ", column " << s;
}

TEST(SCA, PedestalDriftSurvivesInCds)
{
    // Linear drift f(t) = 3 t e-, gain 1, frame time 1 s, CDS over 10 s: read 1 is
    // frame 0 and read N is frame 9. Each row is read at the same point of both frames,
    // so CDS = f(9 + u) − f(u) = 27 ADU for every pixel.
    const size_t W = 128, H = 8;
    PixelConfig cfg{
        {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 0.0},
        {"KTCnoise", 0.0}, {"baseline", 0.0}, {"frame_time", 1.0},
        {"exposureMode", std::string("CDS")}
    };
    nlohmann::json ped;
    ped["drift"]["poly"] = nlohmann::json::array({0.0, 3.0});
    cfg.set("pedestal_drift", ped);

    nisp::SCA sca(W, H);
    sca.initialize(cfg);

    auto result = sca.expose(makeUniformFlux(W, H, 0.0), 10.0);
    for (auto v : *(result->GetData<uint16_t>()))
        ASSERT_EQ(v, 27u);
}

// ── Detector geometry and reference pixels ───────────────────────────────────

TEST(SCA, ReferencePixelsFollowThePhysicalBorder)
{
    // A 2044 × 16 matrix at physical column 2 and row 1016: it keeps 2 reference
    // columns on each side (matrix columns 0–1 and 2042–2043) and no reference row.
    // Reference pixels get QE = 0 and, here, refDark.
    nisp::SCA sca(2044, 16);
    sca.initialize(PixelConfig{{"QE", 0.9}, {"dark", 0.01}, {"refDark", 0.002},
                               {"detectorOrigin", nlohmann::json::array({2, 1016})}});

    for (size_t row = 0; row < 16; ++row)
        for (size_t col = 0; col < 2044; ++col)
        {
            const bool ref = col < 2 || col >= 2042;
            ASSERT_DOUBLE_EQ(sca.at(col, row).getQuantumEfficiency(), ref ? 0.0   : 0.9)  << "column " << col << ", row " << row;
            ASSERT_DOUBLE_EQ(sca.at(col, row).getDarkCurrent(),       ref ? 0.002 : 0.01) << "column " << col << ", row " << row;
        }
}

TEST(SCA, DefaultOriginIsTheDetectorCorner)
{
    // Without "detectorOrigin", a 16 × 16 matrix starts at the first detector pixel:
    // 4 reference columns on the left and 4 reference rows at the bottom.
    nisp::SCA sca(16, 16);
    sca.initialize(PixelConfig{{"QE", 0.9}});
    EXPECT_EQ(sca.getOrigin(), std::make_pair(size_t(0), size_t(0)));

    for (size_t row = 0; row < 16; ++row)
        for (size_t col = 0; col < 16; ++col)
            ASSERT_DOUBLE_EQ(sca.at(col, row).getQuantumEfficiency(), (col < 4 || row < 4) ? 0.0 : 0.9)
                << "column " << col << ", row " << row;
}

TEST(SCA, RefBorderZeroDisablesReferencePixels)
{
    nisp::SCA sca(2048, 4);
    sca.initialize(PixelConfig{{"QE", 0.9}, {"refBorder", 0}});

    for (size_t i = 0; i < sca.size(); ++i)
        ASSERT_DOUBLE_EQ(sca.at(i).getQuantumEfficiency(), 0.9);
}

TEST(SCA, MatrixLargerThanTheDetectorThrows)
{
    nisp::SCA wide(2050, 4);
    EXPECT_THROW(wide.initialize(PixelConfig{{"gain", 1.0}}), std::invalid_argument);

    // ... unless the detector is declared larger.
    nisp::SCA declared(2050, 4);
    EXPECT_NO_THROW(declared.initialize(PixelConfig{{"gain", 1.0}, {"detectorSize", nlohmann::json::array({4096, 4096})}}));
}

TEST(SCA, DetectorSizedMapIsCroppedAndKeepsTheReferencePixels)
{
    // QE map of a 64 × 64 detector (refBorder 4). A 16 × 16 matrix at the default origin
    // takes its corner, whose 4 left columns and 4 bottom rows are reference pixels
    // (QE = 0); at [24, 24] it takes the section [24, 40) × [24, 40), with no reference pixel.
    auto qe = [](size_t x, size_t y) { return static_cast<double>(x + 64 * y + 1) / 4096.0; };
    const std::string path = writeTestMap("sca_qe_64x64.fits", 64, 64, qe);

    nisp::SCA sca(16, 16);
    sca.initialize(PixelConfig{{"QE", path}, {"detectorSize", nlohmann::json::array({64, 64})}});
    for (size_t row = 0; row < 16; ++row)
        for (size_t col = 0; col < 16; ++col)
            ASSERT_DOUBLE_EQ(sca.at(col, row).getQuantumEfficiency(), (col < 4 || row < 4) ? 0.0 : qe(col, row))
                << "origin (0, 0), column " << col << ", row " << row;

    sca.initialize(PixelConfig{{"QE", path}, {"detectorOrigin", nlohmann::json::array({24, 24})}});
    for (size_t row = 0; row < 16; ++row)
        for (size_t col = 0; col < 16; ++col)
            ASSERT_DOUBLE_EQ(sca.at(col, row).getQuantumEfficiency(), qe(24 + col, 24 + row))
                << "origin (24, 24), column " << col << ", row " << row;

    // A map that is neither 16 × 16 nor 64 × 64 is rejected.
    const std::string wrong = writeTestMap("sca_qe_32x32.fits", 32, 32, qe);
    EXPECT_THROW(sca.initialize(PixelConfig{{"QE", wrong}}), std::invalid_argument);

    std::filesystem::remove(path);
    std::filesystem::remove(wrong);
}

TEST(SCA, OutputHeaderRecordsTheInputs)
{
    // A 128 × 8 matrix at the default origin covers the physical columns 1–128 and rows 1–8.
    PixelConfig cfg{
        {"gain", 1.92}, {"dark", 0.01}, {"readNoise", 13.0}, {"QE", 0.94},
        {"KTCnoise", 48.0}, {"baseline", 1000.0}, {"refDark", 0.002}, {"oneOverF", 7.0}
    };
    nlohmann::json ped;
    ped["sigma"]         = 4.0;
    ped["drift"]["poly"] = nlohmann::json::array({0.0, 0.05});
    cfg.set("pedestal_drift", ped);

    nisp::SCA sca(128, 8);
    sca.initialize(cfg);
    auto result = sca.expose(makeUniformFlux(128, 8, 0.0), nisp::macc{2, 1, 0});
    const auto& hdu = result->HDU();

    EXPECT_NEAR(hdu.GetDoubleValueForKey("PIXGAIN"), 1.92,   1e-9);
    EXPECT_NEAR(hdu.GetDoubleValueForKey("PIXDARK"), 0.01,   1e-9);
    EXPECT_NEAR(hdu.GetDoubleValueForKey("PIXRN"),   13.0,   1e-9);
    EXPECT_NEAR(hdu.GetDoubleValueForKey("PIXQE"),   0.94,   1e-9);
    EXPECT_NEAR(hdu.GetDoubleValueForKey("PIXKTC"),  48.0,   1e-9);
    EXPECT_NEAR(hdu.GetDoubleValueForKey("PIXBASE"), 1000.0, 1e-9);
    EXPECT_NEAR(hdu.GetDoubleValueForKey("REFDARK"), 0.002,  1e-9);
    EXPECT_NEAR(hdu.GetDoubleValueForKey("REFBORD"), 4.0,    1e-9);
    EXPECT_EQ  (hdu.GetValueForKey("DETSIZE"), "[1:2048,1:2048]");
    EXPECT_EQ  (hdu.GetValueForKey("DETSEC"),  "[1:128,1:8]");

    EXPECT_NEAR(hdu.GetDoubleValueForKey("OOFAMP"),   7.0, 1e-9);
    EXPECT_NEAR(hdu.GetDoubleValueForKey("PEDSIGMA"), 4.0, 1e-9);
    EXPECT_EQ  (hdu.GetValueForKey("PEDPOLY"), "0, 0.05");
    EXPECT_FALSE(hdu.Exists("BIASAMP"));   // bias voltage noise off

    EXPECT_TRUE(hdu.Exists("PHOT2PIX"));
}

TEST(SCA, ChannelLayoutIsCheckedAtInitialize)
{
    // colPerChannel 100 does not divide the 2048-wide detector.
    PixelConfig cfg = makeQuietPixelCfg();
    nlohmann::json oneOverF;
    oneOverF["amplitude"]     = 5.0;
    oneOverF["colPerChannel"] = 100;
    cfg.set("oneOverF", oneOverF);

    nisp::SCA sca(128, 8);
    EXPECT_THROW(sca.initialize(cfg), std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (e) Combined smoke test — persistence + 1/f + base noise together
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCA, CombinedPersistenceOneOverFAndBaseNoiseSmokeTest)
{
    PixelConfig cfg{
        {"gain", 2.0}, {"dark", 0.5}, {"readNoise", 4.0}, {"QE", 0.85},
        {"KTCnoise", 2.0}, {"baseline", 200.0}, {"oneOverF", 6.0}
    };
    nlohmann::json traps = nlohmann::json::array();
    traps.push_back({{"f", 0.01}, {"t", 80.0}});
    cfg.set("persistence", traps);

    nisp::SCA sca(32, 32);
    sca.initialize(cfg);

    auto flux = makeUniformFlux(32, 32, 30.0);

    std::shared_ptr<DSL::FITScube> first, second;
    EXPECT_NO_THROW(first  = sca.expose(flux, nisp::macc{2, 1, 0}));
    EXPECT_NO_THROW(second = sca.expose(flux, nisp::macc{2, 1, 0}));


    ASSERT_NE(first,  nullptr);
    ASSERT_NE(second, nullptr);

    writeTestOutput(first,  "test_expose_macc_combined01.fits");
    writeTestOutput(second, "test_expose_macc_combined02.fits");

    EXPECT_EQ(first->Size(3),  2u);
    EXPECT_EQ(second->Size(3), 2u);

    auto meanOfLastLayer = [](const std::shared_ptr<DSL::FITScube>& cube) -> double
    {
        const auto& data = *(cube->GetData<uint16_t>());
        const size_t planeSize = cube->Size(1) * cube->Size(2);
        const size_t lastLayerStart = (cube->Size(3) - 1) * planeSize;
        double sum = 0.0;
        for (size_t i = 0; i < planeSize; ++i)
            sum += static_cast<double>(data[lastLayerStart + i]);
        return sum / static_cast<double>(planeSize);
    };

    const double firstMean  = meanOfLastLayer(first);
    const double secondMean = meanOfLastLayer(second);

    // This is an integration smoke test (each effect is verified precisely
    // in isolation above) — loose direction-only check with a wide slack to
    // absorb per-call shot/read/baseline noise (readNoise=4, KTC=2 here can
    // plausibly shift a single realization by several ADU on their own).
    EXPECT_GE(secondMean, firstMean - 30.0);

    EXPECT_EQ(first->HDU().GetValueForKey("PERSIST"),  "ON");
    EXPECT_EQ(second->HDU().GetValueForKey("PERSIST"), "ON");
    EXPECT_EQ(first->HDU().GetValueForKey("MODE"),     "MACC");
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (f) 2048x2048 performance sanity check
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCA, FullDetectorMaccExposureCompletesWithinBudget)
{
    nisp::SCA sca; // default 2048x2048
    sca.initialize(PixelConfig{
        {"gain", 0.2}, {"readNoise", 0.0}, {"dark", 0.001}, {"QE", 0.9}, {"frame_time", 1.45408},
        {"oneOverF", nlohmann::json{{"nPoles", 8}, {"amplitude", 25.0},{"slope", 2.0}}} // well above typical NISP 5-10 e- RMS, so banding is clearly visible
    });

    std::valarray<double> flux(50.0, sca.size());

    const auto t0 = std::chrono::steady_clock::now();
    auto result = sca.expose(flux, nisp::macc{5, 1, 0});
    const auto t1 = std::chrono::steady_clock::now();

    const double seconds = std::chrono::duration<double>(t1 - t0).count();
    RecordProperty("wall_clock_seconds", std::to_string(seconds));
    std::cout << "[SCA 2048x2048 x5 groups] wall clock: " << seconds << " s\n";

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(3), 5u);
    // Deliberately generous: ~4.2M pixels x 5 groups spread across
    // hardware_concurrency() threads should finish in well under a few
    // seconds on any reasonable machine; 30s only guards against a
    // catastrophic regression (e.g. accidental serialization), not against
    // normal machine-to-machine performance variation.
    EXPECT_LT(seconds, 30.0);

    writeTestOutput(result, "test_expose_macc_full_detector.fits");
}

TEST(SCA, FullDetectorMaccExposureCorrectSize)
{
    nisp::SCA sca;
    sca.initialize(PixelConfig{{"gain", 1.0}});
    std::valarray<double> flux(10.0, sca.size());

    auto result = sca.expose(flux, nisp::macc{1, 1, 0});
    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), 2048u);
    EXPECT_EQ(result->Size(2), 2048u);
    EXPECT_EQ(result->GetEqBitPerPixel(), USHORT_IMG);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (g) Exposure parameters: readout::fromConfig() and expose(flux, readout)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCA, ReadoutFromConfigReadsOneTimingAndTheMode)
{
    const nisp::readout t = nisp::readout::fromConfig(PixelConfig{{"exposureTime", 10.0}});
    EXPECT_FALSE(t.maccParams);
    EXPECT_FALSE(t.exposureMode);
    EXPECT_DOUBLE_EQ(t.exposureTime, 10.0);

    const nisp::readout m = nisp::readout::fromConfig(PixelConfig{
        {"macc", nlohmann::json::array({3, 2, 1})}, {"exposureMode", std::string("cds")}});
    ASSERT_TRUE(m.maccParams);
    EXPECT_EQ(m.maccParams->nGroups, 3u);
    EXPECT_EQ(m.maccParams->nFrames, 2u);
    EXPECT_EQ(m.maccParams->nDrop,   1u);
    ASSERT_TRUE(m.exposureMode);
    EXPECT_EQ(*m.exposureMode, nisp::mode::CDS);   // case-insensitive

    const nisp::readout k = nisp::readout::fromConfig(PixelConfig{{"ng", 4}, {"nr", 16}, {"nd", 4}});
    ASSERT_TRUE(k.maccParams);
    EXPECT_EQ(k.maccParams->nGroups, 4u);

    const std::vector<PixelConfig> bad{
        PixelConfig{},                                                                    // no timing
        PixelConfig{{"exposureTime", 10.0}, {"macc", nlohmann::json::array({3, 2, 1})}},  // two timings
        PixelConfig{{"exposureTime", 0.0}},
        PixelConfig{{"exposureTime", -1.0}},
        PixelConfig{{"exposureTime", std::string("10")}},
        PixelConfig{{"exposureTime", 10.0}, {"exposureMode", std::string("FOO")}},
        PixelConfig{{"exposureTime", 10.0}, {"exposureMode", 3.0}},
    };
    for (const PixelConfig& cfg : bad)
        EXPECT_THROW(nisp::readout::fromConfig(cfg), std::invalid_argument) << cfg.toJson().dump();
}

TEST(SCA, ExposeWithReadoutUsesItsModeForThatExposureOnly)
{
    nisp::SCA sca(16, 16);
    sca.initialize(makeQuietPixelCfg());   // MACC mode, frame time 1.45408 s
    const auto flux = makeUniformFlux(16, 16, 10.0);

    // Same shapes as the expose() overloads it calls.
    const auto ramp = sca.expose(flux, nisp::readout{std::nullopt, nisp::macc{3, 2, 1}, 0.0});
    ASSERT_NE(ramp, nullptr);
    EXPECT_EQ(ramp->HDU().GetValueForKey("MODE"), "MACC");
    EXPECT_EQ(ramp->Size(3), 6u);   // every read is a layer: nGroups × nFrames

    const auto single = sca.expose(flux, nisp::readout{nisp::mode::FLUX, std::nullopt, 10.0});
    ASSERT_NE(single, nullptr);
    EXPECT_EQ(single->HDU().GetValueForKey("MODE"), "FLUX");

    // The SCA is back in MACC mode, also after an exposure that throws (CDS needs >= 2 frames).
    auto expectMacc = [&sca, &flux]()
    {
        const auto next = sca.expose(flux, nisp::macc{2, 1, 0});
        ASSERT_NE(next, nullptr);
        EXPECT_EQ(next->HDU().GetValueForKey("MODE"), "MACC");
        EXPECT_EQ(next->Size(3), 2u);
    };
    expectMacc();

    EXPECT_THROW(sca.expose(flux, nisp::readout{nisp::mode::CDS, std::nullopt, 0.5}), std::invalid_argument);
    expectMacc();
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (h) Full well and non-linearity
//
//  The full well is 20 ke- so that it is reached far below the ADC limit (65535 ADU):
//  the NISP one (130 ke-) is above it at the NISP gain and baseline, and never shows in a
//  MACC ramp. Response: P(Q) = Q·(1 + c2·Q + c3·Q²) of Q' = min(Q, 20 ke-), that is −4.4 %
//  at the full well. A read is B + P(Q)/g, of variance (P'(Q)²·Q + RN²)/g².
// ═══════════════════════════════════════════════════════════════════════════════

namespace
{
    constexpr double kFullWell = 20000.0;                  // e-
    constexpr double kC2 = -2.0e-6, kC3 = -1.0e-11;        // 1/e-, 1/e-²
    constexpr double kGain = 2.0, kRN = 10.0, kBase = 1000.0;
    constexpr double kS = 2000.0;                          // e- per frame: the well is full at frame 10
    constexpr size_t kW = 128, kH = 128;
    constexpr double kNPix  = static_cast<double>(kW * kH);
    constexpr double kSigma = 5.0;                         // tolerance in sigma

    double response(double Q, double fullWell = kFullWell)
    {
        const double q = std::min(Q, fullWell);
        return q * (1.0 + q * (kC2 + q * kC3));
    }

    double responseSlope(double Q)   // P'(Q), 0 above the full well
    {
        return (Q < kFullWell) ? 1.0 + Q * (2.0 * kC2 + 3.0 * kC3 * Q) : 0.0;
    }

    // Full well without a model, frame time 1 s.
    PixelConfig makeWellCfg()
    {
        PixelConfig cfg = makeQuietPixelCfg();
        cfg.set("gain",       kGain);
        cfg.set("readNoise",  kRN);
        cfg.set("baseline",   kBase);
        cfg.set("frame_time", 1.0);
        cfg.set("fullWell",   kFullWell);
        return cfg;
    }

    // The same with the polynomial model.
    PixelConfig makeNonLinearCfg()
    {
        PixelConfig cfg = makeWellCfg();
        cfg.set("nlModel", std::string("polynomial"));
        cfg.set("nl0",     kC2);
        cfg.set("nl1",     kC3);
        return cfg;
    }

    const std::valarray<uint16_t>& uint16Data(const std::shared_ptr<DSL::FITScube>& img)
    {
        const std::valarray<uint16_t>* d = img ? img->GetData<uint16_t>() : nullptr;
        if (!d)
            throw std::runtime_error("output is not a uint16 image");
        return *d;
    }

    // Mean and population variance of plane p, of n pixels.
    std::pair<double, double> planeStats(const std::valarray<uint16_t>& d, size_t p, size_t n)
    {
        double sum = 0.0;
        for (size_t i = p * n; i < (p + 1) * n; ++i)
            sum += static_cast<double>(d[i]);
        const double mean = sum / static_cast<double>(n);

        double ss = 0.0;
        for (size_t i = p * n; i < (p + 1) * n; ++i)
        {
            const double x = static_cast<double>(d[i]) - mean;
            ss += x * x;
        }
        return {mean, ss / static_cast<double>(n)};
    }

    // Tolerance on the mean of n pixels: variance not reduced by the non-linearity, uint16
    // rounding, and P(mean Q) instead of the mean of P(Q).
    double meanTol(double Q, double n = kNPix)
    {
        return kSigma * std::sqrt(((Q + 2.0 * kRN * kRN) / (kGain * kGain) + 1.0 / 12.0) / n) + 0.1;
    }
}

TEST(SCA, NoFullWellNorNonLinearityByDefault)
{
    nisp::SCA sca(16, 16);
    sca.initialize(makeQuietPixelCfg());
    for (size_t i = 0; i < sca.size(); ++i)
        ASSERT_TRUE(std::isinf(sca.at(i).getFullWell()));

    const auto  result = sca.expose(makeUniformFlux(16, 16, 10.0), nisp::macc{2, 1, 0});
    const auto& hdu    = result->HDU();
    for (const char* key : {"PIXFW", "NLMODEL", "PIXNL0", "PIXNL1", "PIXNL2"})
        EXPECT_FALSE(hdu.Exists(key)) << key;
}

TEST(SCA, MaccReadsFollowTheNonLinearityThenTheFullWell)
{
    // 15 reads of 2 ke- each: reads 1–9 are below the full well, 11–15 above. Read 10 is
    // at the full well, where the shot noise is cut on one side: not checked.
    constexpr size_t nReads = 15;
    nisp::SCA sca(kW, kH);
    sca.initialize(makeNonLinearCfg());

    const auto result = sca.expose(makeUniformFlux(kW, kH, kS), nisp::macc{nReads, 1, 0});
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->Size(3), nReads);
    const auto& data = uint16Data(result);

    for (size_t k = 1; k <= nReads; ++k)
    {
        const double Q = static_cast<double>(k) * kS;
        if (std::abs(Q - kFullWell) < 10.0 * std::sqrt(Q))
            continue;

        const double expectedMean = kBase + response(Q) / kGain;
        const double slope        = responseSlope(Q);
        const double expectedVar  = (slope * slope * Q + kRN * kRN) / (kGain * kGain) + 1.0 / 12.0;

        const auto [mean, var] = planeStats(data, k - 1, kW * kH);
        EXPECT_NEAR(mean, expectedMean, kSigma * std::sqrt(expectedVar / kNPix) + 0.1) << "read " << k;
        EXPECT_NEAR(var,  expectedVar,  (kSigma * std::sqrt(2.0 / kNPix) + 0.01) * expectedVar)
            << "read " << k << (Q > kFullWell ? " (above the full well: read noise only)" : " (compressed shot noise)");
    }
}

TEST(SCA, FluxCdsAndWindowFollowTheNonLinearityThenTheFullWell)
{
    // One step of N frames gives the same signal as N reads of a ramp: the response only
    // depends on the accumulated charge. 8 frames: 16 ke-, below the full well; 14: above.
    nisp::SCA sca(kW, kH);
    sca.initialize(makeNonLinearCfg());
    const auto flux = makeUniformFlux(kW, kH, kS);

    for (double nFrames : {8.0, 14.0})
    {
        const double Q        = nFrames * kS;
        const double expected = kBase + response(Q) / kGain;

        const auto single = sca.expose(flux, nisp::readout{nisp::mode::FLUX, std::nullopt, nFrames});
        EXPECT_NEAR(planeStats(uint16Data(single), 0, kW * kH).first, expected, meanTol(Q)) << "FLUX, " << nFrames << " frames";

        constexpr size_t w = 64, h = 48;
        const auto window = sca.expose(flux, nFrames, 32, 16, w, h);
        EXPECT_NEAR(planeStats(uint16Data(window), 0, w * h).first, expected, meanTol(Q, static_cast<double>(w * h)))
            << "window, " << nFrames << " frames";

        // Read N − read 1: the baseline cancels, each read has its own non-linearity.
        const auto cds = sca.expose(flux, nisp::readout{nisp::mode::CDS, std::nullopt, nFrames});
        EXPECT_NEAR(planeStats(uint16Data(cds), 0, kW * kH).first, (response(Q) - response(kS)) / kGain, meanTol(Q))
            << "CDS, " << nFrames << " frames";
    }
}

TEST(SCA, SaturatedReadsGiveTheFullWellOfEachPixel)
{
    // Full well map, gain 1, no read noise, 1e5 e- per frame: every read is exactly
    // baseline + full well of its own pixel, and the CDS is 0.
    constexpr size_t W = 16, H = 16;
    auto fullWell = [](size_t x, size_t y) { return 10000.0 + static_cast<double>(x + W * y); };
    const std::string path = writeTestMap("sca_fw_16x16.fits", W, H, fullWell);

    PixelConfig cfg = makeQuietPixelCfg();
    cfg.set("baseline",   kBase);
    cfg.set("frame_time", 1.0);
    cfg.set("fullWell",   path);
    nisp::SCA sca(W, H);
    sca.initialize(cfg);
    const auto flux = makeUniformFlux(W, H, 1.0e5);

    const auto  ramp   = sca.expose(flux, nisp::macc{2, 1, 1});
    const auto  single = sca.expose(flux, nisp::readout{nisp::mode::FLUX, std::nullopt, 3.0});
    const auto  cds    = sca.expose(flux, nisp::readout{nisp::mode::CDS,  std::nullopt, 3.0});
    const auto& rampData = uint16Data(ramp);
    for (size_t y = 0; y < H; ++y)
        for (size_t x = 0; x < W; ++x)
        {
            const size_t i        = y * W + x;
            const auto   expected = static_cast<uint16_t>(kBase + fullWell(x, y));
            ASSERT_FLOAT_EQ(sca.at(x, y).getFullWell(), static_cast<float>(fullWell(x, y)));
            ASSERT_EQ(rampData[i],             expected) << "MACC read 1, x = " << x << ", y = " << y;
            ASSERT_EQ(rampData[W * H + i],     expected) << "MACC read 2, x = " << x << ", y = " << y;
            ASSERT_EQ(uint16Data(single)[i],   expected) << "FLUX, x = " << x << ", y = " << y;
            ASSERT_EQ(uint16Data(cds)[i],      0u)       << "CDS, x = " << x << ", y = " << y;
        }

    constexpr size_t c0 = 3, r0 = 5, w = 8, h = 6;
    const auto window = sca.expose(flux, 3.0, c0, r0, w, h);
    for (size_t y = 0; y < h; ++y)
        for (size_t x = 0; x < w; ++x)
            ASSERT_EQ(uint16Data(window)[y * w + x], static_cast<uint16_t>(kBase + fullWell(c0 + x, r0 + y)))
                << "window, x = " << x << ", y = " << y;

    EXPECT_EQ(ramp->HDU().GetValueForKey("PIXFW"), "sca_fw_16x16.fits");
    std::filesystem::remove(path);
}

TEST(SCA, NonLinearityParametersAreReadFromMaps)
{
    constexpr size_t W = 16, H = 16;
    auto c2 = [](size_t x, size_t y) { return -1.0e-6 * (1.0 + static_cast<double>(x + W * y) / 256.0); };
    const std::string path = writeTestMap("sca_nl0_16x16.fits", W, H, c2);

    PixelConfig cfg = makeNonLinearCfg();
    cfg.set("nl0", path);
    nisp::SCA sca(W, H);
    sca.initialize(cfg);

    for (size_t y = 0; y < H; ++y)
        for (size_t x = 0; x < W; ++x)
        {
            ASSERT_FLOAT_EQ(sca.at(x, y).getFullWell(), static_cast<float>(kFullWell));
            ASSERT_FLOAT_EQ(sca.at(x, y).getNlParam(0), static_cast<float>(c2(x, y))) << "x = " << x << ", y = " << y;
            ASSERT_FLOAT_EQ(sca.at(x, y).getNlParam(1), static_cast<float>(kC3));
            ASSERT_FLOAT_EQ(sca.at(x, y).getNlParam(2), 0.0f);
        }
    std::filesystem::remove(path);
}

TEST(SCA, NonLinearityCanBeSetWithAFunctor)
{
    nisp::SCA sca(kW, kH);
    sca.initialize(makeNonLinearCfg());
    const auto   flux = makeUniformFlux(kW, kH, kS);
    const double Q    = 4.0 * kS;   // 8 ke-
    auto read = [&sca, &flux]() { return sca.expose(flux, nisp::readout{nisp::mode::FLUX, std::nullopt, 4.0}); };

    auto result = read();
    EXPECT_EQ(result->HDU().GetValueForKey("NLMODEL"), "polynomial");
    EXPECT_NEAR(planeStats(uint16Data(result), 0, kW * kH).first, kBase + response(Q) / kGain, meanTol(Q));

    // Captured by value: the SCA keeps its own copy of the table.
    const double a = 30000.0;
    {
        const std::vector<double> table{a};
        sca.setNonLinearity([table](std::span<const float>, double q) { return table[0] * (1.0 - std::exp(-q / table[0])); });
    }
    result = read();
    EXPECT_EQ(result->HDU().GetValueForKey("NLMODEL"), "custom");
    EXPECT_NEAR(planeStats(uint16Data(result), 0, kW * kH).first, kBase + a * (1.0 - std::exp(-Q / a)) / kGain, meanTol(Q));

    // Empty functor: linear.
    sca.setNonLinearity(nonlinearity::functor{});
    result = read();
    EXPECT_EQ(result->HDU().GetValueForKey("NLMODEL"), "linear");
    EXPECT_NEAR(planeStats(uint16Data(result), 0, kW * kH).first, kBase + Q / kGain, meanTol(Q));
}

TEST(SCA, InvalidFullWellOrNonLinearityThrows)
{
    const nonlinearity::functor identity   = [](std::span<const float>, double Q) { return Q; };
    const nonlinearity::functor decreasing = [](std::span<const float>, double Q) { return Q - Q * Q / kFullWell; };
    const nonlinearity::functor notANumber = [](std::span<const float>, double) { return std::numeric_limits<double>::quiet_NaN(); };

    // Full well.
    for (double fullWell : {0.0, -1.0})
    {
        PixelConfig cfg = makeQuietPixelCfg();
        cfg.set("fullWell", fullWell);
        nisp::SCA sca(4, 4);
        EXPECT_THROW(sca.initialize(cfg), std::invalid_argument) << "fullWell = " << fullWell;
    }

    // Model without a full well, from the configuration or from a functor.
    {
        PixelConfig cfg = makeQuietPixelCfg();
        cfg.set("nlModel", std::string("polynomial"));
        nisp::SCA sca(4, 4);
        EXPECT_THROW(sca.initialize(cfg), std::invalid_argument);

        nisp::SCA quiet(4, 4);
        quiet.initialize(makeQuietPixelCfg());
        EXPECT_THROW(quiet.setNonLinearity(identity), std::invalid_argument);
        EXPECT_NO_THROW(quiet.setNonLinearity(nonlinearity::functor{}));
    }

    // Unknown model.
    {
        PixelConfig cfg = makeWellCfg();
        cfg.set("nlModel", std::string("unknown"));
        nisp::SCA sca(4, 4);
        EXPECT_THROW(sca.initialize(cfg), std::invalid_argument);
    }

    // Response that decreases before the full well: P is maximum at 10 ke-.
    {
        PixelConfig cfg = makeNonLinearCfg();
        cfg.set("nl0", -5.0e-5);
        nisp::SCA sca(4, 4);
        EXPECT_THROW(sca.initialize(cfg), std::invalid_argument);
    }

    // A rejected functor leaves the model in place.
    nisp::SCA sca(4, 4);
    sca.initialize(makeNonLinearCfg());
    EXPECT_THROW(sca.setNonLinearity(decreasing), std::invalid_argument);
    EXPECT_THROW(sca.setNonLinearity(notANumber), std::invalid_argument);
    EXPECT_EQ(sca.expose(makeUniformFlux(4, 4, 10.0), nisp::macc{2, 1, 0})->HDU().GetValueForKey("NLMODEL"), "polynomial");

    // New parameters are checked with the model already set.
    EXPECT_THROW(sca.initialize(PixelConfig{{"nl0", -5.0e-5}}), std::invalid_argument);
}

TEST(SCA, OutputHeaderRecordsTheFullWellAndTheNonLinearity)
{
    nisp::SCA sca(4, 4);
    sca.initialize(makeNonLinearCfg());
    const auto  result = sca.expose(makeUniformFlux(4, 4, 10.0), nisp::macc{2, 1, 0});
    const auto& hdu    = result->HDU();

    EXPECT_NEAR (hdu.GetDoubleValueForKey("PIXFW"), kFullWell, 1e-9);
    EXPECT_EQ   (hdu.GetValueForKey("NLMODEL"), "polynomial");
    EXPECT_NEAR (hdu.GetDoubleValueForKey("PIXNL0"), kC2, 1e-6 * std::abs(kC2));   // not written as 0
    EXPECT_NEAR (hdu.GetDoubleValueForKey("PIXNL1"), kC3, 1e-6 * std::abs(kC3));
    EXPECT_FALSE(hdu.Exists("PIXNL2"));   // not given
}

TEST(SCA, DpuSlopeFollowsTheNonLinearity)
{
    // MACC(4,4,4), 200 e- per frame, gain 1 and no read noise: 5.6 ke- at the last read.
    // The DPU fits a straight line, so the slope is that of the chord between the first and
    // the last group: −1.2 % here. It is compared with the slope of the same linear pixels.
    constexpr size_t nG = 4, nR = 4, nD = 4;
    constexpr double S  = 200.0;
    auto dpuMean = [](PixelConfig cfg)
    {
        cfg.set("gain",         1.0);
        cfg.set("readNoise",    0.0);
        cfg.set("exposureMode", std::string("DPU"));
        nisp::SCA sca(kW, kH);
        sca.initialize(cfg);
        const auto result = sca.expose(makeUniformFlux(kW, kH, S), nisp::macc{nG, nR, nD});
        return planeStats(uint16Data(result), 0, kW * kH).first;
    };

    // Average charge and average response of the reads of group g.
    auto group = [](size_t g)
    {
        double Q = 0.0, P = 0.0;
        for (size_t r = 1; r <= nR; ++r)
        {
            const double q = S * static_cast<double>(g * (nR + nD) + r);
            Q += q / static_cast<double>(nR);
            P += response(q) / static_cast<double>(nR);
        }
        return std::make_pair(Q, P);
    };
    const auto [Q0, P0] = group(0);
    const auto [Q3, P3] = group(nG - 1);

    const double linear    = dpuMean(makeWellCfg());
    const double nonLinear = dpuMean(makeNonLinearCfg());
    EXPECT_NEAR(linear, S * static_cast<double>(nR + nD), 2.0);
    EXPECT_NEAR(nonLinear / linear, (P3 - P0) / (Q3 - Q0), 1.5e-3);
}

TEST(SCA, DpuDoesNotFlagAPlateauBelowTheAdcLimit)
{
    // The NISP DPU only tests 65535: a pixel whose full well is below is not flagged, and
    // its plateau is fitted as signal. Full well 3 ke-, 1600 e- per group: the group
    // averages are 500, 2100, 3000 and 3000 e-, that is a slope of √((1600² + 900²) / 3).
    PixelConfig cfg = makeWellCfg();
    cfg.set("gain",         1.0);
    cfg.set("readNoise",    0.0);
    cfg.set("fullWell",     3000.0);
    cfg.set("exposureMode", std::string("DPU"));
    nisp::SCA sca(kW, kH);
    sca.initialize(cfg);

    const auto  result = sca.expose(makeUniformFlux(kW, kH, 200.0), nisp::macc{4, 4, 4});
    const auto& data   = uint16Data(result);
    const auto flagged = std::count(std::begin(data), std::begin(data) + kW * kH, uint16_t(65535));   // slope plane
    EXPECT_EQ(flagged, 0) << "a plateau below the ADC limit must not be flagged";

    const double expected = std::sqrt((1600.0 * 1600.0 + 900.0 * 900.0) / 3.0);   // 1059.9 ADU, instead of 1600
    EXPECT_NEAR(planeStats(data, 0, kW * kH).first, expected, 0.01 * expected);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (i) Photon transfer curve to check the gain and the read noise
// ═══════════════════════════════════════════════════════════════════════════════

namespace
{
    // Gain and read noise of a photon transfer curve.
    struct PtcFit
    {
        double gain,      gainErr;         // e-/ADU
        double readNoise, readNoiseErr;    // e-
        double chi2;                       // (number of points − 2) degrees of freedom
    };

    // Straight line Var = mean / g + 2·RN² / g² + quant fitted to the points (mean, var), in ADU
    // and ADU², each measured over n pixels: the error on a variance is Var·√(2/n).
    // quant is the rounding noise (ADU²). NaN read noise if the intercept is below it.
    PtcFit fitPtc(const std::vector<double>& mean, const std::vector<double>& var, double n, double quant)
    {
        double sw = 0.0, swx = 0.0, swxx = 0.0, swy = 0.0, swxy = 0.0;
        for (size_t i = 0; i < mean.size(); ++i)
        {
            const double w = n / (2.0 * var[i] * var[i]);
            sw   += w;
            swx  += w * mean[i];
            swxx += w * mean[i] * mean[i];
            swy  += w * var[i];
            swxy += w * mean[i] * var[i];
        }

        // Var = slope·mean + intercept, and the covariance matrix of the fit.
        const double det       = sw * swxx - swx * swx;
        const double slope     = (sw * swxy - swx * swy) / det;
        const double intercept = (swxx * swy - swx * swxy) / det;
        const double varSlope  = sw / det;
        const double varInter  = swxx / det;
        const double covar     = -swx / det;

        PtcFit fit{};
        for (size_t i = 0; i < mean.size(); ++i)
        {
            const double r = var[i] - slope * mean[i] - intercept;
            fit.chi2 += r * r * n / (2.0 * var[i] * var[i]);
        }

        const double readVar = intercept - quant;   // 2·RN²/g²
        fit.gain         = 1.0 / slope;
        fit.gainErr      = std::sqrt(varSlope) / (slope * slope);
        fit.readNoise    = fit.gain * std::sqrt(readVar / 2.0);
        fit.readNoiseErr = fit.readNoise * std::sqrt(varInter / (4.0 * readVar * readVar)
                                                     + varSlope / (slope * slope)
                                                     - covar / (slope * readVar));
        return fit;
    }
}

TEST(SCA, PhotonTransferCurveGivesTheGainAndTheReadNoise)
{
    // Use of CDS exposures of increasing length, uniform flux, uniform pixels matrixes without correlated noise.
    // The reset level (baseline + kTC) is canceled in CDS exposure: One obtain
    //   mean signal \$ \bar{S}= \frac{S}{g} \$ with \$S\$ the signal in e-,   Var = (S + 2·RN²) / g² + 1/12 = mean / g + 2·RN² / g² + 1/12
    // (factor 1/12come from uint16 rounding). A straight line fitted to Var=f(mean) gives the gain
    // (slope 1/g) and the read noise as the intercept; the flux, QE and dark value are not used in the PTC estimate.
    // The exposures are independent, while the reads of one ramp share their shot noise.
    constexpr size_t W = 256, H = 256;
    constexpr double n      = static_cast<double>(W * H);
    constexpr double offset = 500.0;         // ADU: keeps the CDS away from the clamp at 0
    constexpr double quant  = 1.0 / 12.0;    // uint16 rounding noise, ADU²
    const auto       flux   = makeUniformFlux(W, H, 100.0);

    for (const double RN : {5.0, 10.0, 13.0, 20.0, 24.0,})   // RN e-1
    {
        std::cout<<std::endl;
        SCOPED_TRACE("read noise " + std::to_string(RN) + " e-");
        
        for (const double g : {0.5, 1.0, 1.92, 4.0})
        {
            SCOPED_TRACE("gain " + std::to_string(g) + " e-/ADU");

            nisp::SCA sca(W, H);
            sca.initialize(PixelConfig{
                {"gain", g}, {"dark", 0.01}, {"readNoise", RN}, {"QE", 0.94},
                {"KTCnoise", 48.0}, {"baseline", 1000.0}, {"frame_time", 1.0},
                {"offset", offset}, {"exposureMode", std::string("CDS")},
                {"detectorOrigin", activeOrigin()}
            });

            // 94 e- to 18.7 ke-: at most 37.9 kADU (gain 0.5, offset included), below the ADC limit.
            std::vector<double> means, vars;
            for (const double nFrames : {2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 200.0})
            {
                const auto result = sca.expose(flux, nFrames);
                ASSERT_NE(result, nullptr);

                const auto [mean, var] = planeStats(uint16Data(result), 0, W * H);
                ASSERT_GT(var, 0.0) << nFrames << " frames";
                means.push_back(mean - offset);
                vars.push_back(var);
            }

            const PtcFit fit = fitPtc(means, vars, n, quant);

            std::cout << "[SCA PTC] gain " << g << " e-/ADU & RN " << RN << " e- : "
                      << fit.gain      << " +/- " << fit.gainErr      << " e-/ADU, read noise "
                      << fit.readNoise << " +/- " << fit.readNoiseErr << " e-\n";

            // Straight line: 8 points, 6 degrees of freedom (P(chi2 > 40) = 5e-7).
            EXPECT_LT(fit.chi2, 40.0) << "the photon transfer curve is not a straight line";

            EXPECT_NEAR(fit.gain,      g,  kSigma * fit.gainErr)      << "gain (e-/ADU)";
            EXPECT_NEAR(fit.readNoise, RN, kSigma * fit.readNoiseErr) << "read noise (e-)";

            // The fit must be precise enough for these checks to mean something (0.3 % and 0.5 %).
            EXPECT_LT(fit.gainErr,      0.005 * g);
            EXPECT_LT(fit.readNoiseErr, 0.010 * RN);
        }
    }
}

TEST(SCA, DpuPhotonTransferCurveGivesTheGainAndTheReadNoise)
{
    // Same curve with the DPU: MACC(ng, 1, 0) with ng increasing, at constant flux.
    // The slope is per frame, so it does not change with ng: the signal of a point is
    // slope × (ng − 1), which is the value of plane 0 with "scale" = ng − 1. At first order it
    // is last read − first read (baseline and kTC cancel), so with S the signal in e-:
    //   mean = S / g,   Var = (S + 2·RN²) / g² + 3/12
    // 3/12: uint16 rounding of both reads by the ADC, then of the output by the DPU.
    //
    // Beyond first order, the estimator also uses the scatter of the differences (and dpu_gain,
    // here the pixel gain): a relative effect of 1 / (2·S per frame) on the variance. With
    // 470 e- per frame, it moves the measured gain by less than 0.1 %, below the fit error.
    // Plane 2 cannot replace the variance over the matrix: it is computed with dpu_gain and
    // would give it back.
    constexpr size_t W = 256, H = 256;
    constexpr double n     = static_cast<double>(W * H);
    constexpr double quant = 3.0 / 12.0;     // ADU²
    const auto       flux  = makeUniformFlux(W, H, 500.0);

    // ng − 1: 470 e- to 18.8 ke-; the last read is at most 39.5 kADU (gain 0.5), below the ADC limit.
    const std::vector<size_t> intervals{1, 2, 3, 5, 8, 12, 20, 30, 40};

    for (const double RN : {5.0, 10.0, 13.0, 20.0, 24.0})   // e-
    {
        std::cout << std::endl;
        SCOPED_TRACE("read noise " + std::to_string(RN) + " e-");

        for (const double g : {0.5, 1.0, 1.92, 4.0})
        {
            SCOPED_TRACE("gain " + std::to_string(g) + " e-/ADU");

            std::vector<double> means, vars;
            for (const size_t nInt : intervals)
            {
                nisp::SCA sca(W, H);
                sca.initialize(PixelConfig{
                    {"gain", g}, {"dark", 0.01}, {"readNoise", RN}, {"QE", 0.94},
                    {"KTCnoise", 48.0}, {"baseline", 1000.0}, {"frame_time", 1.0},
                    {"scale", static_cast<double>(nInt)}, {"offset", 0.0},
                    {"exposureMode", std::string("DPU")}, {"refPixCorr", std::string("OFF")},
                    {"detectorOrigin", activeOrigin()}
                });

                const auto result = sca.expose(flux, nisp::macc{nInt + 1, 1, 0});
                ASSERT_NE(result, nullptr);
                ASSERT_EQ(result->Size(3), 3u);

                const auto [mean, var] = planeStats(uint16Data(result), 0, W * H);
                ASSERT_GT(var, 0.0) << "MACC(" << nInt + 1 << ",1,0)";
                means.push_back(mean);
                vars.push_back(var);
            }

            const PtcFit fit = fitPtc(means, vars, n, quant);

            std::cout << "[SCA DPU PTC] gain " << g << " e-/ADU & RN " << RN << " e- : "
                      << fit.gain      << " +/- " << fit.gainErr      << " e-/ADU, read noise "
                      << fit.readNoise << " +/- " << fit.readNoiseErr << " e-\n";

            // Straight line: 9 points, 7 degrees of freedom (P(chi2 > 45) = 1e-7).
            EXPECT_LT(fit.chi2, 45.0) << "the photon transfer curve is not a straight line";

            EXPECT_NEAR(fit.gain,      g,  kSigma * fit.gainErr)      << "gain (e-/ADU)";
            EXPECT_NEAR(fit.readNoise, RN, kSigma * fit.readNoiseErr) << "read noise (e-)";

            // Precision of the fit: 0.3 % on the gain. The first point is at 470 e-, so the read
            // noise is less constrained than with the CDS: 0.6 % at 24 e-, 3.3 % at 5 e-.
            EXPECT_LT(fit.gainErr,      0.005 * g);
            EXPECT_LT(fit.readNoiseErr, 0.050 * RN);
        }
    }
}
