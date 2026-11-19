//
// Created by gillard on 2026-11-19.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_nispReadoutModes.cxx
//
// Validity of the nisp::SCA readout modes and of the flux each one yields:
// MACC, UTR, CDS, FLUX and DPU (selected with "exposureMode").
//
// Set-up, identical for every mode:
//   - full 2048 x 2048 detector,
//   - every detector effect off: gain 1 e-/ADU, QE 1, read noise 0, dark 0,
//     kTC 0, baseline 0, no persistence, no 1/f noise, DPU offset 0 / scale 1,
//   - uniform flux F = 10 ph/s/pixel, NISP frame time tf = 1.45408 s,
//   - MACC(4,16,4): 4 groups of 16 reads, 4 dropped frames between groups,
//     i.e. 76 frames from reset to last read.
//
// With these settings 1 ADU = 1 e- = 1 photon and each frame adds
// S = F * tf = 14.5408 ADU. The only noise left is photon shot noise, drawn by
// pixsim::integrate() as a Gaussian of variance = signal, so a pixel that
// integrated n frames reads n*S ± sqrt(n*S). Every check is statistical over
// the 4,194,304 pixels: means within 5 sigma of the expected value, variances
// within 5 sigma + 1 % (uint16 rounding adds 1/12 ADU² of quantisation noise).
//
// Expected results:
//   MACC  64 layers, layer 16g+f holds (20g+f+1)*S         flux = (G3-G0) / (3*20*tf)
//   UTR   76 layers, layer k holds (k+1)*S                 flux = least-squares slope
//   CDS   1 frame, last - first read = 75*S = 1090.56      flux = CDS / (75*tf)
//   FLUX  1 frame, 76*S = 1105.10                          flux = value / (76*tf)
//   DPU   1 frame, Kubik et al. 2016 slope in ADU/group,   flux = value / (20*tf)
//         (nr+nd)*S = 290.82
//
// The DPU test is the last one in the file so that, run as a single binary,
// a crash there does not hide the other results (ctest runs each test in its
// own process anyway).
//

#include <gtest/gtest.h>
#include <nisp/nispSCA.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <valarray>
#include <vector>

using namespace nisp;
using namespace phot2pix;

namespace
{
    constexpr size_t kW    = 2048;
    constexpr size_t kH    = 2048;
    constexpr size_t kNPix = kW * kH;

    constexpr double kFlux      = 10.0;      // ph/s/pixel
    constexpr double kFrameTime = 1.45408;   // s, NISP H2RG frame time
    constexpr size_t kNG = 4, kNR = 16, kND = 4;
    constexpr size_t kNFrames = kNG * (kNR + kND) - kND;   // 76 frames from reset to last read
    constexpr double kS       = kFlux * kFrameTime;        // ADU added per frame (gain 1, QE 1)

    constexpr double kK     = 5.0;          // tolerance in sigma
    constexpr double kQuant = 1.0 / 12.0;   // uint16 rounding noise, ADU²

    const nisp::macc kMacc{kNG, kNR, kND};

    // SCA with every detector effect off, in the given exposureMode.
    std::unique_ptr<nisp::SCA> makeSCA(const std::string& mode)
    {
        auto sca = std::make_unique<nisp::SCA>(kW, kH);
        sca->initialize(PixelConfig{
            {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 1.0},
            {"KTCnoise", 0.0}, {"baseline", 0.0},
            {"offset", 0.0}, {"scale", 1.0}, {"frame_time", kFrameTime},
            {"exposureMode", mode}
        });
        return sca;
    }

    struct Stats { double mean; double var; };

    // Mean and population variance of d[offset .. offset+n).
    Stats planeStats(const std::valarray<uint16_t>& d, size_t offset, size_t n)
    {
        double sum = 0.0;
        for (size_t i = offset; i < offset + n; ++i)
            sum += static_cast<double>(d[i]);
        const double mean = sum / static_cast<double>(n);

        double ss = 0.0;
        for (size_t i = offset; i < offset + n; ++i)
        {
            const double x = static_cast<double>(d[i]) - mean;
            ss += x * x;
        }
        return {mean, ss / static_cast<double>(n)};
    }

    // Tolerance on the mean of kNPix pixels of variance var (+0.01 ADU for rounding / clamping at 0).
    double meanTol(double var) { return kK * std::sqrt(var / kNPix) + 0.01; }

    // Tolerance on the sample variance of kNPix Gaussian pixels (+1 % for the quantisation approximation).
    double varTol(double var)  { return kK * var * std::sqrt(2.0 / kNPix) + 0.01 * var; }

    // Check that the plane holds n frames of signal: mean n*S, variance n*S (+ quantisation).
    void expectFrames(const Stats& s, double nFrames, const std::string& what)
    {
        const double signal = nFrames * kS;
        EXPECT_NEAR(s.mean, signal,          meanTol(signal))         << what << ": mean";
        EXPECT_NEAR(s.var,  signal + kQuant, varTol(signal + kQuant)) << what << ": variance (shot noise)";
    }

    const std::valarray<uint16_t>& uint16Data(const std::shared_ptr<DSL::FITScube>& img)
    {
        const std::valarray<uint16_t>* d = img->GetData<uint16_t>();
        if (!d)
            throw std::runtime_error("output is not a uint16 image");
        return *d;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  MACC(4,16,4): 64 cumulative reads, group-average flux estimate
// ═══════════════════════════════════════════════════════════════════════════════

TEST(NispReadoutModes, MACC)
{
    auto sca = makeSCA("MACC");
    const std::valarray<double> flux(kFlux, kNPix);

    auto result = sca->expose(flux, kMacc);
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->Size(1), kW);
    ASSERT_EQ(result->Size(2), kH);
    ASSERT_EQ(result->Size(3), kNG * kNR);
    const auto& data = uint16Data(result);

    // Layer 16g+f is read f+1 frames into group g, i.e. g*(nr+nd) + f + 1 frames after reset.
    std::vector<double> groupMean(kNG, 0.0);
    for (size_t g = 0; g < kNG; ++g)
        for (size_t f = 0; f < kNR; ++f)
        {
            const size_t layer   = g * kNR + f;
            const size_t nFrames = g * (kNR + kND) + f + 1;
            const Stats  s       = planeStats(data, layer * kNPix, kNPix);
            expectFrames(s, static_cast<double>(nFrames), "layer " + std::to_string(layer));
            groupMean[g] += s.mean / static_cast<double>(kNR);
        }

    // Group averages are (nr+nd)*tf apart: G_g = S * (20g + 8.5).
    for (size_t g = 0; g < kNG; ++g)
        EXPECT_NEAR(groupMean[g], kS * (static_cast<double>(g * (kNR + kND)) + 0.5 * (kNR + 1)), 0.05)
            << "group " << g << " average";

    const double fluxEstimate = (groupMean.back() - groupMean.front()) /
                                (static_cast<double>((kNG - 1) * (kNR + kND)) * kFrameTime);
    EXPECT_NEAR(fluxEstimate, kFlux, 1e-3 * kFlux) << "MACC flux estimate (ph/s)";

    EXPECT_EQ  (result->HDU().GetValueForKey("MODE"), "MACC");
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("NG"),      kNG, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("NR"),      kNR, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("ND"),      kND, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("GAIN"),    1.0, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("RN"),      0.0, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("FRTIME"),  kFrameTime, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("INTTIME"), kNFrames * kFrameTime, 1e-6);
    // Reset to last read: the drops after the last group are never executed.
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("EXPTIME"), (kNFrames + kND ) * kFrameTime, 1e-6);
    EXPECT_FALSE(result->HDU().Exists("PERSIST"));
}

// ═══════════════════════════════════════════════════════════════════════════════
//  UTR over the same 76 frames: every frame read, least-squares slope
// ═══════════════════════════════════════════════════════════════════════════════

TEST(NispReadoutModes, UTR)
{
    auto sca = makeSCA("UTR");
    const std::valarray<double> flux(kFlux, kNPix);

    // UTR is selected by an exposure time: one read per frame, no drops.
    auto result = sca->expose(flux, kNFrames * kFrameTime);
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->Size(1), kW);
    ASSERT_EQ(result->Size(2), kH);
    ASSERT_EQ(result->Size(3), kNFrames);
    const auto& data = uint16Data(result);

    // Layer k holds k+1 frames. Least-squares slope of the mean ramp against
    // t_k = (k+1)*tf (= mean of the per-pixel slopes, the fit being linear).
    double st = 0.0, sy = 0.0, stt = 0.0, sty = 0.0;
    for (size_t k = 0; k < kNFrames; ++k)
    {
        const Stats s = planeStats(data, k * kNPix, kNPix);
        expectFrames(s, static_cast<double>(k + 1), "read " + std::to_string(k));

        const double t = static_cast<double>(k + 1) * kFrameTime;
        st  += t;      sy  += s.mean;
        stt += t * t;  sty += t * s.mean;
    }
    const double n            = static_cast<double>(kNFrames);
    const double fluxEstimate = (n * sty - st * sy) / (n * stt - st * st);
    EXPECT_NEAR(fluxEstimate, kFlux, 1e-3 * kFlux) << "UTR flux estimate (ph/s)";

    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("NG"),      kNFrames, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("NR"),      1.0,      1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("ND"),      0.0,      1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("INTTIME"), kNFrames * kFrameTime, 1e-6);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  CDS over the MACC(4,16,4) duration: last read - first read
// ═══════════════════════════════════════════════════════════════════════════════

TEST(NispReadoutModes, CDS)
{
    auto sca = makeSCA("CDS");
    const std::valarray<double> flux(kFlux, kNPix);

    auto result = sca->expose(flux, kMacc);
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->GetDimension(), 2u);
    ASSERT_EQ(result->Size(1), kW);
    ASSERT_EQ(result->Size(2), kH);
    const auto& data = uint16Data(result);

    // Reads at frame 1 and frame 76: the difference holds 75 frames.
    const Stats s = planeStats(data, 0, kNPix);
    expectFrames(s, static_cast<double>(kNFrames - 1), "CDS");

    const double fluxEstimate = s.mean / (static_cast<double>(kNFrames - 1) * kFrameTime);
    EXPECT_NEAR(fluxEstimate, kFlux, 1e-3 * kFlux) << "CDS flux estimate (ph/s)";

    EXPECT_EQ(result->HDU().GetValueForKey("MODE"), "CDS");
}

// ═══════════════════════════════════════════════════════════════════════════════
//  FLUX over the MACC(4,16,4) duration: single frame integrated for 76 frames
// ═══════════════════════════════════════════════════════════════════════════════

TEST(NispReadoutModes, FLUX)
{
    auto sca = makeSCA("FLUX");
    const std::valarray<double> flux(kFlux, kNPix);

    auto result = sca->expose(flux, kMacc);
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->GetDimension(), 2u);
    ASSERT_EQ(result->Size(1), kW);
    ASSERT_EQ(result->Size(2), kH);
    const auto& data = uint16Data(result);

    const Stats s = planeStats(data, 0, kNPix);
    expectFrames(s, static_cast<double>(kNFrames), "FLUX");

    const double fluxEstimate = s.mean / (static_cast<double>(kNFrames) * kFrameTime);
    EXPECT_NEAR(fluxEstimate, kFlux, 1e-3 * kFlux) << "FLUX flux estimate (ph/s)";

    EXPECT_EQ  (result->HDU().GetValueForKey("MODE"), "FLUX");
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("INTTIME"), kNFrames * kFrameTime, 1e-6);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  DPU: onboard slope estimator (Kubik et al. 2016) on MACC(4,16,4)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(NispReadoutModes, DPU)
{
    auto sca = makeSCA("DPU");
    const std::valarray<double> flux(kFlux, kNPix);

    auto result = sca->expose(flux, kMacc);
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->GetDimension(), 2u);
    ASSERT_EQ(result->Size(1), kW);
    ASSERT_EQ(result->Size(2), kH);
    const auto& data = uint16Data(result);

    // No pixel may come out empty or saturated (would flag a NaN / failed fit).
    size_t nBad = 0;
    for (size_t i = 0; i < kNPix; ++i)
        if (data[i] == 0 || data[i] == 65535)
            ++nBad;
    EXPECT_EQ(nBad, 0u) << "pixels at 0 or 65535";

    // Estimator output: signal per group interval, (nr+nd)*S = 290.816 ADU.
    // Its non-linear bias is about -0.1 ADU and its per-pixel scatter about
    // 9.4 ADU (independent re-implementation of Ramp::getSlope() on the same
    // Gaussian ramps); both are well inside the tolerances below.
    const double perGroup = static_cast<double>(kNR + kND) * kS;
    const Stats  s        = planeStats(data, 0, kNPix);
    EXPECT_NEAR(s.mean, perGroup, 5e-3 * perGroup) << "DPU mean (ADU/group)";
    EXPECT_NEAR(std::sqrt(s.var), 9.4, 0.25 * 9.4) << "DPU per-pixel scatter (ADU)";

    const double fluxEstimate = s.mean / (static_cast<double>(kNR + kND) * kFrameTime);
    EXPECT_NEAR(fluxEstimate, kFlux, 5e-3 * kFlux) << "DPU flux estimate (ph/s)";

    EXPECT_EQ  (result->HDU().GetValueForKey("MODE"), "DPU");
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("NG"), kNG, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("NR"), kNR, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("ND"), kND, 1e-9);
}
