//
// Created by gillard on 2026-07-03.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_nispOneOverFNoise.cxx
//

#include <gtest/gtest.h>
#include <nisp/nispOneOverFNoise.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <valarray>
#include <vector>

using namespace nisp;
using namespace phot2pix;

// ── Helpers ──────────────────────────────────────────────────────────────────
//
// NOTE on tolerances used throughout this file: OneOverFNoise::generate()
// builds each channel as a sum of AR(1) poles with corner frequencies
// log-spaced down to f_min = 1/(2*L), L = 2048*64 samples per channel and per
// physical frame — i.e. the slowest pole's correlation time is comparable to
// (or longer than) the *entire* observed series. This is
// the defining feature of 1/f-type noise, but it also means a single finite
// realization converges much more slowly (and with much fatter statistical
// tails) than i.i.d. Gaussian noise of the kind used in gTest_pixels.cxx.
// Tolerances below are therefore deliberately broad (factor-of-3 bands,
// direction-only comparisons) rather than tight k*sigma bounds — they are
// meant to catch gross implementation bugs (wrong scale, wrong sign, noise
// not actually colored), not to pin down precise statistical moments.

// Reconstructs one channel's row-major AR(1) time axis from a 2D noise map
// (t = row*colWidth + col_offset_within_channel), matching _fillChannel()'s
// scatter order exactly — a naive column-stride read would incorrectly cross
// into a neighboring channel's columns.
static std::vector<double> extractChannelSeries(const std::valarray<double>& noise,
                                                 size_t mapWidth, size_t height,
                                                 size_t colStart, size_t colWidth)
{
    std::vector<double> series;
    series.reserve(height * colWidth);
    for (size_t r = 0; r < height; ++r)
        for (size_t c = 0; c < colWidth; ++c)
            series.push_back(noise[r * mapWidth + colStart + c]);
    return series;
}

static double mean(const std::vector<double>& v)
{
    double s = 0.0;
    for (double x : v) s += x;
    return s / static_cast<double>(v.size());
}

static double rms(const std::vector<double>& v)
{
    double s = 0.0;
    for (double x : v) s += x * x;
    return std::sqrt(s / static_cast<double>(v.size()));
}

static double lag1Autocorr(const std::vector<double>& v)
{
    const double m = mean(v);
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i + 1 < v.size(); ++i)
        num += (v[i] - m) * (v[i + 1] - m);
    for (double x : v) den += (x - m) * (x - m);
    return num / den;
}

static double crossCorr(const std::vector<double>& a, const std::vector<double>& b)
{
    const size_t n = std::min(a.size(), b.size());
    const double ma = mean(a), mb = mean(b);
    double num = 0.0, da = 0.0, db = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        num += (a[i] - ma) * (b[i] - mb);
        da  += (a[i] - ma) * (a[i] - ma);
        db  += (b[i] - mb) * (b[i] - mb);
    }
    return num / std::sqrt(da * db);
}

// First difference: v[i+1]-v[i]. A function of one series alone stays
// independent of a function of another independent series alone, so this
// preserves the independence question we're testing — but it suppresses the
// shared slow/low-frequency structure that otherwise makes a naive
// cross-correlation of two 1/f-like series a classic "spurious correlation
// between independent random walks" trap (two genuinely independent but
// strongly autocorrelated series can show large apparent correlation by pure
// chance in any single finite realization).
static std::vector<double> firstDifference(const std::vector<double>& v)
{
    std::vector<double> d;
    if (v.size() < 2) return d;
    d.reserve(v.size() - 1);
    for (size_t i = 1; i < v.size(); ++i)
        d.push_back(v[i] - v[i - 1]);
    return d;
}

// Broad sanity band: RMS should be within a factor of `factor` of the target
// amplitude — loose on purpose (see note above).
static void expectRmsWithinBand(double rmsValue, double amplitude, double factor = 3.0)
{
    EXPECT_GT(rmsValue, amplitude / factor);
    EXPECT_LT(rmsValue, amplitude * factor);
}

// The matrix starts at detector pixel (0, 0) by default. Tests that check where the
// noise lands on the detector (partial channels, row read times) place their W × H
// matrix at the centre of the 2048 × 2048 detector explicitly.
static void placeAtCentre(OneOverFNoise& n, size_t W, size_t H)
{
    nisp::Geometry g;
    g.origin = std::make_pair((g.width - W) / 2, (g.height - H) / 2);
    n.setGeometry(g);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Construction / validation
// ═══════════════════════════════════════════════════════════════════════════════

TEST(OneOverFNoise, DefaultConstructedDisabled)
{
    OneOverFNoise n;
    EXPECT_FALSE(n.isEnabled());
    EXPECT_DOUBLE_EQ(n.getAmplitude(), 0.0);
}

TEST(OneOverFNoise, ExplicitConstructorStoresAllFields)
{
    OneOverFNoise n(7.0, 1.5, 4, 32, 8);
    EXPECT_TRUE(n.isEnabled());
    EXPECT_DOUBLE_EQ(n.getAmplitude(),     7.0);
    EXPECT_DOUBLE_EQ(n.getSlope(),         1.5);
    EXPECT_EQ(n.getNChannels(),            4u);
    EXPECT_EQ(n.getColPerChannel(),        32u);
    EXPECT_EQ(n.getNPoles(),               8u);
}

TEST(OneOverFNoise, DefaultParameterValues)
{
    OneOverFNoise n(5.0);
    EXPECT_DOUBLE_EQ(n.getSlope(),  1.0);
    EXPECT_EQ(n.getNChannels(),     0u);
    EXPECT_EQ(n.getColPerChannel(), 64u);
    EXPECT_EQ(n.getNPoles(),        16u);
}

TEST(OneOverFNoise, InvalidColPerChannelThrows)
{
    EXPECT_THROW(OneOverFNoise(5.0, 1.0, 0, 0, 16), std::invalid_argument);
}

TEST(OneOverFNoise, InvalidNPolesThrows)
{
    EXPECT_THROW(OneOverFNoise(5.0, 1.0, 0, 64, 0), std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  initialize(cfg)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(OneOverFNoise, InitializeNoKeyIsNoop)
{
    OneOverFNoise n;
    PixelConfig cfg;
    n.initialize(cfg);
    EXPECT_FALSE(n.isEnabled());
}

TEST(OneOverFNoise, InitializeScalarAmplitudeOnly)
{
    OneOverFNoise n;
    PixelConfig cfg;
    cfg.set("oneOverF", 7.0);
    n.initialize(cfg);

    EXPECT_DOUBLE_EQ(n.getAmplitude(), 7.0);
    EXPECT_DOUBLE_EQ(n.getSlope(),     1.0);
    EXPECT_EQ(n.getColPerChannel(),    64u);
    EXPECT_EQ(n.getNPoles(),           16u);
}

TEST(OneOverFNoise, InitializeObjectFormAppliesAllFields)
{
    OneOverFNoise n;
    PixelConfig cfg;

    nlohmann::json obj;
    obj["amplitude"]     = 9.0;
    obj["slope"]         = 1.6;
    obj["nChannels"]     = 4;
    obj["colPerChannel"] = 32;
    obj["nPoles"]        = 8;
    cfg.set("oneOverF", obj);

    n.initialize(cfg);

    EXPECT_DOUBLE_EQ(n.getAmplitude(),  9.0);
    EXPECT_DOUBLE_EQ(n.getSlope(),      1.6);
    EXPECT_EQ(n.getNChannels(),         4u);
    EXPECT_EQ(n.getColPerChannel(),     32u);
    EXPECT_EQ(n.getNPoles(),            8u);
}

TEST(OneOverFNoise, InvalidColPerChannelInObjectFormThrows)
{
    OneOverFNoise n;
    PixelConfig cfg;
    nlohmann::json obj;
    obj["amplitude"]     = 5.0;
    obj["colPerChannel"] = 0;
    cfg.set("oneOverF", obj);

    EXPECT_THROW(n.initialize(cfg), std::invalid_argument);
}

TEST(OneOverFNoise, InvalidNPolesInObjectFormThrows)
{
    OneOverFNoise n;
    PixelConfig cfg;
    nlohmann::json obj;
    obj["amplitude"] = 5.0;
    obj["nPoles"]    = 0;
    cfg.set("oneOverF", obj);

    EXPECT_THROW(n.initialize(cfg), std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Channel layout (exercised indirectly through generate())
// ═══════════════════════════════════════════════════════════════════════════════

TEST(OneOverFNoise, UniformLayoutWidthExactMultiple)
{
    // nCh = round(256/64) = 4, 4*64 == 256 -> uniform layout.
    OneOverFNoise n(5.0, 1.0, 0, 64, 16);
    std::valarray<double> out;
    EXPECT_NO_THROW(out = n.generate(256, 10));
    EXPECT_EQ(out.size(), 256u * 10u);
}

TEST(OneOverFNoise, ReferencePixelLayoutWidthPlus8)
{
    // nCh = round(2040/64) = 32, 32*64 = 2048 = 2040+8 -> reference-pixel layout.
    OneOverFNoise n(5.0, 1.0, 0, 64, 16);
    std::valarray<double> out;
    EXPECT_NO_THROW(out = n.generate(2040, 10));
    EXPECT_EQ(out.size(), 2040u * 10u);
}

TEST(OneOverFNoise, AnyMatrixSizeUpToTheDetectorIsValid)
{
    // The layout comes from the 2048-wide detector, not from the matrix: a
    // 1-column, a 2044-wide or a 100-wide matrix is just a part of it.
    OneOverFNoise n(5.0, 1.0, 0, 64, 16);
    for (size_t w : {1u, 100u, 2044u, 2048u})
    {
        std::valarray<double> out;
        EXPECT_NO_THROW(out = n.generate(w, 4)) << "width " << w;
        EXPECT_EQ(out.size(), w * 4u);
    }
}

TEST(OneOverFNoise, MatrixLargerThanDetectorThrows)
{
    OneOverFNoise n(5.0, 1.0, 0, 64, 16);
    EXPECT_THROW(n.generate(2049, 10), std::invalid_argument);
    EXPECT_THROW(n.generate(64, 4096),  std::invalid_argument);
}

TEST(OneOverFNoise, ChannelLayoutMustMatchTheDetector)
{
    // colPerChannel must divide the detector width (2048), and an explicit
    // nChannels must equal width / colPerChannel.
    EXPECT_THROW(OneOverFNoise(5.0, 1.0, 0, 100, 16).generate(64, 4), std::invalid_argument);
    EXPECT_THROW(OneOverFNoise(5.0, 1.0, 5, 64,  16).generate(64, 4), std::invalid_argument);
    EXPECT_NO_THROW(OneOverFNoise(5.0, 1.0, 32, 64, 16).generate(64, 4));
}

TEST(OneOverFNoise, DisabledGenerateNeverBuildsLayoutOrThrows)
{
    // Disabled (amplitude==0): the isEnabled() early-return happens before the
    // matrix is placed, so a matrix that WOULD throw if enabled (see
    // MatrixLargerThanDetectorThrows) must not throw here.
    OneOverFNoise n; // disabled
    std::valarray<double> out;
    EXPECT_NO_THROW(out = n.generate(4096, 10));
    ASSERT_EQ(out.size(), 40960u);
    for (size_t i = 0; i < out.size(); ++i)
        EXPECT_DOUBLE_EQ(out[i], 0.0);
}

TEST(OneOverFNoise, GeometryPlacesTheMatrix)
{
    nisp::Geometry g;   // 2048 × 2048, reference border 4, matrix at (0, 0)
    EXPECT_EQ(g.place(2048, 2048), std::make_pair(size_t(0), size_t(0)));
    EXPECT_EQ(g.place(2040, 2040), std::make_pair(size_t(0), size_t(0)));
    EXPECT_EQ(g.place(2044, 256),  std::make_pair(size_t(0), size_t(0)));
    EXPECT_THROW(g.place(2050, 10), std::invalid_argument);
    EXPECT_THROW(g.place(10, 2049), std::invalid_argument);

    EXPECT_TRUE (g.isReference(0, 1000));
    EXPECT_TRUE (g.isReference(3, 1000));
    EXPECT_FALSE(g.isReference(4, 1000));
    EXPECT_FALSE(g.isReference(2043, 1000));
    EXPECT_TRUE (g.isReference(2044, 1000));
    EXPECT_TRUE (g.isReference(1000, 2047));

    PixelConfig cfg;
    cfg.set("detectorSize",   nlohmann::json::array({128, 64}));
    cfg.set("detectorOrigin", nlohmann::json::array({64, 0}));
    cfg.set("refBorder",      2);
    g.configure(cfg);
    EXPECT_EQ(g.width,  128u);
    EXPECT_EQ(g.height, 64u);
    EXPECT_EQ(g.refBorder, 2u);
    EXPECT_EQ(g.place(64, 64), std::make_pair(size_t(64), size_t(0)));
    EXPECT_THROW(g.place(65, 64), std::invalid_argument);   // does not fit at that origin

    PixelConfig badBorder;
    badBorder.set("refBorder", 100);
    EXPECT_THROW(g.configure(badBorder), std::invalid_argument);
    EXPECT_EQ(g.refBorder, 2u);   // unchanged after the error
}

TEST(OneOverFNoise, PartialChannelsLineUpWithThePhysicalChannels)
{
    // A centred 100-wide matrix starts at physical column 974: its columns 0–49
    // are the end of channel 15 (pixel-clock offsets 14–63) and its columns 50–99
    // the start of channel 16 (offsets 0–49). With bias voltage noise only, the
    // same offset carries the same value: column k of channel 15 (offset 14 + k)
    // matches column 64 + k of channel 16.
    OneOverFNoise n;
    PixelConfig cfg;
    cfg.set("bias_voltage_noise", 10.0);
    n.initialize(cfg);

    const size_t W = 100, H = 8;
    placeAtCentre(n, W, H);
    const auto out = n.generate(W, H);
    for (size_t y = 0; y < H; ++y)
        for (size_t k = 0; k < 36; ++k)
            ASSERT_DOUBLE_EQ(out[y * W + k], out[y * W + 64 + k]) << "row " << y << ", column " << k;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Noise statistics
// ═══════════════════════════════════════════════════════════════════════════════

// A 128-wide matrix at the default origin (0, 0) covers exactly the physical
// channels 0 and 1, so its columns 0–63 are one whole channel.

TEST(OneOverFNoise, GeneratedNoiseRMSWithinBroadBandOfAmplitude)
{
    OneOverFNoise n(10.0, 1.0, 0, 64, 16);
    const auto out = n.generate(128, 800);   // channel samples = 800*64 = 51200
    const auto series = extractChannelSeries(out, 128, 800, 0, 64);

    expectRmsWithinBand(rms(series), 10.0);
}

TEST(OneOverFNoise, PositiveLag1AutocorrelationConfirmsColoredNoise)
{
    OneOverFNoise n(10.0, 1.0, 0, 64, 16);
    const auto out = n.generate(128, 800);
    const auto series = extractChannelSeries(out, 128, 800, 0, 64);

    // Loose threshold, well above the white-noise null (~1/sqrt(N)).
    EXPECT_GT(lag1Autocorr(series), 0.02);
}

TEST(OneOverFNoise, IndependentChannelsUncorrelated)
{
    OneOverFNoise n(10.0, 1.0, 0, 64, 16); // width=128 -> 2 whole channels
    const auto out = n.generate(128, 800);

    const auto ch0 = extractChannelSeries(out, 128, 800, 0,  64);
    const auto ch1 = extractChannelSeries(out, 128, 800, 64, 64);

    // See firstDifference()'s comment: correlating the raw (strongly
    // autocorrelated) series directly is unreliable — this failed
    // intermittently with raw cross-correlations up to ~0.45 despite the
    // channels being genuinely independent (separate RNG per channel).
    const auto d0 = firstDifference(ch0);
    const auto d1 = firstDifference(ch1);

    EXPECT_LT(std::abs(crossCorr(d0, d1)), 0.1);
}

TEST(OneOverFNoise, ScaleWithAmplitudeIsMonotonic)
{
    OneOverFNoise lo(5.0,  1.0, 0, 64, 16);
    OneOverFNoise hi(20.0, 1.0, 0, 64, 16);

    const auto outLo = lo.generate(128, 500);
    const auto outHi = hi.generate(128, 500);

    const double rmsLo = rms(extractChannelSeries(outLo, 128, 500, 0, 64));
    const double rmsHi = rms(extractChannelSeries(outHi, 128, 500, 0, 64));

    expectRmsWithinBand(rmsLo, 5.0);
    expectRmsWithinBand(rmsHi, 20.0);
    EXPECT_GT(rmsHi, rmsLo);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Windowed generate()
// ═══════════════════════════════════════════════════════════════════════════════

TEST(OneOverFNoise, WindowOutOfRangeThrows)
{
    OneOverFNoise n(5.0, 1.0, 0, 64, 16);
    EXPECT_THROW(n.generate(64, 64, 60, 0, 10, 10), std::out_of_range);
}

TEST(OneOverFNoise, DisabledWindowSkipsBoundsCheck)
{
    OneOverFNoise n; // disabled
    std::valarray<double> out;
    EXPECT_NO_THROW(out = n.generate(64, 64, 60, 0, 10, 10));
    ASSERT_EQ(out.size(), 100u);
    for (size_t i = 0; i < out.size(); ++i)
        EXPECT_DOUBLE_EQ(out[i], 0.0);
}

TEST(OneOverFNoise, WindowedOutputWithinBroadBandOfAmplitude)
{
    OneOverFNoise n(10.0, 1.0, 0, 64, 16);
    // Matrix is 128x500; the window is its first whole channel (c0=r0=0, ncols=64, nrows=500).
    const auto out = n.generate(128, 500, 0, 0, 64, 500);
    const auto series = extractChannelSeries(out, 64, 500, 0, 64);

    expectRmsWithinBand(rms(series), 10.0);
}

TEST(OneOverFNoise, WindowNotAlignedToChannelZeroStillWorks)
{
    OneOverFNoise n(10.0, 1.0, 0, 64, 16); // width=128 -> 2 whole channels
    // Matrix is 128x500; window selects its second channel only (cols [64,128)).
    const auto out = n.generate(128, 500, /*c0=*/64, /*r0=*/0, /*ncols=*/64, /*nrows=*/500);
    const auto series = extractChannelSeries(out, 64, 500, 0, 64);

    expectRmsWithinBand(rms(series), 10.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Frame time
// ═══════════════════════════════════════════════════════════════════════════════

TEST(OneOverFNoise, SetFrameTimeRejectsNonPositive)
{
    OneOverFNoise n;
    EXPECT_DOUBLE_EQ(n.getFrameTime(), 1.45408);   // NISP default
    EXPECT_THROW(n.setFrameTime(0.0),  std::invalid_argument);
    EXPECT_THROW(n.setFrameTime(-1.0), std::invalid_argument);
    n.setFrameTime(2.0);
    EXPECT_DOUBLE_EQ(n.getFrameTime(), 2.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  bias_voltage_noise: 1/f noise common to every channel
// ═══════════════════════════════════════════════════════════════════════════════

TEST(OneOverFNoise, InitializeBiasVoltageNoiseScalarAndObject)
{
    OneOverFNoise n;
    PixelConfig scalar;
    scalar.set("bias_voltage_noise", 3.0);
    n.initialize(scalar);

    EXPECT_TRUE(n.isEnabled());                 // enabled on its own, without "oneOverF"
    EXPECT_DOUBLE_EQ(n.getAmplitude(),        0.0);
    EXPECT_DOUBLE_EQ(n.getBiasVoltageNoise(), 3.0);
    EXPECT_DOUBLE_EQ(n.getBiasSlope(),        1.0);
    EXPECT_EQ(n.getBiasNPoles(),              16u);

    nlohmann::json obj;
    obj["amplitude"] = 2.0;
    obj["slope"]     = 1.5;
    obj["nPoles"]    = 8;
    PixelConfig object;
    object.set("bias_voltage_noise", obj);
    n.initialize(object);

    EXPECT_DOUBLE_EQ(n.getBiasVoltageNoise(), 2.0);
    EXPECT_DOUBLE_EQ(n.getBiasSlope(),        1.5);
    EXPECT_EQ(n.getBiasNPoles(),              8u);
}

TEST(OneOverFNoise, InvalidBiasVoltageNoiseThrows)
{
    auto initWith = [](const PixelConfig& cfg) { OneOverFNoise n; n.initialize(cfg); };

    PixelConfig negative;
    negative.set("bias_voltage_noise", -1.0);
    EXPECT_THROW(initWith(negative), std::invalid_argument);

    nlohmann::json typo;
    typo["amplitud"] = 2.0;
    PixelConfig unknownKey;
    unknownKey.set("bias_voltage_noise", typo);
    EXPECT_THROW(initWith(unknownKey), std::invalid_argument);

    nlohmann::json zeroPoles;
    zeroPoles["amplitude"] = 2.0;
    zeroPoles["nPoles"]    = 0;
    PixelConfig badPoles;
    badPoles.set("bias_voltage_noise", zeroPoles);
    EXPECT_THROW(initWith(badPoles), std::invalid_argument);
}

TEST(OneOverFNoise, BiasVoltageNoiseIsCommonToAllChannels)
{
    // Bias voltage noise only: every channel carries the same series, so sample s
    // of a row is identical in every channel (uniform layout, time offset 0).
    OneOverFNoise n;
    PixelConfig cfg;
    cfg.set("bias_voltage_noise", 10.0);
    n.initialize(cfg);

    const size_t W = 256, H = 800;   // 4 channels of 64 columns
    const auto out = n.generate(W, H);

    const auto ch0 = extractChannelSeries(out, W, H, 0, 64);
    for (size_t c = 1; c < 4; ++c)
    {
        const auto chc = extractChannelSeries(out, W, H, c * 64, 64);
        for (size_t i = 0; i < ch0.size(); ++i)
            ASSERT_DOUBLE_EQ(chc[i], ch0[i]) << "channel " << c << ", sample " << i;
    }

    expectRmsWithinBand(rms(ch0), 10.0);
}

TEST(OneOverFNoise, BiasVoltageNoiseFollowsPixelClockInReferenceLayout)
{
    // 2040-wide frame at x0 = 4 (the active area): channel 0 lost its 4 reference
    // columns, which were read first, so its column s is pixel-clock sample s + 4.
    // Channel 1 starts at column 60 and channel 31 at column 1980, both at sample 0.
    OneOverFNoise n;
    PixelConfig cfg;
    cfg.set("bias_voltage_noise", 10.0);
    n.initialize(cfg);

    const size_t W = 2040, H = 4;
    placeAtCentre(n, W, H);
    const auto out = n.generate(W, H);

    for (size_t y = 0; y < H; ++y)
        for (size_t s = 0; s < 60; ++s)
        {
            ASSERT_DOUBLE_EQ(out[y * W + s],        out[y * W + 60 + s + 4]) << "row " << y << ", column " << s;
            ASSERT_DOUBLE_EQ(out[y * W + 1980 + s], out[y * W + 60 + s])     << "row " << y << ", column " << s;
        }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  pedestal_drift: channel offsets, mean f(t) + noise per channel and per frame
// ═══════════════════════════════════════════════════════════════════════════════

TEST(OneOverFNoise, InitializePedestalDriftScalarIsWhiteNoiseOnly)
{
    OneOverFNoise n;
    PixelConfig cfg;
    cfg.set("pedestal_drift", 4.0);
    n.initialize(cfg);

    EXPECT_TRUE(n.isEnabled());
    EXPECT_DOUBLE_EQ(n.getPedestalSigma(), 4.0);
    EXPECT_FALSE(n.hasPedestalDrift());
    EXPECT_DOUBLE_EQ(n.pedestalMean(10.0), 0.0);
}

TEST(OneOverFNoise, PedestalDriftEquationSumsEveryTerm)
{
    nlohmann::json expTerm;
    expTerm["amplitude"] = 20.0;
    expTerm["tau"]       = 30.0;

    nlohmann::json sinTerm;               // a single object is accepted as well as an array
    sinTerm["amplitude"] = 2.0;
    sinTerm["period"]    = 300.0;
    sinTerm["phase"]     = 0.25;

    nlohmann::json ped;
    ped["sigma"]         = 1.5;
    ped["t0"]            = 100.0;
    ped["drift"]["poly"] = nlohmann::json::array({1.0, 0.5, 0.01});   // 1 + 0.5 t + 0.01 t²
    ped["drift"]["exp"]  = nlohmann::json::array({expTerm});
    ped["drift"]["sin"]  = sinTerm;

    OneOverFNoise n;
    PixelConfig cfg;
    cfg.set("pedestal_drift", ped);
    n.initialize(cfg);

    EXPECT_DOUBLE_EQ(n.getPedestalSigma(), 1.5);
    EXPECT_DOUBLE_EQ(n.getPedestalT0(),    100.0);
    EXPECT_TRUE(n.hasPedestalDrift());

    const double t  = 7.0;
    const double tt = 100.0 + t;          // t0 + time since reset
    const double expected = 1.0 + 0.5 * tt + 0.01 * tt * tt
                          + 20.0 * std::exp(-tt / 30.0)
                          + 2.0 * std::sin(2.0 * std::numbers::pi * tt / 300.0 + 0.25);
    EXPECT_NEAR(n.pedestalMean(t), expected, 1e-9);
}

TEST(OneOverFNoise, InvalidPedestalDriftThrows)
{
    auto initWith = [](const nlohmann::json& ped)
    {
        OneOverFNoise n;
        PixelConfig cfg;
        cfg.set("pedestal_drift", ped);
        n.initialize(cfg);
    };

    nlohmann::json negativeSigma;
    negativeSigma["sigma"] = -1.0;
    EXPECT_THROW(initWith(negativeSigma), std::invalid_argument);

    nlohmann::json typo;
    typo["sigma"] = 1.0;
    typo["drfit"] = nlohmann::json::object();
    EXPECT_THROW(initWith(typo), std::invalid_argument);

    nlohmann::json polyNotArray;
    polyNotArray["drift"]["poly"] = 3.0;
    EXPECT_THROW(initWith(polyNotArray), std::invalid_argument);

    nlohmann::json zeroTau;
    zeroTau["amplitude"] = 1.0;
    zeroTau["tau"]       = 0.0;
    nlohmann::json badTau;
    badTau["drift"]["exp"] = zeroTau;
    EXPECT_THROW(initWith(badTau), std::invalid_argument);

    nlohmann::json sinNoPeriod;
    sinNoPeriod["amplitude"] = 1.0;
    nlohmann::json noPeriod;
    noPeriod["drift"]["sin"] = sinNoPeriod;
    EXPECT_THROW(initWith(noPeriod), std::invalid_argument);
}

TEST(OneOverFNoise, PedestalNoiseIsOneOffsetPerChannelPerFrame)
{
    // Pedestal noise only: within a frame every pixel of a channel carries the same
    // offset; offsets are independent between channels and frames, with RMS sigma.
    OneOverFNoise n;
    PixelConfig cfg;
    cfg.set("pedestal_drift", 4.0);
    n.initialize(cfg);

    const size_t W = 2048, H = 2, nCh = 32, nFrames = 100;
    std::vector<double> offsets, ch0, ch1;
    for (size_t frame = 0; frame < nFrames; ++frame)
    {
        const auto out = n.generate(W, H);
        for (size_t ch = 0; ch < nCh; ++ch)
        {
            const double v = out[ch * 64];
            for (size_t y = 0; y < H; ++y)
                for (size_t s = 0; s < 64; ++s)
                    ASSERT_DOUBLE_EQ(out[y * W + ch * 64 + s], v) << "frame " << frame << ", channel " << ch;
            offsets.push_back(v);
        }
        ch0.push_back(out[0]);
        ch1.push_back(out[64]);
    }

    // 3200 draws: the RMS is known to ~1.3 % and the mean to σ/√3200 ≈ 0.07,
    // so these bands are > 5σ margins.
    EXPECT_NEAR(rms(offsets),  4.0, 0.4);
    EXPECT_NEAR(mean(offsets), 0.0, 0.4);
    // 100 frames: σ(r) ≈ 0.1 for independent channels.
    EXPECT_LT(std::abs(crossCorr(ch0, ch1)), 0.4);
}

TEST(OneOverFNoise, PedestalDriftEvaluatedAtEachRowReadTime)
{
    // f(t) = 2 t, t0 = 1 s, frame time 1 s. A centred 4-row matrix covers the
    // physical rows 1022–1025 of the 2048-row detector: its row y of the frame
    // read from tStart = 3 s is read at 3 + (1022 + y + 0.5)/2048 s.
    nlohmann::json ped;
    ped["t0"]            = 1.0;
    ped["drift"]["poly"] = nlohmann::json::array({0.0, 2.0});

    OneOverFNoise n;
    n.setFrameTime(1.0);
    PixelConfig cfg;
    cfg.set("pedestal_drift", ped);
    n.initialize(cfg);

    const size_t W = 128, H = 4;
    placeAtCentre(n, W, H);
    const auto out = n.generate(W, H, 3.0);
    for (size_t y = 0; y < H; ++y)
    {
        const double expected = 2.0 * (1.0 + 3.0 + (1022.0 + static_cast<double>(y) + 0.5) / 2048.0);
        for (size_t x = 0; x < W; ++x)
            ASSERT_NEAR(out[y * W + x], expected, 1e-12) << "row " << y << ", column " << x;
    }
}

TEST(OneOverFNoise, WindowedPedestalDriftUsesDetectorRows)
{
    // f(t) = t, frame time 1 s. A centred 8-row matrix covers the physical rows
    // 1020–1027; a window on its rows 2–3 is read at (1022.5/2048) s and (1023.5/2048) s.
    nlohmann::json ped;
    ped["drift"]["poly"] = nlohmann::json::array({0.0, 1.0});

    OneOverFNoise n;
    n.setFrameTime(1.0);
    PixelConfig cfg;
    cfg.set("pedestal_drift", ped);
    n.initialize(cfg);

    placeAtCentre(n, 128, 8);
    const auto out = n.generate(128, 8, /*c0=*/64, /*r0=*/2, /*ncols=*/32, /*nrows=*/2, /*tStart=*/0.0);
    for (size_t r = 0; r < 2; ++r)
        for (size_t c = 0; c < 32; ++c)
            ASSERT_NEAR(out[r * 32 + c], (1022.0 + static_cast<double>(r) + 0.5) / 2048.0, 1e-12);
}

TEST(OneOverFNoise, StreamCountsFramesForPedestalDrift)
{
    // f(t) = t, frame time 2 s. A centred 1-row matrix is physical row 1023:
    // in frame k of the stream it is read at (k + 1023.5/2048)·2 s.
    nlohmann::json ped;
    ped["drift"]["poly"] = nlohmann::json::array({0.0, 1.0});

    OneOverFNoise n;
    n.setFrameTime(2.0);
    PixelConfig cfg;
    cfg.set("pedestal_drift", ped);
    n.initialize(cfg);

    placeAtCentre(n, 64, 1);
    auto stream = n.beginStream(64, 1);
    for (size_t k = 0; k < 5; ++k)
    {
        const auto out = stream.next();
        for (size_t x = 0; x < 64; ++x)
            ASSERT_NEAR(out[x], (static_cast<double>(k) + 1023.5 / 2048.0) * 2.0, 1e-12) << "frame " << k;
    }
}

TEST(OneOverFNoise, SlopeIsTheSpectralIndex)
{
    // Lag-1 autocorrelation of the stationary series: Σ Var_k ρ_k / Σ Var_k.
    //   slope 1 (pink): the same variance for every pole → ≈ 0.83 with 16 poles over
    //                   a 2048 × 64 frame (it was ≈ 0.99 when slope 1 gave 1/f²);
    //   slope 0 (variance ∝ f_k, nearly white) → ≈ 0.2.
    // One full channel (W = 128 is two whole channels at x0 = 0), 51 200 samples.
    auto lag1 = [](double slope)
    {
        OneOverFNoise n(10.0, slope, 0, 64, 16);
        const auto out = n.generate(128, 800);
        return lag1Autocorr(extractChannelSeries(out, 128, 800, 0, 64));
    };

    const double pink  = lag1(1.0);
    const double white = lag1(0.0);
    EXPECT_GT(pink,  0.6);
    EXPECT_LT(pink,  0.93);
    EXPECT_LT(white, 0.4);
}

TEST(OneOverFNoise, SmallMatrixStreamRunsOverPartialChannels)
{
    // Every component on, on a centred 100 × 16 matrix (two partial channels, rows
    // skipped between frames): the stream must produce full-size frames with noise in them.
    nlohmann::json bias;
    bias["amplitude"] = 5.0;
    PixelConfig cfg;
    cfg.set("oneOverF", 10.0);
    cfg.set("bias_voltage_noise", bias);
    cfg.set("pedestal_drift", 2.0);

    OneOverFNoise n;
    n.initialize(cfg);

    placeAtCentre(n, 100, 16);
    auto stream = n.beginStream(100, 16);
    for (size_t k = 0; k < 3; ++k)
    {
        const auto out = stream.next();
        ASSERT_EQ(out.size(), 100u * 16u);
        EXPECT_GT(rms(std::vector<double>(std::begin(out), std::end(out))), 1.0) << "frame " << k;
    }
}
