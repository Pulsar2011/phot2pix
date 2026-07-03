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
#include <valarray>
#include <vector>

using namespace nisp;
using namespace phot2pix;

// ── Helpers ──────────────────────────────────────────────────────────────────
//
// NOTE on tolerances used throughout this file: OneOverFNoise::generate()
// builds each channel as a sum of AR(1) poles with corner frequencies
// log-spaced down to f_min = 1/(2*channelLength) — i.e. the slowest pole's
// correlation time is comparable to the *entire* observed series. This is
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

TEST(OneOverFNoise, ExplicitChannelCountWithRemainder)
{
    // nCh=5, colPerChannel=64: nominal_total=320, matches neither width(101)
    // nor width+8(109) -> fallback layout: col_w=20, rem=1 (channel 0 gets 21).
    OneOverFNoise n(5.0, 1.0, /*nChannels=*/5, /*colPerChannel=*/64, 16);
    const auto out = n.generate(101, 8);
    EXPECT_EQ(out.size(), 101u * 8u);

    const auto ch0 = extractChannelSeries(out, 101, 8, 0,  21);
    const auto ch1 = extractChannelSeries(out, 101, 8, 21, 20);

    // Independent per-channel RNG draws -> not element-wise identical over
    // their common length.
    const size_t n0 = std::min(ch0.size(), ch1.size());
    bool anyDifferent = false;
    for (size_t i = 0; i < n0; ++i)
        if (ch0[i] != ch1[i]) { anyDifferent = true; break; }
    EXPECT_TRUE(anyDifferent);
}

TEST(OneOverFNoise, ZeroDerivedChannelsThrows)
{
    // nCh = round(1/64) = 0.
    OneOverFNoise n(5.0, 1.0, 0, 64, 16);
    EXPECT_THROW(n.generate(1, 10), std::invalid_argument);
}

TEST(OneOverFNoise, DisabledGenerateNeverBuildsLayoutOrThrows)
{
    // Disabled (amplitude==0): the isEnabled() early-return happens before
    // _buildLayout() is ever called, so a width that WOULD throw if enabled
    // (see ZeroDerivedChannelsThrows) must not throw here.
    OneOverFNoise n; // disabled
    std::valarray<double> out;
    EXPECT_NO_THROW(out = n.generate(1, 10));
    ASSERT_EQ(out.size(), 10u);
    for (size_t i = 0; i < out.size(); ++i)
        EXPECT_DOUBLE_EQ(out[i], 0.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Noise statistics
// ═══════════════════════════════════════════════════════════════════════════════

TEST(OneOverFNoise, GeneratedNoiseRMSWithinBroadBandOfAmplitude)
{
    OneOverFNoise n(10.0, 1.0, 0, 64, 16); // single channel (width==colPerChannel)
    const auto out = n.generate(64, 800);  // channel length = 800*64 = 51200
    const auto series = extractChannelSeries(out, 64, 800, 0, 64);

    expectRmsWithinBand(rms(series), 10.0);
}

TEST(OneOverFNoise, PositiveLag1AutocorrelationConfirmsColoredNoise)
{
    OneOverFNoise n(10.0, 1.0, 0, 64, 16);
    const auto out = n.generate(64, 800);
    const auto series = extractChannelSeries(out, 64, 800, 0, 64);

    // Loose threshold, well above the white-noise null (~1/sqrt(N)).
    EXPECT_GT(lag1Autocorr(series), 0.02);
}

TEST(OneOverFNoise, IndependentChannelsUncorrelated)
{
    OneOverFNoise n(10.0, 1.0, 0, 64, 16); // width=128 -> 2 channels
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

    const auto outLo = lo.generate(64, 500);
    const auto outHi = hi.generate(64, 500);

    const double rmsLo = rms(extractChannelSeries(outLo, 64, 500, 0, 64));
    const double rmsHi = rms(extractChannelSeries(outHi, 64, 500, 0, 64));

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
    // Detector is 64x500; window covers it entirely (c0=r0=0, ncols=64, nrows=500).
    const auto out = n.generate(64, 500, 0, 0, 64, 500);
    const auto series = extractChannelSeries(out, 64, 500, 0, 64);

    expectRmsWithinBand(rms(series), 10.0);
}

TEST(OneOverFNoise, WindowNotAlignedToChannelZeroStillWorks)
{
    OneOverFNoise n(10.0, 1.0, 0, 64, 16); // width=128 -> 2 channels
    // Detector is 128x500; window selects channel 1 only (cols [64,128)).
    const auto out = n.generate(128, 500, /*c0=*/64, /*r0=*/0, /*ncols=*/64, /*nrows=*/500);
    const auto series = extractChannelSeries(out, 64, 500, 0, 64);

    expectRmsWithinBand(rms(series), 10.0);
}
