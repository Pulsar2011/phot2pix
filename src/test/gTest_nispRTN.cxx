//
// Created by gillard on 2026-11-27.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_nispRTN.cxx
//
// Random telegraph noise of the NISP pixels (nisp::RTN): configuration, draw of the RTN pixels,
// time evolution of their levels, and its use by nisp::SCA at each read.
//
// The statistical tolerances are about 6 sigma of a binomial count: they catch a wrong
// probability, not a fluctuation.

#include <gtest/gtest.h>
#include <nisp/nispRTN.h>
#include <nisp/nispSCA.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <valarray>
#include <vector>

using namespace nisp;
using namespace phot2pix;

namespace
{
    // "rtn" object with the four required keys.
    nlohmann::json rtnObject(double fraction, double tauMin, double tauMax, double amplitude)
    {
        return nlohmann::json{{"fraction", fraction}, {"tauMin", tauMin}, {"tauMax", tauMax}, {"amplitude", amplitude}};
    }

    // Configuration holding only "rtn".
    PixelConfig rtnConfig(const nlohmann::json& rtn)
    {
        return PixelConfig{{"rtn", rtn}};
    }

    // RTN of a width × height matrix.
    RTN makeRTN(const nlohmann::json& rtn, size_t width, size_t height)
    {
        RTN model;
        model.initialize(rtnConfig(rtn), width, height);
        return model;
    }

    // Fraction of the RTN pixels in the shifted level.
    double shiftedFraction(const RTN::Stream& stream)
    {
        size_t n = 0;
        for (uint8_t s : stream.states())
            n += s;
        return static_cast<double>(n) / static_cast<double>(stream.states().size());
    }

    // Offset (ADU) of RTN pixel k of the stream: s·A when shifted, 0 at rest.
    double offset(const RTN& model, const RTN::Stream& stream, size_t k)
    {
        return stream.states()[k] ? model.getPixels()[k].sign * model.getPixels()[k].amplitude : 0.0;
    }
}

// ── Configuration ────────────────────────────────────────────────────────────

TEST(RTN, OffByDefault)
{
    const RTN model;
    EXPECT_FALSE(model.isEnabled());
    EXPECT_TRUE(model.getPixels().empty());
    EXPECT_DOUBLE_EQ(model.getPositiveFraction(), 0.4);
    EXPECT_FALSE(model.getSeed().has_value());
}

TEST(RTN, NothingChangesWithoutTheKey)
{
    RTN model = makeRTN(rtnObject(0.5, 1.0, 2.0, 10.0), 10, 10);
    ASSERT_TRUE(model.isEnabled());

    model.initialize(PixelConfig{{"gain", 2.0}}, 10, 10);
    EXPECT_TRUE(model.isEnabled());
    EXPECT_EQ(model.getPixels().size(), 50u);
}

TEST(RTN, FractionZeroOrAmplitudeZeroIsOff)
{
    EXPECT_FALSE(makeRTN(rtnObject(0.0, 1.0, 2.0, 10.0), 10, 10).isEnabled());
    EXPECT_THROW(makeRTN(rtnObject(0.5, 1.0, 2.0, 0.0),  10, 10).isEnabled(),std::invalid_argument);
    EXPECT_THROW(makeRTN(rtnObject(0.5, 1.0, 2.0, -5.0),  10, 10).isEnabled(),std::invalid_argument);
}

TEST(RTN, DrawsDistinctPixelsWithinTheRanges)
{
    const size_t W = 100, H = 50;
    const RTN model = makeRTN(rtnObject(0.1, 0.5, 20.0, 30.0), W, H);

    const auto& pixels = model.getPixels();
    ASSERT_EQ(pixels.size(), 500u);   // round(0.1 · 5000)
    EXPECT_EQ(model.getWidth(),  W);
    EXPECT_EQ(model.getHeight(), H);

    for (size_t k = 0; k < pixels.size(); ++k)
    {
        EXPECT_LT(pixels[k].idx, W * H);
        if (k > 0)
            EXPECT_GT(pixels[k].idx, pixels[k - 1].idx);   // distinct, in increasing order
        EXPECT_GE(pixels[k].tau0, 0.5);
        EXPECT_LE(pixels[k].tau0, 20.0);
        EXPECT_GE(pixels[k].tau1, 0.5);
        EXPECT_LE(pixels[k].tau1, 20.0);
        EXPECT_GE(pixels[k].amplitude, 1.0);   // uniform in [1 ADU, amplitude]
        EXPECT_LE(pixels[k].amplitude, 30.0);
        EXPECT_TRUE(pixels[k].sign == -1 || pixels[k].sign == +1);
    }
}

TEST(RTN, EachPixelHasItsOwnAmplitude)
{
    const RTN model = makeRTN(rtnObject(1.0, 1.0, 1.0, 31.0), 100, 100);

    double sum = 0.0;
    bool   differs = false;
    for (const auto& px : model.getPixels())
    {
        sum     += px.amplitude;
        differs |= (px.amplitude != model.getPixels().front().amplitude);
    }
    EXPECT_TRUE(differs);
    EXPECT_NEAR(sum / 10000.0, 16.0, 0.6);   // mean of U[1, 31]: 16, sigma of the mean 0.087

    // amplitude 1: every pixel has the smallest amplitude.
    for (const auto& px : makeRTN(rtnObject(0.5, 1.0, 1.0, 1.0), 10, 10).getPixels())
        EXPECT_DOUBLE_EQ(px.amplitude, 1.0);
}

TEST(RTN, RoundsTheNumberOfPixels)
{
    EXPECT_EQ(makeRTN(rtnObject(0.25, 1.0, 1.0, 1.0), 10, 1).getPixels().size(), 3u);   // 2.5 → 3

    const RTN all = makeRTN(rtnObject(1.0, 1.0, 1.0, 1.0), 4, 3);
    ASSERT_EQ(all.getPixels().size(), 12u);
    for (size_t k = 0; k < 12; ++k)
        EXPECT_EQ(all.getPixels()[k].idx, k);
}

TEST(RTN, EqualTauMinAndTauMaxFixTheTimeConstants)
{
    const RTN model = makeRTN(rtnObject(0.5, 3.0, 3.0, 1.0), 10, 10);
    for (const auto& px : model.getPixels())
    {
        EXPECT_DOUBLE_EQ(px.tau0, 3.0);
        EXPECT_DOUBLE_EQ(px.tau1, 3.0);
        EXPECT_DOUBLE_EQ(px.shiftedProbability(), 0.5);
    }
}

TEST(RTN, PositiveFractionSetsTheSigns)
{
    auto signs = [](double positiveFraction) {
        nlohmann::json j = rtnObject(1.0, 1.0, 1.0, 1.0);
        if (positiveFraction >= 0.0)
            j["positiveFraction"] = positiveFraction;
        const RTN model = makeRTN(j, 100, 100);

        size_t positive = 0;
        for (const auto& px : model.getPixels())
            positive += (px.sign > 0);
        return static_cast<double>(positive) / static_cast<double>(model.getPixels().size());
    };

    EXPECT_DOUBLE_EQ(signs(0.0), 0.0);
    EXPECT_DOUBLE_EQ(signs(1.0), 1.0);
    EXPECT_NEAR(signs(-1.0), 0.4, 0.03);   // absent: 40 % up, 60 % down
}

TEST(RTN, SameSeedSameDraw)
{
    nlohmann::json j = rtnObject(0.01, 0.5, 20.0, 30.0);
    j["seed"] = 42;
    const RTN a = makeRTN(j, 100, 100);
    const RTN b = makeRTN(j, 100, 100);
    j["seed"] = 43;
    const RTN c = makeRTN(j, 100, 100);

    ASSERT_EQ(a.getSeed().value_or(0), 42u);
    ASSERT_EQ(a.getPixels().size(), 100u);
    ASSERT_EQ(b.getPixels().size(), 100u);
    ASSERT_EQ(c.getPixels().size(), 100u);

    bool differs = false;
    for (size_t k = 0; k < a.getPixels().size(); ++k)
    {
        EXPECT_EQ(a.getPixels()[k].idx,  b.getPixels()[k].idx);
        EXPECT_EQ(a.getPixels()[k].tau0, b.getPixels()[k].tau0);
        EXPECT_EQ(a.getPixels()[k].tau1, b.getPixels()[k].tau1);
        EXPECT_EQ(a.getPixels()[k].amplitude, b.getPixels()[k].amplitude);
        EXPECT_EQ(a.getPixels()[k].sign, b.getPixels()[k].sign);
        differs |= (a.getPixels()[k].idx != c.getPixels()[k].idx);
    }
    EXPECT_TRUE(differs);
}

TEST(RTN, RejectsInvalidConfigurations)
{
    auto with = [](const char* key, const nlohmann::json& value) {
        nlohmann::json j = rtnObject(0.1, 1.0, 2.0, 10.0);
        j[key] = value;
        return j;
    };
    auto without = [](const char* key) {
        nlohmann::json j = rtnObject(0.1, 1.0, 2.0, 10.0);
        j.erase(key);
        return j;
    };

    const std::vector<nlohmann::json> invalid = {
        without("fraction"), without("tauMin"), without("tauMax"), without("amplitude"),
        with("tau", 1.0),                                              // unknown key
        with("fraction", -0.1), with("fraction", 1.5), with("fraction", "a"),
        with("tauMin", 0.0), with("tauMin", -1.0), with("tauMax", 0.5),   // tauMax < tauMin
        with("amplitude", -1.0), with("amplitude", 0.5),                 // amplitude < 1 ADU
        with("positiveFraction", -0.1), with("positiveFraction", 1.2),
        with("seed", -1), with("seed", 1.5), with("seed", "1"),
        nlohmann::json::array({0.1, 1.0, 2.0, 10.0}),                     // not an object
    };

    RTN model = makeRTN(rtnObject(0.5, 1.0, 2.0, 10.0), 10, 10);
    for (const auto& j : invalid)
    {
        EXPECT_THROW(model.initialize(rtnConfig(j), 10, 10), std::invalid_argument) << j.dump();
        EXPECT_EQ(model.getPixels().size(), 50u) << j.dump();   // unchanged
        EXPECT_DOUBLE_EQ(model.getFraction(), 0.5) << j.dump();
    }

    const PixelConfig number{{"rtn", 1.0}};
    EXPECT_THROW(model.initialize(number, 10, 10), std::invalid_argument);   // a number
}

// ── Stream ───────────────────────────────────────────────────────────────────

TEST(RTNStream, StartsInTheStationaryState)
{
    const RTN model = makeRTN(rtnObject(1.0, 1.0, 1.0, 1.0), 100, 100);   // p = 0.5
    const RTN::Stream stream = model.beginStream();
    ASSERT_EQ(stream.states().size(), 10000u);
    EXPECT_NEAR(shiftedFraction(stream), 0.5, 0.03);
}

TEST(RTNStream, AddToAndWriteToGiveTheOffsets)
{
    const size_t W = 8, H = 6;
    nlohmann::json j = rtnObject(1.0, 1.0, 1.0, 7.0);
    j["positiveFraction"] = 0.5;
    const RTN model = makeRTN(j, W, H);
    const RTN::Stream stream = model.beginStream();

    // Whole matrix: every pixel is an RTN pixel (idx = k).
    std::valarray<double> values(1.0, W * H);
    stream.addTo(values);
    std::valarray<double> map(99.0, W * H);
    stream.writeTo(map);
    for (size_t k = 0; k < W * H; ++k)
    {
        EXPECT_DOUBLE_EQ(values[k], 1.0 + offset(model, stream, k)) << k;
        EXPECT_DOUBLE_EQ(map[k], offset(model, stream, k)) << k;
    }

    // Window [2, 5) × [1, 5).
    const size_t c0 = 2, r0 = 1, ncols = 3, nrows = 4;
    std::valarray<double> window(0.0, ncols * nrows);
    stream.addTo(window, c0, r0, ncols, nrows);
    for (size_t r = 0; r < nrows; ++r)
        for (size_t c = 0; c < ncols; ++c)
            EXPECT_DOUBLE_EQ(window[r * ncols + c], offset(model, stream, (r0 + r) * W + c0 + c)) << r << "," << c;
}

TEST(RTNStream, WriteToLeavesTheOtherPixels)
{
    const size_t W = 8, H = 6;
    const RTN model = makeRTN(rtnObject(0.25, 1.0, 1.0, 7.0), W, H);
    const RTN::Stream stream = model.beginStream();

    std::valarray<double> map(99.0, W * H);
    stream.writeTo(map);

    std::vector<bool> isRTN(W * H, false);
    for (size_t k = 0; k < model.getPixels().size(); ++k)
    {
        isRTN[model.getPixels()[k].idx] = true;
        EXPECT_DOUBLE_EQ(map[model.getPixels()[k].idx], offset(model, stream, k));
    }
    for (size_t i = 0; i < W * H; ++i)
        if (!isRTN[i])
            EXPECT_DOUBLE_EQ(map[i], 99.0) << i;
}

TEST(RTNStream, RejectsInvalidArguments)
{
    const RTN model = makeRTN(rtnObject(1.0, 1.0, 1.0, 7.0), 8, 6);
    RTN::Stream stream = model.beginStream();

    std::valarray<double> small(0.0, 47);
    EXPECT_THROW(stream.addTo(small), std::invalid_argument);
    EXPECT_THROW(stream.writeTo(small), std::invalid_argument);

    std::valarray<double> window(0.0, 12);
    EXPECT_THROW(stream.addTo(window, 6, 0, 3, 4), std::out_of_range);   // columns 6..8 of 8
    EXPECT_THROW(stream.addTo(window, 0, 3, 3, 4), std::out_of_range);   // rows 3..6 of 6
    EXPECT_THROW(stream.addTo(window, 0, 0, 4, 4), std::invalid_argument);   // 16 values needed

    EXPECT_THROW(stream.advance(-1.0), std::invalid_argument);
    EXPECT_THROW(stream.advance(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
}

TEST(RTNStream, NoTimeKeepsTheLevels)
{
    const RTN model = makeRTN(rtnObject(1.0, 1.0, 1.0, 1.0), 100, 100);
    RTN::Stream stream = model.beginStream();
    const std::vector<uint8_t> before = stream.states();

    stream.advance(0.0);
    EXPECT_EQ(stream.states(), before);

    stream.advance(1e-12);   // flip probability 1e-12 per pixel
    EXPECT_EQ(stream.states(), before);
}

TEST(RTNStream, FlipRateFollowsTheTimeConstants)
{
    // τ_0 = τ_1 = 10 s: over dt = 0.1 s a pixel flips with the probability ½ (1 − exp(−0.02)),
    // whatever its level, i.e. 99.0 flips among 10 000 pixels (σ ≈ 10).
    const RTN model = makeRTN(rtnObject(1.0, 10.0, 10.0, 1.0), 100, 100);
    RTN::Stream stream = model.beginStream();
    const std::vector<uint8_t> before = stream.states();

    stream.advance(0.1);

    size_t flips = 0;
    for (size_t k = 0; k < before.size(); ++k)
        flips += (stream.states()[k] != before[k]);
    EXPECT_NEAR(static_cast<double>(flips), 10000.0 * 0.5 * (1.0 - std::exp(-0.02)), 50.0);
}

TEST(RTNStream, LongTimeForgetsTheLevels)
{
    const RTN model = makeRTN(rtnObject(1.0, 1.0, 1.0, 1.0), 100, 100);
    RTN::Stream stream = model.beginStream();
    const std::vector<uint8_t> before = stream.states();

    stream.advance(1e3);   // λ·dt = 2000

    size_t kept = 0;
    for (size_t k = 0; k < before.size(); ++k)
        kept += (stream.states()[k] == before[k]);
    EXPECT_NEAR(shiftedFraction(stream), 0.5, 0.03);
    EXPECT_NEAR(static_cast<double>(kept) / 10000.0, 0.5, 0.03);   // independent of the past level
}

// ── SCA ──────────────────────────────────────────────────────────────────────

namespace
{
    constexpr size_t kW = 16, kH = 16, kNPix = kW * kH;
    constexpr double kBaseline  = 1000.0;   // ADU
    constexpr double kAmplitude = 50.0;     // ADU

    // SCA with every other effect off, gain 2 (the RTN offset must not be divided by it), every pixel
    // an RTN pixel with time constants tau (s).
    std::unique_ptr<nisp::SCA> makeSCA(const std::string& mode, double tau)
    {
        nlohmann::json rtn = rtnObject(1.0, tau, tau, kAmplitude);
        rtn["seed"] = 3;

        auto sca = std::make_unique<nisp::SCA>(kW, kH);
        sca->initialize(PixelConfig{
            {"gain", 2.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 1.0},
            {"KTCnoise", 0.0}, {"baseline", kBaseline},
            {"offset", 0.0}, {"scale", 1.0}, {"frame_time", 1.0},
            {"exposureMode", mode}, {"rtn", rtn}
        });
        return sca;
    }

    // True if v is the read of RTN pixel k at rest or shifted (rounded to the uint16 output).
    bool isLevel(const nisp::SCA& sca, size_t k, double v)
    {
        const nisp::RTN::Pixel& px = sca.getRTN().getPixels()[k];
        const double shifted = std::floor(kBaseline + px.sign * px.amplitude + 0.5);
        return v == kBaseline || v == shifted;
    }
}

TEST(RTNSCA, FluxReadIsTheBaselineOrTheShiftedLevel)
{
    auto sca = makeSCA("FLUX", 1e9);
    ASSERT_TRUE(sca->getRTN().isEnabled());
    ASSERT_EQ(sca->getRTN().getPixels().size(), kNPix);   // idx = k

    const std::valarray<double> flux(0.0, kNPix);
    auto result = sca->expose(flux, 10.0);
    const std::valarray<uint16_t>* d = result->GetData<uint16_t>();
    ASSERT_NE(d, nullptr);
    ASSERT_EQ(d->size(), kNPix);

    size_t shifted = 0;
    for (size_t k = 0; k < kNPix; ++k)
    {
        EXPECT_TRUE(isLevel(*sca, k, (*d)[k])) << k << ": " << (*d)[k];
        shifted += ((*d)[k] != kBaseline);
    }
    EXPECT_GT(shifted, 0u);       // p = 0.5 over 256 pixels
    EXPECT_LT(shifted, kNPix);

    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("RTNAMP"),  kAmplitude, 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("RTNFRAC"), 1.0,        1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("RTNNPIX"), static_cast<double>(kNPix), 1e-9);
    EXPECT_NEAR(result->HDU().GetDoubleValueForKey("RTNSEED"), 3.0,        1e-9);
}

TEST(RTNSCA, SlowRTNKeepsItsLevelAlongTheRamp)
{
    auto sca = makeSCA("MACC", 1e9);   // flip probability ~1e-8 per pixel over the ramp

    const std::valarray<double> flux(0.0, kNPix);
    auto result = sca->expose(flux, nisp::macc{3, 2, 1});
    ASSERT_EQ(result->Size(3), 6u);
    const std::valarray<uint16_t>* d = result->GetData<uint16_t>();
    ASSERT_NE(d, nullptr);

    for (size_t k = 0; k < kNPix; ++k)
    {
        EXPECT_TRUE(isLevel(*sca, k, (*d)[k])) << k;
        for (size_t layer = 1; layer < 6; ++layer)
            EXPECT_EQ((*d)[layer * kNPix + k], (*d)[k]) << k << ", layer " << layer;
    }
}

TEST(RTNSCA, FastRTNChangesFromReadToRead)
{
    auto sca = makeSCA("MACC", 1e-6);   // independent levels at every read

    const std::valarray<double> flux(0.0, kNPix);
    auto result = sca->expose(flux, nisp::macc{3, 2, 1});
    const std::valarray<uint16_t>* d = result->GetData<uint16_t>();
    ASSERT_NE(d, nullptr);

    size_t changing = 0;
    for (size_t k = 0; k < kNPix; ++k)
    {
        bool changes = false;
        for (size_t layer = 0; layer < 6; ++layer)
        {
            EXPECT_TRUE(isLevel(*sca, k, (*d)[layer * kNPix + k])) << k << ", layer " << layer;
            changes |= ((*d)[layer * kNPix + k] != (*d)[k]);
        }
        changing += changes;
    }
    EXPECT_GT(changing, kNPix / 2);   // a pixel keeps its level over 6 reads with the probability 1/32
}

TEST(RTNSCA, CDSCancelsASlowRTN)
{
    // Read 1 and read N see the same level: the RTN cancels in the difference.
    auto sca = makeSCA("CDS", 1e9);

    const std::valarray<double> flux(0.0, kNPix);
    auto result = sca->expose(flux, 10.0);
    const std::valarray<uint16_t>* d = result->GetData<uint16_t>();
    ASSERT_NE(d, nullptr);

    for (size_t k = 0; k < kNPix; ++k)
        EXPECT_EQ((*d)[k], uint16_t{0}) << k;
}
