//
// Created by gillard on 2026-07-03.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_nispPixel.cxx
//

#include <gtest/gtest.h>
#include <nisp/nispPixel.h>
#include <DSTmath/DSTstatistic.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <valarray>
#include <vector>

using namespace nisp;
using namespace phot2pix;

// ── Helpers ──────────────────────────────────────────────────────────────────

// Deterministic, noise-free UTR ramp: frame[i] = s*i.
static std::valarray<double> makeLinearRamp(double s, size_t nFrames)
{
    std::valarray<double> ramp(nFrames);
    for (size_t i = 0; i < nFrames; ++i)
        ramp[i] = s * static_cast<double>(i);
    return ramp;
}

// Reproduces Ramp::getSlope()'s closed-form flux for a perfectly constant CDS
// step `c` with zero read noise (beta == 0) — used to check the Kubic et al.
// 2016 formula is applied consistently, not to independently re-derive it.
static double expectedKubicFlux(double c, double gain, size_t nf, size_t nd)
{
    const double alpha = (1.0 - static_cast<double>(nf) * static_cast<double>(nf))
                        / (3.0 * static_cast<double>(nf) * (static_cast<double>(nf) + static_cast<double>(nd)));
    const double fe    = gain;   // e-/ADU
    const double gamma = (2.0 * fe) / (alpha + 1.0);
    const double ffit  = c * c;
    return (1.0 / gamma) * (std::sqrt(1.0 + gamma * gamma * ffit) - 1.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  nisp::pixel — construction
// ═══════════════════════════════════════════════════════════════════════════════

TEST(nispPixel, DefaultConstruction)
{
    nisp::pixel p;
    EXPECT_DOUBLE_EQ(p.getDarkCurrent(),       0.0);
    EXPECT_DOUBLE_EQ(p.getReadNoise(),         0.0);
    EXPECT_DOUBLE_EQ(p.getGain(),              1.0);
    EXPECT_DOUBLE_EQ(p.getQuantumEfficiency(), 1.0);
    EXPECT_DOUBLE_EQ(p.getKTCnoise(),          0.0);
    EXPECT_DOUBLE_EQ(p.getBaseline(),          0.0);
}

TEST(nispPixel, ParameterizedConstruction)
{
    nisp::pixel p(1.0, 2.0, 3.0, 0.4, 5.0, 6.0);
    EXPECT_DOUBLE_EQ(p.getDarkCurrent(),       1.0);
    EXPECT_DOUBLE_EQ(p.getReadNoise(),         2.0);
    EXPECT_DOUBLE_EQ(p.getGain(),              3.0);
    EXPECT_DOUBLE_EQ(p.getQuantumEfficiency(), 0.4);
    EXPECT_DOUBLE_EQ(p.getKTCnoise(),          5.0);
    EXPECT_DOUBLE_EQ(p.getBaseline(),          6.0);
}

TEST(nispPixel, CopyConstructorPreservesBehaviour)
{
    // All noise sources zero except a fixed baseline ⇒ deterministic output.
    nisp::pixel original(0.0, 0.0, 1.0, 1.0, 0.0, 42.0);
    nisp::pixel copy(original);

    EXPECT_DOUBLE_EQ(copy.getBaseline(), original.getBaseline());
    EXPECT_DOUBLE_EQ(original.pixelresponse(0.0, 1.0), 42.0);
    EXPECT_DOUBLE_EQ(copy.pixelresponse(0.0, 1.0),     42.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  nisp::pixel — pixel stages (resetLevel / integrate / read) and pixelresponse()
// ═══════════════════════════════════════════════════════════════════════════════

TEST(nispPixel, PixelresponseIncludesResetLevel)
{
    // flux=0, dark=0, readNoise=0 ⇒ integrate and read add nothing.
    // KTC=5, baseline=100 ⇒ single read ~ N(100, 5).
    nisp::pixel p(0.0, 0.0, 1.0, 1.0, 5.0, 100.0);

    constexpr double k        = 5.0;
    constexpr size_t nSample  = 10'000;
    const double     n        = static_cast<double>(nSample);
    constexpr double mu       = 100.0;
    constexpr double var      = 25.0;

    double sigMean = 0.0, sigVar = 0.0;
    for (size_t i = 0; i < nSample; ++i)
    {
        const double x = p.pixelresponse(0.0, 1.0);
        sigMean += x;
        sigVar  += x * x;
    }
    sigMean /= n;
    sigVar   = sigVar / n - sigMean * sigMean;

    EXPECT_NEAR(sigMean, mu,  k * std::sqrt(var / n));
    EXPECT_NEAR(sigVar,  var, k * std::sqrt(2.0) * var / std::sqrt(n));
}

TEST(nispPixel, CfgOverloadMatchesTwoArgOverload)
{
    // The base pixelresponse(flux, t, cfg) ignores cfg: same single read as the
    // two-argument form, whatever the keys (the old per-pixel "baseline" switch
    // is gone, the reset level is its own stage: resetLevel()).
    nisp::pixel p(0.0, 0.0, 1.0, 1.0, 5.0, 100.0);
    PixelConfig cfg;
    cfg.set("baseline", std::string("off"));

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;
    const double     n       = static_cast<double>(nSample);

    double sigMean = 0.0, sigVar = 0.0;
    for (size_t i = 0; i < nSample; ++i)
    {
        const double x = p.pixelresponse(0.0, 1.0, cfg);
        sigMean += x;
        sigVar  += x * x;
    }
    sigMean /= n;
    sigVar   = sigVar / n - sigMean * sigMean;

    EXPECT_NEAR(sigMean, 100.0, k * std::sqrt(25.0 / n));
    EXPECT_NEAR(sigVar,  25.0,  k * std::sqrt(2.0) * 25.0 / std::sqrt(n));
}

TEST(nispPixel, ResetLevelIsBaselinePlusKtc)
{
    // resetLevel() = N(baseline, KTC/gain); exactly the baseline without kTC.
    nisp::pixel noKtc(0.0, 5.0, 2.0, 1.0, 0.0, 77.0);
    EXPECT_DOUBLE_EQ(noKtc.resetLevel(), 77.0);

    constexpr double KTC = 10.0, gain = 2.0, baseline = 100.0;
    constexpr double var = (KTC / gain) * (KTC / gain);   // 25 ADU²
    nisp::pixel p(0.0, 5.0, gain, 1.0, KTC, baseline);    // read noise must not appear

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;
    const double     n       = static_cast<double>(nSample);

    double sigMean = 0.0, sigVar = 0.0;
    for (size_t i = 0; i < nSample; ++i)
    {
        const double x = p.resetLevel();
        sigMean += x;
        sigVar  += x * x;
    }
    sigMean /= n;
    sigVar   = sigVar / n - sigMean * sigMean;

    EXPECT_NEAR(sigMean, baseline, k * std::sqrt(var / n));
    EXPECT_NEAR(sigVar,  var,      k * std::sqrt(2.0) * var / std::sqrt(n));
}

TEST(nispPixel, IntegrateHasNoResetLevelNorReadNoise)
{
    // No flux, no dark ⇒ no charge collected. The reset level (baseline 100,
    // kTC 5) and the read noise (5) belong to other stages: integrate() is
    // exactly 0, through both overloads.
    nisp::pixel p(0.0, 5.0, 1.0, 1.0, 5.0, 100.0);
    PixelConfig cfg;

    for (size_t i = 0; i < 1000; ++i)
    {
        EXPECT_DOUBLE_EQ(p.integrate(0.0, 1.0),      0.0);
        EXPECT_DOUBLE_EQ(p.integrate(0.0, 1.0, cfg), 0.0);
    }
}

TEST(nispPixel, FullNoiseBudgetViaCfgOverload)
{
    // Same setup as PixelSimulator.FullNoiseBudget in gTest_pixels.cxx, but
    // driven through the 3-arg cfg overload with an empty cfg, to confirm it
    // matches the 2-arg overload's behaviour exactly: one complete single read
    // = read(resetLevel() + integrate()).
    constexpr double flux      = 10.0;
    constexpr double dark      =  1.0;
    constexpr double readNoise =  5.0;
    constexpr double gain      =  2.0;
    constexpr double QE        =  0.75;
    constexpr double KTC       =  3.0;
    constexpr double baseline  = 100.0;
    constexpr double t         = 100.0;

    constexpr double signal_e = (flux * QE + dark) * t;                                    // 850 e
    constexpr double mu_ADU   = signal_e / gain + baseline;                                // 525 ADU
    constexpr double var_ADU  = (signal_e + readNoise * readNoise + KTC * KTC) / (gain * gain); // 221 ADU^2

    nisp::pixel p(dark, readNoise, gain, QE, KTC, baseline);
    PixelConfig cfg; // empty

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;
    const double     n       = static_cast<double>(nSample);

    double sigMean = 0.0, sigVar = 0.0;
    for (size_t i = 0; i < nSample; ++i)
    {
        const double x = p.pixelresponse(flux, t, cfg);
        sigMean += x;
        sigVar  += x * x;
    }
    sigMean /= n;
    sigVar   = sigVar / n - sigMean * sigMean;

    EXPECT_NEAR(sigMean, mu_ADU,  k * std::sqrt(var_ADU / n));
    EXPECT_NEAR(sigVar,  var_ADU, k * std::sqrt(2.0) * var_ADU / std::sqrt(n));
}

TEST(nispPixel, ReadAddsOnlyReadNoiseInADU)
{
    // read(charge) = charge + N(0, readNoise)/gain: mean charge, variance (RN/gain)².
    constexpr double readNoise = 5.0;
    constexpr double gain      = 2.0;
    constexpr double charge    = 1000.0;
    constexpr double var       = (readNoise / gain) * (readNoise / gain);   // 6.25 ADU²

    nisp::pixel p(0.0, readNoise, gain, 1.0, 0.0, 0.0);

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;
    const double     n       = static_cast<double>(nSample);

    double sigMean = 0.0, sigVar = 0.0;
    for (size_t i = 0; i < nSample; ++i)
    {
        const double x = p.read(charge);
        sigMean += x;
        sigVar  += x * x;
    }
    sigMean /= n;
    sigVar   = sigVar / n - sigMean * sigMean;

    EXPECT_NEAR(sigMean, charge, k * std::sqrt(var / n));
    EXPECT_NEAR(sigVar,  var,    k * std::sqrt(2.0) * var / std::sqrt(n));
}

TEST(nispPixel, ReadWithoutReadNoiseIsExact)
{
    nisp::pixel p(0.0, 0.0, 2.0, 1.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(p.read(123.25), 123.25);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  nisp::pixel — full well and non-linearity: read(charge, nl)
// ═══════════════════════════════════════════════════════════════════════════════

namespace
{
    constexpr float  kFullWell = 130000.0f;                  // e-
    constexpr float  kC2 = -6.0e-7f, kC3 = -3.1e-13f;        // 1/e-, 1/e-²
    constexpr double kGain = 1.92, kBaseline = 10000.0;      // e-/ADU, ADU

    // Noiseless pixel with a full well and the coefficients of the polynomial model.
    nisp::pixel makeWellPixel(double readNoise = 0.0)
    {
        nisp::pixel p(0.0, readNoise, kGain, 1.0, 0.0, kBaseline);
        p.setFullWell(kFullWell);
        p.setNlParam(0, kC2);
        p.setNlParam(1, kC3);
        return p;
    }

    // Accumulated charge (ADU, reset level included) of Q electrons.
    double charge(double Q) { return kBaseline + Q / kGain; }

    // Q·(1 + c2·Q + c3·Q²), with the coefficients the pixel holds (float).
    double polynomialResponse(double Q)
    {
        const double c2 = kC2, c3 = kC3;
        return Q + c2 * Q * Q + c3 * Q * Q * Q;
    }
}

TEST(nispPixel, NoFullWellByDefault)
{
    // Infinite full well: read(charge, nl) is read(charge), whatever the model.
    nisp::pixel p(0.0, 0.0, kGain, 1.0, 0.0, kBaseline);
    EXPECT_TRUE(std::isinf(p.getFullWell()));
    for (size_t k = 0; k < nonlinearity::nParam; ++k)
        EXPECT_FLOAT_EQ(p.getNlParam(k), 0.0f) << "parameter " << k;

    const nonlinearity::functor half = [](std::span<const float>, double Q) { return 0.5 * Q; };
    for (double Q : {0.0, 1.0e3, 1.0e6, 1.0e9})
        EXPECT_DOUBLE_EQ(p.read(charge(Q), half), p.read(charge(Q))) << "Q = " << Q;
}

TEST(nispPixel, ReadAppliesTheNonLinearityBelowTheFullWell)
{
    const nisp::pixel           p  = makeWellPixel();
    const nonlinearity::functor nl = nonlinearity::polynomial;

    for (double Q : {0.0, 1.0, 1.0e3, 2.0e4, 8.0e4, 1.299e5})
        EXPECT_NEAR(p.read(charge(Q), nl), charge(polynomialResponse(Q)), 1e-6) << "Q = " << Q;

    // These coefficients are those of a 5 % deviation at 80 ke-.
    const double measured = (p.read(charge(8.0e4), nl) - kBaseline) * kGain;
    EXPECT_NEAR(measured / 8.0e4 - 1.0, -0.05, 1e-4);
}

TEST(nispPixel, ReadIsConstantAboveTheFullWell)
{
    const nisp::pixel           p  = makeWellPixel();
    const nonlinearity::functor nl = nonlinearity::polynomial;

    const double plateau = charge(polynomialResponse(kFullWell));
    for (double Q : {1.0 * kFullWell, 1.001 * kFullWell, 2.0 * kFullWell, 1.0e7})
        EXPECT_NEAR(p.read(charge(Q), nl), plateau, 1e-6) << "Q = " << Q;
}

TEST(nispPixel, EmptyModelIsLinearUpToTheFullWell)
{
    const nisp::pixel           p = makeWellPixel();
    const nonlinearity::functor none;

    for (double Q : {1.0e3, 8.0e4})
        EXPECT_NEAR(p.read(charge(Q), none), charge(Q), 1e-6) << "Q = " << Q;
    for (double Q : {1.4e5, 1.0e7})
        EXPECT_NEAR(p.read(charge(Q), none), charge(kFullWell), 1e-6) << "Q = " << Q;
}

TEST(nispPixel, BaselineDoesNotFillTheWell)
{
    // The baseline is an electronic offset: a baseline larger than the full well
    // changes neither the non-linearity nor the plateau.
    const nonlinearity::functor nl = nonlinearity::polynomial;
    nisp::pixel low (0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
    nisp::pixel high(0.0, 0.0, 1.0, 1.0, 0.0, 50000.0);
    for (nisp::pixel* p : {&low, &high})
    {
        p->setFullWell(20000.0f);
        p->setNlParam(0, -2.0e-6f);
    }

    EXPECT_DOUBLE_EQ(high.read(50000.0, nl), 50000.0);
    for (double Q : {1.0e3, 1.5e4, 2.0e4, 1.0e5})
        EXPECT_NEAR(high.read(50000.0 + Q, nl) - 50000.0, low.read(Q, nl), 1e-6) << "Q = " << Q;
}

TEST(nispPixel, ChargeBelowTheBaselineIsLinear)
{
    // kTC noise can put the charge below the baseline. A model is not called there:
    // this one would return NaN.
    const nisp::pixel           p  = makeWellPixel();
    const nonlinearity::functor nl = [](std::span<const float>, double Q) { return std::exp(std::log(Q)); };

    EXPECT_NEAR(p.read(kBaseline - 50.0, nl), kBaseline - 50.0, 1e-9);
    EXPECT_NEAR(p.read(kBaseline,        nl), kBaseline,        1e-9);
}

TEST(nispPixel, ModelReceivesTheParametersOfThePixelAndItsCaptures)
{
    nisp::pixel p;
    p.setFullWell(20000.0f);
    const std::vector<float> values{0.125f, 0.5f, 0.0625f};
    ASSERT_EQ(values.size(), nonlinearity::nParam);
    for (size_t k = 0; k < values.size(); ++k)
        p.setNlParam(k, values[k]);

    const double captured = 0.25;
    size_t       nParam   = 0;
    const nonlinearity::functor f = [captured, &nParam](std::span<const float> param, double Q)
    {
        nParam = param.size();
        return Q * (captured + param[0] + param[1] + param[2]);   // Q × 0.9375
    };

    EXPECT_DOUBLE_EQ(p.nlResponse(1000.0, f), 937.5);
    EXPECT_EQ(nParam, nonlinearity::nParam);
}

TEST(nispPixel, CopiesKeepTheFullWellAndTheNonLinearity)
{
    const nisp::pixel           original = makeWellPixel();
    const nonlinearity::functor nl       = nonlinearity::polynomial;

    const nisp::pixel copy(original);
    nisp::pixel       assigned;
    assigned = original;

    for (const nisp::pixel* p : {&copy, static_cast<const nisp::pixel*>(&assigned)})
    {
        EXPECT_FLOAT_EQ(p->getFullWell(), kFullWell);
        EXPECT_FLOAT_EQ(p->getNlParam(0), kC2);
        EXPECT_FLOAT_EQ(p->getNlParam(1), kC3);
        for (double Q : {8.0e4, 1.0e7})
            EXPECT_DOUBLE_EQ(p->read(charge(Q), nl), original.read(charge(Q), nl)) << "Q = " << Q;
    }
}

TEST(nispPixel, InvalidFullWellOrParameterThrows)
{
    nisp::pixel p = makeWellPixel();
    EXPECT_THROW(p.setFullWell(0.0f),                                    std::invalid_argument);
    EXPECT_THROW(p.setFullWell(-1.0f),                                   std::invalid_argument);
    EXPECT_THROW(p.setFullWell(std::numeric_limits<float>::quiet_NaN()), std::invalid_argument);
    EXPECT_FLOAT_EQ(p.getFullWell(), kFullWell) << "a rejected value must not be kept";

    EXPECT_THROW(p.setNlParam(nonlinearity::nParam, 0.0f), std::out_of_range);
}

TEST(nispPixel, StoredChargeIsLimitedByTheFullWell)
{
    const nisp::pixel p = makeWellPixel();
    EXPECT_DOUBLE_EQ(p.storedCharge(1000.0), 1000.0);
    EXPECT_DOUBLE_EQ(p.storedCharge(1.0e7),  static_cast<double>(kFullWell));

    const nisp::pixel none;
    EXPECT_DOUBLE_EQ(none.storedCharge(1.0e7), 1.0e7);
}

TEST(nispPixel, ReadNoiseIsNotCompressed)
{
    // The read noise is added after the full well and the non-linearity: at any charge,
    // below or above the full well, the read variance is (RN/gain)².
    constexpr double readNoise = 13.0;
    constexpr double var       = (readNoise / kGain) * (readNoise / kGain);
    const nisp::pixel           p  = makeWellPixel(readNoise);
    const nonlinearity::functor nl = nonlinearity::polynomial;

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;
    const double     n       = static_cast<double>(nSample);

    for (double Q : {8.0e4, 1.0e7})
    {
        const double expected = charge(polynomialResponse(std::min<double>(Q, kFullWell)));

        double sigMean = 0.0, sigVar = 0.0;
        for (size_t i = 0; i < nSample; ++i)
        {
            const double x = p.read(charge(Q), nl) - expected;
            sigMean += x;
            sigVar  += x * x;
        }
        sigMean /= n;
        sigVar   = sigVar / n - sigMean * sigMean;

        EXPECT_NEAR(sigMean, 0.0, k * std::sqrt(var / n))                  << "Q = " << Q;
        EXPECT_NEAR(sigVar,  var, k * std::sqrt(2.0) * var / std::sqrt(n)) << "Q = " << Q;
    }
}

TEST(nispPixel, NonLinearityModelFromItsName)
{
    const nonlinearity::functor f = nonlinearity::fromString("polynomial");
    ASSERT_TRUE(f);

    const std::vector<float> param{kC2, kC3, 0.0f};
    EXPECT_DOUBLE_EQ(f(param, 8.0e4), nonlinearity::polynomial(param, 8.0e4));
    EXPECT_NEAR     (f(param, 8.0e4), polynomialResponse(8.0e4), 1e-6);

    EXPECT_THROW(nonlinearity::fromString("unknown"), std::invalid_argument);
    EXPECT_THROW(nonlinearity::fromString(""),        std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  nisp::Ramp — coadd() / getSlope() / getFowler() on noise-free linear ramps
// ═══════════════════════════════════════════════════════════════════════════════

TEST(nispRamp, CoaddNoiselessLinearRampMatchesKubicFormula)
{
    // For a perfectly linear, noise-free UTR ramp (frame[i] = s*i), every
    // coadded CDS difference equals the same constant c = s*(nreads+ndrops).
    // With fsigma=0 (beta=0), Ramp::getSlope() must then report chi2 == 0
    // exactly (size*pseudo_flux == lf identically, regardless of alpha/gamma),
    // and its "flux" output must match the Kubic-formula closed form for a
    // constant CDS step, which does exercise alpha/gamma.
    struct RampCase { double s; size_t nreads; size_t ndrops; size_t ngr; double gain; };
    const std::vector<RampCase> cases = {
        {1.0,  2, 0, 5,  1.0},
        {3.7,  4, 1, 8,  2.5},
        {0.5,  8, 3, 6,  1.0},
        {10.0, 1, 0, 10, 4.0},
    };

    for (const auto& c : cases)
    {
        const size_t N = c.ngr * (c.nreads + c.ndrops);
        std::valarray<double> utr = makeLinearRamp(c.s, N);

        nisp::Ramp  raw(&utr, std::slice(0, N, 1), 0.0, c.gain);
        nisp::Ramp& macc = raw.coadd(c.ngr, c.nreads, c.ndrops);

        double flux = 0.0, chi2 = 0.0, err = 0.0;
        macc.getSlope(flux, chi2, err);

        EXPECT_NEAR(chi2, 0.0, 1e-6)
            << "s=" << c.s << " nreads=" << c.nreads << " ndrops=" << c.ndrops
            << " ngr=" << c.ngr << " gain=" << c.gain;

        const double step         = c.s * static_cast<double>(c.nreads + c.ndrops);
        const double expectedFlux = expectedKubicFlux(step, c.gain, c.nreads, c.ndrops);
        EXPECT_NEAR(flux, expectedFlux, 1e-6 * std::max(1.0, std::abs(expectedFlux)))
            << "s=" << c.s << " nreads=" << c.nreads << " ndrops=" << c.ndrops
            << " ngr=" << c.ngr << " gain=" << c.gain;

        delete &macc;
    }
}

TEST(nispRamp, CoaddIntoBufferMatchesHeapCoadd)
{
    // The allocation-free coadd(groups, ...) used per pixel by SCA::expose_dpu must
    // give exactly the same groups, settings and slope as the heap-allocating coadd().
    constexpr size_t nreads = 4, ndrops = 1, ngr = 6, N = ngr * (nreads + ndrops);
    std::valarray<double> utr(N);
    for (size_t i = 0; i < N; ++i)
        utr[i] = 3.0 * static_cast<double>(i) + 0.1 * static_cast<double>(i * i);   // non-linear

    nisp::Ramp raw(&utr, std::slice(0, N, 1), 2.0, 1.5);

    nisp::Ramp& heap = raw.coadd(ngr, nreads, ndrops);

    nisp::rawLines groups(0.0, 2);   // too small on purpose: coadd must resize it
    nisp::Ramp     buffered = raw.coadd(groups, ngr, nreads, ndrops);

    ASSERT_EQ(groups.size(), ngr);
    const nisp::rawLines heapGroups = heap.getRamp();
    for (size_t g = 0; g < ngr; ++g)
        EXPECT_DOUBLE_EQ(groups[g], heapGroups[g]) << "group " << g;

    EXPECT_EQ(buffered.get_reads_per_group(), nreads);
    EXPECT_EQ(buffered.get_drops_per_group(), ndrops);

    double f1 = 0., c1 = 0., e1 = 0., f2 = 0., c2 = 0., e2 = 0.;
    heap.getSlope(f1, c1, e1);
    buffered.getSlope(f2, c2, e2);
    EXPECT_DOUBLE_EQ(f1, f2);
    EXPECT_DOUBLE_EQ(c1, c2);
    EXPECT_DOUBLE_EQ(e1, e2);

    delete &heap;
}

TEST(nispRamp, GetSlopeThrowsOnSingleGroup)
{
    std::valarray<double> data = {5.0};
    nisp::Ramp r(&data, std::slice(0, 1, 1));
    r.reads_per_group(1);

    double flux = 0., chi2 = 0., err = 0.;
    EXPECT_THROW(r.getSlope(flux, chi2, err), std::runtime_error);
}

TEST(nispRamp, GetFowlerOnRawRamp)
{
    std::valarray<double> data = {2.0, 5.0, 3.0, 9.0};
    nisp::Ramp r(&data, std::slice(0, 4, 1));

    double flux = 0.0;
    r.getFowler(flux);
    EXPECT_DOUBLE_EQ(flux, 9.0 - 2.0);
}

TEST(nispRamp, GetFowlerOnCoaddedRamp)
{
    constexpr double s       = 2.0;
    constexpr size_t nreads  = 3;
    constexpr size_t ndrops  = 1;
    constexpr size_t ngr     = 4;
    constexpr size_t N       = ngr * (nreads + ndrops);

    std::valarray<double> utr = makeLinearRamp(s, N);
    nisp::Ramp  raw(&utr, std::slice(0, N, 1));
    nisp::Ramp& macc = raw.coadd(ngr, nreads, ndrops);

    double flux = 0.0;
    macc.getFowler(flux);

    // Telescoping sum of (ngr-1) equal steps of size s*(nreads+ndrops).
    const double expected = static_cast<double>(ngr - 1) * s * static_cast<double>(nreads + ndrops);
    EXPECT_DOUBLE_EQ(flux, expected);

    delete &macc;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  nisp::Ramp — getCDS()
// ═══════════════════════════════════════════════════════════════════════════════

TEST(nispRamp, GetCDSStandardBranchMatchesHandComputation)
{
    // Non-linear ramp (squares), default ndrops=0 <= 1 ⇒ standard whole-ramp
    // consecutive-difference branch.
    std::valarray<double> data = {1.0, 4.0, 9.0, 16.0, 25.0}; // diffs: 3,5,7,9
    nisp::Ramp r(&data, std::slice(0, 5, 1));

    double cds = 0.0, noise = 0.0;
    r.getCDS(cds, noise);

    EXPECT_NEAR(cds,   6.0,            1e-9); // mean(3,5,7,9)
    EXPECT_NEAR(noise, std::sqrt(20.0/3.0), 1e-9); // stddev(ddof=1)
}

TEST(nispRamp, GetCDSThrowsOnTooShortRamp)
{
    std::valarray<double> data = {5.0};
    nisp::Ramp r(&data, std::slice(0, 1, 1));

    double cds = 0.0, noise = 0.0;
    EXPECT_THROW(r.getCDS(cds, noise), std::runtime_error);
}

TEST(nispRamp, GetCDSMaccBranchCharacterization)
{
    // Characterization test: this locks in CURRENT behaviour of the
    // `else if (nreads > 1)` branch, which is reached here because ndrops=2
    // fails the `ndrops <= 1` guard. NOTE: the loop advances each group-start
    // offset by 1 per iteration rather than by `nreads`, which contradicts
    // its own doc comment ("within each group of nreads frames..."). This is
    // NOT asserted to be "correct" — only that it doesn't silently change.
    std::valarray<double> data(12);
    for (size_t i = 0; i < 12; ++i)
        data[i] = static_cast<double>(i * i); // 0,1,4,9,16,25,36,49,64,81,100,121

    nisp::Ramp r(&data, std::slice(0, 12, 1));
    r.reads_per_group(4);
    r.drops_per_group(2);

    double cds = 0.0, noise = 0.0;
    r.getCDS(cds, noise);

    // Hand-derived from the current loop (iS=0..2, k=0..2 per iS):
    // diffs = [1,3,5, 3,5,7, 5,7,9] -> mean=5.0, stddev(ddof=1)=sqrt(6.0)
    EXPECT_NEAR(cds,   5.0,            1e-9);
    EXPECT_NEAR(noise, std::sqrt(6.0), 1e-9);
}

TEST(nispRamp, AccessorsAndSliceRoundTrip)
{
    std::valarray<double> data = {10.0, 20.0, 30.0, 40.0, 50.0};
    std::slice sl(1, 3, 1); // elements 20,30,40

    nisp::Ramp r(&data, sl, 2.5, 1.5);
    r.reads_per_group(4);
    r.drops_per_group(2);

    EXPECT_EQ(r.get_reads_per_group(), 4u);
    EXPECT_EQ(r.get_drops_per_group(), 2u);

    EXPECT_EQ(r.get_slice().start(),  sl.start());
    EXPECT_EQ(r.get_slice().size(),   sl.size());
    EXPECT_EQ(r.get_slice().stride(), sl.stride());

    nisp::rawLines got = r.getRamp();
    ASSERT_EQ(got.size(), 3u);
    EXPECT_DOUBLE_EQ(got[0], 20.0);
    EXPECT_DOUBLE_EQ(got[1], 30.0);
    EXPECT_DOUBLE_EQ(got[2], 40.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  nisp::Ramp — statistical consistency of the reported flux uncertainty
// ═══════════════════════════════════════════════════════════════════════════════

TEST(nispRamp, CoaddStatisticalErrorMatchesObservedScatter)
{
    // Build many independent noisy realizations of the same underlying ramp
    // (constant slope s, per-frame read noise sigma) and check that the
    // *reported* uncertainty `err` from getSlope() is consistent, in order of
    // magnitude, with the *observed* trial-to-trial scatter of `flux`. This
    // validates the alpha/beta/gamma error-propagation path without needing
    // to independently re-derive the Kubic et al. bias-corrected flux formula
    // in closed form.
    constexpr double s       = 8.0;
    constexpr double gain    = 1.0;
    constexpr double sigma   = 3.0;
    constexpr size_t nreads  = 4;
    constexpr size_t ndrops  = 1;
    constexpr size_t ngr     = 30;
    constexpr size_t nFrames = ngr * (nreads + ndrops);
    constexpr size_t nTrials = 1500;

    std::vector<double> fluxSamples;
    std::vector<double> errSamples;
    fluxSamples.reserve(nTrials);
    errSamples.reserve(nTrials);

    for (size_t trial = 0; trial < nTrials; ++trial)
    {
        std::valarray<double> utr(nFrames);
        double cumulative = 0.0;
        for (size_t i = 0; i < nFrames; ++i)
        {
            cumulative += s + DST::Math::sample_normal_distribution(0.0, sigma);
            utr[i] = cumulative;
        }

        nisp::Ramp  raw(&utr, std::slice(0, nFrames, 1), sigma, gain);
        nisp::Ramp& macc = raw.coadd(ngr, nreads, ndrops);

        double flux = 0.0, chi2 = 0.0, err2 = 0.0;
        macc.getSlope(flux, chi2, err2);

        fluxSamples.push_back(flux);
        errSamples.push_back(err2);

        delete &macc;
    }

    const double n = static_cast<double>(nTrials);
    double meanFlux = 0.0, meanErr = 0.0;
    for (size_t i = 0; i < nTrials; ++i) { meanFlux += fluxSamples[i]; meanErr += errSamples[i]; }
    meanFlux /= n;
    meanErr  /= n;

    double varFlux = 0.0;
    for (double f : fluxSamples) varFlux += (f - meanFlux) * (f - meanFlux);
    varFlux /= (n - 1.0);
    const double empiricalStd = varFlux;

    // Broad factor-of-3 band rather than a tight relative tolerance: this
    // checks `err` is in the right ballpark (not zero, not orders of
    // magnitude off), without assuming an exact calibration constant that
    // would need to be verified empirically against a real build.
    EXPECT_GT(meanErr, 0.0);
    EXPECT_GT(empiricalStd, 0.3 * meanErr)
        << "empiricalStd=" << empiricalStd << " meanErr=" << meanErr;
    EXPECT_LT(empiricalStd, 3.0 * meanErr)
        << "empiricalStd=" << empiricalStd << " meanErr=" << meanErr;
}
