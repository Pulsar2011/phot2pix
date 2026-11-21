//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_persistence.cxx
//

#include <gtest/gtest.h>
#include <nisp/nispSCA.h>
#include <DSTfits/FITSimg.h>

#include <cmath>
#include <memory>
#include <valarray>

using namespace nisp::persistence;
using namespace phot2pix;

// ── Helpers ──────────────────────────────────────────────────────────────────

// Build a PixelConfig carrying a single trap component in array format.
static PixelConfig makeTrapConfig(double f, double tau)
{
    PixelConfig cfg;
    nlohmann::json traps = nlohmann::json::array();
    traps.push_back({{"f", f}, {"t", tau}});
    cfg.set("persistence", traps);
    return cfg;
}

// Create a uniform DOUBLE_IMG FITScube of given dimensions.
static std::shared_ptr<DSL::FITScube> makeFluxImage(size_t w, size_t h, double value)
{
    std::shared_ptr<DSL::FITScube> img(new DSL::FITSimg<double>({w, h}));
    *(img) += std::valarray<double>(value, w * h);
    return img;
}

// Expected persistence signal for a single trap component: the charge released over
// the step [dT, dT + T] by a charge `flux` stored at t = 0, at the rate f·e^(−t/τ).
static double expectedSignal(double f, double tau, double dT, double T, double flux)
{
    return f * tau * std::exp(-dT / tau) * (1.0 - std::exp(-T / tau)) * flux;
}

// ═════════════════════════════════════════════════════════════════════════════
//  persistence::model
// ═════════════════════════════════════════════════════════════════════════════

TEST(PersistenceModel, DefaultConstruction)
{
    model m;
    EXPECT_TRUE(m.getFactor().empty());
    EXPECT_TRUE(m.getDecayTime().empty());
    EXPECT_DOUBLE_EQ(m.getTimeSinceLastExposure(), 0.0);
}

TEST(PersistenceModel, InsertAndClearTraps)
{
    model m;
    m.insertParameter(0.01, 100.0);
    m.insertParameter(0.005, 500.0);

    ASSERT_EQ(m.getFactor().size(),    2u);
    ASSERT_EQ(m.getDecayTime().size(), 2u);
    EXPECT_DOUBLE_EQ(m.getFactor().at(0),    0.01);
    EXPECT_DOUBLE_EQ(m.getDecayTime().at(0), 100.0);
    EXPECT_DOUBLE_EQ(m.getFactor().at(1),    0.005);
    EXPECT_DOUBLE_EQ(m.getDecayTime().at(1), 500.0);

    m.clearTraps();
    EXPECT_TRUE(m.getFactor().empty());
    EXPECT_TRUE(m.getDecayTime().empty());
}

TEST(PersistenceModel, IncrementTime)
{
    model m;
    m.incrementTime(50.0);
    EXPECT_DOUBLE_EQ(m.getTimeSinceLastExposure(), 50.0);
    m.incrementTime(30.0);
    EXPECT_DOUBLE_EQ(m.getTimeSinceLastExposure(), 80.0);
}

TEST(PersistenceModel, ZeroResponseNoTraps)
{
    model m;
    EXPECT_DOUBLE_EQ(m.pixelresponse(1000.0, 10.0), 0.0);
}

TEST(PersistenceModel, ZeroResponseZeroFlux)
{
    model m;
    m.insertParameter(0.01, 100.0);
    m.incrementTime(50.0);
    EXPECT_DOUBLE_EQ(m.pixelresponse(0.0, 10.0), 0.0);
}

TEST(PersistenceModel, SingleTrapResponse)
{
    const double f    = 0.01;
    const double tau  = 100.0;
    const double dT   = 50.0;
    const double T    = 10.0;
    const double flux = 1000.0;

    model m;
    m.insertParameter(f, tau);
    m.incrementTime(dT);

    const double exp = expectedSignal(f, tau, dT, T, flux);
    EXPECT_NEAR(m.pixelresponse(flux, T), exp, 1e-9 * exp);
}

TEST(PersistenceModel, CfgDTOverridesMemberDT)
{
    const double f        = 0.01;
    const double tau      = 100.0;
    const double memberDT = 50.0;
    const double cfgDT    = 200.0;
    const double T        = 10.0;
    const double flux     = 1000.0;

    model m;
    m.insertParameter(f, tau);
    m.incrementTime(memberDT);

    PixelConfig cfg;
    cfg.set("dT", cfgDT);

    const double exp = expectedSignal(f, tau, cfgDT, T, flux);
    EXPECT_NEAR(m.pixelresponse(flux, T, cfg), exp, 1e-9 * exp);
}

TEST(PersistenceModel, TwoTrapResponse)
{
    const double f1 = 0.01,  tau1 = 100.0;
    const double f2 = 0.005, tau2 = 500.0;
    const double dT   = 50.0;
    const double T    = 10.0;
    const double flux = 1000.0;

    model m;
    m.insertParameter(f1, tau1);
    m.insertParameter(f2, tau2);
    m.incrementTime(dT);

    const double exp = expectedSignal(f1, tau1, dT, T, flux)
                     + expectedSignal(f2, tau2, dT, T, flux);
    EXPECT_NEAR(m.pixelresponse(flux, T), exp, 1e-9 * std::abs(exp));
}

TEST(PersistenceModel, ConsecutiveStepsAddUpToOneStep)
{
    // The release over [Δ, Δ + T] is the sum of the releases over [Δ, Δ + T/2] and
    // [Δ + T/2, Δ + T]: splitting a step does not change the charge released.
    const double f = 0.01, tau = 100.0, flux = 1000.0, delta = 20.0, T = 30.0;

    model m;
    m.insertParameter(f, tau);

    auto release = [&m, flux](double start, double length)
    {
        PixelConfig cfg;
        cfg.set("dT", start);
        return m.pixelresponse(flux, length, cfg);
    };

    const double whole = release(delta, T);
    EXPECT_NEAR(whole, release(delta, T / 2) + release(delta + T / 2, T / 2), 1e-9 * whole);
}

TEST(PersistenceModel, TotalReleaseIsFTau)
{
    // Over all the time after the exposure, a trap releases f·τ of the stored charge.
    model m;
    m.insertParameter(1e-3, 80.0);

    PixelConfig cfg;
    cfg.set("dT", 0.0);
    EXPECT_NEAR(m.pixelresponse(1.0, 1e6, cfg), 1e-3 * 80.0, 1e-12);
}

// ═════════════════════════════════════════════════════════════════════════════
//  persistence::pixelArray
// ═════════════════════════════════════════════════════════════════════════════

TEST(PersistencePixelArray, Construction)
{
    pixelArray arr(8, 4);
    EXPECT_EQ(arr.getWidth(),  8u);
    EXPECT_EQ(arr.getHeight(), 4u);
    EXPECT_EQ(arr.size(),      32u);
}

TEST(PersistencePixelArray, InitializeArrayFormat)
{
    const size_t W = 4, H = 4;
    pixelArray arr(W, H);
    arr.initialize(makeTrapConfig(0.01, 100.0));

    PixelConfig simCfg;
    simCfg.set("dT", 50.0);
    std::valarray<double> flux(1000.0, W * H);

    std::shared_ptr<DSL::FITScube> result = arr.simulate(flux, 10.0, simCfg);
    ASSERT_NE(result, nullptr);

    const std::valarray<double>& data = *(result->GetData<double>());
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_GT(data[i], 0.0) << "Pixel " << i << " should have non-zero persistence";
}

TEST(PersistencePixelArray, InitializeObjectFormat)
{
    const size_t W = 4, H = 4;
    pixelArray arr(W, H);

    PixelConfig cfg;
    nlohmann::json persistObj;
    nlohmann::json trapsArr = nlohmann::json::array();
    trapsArr.push_back({{"f", 0.01}, {"t", 100.0}});
    persistObj["traps"] = trapsArr;
    cfg.set("persistence", persistObj);
    arr.initialize(cfg);

    PixelConfig simCfg;
    simCfg.set("dT", 50.0);
    std::valarray<double> flux(1000.0, W * H);

    std::shared_ptr<DSL::FITScube> result = arr.simulate(flux, 10.0, simCfg);
    ASSERT_NE(result, nullptr);

    const std::valarray<double>& data = *(result->GetData<double>());
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_GT(data[i], 0.0) << "Pixel " << i;
}

TEST(PersistencePixelArray, ArrayAndObjectFormatAgreement)
{
    const size_t W = 4, H = 4;
    const double f = 0.01, tau = 100.0;

    pixelArray arrA(W, H);
    arrA.initialize(makeTrapConfig(f, tau));

    PixelConfig cfgObj;
    nlohmann::json persistObj;
    nlohmann::json trapsArr = nlohmann::json::array();
    trapsArr.push_back({{"f", f}, {"t", tau}});
    persistObj["traps"] = trapsArr;
    cfgObj.set("persistence", persistObj);
    pixelArray arrB(W, H);
    arrB.initialize(cfgObj);

    PixelConfig simCfg;
    simCfg.set("dT", 50.0);
    std::valarray<double> flux(1000.0, W * H);

    std::shared_ptr<DSL::FITScube> rA = arrA.simulate(flux, 10.0, simCfg);
    std::shared_ptr<DSL::FITScube> rB = arrB.simulate(flux, 10.0, simCfg);

    const std::valarray<double>& dA = *(rA->GetData<double>());
    const std::valarray<double>& dB = *(rB->GetData<double>());
    for (size_t i = 0; i < dA.size(); ++i)
        EXPECT_NEAR(dA[i], dB[i], 1e-12) << "Array and object formats disagree at pixel " << i;
}

TEST(PersistencePixelArray, IncrementTimeDecaysSignal)
{
    const size_t W = 2, H = 2;
    pixelArray arr(W, H);
    arr.initialize(makeTrapConfig(0.01, 100.0));

    std::valarray<double> flux(1000.0, W * H);

    // Without incrementTime: member dT = 0
    std::shared_ptr<DSL::FITScube> r0 = arr.simulate(flux, 10.0);

    arr.incrementTime(50.0);

    // After incrementTime(50): member dT = 50
    std::shared_ptr<DSL::FITScube> r50 = arr.simulate(flux, 10.0);

    const double s0  = (*r0 ->GetData<double>())[0];
    const double s50 = (*r50->GetData<double>())[0];
    EXPECT_GT(s0, s50) << "More elapsed time should reduce persistence signal";
}

TEST(PersistencePixelArray, SimulateMatchesModelFormula)
{
    const double f    = 0.01;
    const double tau  = 100.0;
    const double dT   = 50.0;
    const double T    = 10.0;
    const double flux = 1000.0;
    const size_t W = 4, H = 4;

    pixelArray arr(W, H);
    arr.initialize(makeTrapConfig(f, tau));

    PixelConfig simCfg;
    simCfg.set("dT", dT);
    std::valarray<double> fluxArr(flux, W * H);

    std::shared_ptr<DSL::FITScube> result = arr.simulate(fluxArr, T, simCfg);
    ASSERT_NE(result, nullptr);

    const double exp = expectedSignal(f, tau, dT, T, flux);
    const std::valarray<double>& data = *(result->GetData<double>());
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_NEAR(data[i], exp, 1e-9 * exp) << "Pixel " << i;
}

// ═════════════════════════════════════════════════════════════════════════════
//  persistence::scheduler
// ═════════════════════════════════════════════════════════════════════════════

TEST(PersistenceScheduler, EmptyReturnsNullptr)
{
    PixelConfig cfg = makeTrapConfig(0.01, 100.0);
    scheduler sched(cfg);
    EXPECT_EQ(sched.computePersistenceSignal(50.0, 10.0), nullptr);
}

TEST(PersistenceScheduler, DimensionMismatchThrows)
{
    PixelConfig cfg = makeTrapConfig(0.01, 100.0);
    scheduler sched(cfg);
    sched.insertCharges(makeFluxImage(4, 4, 1000.0), 0.0);
    EXPECT_THROW(sched.insertCharges(makeFluxImage(8, 4, 1000.0), 0.0),
                 std::invalid_argument);
}

TEST(PersistenceScheduler, SingleExposureMatchesFormula)
{
    const double f    = 0.01;
    const double tau  = 100.0;
    const double dT   = 50.0;
    const double T    = 10.0;
    const double flux = 1000.0;
    const size_t W = 4, H = 4;

    PixelConfig cfg = makeTrapConfig(f, tau);
    scheduler sched(cfg);
    sched.insertCharges(makeFluxImage(W, H, flux), 0.0);

    std::shared_ptr<DSL::FITScube> result = sched.computePersistenceSignal(dT, T);
    ASSERT_NE(result, nullptr);

    const double exp = expectedSignal(f, tau, dT, T, flux);
    const std::valarray<double>& data = *(result->GetData<double>());
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_NEAR(data[i], exp, 1e-9 * exp) << "Pixel " << i;
}

TEST(PersistenceScheduler, T0AccumulatesAfterCompute)
{
    const size_t W = 4, H = 4;
    PixelConfig cfg = makeTrapConfig(0.01, 100.0);
    scheduler sched(cfg);
    sched.insertCharges(makeFluxImage(W, H, 1000.0), 0.0);
    EXPECT_NEAR(sched.getHistory()[0].t0, 0.0, 1e-9);

    // First compute: record.t0 starts at 0, dT seen = 0+50 = 50, then t0 += 10 → 10
    std::shared_ptr<DSL::FITScube> r1 = sched.computePersistenceSignal(50.0, 10.0);
    EXPECT_NEAR(sched.getHistory()[0].t0, 60.0, 1e-9);

    // Second compute: record.t0 = 10, dT seen = 10+50 = 60, then t0 += 10 → 20
    std::shared_ptr<DSL::FITScube> r2 = sched.computePersistenceSignal(50.0, 10.0);
    EXPECT_NEAR(sched.getHistory()[0].t0, 120.0, 1e-9);

    ASSERT_NE(r1, nullptr);
    ASSERT_NE(r2, nullptr);

    const double s1 = (*r1->GetData<double>())[0];
    const double s2 = (*r2->GetData<double>())[0];
    EXPECT_LT(s2, s1) << "Signal should decrease as elapsed time accumulates";
}

TEST(PersistenceScheduler, TwoIdenticalExposuresDoubleSignal)
{
    const double f    = 0.01;
    const double tau  = 100.0;
    const double dT   = 50.0;
    const double T    = 10.0;
    const double flux = 1000.0;
    const size_t W = 4, H = 4;

    PixelConfig cfg = makeTrapConfig(f, tau);

    scheduler schedOne(cfg);
    schedOne.insertCharges(makeFluxImage(W, H, flux), 0.0);

    scheduler schedTwo(cfg);
    schedTwo.insertCharges(makeFluxImage(W, H, flux), 0.0);
    schedTwo.insertCharges(makeFluxImage(W, H, flux), 0.0);

    std::shared_ptr<DSL::FITScube> rOne = schedOne.computePersistenceSignal(dT, T);
    std::shared_ptr<DSL::FITScube> rTwo = schedTwo.computePersistenceSignal(dT, T);

    ASSERT_NE(rOne, nullptr);
    ASSERT_NE(rTwo, nullptr);

    const double sOne = (*rOne->GetData<double>())[0];
    const double sTwo = (*rTwo->GetData<double>())[0];
    EXPECT_NEAR(sTwo, 2.0 * sOne, 1e-9 * sTwo)
        << "Two identical exposures should sum to twice the single-exposure signal";
}

TEST(PersistenceScheduler, RemoveCharges)
{
    const size_t W = 4, H = 4;
    PixelConfig cfg = makeTrapConfig(0.01, 100.0);

    scheduler schedFull(cfg);
    schedFull.insertCharges(makeFluxImage(W, H, 1000.0), 0.0);
    schedFull.insertCharges(makeFluxImage(W, H, 1000.0), 5.0);

    // Remove the first entry (t0=0)
    schedFull.removeCharges(0.0);

    // Reference: only the second entry (t0=5)
    scheduler schedRef(cfg);
    schedRef.insertCharges(makeFluxImage(W, H, 1000.0), 5.0);

    std::shared_ptr<DSL::FITScube> rFull = schedFull.computePersistenceSignal(50.0, 10.0);
    std::shared_ptr<DSL::FITScube> rRef  = schedRef .computePersistenceSignal(50.0, 10.0);

    ASSERT_NE(rFull, nullptr);
    ASSERT_NE(rRef,  nullptr);

    const double sFull = (*rFull->GetData<double>())[0];
    const double sRef  = (*rRef ->GetData<double>())[0];
    EXPECT_NEAR(sFull, sRef, 1e-9 * sRef)
        << "After removal, remaining signal should match single-entry reference";
}

// Persistence configuration in object form, with a delay between exposures.
static PixelConfig makeDelayedTrapConfig(double f, double tau, double delay)
{
    PixelConfig cfg;
    nlohmann::json p;
    p["traps"] = nlohmann::json::array({ {{"f", f}, {"t", tau}} });
    p["dT"]    = delay;
    cfg.set("persistence", p);
    return cfg;
}

TEST(PersistenceScheduler, DelayBetweenExposuresAgesEveryStoredExposure)
{
    // "dT": 30 → the next exposure starts 30 s after the end of each stored one.
    const double f = 0.01, tau = 100.0, flux = 1000.0, T = 10.0, delay = 30.0;
    const size_t W = 4, H = 4;

    PixelConfig cfg = makeDelayedTrapConfig(f, tau, delay);
    scheduler sched(cfg);
    EXPECT_DOUBLE_EQ(sched.getDelay(), delay);

    sched.insertCharges(makeFluxImage(W, H, flux), 0.0);
    EXPECT_NEAR(sched.getHistory()[0].t0, delay, 1e-9);

    // The first step of the next exposure starts `delay` after the stored one.
    std::shared_ptr<DSL::FITScube> r = sched.computePersistenceSignal(0.0, T);
    ASSERT_NE(r, nullptr);
    const double exp = expectedSignal(f, tau, delay, T, flux);
    EXPECT_NEAR((*r->GetData<double>())[0], exp, 1e-9 * exp);

    // A second exposure is stored: both exposures age by the delay before the next one.
    sched.insertCharges(makeFluxImage(W, H, flux), 0.0);
    EXPECT_NEAR(sched.getHistory()[0].t0, delay + T + delay, 1e-9);
    EXPECT_NEAR(sched.getHistory()[1].t0, delay,             1e-9);
}

TEST(PersistenceScheduler, NoDelayByDefault)
{
    PixelConfig cfg = makeTrapConfig(0.01, 100.0);   // array form: no "dT"
    scheduler sched(cfg);                            // keeps a reference to cfg
    EXPECT_DOUBLE_EQ(sched.getDelay(), 0.0);
}

TEST(PersistenceScheduler, NegativeDelayThrows)
{
    PixelConfig cfg = makeDelayedTrapConfig(0.01, 100.0, -1.0);
    EXPECT_THROW(scheduler{cfg}, std::invalid_argument);
}
