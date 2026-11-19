//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_pixels.cxx
//

#include <gtest/gtest.h>
#include <algorithm>              // std::min
#include <phot2pix/PixelSimulator.h>

using namespace phot2pix;

// ═══════════════════════════════════════════════════════════════════════════════
//  pixsim — construction, copy semantics, parameter validation
// ═══════════════════════════════════════════════════════════════════════════════

TEST(pixsim, DefaultConstructor)
{
    pixsim sim = pixsim();
    EXPECT_DOUBLE_EQ(sim.getDarkCurrent(),       0.0);
    EXPECT_DOUBLE_EQ(sim.getReadNoise(),         0.0);
    EXPECT_DOUBLE_EQ(sim.getGain(),              1.0);
    EXPECT_DOUBLE_EQ(sim.getKTCnoise(),          0.0);
    EXPECT_DOUBLE_EQ(sim.getQuantumEfficiency(), 1.0);
    EXPECT_DOUBLE_EQ(sim.getBaseline(),          0.0);
}

TEST(pixsim, Constructor)
{
    pixsim sim(1.0, 2.0, 3.0, 0.4, 5.0, 6.0);
    EXPECT_DOUBLE_EQ(sim.getDarkCurrent(),       1.0);
    EXPECT_DOUBLE_EQ(sim.getReadNoise(),         2.0);
    EXPECT_DOUBLE_EQ(sim.getGain(),              3.0);
    EXPECT_DOUBLE_EQ(sim.getQuantumEfficiency(), 0.4);
    EXPECT_DOUBLE_EQ(sim.getKTCnoise(),          5.0);
    EXPECT_DOUBLE_EQ(sim.getBaseline(),          6.0);
}

TEST(pixsim, InvalidParameters)
{
    EXPECT_THROW(pixsim(-1.0, 0.0, 1.0, 1.0, 0.0, 0.0), std::invalid_argument); // dark < 0
    EXPECT_THROW(pixsim(0.0, -1.0, 1.0, 1.0, 0.0, 0.0), std::invalid_argument); // readNoise < 0
    EXPECT_THROW(pixsim(0.0,  0.0, 0.0, 1.0, 0.0, 0.0), std::invalid_argument); // gain = 0
    EXPECT_THROW(pixsim(0.0,  0.0,-1.0, 1.0, 0.0, 0.0), std::invalid_argument); // gain < 0
    EXPECT_THROW(pixsim(0.0,  0.0, 1.0,-1.0, 0.0, 0.0), std::invalid_argument); // QE < 0
    EXPECT_THROW(pixsim(0.0,  0.0, 1.0, 2.0, 0.0, 0.0), std::invalid_argument); // QE > 1
    EXPECT_THROW(pixsim(0.0,  0.0, 1.0, 1.0,-1.0, 0.0), std::invalid_argument); // KTC < 0
    EXPECT_THROW(pixsim(0.0,  0.0, 1.0, 1.0, 0.0,-1.0), std::invalid_argument); // baseline < 0
}

TEST(pixsim, CopyConstructor)
{
    pixsim sim1(1.0, 2.0, 3.0, 0.4, 5.0, 6.0);
    pixsim sim2(sim1);
    EXPECT_DOUBLE_EQ(sim2.getDarkCurrent(),       sim1.getDarkCurrent());
    EXPECT_DOUBLE_EQ(sim2.getReadNoise(),         sim1.getReadNoise());
    EXPECT_DOUBLE_EQ(sim2.getGain(),              sim1.getGain());
    EXPECT_DOUBLE_EQ(sim2.getKTCnoise(),          sim1.getKTCnoise());
    EXPECT_DOUBLE_EQ(sim2.getQuantumEfficiency(), sim1.getQuantumEfficiency());
    EXPECT_DOUBLE_EQ(sim2.getBaseline(),          sim1.getBaseline());
}

TEST(pixsim, AssignmentOperator)
{
    pixsim sim1(1.0, 2.0, 3.0, 0.4, 5.0, 6.0);
    pixsim sim2;
    sim2 = sim1;
    EXPECT_DOUBLE_EQ(sim2.getDarkCurrent(),       sim1.getDarkCurrent());
    EXPECT_DOUBLE_EQ(sim2.getReadNoise(),         sim1.getReadNoise());
    EXPECT_DOUBLE_EQ(sim2.getGain(),              sim1.getGain());
    EXPECT_DOUBLE_EQ(sim2.getKTCnoise(),          sim1.getKTCnoise());
    EXPECT_DOUBLE_EQ(sim2.getQuantumEfficiency(), sim1.getQuantumEfficiency());
    EXPECT_DOUBLE_EQ(sim2.getBaseline(),          sim1.getBaseline());
}

TEST(pixsim, SelfAssignment)
{
    pixsim sim(1.0, 2.0, 3.0, 0.4, 5.0, 6.0);
    sim = sim;
    EXPECT_DOUBLE_EQ(sim.getDarkCurrent(),       1.0);
    EXPECT_DOUBLE_EQ(sim.getReadNoise(),         2.0);
    EXPECT_DOUBLE_EQ(sim.getGain(),              3.0);
    EXPECT_DOUBLE_EQ(sim.getQuantumEfficiency(), 0.4);
    EXPECT_DOUBLE_EQ(sim.getKTCnoise(),          5.0);
    EXPECT_DOUBLE_EQ(sim.getBaseline(),          6.0);
}

TEST(pixsim, ValidateState)
{
    pixsim sim(1.0, 2.0, 3.0, 0.4, 5.0, 6.0);
    EXPECT_THROW(sim.setDarkCurrent(-1.0),       std::invalid_argument);
    EXPECT_THROW(sim.setReadNoise(-1.0),         std::invalid_argument);
    EXPECT_THROW(sim.setGain(-1.0),              std::invalid_argument);
    EXPECT_THROW(sim.setGain(0.0),               std::invalid_argument);
    EXPECT_THROW(sim.setQuantumEfficiency(-1.0), std::invalid_argument);
    EXPECT_THROW(sim.setQuantumEfficiency(1.1),  std::invalid_argument);
    EXPECT_THROW(sim.setBaseline(-2.0),          std::invalid_argument);
    EXPECT_THROW(sim.setKTCnoise(-10.0),         std::invalid_argument);
}

TEST(pixsim, Modifiers)
{
    pixsim sim = pixsim();
    sim.setDarkCurrent(1.0);
    sim.setReadNoise(2.0);
    sim.setGain(3.0);
    sim.setQuantumEfficiency(0.4);
    sim.setKTCnoise(5.0);
    sim.setBaseline(6.0);
    EXPECT_DOUBLE_EQ(sim.getDarkCurrent(),       1.0);
    EXPECT_DOUBLE_EQ(sim.getReadNoise(),         2.0);
    EXPECT_DOUBLE_EQ(sim.getGain(),              3.0);
    EXPECT_DOUBLE_EQ(sim.getQuantumEfficiency(), 0.4);
    EXPECT_DOUBLE_EQ(sim.getKTCnoise(),          5.0);
    EXPECT_DOUBLE_EQ(sim.getBaseline(),          6.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  PixelSimulator — physics tests
//
//  The model draws:
//    output = N(signal_e/gain, sqrt(signal_e + RN²)/gain) + N(baseline, KTC/gain)
//  giving:
//    E[output]   = signal_e/gain + baseline
//    Var[output] = (signal_e + RN² + KTC²) / gain²
//  where signal_e = (flux*QE + dark)*t
//
//  Statistical assertions use 5-sigma bounds (k=5, n=10 000):
//    P(spurious failure) < 6×10⁻⁷ per assertion.
//    SEM of mean       : sqrt(Var/n)
//    std of sample Var : sqrt(2)*Var/sqrt(n)   [N(μ,σ²) exact]
// ═══════════════════════════════════════════════════════════════════════════════

TEST(PixelSimulator, LinearResponse)
{
    // E[output] = flux * t  (QE=1, gain=1, all noise=0, dark=0, baseline=0).
    // Verified via sample mean — a single draw is insufficient for nonzero flux.
    pixsim sim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0);

    const std::vector<double> fluxes        = {0.0, 10.0, 50.0, 100.0, 200.0};  // ph/s
    const std::vector<double> exposureTimes = {1.0, 10.0, 100.0, 500.0, 1000.0, 3600.0}; // s

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;

    for (const double t : exposureTimes)
    {
        for (const double flux : fluxes)
        {
            const double mu = flux * t;  // expected mean [ADU], Var[output] = mu
            const double n  = static_cast<double>(nSample);

            double sigMean = 0.0;
            for (size_t i = 0; i < nSample; ++i)
                sigMean += sim.pixelresponse(flux, t);
            sigMean /= n;

            // flux=0 → output is deterministic (std=0), EXPECT_NEAR(0,0,0) passes trivially
            EXPECT_NEAR(sigMean, mu, k * std::sqrt(mu / n))
                << "flux=" << flux << " t=" << t;
        }
    }
}

TEST(PixelSimulator, Gain)
{
    // E[output]   = signal_e / gain          (mean scales as 1/gain)
    // Var[output] = signal_e / gain²         (variance scales as 1/gain²)
    // where signal_e = flux * t  (QE=1, dark=0, no noise, baseline=0).
    pixsim sim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0);

    const std::vector<double> gains        = {0.5, 0.75, 1.0, 1.5, 2.0, 4.0};
    const std::vector<double> fluxes       = {10.0, 50.0, 100.0, 200.0};
    const std::vector<double> exposureTimes = {100.0};

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;

    for (const double t : exposureTimes)
    {
        for (const double gain : gains)
        {
            sim.setGain(gain);

            for (const double flux : fluxes)
            {
                const double signal_e = flux * t;
                const double mu_ADU   = signal_e / gain;
                const double var_ADU  = signal_e / (gain * gain);
                const double n        = static_cast<double>(nSample);

                double sigMean = 0.0, sigVar = 0.0;
                for (size_t i = 0; i < nSample; ++i)
                {
                    const double x = sim.pixelresponse(flux, t);
                    sigMean += x;
                    sigVar  += x * x;
                }
                sigMean /= n;
                sigVar   = sigVar / n - sigMean * sigMean;

                EXPECT_NEAR(sigMean, mu_ADU,  k * std::sqrt(var_ADU / n))
                    << "gain=" << gain << " flux=" << flux << " t=" << t;
                EXPECT_NEAR(sigVar,  var_ADU, k * std::sqrt(2.0) * var_ADU / std::sqrt(n))
                    << "gain=" << gain << " flux=" << flux << " t=" << t;
            }
        }
    }
}

TEST(PixelSimulator, PTC)
{
    // Photon Transfer Curve: Var/Mean = 1/gain, independent of signal level.
    // By the delta method (large-signal regime):
    //   std(Var/Mean) ≈ sqrt(2) / (gain * sqrt(n))
    pixsim sim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0);

    const std::vector<double> gains        = {0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0};
    const std::vector<double> fluxes       = {10.0, 50.0, 100.0, 200.0};
    const std::vector<double> exposureTimes = {100.0};

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;

    for (const double t : exposureTimes)
    {
        for (const double gain : gains)
        {
            sim.setGain(gain);

            for (const double flux : fluxes)
            {
                const double n = static_cast<double>(nSample);

                double sigMean = 0.0, sigVar = 0.0;
                for (size_t i = 0; i < nSample; ++i)
                {
                    const double x = sim.pixelresponse(flux, t);
                    sigMean += x;
                    sigVar  += x * x;
                }
                sigMean /= n;
                sigVar   = sigVar / n - sigMean * sigMean;

                EXPECT_NEAR(sigVar / sigMean, 1.0 / gain,
                            k * std::sqrt(2.0) / (gain * std::sqrt(n)))
                    << "gain=" << gain << " flux=" << flux << " t=" << t;
            }
        }
    }
}

TEST(PixelSimulator, QE)
{
    // QE scales the number of detected electrons: signal_e = flux * t * QE.
    //   E[output]   = flux * t * QE    (absolute mean level scales with QE)
    //   Var[output] = flux * t * QE    (Poisson, gain=1)
    // NOTE: Var/Mean = 1/gain = 1 regardless of QE — QE does not change SNR ratio.
    pixsim sim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0);

    const std::vector<double> inputQE      = {0.25, 0.5, 0.75, 1.0};
    const std::vector<double> fluxes       = {10.0, 50.0, 100.0, 200.0};
    const std::vector<double> exposureTimes = {100.0};

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;

    for (const double t : exposureTimes)
    {
        for (const double qe : inputQE)
        {
            sim.setQuantumEfficiency(qe);

            for (const double flux : fluxes)
            {
                const double mu = flux * t * qe;
                const double n  = static_cast<double>(nSample);

                double sigMean = 0.0, sigVar = 0.0;
                for (size_t i = 0; i < nSample; ++i)
                {
                    const double x = sim.pixelresponse(flux, t);
                    sigMean += x;
                    sigVar  += x * x;
                }
                sigMean /= n;
                sigVar   = sigVar / n - sigMean * sigMean;

                EXPECT_NEAR(sigMean, mu, k * std::sqrt(mu / n))
                    << "QE=" << qe << " flux=" << flux << " t=" << t;
                EXPECT_NEAR(sigVar,  mu, k * std::sqrt(2.0) * mu / std::sqrt(n))
                    << "QE=" << qe << " flux=" << flux << " t=" << t;
            }
        }
    }
}

TEST(PixelSimulator, DarkCurent)
{
    // Dark current adds signal_e = dark * t electrons, raising mean and Poisson variance.
    // With flux=0: E[output] = dark*t, Var[output] = dark*t  (gain=1).
    pixsim sim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0);

    const std::vector<double> darkRates    = {1e-3, 1e-2, 1e-1, 1.0, 10.0}; // e/s
    const std::vector<double> exposureTimes = {1., 10., 50., 100.0};

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;

    for (const double t : exposureTimes)
    {
        for (const double dark : darkRates)
        {
            sim.setDarkCurrent(dark);
            const double mu = dark * t;
            const double n  = static_cast<double>(nSample);

            double sigMean = 0.0, sigVar = 0.0;
            for (size_t i = 0; i < nSample; ++i)
            {
                const double x = sim.pixelresponse(0.0, t);
                sigMean += x;
                sigVar  += x * x;
            }
            sigMean /= n;
            sigVar   = sigVar / n - sigMean * sigMean;

            EXPECT_NEAR(sigMean, mu, k * std::sqrt(mu / n))
                << "dark=" << dark << " t=" << t;
            EXPECT_NEAR(sigVar,  mu, k * std::sqrt(2.0) * mu / std::sqrt(n))
                << "dark=" << dark << " t=" << t;
        }
    }
}

TEST(PixelSimulator, ReadNoise)
{
    // Read noise is a constant variance floor: Var[output] += RN²/gain².
    // Independent of flux and exposure time.
    // With flux=0, dark=0, gain=1: output ~ N(0, RN), E=0, Var=RN².
    pixsim sim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0);

    const std::vector<double> readNoises   = {1e-2, 1e-1, 1.0, 5.0, 10.0}; // e rms
    const std::vector<double> exposureTimes = {1., 10., 50., 100.0};

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;

    for (const double t : exposureTimes)
    {
        for (const double rn : readNoises)
        {
            sim.setReadNoise(rn);
            const double var = rn * rn;
            const double n   = static_cast<double>(nSample);

            double sigMean = 0.0, sigVar = 0.0;
            for (size_t i = 0; i < nSample; ++i)
            {
                const double x = sim.pixelresponse(0.0, t);
                sigMean += x;
                sigVar  += x * x;
            }
            sigMean /= n;
            sigVar   = sigVar / n - sigMean * sigMean;

            // Mean stays at 0; read noise is zero-mean
            EXPECT_NEAR(sigMean, 0.0, k * std::sqrt(var / n))
                << "readNoise=" << rn << " t=" << t;
            // Variance is constant (no t dependence)
            EXPECT_NEAR(sigVar, var, k * std::sqrt(2.0) * var / std::sqrt(n))
                << "readNoise=" << rn << " t=" << t;
        }
    }
}

TEST(PixelSimulator, Baseline)
{
    // Baseline is a deterministic DC offset on every output.
    // With all noise sources = 0 and flux=0: output == baseline exactly.
    // (sample_normal_distribution(x, 0) must return x)
    pixsim sim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0);

    const std::vector<double> baselines    = {10.0, 100.0, 500.0, 1000.0}; // ADU
    const std::vector<double> exposureTimes = {1., 10., 50., 100.0};

    for (const double t : exposureTimes)
    {
        for (const double baseline : baselines)
        {
            sim.setBaseline(baseline);
            const double output = sim.pixelresponse(0.0, t);
            EXPECT_NEAR(output, baseline, 1e-9)
                << "baseline=" << baseline << " t=" << t;
        }
    }
}

TEST(PixelSimulator, KTCNoise)
{
    // KTC (kT/C) noise is reset noise on the baseline term: Var += (KTC/gain)².
    // Independent of flux and exposure time.
    // With flux=0, dark=0, baseline=0, gain=1: output ~ N(0, KTC), E=0, Var=KTC².
    pixsim sim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0);

    const std::vector<double> ktcNoises    = {1e-2, 1e-1, 1.0, 3.0, 10.0}; // e rms
    const std::vector<double> exposureTimes = {1., 10., 50., 100.0};

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;

    for (const double t : exposureTimes)
    {
        for (const double ktc : ktcNoises)
        {
            sim.setKTCnoise(ktc);
            const double var = ktc * ktc;
            const double n   = static_cast<double>(nSample);

            double sigMean = 0.0, sigVar = 0.0;
            for (size_t i = 0; i < nSample; ++i)
            {
                const double x = sim.pixelresponse(0.0, t);
                sigMean += x;
                sigVar  += x * x;
            }
            sigMean /= n;
            sigVar   = sigVar / n - sigMean * sigMean;

            EXPECT_NEAR(sigMean, 0.0, k * std::sqrt(var / n))
                << "KTCnoise=" << ktc << " t=" << t;
            EXPECT_NEAR(sigVar, var, k * std::sqrt(2.0) * var / std::sqrt(n))
                << "KTCnoise=" << ktc << " t=" << t;
        }
    }
}

TEST(PixelSimulator, FullNoiseBudget)
{
    // Verifies the complete noise budget with ALL sources active simultaneously.
    //
    //  Setup:  flux=10 ph/s,  dark=1 e/s,  readNoise=5 e,  gain=2 e/ADU,
    //          QE=0.75,  KTCnoise=3 e,  baseline=100 ADU,  t=100 s
    //
    //  signal_e = (10×0.75 + 1)×100                       =  850.0  e
    //  E[output]= signal_e/gain + baseline = 425 + 100     =  525.0  ADU
    //  Var[out] = (signal_e + RN² + KTC²) / gain²
    //           = (850 + 25 + 9) / 4                       =  221.0  ADU²
    //
    constexpr double flux      = 10.0;   // ph/s
    constexpr double dark      =  1.0;   // e/s
    constexpr double readNoise =  5.0;   // e rms
    constexpr double gain      =  2.0;   // e/ADU
    constexpr double QE        =  0.75;
    constexpr double KTC       =  3.0;   // e rms
    constexpr double baseline  = 100.0;  // ADU
    constexpr double t         = 100.0;  // s

    constexpr double signal_e = (flux * QE + dark) * t;                              // 850 e
    constexpr double mu_ADU   = signal_e / gain + baseline;                          // 525 ADU
    constexpr double var_ADU  = (signal_e + readNoise*readNoise + KTC*KTC) / (gain*gain); // 221 ADU²

    constexpr double k       = 5.0;
    constexpr size_t nSample = 10'000;
    const     double n       = static_cast<double>(nSample);

    pixsim sim(dark, readNoise, gain, QE, KTC, baseline);

    double sigMean = 0.0, sigVar = 0.0;
    for (size_t i = 0; i < nSample; ++i)
    {
        const double x = sim.pixelresponse(flux, t);
        sigMean += x;
        sigVar  += x * x;
    }
    sigMean /= n;
    sigVar   = sigVar / n - sigMean * sigMean;

    EXPECT_NEAR(sigMean, mu_ADU,  k * std::sqrt(var_ADU / n));
    EXPECT_NEAR(sigVar,  var_ADU, k * std::sqrt(2.0) * var_ADU / std::sqrt(n));
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Inheritance — design contract of the class hierarchy
// ═══════════════════════════════════════════════════════════════════════════════

TEST(Inheritance, SaturationOverride)
{
    // A derived pixel that clips output to a full-well capacity.
    // Validates virtual dispatch on pixelresponse and access to protected members.
    class SaturatedPixel : public pixsim
    {
    public:
        explicit SaturatedPixel(double fullWell)
            : pixsim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0), m_fullWell(fullWell) {}

        double pixelresponse(const double& flux, const double& t) const override
        {
            return std::min(pixsim::pixelresponse(flux, t), m_fullWell);
        }
    private:
        double m_fullWell;
    };

    constexpr double fullWell = 1000.0; // ADU

    SaturatedPixel pixel(fullWell);

    // Below saturation: output equals flux*t (no noise, deterministic)
    EXPECT_NEAR(pixel.pixelresponse(5.0, 10.0), 50.0, 4*std::sqrt(50.));    // 50 ADU  << 1000

    // Above saturation: output is clipped to full-well
    EXPECT_DOUBLE_EQ(pixel.pixelresponse(100.0, 100.0), fullWell); // 10000 >> 1000

    // Virtual dispatch works through a base-class pointer
    const pixsim* base = &pixel;
    EXPECT_DOUBLE_EQ(base->pixelresponse(100.0, 100.0), fullWell);
}

TEST(Inheritance, StageOverrideUsedByPixelresponse)
{
    // A derived pixel that only replaces the integration stage (here a
    // deterministic charge, no shot noise). The default pixelresponse() still
    // chains it with the base reset level and read stages:
    // pixelresponse = read(resetLevel() + integrate()).
    class NoShotNoisePixel : public pixsim
    {
    public:
        NoShotNoisePixel() : pixsim(0.0, 0.0, 2.0, 1.0, 0.0, 10.0) {} // gain 2, baseline 10

        using pixsim::integrate;   // keep the cfg overload visible
        double integrate(const double& flux, const double& t) const override
        {
            return flux * t * m_quantumEfficiency / m_gain;   // ADU
        }
    };

    NoShotNoisePixel pixel;

    EXPECT_DOUBLE_EQ(pixel.integrate(5.0, 4.0), 10.0);              // 20 e- / 2
    EXPECT_DOUBLE_EQ(pixel.resetLevel(),        10.0);              // baseline, no kTC
    EXPECT_DOUBLE_EQ(pixel.pixelresponse(5.0, 4.0), 20.0);          // read(10 + 10), no read noise

    // Through the base class and the cfg overload as well.
    const pixsim* base = &pixel;
    EXPECT_DOUBLE_EQ(base->pixelresponse(5.0, 4.0), 20.0);
    EXPECT_DOUBLE_EQ(base->pixelresponse(5.0, 4.0, PixelConfig{}), 20.0);
}

TEST(Inheritance, UserValidationOverride)
{
    // A derived pixel that constrains dark current to <= 0.01 e/s ("cold pixel").
    //
    // _validateState() is protected so derived constructors can call it explicitly
    // AFTER base construction, when the vtable is fully set up and _userValidation()
    // dispatches to the derived override.
    class ColdPixel : public pixsim
    {
    public:
        // 1-arg ctor: setDarkCurrent() triggers _validateState() with ColdPixel vtable
        explicit ColdPixel(double dark = 0.0) : pixsim()
        {
            setDarkCurrent(dark);
        }

        // Full ctor: base ctor validates base params (pixsim vtable);
        //            explicit call then validates derived params (ColdPixel vtable).
        ColdPixel(double dark, double readNoise, double gain, double QE, double KTC, double baseline)
            : pixsim(dark, readNoise, gain, QE, KTC, baseline)
        {
            _userValidation(); // fires ColdPixel::_userValidation() with derived vtable
        }

    protected:
        void _userValidation() override
        {
            if (m_darkCurrent > 0.01)
                throw std::invalid_argument("Dark current exceeds cold-pixel limit");
        }
    };

    // 1-arg constructor: derived constraint enforced at construction time
    EXPECT_NO_THROW(ColdPixel(0.001));
    EXPECT_THROW(ColdPixel(1.0), std::invalid_argument);           // dark > 0.01

    // Full constructor: base catches invalid base params (pixsim vtable)
    // Note: lambda must be CALLED with () — just defining it does not execute the body
    EXPECT_THROW([]() { ColdPixel p(0.001, 0.0, 1.0, -1.0, 0.0, 0.0); }(),
                 std::invalid_argument);                            // QE < 0, caught by base

    // Full constructor: derived catches its own constraint (ColdPixel vtable)
    EXPECT_THROW([]() { ColdPixel p(1.0, 0.0, 1.0, 1.0, 0.0, 0.0); }(),
                 std::invalid_argument);                            // dark > 0.01, caught by derived

    // Setters go through the same _validateState() path
    ColdPixel pixel(0.0);
    EXPECT_NO_THROW(pixel.setDarkCurrent(0.01));
    EXPECT_THROW(pixel.setDarkCurrent(1.0), std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — getAsDouble()
// ═══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, GetAsDoubleFromDouble)
{
    PixelConfig cfg{ {"gain", 2.5} };
    EXPECT_DOUBLE_EQ(cfg.getAsDouble("gain"), 2.5);
}

TEST(PixelConfig, GetAsDoubleFromInt)
{
    PixelConfig cfg;
    cfg.set<int64_t>("dark", 3);
    EXPECT_DOUBLE_EQ(cfg.getAsDouble("dark"), 3.0);
}

TEST(PixelConfig, GetAsDoubleFromUInt)
{
    PixelConfig cfg;
    cfg.set<uint64_t>("readNoise", 5u);
    EXPECT_DOUBLE_EQ(cfg.getAsDouble("readNoise"), 5.0);
}

TEST(PixelConfig, GetAsDoubleFromFloat)
{
    PixelConfig cfg;
    cfg.set<float>("QE", 0.8f);
    EXPECT_NEAR(cfg.getAsDouble("QE"), 0.8, 1e-6);
}

TEST(PixelConfig, GetAsDoubleStringThrows)
{
    PixelConfig cfg{ {"path", std::string("/some/file.fits")} };
    EXPECT_THROW(cfg.getAsDouble("path"), std::invalid_argument);
}

TEST(PixelConfig, GetAsDoubleMissingKey)
{
    PixelConfig cfg;
    EXPECT_THROW(cfg.getAsDouble("nonexistent"), std::out_of_range);
}

TEST(PixelConfig, GetAsDoubleDefault)
{
    PixelConfig cfg;
    EXPECT_DOUBLE_EQ(cfg.getAsDouble("missing", 42.0), 42.0); // absent → default
    cfg.set<std::string>("path", std::string("/file.fits"));
    EXPECT_DOUBLE_EQ(cfg.getAsDouble("path", 99.0), 99.0);    // non-numeric → default
}
