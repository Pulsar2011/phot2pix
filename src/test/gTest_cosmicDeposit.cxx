//
// Created by gillard on 2026-11-26.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_cosmicDeposit.cxx
//
// Unit tests of the cosmic-ray source phot2pix::CosmicRay::CosmicDeposit:
// (a) species: energy per nucleon or per particle, Lorentz factor;
// (b) CosmicSpectrum: tables, closed-form flux of a bin, probabilities, sampling;
// (c) CosmicSpectrum::load(): CSV file and FITS table, filter, required columns, errors;
// (d) sensor materials: compound(), HgCdTe(x);
// (e) CosmicDeposit configuration: constructor, load(), configure(), all-or-nothing;
// (f) probabilities: rate, Poisson, truncated share;
// (g) generation: number of hits, entry points, directions, reproducibility;
// (h) footprint: charge in each pixel along the track, lateral spread, detector edges, per-species core and halo,
//     delta-ray halo (share from the velocity, charge kept, symmetry, cut, configuration);
// (i) energy loss: Landau parameters, deposit capped at the kinetic energy;
// (j) outputs: table of hits and header cards;
// (k) stopping power and range: StoppingTable, Bethe against the PDG muon tables (Si, Ge, CdTe), e+ and e-;
// (l) range cut, thick crossings, Bragg profile;
// (m) angular models: geometric factors, truncation, directions, faces, tilted sensor, configuration;
// (n) geometry helpers without a source: incidence(), effectiveArea(), geometricFactor(), etendue().
//
// The spectra of the tests are power laws J = J0 (E/10)^-2 on [10, 100] MeV/nucleon, whose
// integral is 9 J0 m-2 s-1 sr-1. Statistical tests use a KS test at alpha = 1e-6.
//

#include <gtest/gtest.h>
#include <phot2pix/CosmicDeposit.h>
#include <DSTmath/DSTstatistic.h>   // DST::Math::ks_uniform_pvalue
#include <DSTfits/FITSmanager.h>    // DSL::FITStable, DSL::FITScolumn

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <variant>
#include <vector>

namespace
{
    using phot2pix::Deposit;
    using phot2pix::Geometry;
    using phot2pix::PixelConfig;
    using phot2pix::json;
    using phot2pix::CosmicRay::CosmicDeposit;
    using phot2pix::CosmicRay::AngularModel;
    using phot2pix::CosmicRay::CosmicSpectrum;
    using phot2pix::CosmicRay::StoppingTable;
    using phot2pix::CosmicRay::specie;
    using phot2pix::CosmicRay::MeV_per_u;
    using material = phot2pix::CosmicRay::target::material;
    using element  = phot2pix::CosmicRay::target::element;

    const specie& proton   = phot2pix::CosmicRay::proton;
    const specie& helium   = phot2pix::CosmicRay::helium;
    const specie& electron = phot2pix::CosmicRay::electron;
    const specie& positron = phot2pix::CosmicRay::positron;

    constexpr double kAlpha = 1e-6;   // KS significance

    // ── Spectra and parameters ───────────────────────────────────────────────────

    // Protons, J = J0 (E/10)^-2 on [10, 100] MeV/nucleon: ∫ J dE = 9 J0.
    CosmicSpectrum protonSpectrum(double J0 = 1.0)
    {
        CosmicSpectrum s;
        s.add(proton, {10.0, 100.0}, {J0, 0.01 * J0});
        return s;
    }

    // 10 um of silicon, 18 um pixels (the material is silicon by default).
    CosmicDeposit::Params siliconParams()
    {
        CosmicDeposit::Params p;
        p.depth_um      = 10.0;
        p.pitch_um      = 18.0;
        p.eVPerElectron = 3.6;
        return p;
    }

    // Duration of a step with `hits` drawn hits on average, on geometry g.
    double stepFor(const CosmicDeposit& d, const Geometry& g, double hits)
    {
        const double area_m2 = static_cast<double>(g.width * g.height) * d.params().pitch_um * d.params().pitch_um * 1.0e-12;
        return hits / (d.eventRate(area_m2) * (1.0 - d.truncatedFraction()));
    }

    // Hit built by hand, for the footprint tests.
    CosmicDeposit::Event hit(double x, double y, double dx, double dy, double electrons)
    {
        CosmicDeposit::Event e{};
        e.species   = 0;
        e.energy    = 100.0;
        e.mu        = 1.0;
        e.phi       = 0.0;
        e.x         = x;
        e.y         = y;
        e.dx        = dx;
        e.dy        = dy;
        e.electrons = electrons;
        return e;
    }

    double total(const std::vector<Deposit>& out)
    {
        double sum = 0.0;
        for (const Deposit& d : out)
            sum += d.electrons;
        return sum;
    }

    // Lateral spread of a species: core sigma (NaN: that of the source), halo sigma and share.
    CosmicSpectrum::Footprint footprint(double core, double haloSigma = 0.0, double haloFraction = 0.0)
    {
        CosmicSpectrum::Footprint f;
        f.sigma_px     = core;
        f.haloSigma_px = haloSigma;
        f.haloFraction = haloFraction;
        return f;
    }

    // Lateral spread of the source, with the delta-ray halo.
    CosmicSpectrum::Footprint withDeltaRays()
    {
        CosmicSpectrum::Footprint f;
        f.deltaRays = true;
        return f;
    }

    // Charge of each pixel of the deposits.
    std::map<size_t, double> perPixel(const std::vector<Deposit>& out)
    {
        std::map<size_t, double> m;
        for (const Deposit& d : out)
            m[d.idx] += d.electrons;
        return m;
    }

    // ── Files ────────────────────────────────────────────────────────────────────

    // File in the temporary directory, removed at the end of the test.
    struct TempFile
    {
        std::filesystem::path path;
        explicit TempFile(const std::string& name) : path(std::filesystem::temp_directory_path() / name) {}
        ~TempFile() { std::error_code ec; std::filesystem::remove(path, ec); }
        std::string str() const { return path.string(); }
    };

    // Protons and helium in one CSV table, in the layout of DSL::FITStable::toCSV() (units after the
    // names), with a comment and a blank line.
    void writeSpectrumCsv(const TempFile& f)
    {
        std::ofstream(f.path) << "# protons and helium, J = J0 (E/10)^-2\n"
                                 "NAME,ENERGY [MeV/nucleon],FLUX [m-2 s-1 sr-1 (MeV/nucleon)-1]\n"
                                 "p,10,1\n"
                                 "p,100,0.01\n"
                                 "\n"
                                 "He,10,0.5\n"
                                 "He,100,0.005\n";
    }

    // The same table as a FITS binary table named GCR, in HDU 2.
    void writeSpectrumFits(const TempFile& f)
    {
        const std::vector<std::string> names{"p", "p", "He", "He"};
        const std::vector<double>      E{10.0, 100.0, 10.0, 100.0};
        const std::vector<double>      J{1.0, 0.01, 0.5, 0.005};

        auto name   = std::make_shared<DSL::FITScolumn<std::string>>("NAME",   DSL::tstring, "",                             1);
        auto energy = std::make_shared<DSL::FITScolumn<double>>     ("ENERGY", DSL::tdouble, "MeV/nucleon",                  2);
        auto flux   = std::make_shared<DSL::FITScolumn<double>>     ("FLUX",   DSL::tdouble, "m-2 s-1 sr-1 (MeV/nucleon)-1", 3);
        for (size_t i = 0; i < names.size(); ++i)
        {
            name->push_back(names[i]);
            energy->push_back(E[i]);
            flux->push_back(J[i]);
        }

        DSL::FITStable table;
        table.InsertColumn(name);
        table.InsertColumn(energy);
        table.InsertColumn(flux);
        table.setName("GCR");
        table.write(f.str(), 0, true);
    }

    // Spectrum entry of a species, with the column names.
    json entry(const std::string& species, double Z, double A, const std::string& file)
    {
        return json{{"species", species}, {"Z", Z}, {"A", A}, {"file", file}, {"colEnergy", "ENERGY"}, {"colFlux", "FLUX"}};
    }

    json filtered(json e, const std::string& value)
    {
        e["filterColumn"] = "NAME";
        e["filterValue"]  = value;
        return e;
    }

    PixelConfig spectrumConfig(const std::vector<json>& entries)
    {
        json cfg;
        cfg["spectrum"] = json::array();
        for (const json& e : entries)
            cfg["spectrum"].push_back(e);
        return PixelConfig(cfg);
    }
}

// ── (a) Species ──────────────────────────────────────────────────────────────────

TEST(CosmicSpecieTest, EnergyPerNucleonOrPerParticle)
{
    EXPECT_TRUE(proton.perNucleon());
    EXPECT_TRUE(helium.perNucleon());
    EXPECT_FALSE(electron.perNucleon());

    EXPECT_DOUBLE_EQ(helium.kineticEnergy(100.0), 400.0);   // 100 MeV/nucleon x 4 nucleons
    EXPECT_DOUBLE_EQ(electron.kineticEnergy(1.0), 1.0);      // 1 MeV per electron
}

TEST(CosmicSpecieTest, LorentzFactor)
{
    EXPECT_NEAR(proton.gamma(1000.0), 1.0 + 1000.0 / MeV_per_u, 1e-12);
    EXPECT_NEAR(proton.gamma(1000.0), 2.0735, 1e-4);
    EXPECT_NEAR(helium.gamma(1000.0), proton.gamma(1000.0), 1e-12);   // same speed per nucleon
    EXPECT_NEAR(electron.gamma(1.0), 2.957, 1e-3);                    // 1 + 1 MeV / 0.511 MeV
}

TEST(CosmicSpecieTest, RejectsAMassNumberNotPositive)
{
    EXPECT_THROW(specie("X", 1.0, 0.0), std::invalid_argument);
    EXPECT_THROW(specie("X", 1.0, -1.0), std::invalid_argument);
}

// ── (b) CosmicSpectrum ───────────────────────────────────────────────────────────

TEST(CosmicSpectrumTest, EmptyByDefault)
{
    const CosmicSpectrum s;
    EXPECT_TRUE(s.empty());
    EXPECT_EQ(s.integralIntensity(), 0.0);
}

TEST(CosmicSpectrumTest, FluxOfABinIsTheIntegralOfThePowerLaw)
{
    // g = 2: ∫_10^100 (E/10)^-2 dE = 10 (1 - 1/10) = 9
    EXPECT_NEAR(protonSpectrum().integralIntensity(), 9.0, 1e-12);

    // g = 1: ∫_10^100 (E/10)^-1 dE = 10 ln(10)
    CosmicSpectrum s1;
    s1.add(proton, {10.0, 100.0}, {1.0, 0.1});
    EXPECT_NEAR(s1.integralIntensity(), 10.0 * std::log(10.0), 1e-9);

    // g = 0: flat, J0 (E1 - E0)
    CosmicSpectrum s0;
    s0.add(proton, {10.0, 100.0}, {2.0, 2.0});
    EXPECT_NEAR(s0.integralIntensity(), 180.0, 1e-9);
}

TEST(CosmicSpectrumTest, AddRejectsBadTablesAndKeepsTheSpectrum)
{
    CosmicSpectrum s = protonSpectrum();
    EXPECT_THROW(s.add(helium, {10.0}, {1.0}), std::invalid_argument);                  // one point
    EXPECT_THROW(s.add(helium, {10.0, 100.0}, {1.0}), std::invalid_argument);           // sizes differ
    EXPECT_THROW(s.add(helium, {100.0, 10.0}, {1.0, 0.1}), std::invalid_argument);      // not ascending
    EXPECT_THROW(s.add(helium, {10.0, 10.0}, {1.0, 0.1}), std::invalid_argument);       // not strictly
    EXPECT_THROW(s.add(helium, {10.0, 100.0}, {1.0, 0.0}), std::invalid_argument);      // J = 0
    EXPECT_THROW(s.add(helium, {0.0, 100.0}, {1.0, 0.1}), std::invalid_argument);       // E = 0

    EXPECT_EQ(s.species().size(), 1u);
    EXPECT_EQ(s.bins().size(), 1u);
    EXPECT_NEAR(s.integralIntensity(), 9.0, 1e-12);
}

TEST(CosmicSpectrumTest, FootprintOfEachSpecies)
{
    CosmicSpectrum s = protonSpectrum();   // default footprint: the source's sigma, no halo
    s.add(helium, {10.0, 100.0}, {1.0, 0.01}, AngularModel(), footprint(0.4, 3.0, 0.2));

    ASSERT_EQ(s.footprints().size(), 2u);
    EXPECT_TRUE(std::isnan(s.footprints()[0].sigma_px));
    EXPECT_DOUBLE_EQ(s.footprints()[0].haloFraction, 0.0);
    EXPECT_DOUBLE_EQ(s.footprints()[1].sigma_px, 0.4);
    EXPECT_DOUBLE_EQ(s.footprints()[1].haloSigma_px, 3.0);
    EXPECT_DOUBLE_EQ(s.footprints()[1].haloFraction, 0.2);

    const std::vector<double> E{10.0, 100.0}, J{1.0, 0.01};
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(s.add(helium, E, J, AngularModel(), footprint(-1.0)), std::invalid_argument);            // core sigma < 0
    EXPECT_THROW(s.add(helium, E, J, AngularModel(), footprint(inf)), std::invalid_argument);             // not finite
    EXPECT_THROW(s.add(helium, E, J, AngularModel(), footprint(0.5, 2.0, 1.5)), std::invalid_argument);   // share > 1
    EXPECT_THROW(s.add(helium, E, J, AngularModel(), footprint(0.5, 2.0, -0.1)), std::invalid_argument);  // share < 0
    EXPECT_THROW(s.add(helium, E, J, AngularModel(), footprint(0.5, 0.0, 0.2)), std::invalid_argument);   // halo without sigma
    EXPECT_THROW(s.add(helium, E, J, AngularModel(), footprint(0.5, -1.0)), std::invalid_argument);       // halo sigma < 0
    EXPECT_EQ(s.species().size(), 2u);
    EXPECT_EQ(s.footprints().size(), 2u);
}

TEST(CosmicSpectrumTest, ProbabilitiesSumToOne)
{
    CosmicSpectrum s = protonSpectrum();                                  // 9
    s.add(helium, {10.0, 100.0, 1000.0}, {0.5, 0.005, 0.00005});         // two bins: 4.5 + 0.45

    ASSERT_EQ(s.bins().size(), 3u);
    EXPECT_NEAR(s.integralIntensity(1), 4.95, 1e-12);

    const std::vector<double> pb = s.binProbabilities();
    ASSERT_EQ(pb.size(), 3u);
    EXPECT_NEAR(std::accumulate(pb.begin(), pb.end(), 0.0), 1.0, 1e-12);

    const std::vector<double> ps = s.speciesProbabilities();
    ASSERT_EQ(ps.size(), 2u);
    EXPECT_NEAR(ps[0], 9.0 / 13.95, 1e-12);
    EXPECT_NEAR(ps[0] + ps[1], 1.0, 1e-12);
}

TEST(CosmicSpectrumTest, SampleFollowsThePowerLaw)
{
    const CosmicSpectrum s = protonSpectrum();   // g = 2 on [10, 100]
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    const size_t N = 20000;
    std::vector<double> u(N);
    for (size_t i = 0; i < N; ++i)
    {
        const double uBin = U(rng), uEnergy = U(rng);
        const CosmicSpectrum::Sample smp = s.sample(uBin, uEnergy);
        ASSERT_EQ(smp.species, 0u);
        ASSERT_GE(smp.E, 10.0);
        ASSERT_LE(smp.E, 100.0);
        u[i] = (1.0 - 10.0 / smp.E) / 0.9;   // CDF of E^-2 on [10, 100]
    }
    EXPECT_GT(DST::Math::ks_uniform_pvalue(u, 0.0, 1.0), kAlpha);

    EXPECT_NEAR(s.sample(0.5, 0.0).E, 10.0, 1e-9);
    EXPECT_NEAR(s.sample(0.5, 1.0).E, 100.0, 1e-9);
}

TEST(CosmicSpectrumTest, SampleChoosesTheSpeciesByItsFlux)
{
    CosmicSpectrum s = protonSpectrum();             // 9
    s.add(helium, {10.0, 100.0}, {1.0, 0.01});      // 9
    EXPECT_EQ(s.sample(0.25, 0.5).species, 0u);
    EXPECT_EQ(s.sample(0.75, 0.5).species, 1u);
}

// ── (c) CosmicSpectrum::load() ───────────────────────────────────────────────────

TEST(CosmicSpectrumLoadTest, CsvFileFilteredBySpecies)
{
    const TempFile csv("p2p_gTest_cosmic_load_csv.csv");
    writeSpectrumCsv(csv);

    const CosmicSpectrum s = CosmicSpectrum::load(spectrumConfig({filtered(entry("p", 1, 1, csv.str()), "p"),
                                                                  filtered(entry("He", 2, 4, csv.str()), "He")}));
    ASSERT_EQ(s.species().size(), 2u);
    EXPECT_EQ(s.species()[0].name(), "p");
    EXPECT_EQ(s.species()[1].name(), "He");
    EXPECT_DOUBLE_EQ(s.species()[1].Z(), 2.0);
    EXPECT_DOUBLE_EQ(s.species()[1].A(), 4.0);
    EXPECT_NEAR(s.integralIntensity(0), 9.0, 1e-12);
    EXPECT_NEAR(s.integralIntensity(1), 4.5, 1e-12);
}

TEST(CosmicSpectrumLoadTest, ParticlesFromTheirRestEnergyHaveAnEnergyPerParticle)
{
    const TempFile csv("p2p_gTest_cosmic_load_leptons.csv");
    std::ofstream(csv.path) << "NAME,E,J\n"
                               "e-,1,100\n"
                               "e-,10,1\n"
                               "mu-,100,10\n"
                               "mu-,1000,0.1\n";

    const auto lepton = [&csv](const std::string& name, double mc2)
    {
        return json{{"species", name}, {"Z", -1.0}, {"mc2", mc2}, {"file", csv.str()},
                    {"filterColumn", "NAME"}, {"filterValue", name}, {"colEnergy", "E"}, {"colFlux", "J"}};
    };
    const CosmicSpectrum s = CosmicSpectrum::load(spectrumConfig({lepton("e-", 0.51099895), lepton("mu-", 105.6584)}));

    ASSERT_EQ(s.species().size(), 2u);
    const specie& e  = s.species()[0];
    const specie& mu = s.species()[1];
    EXPECT_FALSE(e.perNucleon());
    EXPECT_FALSE(mu.perNucleon());
    EXPECT_DOUBLE_EQ(e.Z(), -1.0);
    EXPECT_NEAR(e.mc2(), 0.51099895, 1e-12);
    EXPECT_NEAR(mu.mc2(), 105.6584, 1e-9);
    EXPECT_DOUBLE_EQ(e.kineticEnergy(10.0), 10.0);              // per particle
    EXPECT_NEAR(e.gamma(1.0), 1.0 + 1.0 / 0.51099895, 1e-12);
    EXPECT_NEAR(mu.gamma(1000.0), 1.0 + 1000.0 / 105.6584, 1e-12);

    // g = 2 on a decade: ∫ J dE = J0 E0 (1 - 1/10)
    EXPECT_NEAR(s.integralIntensity(0), 100.0 * 1.0 * 0.9, 1e-9);
    EXPECT_NEAR(s.integralIntensity(1), 10.0 * 100.0 * 0.9, 1e-9);
}

TEST(CosmicSpectrumLoadTest, MassNumberOrRestEnergyNotBoth)
{
    const TempFile csv("p2p_gTest_cosmic_load_mass.csv");
    writeSpectrumCsv(csv);
    const auto load = [](const json& e) { return CosmicSpectrum::load(spectrumConfig({e})); };

    json both = filtered(entry("p", 1, 1, csv.str()), "p");
    both["mc2"] = 938.272;
    EXPECT_THROW(load(both), std::invalid_argument);

    json neither = filtered(entry("p", 1, 1, csv.str()), "p");
    neither.erase("A");
    EXPECT_THROW(load(neither), std::invalid_argument);

    json negative = neither;
    negative["mc2"] = -1.0;
    EXPECT_THROW(load(negative), std::invalid_argument);   // A = mc2 / u <= 0

    json restEnergy = neither;
    restEnergy["mc2"] = 938.272;                            // a proton given by its rest energy: per particle
    const CosmicSpectrum s = load(restEnergy);
    EXPECT_FALSE(s.species()[0].perNucleon());
}

TEST(CosmicSpectrumLoadTest, FitsTableByIndexOrName)
{
    const TempFile fits("p2p_gTest_cosmic_load_fits.fits");
    writeSpectrumFits(fits);

    json p = filtered(entry("p", 1, 1, fits.str()), "p");
    p["hdu"] = 2;        // first extension
    json he = filtered(entry("He", 2, 4, fits.str()), "He");
    he["hdu"] = "GCR";   // EXTNAME
    const CosmicSpectrum s = CosmicSpectrum::load(spectrumConfig({p, he}));
    ASSERT_EQ(s.species().size(), 2u);
    EXPECT_NEAR(s.integralIntensity(0), 9.0, 1e-12);
    EXPECT_NEAR(s.integralIntensity(1), 4.5, 1e-12);

    json noHdu = filtered(entry("p", 1, 1, fits.str()), "p");
    EXPECT_THROW(CosmicSpectrum::load(spectrumConfig({noHdu})), std::invalid_argument);     // "hdu" is required

    json primary = noHdu;
    primary["hdu"] = 1;   // the primary HDU holds no table
    EXPECT_THROW(CosmicSpectrum::load(spectrumConfig({primary})), std::invalid_argument);
}

TEST(CosmicSpectrumLoadTest, ErrorsOfAnEntry)
{
    const TempFile csv("p2p_gTest_cosmic_load_errors.csv");
    writeSpectrumCsv(csv);
    const auto load = [](const json& e) { return CosmicSpectrum::load(spectrumConfig({e})); };

    // Without a filter, the rows of both species make E non ascending
    EXPECT_THROW(load(entry("p", 1, 1, csv.str())), std::invalid_argument);

    json noEnergy = filtered(entry("p", 1, 1, csv.str()), "p");
    noEnergy.erase("colEnergy");
    EXPECT_THROW(load(noEnergy), std::invalid_argument);   // colEnergy is required

    json noFlux = filtered(entry("p", 1, 1, csv.str()), "p");
    noFlux.erase("colFlux");
    EXPECT_THROW(load(noFlux), std::invalid_argument);     // colFlux is required

    json badColumn = filtered(entry("p", 1, 1, csv.str()), "p");
    badColumn["colEnergy"] = "E_KIN";
    EXPECT_THROW(load(badColumn), std::invalid_argument);  // no such column

    EXPECT_THROW(load(filtered(entry("Li", 3, 7, csv.str()), "Li")), std::invalid_argument);   // no row

    json halfFilter = entry("p", 1, 1, csv.str());
    halfFilter["filterColumn"] = "NAME";
    EXPECT_THROW(load(halfFilter), std::invalid_argument); // filterColumn without filterValue

    json hduCsv = filtered(entry("p", 1, 1, csv.str()), "p");
    hduCsv["hdu"] = 2;
    EXPECT_THROW(load(hduCsv), std::invalid_argument);     // "hdu" goes with a FITS file

    EXPECT_THROW(load(filtered(entry("p", 1, 1, "/no/such/spectrum.csv"), "p")), std::invalid_argument);

    EXPECT_THROW(CosmicSpectrum::load(PixelConfig(json{{"spectrum", json::array()}})), std::invalid_argument);
}

TEST(CosmicSpectrumLoadTest, FootprintOfAnEntry)
{
    const TempFile csv("p2p_gTest_cosmic_load_footprint.csv");
    writeSpectrumCsv(csv);

    json p = filtered(entry("p", 1, 1, csv.str()), "p");
    p["sigmaDiffusion_px"] = 0.3;
    p["halo"] = {{"sigma_px", 2.5}, {"fraction", 0.1}};
    const json he = filtered(entry("He", 2, 4, csv.str()), "He");   // nothing given

    const CosmicSpectrum s = CosmicSpectrum::load(spectrumConfig({p, he}));
    ASSERT_EQ(s.footprints().size(), 2u);
    EXPECT_DOUBLE_EQ(s.footprints()[0].sigma_px, 0.3);
    EXPECT_DOUBLE_EQ(s.footprints()[0].haloSigma_px, 2.5);
    EXPECT_DOUBLE_EQ(s.footprints()[0].haloFraction, 0.1);
    EXPECT_TRUE(std::isnan(s.footprints()[1].sigma_px));
    EXPECT_DOUBLE_EQ(s.footprints()[1].haloFraction, 0.0);

    const auto load = [](const json& e) { return CosmicSpectrum::load(spectrumConfig({e})); };
    const auto with = [&he](const char* key, const json& value) { json e = he; e[key] = value; return e; };

    EXPECT_THROW(load(with("sigmaDiffusion_px", -1.0)), std::invalid_argument);
    EXPECT_THROW(load(with("sigmaDiffusion_px", "wide")), std::invalid_argument);
    EXPECT_THROW(load(with("halo", json::array({2.0, 0.1}))), std::invalid_argument);                  // not an object
    EXPECT_THROW(load(with("halo", json{{"sigma_px", 2.0}})), std::invalid_argument);                  // no fraction
    EXPECT_THROW(load(with("halo", json{{"fraction", 0.1}})), std::invalid_argument);                  // no sigma_px
    EXPECT_THROW(load(with("halo", json{{"sigma_px", 2.0}, {"fraction", 1.5}})), std::invalid_argument);
    EXPECT_THROW(load(with("halo", json{{"sigma_px", 0.0}, {"fraction", 0.1}})), std::invalid_argument);
}

// ── (d) Sensor materials ─────────────────────────────────────────────────────────

TEST(CosmicMaterialTest, CompoundOfOneElementIsTheElement)
{
    // Silicon: Z/A, I and density of PDG; Cbar = 2 ln(I / hbar omega_p) + 1 = 4.4355 (PDG)
    const material si = material::compound({{14.0, 28.0855, 173.0, 1.0}}, 2.329);
    EXPECT_NEAR(si.ZoverA, 0.49848, 1e-5);
    EXPECT_NEAR(si.meanExcitation_eV, 173.0, 1e-9);
    EXPECT_NEAR(si.sternheimerC, 4.4355, 1e-3);

    // Sternheimer & Peierls, solid with I >= 100 eV and Cbar < 5.215
    EXPECT_DOUBLE_EQ(si.sternheimerX0, 0.2);
    EXPECT_DOUBLE_EQ(si.sternheimerX1, 3.0);
    EXPECT_DOUBLE_EQ(si.sternheimerM, 3.0);
    EXPECT_DOUBLE_EQ(si.sternheimerDelta0, 0.0);
    EXPECT_TRUE(std::isnan(si.sternheimerDeltaMax));
}

TEST(CosmicMaterialTest, HgCdTeOfCadmiumFractionOneIsCdTe)
{
    // PDG 2024, cadmium telluride
    const material cdte = material::HgCdTe(1.0);
    EXPECT_NEAR(cdte.ZoverA,            0.41665, 1e-3 * 0.41665);
    EXPECT_NEAR(cdte.meanExcitation_eV, 539.3,   1e-3 * 539.3);
    EXPECT_NEAR(cdte.density_gcm3,      5.850,   2e-3 * 5.850);   // lattice constant of Higgins et al.: 5.857
    EXPECT_NEAR(cdte.sternheimerC,      5.9677,  1e-3 * 5.9677);
}

TEST(CosmicMaterialTest, HgCdTeDensityBetweenHgTeAndCdTe)
{
    const double rhoHgTe = material::HgCdTe(0.0).density_gcm3;
    const double rhoMid  = material::HgCdTe(0.45).density_gcm3;
    const double rhoCdTe = material::HgCdTe(1.0).density_gcm3;
    EXPECT_NEAR(rhoHgTe, 8.08, 0.01);
    EXPECT_GT(rhoMid, rhoCdTe);
    EXPECT_LT(rhoMid, rhoHgTe);
}

TEST(CosmicMaterialTest, RejectsBadInput)
{
    EXPECT_THROW(material::HgCdTe(-0.1), std::invalid_argument);
    EXPECT_THROW(material::HgCdTe(1.1), std::invalid_argument);
    EXPECT_THROW(material::compound({}, 2.0), std::invalid_argument);
    EXPECT_THROW(material::compound({{14.0, 28.0855, 173.0, 1.0}}, 0.0), std::invalid_argument);
    EXPECT_THROW(material::compound({{0.0, 28.0855, 173.0, 1.0}}, 2.0), std::invalid_argument);
    EXPECT_THROW(material::compound({{14.0, 28.0855, 173.0, 0.0}}, 2.0), std::invalid_argument);   // no atom
}

// ── (e) Configuration ────────────────────────────────────────────────────────────

TEST(CosmicDepositConfigTest, NotConfiguredByDefault)
{
    const CosmicDeposit d;
    EXPECT_FALSE(d.isConfigured());

    std::vector<Deposit> out;
    EXPECT_THROW(d.generate(1.0, Geometry(8, 8), out), std::logic_error);
}

TEST(CosmicDepositConfigTest, ConstructorChecksItsArguments)
{
    EXPECT_TRUE(CosmicDeposit(protonSpectrum(), siliconParams()).isConfigured());

    CosmicDeposit::Params noDepth = siliconParams();
    noDepth.depth_um = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW((void)CosmicDeposit(protonSpectrum(), noDepth), std::invalid_argument);

    CosmicDeposit::Params negativeSpread = siliconParams();
    negativeSpread.sigmaDiffusion_px = -1.0;
    EXPECT_THROW((void)CosmicDeposit(protonSpectrum(), negativeSpread), std::invalid_argument);

    EXPECT_THROW((void)CosmicDeposit(CosmicSpectrum(), siliconParams()), std::invalid_argument);
}

TEST(CosmicDepositConfigTest, LoadReadsEveryKey)
{
    const TempFile csv("p2p_gTest_cosmic_config.csv");
    writeSpectrumCsv(csv);

    json cosmic = {{"depth_um", 10.0}, {"pitch_um", 18.0}, {"eVPerElectron", 3.6}, {"tilt_deg", 30.0},
                   {"sigmaDiffusion_px", 0.5}, {"maxTrackPixels", 32.0}, {"seed", 17}};
    cosmic["spectrum"] = json::array({filtered(entry("p", 1, 1, csv.str()), "p")});

    // configure() reads the "cosmic" key of a detector configuration
    CosmicDeposit d;
    d.configure(PixelConfig(json{{"cosmic", cosmic}}));
    ASSERT_TRUE(d.isConfigured());
    EXPECT_DOUBLE_EQ(d.params().depth_um, 10.0);
    EXPECT_DOUBLE_EQ(d.params().pitch_um, 18.0);
    EXPECT_DOUBLE_EQ(d.params().eVPerElectron, 3.6);
    EXPECT_DOUBLE_EQ(d.params().tilt_deg, 30.0);
    EXPECT_DOUBLE_EQ(d.params().sigmaDiffusion_px, 0.5);
    EXPECT_DOUBLE_EQ(d.params().maxTrackPixels, 32.0);
    EXPECT_EQ(d.seed(), 17u);
    EXPECT_NEAR(d.spectrum().integralIntensity(), 9.0, 1e-12);

    // and nothing without it
    CosmicDeposit untouched;
    untouched.configure(PixelConfig(json{{"gain", 1.5}}));
    EXPECT_FALSE(untouched.isConfigured());
}

TEST(CosmicDepositConfigTest, ABadLoadLeavesTheSourceUnchanged)
{
    // A detector configuration whose "cosmic" object is j
    const auto cosmic = [](const json& j) { return PixelConfig(json{{"cosmic", j}}); };

    CosmicDeposit d(protonSpectrum(), siliconParams());
    d.setSeed(3);

    EXPECT_THROW(d.configure(cosmic(json{{"depth_um", "ten"}})), std::invalid_argument);                   // not a number
    EXPECT_THROW(d.configure(cosmic(json{{"depth_um", -1.0}})), std::invalid_argument);                    // not > 0
    EXPECT_THROW(d.configure(cosmic(json{{"seed", -1}})), std::invalid_argument);                          // not >= 0
    EXPECT_THROW(d.configure(cosmic(json{{"tilt_deg", 200.0}})), std::invalid_argument);                   // not in [0, 180]
    EXPECT_THROW(d.configure(cosmic(json{{"acceptance_sr", 3.0}})), std::invalid_argument);                // no longer read
    EXPECT_THROW(d.configure(cosmic(json{{"depth_um", 20.0}, {"seed", -1}})), std::invalid_argument);      // the good value isn't kept

    EXPECT_DOUBLE_EQ(d.params().depth_um, 10.0);
    EXPECT_EQ(d.seed(), 3u);
    EXPECT_TRUE(d.isConfigured());

    CosmicDeposit empty;
    EXPECT_THROW(empty.configure(cosmic(json{{"depth_um", 10.0}, {"pitch_um", 18.0}})), std::invalid_argument);   // no spectrum
    EXPECT_FALSE(empty.isConfigured());
    EXPECT_THROW(empty.configure(PixelConfig(json{{"cosmic", json::array()}})), std::invalid_argument);           // not an object
}

// ── (f) Probabilities ────────────────────────────────────────────────────────────

TEST(CosmicDepositProbabilityTest, RateAndPoisson)
{
    const CosmicDeposit d(protonSpectrum(), siliconParams());   // ∫ J dE = 9, acceptance pi
    const double area = 1.0e-4;   // m2
    EXPECT_NEAR(d.eventRate(area), std::numbers::pi * 9.0 * area, 1e-15);

    const double dt = 1000.0;
    EXPECT_NEAR(d.expectedEvents(area, dt), d.eventRate(area) * dt, 1e-12);

    double sum = 0.0;
    for (size_t n = 0; n < 200; ++n)
        sum += d.probability(n, area, dt);
    EXPECT_NEAR(sum, 1.0, 1e-12);
    EXPECT_NEAR(d.probabilityAtLeastOne(area, dt), 1.0 - d.probability(0, area, dt), 1e-12);
}

TEST(CosmicDepositProbabilityTest, TruncatedFraction)
{
    const CosmicDeposit d(protonSpectrum(), siliconParams());
    const CosmicDeposit::Params& p = d.params();
    const double q = p.maxTrackPixels * p.pitch_um / p.depth_um;
    EXPECT_DOUBLE_EQ(d.truncatedFraction(), 1.0 / (1.0 + q * q));
}

// ── (g) Generation ───────────────────────────────────────────────────────────────

TEST(CosmicDepositGenerationTest, NumberEntryPointsAndDirections)
{
    CosmicDeposit::Params p = siliconParams();
    p.maxTrackPixels = 8.0;
    const CosmicDeposit d(protonSpectrum(1.0e6), p);
    const Geometry g(256, 256);
    const double dt = stepFor(d, g, 50.0);
    const double f  = d.truncatedFraction();

    std::mt19937_64 rng(2026);
    std::vector<double> mu2, phi;
    double count = 0.0;
    const size_t steps = 400;
    for (size_t s = 0; s < steps; ++s)
    {
        const std::vector<CosmicDeposit::Event> events = d.drawEvents(dt, g, rng);
        count += static_cast<double>(events.size());
        for (const CosmicDeposit::Event& e : events)
        {
            ASSERT_GE(e.x, 0.0);
            ASSERT_LT(e.x, 256.0);
            ASSERT_GE(e.y, 0.0);
            ASSERT_LT(e.y, 256.0);
            mu2.push_back((e.mu * e.mu - f) / (1.0 - f));
            phi.push_back(e.phi / (2.0 * std::numbers::pi));

            // Projected track: (d / pitch) tan(theta) pixels
            ASSERT_NEAR(std::hypot(e.dx, e.dy), p.depth_um / p.pitch_um * std::sqrt(1.0 - e.mu * e.mu) / e.mu, 1e-9);
        }
    }

    // Poisson, 50 per step: standard error of the mean sqrt(50 / steps)
    EXPECT_NEAR(count / static_cast<double>(steps), 50.0, 4.0 * std::sqrt(50.0 / static_cast<double>(steps)));
    // mu^2 uniform above the truncated share, phi uniform
    EXPECT_GT(DST::Math::ks_uniform_pvalue(mu2, 0.0, 1.0), kAlpha);
    EXPECT_GT(DST::Math::ks_uniform_pvalue(phi, 0.0, 1.0), kAlpha);
}

TEST(CosmicDepositGenerationTest, NoStepNoHit)
{
    const CosmicDeposit d(protonSpectrum(1.0e6), siliconParams());
    std::mt19937_64 rng(1);
    EXPECT_TRUE(d.drawEvents(0.0, Geometry(64, 64), rng).empty());
}

TEST(CosmicDepositGenerationTest, SameSeedSameHitsAndDeposits)
{
    CosmicDeposit a(protonSpectrum(1.0e6), siliconParams());
    CosmicDeposit b(protonSpectrum(1.0e6), siliconParams());
    a.setSeed(77);
    b.setSeed(77);
    const Geometry g(128, 128);
    const double dt = stepFor(a, g, 20.0);

    std::vector<Deposit> outA, outB;
    std::vector<CosmicDeposit::Event> evA, evB;
    a.generate(dt, g, outA, &evA);
    b.generate(dt, g, outB, &evB);

    ASSERT_FALSE(evA.empty());
    ASSERT_EQ(evA.size(), evB.size());
    for (size_t i = 0; i < evA.size(); ++i)
    {
        EXPECT_EQ(evA[i].species, evB[i].species);
        EXPECT_EQ(evA[i].energy, evB[i].energy);
        EXPECT_EQ(evA[i].mu, evB[i].mu);
        EXPECT_EQ(evA[i].phi, evB[i].phi);
        EXPECT_EQ(evA[i].x, evB[i].x);
        EXPECT_EQ(evA[i].y, evB[i].y);
        EXPECT_EQ(evA[i].electrons, evB[i].electrons);
    }
    ASSERT_EQ(outA.size(), outB.size());
    for (size_t i = 0; i < outA.size(); ++i)
    {
        EXPECT_EQ(outA[i].idx, outB[i].idx);
        EXPECT_EQ(outA[i].electrons, outB[i].electrons);
    }

    // The next step draws other hits
    std::vector<Deposit> next;
    std::vector<CosmicDeposit::Event> evNext;
    a.generate(dt, g, next, &evNext);
    ASSERT_FALSE(evNext.empty());
    EXPECT_NE(evNext[0].x, evA[0].x);
}

// ── (h) Footprint ────────────────────────────────────────────────────────────────

TEST(CosmicDepositFootprintTest, NormalIncidenceIsOnePixel)
{
    const CosmicDeposit d(protonSpectrum(), siliconParams());
    const Geometry g(32, 32);
    std::vector<Deposit> out;
    d.depositEvent(hit(10.3, 20.7, 0.0, 0.0, 1000.0), g, out);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].idx, 20u * 32u + 10u);
    EXPECT_DOUBLE_EQ(out[0].electrons, 1000.0);
}

TEST(CosmicDepositFootprintTest, ChargeFollowsTheLengthInEachPixel)
{
    // From x = 10.5 over 3 pixels: 1/6, 1/3, 1/3, 1/6 of the track in pixels 10 to 13
    const CosmicDeposit d(protonSpectrum(), siliconParams());
    const Geometry g(32, 32);
    std::vector<Deposit> out;
    d.depositEvent(hit(10.5, 5.5, 3.0, 0.0, 600.0), g, out);

    ASSERT_EQ(out.size(), 4u);
    const std::vector<double> expected{100.0, 200.0, 200.0, 100.0};
    for (size_t i = 0; i < 4; ++i)
    {
        EXPECT_EQ(out[i].idx, 5u * 32u + 10u + i);
        EXPECT_NEAR(out[i].electrons, expected[i], 1e-9);
    }
}

TEST(CosmicDepositFootprintTest, ChargeIsConservedInsideTheDetector)
{
    const Geometry g(64, 64);
    const CosmicDeposit::Event e = hit(30.2, 40.6, 5.3, -3.7, 1234.0);

    const CosmicDeposit sharp(protonSpectrum(), siliconParams());
    std::vector<Deposit> out;
    sharp.depositEvent(e, g, out);
    EXPECT_NEAR(total(out), 1234.0, 1e-9);

    // With a lateral spread: more pixels, the same charge
    CosmicDeposit::Params p = siliconParams();
    p.sigmaDiffusion_px = 0.7;
    const CosmicDeposit spread(protonSpectrum(), p);
    std::vector<Deposit> outSpread;
    spread.depositEvent(e, g, outSpread);
    EXPECT_GT(outSpread.size(), out.size());
    EXPECT_NEAR(total(outSpread), 1234.0, 1e-9);
}

TEST(CosmicDepositFootprintTest, ChargeOutsideTheDetectorIsLost)
{
    // From x = 1.5 towards -2.5: 1/8 in pixel 1, 1/4 in pixel 0, the rest outside
    const CosmicDeposit d(protonSpectrum(), siliconParams());
    const Geometry g(16, 16);
    std::vector<Deposit> out;
    d.depositEvent(hit(1.5, 1.5, -4.0, 0.0, 400.0), g, out);
    ASSERT_EQ(out.size(), 2u);
    EXPECT_NEAR(total(out), 150.0, 1e-9);
}

TEST(CosmicDepositFootprintTest, EachSpeciesHasItsOwnSpread)
{
    // Source sigma 0.7 (R = 3, 7 x 7 pixels); helium sharp
    CosmicSpectrum s = protonSpectrum();
    s.add(helium, {10.0, 100.0}, {1.0, 0.01}, AngularModel(), footprint(0.0));
    CosmicDeposit::Params p = siliconParams();
    p.sigmaDiffusion_px = 0.7;
    const CosmicDeposit d(s, p);
    const Geometry g(32, 32);

    CosmicDeposit::Event e = hit(15.5, 15.5, 0.0, 0.0, 1000.0);
    std::vector<Deposit> outP;
    d.depositEvent(e, g, outP);
    EXPECT_EQ(outP.size(), 49u);
    EXPECT_NEAR(total(outP), 1000.0, 1e-9);

    e.species = 1;
    std::vector<Deposit> outHe;
    d.depositEvent(e, g, outHe);
    ASSERT_EQ(outHe.size(), 1u);
    EXPECT_DOUBLE_EQ(outHe[0].electrons, 1000.0);

    // A species outside the spectrum gets the source's sigma
    e.species = 7;
    std::vector<Deposit> outOther;
    d.depositEvent(e, g, outOther);
    const std::map<size_t, double> a = perPixel(outP), b = perPixel(outOther);
    ASSERT_EQ(a.size(), b.size());
    for (const auto& [idx, electrons] : a)
        EXPECT_NEAR(b.at(idx), electrons, 1e-9);
}

TEST(CosmicDepositFootprintTest, HaloHoldsItsShareOfTheCharge)
{
    // Proton: sharp core and a halo of sigma 2 holding 25 %; helium: that halo alone (sigma 2, no halo)
    CosmicSpectrum s;
    s.add(proton, {10.0, 100.0}, {1.0, 0.01}, AngularModel(), footprint(0.0, 2.0, 0.25));
    s.add(helium, {10.0, 100.0}, {1.0, 0.01}, AngularModel(), footprint(2.0));
    const CosmicDeposit d(s, siliconParams());
    const Geometry g(64, 64);

    CosmicDeposit::Event e = hit(32.5, 32.5, 0.0, 0.0, 1000.0);
    std::vector<Deposit> outHalo;
    d.depositEvent(e, g, outHalo);
    e.species = 1;
    std::vector<Deposit> outWide;
    d.depositEvent(e, g, outWide);

    EXPECT_EQ(outHalo.size(), 13u * 13u);   // R = ceil(3 sigma) = 6
    EXPECT_NEAR(total(outHalo), 1000.0, 1e-9);

    // 75 % in the pixel hit, 25 % spread as the wide Gaussian
    const size_t centre = 32u * 64u + 32u;
    const std::map<size_t, double> halo = perPixel(outHalo), wide = perPixel(outWide);
    ASSERT_EQ(halo.size(), wide.size());
    for (const auto& [idx, electrons] : wide)
        EXPECT_NEAR(halo.at(idx), 0.25 * electrons + (idx == centre ? 750.0 : 0.0), 1e-9);

    // A thin track: the same profile around each pixel, the charge kept
    e = hit(20.2, 30.7, 12.4, 3.1, 5000.0);   // species 0
    std::vector<Deposit> outTrack;
    d.depositEvent(e, g, outTrack);
    EXPECT_NEAR(total(outTrack), 5000.0, 1e-9);
}

TEST(CosmicDepositFootprintTest, FromTheConfiguration)
{
    const TempFile csv("p2p_gTest_cosmic_footprint.csv");
    writeSpectrumCsv(csv);

    json p = filtered(entry("p", 1, 1, csv.str()), "p");
    p["sigmaDiffusion_px"] = 1.0;
    json cosmic = {{"depth_um", 10.0}, {"pitch_um", 18.0}, {"eVPerElectron", 3.6}};
    cosmic["spectrum"] = json::array({p, filtered(entry("He", 2, 4, csv.str()), "He")});

    CosmicDeposit d;
    d.configure(PixelConfig(json{{"cosmic", cosmic}}));
    ASSERT_TRUE(d.isConfigured());
    const Geometry g(32, 32);
    const auto pixels = [&d, &g](size_t species)
    {
        CosmicDeposit::Event e = hit(15.5, 15.5, 0.0, 0.0, 1000.0);
        e.species = species;
        std::vector<Deposit> out;
        d.depositEvent(e, g, out);
        return out.size();
    };
    EXPECT_EQ(pixels(0), 49u);   // its own sigma 1: R = 3
    EXPECT_EQ(pixels(1), 1u);    // the source's sigma: 0

    // A new source sigma reaches the species that give none, not the others
    d.configure(PixelConfig(json{{"cosmic", json{{"sigmaDiffusion_px", 0.5}}}}));
    EXPECT_EQ(pixels(0), 49u);
    EXPECT_EQ(pixels(1), 25u);   // sigma 0.5: R = ceil(1.5) = 2, 5 x 5
}

TEST(CosmicDeltaRayTest, ShareFromTheVelocity)
{
    CosmicSpectrum s;
    s.add(proton, {10.0, 1.0e4}, {1.0, 1.0e-6}, AngularModel(), withDeltaRays());
    s.add(helium, {10.0, 1.0e4}, {1.0, 1.0e-6}, AngularModel(), withDeltaRays());
    s.add(helium, {10.0, 1.0e4}, {1.0, 1.0e-6});   // without delta rays
    const CosmicDeposit d(s, siliconParams());

    // A 5 MeV proton gives at most Wmax = 11 keV, below T0 (CSDA range of half a pixel, about 30 keV): no halo
    EXPECT_DOUBLE_EQ(d.deltaRayShare(0, 5.0), 0.0);

    for (double E : {100.0, 1000.0, 1.0e4})
    {
        const double h = d.deltaRayShare(0, E);
        EXPECT_GT(h, 0.0) << E;
        EXPECT_LT(h, 0.5) << E;
        // h follows the velocity, not the charge: helium at the same energy per nucleon (Wmax within 1 %)
        EXPECT_NEAR(d.deltaRayShare(1, E), h, 1e-2 * h) << E;
        EXPECT_DOUBLE_EQ(d.deltaRayShare(2, E), 0.0) << E;
    }
    EXPECT_THROW((void)d.deltaRayShare(3, 100.0), std::out_of_range);

    // Heavy particles only
    CosmicSpectrum leptons;
    EXPECT_THROW(leptons.add(electron, {0.1, 1.0}, {1.0, 0.1}, AngularModel(), withDeltaRays()), std::invalid_argument);
    EXPECT_THROW(leptons.add(positron, {0.1, 1.0}, {1.0, 0.1}, AngularModel(), withDeltaRays()), std::invalid_argument);
    EXPECT_TRUE(leptons.empty());
}

TEST(CosmicDeltaRayTest, HaloKeepsTheChargeAndTheSymmetry)
{
    CosmicSpectrum s;
    s.add(proton, {10.0, 1.0e4}, {1.0, 1.0e-6}, AngularModel(), withDeltaRays());
    CosmicDeposit::Params p = siliconParams();
    p.deltaRayCut_e = 0.01;
    const CosmicDeposit d(s, p);
    const Geometry g(128, 128);

    CosmicDeposit::Event e = hit(64.5, 64.5, 0.0, 0.0, 1.0e5);   // normal incidence, 1 GeV
    e.energy = 1000.0;
    std::vector<Deposit> out;
    d.depositEvent(e, g, out);
    EXPECT_NEAR(total(out), 1.0e5, 1e-6);

    const double h = d.deltaRayShare(0, 1000.0);
    const std::map<size_t, double> q = perPixel(out);
    const auto at = [&q](long a, long b)
    {
        const auto it = q.find(static_cast<size_t>((64 + b) * 128 + 64 + a));
        return (it == q.end()) ? 0.0 : it->second;
    };

    // The track keeps 1 - h, the halo spreads h around it, decreasing with the distance
    EXPECT_GT(at(0, 0), (1.0 - h) * 1.0e5);
    EXPECT_LT(at(0, 0), 1.0e5);
    EXPECT_GT(q.size(), 9u);
    EXPECT_GT(at(1, 0), at(2, 0));
    EXPECT_GT(at(2, 0), at(4, 0));
    EXPECT_GT(at(4, 0), 0.0);

    // Symmetries of the pixel grid
    for (long r : {1L, 2L, 3L, 5L})
    {
        EXPECT_NEAR(at(-r, 0),  at(r, 0), 1e-9 * at(r, 0)) << r;
        EXPECT_NEAR(at(0, r),   at(r, 0), 1e-9 * at(r, 0)) << r;
        EXPECT_NEAR(at(0, -r),  at(r, 0), 1e-9 * at(r, 0)) << r;
        EXPECT_NEAR(at(-r, -r), at(r, r), 1e-9 * at(r, r)) << r;
        EXPECT_NEAR(at(r, -r),  at(r, r), 1e-9 * at(r, r)) << r;
    }
}

TEST(CosmicDeltaRayTest, CutAndChargeSetTheRadius)
{
    CosmicSpectrum s;
    s.add(proton, {10.0, 1.0e4}, {1.0, 1.0e-6}, AngularModel(), withDeltaRays());
    const CosmicDeposit d(s, siliconParams());   // deltaRayCut_e = 1 e-
    const Geometry g(256, 256);
    const auto pixels = [&d, &g](double energy, double electrons, double dx)
    {
        CosmicDeposit::Event e = hit(128.5, 128.5, dx, 0.0, electrons);
        e.energy = energy;
        std::vector<Deposit> out;
        d.depositEvent(e, g, out);
        EXPECT_NEAR(total(out), electrons, 1e-9 * electrons);   // the charge is kept
        return perPixel(out).size();
    };

    EXPECT_EQ(pixels(1000.0, 10.0, 0.0), 1u);    // the halo puts less than 1 e- on every pixel
    EXPECT_EQ(pixels(5.0, 1.0e6, 0.0), 1u);      // slow: Wmax below T0
    const size_t small = pixels(1000.0, 1.0e4, 0.0);
    const size_t large = pixels(1000.0, 1.0e7, 0.0);
    EXPECT_GT(small, 1u);
    EXPECT_GT(large, small);                     // a larger charge reaches farther
    EXPECT_GT(pixels(1000.0, 1.0e6, 10.0), 11u); // an inclined track: a halo around each of its pixels
}

TEST(CosmicDeltaRayTest, FromTheConfiguration)
{
    const TempFile csv("p2p_gTest_cosmic_delta.csv");
    writeSpectrumCsv(csv);

    json p = filtered(entry("p", 1, 1, csv.str()), "p");
    p["deltaRays"] = true;
    json cosmic = {{"depth_um", 10.0}, {"pitch_um", 18.0}, {"eVPerElectron", 3.6}, {"deltaRayCut_e", 0.5}};
    cosmic["spectrum"] = json::array({p, filtered(entry("He", 2, 4, csv.str()), "He")});

    CosmicDeposit d;
    d.configure(PixelConfig(json{{"cosmic", cosmic}}));
    ASSERT_TRUE(d.isConfigured());
    EXPECT_DOUBLE_EQ(d.params().deltaRayCut_e, 0.5);
    EXPECT_TRUE(d.spectrum().footprints()[0].deltaRays);
    EXPECT_FALSE(d.spectrum().footprints()[1].deltaRays);
    EXPECT_GT(d.deltaRayShare(0, 100.0), 0.0);
    EXPECT_DOUBLE_EQ(d.deltaRayShare(1, 100.0), 0.0);

    const auto cfg = [](const json& j) { return PixelConfig(json{{"cosmic", j}}); };
    EXPECT_THROW(d.configure(cfg(json{{"deltaRayCut_e", 0.0}})), std::invalid_argument);   // not > 0
    json notBool = cosmic;
    notBool["spectrum"][0]["deltaRays"] = 1;
    EXPECT_THROW(d.configure(cfg(notBool)), std::invalid_argument);                       // not a boolean
    json lepton = cosmic;
    lepton["spectrum"] = json::array({json{{"species", "e-"}, {"Z", -1}, {"mc2", 0.51099895}, {"file", csv.str()},
                                           {"filterColumn", "NAME"}, {"filterValue", "p"}, {"colEnergy", "ENERGY"},
                                           {"colFlux", "FLUX"}, {"deltaRays", true}}});
    EXPECT_THROW(d.configure(cfg(lepton)), std::invalid_argument);                         // heavy particles only
    EXPECT_DOUBLE_EQ(d.params().deltaRayCut_e, 0.5);                                       // unchanged
    EXPECT_TRUE(d.spectrum().footprints()[0].deltaRays);
}

// ── (i) Energy loss ──────────────────────────────────────────────────────────────

TEST(CosmicDepositEnergyLossTest, LandauWidthScalesAsChargeSquaredAndPath)
{
    const CosmicDeposit d(protonSpectrum(), siliconParams());
    const CosmicDeposit::EnergyLoss p  = d.energyLoss(proton, 100.0, 10.0);
    const CosmicDeposit::EnergyLoss he = d.energyLoss(helium, 100.0, 10.0);
    EXPECT_NEAR(he.xi_eV / p.xi_eV, 4.0, 1e-12);
    EXPECT_NEAR(d.energyLoss(proton, 100.0, 20.0).xi_eV / p.xi_eV, 2.0, 1e-12);
    EXPECT_GT(p.mpv_eV, 0.0);
}

TEST(CosmicDepositEnergyLossTest, ElectronEnergyIsPerParticle)
{
    // A 10 MeV electron is close to minimum ionising; a 10 MeV proton is slow and loses far more.
    // (With an energy per nucleon, the electron would be as slow as the proton.)
    const CosmicDeposit d(protonSpectrum(), siliconParams());
    EXPECT_GT(d.energyLoss(proton, 10.0, 10.0).mpv_eV, 5.0 * d.energyLoss(electron, 10.0, 10.0).mpv_eV);
}

TEST(CosmicDepositEnergyLossTest, DepositNeverAboveTheKineticEnergy)
{
    // 10 to 20 keV protons: the Landau loss on 10 um is far above their energy
    CosmicSpectrum s;
    s.add(proton, {0.01, 0.02}, {1.0e6, 1.0e6});
    const CosmicDeposit d(s, siliconParams());
    const Geometry g(64, 64);

    std::mt19937_64 rng(5);
    const std::vector<CosmicDeposit::Event> events = d.drawEvents(stepFor(d, g, 100.0), g, rng);
    ASSERT_FALSE(events.empty());
    for (const CosmicDeposit::Event& e : events)
        EXPECT_LE(e.electrons, proton.kineticEnergy(e.energy) * 1.0e6 / 3.6 * (1.0 + 1e-12));
}

// ── (j) Outputs ──────────────────────────────────────────────────────────────────

TEST(CosmicDepositOutputTest, TableOfHits)
{
    CosmicSpectrum s = protonSpectrum(1.0e6);
    s.add(helium, {10.0, 100.0}, {1.0e5, 1.0e3});
    CosmicDeposit d(s, siliconParams());
    d.setSeed(11);
    const Geometry g(64, 64);

    std::vector<Deposit> out;
    std::vector<CosmicDeposit::Event> events;
    d.generate(stepFor(d, g, 30.0), g, out, &events);
    ASSERT_FALSE(events.empty());

    const std::vector<std::string> names{"SPECIES", "Z", "A", "ENERGY", "MU", "PHI", "FACE", "X", "Y", "DX", "DY", "PATH", "NE"};
    const std::vector<CosmicDeposit::ColumnInfo>& info = CosmicDeposit::columnInfo();
    ASSERT_EQ(info.size(), names.size());
    for (size_t i = 0; i < names.size(); ++i)
        EXPECT_EQ(info[i].name, names[i]);

    const CosmicDeposit::EventColumns c = d.toColumns(events);
    EXPECT_EQ(c.size(), events.size());
    EXPECT_EQ(c.species.size(), events.size());
    for (const std::vector<double>* v : {&c.Z, &c.A, &c.energy, &c.mu, &c.phi, &c.face, &c.x, &c.y, &c.dx, &c.dy, &c.path, &c.electrons})
        EXPECT_EQ(v->size(), events.size());

    for (size_t i = 0; i < events.size(); ++i)
    {
        const specie& sp = d.spectrum().species()[events[i].species];
        EXPECT_EQ(c.species[i], sp.name());
        EXPECT_EQ(c.Z[i], sp.Z());
        EXPECT_EQ(c.A[i], sp.A());
        EXPECT_EQ(c.energy[i], events[i].energy);
        EXPECT_EQ(c.path[i], events[i].path_um);
        EXPECT_EQ(c.face[i], static_cast<double>(events[i].face));
        EXPECT_EQ(c.electrons[i], events[i].electrons);
    }
}

TEST(CosmicDepositOutputTest, HeaderCards)
{
    CosmicDeposit d(protonSpectrum(), siliconParams());
    d.setSeed(8);

    std::set<std::string> keys;
    for (const phot2pix::HeaderCard& card : d.inputCards())
    {
        keys.insert(card.key);
        EXPECT_LE(card.key.size(), 8u) << card.key;
        if (card.key == "CRFLUX")
            EXPECT_NEAR(std::get<double>(card.value), 9.0, 1e-12);
        if (card.key == "CRACC")
            EXPECT_NEAR(std::get<double>(card.value), std::numbers::pi, 1e-12);   // isotropic, one face
        if (card.key == "CRTILT")
            EXPECT_DOUBLE_EQ(std::get<double>(card.value), 0.0);
        if (card.key == "CRSEED")
            EXPECT_EQ(std::get<size_t>(card.value), 8u);
    }
    for (const char* k : {"CRFLUX", "CRACC", "CRTILT", "CRDEPTH", "CRPITCH", "CREVPE", "CRSIGD", "CRMAXTRK", "CRDCUT", "CRTRUNC", "CRSEED"})
        EXPECT_EQ(keys.count(k), 1u) << k;
}

// ── (k) Stopping power and range ─────────────────────────────────────────────────

namespace
{
    // Material with the values of a PDG muon-table header (Atomic and Nuclear Properties, 2024).
    material pdgMaterial(double ZoverA, double I_eV, double rho, double a, double k, double x0, double x1, double Cbar, double delta0)
    {
        material m;
        m.ZoverA            = ZoverA;
        m.meanExcitation_eV = I_eV;
        m.density_gcm3      = rho;
        m.sternheimerA      = a;
        m.sternheimerM      = k;
        m.sternheimerX0     = x0;
        m.sternheimerX1     = x1;
        m.sternheimerC      = Cbar;
        m.sternheimerDelta0 = delta0;
        return m;
    }

    const specie muon = specie::fromRestEnergy("mu-", -1.0, 105.65839);

    // Source of muons from 10 MeV to 2 GeV in the material m: its range table covers the PDG rows.
    CosmicDeposit muonSource(const material& m)
    {
        CosmicSpectrum s;
        s.add(muon, {10.0, 2000.0}, {1.0, 1.0e-4});
        CosmicDeposit::Params p = siliconParams();
        p.material = m;
        return CosmicDeposit(s, p);
    }

    // Row of a PDG muon table: kinetic energy (MeV), ionization (MeV cm2/g), CSDA range (g/cm2).
    struct PdgRow { double T, ionization, csdaRange; };
}

TEST(StoppingTableTest, ConstantStoppingPower)
{
    // S = c above Tc, c sqrt(T / Tc) below: R = 2 sqrt(T Tc) / c below Tc, (T + Tc) / c above
    const double c = 2.0, Tc = 1.0;
    const StoppingTable t([c](double) { return c; }, Tc, 1000.0);
    ASSERT_FALSE(t.empty());
    EXPECT_DOUBLE_EQ(t.lowEnergyLimit(), Tc);

    EXPECT_NEAR(t.range(0.25), 2.0 * std::sqrt(0.25 * Tc) / c, 1e-12);
    for (double T : {1.0, 10.0, 37.3, 500.0, 999.0})
        EXPECT_NEAR(t.range(T), (T + Tc) / c, 1e-5 * (T + Tc) / c) << "T = " << T;

    // energy() is the exact inverse of range(), on both parts
    for (double T : {0.01, 0.5, 1.0, 3.7, 123.0, 999.0})
        EXPECT_NEAR(t.energy(t.range(T)), T, 1e-12 * T) << "T = " << T;

    EXPECT_DOUBLE_EQ(t.energyAfter(50.0, t.range(50.0)), 0.0);   // stopped
    EXPECT_NEAR(t.energyAfter(50.0, 0.0), 50.0, 1e-9);
    EXPECT_NEAR(t.energyAfter(50.0, 10.0 / c), 40.0, 1e-3);      // 10 MeV lost over 10 / c g/cm2

    EXPECT_EQ(t.range(0.0), 0.0);
    EXPECT_EQ(t.energy(0.0), 0.0);
    EXPECT_TRUE(StoppingTable().empty());
}

TEST(StoppingTableTest, RejectsBadInput)
{
    EXPECT_THROW((void)StoppingTable([](double) { return 1.0; }, 0.0, 10.0), std::invalid_argument);       // Tc = 0
    EXPECT_THROW((void)StoppingTable([](double) { return 1.0; }, 1.0, 10.0, 0), std::invalid_argument);    // no node
    EXPECT_THROW((void)StoppingTable([](double T) { return 2.0 - T; }, 1.0, 10.0), std::invalid_argument); // S <= 0 above 2 MeV
}

TEST(CosmicStoppingTest, BetheAgainstPdgMuonTables)
{
    struct Case
    {
        const char*         name;
        material            m;
        std::vector<PdgRow> rows;
    };
    // Rows copied from pdg.lbl.gov/2024/AtomicNuclearProperties/MUE/muE_{silicon_Si,germanium_Ge,cadmium_telluride_CdTe}.txt
    const std::vector<Case> cases = {
        {"Si",   pdgMaterial(14.0 / 28.0855, 173.0, 2.329, 0.1492, 3.2546, 0.2015, 2.8716, 4.4355, 0.14),
                 {{10.0, 6.363, 0.8779}, {100.0, 1.849, 38.22}, {1000.0, 1.803, 559.5}}},
        {"Ge",   pdgMaterial(32.0 / 72.630, 350.0, 5.323, 0.0719, 3.3306, 0.3376, 3.6096, 5.1411, 0.14),
                 {{10.0, 5.033, 1.124}, {100.0, 1.510, 47.20}, {1000.0, 1.504, 676.7}}},
        {"CdTe", pdgMaterial(0.41665, 539.3, 5.850, 0.2484, 2.6665, 0.0564, 3.2962, 5.9677, 0.00),
                 {{10.0, 4.453, 1.273}, {100.0, 1.359, 52.72}, {1000.0, 1.387, 741.2}}},
    };

    for (const Case& c : cases)
    {
        const CosmicDeposit d = muonSource(c.m);
        for (const PdgRow& r : c.rows)
        {
            EXPECT_NEAR(d.stoppingPower(muon, r.T), r.ionization, 0.01 * r.ionization) << c.name << ", T = " << r.T << " MeV";

            // At 10 MeV the low-energy part (no shell corrections, sqrt(T) below the peak) weighs more in the range
            const double tolerance = (r.T < 50.0) ? 0.03 : 0.02;
            EXPECT_NEAR(d.stoppingTable(0).range(r.T), r.csdaRange, tolerance * r.csdaRange) << c.name << ", T = " << r.T << " MeV";
        }
    }
}

TEST(CosmicStoppingTest, HgCdTeOfCadmiumFractionOneStopsLikeCdTe)
{
    // PDG CdTe; at 1 GeV the Sternheimer-Peierls delta of compound() differs from the PDG fit by about 0.4
    const CosmicDeposit d = muonSource(material::HgCdTe(1.0));
    EXPECT_NEAR(d.stoppingPower(muon, 10.0),   4.453, 0.01 * 4.453);
    EXPECT_NEAR(d.stoppingPower(muon, 100.0),  1.359, 0.01 * 1.359);
    EXPECT_NEAR(d.stoppingPower(muon, 1000.0), 1.387, 0.03 * 1.387);
}

TEST(CosmicStoppingTest, OneTablePerSpeciesAndRangeScaling)
{
    CosmicSpectrum s = protonSpectrum();
    s.add(helium, {10.0, 100.0}, {1.0, 0.01});
    const CosmicDeposit d(s, siliconParams());

    EXPECT_FALSE(d.stoppingTable(0).empty());
    EXPECT_FALSE(d.stoppingTable(1).empty());
    EXPECT_THROW((void)d.stoppingTable(2), std::out_of_range);

    // At the same speed, the range scales as M / z^2: 4 / 4 = 1 for helium against the proton
    const double Rp = d.stoppingTable(0).range(50.0);
    EXPECT_NEAR(d.stoppingTable(1).range(4.0 * 50.0), Rp, 1e-3 * Rp);
}

TEST(CosmicStoppingTest, ElectronsAndPositrons)
{
    // No independent table for electrons: consistency only
    CosmicSpectrum s;
    s.add(electron, {0.01, 10.0}, {1.0, 1.0e-3});
    s.add(positron, {0.01, 10.0}, {1.0, 1.0e-3});
    const CosmicDeposit d(s, siliconParams());

    for (double T : {0.01, 0.1, 1.0, 10.0})
    {
        EXPECT_GT(d.stoppingPower(electron, T), 0.0) << "T = " << T;
        EXPECT_GT(d.stoppingPower(positron, T), 0.0) << "T = " << T;
    }
    EXPECT_GT(d.stoppingPower(electron, 0.01), d.stoppingPower(electron, 0.1));   // falls up to minimum ionisation
    EXPECT_GT(d.stoppingPower(electron, 0.1), d.stoppingPower(electron, 1.0));
    EXPECT_NE(d.stoppingPower(electron, 0.1), d.stoppingPower(positron, 0.1));

    for (size_t i : {0u, 1u})
    {
        double previous = 0.0;
        for (double T : {0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0})
        {
            const double R = d.stoppingTable(i).range(T);
            EXPECT_GT(R, previous) << "species " << i << ", T = " << T;
            previous = R;
        }
    }
}

// ── (l) Range cut, thick crossings, Bragg profile ────────────────────────────────

namespace
{
    // Protons of one kinetic energy T (MeV) in silicon of the given depth.
    CosmicDeposit monoProtons(double T, double depth_um)
    {
        CosmicSpectrum s;
        s.add(proton, {T, T * 1.0001}, {1.0e6, 1.0e6});
        CosmicDeposit::Params p = siliconParams();
        p.depth_um = depth_um;
        return CosmicDeposit(s, p);
    }

    // CSDA range (um) of the protons of d.
    double rangeUm(const CosmicDeposit& d, double T)
    {
        return d.stoppingTable(0).range(T) / (1.0e-4 * d.params().material.density_gcm3);
    }
}

TEST(CosmicBraggTest, StoppingParticlesDepositAllTheirEnergy)
{
    const CosmicDeposit d = monoProtons(0.3, 10.0);   // a range of a few um: every proton stops in 10 um
    const Geometry g(64, 64);
    std::mt19937_64 rng(31);
    const std::vector<CosmicDeposit::Event> events = d.drawEvents(stepFor(d, g, 200.0), g, rng);
    ASSERT_FALSE(events.empty());

    for (const CosmicDeposit::Event& e : events)
    {
        const double T = proton.kineticEnergy(e.energy);
        const double R = rangeUm(d, T);
        ASSERT_LT(R, d.params().depth_um);
        EXPECT_NEAR(e.path_um, R, 1e-9 * R);
        EXPECT_NEAR(e.electrons, T * 1.0e6 / 3.6, 1e-9 * e.electrons);
        // Projected track of the path R: R sin(theta) / pitch pixels
        EXPECT_NEAR(std::hypot(e.dx, e.dy), R * std::sqrt(1.0 - e.mu * e.mu) / d.params().pitch_um, 1e-9);
    }
}

TEST(CosmicBraggTest, FastParticlesCrossTheDepth)
{
    const CosmicDeposit d = monoProtons(100.0, 10.0);
    const Geometry g(64, 64);
    std::mt19937_64 rng(32);
    const std::vector<CosmicDeposit::Event> events = d.drawEvents(stepFor(d, g, 200.0), g, rng);
    ASSERT_FALSE(events.empty());

    for (const CosmicDeposit::Event& e : events)
    {
        const double L = d.params().depth_um / e.mu;
        EXPECT_NEAR(e.path_um, L, 1e-12 * L);
        EXPECT_LT(e.electrons, 100.0 * 1.0e6 / 3.6);
    }
}

TEST(CosmicBraggTest, ThickCrossingsAverageTheCsdaLoss)
{
    // 3 MeV protons in 50 um of silicon: they stop, or lose far more than 10 % of their energy
    const CosmicDeposit d = monoProtons(3.0, 50.0);
    const Geometry g(256, 256);
    const double um_to_gcm2 = 1.0e-4 * d.params().material.density_gcm3;
    std::mt19937_64 rng(7);

    double sumDeposit = 0.0, sumMean = 0.0, sumSquares = 0.0;
    size_t nThick = 0;
    for (int step = 0; step < 20 && nThick < 2000; ++step)
        for (const CosmicDeposit::Event& e : d.drawEvents(stepFor(d, g, 500.0), g, rng))
        {
            const double T = proton.kineticEnergy(e.energy);
            const double L = d.params().depth_um / e.mu;
            if (rangeUm(d, T) <= L)
                continue;   // stops
            const double mean = T - d.stoppingTable(0).energyAfter(T, L * um_to_gcm2);
            if (!(mean > 0.1 * T))
                continue;   // thin
            const double deposit = e.electrons * 3.6e-6;   // MeV
            sumDeposit += deposit;
            sumMean    += mean;
            sumSquares += (deposit - mean) * (deposit - mean);
            ++nThick;
        }

    ASSERT_GT(nThick, 500u);
    EXPECT_NEAR(sumDeposit / sumMean, 1.0, 1e-2);   // Gaussian around the CSDA mean loss
    EXPECT_GT(sumSquares, 0.0);                     // with straggling
}

TEST(CosmicBraggTest, ChargeRisesToABraggPeak)
{
    const double T = 5.0;
    const CosmicDeposit d = monoProtons(T, 300.0);
    const double R     = rangeUm(d, T);
    const double pitch = d.params().pitch_um;
    ASSERT_GT(R, 5.0 * pitch);   // several pixels

    // Horizontal track from x = 0.25 along +x, over the whole range
    CosmicDeposit::Event e = hit(0.25, 3.5, R / pitch, 0.0, 1.0e5);
    e.energy  = T;
    e.path_um = R;
    std::vector<Deposit> out;
    d.depositEvent(e, Geometry(64, 8), out);

    ASSERT_GE(out.size(), 6u);
    EXPECT_NEAR(total(out), 1.0e5, 1e-6);
    // The first and the last pixels are crossed in part; over the full ones the charge rises, up to the last
    for (size_t i = 2; i + 1 < out.size(); ++i)
        EXPECT_GT(out[i].electrons, out[i - 1].electrons) << "pixel " << i;
}

TEST(CosmicBraggTest, FastParticleChargeIsUniform)
{
    const CosmicDeposit d = monoProtons(1000.0, 300.0);
    CosmicDeposit::Event e = hit(10.0, 3.5, 4.0, 0.0, 4000.0);   // four full pixels
    e.energy  = 1000.0;
    e.path_um = 4.0 * d.params().pitch_um;
    std::vector<Deposit> out;
    d.depositEvent(e, Geometry(64, 8), out);

    ASSERT_EQ(out.size(), 4u);
    for (const Deposit& p : out)
        EXPECT_NEAR(p.electrons, 1000.0, 5.0);
}

TEST(CosmicBraggTest, TablesMustMatchTheSpectrum)
{
    // A source whose tables were never built (as a derived class that fills m_spectrum itself would be)
    struct Raw : public CosmicDeposit
    {
        Raw()
        {
            m_params = siliconParams();
            m_spectrum = protonSpectrum();
        }
    };
    const Raw raw;
    std::mt19937_64 rng(1);
    EXPECT_THROW(raw.drawEvents(1.0, Geometry(8, 8), rng), std::logic_error);
}

// ── (m) Angular models ───────────────────────────────────────────────────────────

namespace
{
    // Fast protons (no stop) with an angular model, on a sensor tilted by tilt_deg; maxTrackPixels large
    // enough to make the truncation negligible unless given.
    CosmicDeposit angularSource(const AngularModel& a, double tilt_deg = 0.0, double maxTrackPixels = 1.0e4)
    {
        CosmicSpectrum s;
        s.add(proton, {1000.0, 2000.0}, {1.0e6, 0.5e6}, a);
        CosmicDeposit::Params p = siliconParams();
        p.tilt_deg       = tilt_deg;
        p.maxTrackPixels = maxTrackPixels;
        return CosmicDeposit(s, p);
    }

    // Hits of about n steps of 2000 hits each.
    std::vector<CosmicDeposit::Event> drawMany(const CosmicDeposit& d, size_t steps, uint64_t seed)
    {
        const Geometry g(256, 256);
        std::mt19937_64 rng(seed);
        std::vector<CosmicDeposit::Event> all;
        for (size_t k = 0; k < steps; ++k)
        {
            const std::vector<CosmicDeposit::Event> ev = d.drawEvents(stepFor(d, g, 2000.0), g, rng);
            all.insert(all.end(), ev.begin(), ev.end());
        }
        return all;
    }

    constexpr double kPi = std::numbers::pi;
}

TEST(CosmicAngularTest, ModelsAndTheirChecks)
{
    EXPECT_EQ(AngularModel().kind(), AngularModel::Kind::isotropic);
    EXPECT_EQ(AngularModel().faces(), AngularModel::Faces::front);
    EXPECT_FALSE(AngularModel::isotropic().zenith());
    EXPECT_TRUE(AngularModel::cosn(2.0).zenith());

    const AngularModel c2 = AngularModel::cosn(2.0);
    EXPECT_NEAR(c2.relativeIntensity(1.0), 1.0, 1e-15);
    EXPECT_NEAR(c2.relativeIntensity(0.5), 0.25, 1e-15);
    EXPECT_EQ(c2.relativeIntensity(-0.5), 0.0);   // upward-going

    const AngularModel t = AngularModel::table({0.0, 60.0, 90.0}, {1.0, 0.5, 0.0});
    EXPECT_NEAR(t.relativeIntensity(1.0), 1.0, 1e-12);
    EXPECT_NEAR(t.relativeIntensity(0.5), 0.5, 1e-9);            // 60 deg
    EXPECT_NEAR(t.relativeIntensity(std::cos(30.0 * kPi / 180.0)), 0.75, 1e-9);
    EXPECT_DOUBLE_EQ(t.maxRelativeIntensity(), 1.0);

    EXPECT_THROW(AngularModel::cosn(-1.0), std::invalid_argument);
    EXPECT_THROW(AngularModel::table({0.0}, {1.0}), std::invalid_argument);                    // one point
    EXPECT_THROW(AngularModel::table({0.0, 95.0}, {1.0, 1.0}), std::invalid_argument);         // above 90 deg
    EXPECT_THROW(AngularModel::table({10.0, 0.0}, {1.0, 1.0}), std::invalid_argument);         // not ascending
    EXPECT_THROW(AngularModel::table({0.0, 90.0}, {1.0, -1.0}), std::invalid_argument);        // negative
    EXPECT_THROW(AngularModel::table({0.0, 90.0}, {0.0, 0.0}), std::invalid_argument);         // all 0
}

TEST(CosmicAngularTest, GeometricFactors)
{
    using F = AngularModel::Faces;
    EXPECT_NEAR(angularSource(AngularModel::isotropic()).geometricFactor(0),          kPi,       1e-12);
    EXPECT_NEAR(angularSource(AngularModel::isotropic(F::back)).geometricFactor(0),   kPi,       1e-12);
    EXPECT_NEAR(angularSource(AngularModel::isotropic(F::both)).geometricFactor(0),   2.0 * kPi, 1e-12);
    EXPECT_NEAR(angularSource(AngularModel::surface()).geometricFactor(0),            2.0 * kPi, 1e-12);
    EXPECT_NEAR(angularSource(AngularModel::cosn(2.0)).geometricFactor(0),            kPi / 2.0, 1e-12);
    EXPECT_NEAR(angularSource(AngularModel::cosn(2.0, F::back)).geometricFactor(0),   0.0,       1e-12);   // nothing from below

    // Numerical: cos^2 tabulated on a horizontal sensor, and cos^2 on a vertical sensor (front: pi / 8)
    std::vector<double> theta, f;
    for (int k = 0; k <= 90; ++k)
    {
        theta.push_back(k);
        f.push_back(std::pow(std::cos(k * kPi / 180.0), 2));
    }
    EXPECT_NEAR(angularSource(AngularModel::table(theta, f)).geometricFactor(0), kPi / 2.0, 1e-3 * kPi / 2.0);
    EXPECT_NEAR(angularSource(AngularModel::cosn(2.0), 90.0).geometricFactor(0), kPi / 8.0, 1e-3 * kPi / 8.0);
    EXPECT_NEAR(angularSource(AngularModel::cosn(2.0, F::both), 90.0).geometricFactor(0), kPi / 4.0, 1e-3 * kPi / 4.0);

    // The rate follows G: G * ∫ J dE * area
    const CosmicDeposit d = angularSource(AngularModel::cosn(2.0));
    EXPECT_NEAR(d.eventRate(1.0e-4), kPi / 2.0 * d.spectrum().integralIntensity() * 1.0e-4, 1e-9 * d.eventRate(1.0e-4));
}

TEST(CosmicAngularTest, TruncatedShares)
{
    // maxTrackPixels = 8, depth 10 um, pitch 18 um: q = 14.4, muMin = 1 / sqrt(1 + q^2)
    const double q = 8.0 * 18.0 / 10.0, mu2 = 1.0 / (1.0 + q * q), muMin = std::sqrt(mu2);
    EXPECT_NEAR(angularSource(AngularModel::isotropic(), 0.0, 8.0).truncatedFraction(0), mu2, 1e-15);
    EXPECT_NEAR(angularSource(AngularModel::surface(), 0.0, 8.0).truncatedFraction(0), muMin, 1e-15);
    EXPECT_NEAR(angularSource(AngularModel::cosn(2.0), 0.0, 8.0).truncatedFraction(0), std::pow(muMin, 4.0), 1e-15);
}

TEST(CosmicAngularTest, DirectionsFollowTheModel)
{
    // mu above muMin: surface mu uniform, cos^2 on a horizontal sensor mu^4 uniform
    const double q = 8.0 * 18.0 / 10.0, muMin = 1.0 / std::sqrt(1.0 + q * q);

    std::vector<double> u;
    for (const CosmicDeposit::Event& e : drawMany(angularSource(AngularModel::surface(), 0.0, 8.0), 5, 41))
        u.push_back((e.mu - muMin) / (1.0 - muMin));
    EXPECT_GT(DST::Math::ks_uniform_pvalue(u, 0.0, 1.0), kAlpha);

    u.clear();
    const double low = std::pow(muMin, 4.0);
    for (const CosmicDeposit::Event& e : drawMany(angularSource(AngularModel::cosn(2.0), 0.0, 8.0), 5, 42))
    {
        u.push_back((std::pow(e.mu, 4.0) - low) / (1.0 - low));
        ASSERT_EQ(e.face, 1);
    }
    EXPECT_GT(DST::Math::ks_uniform_pvalue(u, 0.0, 1.0), kAlpha);
}

TEST(CosmicAngularTest, Faces)
{
    using F = AngularModel::Faces;
    for (const CosmicDeposit::Event& e : drawMany(angularSource(AngularModel::isotropic(F::front)), 1, 51))
        ASSERT_EQ(e.face, 1);
    for (const CosmicDeposit::Event& e : drawMany(angularSource(AngularModel::isotropic(F::back)), 1, 52))
        ASSERT_EQ(e.face, -1);

    // Both faces alike: binomial share of the back face
    const std::vector<CosmicDeposit::Event> both = drawMany(angularSource(AngularModel::isotropic(F::both)), 5, 53);
    double back = 0.0;
    for (const CosmicDeposit::Event& e : both)
        back += (e.face == -1) ? 1.0 : 0.0;
    const double n = static_cast<double>(both.size());
    EXPECT_NEAR(back / n, 0.5, 4.0 * 0.5 / std::sqrt(n));
}

TEST(CosmicAngularTest, VerticalSensorUnderCosSquare)
{
    // Sensor normal horizontal (tilt 90 deg), I = J cos^2(theta_z): E[mu] = (pi / 15) / (pi / 8) = 8 / 15,
    // E[mu^2] = 1 / 3, so sigma(mu) = sqrt(1/3 - (8/15)^2); both faces alike
    const std::vector<CosmicDeposit::Event> ev = drawMany(angularSource(AngularModel::cosn(2.0, AngularModel::Faces::both), 90.0), 10, 61);
    ASSERT_GT(ev.size(), 5000u);

    double sum = 0.0, back = 0.0;
    for (const CosmicDeposit::Event& e : ev)
    {
        sum  += e.mu;
        back += (e.face == -1) ? 1.0 : 0.0;
    }
    const double n = static_cast<double>(ev.size());
    const double sigma = std::sqrt(1.0 / 3.0 - (8.0 / 15.0) * (8.0 / 15.0));
    EXPECT_NEAR(sum / n, 8.0 / 15.0, 4.0 * sigma / std::sqrt(n));
    EXPECT_NEAR(back / n, 0.5, 4.0 * 0.5 / std::sqrt(n));
}

TEST(CosmicAngularTest, FromTheConfiguration)
{
    const TempFile csv("p2p_gTest_cosmic_angular.csv");
    writeSpectrumCsv(csv);
    const TempFile table("p2p_gTest_cosmic_angular_table.csv");
    std::ofstream(table.path) << "THETA [deg],I\n0,1\n45,0.5\n90,0\n";

    json cosn = filtered(entry("p", 1, 1, csv.str()), "p");
    cosn["angular"] = {{"model", "cosn"}, {"n", 2.0}};
    cosn["faces"]   = "both";
    json tab = filtered(entry("He", 2, 4, csv.str()), "He");
    tab["angular"] = {{"model", "table"}, {"file", table.str()}, {"colTheta", "THETA"}, {"colIntensity", "I"}};

    const CosmicSpectrum s = CosmicSpectrum::load(spectrumConfig({cosn, tab}));
    ASSERT_EQ(s.angular().size(), 2u);
    EXPECT_EQ(s.angular()[0].kind(), AngularModel::Kind::cosn);
    EXPECT_DOUBLE_EQ(s.angular()[0].n(), 2.0);
    EXPECT_EQ(s.angular()[0].faces(), AngularModel::Faces::both);
    EXPECT_EQ(s.angular()[1].kind(), AngularModel::Kind::table);
    EXPECT_EQ(s.angular()[1].faces(), AngularModel::Faces::front);
    EXPECT_NEAR(s.angular()[1].relativeIntensity(std::cos(45.0 * kPi / 180.0)), 0.5, 1e-9);

    const auto load = [](const json& e) { return CosmicSpectrum::load(spectrumConfig({e})); };
    json bad = filtered(entry("p", 1, 1, csv.str()), "p");
    bad["angular"] = {{"model", "beam"}};
    EXPECT_THROW(load(bad), std::invalid_argument);   // unknown model
    bad["angular"] = {{"model", "cosn"}};
    EXPECT_THROW(load(bad), std::invalid_argument);   // no n
    bad.erase("angular");
    bad["faces"] = "top";
    EXPECT_THROW(load(bad), std::invalid_argument);   // unknown face
}

// ── (n) Geometry helpers ─────────────────────────────────────────────────────────

TEST(CosmicGeometryTest, IncidenceOfABeam)
{
    using phot2pix::CosmicRay::incidence;
    using phot2pix::CosmicRay::Incidence;

    // Vertical beam on a horizontal sensor: normal incidence, front face
    Incidence in = incidence(0.0, 0.0);
    EXPECT_NEAR(in.mu, 1.0, 1e-15);
    EXPECT_EQ(in.face, 1);

    // From 60 deg on the +x side: mu = 1/2, the track runs towards -x
    in = incidence(60.0, 0.0);
    EXPECT_NEAR(in.mu, 0.5, 1e-12);
    EXPECT_NEAR(in.phi, kPi, 1e-12);
    EXPECT_EQ(in.face, 1);

    // From below: back face
    in = incidence(120.0, 0.0);
    EXPECT_NEAR(in.mu, 0.5, 1e-12);
    EXPECT_EQ(in.face, -1);

    // Sensor tilted by 30 deg towards +x, beam along its normal
    EXPECT_NEAR(incidence(30.0, 0.0, 30.0).mu, 1.0, 1e-12);

    // Vertical sensor (tilt 90): a vertical beam grazes it, a horizontal one from +x hits it at normal incidence
    EXPECT_NEAR(incidence(0.0, 0.0, 90.0).mu, 0.0, 1e-12);
    EXPECT_NEAR(incidence(90.0, 0.0, 90.0).mu, 1.0, 1e-12);
    EXPECT_EQ(incidence(90.0, 180.0, 90.0).face, -1);   // from -x: back face

    EXPECT_THROW(incidence(190.0, 0.0), std::invalid_argument);
    EXPECT_THROW(incidence(10.0, 0.0, -5.0), std::invalid_argument);
}

TEST(CosmicGeometryTest, EffectiveArea)
{
    using phot2pix::CosmicRay::effectiveArea;
    EXPECT_NEAR(effectiveArea(2.0e-3, 0.0, 0.0), 2.0e-3, 1e-18);
    EXPECT_NEAR(effectiveArea(2.0e-3, 60.0, 45.0), 1.0e-3, 1e-15);           // A cos(60)
    EXPECT_NEAR(effectiveArea(2.0e-3, 0.0, 0.0, 60.0), 1.0e-3, 1e-15);       // tilted sensor, vertical beam
    EXPECT_THROW(effectiveArea(0.0, 0.0, 0.0), std::invalid_argument);
}

TEST(CosmicGeometryTest, GeometricFactorOfTheModels)
{
    using phot2pix::CosmicRay::geometricFactor;
    using F = AngularModel::Faces;

    // Full sky, closed forms
    EXPECT_NEAR(geometricFactor(AngularModel::isotropic()),          kPi,       1e-12);
    EXPECT_NEAR(geometricFactor(AngularModel::isotropic(F::both)),   2.0 * kPi, 1e-12);
    EXPECT_NEAR(geometricFactor(AngularModel::isotropic(), 37.0),    kPi,       1e-12);   // sensor frame: no tilt
    EXPECT_NEAR(geometricFactor(AngularModel::surface()),            2.0 * kPi, 1e-12);
    EXPECT_NEAR(geometricFactor(AngularModel::cosn(2.0)),            kPi / 2.0, 1e-12);

    // Window of incidence angles
    EXPECT_NEAR(geometricFactor(AngularModel::isotropic(), 0.0, 0.0, 30.0),  kPi / 4.0,           1e-12);   // pi (1 - 3/4)
    EXPECT_NEAR(geometricFactor(AngularModel::surface(), 0.0, 60.0, 90.0),   kPi,                 1e-12);   // 2 pi cos(60)
    EXPECT_NEAR(geometricFactor(AngularModel::cosn(2.0), 0.0, 0.0, 60.0),    15.0 * kPi / 32.0,   1e-12);   // pi/2 (1 - 1/16)

    // Numerical: tilted cos^2 (vertical sensor: pi / 8 per face), and the window splits the full sky
    EXPECT_NEAR(geometricFactor(AngularModel::cosn(2.0), 90.0),          kPi / 8.0, 1e-3 * kPi / 8.0);
    EXPECT_NEAR(geometricFactor(AngularModel::cosn(2.0, F::both), 90.0), kPi / 4.0, 1e-3 * kPi / 4.0);
    const double inner = geometricFactor(AngularModel::cosn(2.0), 45.0, 0.0, 40.0);
    const double outer = geometricFactor(AngularModel::cosn(2.0), 45.0, 40.0, 90.0);
    EXPECT_NEAR(inner + outer, geometricFactor(AngularModel::cosn(2.0), 45.0), 2e-3 * (inner + outer));

    // The source uses the same geometric factor
    EXPECT_NEAR(angularSource(AngularModel::cosn(2.0), 90.0).geometricFactor(0), geometricFactor(AngularModel::cosn(2.0), 90.0), 1e-12);

    EXPECT_THROW(geometricFactor(AngularModel::isotropic(), 0.0, 30.0, 30.0), std::invalid_argument);   // empty window
    EXPECT_THROW(geometricFactor(AngularModel::isotropic(), 0.0, 0.0, 95.0), std::invalid_argument);
    EXPECT_THROW(geometricFactor(AngularModel::isotropic(), 200.0), std::invalid_argument);
}

TEST(CosmicGeometryTest, Etendue)
{
    using phot2pix::CosmicRay::etendue;
    using phot2pix::CosmicRay::geometricFactor;

    // A 2048 x 2048 matrix of 18 um pixels, isotropic field on one face
    const double area = std::pow(2048.0 * 18.0e-6, 2);
    EXPECT_NEAR(etendue(AngularModel::isotropic(), area), kPi * area, 1e-15);
    EXPECT_NEAR(etendue(AngularModel::cosn(2.0), area, 90.0), area * geometricFactor(AngularModel::cosn(2.0), 90.0), 1e-15);
    EXPECT_THROW(etendue(AngularModel::isotropic(), -1.0), std::invalid_argument);
}
