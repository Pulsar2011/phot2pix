//
// Created by gillard on 2026-11-27.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_detectorExpose.cxx
//
// Unit tests of DetectorSimulator::expose(): simulate() plus the charge deposits of the source that
// initialize() creates from "cosmic" (_makeChargeDeposit()), added by _addDeposits().

#include <gtest/gtest.h>
#include <phot2pix/DetectorSimulator.h>
#include <phot2pix/CosmicDeposit.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <valarray>
#include <vector>

using namespace phot2pix;

namespace
{
    // Pixel without noise: its response is flux * t, so that every image is known exactly.
    class QuietPixel : public pixsim
    {
    public:
        double pixelresponse(const double& flux, const double& t) const override { return flux * t; }
        double pixelresponse(const double& flux, const double& t, const PixelConfig&) const override { return flux * t; }
    };

    // Source that deposits the same charges at every call, once a "cosmic" key has configured it.
    class FixedSource : public ChargeDeposit
    {
    public:
        explicit FixedSource(std::vector<Deposit> deposits) : m_deposits(std::move(deposits)) {}

        std::shared_ptr<ChargeDeposit> clone() const override { return std::make_shared<FixedSource>(*this); }
        void configure(const PixelConfig& cfg) override { m_configured = m_configured || cfg.has("cosmic"); }
        bool isConfigured() const override { return m_configured; }
        std::vector<HeaderCard> inputCards() const override { return {{"CRTEST", 1.0, "Test source"}}; }

        void generate(double, const Geometry&, std::vector<Deposit>& out) const override
        {
            out.insert(out.end(), m_deposits.begin(), m_deposits.end());
        }

    private:
        std::vector<Deposit> m_deposits;
        bool                 m_configured = false;
    };

    // Detector of quiet pixels whose charge source is a FixedSource.
    class FixedDetector : public DetectorSimulator<QuietPixel>
    {
    public:
        FixedDetector(size_t w, size_t h, std::vector<Deposit> deposits)
            : DetectorSimulator<QuietPixel>(w, h), m_deposits(std::move(deposits)) {}

    protected:
        std::shared_ptr<ChargeDeposit> _makeChargeDeposit() const override { return std::make_shared<FixedSource>(m_deposits); }

    private:
        std::vector<Deposit> m_deposits;
    };

    // Detector of quiet pixels without a charge source, as nisp::SCA and nisp::persistence::pixelArray.
    class NoSourceDetector : public DetectorSimulator<QuietPixel>
    {
    public:
        NoSourceDetector(size_t w, size_t h) : DetectorSimulator<QuietPixel>(w, h) {}

    protected:
        std::shared_ptr<ChargeDeposit> _makeChargeDeposit() const override { return nullptr; }
    };

    // Values of a double-precision image.
    const std::valarray<double>& values(const std::shared_ptr<DSL::FITScube>& img) { return *img->GetData<double>(); }

    // True if a card has this keyword.
    bool hasCard(const std::vector<HeaderCard>& cards, const std::string& key)
    {
        return std::any_of(cards.begin(), cards.end(), [&key](const HeaderCard& c) { return c.key == key; });
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  expose() with a fixed source
// ═══════════════════════════════════════════════════════════════════════════════

TEST(DetectorExposeTest, WithoutCosmicKeyExposeIsSimulate)
{
    const std::vector<Deposit> hits{{1 * 8 + 2, 100.0}};
    FixedDetector d(8, 4, hits);
    d.initialize(PixelConfig(json{{"gain", 2.0}}));
    EXPECT_EQ(d.getChargeDeposit(), nullptr);
    EXPECT_FALSE(hasCard(d.getInputCards(), "CRTEST"));

    const auto exposed   = d.expose(3.0, 2.0);
    const auto simulated = d.simulate(3.0, 2.0);
    ASSERT_EQ(values(exposed).size(), d.size());
    for (size_t i = 0; i < d.size(); ++i)
        EXPECT_DOUBLE_EQ(values(exposed)[i], values(simulated)[i]);
}

TEST(DetectorExposeTest, AddsTheDepositsDividedByTheGain)
{
    // Two deposits on matrix pixel (2, 1) and one on (5, 3); the matrix is the physical detector.
    const std::vector<Deposit> hits{{1 * 8 + 2, 100.0}, {1 * 8 + 2, 50.0}, {3 * 8 + 5, 10.0}};
    FixedDetector d(8, 4, hits);
    d.initialize(PixelConfig(json{{"gain", 2.0}, {"cosmic", json::object()}}));
    ASSERT_NE(d.getChargeDeposit(), nullptr);
    EXPECT_TRUE(hasCard(d.getInputCards(), "CRTEST"));

    const auto exposed   = d.expose(3.0, 2.0);   // 6 ADU per pixel, plus the deposits
    const auto simulated = d.simulate(3.0, 2.0);  // 6 ADU per pixel: simulate() adds no deposit
    for (size_t row = 0; row < 4; ++row)
        for (size_t col = 0; col < 8; ++col)
        {
            double expected = 6.0;
            if (col == 2 && row == 1) expected += (100.0 + 50.0) / 2.0;
            if (col == 5 && row == 3) expected += 10.0 / 2.0;
            EXPECT_DOUBLE_EQ(values(exposed)[row * 8 + col], expected) << "pixel (" << col << ", " << row << ")";
            EXPECT_DOUBLE_EQ(values(simulated)[row * 8 + col], 6.0);
        }

    // The overload with runtime parameters adds the same deposits
    const auto withCfg = d.expose(std::valarray<double>(3.0, d.size()), 2.0, PixelConfig());
    for (size_t i = 0; i < d.size(); ++i)
        EXPECT_DOUBLE_EQ(values(withCfg)[i], values(exposed)[i]);
}

TEST(DetectorExposeTest, AWindowKeepsOnlyItsDeposits)
{
    const std::vector<Deposit> hits{{1 * 8 + 2, 100.0}, {3 * 8 + 5, 10.0}};
    FixedDetector d(8, 4, hits);
    d.initialize(PixelConfig(json{{"gain", 2.0}, {"cosmic", json::object()}}));

    // Window of 3 × 2 pixels from (2, 1): it holds matrix pixel (2, 1), not (5, 3)
    const std::valarray<double> flux(3.0, 3 * 2);
    const auto img     = d.expose(flux, 2.0, 2, 1, 3, 2);
    const auto withCfg = d.expose(flux, 2.0, 2, 1, 3, 2, PixelConfig());
    ASSERT_EQ(values(img).size(), 6u);
    EXPECT_DOUBLE_EQ(values(img)[0], 6.0 + 100.0 / 2.0);
    EXPECT_DOUBLE_EQ(values(withCfg)[0], 6.0 + 100.0 / 2.0);
    for (size_t i = 1; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(values(img)[i], 6.0);
        EXPECT_DOUBLE_EQ(values(withCfg)[i], 6.0);
    }

    // The window is checked as by simulate()
    EXPECT_THROW(d.expose(flux, 2.0, 6, 3, 3, 2), std::out_of_range);
    EXPECT_THROW(d.expose(std::valarray<double>(3.0, 5), 2.0, 2, 1, 3, 2), std::invalid_argument);
}

TEST(DetectorExposeTest, UsesThePlacementOfTheMatrixOnTheDetector)
{
    // Matrix of 8 × 4 pixels at (4, 2) on a physical detector of 16 × 8 pixels
    const size_t W = 16;
    const std::vector<Deposit> hits{{(2 + 1) * W + (4 + 2), 100.0},   // matrix pixel (2, 1)
                                    {1 * W + 1, 30.0},                // left of the matrix and above it
                                    {7 * W + 15, 20.0}};              // right of the matrix and below it
    FixedDetector d(8, 4, hits);
    d.initialize(PixelConfig(json{{"detectorSize", {16, 8}}, {"detectorOrigin", {4, 2}},
                                  {"gain", 2.0}, {"cosmic", json::object()}}));

    const auto img = d.expose(3.0, 2.0);
    for (size_t row = 0; row < 4; ++row)
        for (size_t col = 0; col < 8; ++col)
        {
            const double expected = (col == 2 && row == 1) ? 6.0 + 100.0 / 2.0 : 6.0;
            EXPECT_DOUBLE_EQ(values(img)[row * 8 + col], expected) << "pixel (" << col << ", " << row << ")";
        }
}

TEST(DetectorExposeTest, NoSourceWhenTheFamilyReturnsNone)
{
    NoSourceDetector d(8, 4);
    d.initialize(PixelConfig(json{{"gain", 2.0}, {"cosmic", json::object()}}));
    EXPECT_EQ(d.getChargeDeposit(), nullptr);

    const auto exposed   = d.expose(3.0, 2.0);
    const auto simulated = d.simulate(3.0, 2.0);
    for (size_t i = 0; i < d.size(); ++i)
        EXPECT_DOUBLE_EQ(values(exposed)[i], values(simulated)[i]);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  The default source: CosmicRay::CosmicDeposit
// ═══════════════════════════════════════════════════════════════════════════════

TEST(DetectorExposeTest, TheDefaultSourceIsACosmicDeposit)
{
    // Protons with J = 1e9 (E/10)^-2 between 10 and 100 MeV: about 94 hits in 0.01 s on 32 × 32 pixels of 18 µm
    const std::filesystem::path csv = std::filesystem::temp_directory_path() / "p2p_gTest_detectorExpose.csv";
    std::ofstream(csv) << "E,J\n10,1e9\n100,1e7\n";

    json cosmic = {{"depth_um", 10.0}, {"pitch_um", 18.0}, {"eVPerElectron", 3.6}, {"seed", 7}};
    cosmic["spectrum"] = json::array({json{{"species", "p"}, {"Z", 1}, {"A", 1}, {"file", csv.string()},
                                           {"colEnergy", "E"}, {"colFlux", "J"}}});
    const PixelConfig cfg(json{{"gain", 2.0}, {"cosmic", cosmic}});

    DetectorSimulator<QuietPixel> d(32, 32);
    d.initialize(cfg);
    ASSERT_NE(d.getChargeDeposit(), nullptr);
    EXPECT_TRUE(d.getChargeDeposit()->isConfigured());
    EXPECT_NE(dynamic_cast<const CosmicRay::CosmicDeposit*>(d.getChargeDeposit()), nullptr);
    EXPECT_TRUE(hasCard(d.getInputCards(), "CRFLUX"));

    // A copy of the detector gets its own copy of the source, at the same point of the random sequence
    const DetectorSimulator<QuietPixel> copy(d);

    // The same source on its own: its first call draws the hits of the first expose()
    CosmicRay::CosmicDeposit reference;
    reference.configure(cfg);
    std::vector<Deposit> hits;
    reference.generate(0.01, d.getGeometry(), hits);
    ASSERT_FALSE(hits.empty());

    std::valarray<double> expected(0.0, d.size());
    for (const Deposit& h : hits)
        expected[h.idx] += h.electrons / 2.0;   // the matrix is the physical detector

    const auto img       = d.expose(0.0, 0.01);
    const auto imgOfCopy = copy.expose(0.0, 0.01);
    for (size_t i = 0; i < d.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(values(img)[i], expected[i]) << "pixel " << i;
        EXPECT_DOUBLE_EQ(values(imgOfCopy)[i], expected[i]) << "pixel " << i;
    }

    std::error_code ec;
    std::filesystem::remove(csv, ec);
}
