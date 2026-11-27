//
// Created by gillard on 2026-11-22.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_scenerySequencer.cxx
//
// phot2pix::ScenerySequencer: the "flux" array (parsing and checks), the hook of a detector
// family (_checkExposure()), and run() on a plain DetectorSimulator: order, saved outputs
// and their header cards. The NISP sequencer is tested in gTest_nispScenerySequencer.cxx.
//

#include <gtest/gtest.h>
#include <phot2pix/ScenerySequencer.h>
#include <DSTfits/FITSmanager.h>

#include <cmath>
#include <filesystem>
#include <string>
#include <valarray>
#include <vector>

using namespace phot2pix;

namespace
{
    // Write a W × H double map holding value in every pixel to test_data/<name> and return its
    // path. CFITSIO does not overwrite: a file left by a previous run is removed first.
    std::string writeUniformMap(const std::string& name, size_t W, size_t H, double value)
    {
        const std::filesystem::path dir  = "test_data";
        const std::filesystem::path file = dir / name;

        std::filesystem::create_directories(dir);
        std::filesystem::remove(file);

        DSL::FITSimg<double> map({W, H});
        for (size_t i = 0; i < W * H; ++i)
            map.SetPixelValue(value, i);
        map.Write(file.string());
        return file.string();
    }

    // Path of an output in test_data/, the file left by a previous run removed.
    std::string outputPath(const std::string& name)
    {
        const std::filesystem::path dir  = "test_data";
        const std::filesystem::path file = dir / name;

        std::filesystem::create_directories(dir);
        std::filesystem::remove(file);
        return file.string();
    }

    // Primary HDU of a FITS file (the image keeps its data once the file is closed).
    std::shared_ptr<DSL::FITScube> readImage(const std::string& path)
    {
        DSL::FITSmanager file(path, true);
        return file.GetImageAtIndex(1);
    }

    // Mean pixel value of an image, whatever its storage type.
    double meanOf(const DSL::FITScube& image)
    {
        const std::valarray<double> data =
            (image.GetEqBitPerPixel() == DOUBLE_IMG && image.GetBitPerPixel() == DOUBLE_IMG)
                ? *image.GetData<double>()
                : *image.cast<double>()->GetData<double>();
        return data.sum() / static_cast<double>(data.size());
    }

    PixelConfig sequence(const json& entries) { return PixelConfig(json{{"flux", entries}}); }

    // No noise but the shot noise: gain 1, QE 1, no dark current, read noise, kTC or baseline.
    PixelConfig quietPixels()
    {
        return PixelConfig{{"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 1.0},
                           {"KTCnoise", 0.0}, {"baseline", 0.0}};
    }

    // A detector family whose entries give "nFrames" instead of "exposureTime".
    class FramesSequencer : public ScenerySequencer
    {
    protected:
        void _checkExposure(const Exposure& e) const override
        {
            if (!e.cfg.has("nFrames"))
                throw std::invalid_argument("no \"nFrames\".");
        }
    };

    // Charge source that deposits the same charge on one pixel at every call (a fixed cosmic-ray hit).
    class OneHitSource : public ChargeDeposit
    {
    public:
        OneHitSource(size_t idx, double electrons) : m_hit{idx, electrons} {}

        std::shared_ptr<ChargeDeposit> clone() const override { return std::make_shared<OneHitSource>(*this); }
        void configure(const PixelConfig&) override {}
        void generate(double, const Geometry&, std::vector<Deposit>& out) const override { out.push_back(m_hit); }

    private:
        Deposit m_hit;
    };

    // Detector whose charge source, created from the "cosmic" key, is a OneHitSource.
    class OneHitDetector : public DetectorSimulator<>
    {
    public:
        OneHitDetector(size_t w, size_t h, size_t idx, double electrons)
            : DetectorSimulator<>(w, h), m_idx(idx), m_electrons(electrons) {}

    protected:
        std::shared_ptr<ChargeDeposit> _makeChargeDeposit() const override { return std::make_shared<OneHitSource>(m_idx, m_electrons); }

    private:
        size_t m_idx;
        double m_electrons;
    };
}

// ═══════════════════════════════════════════════════════════════════════════════
//  The "flux" array
// ═══════════════════════════════════════════════════════════════════════════════

TEST(ScenerySequencer, InConfigOnlyForAnArray)
{
    EXPECT_TRUE (ScenerySequencer::inConfig(sequence(json::array({{{"file", "a.fits"}, {"exposureTime", 1.0}}}))));
    EXPECT_TRUE (ScenerySequencer::inConfig(sequence(json::array())));   // empty: rejected by initialize()

    EXPECT_FALSE(ScenerySequencer::inConfig(PixelConfig{}));
    EXPECT_FALSE(ScenerySequencer::inConfig(PixelConfig{{"flux", 10.0}}));
    EXPECT_FALSE(ScenerySequencer::inConfig(PixelConfig{{"flux", std::string("a.fits")}}));
    EXPECT_FALSE(ScenerySequencer::inConfig(PixelConfig(json{{"flux", {{"file", "a.fits"}, {"hdu", 2}}}})));
}

TEST(ScenerySequencer, InitializeReadsEveryEntry)
{
    // The files are not opened by initialize(): they need not exist.
    ScenerySequencer seq;
    seq.initialize(sequence(json::array({
        {{"file", "a.fits"}, {"hdu", 2},           {"exposureTime", 10.0}, {"output", "out_a.fits"}},
        {{"file", "b.fits"}, {"hdu", "DET13.SCI"}, {"exposureTime", 5},    {"custom", 3}},
        {{"file", "c.fits"},                       {"exposureTime", 2.5}}
    })));

    ASSERT_EQ(seq.size(), 3u);
    const auto& e = seq.getExposures();

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_EQ(e[i].index, i);
        EXPECT_FALSE(e[i].cfg.has("file"));     // moved into "flux"
        EXPECT_FALSE(e[i].cfg.has("hdu"));
        EXPECT_FALSE(e[i].cfg.has("output"));
    }

    EXPECT_EQ(e[0].cfg.getJson("flux"), json({{"file", "a.fits"}, {"hdu", 2}}));
    EXPECT_EQ(e[1].cfg.getJson("flux"), json({{"file", "b.fits"}, {"hdu", "DET13.SCI"}}));
    EXPECT_EQ(e[2].cfg.getJson("flux"), json({{"file", "c.fits"}}));

    EXPECT_DOUBLE_EQ(e[0].cfg.getAsDouble("exposureTime"), 10.0);
    EXPECT_DOUBLE_EQ(e[1].cfg.getAsDouble("exposureTime"), 5.0);   // a JSON integer
    EXPECT_DOUBLE_EQ(e[1].cfg.getAsDouble("custom"),       3.0);   // other keys are kept, for a derived class

    EXPECT_EQ(e[0].output, "out_a.fits");
    EXPECT_TRUE(e[1].output.empty());
    EXPECT_TRUE(e[2].output.empty());

    // The FLUX header card of each exposure.
    EXPECT_EQ(std::get<std::string>(*e[0].cfg.card("flux")), "a.fits[2]");
    EXPECT_EQ(std::get<std::string>(*e[1].cfg.card("flux")), "b.fits[DET13.SCI]");
    EXPECT_EQ(std::get<std::string>(*e[2].cfg.card("flux")), "c.fits");
}

TEST(ScenerySequencer, InitializeRejectsMalformedEntries)
{
    const json ok = {{"file", "a.fits"}, {"exposureTime", 1.0}};
    auto with = [&ok](const std::string& key, const json& value) { json e = ok; e[key] = value; return e; };
    auto without = [&ok](const std::string& key) { json e = ok; e.erase(key); return e; };

    const std::vector<PixelConfig> bad{
        PixelConfig{},                                         // no "flux"
        PixelConfig{{"flux", 10.0}},                           // not an array
        PixelConfig(json{{"flux", {{"file", "a.fits"}}}}),     // an object, not an array
        sequence(json::array()),                               // empty
        sequence(json::array({"a.fits"})),                     // an entry that is not an object
        sequence(json::array({without("file")})),
        sequence(json::array({with("file", "")})),
        sequence(json::array({with("file", 3)})),
        sequence(json::array({with("hdu", 0)})),
        sequence(json::array({with("hdu", -1)})),
        sequence(json::array({with("hdu", 1.5)})),
        sequence(json::array({with("hdu", "")})),
        sequence(json::array({with("flux", 10.0)})),                       // both "file" and "flux"
        sequence(json::array({with("flux", -1.0)})),
        sequence(json::array({json{{"flux", -1.0}, {"exposureTime", 1.0}}})),
        sequence(json::array({json{{"flux", "10"}, {"exposureTime", 1.0}}})),
        sequence(json::array({json{{"flux", 10.0}, {"hdu", 2}, {"exposureTime", 1.0}}})),
        sequence(json::array({without("exposureTime")})),     // required by the base class
        sequence(json::array({with("exposureTime", 0.0)})),
        sequence(json::array({with("exposureTime", -1.0)})),
        sequence(json::array({with("exposureTime", "10")})),
        sequence(json::array({with("output", "")})),
        sequence(json::array({with("output", 1)})),
        sequence(json::array({with("output", "same.fits"), with("output", "same.fits")})),
        sequence(json::array({with("output", "same.fits"), with("output", "./same.fits")})),
    };

    for (const PixelConfig& cfg : bad)
    {
        ScenerySequencer seq;
        EXPECT_THROW(seq.initialize(cfg), std::invalid_argument) << cfg.toJson().dump();
    }
}

TEST(ScenerySequencer, ErrorNamesTheEntry)
{
    ScenerySequencer seq;
    try
    {
        seq.initialize(sequence(json::array({
            {{"file", "a.fits"}, {"exposureTime", 1.0}},
            {{"file", "b.fits"}}
        })));
        FAIL() << "the second entry has no exposureTime";
    }
    catch (const std::invalid_argument& err)
    {
        EXPECT_NE(std::string(err.what()).find("'flux'[1]"), std::string::npos) << err.what();
    }
}

TEST(ScenerySequencer, FailedInitializeKeepsThePreviousSequence)
{
    ScenerySequencer seq;
    seq.initialize(sequence(json::array({{{"file", "a.fits"}, {"exposureTime", 1.0}},
                                         {{"file", "b.fits"}, {"exposureTime", 2.0}}})));
    ASSERT_EQ(seq.size(), 2u);

    EXPECT_THROW(seq.initialize(sequence(json::array({{{"file", "c.fits"}, {"exposureTime", 1.0}},
                                                      {{"file", "d.fits"}}}))),
                 std::invalid_argument);

    ASSERT_EQ(seq.size(), 2u);
    EXPECT_EQ(seq.getExposures()[1].cfg.getJson("flux"), json({{"file", "b.fits"}}));
}

TEST(ScenerySequencer, DerivedClassReplacesTheTimingCheck)
{
    FramesSequencer seq;
    EXPECT_NO_THROW(seq.initialize(sequence(json::array({{{"file", "a.fits"}, {"nFrames", 4}}}))));
    EXPECT_EQ(seq.size(), 1u);

    // Its own check replaces the base one: "exposureTime" alone is no longer enough.
    EXPECT_THROW(seq.initialize(sequence(json::array({{{"file", "a.fits"}, {"exposureTime", 1.0}}}))),
                 std::invalid_argument);
    // The base checks of the common keys still apply.
    EXPECT_THROW(seq.initialize(sequence(json::array({{{"file", "a.fits"}, {"nFrames", 4}, {"exposureTime", -1.0}}}))),
                 std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  run() on a plain DetectorSimulator
// ═══════════════════════════════════════════════════════════════════════════════

TEST(ScenerySequencer, RunBeforeInitializeThrows)
{
    const DetectorSimulator<> det(4, 4);
    ScenerySequencer seq;
    EXPECT_THROW(seq.run(det), std::logic_error);
}

TEST(ScenerySequencer, RunSavesTheExposuresThatHaveAnOutput)
{
    const size_t W = 16, H = 16;
    const std::string bright = writeUniformMap("seq_flux_100.fits", W, H, 100.0);
    const std::string faint  = writeUniformMap("seq_flux_50.fits",  W, H, 50.0);
    const std::string outA   = outputPath("seq_out_a.fits");
    const std::string outC   = outputPath("seq_out_c.fits");

    DetectorSimulator<> det(W, H);
    det.initialize(quietPixels());

    ScenerySequencer seq;
    seq.initialize(sequence(json::array({
        {{"file", bright}, {"exposureTime", 10.0}, {"output", outA}},
        {{"file", faint},  {"exposureTime", 2.0}},                     // simulated, not saved
        {{"file", faint},  {"hdu", 1}, {"exposureTime", 4.0}, {"output", outC}}
    })));
    seq.run(det);

    ASSERT_TRUE(std::filesystem::exists(outA));
    ASSERT_TRUE(std::filesystem::exists(outC));

    // Mean S = flux × t; shot noise only, so the mean of N pixels has variance S / N.
    const double n = static_cast<double>(W * H);

    const auto a = readImage(outA);
    ASSERT_TRUE(a);
    EXPECT_NEAR(meanOf(*a), 1000.0, 5.0 * std::sqrt(1000.0 / n));
    EXPECT_EQ  (a->HDU().GetValueForKey("FLUX"), "seq_flux_100.fits");
    EXPECT_NEAR(a->HDU().GetDoubleValueForKey("SEQIDX"), 1.0, 1e-9);
    EXPECT_NEAR(a->HDU().GetDoubleValueForKey("SEQNUM"), 3.0, 1e-9);

    const auto c = readImage(outC);
    ASSERT_TRUE(c);
    EXPECT_NEAR(meanOf(*c), 200.0, 5.0 * std::sqrt(200.0 / n));
    EXPECT_EQ  (c->HDU().GetValueForKey("FLUX"), "seq_flux_50.fits[1]");
    EXPECT_NEAR(c->HDU().GetDoubleValueForKey("SEQIDX"), 3.0, 1e-9);
    EXPECT_NEAR(c->HDU().GetDoubleValueForKey("SEQNUM"), 3.0, 1e-9);

    for (const std::string& path : {bright, faint, outA, outC})
        std::filesystem::remove(path);
}

TEST(ScenerySequencer, OutputsGoIntoTheOutputDirectory)
{
    const size_t W = 4, H = 4;
    const std::string           flux = writeUniformMap("seq_dir_flux.fits", W, H, 10.0);
    const std::filesystem::path dir  = "test_data/seq_output_dir";
    std::filesystem::remove_all(dir);   // created by run()

    DetectorSimulator<> det(W, H);
    det.initialize(quietPixels());

    ScenerySequencer seq;
    seq.initialize(sequence(json::array({
        {{"file", flux}, {"exposureTime", 1.0}, {"output", "a.fits"}},
        {{"file", flux}, {"exposureTime", 1.0}}                         // not saved
    })));
    seq.setOutputDirectory(dir.string());

    EXPECT_EQ(seq.outputPath(seq.getExposures()[0]), (dir / "a.fits").string());
    EXPECT_TRUE(seq.outputPath(seq.getExposures()[1]).empty());

    seq.run(det);
    EXPECT_TRUE(std::filesystem::exists(dir / "a.fits"));

    std::filesystem::remove_all(dir);
    std::filesystem::remove(flux);
}

TEST(ScenerySequencer, UniformFluxInsteadOfAFile)
{
    const size_t W = 16, H = 16;
    const std::string out  = outputPath("seq_uniform.fits");
    const std::string dark = outputPath("seq_dark.fits");

    ScenerySequencer seq;
    seq.initialize(sequence(json::array({
        {{"flux", 100.0}, {"exposureTime", 10.0}, {"output", out}},
        {{"flux", 0},     {"exposureTime", 5.0},  {"output", dark}}   // a JSON integer: a dark exposure
    })));

    ASSERT_EQ(seq.size(), 2u);
    const auto& e = seq.getExposures();
    EXPECT_DOUBLE_EQ(e[0].cfg.getAsDouble("flux"), 100.0);
    EXPECT_DOUBLE_EQ(std::get<double>(*e[0].cfg.card("flux")), 100.0);   // FLUX card: the value

    DetectorSimulator<> det(W, H);
    det.initialize(quietPixels());
    seq.run(det);

    const double n = static_cast<double>(W * H);

    const auto a = readImage(out);
    ASSERT_TRUE(a);
    EXPECT_NEAR(meanOf(*a), 1000.0, 5.0 * std::sqrt(1000.0 / n));
    EXPECT_NEAR(a->HDU().GetDoubleValueForKey("FLUX"), 100.0, 1e-9);

    const auto d = readImage(dark);
    ASSERT_TRUE(d);
    EXPECT_NEAR(meanOf(*d), 0.0, 1e-9);   // no flux, no dark current, no noise
    EXPECT_NEAR(d->HDU().GetDoubleValueForKey("FLUX"), 0.0, 1e-9);

    std::filesystem::remove(out);
    std::filesystem::remove(dark);
}

TEST(ScenerySequencer, RunAddsTheChargeDepositsOfTheDetector)
{
    // Dark exposures of quiet pixels (gain 1): each image is 0 except the hit of the source, 150 e- on pixel 5
    const size_t W = 4, H = 4;
    const std::string first  = outputPath("seq_hit_1.fits");
    const std::string second = outputPath("seq_hit_2.fits");

    OneHitDetector det(W, H, 5, 150.0);
    PixelConfig cfg = quietPixels();
    cfg.set("cosmic", json::object());   // creates the source of the detector
    det.initialize(cfg);
    ASSERT_NE(det.getChargeDeposit(), nullptr);

    ScenerySequencer seq;
    seq.initialize(sequence(json::array({
        {{"flux", 0.0}, {"exposureTime", 5.0}, {"output", first}},
        {{"flux", 0.0}, {"exposureTime", 5.0}, {"output", second}}
    })));
    seq.run(det);

    for (const std::string& path : {first, second})
    {
        const auto img = readImage(path);
        ASSERT_TRUE(img) << path;
        const std::valarray<double> data =
            (img->GetEqBitPerPixel() == DOUBLE_IMG && img->GetBitPerPixel() == DOUBLE_IMG)
                ? *img->GetData<double>()
                : *img->cast<double>()->GetData<double>();
        ASSERT_EQ(data.size(), W * H);
        for (size_t i = 0; i < W * H; ++i)
            EXPECT_NEAR(data[i], (i == 5) ? 150.0 : 0.0, 1e-9) << path << ", pixel " << i;
        std::filesystem::remove(path);
    }
}
