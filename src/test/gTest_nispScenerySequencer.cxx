//
// Created by gillard on 2026-11-22.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_nispScenerySequencer.cxx
//
// nisp::ScenerySequencer: the NISP keys of an entry (timing as "exposureTime" or MACC
// parameters, optional "exposureMode"), and run() on an SCA: each exposure in its own
// readout mode, the SCA's mode kept, the flux image cropped at the matrix origin.
// The base class (the "flux" array, outputs) is tested in gTest_scenerySequencer.cxx.
//

#include <gtest/gtest.h>
#include <nisp/nispScenerySequencer.h>
#include <DSTfits/FITSmanager.h>

#include <cmath>
#include <filesystem>
#include <string>
#include <valarray>
#include <vector>

using namespace phot2pix;

namespace
{
    // Write a W × H double map whose pixel (x, y) holds value(x, y) to test_data/<name> and
    // return its path. CFITSIO does not overwrite: a file left by a previous run is removed first.
    template<typename F>
    std::string writeTestMap(const std::string& name, size_t W, size_t H, F value)
    {
        const std::filesystem::path dir  = "test_data";
        const std::filesystem::path file = dir / name;

        std::filesystem::create_directories(dir);
        std::filesystem::remove(file);

        DSL::FITSimg<double> map({W, H});
        for (size_t y = 0; y < H; ++y)
            for (size_t x = 0; x < W; ++x)
                map.SetPixelValue(value(x, y), y * W + x);
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

    // Pixel values of an image, whatever its storage type.
    std::valarray<double> pixels(const DSL::FITScube& image)
    {
        if (image.GetEqBitPerPixel() == DOUBLE_IMG && image.GetBitPerPixel() == DOUBLE_IMG)
            return *image.GetData<double>();
        return *image.cast<double>()->GetData<double>();
    }

    PixelConfig sequence(const json& entries) { return PixelConfig(json{{"flux", entries}}); }

    // SCA with no noise but the shot noise (gain 1, QE 1), a 1 s frame time, in MACC mode.
    // The matrix is on the active area (origin [4, 4]) unless cfg says otherwise.
    void initializeQuiet(nisp::SCA& sca, PixelConfig cfg = PixelConfig{})
    {
        cfg.mergeDefaults(PixelConfig{
            {"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 1.0},
            {"KTCnoise", 0.0}, {"baseline", 0.0}, {"frame_time", 1.0},
            {"exposureMode", std::string("MACC")}, {"detectorOrigin", json::array({4, 4})}
        });
        sca.initialize(cfg);
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  NISP keys of an entry
// ═══════════════════════════════════════════════════════════════════════════════

TEST(NispScenerySequencer, InitializeAcceptsAnExposureTimeOrMaccParameters)
{
    // The files are not opened by initialize(): they need not exist.
    nisp::ScenerySequencer seq;
    EXPECT_NO_THROW(seq.initialize(sequence(json::array({
        {{"file", "a.fits"}, {"exposureTime", 10.0}},
        {{"file", "a.fits"}, {"macc", {4, 16, 4}}},
        {{"file", "a.fits"}, {"ng", 3}, {"nr", 1}, {"nd", 0}, {"exposureMode", "utr"}},
        {{"file", "a.fits"}, {"exposureTime", 10.0}, {"exposureMode", "CDS"}}
    }))));
    EXPECT_EQ(seq.size(), 4u);
}

TEST(NispScenerySequencer, InitializeRejectsInvalidNispKeys)
{
    const std::vector<json> bad{
        {{"file", "a.fits"}},                                                 // no timing
        {{"file", "a.fits"}, {"exposureTime", 10.0}, {"macc", {4, 16, 4}}},   // two timings
        {{"file", "a.fits"}, {"exposureTime", 10.0}, {"ng", 4}},              // two timings
        {{"file", "a.fits"}, {"macc", {0, 16, 4}}},                           // ng = 0
        {{"file", "a.fits"}, {"macc", {4, 16}}},                              // not [ng, nr, nd]
        {{"file", "a.fits"}, {"exposureTime", 10.0}, {"exposureMode", "FOO"}},
        {{"file", "a.fits"}, {"exposureTime", 10.0}, {"exposureMode", 3}},
    };

    for (const json& entry : bad)
    {
        nisp::ScenerySequencer seq;
        try
        {
            seq.initialize(sequence(json::array({entry})));
            ADD_FAILURE() << "accepted " << entry.dump();
        }
        catch (const std::invalid_argument& err)
        {
            EXPECT_NE(std::string(err.what()).find("'flux'[0]"), std::string::npos) << err.what();
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  run() on an SCA
// ═══════════════════════════════════════════════════════════════════════════════

TEST(NispScenerySequencer, RunBeforeInitializeThrows)
{
    nisp::SCA sca(4, 4);
    initializeQuiet(sca);
    nisp::ScenerySequencer seq;
    EXPECT_THROW(seq.run(sca), std::logic_error);
}

TEST(NispScenerySequencer, EachExposureUsesItsOwnReadoutMode)
{
    const size_t W = 16, H = 16;
    const std::string flux = writeTestMap("nseq_flux.fits", W, H, [](size_t, size_t) { return 10.0; });
    const std::string outF = outputPath("nseq_flux_mode.fits");
    const std::string outM = outputPath("nseq_macc_mode.fits");
    const std::string outC = outputPath("nseq_cds_mode.fits");

    nisp::SCA sca(W, H);
    initializeQuiet(sca);   // MACC by default

    nisp::ScenerySequencer seq;
    seq.initialize(sequence(json::array({
        {{"file", flux}, {"exposureMode", "FLUX"}, {"exposureTime", 10.0}, {"output", outF}},
        {{"file", flux}, {"macc", {3, 2, 1}},                              {"output", outM}},   // the SCA's mode
        {{"file", flux}, {"exposureMode", "CDS"},  {"exposureTime", 10.0}, {"output", outC}}
    })));
    seq.run(sca);

    const auto f = readImage(outF);
    ASSERT_TRUE(f);
    EXPECT_EQ(f->HDU().GetValueForKey("MODE"), "FLUX");

    const auto m = readImage(outM);
    ASSERT_TRUE(m);
    EXPECT_EQ(m->HDU().GetValueForKey("MODE"), "MACC");
    EXPECT_EQ(m->Size(3), 6u);   // every read is a layer: nGroups × nFrames
    EXPECT_NEAR(m->HDU().GetDoubleValueForKey("SEQIDX"), 2.0, 1e-9);

    const auto c = readImage(outC);
    ASSERT_TRUE(c);
    EXPECT_EQ(c->HDU().GetValueForKey("MODE"), "CDS");

    // The mode of an entry was for its exposure only: the SCA is still in MACC mode.
    const auto after = sca.expose(std::valarray<double>(10.0, sca.size()), nisp::macc{2, 1, 0});
    ASSERT_TRUE(after);
    EXPECT_EQ(after->HDU().GetValueForKey("MODE"), "MACC");
    EXPECT_EQ(after->Size(3), 2u);

    for (const std::string& path : {flux, outF, outM, outC})
        std::filesystem::remove(path);
}

TEST(NispScenerySequencer, DetectorSizedFluxIsCroppedAtTheOrigin)
{
    // 64 × 64 detector, 8 × 8 matrix at [24, 24]. The flux image is 100 (c + 1) ph/s on the
    // section of the matrix (c its column) and 0 elsewhere: a section taken at the wrong
    // place shows a wrong column level or a row of zeros.
    const size_t W = 8, H = 8, x0 = 24, y0 = 24;
    auto value = [=](size_t x, size_t y)
    {
        const bool inside = x >= x0 && x < x0 + W && y >= y0 && y < y0 + H;
        return inside ? 100.0 * static_cast<double>(x - x0 + 1) : 0.0;
    };
    const std::string flux = writeTestMap("nseq_flux_64x64.fits", 64, 64, value);
    const std::string out  = outputPath("nseq_cropped.fits");

    nisp::SCA sca(W, H);
    initializeQuiet(sca, PixelConfig{{"detectorSize",   json::array({64, 64})},
                                     {"detectorOrigin", json::array({x0, y0})}});

    nisp::ScenerySequencer seq;
    seq.initialize(sequence(json::array({
        {{"file", flux}, {"exposureMode", "FLUX"}, {"exposureTime", 10.0}, {"output", out}}
    })));
    seq.run(sca);

    const auto image = readImage(out);
    ASSERT_TRUE(image);
    ASSERT_EQ(image->Size(1), W);
    ASSERT_EQ(image->Size(2), H);

    // S = 1000 (c + 1) electrons = ADU (gain 1), shot noise only.
    const std::valarray<double> data = pixels(*image);
    ASSERT_EQ(data.size(), W * H);
    for (size_t r = 0; r < H; ++r)
        for (size_t col = 0; col < W; ++col)
        {
            const double S = 1000.0 * static_cast<double>(col + 1);
            EXPECT_NEAR(data[r * W + col], S, 6.0 * std::sqrt(S) + 1.0)
                << "column " << col << ", row " << r;
        }

    std::filesystem::remove(flux);
    std::filesystem::remove(out);
}
