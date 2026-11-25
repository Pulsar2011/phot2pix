//
// Created by gillard on 2026-11-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_nispCrosstalk.cxx
//
// Crosstalk between the output channels of nisp::SCA:
// (a) kernel built from the configuration key "crosstalk", and invalid values;
// (b) ghosts, full frame and window;
// (c) baseline, gain and full well of the pixels;
// (d) exposures: MACC, FLUX, CDS and window.
//
// Input image: in every channel (64 columns), one column holds 1 in one row out of 64. This
// comb is shifted up by 10 rows from a channel to the next one. 10 c modulo 64 is different for
// the 32 channels: a row holds at most one source, so the ghosts of two sources never overlap.
// Kernel: 1 for the source (no loss), 0.5 for the other channels. Expected image: in the row of a
// source, 1 at its place and 0.5 at the same position of every other channel.
//
// To look at the images, give a directory in the environment variable NISP_IPC_OUTPUT (the same as
// for gTest_nispIPC): the input and the output are saved in it, as
// SCACrosstalk_<test>_<image>_input.fits and SCACrosstalk_<test>_<image>_output.fits. Nothing is
// saved when it is not set. From the build directory:
//   NISP_IPC_OUTPUT=crosstalk_images ./gTest_nispCrosstalk --gtest_filter='SCACrosstalk.*'
//

#include <gtest/gtest.h>
#include <nisp/nispSCA.h>
#include <DSTfits/FITSimg.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <valarray>
#include <variant>
#include <vector>

namespace
{
    using phot2pix::PixelConfig;

    constexpr size_t kN        = 2048;   // detector size
    constexpr size_t kChannels = 32;
    constexpr size_t kColumns  = 64;     // columns per channel
    constexpr size_t kStep     = 64;     // rows between two sources of a channel
    constexpr size_t kShift    = 10;     // rows between the combs of two consecutive channels
    constexpr size_t kPosition = 20;     // column of the sources in their channel (not its centre: a mirror is seen)

    // ── Kernels ──────────────────────────────────────────────────────────────────

    // 1 for the source, `coupling` for every other channel.
    std::vector<double> flatKernel(double coupling, size_t channels = kChannels)
    {
        std::vector<double> kernel(channels, coupling);
        kernel[0] = 1.0;
        return kernel;
    }

    // 1 for the source, coupling / d at a distance of d channels.
    std::vector<double> decreasingKernel(double coupling, size_t channels = kChannels)
    {
        std::vector<double> kernel(channels, 1.0);
        for (size_t d = 1; d < channels; ++d)
            kernel[d] = coupling / static_cast<double>(d);
        return kernel;
    }

    // ── Configurations ───────────────────────────────────────────────────────────

    // No noise, baseline 0, gain 1, no full well.
    PixelConfig quietCfg()
    {
        return PixelConfig{{"gain", 1.0}, {"dark", 0.0}, {"readNoise", 0.0}, {"QE", 1.0}, {"KTCnoise", 0.0},
                           {"baseline", 0.0}, {"frame_time", 1.0}};
    }

    // "crosstalk": [...]
    PixelConfig arrayCfg(const std::vector<double>& kernel)
    {
        PixelConfig cfg = quietCfg();
        cfg.set("crosstalk", nlohmann::json(kernel));
        return cfg;
    }

    // "crosstalk": {"kernel": [...]}, on an instrument of nChannels channels of colPerChannel columns.
    PixelConfig objectCfg(const std::vector<double>& kernel, int nChannels, int colPerChannel)
    {
        nlohmann::json crosstalk;
        crosstalk["kernel"] = kernel;

        PixelConfig cfg = quietCfg();
        cfg.set("crosstalk",     crosstalk);
        cfg.set("nChannels",     nChannels);
        cfg.set("colPerChannel", colPerChannel);
        return cfg;
    }

    PixelConfig jsonCfg(const std::string& text)
    {
        PixelConfig cfg = quietCfg();
        cfg.set("crosstalk", nlohmann::json::parse(text));
        return cfg;
    }

    // The SCA holds the crosstalk `kernel`, for channels of `columns` columns: one row, the
    // couplings at their distance, 0 elsewhere and in the centre (ghosts alone).
    void expectKernel(const nisp::SCA& sca, const std::vector<double>& kernel, size_t columns)
    {
        ASSERT_TRUE(sca.isCrosstalkEnabled());

        const auto&  k      = sca.getCrosstalkKernel();
        const size_t centre = (kernel.size() - 1) * columns;
        ASSERT_EQ(k.size(0), 2 * centre + 1) << "columns";
        ASSERT_EQ(k.size(1), 1u)             << "rows";

        std::valarray<double> expected(0.0, 2 * centre + 1);
        for (size_t d = 1; d < kernel.size(); ++d)
            expected[centre - d * columns] = expected[centre + d * columns] = kernel[d];
        for (size_t i = 0; i < expected.size(); ++i)
            ASSERT_DOUBLE_EQ(k.data()[i], expected[i]) << "weight " << i;
    }

    // ── Images ───────────────────────────────────────────────────────────────────

    // Row of the first source of a channel.
    size_t firstRow(size_t channel) { return (channel * kShift) % kStep; }

    // Input image of W × H pixels (see the top of the file).
    std::valarray<double> combs(size_t W, size_t H)
    {
        std::valarray<double> image(0.0, W * H);
        for (size_t ch = 0; ch * kColumns + kPosition < W; ++ch)
            for (size_t y = firstRow(ch); y < H; y += kStep)
                image[y * W + ch * kColumns + kPosition] = 1.0;
        return image;
    }

    struct Window { size_t c0, r0, ncols, nrows; };

    std::valarray<double> crop(const std::valarray<double>& image, size_t W, const Window& w)
    {
        return image[std::gslice(w.r0 * W + w.c0, std::valarray<size_t>{w.nrows, w.ncols}, std::valarray<size_t>{W, 1})];
    }

    // Expected image: every pixel keeps its signal, and adds kernel[d] of it to the pixels of the
    // same row that are d channels away, on each side.
    std::valarray<double> expectedCrosstalk(const std::valarray<double>& image, size_t ncols, const std::vector<double>& kernel,
                                            size_t columns = kColumns)
    {
        std::valarray<double> out = image;
        for (size_t i = 0; i < image.size(); ++i)
        {
            if (image[i] == 0.0)
                continue;
            const size_t x = i % ncols;
            for (size_t d = 1; d < kernel.size(); ++d)
            {
                if (x >= d * columns)
                    out[i - d * columns] += image[i] * kernel[d];
                if (x + d * columns < ncols)
                    out[i + d * columns] += image[i] * kernel[d];
            }
        }
        return out;
    }

    // Compare two images pixel by pixel: one failure with the first wrong pixel, not one per pixel.
    void expectSameImage(const std::valarray<double>& image, const std::valarray<double>& expected,
                         size_t ncols, double tolerance, const std::string& what)
    {
        ASSERT_EQ(image.size(), expected.size()) << what;

        size_t nWrong = 0, first = 0;
        for (size_t i = 0; i < image.size(); ++i)
            if (!(std::abs(image[i] - expected[i]) <= tolerance))   // also NaN
            {
                if (nWrong == 0)
                    first = i;
                ++nWrong;
            }

        EXPECT_EQ(nWrong, 0u) << what << ": first wrong pixel x = " << first % ncols << ", y = " << first / ncols
                              << ": " << image[first] << " instead of " << expected[first];
    }

    std::string label(const Window& w)
    {
        return "window_" + std::to_string(w.c0) + "_" + std::to_string(w.r0);
    }

    // Save the input and the output, images of ncols columns, in the directory
    // $NISP_IPC_OUTPUT. Nothing is done when it is not set.
    void saveImages(const std::valarray<double>& input, const std::valarray<double>& output,
                    size_t ncols, const std::string& name)
    {
        const char* directory = std::getenv("NISP_IPC_OUTPUT");
        if (!directory || !*directory)
            return;

        const testing::TestInfo*    test = testing::UnitTest::GetInstance()->current_test_info();
        const std::filesystem::path dir  = directory;
        const std::string           stem = std::string(test->test_suite_name()) + "_" + test->name() + "_" + name;
        std::filesystem::create_directories(dir);

        const auto save = [&](const std::valarray<double>& image, const std::string& suffix)
        {
            const std::filesystem::path file = dir / (stem + suffix);
            std::filesystem::remove(file);   // CFITSIO does not overwrite; no-op when absent

            DSL::FITSimg<double> fits(std::vector<size_t>{ncols, image.size() / ncols});
            *(fits.GetData<double>()) = image;
            fits.Write(file.string());
            std::cout << "[ XTALK    ] saved " << file.string() << std::endl;
        };
        save(input,  "_input.fits");
        save(output, "_output.fits");
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (a) Kernel read from the configuration
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCACrosstalkConfig, OffByDefault)
{
    nisp::SCA sca(kN, 4);
    sca.initialize(quietCfg());
    EXPECT_FALSE(sca.isCrosstalkEnabled());

    // Nothing is done then.
    const std::valarray<double> in = combs(kN, 4);
    std::valarray<double> image = in;
    sca.applyCrosstalk(image, 0, 0, kN, 4);
    expectSameImage(image, in, kN, 0.0, "crosstalk off");
}

TEST(SCACrosstalkConfig, ArrayUsesTheNispChannels)
{
    // 32 channels of 64 columns.
    nisp::SCA sca(kN, 4);
    sca.initialize(arrayCfg(decreasingKernel(0.5)));
    expectKernel(sca, decreasingKernel(0.5), 64);
}

TEST(SCACrosstalkConfig, NumberAndWidthOfTheChannelsAreConfigurable)
{
    nisp::SCA sca(kN, 4);
    sca.initialize(objectCfg(decreasingKernel(0.5, 32), 32, 64));
    expectKernel(sca, decreasingKernel(0.5, 32), 64);

    sca.initialize(objectCfg(decreasingKernel(0.5, 4), 4, 512));
    expectKernel(sca, decreasingKernel(0.5, 4), 512);

    // The channels change without a new "crosstalk" key: same couplings, at the new distances.
    PixelConfig channels = quietCfg();
    channels.set("colPerChannel", 256);
    sca.initialize(channels);
    expectKernel(sca, decreasingKernel(0.5, 4), 256);
}

TEST(SCACrosstalkConfig, KernelMayBeShorterThanTheNumberOfChannels)
{
    // [1, 0.5, 0.1]: 2 × 2 × 64 + 1 = 257 columns, 0.1 at ± 128, 0.5 at ± 64, 0 elsewhere.
    const std::vector<double> kernel = {1.0, 0.5, 0.1};
    nisp::SCA sca(kN, 4);
    sca.initialize(arrayCfg(kernel));
    expectKernel(sca, kernel, 64);
    EXPECT_EQ(sca.getCrosstalkKernel().size(), 257u);

    // The header gives the channels of the detector, and the couplings of the kernel.
    size_t nCouplings = 0;
    for (const auto& card : sca.getInputCards())
    {
        if (card.key == "XTALKNCH")
            EXPECT_EQ(std::get<size_t>(card.value), 32u);
        else if (card.key == "XTALKCOL")
            EXPECT_EQ(std::get<size_t>(card.value), 64u);
        else if (card.key.rfind("XTALK", 0) == 0)
            ++nCouplings;
    }
    EXPECT_EQ(nCouplings, 2u) << "XTALK1 and XTALK2";

    sca.initialize(objectCfg(kernel, 32, 64));
    expectKernel(sca, kernel, 64);
    sca.initialize(objectCfg(kernel, 4, 512));
    expectKernel(sca, kernel, 512);

    // The source alone: no crosstalk.
    sca.initialize(arrayCfg({1.0}));
    EXPECT_FALSE(sca.isCrosstalkEnabled());
}

TEST(SCACrosstalkConfig, OffAndAbsentKey)
{
    nisp::SCA sca(kN, 4);
    sca.initialize(arrayCfg(flatKernel(0.5)));
    expectKernel(sca, flatKernel(0.5), 64);

    sca.initialize(quietCfg());   // no "crosstalk": unchanged
    expectKernel(sca, flatKernel(0.5), 64);

    sca.initialize(jsonCfg("[]"));
    EXPECT_FALSE(sca.isCrosstalkEnabled());

    sca.initialize(arrayCfg(flatKernel(0.0)));   // no coupling
    EXPECT_FALSE(sca.isCrosstalkEnabled());
}

TEST(SCACrosstalkConfig, InvalidValuesThrowAndKeepThePreviousOne)
{
    const std::vector<double> previous = flatKernel(0.5);
    nisp::SCA sca(kN, 4);
    sca.initialize(arrayCfg(previous));

    auto expectThrow = [&sca, &previous](const PixelConfig& cfg, const std::string& what)
    {
        SCOPED_TRACE(what);
        EXPECT_THROW(sca.initialize(cfg), std::invalid_argument);
        expectKernel(sca, previous, 64);
    };

    // 0.5 at every distance but `bad` at a distance of d channels.
    auto with = [](size_t d, double bad)
    {
        std::vector<double> kernel = flatKernel(0.5);
        kernel[d] = bad;
        return arrayCfg(kernel);
    };

    expectThrow(with(0, 0.0), "the source must keep its signal: 0");
    expectThrow(with(0, 0.9), "the source must keep its signal: 0.9");
    expectThrow(with(3, -0.1), "negative coupling");
    expectThrow(with(3, std::numeric_limits<double>::quiet_NaN()), "NaN coupling");
    expectThrow(with(3, std::numeric_limits<double>::infinity()),  "infinite coupling");
    expectThrow(arrayCfg(flatKernel(0.5, 33)), "33 values for 32 channels");

    expectThrow(objectCfg(flatKernel(0.5, 32), 32, 32), "32 channels of 32 columns");
    expectThrow(objectCfg(flatKernel(0.5, 16), 16, 64), "16 channels of 64 columns");
    expectThrow(objectCfg(flatKernel(0.5, 32),  0, 64), "no channel");
    expectThrow(objectCfg(flatKernel(0.5, 32), 32,  0), "no column");
    expectThrow(objectCfg(flatKernel(0.5, 32), 32, -64), "negative columns");

    expectThrow(objectCfg(flatKernel(0.5, 5), 4, 512), "5 values for 4 channels");

    for (const char* text : {R"({})",
                             R"({"kernel": 0.5})",
                             R"({"kernel": []})",
                             R"({"kernel": [1, "a", 0.25]})",
                             R"({"kernel": [1, 0.5], "dynamic": true})",
                             R"({"kernel": [1, 0.5], "nChannels": 32})",
                             R"({"kernel": [1, 0.5], "colPerChannel": 64})"})
        expectThrow(jsonCfg(text), text);

    PixelConfig number = quietCfg();
    number.set("crosstalk", 0.5);
    expectThrow(number, "a number");

    PixelConfig path = quietCfg();
    path.set("crosstalk", std::string("kernel.fits"));
    expectThrow(path, "a file: there is no dynamic crosstalk");
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (b) Ghosts
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCACrosstalk, InputHoldsAtMostOneSourcePerRow)
{
    // Each of the 1024 even rows holds one source, the odd rows none.
    const std::valarray<double> image = combs(kN, kN);
    for (size_t y = 0; y < kN; ++y)
    {
        const std::valarray<double> row = image[std::slice(y * kN, kN, 1)];
        ASSERT_EQ(row.sum(), (y % 2 == 0) ? 1.0 : 0.0) << "row " << y;
    }

    // The comb of channel 1 is 10 rows above the one of channel 0.
    EXPECT_EQ(image[ 0 * kN + kPosition],            1.0);
    EXPECT_EQ(image[64 * kN + kPosition],            1.0);
    EXPECT_EQ(image[10 * kN + kColumns + kPosition], 1.0);
    EXPECT_EQ(image[74 * kN + kColumns + kPosition], 1.0);
}

TEST(SCACrosstalk, FullFrame)
{
    nisp::SCA sca(kN, kN);
    sca.initialize(arrayCfg(flatKernel(0.5)));

    const std::valarray<double> in = combs(kN, kN);
    std::valarray<double> image = in;
    sca.applyCrosstalk(image, 0, 0, kN, kN);
    saveImages(in, image, kN, "frame");

    // Source of channel 1, row 10: 1 at its place, 0.5 at the same position of the 31 others.
    for (size_t ch = 0; ch < kChannels; ++ch)
        EXPECT_NEAR(image[10 * kN + ch * kColumns + kPosition], (ch == 1) ? 1.0 : 0.5, 1e-12) << "channel " << ch;
    EXPECT_NEAR(image[10 * kN + kColumns + kPosition - 1], 0.0, 1e-12) << "next to the source";
    EXPECT_NEAR(image[10 * kN + kColumns + kPosition + 1], 0.0, 1e-12) << "next to the source";
    EXPECT_NEAR(image[10 * kN + kColumns - 1 - kPosition], 0.0, 1e-12) << "mirror position";
    EXPECT_NEAR(image[11 * kN + kPosition],                0.0, 1e-12) << "row without source";

    // Every pixel: in an even row, 1 for the channel that holds the source and 0.5 for the others.
    std::valarray<double> expected(0.0, kN * kN);
    for (size_t source = 0; source < kChannels; ++source)
        for (size_t y = firstRow(source); y < kN; y += kStep)
            for (size_t ch = 0; ch < kChannels; ++ch)
                expected[y * kN + ch * kColumns + kPosition] = (ch == source) ? 1.0 : 0.5;
    expectSameImage(image, expected, kN, 1e-12, "full frame");

    // No loss: 1024 sources of 1, each with 31 ghosts of 0.5.
    EXPECT_NEAR(image.sum(), 1024.0 * (1.0 + 31.0 * 0.5), 1e-6);
}

TEST(SCACrosstalk, CouplingDependsOnTheDistanceBetweenChannels)
{
    constexpr size_t H = 128;
    const std::vector<double> kernel = decreasingKernel(0.5);
    nisp::SCA sca(kN, H);
    sca.initialize(arrayCfg(kernel));

    const std::valarray<double> in = combs(kN, H);
    std::valarray<double> image = in;
    sca.applyCrosstalk(image, 0, 0, kN, H);
    saveImages(in, image, kN, "frame");

    // Source of channel 5, row 50: the same ghost on both sides.
    ASSERT_EQ(firstRow(5), 50u);
    EXPECT_NEAR(image[50 * kN + 5 * kColumns + kPosition], 1.0,        1e-12);
    EXPECT_NEAR(image[50 * kN + 4 * kColumns + kPosition], 0.5,        1e-12);
    EXPECT_NEAR(image[50 * kN + 6 * kColumns + kPosition], 0.5,        1e-12);
    EXPECT_NEAR(image[50 * kN + 3 * kColumns + kPosition], 0.5 / 2.0,  1e-12);
    EXPECT_NEAR(image[50 * kN + 7 * kColumns + kPosition], 0.5 / 2.0,  1e-12);
    EXPECT_NEAR(image[50 * kN + 31 * kColumns + kPosition], 0.5 / 26.0, 1e-12);

    expectSameImage(image, expectedCrosstalk(in, kN, kernel), kN, 1e-12, "decreasing kernel");
}

TEST(SCACrosstalk, ShortKernelOnlyReachesTheNearestChannels)
{
    // [1, 0.5, 0.1]: ghosts up to a distance of 2 channels.
    constexpr size_t H = 128;
    const std::vector<double> kernel = {1.0, 0.5, 0.1};
    nisp::SCA sca(kN, H);
    sca.initialize(arrayCfg(kernel));

    const std::valarray<double> in = combs(kN, H);
    std::valarray<double> image = in;
    sca.applyCrosstalk(image, 0, 0, kN, H);
    saveImages(in, image, kN, "frame");

    // Source of channel 5, row 50.
    ASSERT_EQ(firstRow(5), 50u);
    const std::vector<double> expected = {0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.5, 0.1, 0.0};   // channels 0 to 8
    for (size_t ch = 0; ch < expected.size(); ++ch)
        EXPECT_NEAR(image[50 * kN + ch * kColumns + kPosition], expected[ch], 1e-12) << "channel " << ch;

    expectSameImage(image, expectedCrosstalk(in, kN, kernel), kN, 1e-12, "short kernel");

    // The same as the kernel of 32 values completed with 0.
    std::vector<double> full = kernel;
    full.resize(kChannels, 0.0);
    sca.initialize(arrayCfg(full));
    std::valarray<double> completed = in;
    sca.applyCrosstalk(completed, 0, 0, kN, H);
    expectSameImage(completed, image, kN, 1e-15, "completed with 0");
}

TEST(SCACrosstalk, GaussianSource)
{
    // One large Gaussian source (sigma 15 pixels: 6 sigma are wider than a channel of 64 columns).
    // Its ghosts are copies of it, 0.5 and 0.1 times as bright, 1 and 2 channels away on each side:
    // 64 and 128 columns, or 128 and 256 columns with channels of 128 columns.
    constexpr double amplitude = 1000.0, x0 = 830.3, y0 = 1210.6, sigma = 15.0;
    const auto gaussian = [](double x, double y)
    {
        return amplitude * std::exp(-0.5 * ((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (sigma * sigma));
    };
    const std::vector<double> kernel = {1.0, 0.5, 0.1};

    for (size_t columns : {64, 128})
    {
        const std::string what = std::to_string(columns) + "cols";
        PixelConfig cfg = arrayCfg(kernel);
        cfg.set("colPerChannel", static_cast<int>(columns));

        nisp::SCA sca(kN, kN);
        sca.initialize(cfg);

        // Expected: the source and its shifted copies, from the formula and not from the image.
        std::valarray<double> in(0.0, kN * kN), expected(0.0, kN * kN);
        for (size_t y = 0; y < kN; ++y)
            for (size_t x = 0; x < kN; ++x)
            {
                const double X = static_cast<double>(x), Y = static_cast<double>(y);
                in[y * kN + x] = expected[y * kN + x] = gaussian(X, Y);
                for (size_t d = 1; d < kernel.size(); ++d)
                {
                    const double shift = static_cast<double>(d * columns);
                    expected[y * kN + x] += kernel[d] * (gaussian(X - shift, Y) + gaussian(X + shift, Y));
                }
            }

        std::valarray<double> image = in;
        sca.applyCrosstalk(image, 0, 0, kN, kN);
        saveImages(in, image, kN, "frame_" + what);

        expectSameImage(image, expected, kN, 1e-9, what);
        EXPECT_NEAR(image.sum() / in.sum(), 1.0 + 2.0 * 0.5 + 2.0 * 0.1, 1e-9) << what << ": no loss, 4 ghosts";
    }
}

TEST(SCACrosstalk, OtherChannelWidth)
{
    // 8 channels of 256 columns: ghosts every 256 columns.
    constexpr size_t H = 128;
    const std::vector<double> kernel = decreasingKernel(0.5, 8);
    nisp::SCA sca(kN, H);
    sca.initialize(objectCfg(kernel, 8, 256));

    const std::valarray<double> in = combs(kN, H);
    std::valarray<double> image = in;
    sca.applyCrosstalk(image, 0, 0, kN, H);
    saveImages(in, image, kN, "frame");
    expectSameImage(image, expectedCrosstalk(in, kN, kernel, 256), kN, 1e-12, "8 channels");
}

TEST(SCACrosstalk, Window)
{
    constexpr size_t H = 256;
    const std::vector<double> kernel = decreasingKernel(0.5);
    nisp::SCA sca(kN, H);
    sca.initialize(arrayCfg(kernel));
    const std::valarray<double> frame = combs(kN, H);

    // Whole channels; from the middle of a channel to the middle of another one; last pixels.
    for (const Window& w : {Window{4 * kColumns, 64, 6 * kColumns, 128}, Window{100, 37, 500, 90}, Window{kN - 300, H - 70, 300, 70}})
    {
        const std::valarray<double> in = crop(frame, kN, w);
        std::valarray<double> image = in;
        sca.applyCrosstalk(image, w.c0, w.r0, w.ncols, w.nrows);
        saveImages(in, image, w.ncols, label(w));
        ASSERT_GT(in.sum(), 0.0) << label(w);
        expectSameImage(image, expectedCrosstalk(in, w.ncols, kernel), w.ncols, 1e-12, label(w));
    }

    // At most one channel wide: no ghost.
    for (const Window& w : {Window{3 * kColumns, 0, kColumns, H}, Window{100, 0, kColumns, H}, Window{100, 0, 10, H}})
    {
        const std::valarray<double> in = crop(frame, kN, w);
        std::valarray<double> image = in;
        sca.applyCrosstalk(image, w.c0, w.r0, w.ncols, w.nrows);
        saveImages(in, image, w.ncols, label(w) + "_" + std::to_string(w.ncols) + "cols");
        expectSameImage(image, in, w.ncols, 0.0, label(w));
    }
}

TEST(SCACrosstalk, MatrixOnTheActiveArea)
{
    // 2040 columns at detector column 4: the first channel has 60 columns. The ghosts stay 64
    // columns away from their source.
    constexpr size_t W = 2040, H = 128;
    const std::vector<double> kernel = decreasingKernel(0.5);
    PixelConfig cfg = arrayCfg(kernel);
    cfg.set("detectorOrigin", nlohmann::json::array({4, 4}));

    nisp::SCA sca(W, H);
    sca.initialize(cfg);

    const std::valarray<double> in = combs(W, H);
    std::valarray<double> image = in;
    sca.applyCrosstalk(image, 0, 0, W, H);
    saveImages(in, image, W, "frame");
    expectSameImage(image, expectedCrosstalk(in, W, kernel), W, 1e-12, "active area");
}

TEST(SCACrosstalk, InvalidWindowThrowsAndKeepsTheImage)
{
    constexpr size_t H = 64;
    nisp::SCA sca(kN, H);
    sca.initialize(arrayCfg(flatKernel(0.5)));

    const std::valarray<double> in = combs(kN, H);
    std::valarray<double> image = in;
    EXPECT_THROW(sca.applyCrosstalk(image, 1, 0, kN, H), std::out_of_range);
    EXPECT_THROW(sca.applyCrosstalk(image, 0, 1, kN, H), std::out_of_range);
    expectSameImage(image, in, kN, 0.0, "unchanged");

    std::valarray<double> tooSmall(0.0, kN * H - 1);
    EXPECT_THROW(sca.applyCrosstalk(tooSmall, 0, 0, kN, H), std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (c) Baseline, gain and full well
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCACrosstalk, BaselineGivesNoGhost)
{
    // Baseline of 1000 ADU everywhere, gain 2 e-/ADU: the same ghosts, above the baseline.
    constexpr size_t H = 128;
    const std::vector<double> kernel = decreasingKernel(0.5);
    PixelConfig cfg = arrayCfg(kernel);
    cfg.set("baseline", 1000.0);
    cfg.set("gain",     2.0);

    nisp::SCA sca(kN, H);
    sca.initialize(cfg);

    const std::valarray<double> sources = combs(kN, H);
    const std::valarray<double> in      = 1000.0 + sources;
    std::valarray<double>       image   = in;
    sca.applyCrosstalk(image, 0, 0, kN, H);
    saveImages(in, image, kN, "frame");

    const std::valarray<double> expected = 1000.0 + expectedCrosstalk(sources, kN, kernel);
    expectSameImage(image, expected, kN, 1e-9, "baseline");
}

TEST(SCACrosstalk, GhostsAreLimitedByTheFullWellOfTheSource)
{
    // Full well of 100 e-, gain 2 e-/ADU: a source of 1e6 ADU outputs 50 ADU, its ghosts are 25 ADU.
    // The charge of the source is not changed here: pixel::read() limits it.
    constexpr size_t H = 128;
    PixelConfig cfg = arrayCfg(flatKernel(0.5));
    cfg.set("gain",     2.0);
    cfg.set("fullWell", 100.0);

    nisp::SCA sca(kN, H);
    sca.initialize(cfg);

    const std::valarray<double> sources = combs(kN, H);
    const std::valarray<double> in      = 1.0e6 * sources;
    std::valarray<double>       image   = in;
    sca.applyCrosstalk(image, 0, 0, kN, H);
    saveImages(in, image, kN, "frame");

    const std::valarray<double> ghosts   = expectedCrosstalk(sources, kN, flatKernel(0.5)) - sources;   // 0.5
    const std::valarray<double> expected = 1.0e6 * sources + 50.0 * ghosts;
    expectSameImage(image, expected, kN, 1e-6, "full well");
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (d) Exposures
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCACrosstalk, EveryReadOfAnExposureHasTheGhosts)
{
    // Sources far above a full well of 1000 e- from the first frame (no shot noise left), no flux
    // elsewhere, baseline 100 ADU, gain 1: every read holds 100 + 1000 at the sources, 100 + 500
    // at their ghosts and 100 elsewhere. The ghosts do not build up from a read to the next one.
    constexpr size_t H = 128;
    PixelConfig cfg = arrayCfg(flatKernel(0.5));
    cfg.set("baseline",  100.0);
    cfg.set("fullWell",  1000.0);
    cfg.set("refBorder", 0);   // every pixel is active

    nisp::SCA sca(kN, H);
    sca.initialize(cfg);

    const std::valarray<double> sources = combs(kN, H);
    const std::valarray<double> flux    = 1.0e6 * sources;

    // Saved images: the flux (photons/s) and the read (ADU).
    auto expectRead = [](const std::shared_ptr<DSL::FITScube>& result, size_t layer, const std::valarray<double>& expected,
                         const std::valarray<double>& input, size_t ncols, const std::string& what)
    {
        ASSERT_NE(result, nullptr) << what;
        const std::valarray<uint16_t>* data = result->GetData<uint16_t>();
        ASSERT_NE(data, nullptr) << what;
        ASSERT_GE(data->size(), (layer + 1) * expected.size()) << what;

        std::valarray<double> image(0.0, expected.size());
        for (size_t i = 0; i < image.size(); ++i)
            image[i] = static_cast<double>((*data)[layer * expected.size() + i]);
        saveImages(input, image, ncols, what);
        expectSameImage(image, expected, ncols, 0.0, what);
    };

    const std::valarray<double> expected = 100.0 + 1000.0 * expectedCrosstalk(sources, kN, flatKernel(0.5));

    const auto ramp = sca.expose(flux, nisp::macc{2, 2, 1});
    for (size_t layer = 0; layer < 4; ++layer)
        expectRead(ramp, layer, expected, flux, kN, "MACC_read_" + std::to_string(layer));

    const auto single = sca.expose(flux, nisp::readout{nisp::mode::FLUX, std::nullopt, 3.0});
    expectRead(single, 0, expected, flux, kN, "FLUX");

    // Read N − read 1: the sources and their ghosts are the same in both.
    const auto cds = sca.expose(flux, nisp::readout{nisp::mode::CDS, std::nullopt, 3.0});
    expectRead(cds, 0, std::valarray<double>(0.0, kN * H), flux, kN, "CDS");

    // Window: ghosts of its own sources only.
    const Window w{100, 37, 500, 90};
    const auto window = sca.expose(flux, 3.0, w.c0, w.r0, w.ncols, w.nrows);
    const std::valarray<double> expectedWindow = 100.0 + 1000.0 * expectedCrosstalk(crop(sources, kN, w), w.ncols, flatKernel(0.5));
    expectRead(window, 0, expectedWindow, crop(flux, kN, w), w.ncols, label(w));

    // Output header.
    ASSERT_NE(ramp, nullptr);
    EXPECT_NEAR(ramp->HDU().GetDoubleValueForKey("XTALKCOL"), 64.0, 1e-9);
    EXPECT_NEAR(ramp->HDU().GetDoubleValueForKey("XTALKNCH"), 32.0, 1e-9);
    EXPECT_NEAR(ramp->HDU().GetDoubleValueForKey("XTALK1"),    0.5, 1e-9);
    EXPECT_NEAR(ramp->HDU().GetDoubleValueForKey("XTALK31"),   0.5, 1e-9);
}
