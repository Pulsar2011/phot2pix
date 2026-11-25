//
// Created by gillard on 2026-11-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_nispIPC.cxx
//
// Inter-pixel capacitance (IPC) of nisp::SCA:
// (a) kernel read from the configuration key "ipc": square, rectangular, from a FITS cube
//     (sums of the layers below and above 1), off, and invalid values;
// (b) convolution by a static kernel, full frame and window;
// (c) convolution by the kernels of a FITS cube, full frame and window.
//
// The convolution is checked on a 2048 x 2048 image holding 1 in one pixel out of 64 along
// each axis and 0 elsewhere: each of these pixels gives an image of its own kernel.
//
// To look at the images, give a directory in the environment variable NISP_IPC_OUTPUT: the
// input and the output of every convolution are saved in it, as <test>_<image>_input.fits and
// <test>_<image>_output.fits. Nothing is saved when it is not set. From the build directory:
//   NISP_IPC_OUTPUT=ipc_images ./gTest_nispIPC --gtest_filter='SCAIpcConvolution.*'
//

#include <gtest/gtest.h>
#include <nisp/nispSCA.h>
#include <DSTfits/FITSimg.h>
#include <DSTfits/FITSmanager.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <valarray>
#include <vector>

namespace
{
    using phot2pix::PixelConfig;

    // [[0, 0.1, 0], [0.1, 0.6, 0.1], [0, 0.1, 0]], row by row; sum = 1.
    const std::vector<double> kSimple = {0.025, 0.1, 0.025,
                                         0.1, 0.5, 0.1,
                                         0.025, 0.1, 0.025};

    // ── Configurations ───────────────────────────────────────────────────────────

    // "ipc": [...]
    PixelConfig squareCfg(const std::vector<double>& weights)
    {
        PixelConfig cfg;
        cfg.set("ipc", nlohmann::json(weights));
        return cfg;
    }

    // "ipc": {"kernel": [...], "width": w, "height": h}
    PixelConfig rectangularCfg(const std::vector<double>& weights, int w, int h)
    {
        nlohmann::json ipc;
        ipc["kernel"] = weights;
        ipc["width"]  = w;
        ipc["height"] = h;

        PixelConfig cfg;
        cfg.set("ipc", ipc);
        return cfg;
    }

    // "ipc": "<path>"
    PixelConfig fileCfg(const std::string& path)
    {
        PixelConfig cfg;
        cfg.set("ipc", path);
        return cfg;
    }

    // ── FITS kernels ─────────────────────────────────────────────────────────────

    // Path test_data/<name>, without the file a previous run may have left.
    std::string testFile(const std::string& name)
    {
        const std::filesystem::path dir  = "test_data";
        const std::filesystem::path file = dir / name;

        std::filesystem::create_directories(dir);
        std::filesystem::remove(file);   // no-op when absent
        return file.string();
    }

    // Cube of W × H × nK values of type T: value(x, y, k).
    template<typename T, typename F>
    std::unique_ptr<DSL::FITSimg<T>> makeCube(size_t W, size_t H, size_t nK, F value)
    {
        auto cube = std::make_unique<DSL::FITSimg<T>>(std::vector<size_t>{W, H, nK});
        std::valarray<T>& data = *(cube->template GetData<T>());
        for (size_t k = 0; k < nK; ++k)
            for (size_t y = 0; y < H; ++y)
                for (size_t x = 0; x < W; ++x)
                    data[(k * H + y) * W + x] = static_cast<T>(value(x, y, k));
        return cube;
    }

    // Write that cube to test_data/<name> and return its path.
    template<typename T, typename F>
    std::string writeCube(const std::string& name, size_t W, size_t H, size_t nK, F value)
    {
        const std::string path = testFile(name);
        makeCube<T>(W, H, nK, value)->Write(path);
        return path;
    }

    // Removes a file at the end of the test, also when an ASSERT stops it.
    struct FileRemover
    {
        std::string path;
        ~FileRemover()
        {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }
    };

    // Kernel of pixel (x, y), different for each of its neighbours and from pixel to pixel: the
    // right, left, upper and lower neighbours receive a, a/2, 2a and a/4. Sum = 1.
    double pixelWeight(size_t x, size_t y, size_t k)
    {
        const double a = 0.01 * static_cast<double>(1 + (x + 3 * y) % 7);
        switch (k)
        {
            case 1:  return 0.25 * a;          // ( 0, −1)
            case 3:  return 0.50 * a;          // (−1,  0)
            case 4:  return 1.0 - 3.75 * a;    // kept
            case 5:  return a;                 // (+1,  0)
            case 7:  return 2.00 * a;          // ( 0, +1)
            default: return 0.0;
        }
    }

    // Number of weights of the configured cube that differ from expected(x, y, k).
    template<typename F>
    size_t countWrongWeights(const nisp::SCA& sca, size_t nK, F expected)
    {
        const std::shared_ptr<const DSL::FITScube> cube = sca.getIPCKernelCube();
        const std::valarray<float>*                data = cube ? cube->GetData<float>() : nullptr;   // stored as float
        const size_t W = sca.getWidth(), H = sca.getHeight();
        if (!data || data->size() != W * H * nK)
            return std::numeric_limits<size_t>::max();

        size_t nWrong = 0;
        for (size_t k = 0; k < nK; ++k)
            for (size_t y = 0; y < H; ++y)
                for (size_t x = 0; x < W; ++x)
                    if (!(std::abs(static_cast<double>((*data)[(k * H + y) * W + x]) - expected(x, y, k)) <= 1e-6))
                        ++nWrong;
        return nWrong;
    }

    // The SCA holds no kernel.
    void expectNoKernel(const nisp::SCA& sca)
    {
        EXPECT_FALSE(sca.isIPCEnabled());
        EXPECT_FALSE(sca.isIPCDynamic());
        EXPECT_FALSE(sca.getIPCKernelCube());
    }

    // The SCA holds the static kernel `weights` of w × h values.
    void expectStaticKernel(const nisp::SCA& sca, const std::vector<double>& weights, size_t w, size_t h)
    {
        EXPECT_TRUE (sca.isIPCEnabled());
        EXPECT_FALSE(sca.isIPCDynamic());
        EXPECT_FALSE(sca.getIPCKernelCube());

        const auto& kernel = sca.getIPCKernel();
        ASSERT_EQ(kernel.sizes().size(), 2u);
        EXPECT_EQ(kernel.size(0), w) << "columns";
        EXPECT_EQ(kernel.size(1), h) << "rows";
        ASSERT_EQ(kernel.size(), weights.size());
        for (size_t k = 0; k < weights.size(); ++k)
            EXPECT_NEAR(kernel.data()[k], weights[k], 1e-12) << "weight " << k;
    }

    // The SCA holds the kernels expected(x, y, k) of a cube of 9 layers.
    template<typename F>
    void expectDynamicKernel(const nisp::SCA& sca, F expected, const std::string& what = "")
    {
        EXPECT_TRUE(sca.isIPCEnabled()) << what;
        EXPECT_TRUE(sca.isIPCDynamic()) << what;
        EXPECT_EQ(countWrongWeights(sca, 9, expected), 0u) << what;
    }

    // ── Images ───────────────────────────────────────────────────────────────────

    constexpr size_t kN    = 2048;   // image size
    constexpr size_t kStep = 64;     // one pixel out of 64 holds 1 …
    constexpr size_t kHalf = 32;     // … from this one: 32, 96, 160, …, never on the border

    std::valarray<double> pointSources(size_t W, size_t H)
    {
        std::valarray<double> image(0.0, W * H);
        for (size_t y = kHalf; y < H; y += kStep)
            for (size_t x = kHalf; x < W; x += kStep)
                image[y * W + x] = 1.0;
        return image;
    }

    struct Window { size_t c0, r0, ncols, nrows; };

    // Its sources are inside, far from its border.
    constexpr Window kInside{64, 128, 256, 192};
    // Its first source is its corner pixel (0, 0), the next ones its pixels (64, 0), (0, 64), …
    constexpr Window kOnSource{96, 160, 200, 100};

    // Section of an image of W columns.
    std::valarray<double> crop(const std::valarray<double>& image, size_t W, const Window& w)
    {
        return image[std::gslice(w.r0 * W + w.c0, std::valarray<size_t>{w.nrows, w.ncols}, std::valarray<size_t>{W, 1})];
    }

    // Expected convolution of an image of ncols × nrows by kernels of kw × kh weights,
    // weight(x, y, k) being weight k of the kernel of pixel (x, y) of the image. Each pixel
    // distributes its charge with its own kernel; a pixel of the border, which misses neighbours,
    // is divided by the sum of its weights for the neighbours it has.
    template<typename Weight>
    std::valarray<double> expectedIpc(const std::valarray<double>& image, size_t ncols, size_t nrows,
                                      size_t kw, size_t kh, Weight weight)
    {
        using idx = std::ptrdiff_t;
        const idx W  = static_cast<idx>(ncols), H  = static_cast<idx>(nrows);
        const idx Kw = static_cast<idx>(kw),    Kh = static_cast<idx>(kh);
        const auto inside = [W, H](idx x, idx y) { return x >= 0 && x < W && y >= 0 && y < H; };

        std::valarray<double> out(0.0, ncols * nrows);
        for (idx y = 0; y < H; ++y)
            for (idx x = 0; x < W; ++x)
            {
                const double charge = image[static_cast<size_t>(y * W + x)];
                if (charge == 0.0)
                    continue;
                for (idx k = 0; k < Kw * Kh; ++k)
                {
                    const idx xp = x + k % Kw - Kw / 2, yp = y + k / Kw - Kh / 2;
                    if (inside(xp, yp))
                        out[static_cast<size_t>(yp * W + xp)] += charge * weight(static_cast<size_t>(x), static_cast<size_t>(y), static_cast<size_t>(k));
                }
            }

        for (idx y = 0; y < H; ++y)
            for (idx x = 0; x < W; ++x)
            {
                if (inside(x - Kw / 2, y - Kh / 2) && inside(x + Kw / 2, y + Kh / 2))
                    continue;
                double norm = 0.0;
                for (idx k = 0; k < Kw * Kh; ++k)
                    if (inside(x - (k % Kw - Kw / 2), y - (k / Kw - Kh / 2)))
                        norm += weight(static_cast<size_t>(x), static_cast<size_t>(y), static_cast<size_t>(k));
                out[static_cast<size_t>(y * W + x)] /= norm;
            }
        return out;
    }

    // The same for a static kernel.
    std::valarray<double> expectedIpc(const std::valarray<double>& image, size_t ncols, size_t nrows,
                                      const std::vector<double>& kernel, size_t kw, size_t kh)
    {
        return expectedIpc(image, ncols, nrows, kw, kh, [&kernel](size_t, size_t, size_t k) { return kernel[k]; });
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

    std::string name(const Window& w)
    {
        return "window at (" + std::to_string(w.c0) + ", " + std::to_string(w.r0) + ")";
    }

    // The same for a file name.
    std::string fileLabel(const Window& w)
    {
        return "window_" + std::to_string(w.c0) + "_" + std::to_string(w.r0);
    }

    // Save the input and the output of a convolution, images of ncols columns, in the directory
    // $NISP_IPC_OUTPUT. Nothing is done when it is not set.
    void saveImages(const std::valarray<double>& input, const std::valarray<double>& output,
                    size_t ncols, const std::string& label)
    {
        const char* directory = std::getenv("NISP_IPC_OUTPUT");
        if (!directory || !*directory)
            return;

        const std::filesystem::path dir  = directory;
        const std::string           stem = std::string(testing::UnitTest::GetInstance()->current_test_info()->name()) + "_" + label;
        std::filesystem::create_directories(dir);

        const auto save = [&](const std::valarray<double>& image, const std::string& suffix)
        {
            const std::filesystem::path file = dir / (stem + suffix);
            std::filesystem::remove(file);   // CFITSIO does not overwrite; no-op when absent

            DSL::FITSimg<double> fits(std::vector<size_t>{ncols, image.size() / ncols});
            *(fits.GetData<double>()) = image;
            fits.Write(file.string());
            std::cout << "[ IPC      ] saved " << file.string() << std::endl;
        };
        save(input,  "_input.fits");
        save(output, "_output.fits");
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (a) Kernel read from the configuration
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCAIpcConfig, OffByDefault)
{
    nisp::SCA sca(8, 6);
    sca.initialize(PixelConfig{{"gain", 1.0}});
    expectNoKernel(sca);

    // The convolution then does nothing.
    std::valarray<double> image(0.0, sca.size());
    image[3 * 8 + 4] = 1.0;
    const std::valarray<double> before = image;
    sca.applyIPC(image, 0, 0, 8, 6);
    expectSameImage(image, before, 8, 0.0, "IPC off");
}

TEST(SCAIpcConfig, SquareKernelFromAnArray)
{
    nisp::SCA sca(8, 6);
    sca.initialize(squareCfg(kSimple));
    expectStaticKernel(sca, kSimple, 3, 3);

    // Read from JSON text, as in a configuration file.
    nisp::SCA fromText(8, 6);
    fromText.initialize(PixelConfig(nlohmann::json::parse(R"({"ipc": [0.025, 0.1, 0.025, 0.1, 0.5, 0.1, 0.025, 0.1, 0.025]})")));
    expectStaticKernel(fromText, kSimple, 3, 3);
}

TEST(SCAIpcConfig, AsymmetricWeightsKeepTheirOrder)
{
    const std::vector<double> weights = {0.01, 0.02, 0.03,
                                         0.04, 0.73, 0.05,
                                         0.06, 0.02, 0.04};
    nisp::SCA sca(8, 6);
    sca.initialize(squareCfg(weights));
    expectStaticKernel(sca, weights, 3, 3);
}

TEST(SCAIpcConfig, RectangularKernelFromAnObject)
{
    // 5 columns × 3 rows, then 3 columns × 5 rows: same weights, sum = 1.
    const std::vector<double> weights = {0.00, 0.01, 0.05, 0.02, 0.00,
                                         0.01, 0.06, 0.70, 0.07, 0.02,
                                         0.00, 0.01, 0.04, 0.01, 0.00};
    nisp::SCA sca(8, 6);
    sca.initialize(rectangularCfg(weights, 5, 3));
    expectStaticKernel(sca, weights, 5, 3);

    sca.initialize(rectangularCfg(weights, 3, 5));
    expectStaticKernel(sca, weights, 3, 5);

    // The object form of a square kernel.
    sca.initialize(rectangularCfg(kSimple, 3, 3));
    expectStaticKernel(sca, kSimple, 3, 3);
}

TEST(SCAIpcConfig, StaticKernelIsNormalized)
{
    nisp::SCA sca(8, 6);
    for (double scale : {0.5, 1.0, 10.0})   // sum of the weights
    {
        std::vector<double> weights = kSimple;
        for (double& w : weights)
            w *= scale;
        sca.initialize(squareCfg(weights));
        expectStaticKernel(sca, kSimple, 3, 3);
    }
}

TEST(SCAIpcConfig, EmptyArrayTurnsTheIpcOffAndAnAbsentKeyChangesNothing)
{
    nisp::SCA sca(8, 6);
    sca.initialize(squareCfg(kSimple));

    sca.initialize(PixelConfig{{"gain", 2.0}});   // no "ipc"
    expectStaticKernel(sca, kSimple, 3, 3);

    PixelConfig off;
    off.set("ipc", nlohmann::json::array());
    sca.initialize(off);
    expectNoKernel(sca);
}

// Case (a) of a kernel from a file: no pixel above 1.
TEST(SCAIpcConfig, CubeWithSumsUpToOne)
{
    // Sums of 1 and of 0.9, in float and in double: normalized to 1, without warning.
    constexpr size_t W = 8, H = 6;
    auto weight = [](size_t x, size_t y, size_t k) { return pixelWeight(x, y, k) * ((x + y) % 2 ? 0.9 : 1.0); };

    const FileRemover f{writeCube<float> ("ipc_cube_f.fits", W, H, 9, weight)};
    const FileRemover d{writeCube<double>("ipc_cube_d.fits", W, H, 9, weight)};

    for (const std::string& path : {f.path, d.path})
    {
        nisp::SCA sca(W, H);
        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(sca.initialize(fileCfg(path))) << path;
        const std::string messages = testing::internal::GetCapturedStderr();

        EXPECT_EQ(messages.find("sum to more than 1"), std::string::npos) << path << ": " << messages;
        expectDynamicKernel(sca, pixelWeight, path);
    }
}

// Case (b) of a kernel from a file: some pixels above 1.
TEST(SCAIpcConfig, CubeWithSumsAboveOneIsNormalizedWithAWarning)
{
    // Sums of 1.05 (one pixel out of 3), 0.95 and 1: all normalized to 1.
    constexpr size_t W = 8, H = 6;
    auto scale  = [](size_t x, size_t y) { return (x + y) % 3 == 0 ? 1.05 : ((x + y) % 3 == 1 ? 0.95 : 1.0); };
    auto weight = [scale](size_t x, size_t y, size_t k) { return pixelWeight(x, y, k) * scale(x, y); };

    size_t nAbove = 0;
    for (size_t y = 0; y < H; ++y)
        for (size_t x = 0; x < W; ++x)
            if (scale(x, y) > 1.0)
                ++nAbove;

    const FileRemover f{writeCube<float>("ipc_cube_above.fits", W, H, 9, weight)};
    nisp::SCA sca(W, H);

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(sca.initialize(fileCfg(f.path)));
    const std::string messages = testing::internal::GetCapturedStderr();

    EXPECT_NE(messages.find("sum to more than 1 in " + std::to_string(nAbove) + " of " + std::to_string(W * H) + " pixels"),
              std::string::npos) << "warning: " << messages;
    EXPECT_NE(messages.find("largest sum: 1.05"), std::string::npos) << "warning: " << messages;

    expectDynamicKernel(sca, pixelWeight);
}

TEST(SCAIpcConfig, EveryWayToNameTheFile)
{
    // Primary HDU: a 2-D image; HDU 2, EXTNAME "IPC": the cube.
    constexpr size_t W = 8, H = 6;
    const FileRemover f{testFile("ipc_two_hdu.fits")};
    {
        DSL::FITSimg<double> primary(std::vector<size_t>{W, H});
        primary.Write(f.path);

        auto cube = makeCube<float>(W, H, 9, pixelWeight);
        cube->SetName("IPC");
        DSL::FITSmanager file(f.path, false);
        file.AppendImage(*cube);
        file.Close();
    }

    using nlohmann::json;
    const std::vector<PixelConfig> forms{
        PixelConfig(json{{"ipc", {{"file", f.path}, {"hdu", 2}}}}),
        PixelConfig(json{{"ipc", {{"file", f.path}, {"hdu", "IPC"}}}}),
        PixelConfig(json{{"ipc", f.path}, {"ipc.hdu", 2}}),
        PixelConfig(json{{"ipc.file", f.path}, {"ipc.hdu", "IPC"}}),
    };
    for (size_t i = 0; i < forms.size(); ++i)
    {
        nisp::SCA sca(W, H);
        ASSERT_NO_THROW(sca.initialize(forms[i])) << "form " << i;
        expectDynamicKernel(sca, pixelWeight, "form " + std::to_string(i));
    }

    // The primary HDU is not a cube.
    nisp::SCA sca(W, H);
    EXPECT_THROW(sca.initialize(fileCfg(f.path)), std::invalid_argument);
    expectNoKernel(sca);
}

TEST(SCAIpcConfig, DetectorSizedCubeIsCropped)
{
    // 16 × 12 detector, 8 × 6 matrix at (3, 2).
    constexpr size_t W = 8, H = 6, x0 = 3, y0 = 2;
    const FileRemover f{writeCube<float>("ipc_cube_detector.fits", 16, 12, 9, pixelWeight)};

    PixelConfig cfg = fileCfg(f.path);
    cfg.set("detectorSize",   nlohmann::json::array({16, 12}));
    cfg.set("detectorOrigin", nlohmann::json::array({3, 2}));
    cfg.set("refBorder",      0);

    nisp::SCA sca(W, H);
    sca.initialize(cfg);
    expectDynamicKernel(sca, [](size_t x, size_t y, size_t k) { return pixelWeight(x0 + x, y0 + y, k); });
}

TEST(SCAIpcConfig, StaticAndDynamicKernelsReplaceEachOther)
{
    constexpr size_t W = 8, H = 6;
    const FileRemover f{writeCube<float>("ipc_cube_replace.fits", W, H, 9, pixelWeight)};

    nisp::SCA sca(W, H);
    sca.initialize(squareCfg(kSimple));
    expectStaticKernel(sca, kSimple, 3, 3);

    sca.initialize(fileCfg(f.path));
    expectDynamicKernel(sca, pixelWeight);

    sca.initialize(squareCfg(kSimple));
    expectStaticKernel(sca, kSimple, 3, 3);
}

TEST(SCAIpcConfig, InvalidKernelThrowsAndKeepsThePreviousOne)
{
    constexpr size_t W = 8, H = 6;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    nisp::SCA sca(W, H);
    sca.initialize(squareCfg(kSimple));

    auto expectThrow = [&sca](const PixelConfig& cfg, const std::string& what)
    {
        SCOPED_TRACE(what);
        EXPECT_THROW(sca.initialize(cfg), std::invalid_argument);
        expectStaticKernel(sca, kSimple, 3, 3);
    };

    // ── Not a kernel ──
    expectThrow(PixelConfig{{"ipc", 0.1}},  "a number");
    expectThrow(PixelConfig{{"ipc", true}}, "a boolean");

    // ── Square kernel ──
    expectThrow(squareCfg({0.1, 0.8, 0.1}),                                "3 weights");
    expectThrow(squareCfg({0.1, 0.4, 0.4, 0.1}),                           "2 × 2");
    expectThrow(squareCfg({0.0, -0.1, 0.0, 0.1, 0.8, 0.1, 0.0, 0.1, 0.0}), "negative weight");
    expectThrow(squareCfg({0.0, nan, 0.0, 0.1, 0.6, 0.1, 0.0, 0.1, 0.0}),  "NaN weight");
    expectThrow(squareCfg(std::vector<double>(9, 0.0)),                    "null kernel");
    {
        PixelConfig cfg;
        cfg.set("ipc", nlohmann::json::parse(R"([0, "a", 0, 0.1, 0.6, 0.1, 0, 0.1, 0])"));
        expectThrow(cfg, "a weight that is not a number");
    }

    // ── Rectangular kernel ──
    expectThrow(rectangularCfg(kSimple, 5, 3),                       "9 weights for 5 × 3");
    expectThrow(rectangularCfg(std::vector<double>(8, 0.125), 2, 4), "even sizes");
    expectThrow(rectangularCfg(kSimple, -3, -3),                     "negative sizes");
    for (const char* text : {R"({"kernel": [0.1, 0.8, 0.1], "width": 3})",
                             R"({"kernel": [0.1, 0.8, 0.1], "width": 3, "height": "1"})",
                             R"({"kernel": [0.1, 0.8, 0.1], "width": 3, "height": 1, "dynamic": true})",
                             R"({"kernel": 1.0, "width": 1, "height": 1})",
                             R"({"width": 3, "height": 3})"})
    {
        PixelConfig cfg;
        cfg.set("ipc", nlohmann::json::parse(text));
        expectThrow(cfg, text);
    }

    // ── FITS cube ──
    const FileRemover badLayers{writeCube<float>("ipc_bad_layers.fits", W, H, 8, [](size_t, size_t, size_t) { return 0.125; })};
    const FileRemover badSize  {writeCube<float>("ipc_bad_size.fits", W + 1, H, 9, pixelWeight)};
    const FileRemover badNull  {writeCube<float>("ipc_bad_null.fits", W, H, 9, [](size_t x, size_t y, size_t k)
                                                 { return (x == 2 && y == 3) ? 0.0 : pixelWeight(x, y, k); })};
    const FileRemover badNaN   {writeCube<float>("ipc_bad_nan.fits", W, H, 9, [nan](size_t x, size_t y, size_t k)
                                                 { return (x == 2 && y == 3 && k == 5) ? nan : pixelWeight(x, y, k); })};
    expectThrow(fileCfg(badLayers.path), "8 layers");
    expectThrow(fileCfg(badSize.path),   "neither matrix- nor detector-sized");
    expectThrow(fileCfg(badNull.path),   "a pixel with a null kernel");
    expectThrow(fileCfg(badNaN.path),    "a pixel with a NaN weight");
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (b) Convolution by a static kernel
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCAIpcConvolution, StaticFullFrame)
{
    nisp::SCA sca(kN, kN);
    sca.initialize(squareCfg(kSimple));

    const std::valarray<double> sources = pointSources(kN, kN);
    std::valarray<double>       image   = sources;
    sca.applyIPC(image, 0, 0, kN, kN);
    saveImages(sources, image, kN, "frame");

    // The first source, pixel by pixel.
    const auto at = [&image](size_t x, size_t y) { return image[y * kN + x]; };
    EXPECT_NEAR(at(32, 32), 0.5, 1e-12);
    EXPECT_NEAR(at(31, 32), 0.1, 1e-12);
    EXPECT_NEAR(at(33, 32), 0.1, 1e-12);
    EXPECT_NEAR(at(32, 31), 0.1, 1e-12);
    EXPECT_NEAR(at(32, 33), 0.1, 1e-12);
    for (size_t y : {31, 33})
        for (size_t x : {31, 33})
            EXPECT_NEAR(at(x, y), 0.025, 1e-12) << "diagonal";
    EXPECT_NEAR(at(34, 32), 0.0, 1e-12);

    expectSameImage(image, expectedIpc(sources, kN, kN, kSimple, 3, 3), kN, 1e-12, "full frame");
    EXPECT_NEAR(image.sum(), sources.sum(), 1e-9) << "charge";   // 32 × 32 sources
}

TEST(SCAIpcConvolution, StaticWindow)
{
    nisp::SCA sca(kN, kN);
    sca.initialize(squareCfg(kSimple));

    const std::valarray<double> sources = pointSources(kN, kN);
    std::valarray<double>       frame   = sources;
    sca.applyIPC(frame, 0, 0, kN, kN);

    for (const Window& w : {kInside, kOnSource})
    {
        const std::valarray<double> in = crop(sources, kN, w);
        std::valarray<double>       image = in;
        sca.applyIPC(image, w.c0, w.r0, w.ncols, w.nrows);
        saveImages(in, image, w.ncols, fileLabel(w));
        expectSameImage(image, expectedIpc(in, w.ncols, w.nrows, kSimple, 3, 3), w.ncols, 1e-12, name(w));

        // No source near the border of this window: the same as in the full frame.
        if (w.c0 == kInside.c0)
            expectSameImage(image, crop(frame, kN, w), w.ncols, 1e-12, name(w) + " and full frame");
    }

    // A source in the corner of the window: its missing neighbours are ignored. Each pixel is
    // divided by the sum of the weights of the neighbours it has (the kernel is symmetric).
    const double atCorner = kSimple[4] + kSimple[5] + kSimple[7] + kSimple[8];   // 2 × 2 pixels
    const double onEdge   = 1.0 - kSimple[0] - kSimple[1] - kSimple[2];          // one row or column missing

    std::valarray<double> corner = crop(sources, kN, kOnSource);
    sca.applyIPC(corner, kOnSource.c0, kOnSource.r0, kOnSource.ncols, kOnSource.nrows);
    EXPECT_NEAR(corner[0],               kSimple[4] / atCorner, 1e-12);
    EXPECT_NEAR(corner[1],               kSimple[5] / onEdge,   1e-12);
    EXPECT_NEAR(corner[kOnSource.ncols], kSimple[7] / onEdge,   1e-12);
}

TEST(SCAIpcConvolution, StaticAsymmetricKernelDirection)
{
    // The kernel is the response to a charge in its central pixel: 0.2 goes to the right
    // neighbour and 0.1 to the upper one (next row).
    const std::vector<double> oriented = {0.0, 0.0, 0.0,
                                          0.0, 0.7, 0.2,
                                          0.0, 0.1, 0.0};
    nisp::SCA sca(kN, kN);
    sca.initialize(squareCfg(oriented));

    const std::valarray<double> sources = pointSources(kN, kN);
    std::valarray<double>       image   = sources;
    sca.applyIPC(image, 0, 0, kN, kN);
    saveImages(sources, image, kN, "frame_oriented");
    EXPECT_NEAR(image[32 * kN + 32], 0.7, 1e-12);
    EXPECT_NEAR(image[32 * kN + 33], 0.2, 1e-12) << "right";
    EXPECT_NEAR(image[33 * kN + 32], 0.1, 1e-12) << "up";
    EXPECT_NEAR(image[32 * kN + 31], 0.0, 1e-12) << "left";
    EXPECT_NEAR(image[31 * kN + 32], 0.0, 1e-12) << "down";
    expectSameImage(image, expectedIpc(sources, kN, kN, oriented, 3, 3), kN, 1e-12, "oriented kernel");

    // 5 columns × 3 rows.
    const std::vector<double> wide = {0.00, 0.01, 0.05, 0.02, 0.00,
                                      0.01, 0.06, 0.70, 0.07, 0.02,
                                      0.00, 0.01, 0.04, 0.01, 0.00};
    sca.initialize(rectangularCfg(wide, 5, 3));
    image = sources;
    sca.applyIPC(image, 0, 0, kN, kN);
    saveImages(sources, image, kN, "frame_5x3");
    expectSameImage(image, expectedIpc(sources, kN, kN, wide, 5, 3), kN, 1e-12, "5 × 3");

    for (const Window& w : {kInside, kOnSource})
    {
        const std::valarray<double> in = crop(sources, kN, w);
        std::valarray<double>       window = in;
        sca.applyIPC(window, w.c0, w.r0, w.ncols, w.nrows);
        saveImages(in, window, w.ncols, "5x3_" + fileLabel(w));
        expectSameImage(window, expectedIpc(in, w.ncols, w.nrows, wide, 5, 3), w.ncols, 1e-12, "5 × 3, " + name(w));
    }
}

TEST(SCAIpcConvolution, FlatImageStaysFlat)
{
    // Same kernel everywhere, static or from a file: a flat image (the baseline) is unchanged,
    // border included.
    constexpr size_t W = 96, H = 64;
    const FileRemover f{writeCube<float>("ipc_cube_flat.fits", W, H, 9, [](size_t, size_t, size_t k) { return kSimple[k]; })};

    for (const PixelConfig& cfg : {squareCfg(kSimple), fileCfg(f.path)})
    {
        nisp::SCA sca(W, H);
        sca.initialize(cfg);
        const std::valarray<double> flat(1000.0, W * H);
        std::valarray<double>       image = flat;
        sca.applyIPC(image, 0, 0, W, H);
        saveImages(flat, image, W, sca.isIPCDynamic() ? "dynamic" : "static");
        expectSameImage(image, flat, W, 1e-3, sca.isIPCDynamic() ? "dynamic" : "static");
    }
}

TEST(SCAIpcConvolution, InvalidWindowThrows)
{
    constexpr size_t W = 8, H = 6;
    const FileRemover f{writeCube<float>("ipc_cube_window.fits", W, H, 9, pixelWeight)};

    for (const PixelConfig& cfg : {squareCfg(kSimple), fileCfg(f.path)})
    {
        nisp::SCA sca(W, H);
        sca.initialize(cfg);

        std::valarray<double> tooSmall(0.0, W * H - 1);
        EXPECT_THROW(sca.applyIPC(tooSmall, 0, 0, W, H), std::invalid_argument);

        std::valarray<double> image(0.0, W * H);
        EXPECT_THROW(sca.applyIPC(image, 1, 0, W, H), std::out_of_range);
        EXPECT_THROW(sca.applyIPC(image, 0, 1, W, H), std::out_of_range);
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (c) Convolution by the kernels of a FITS cube
// ═══════════════════════════════════════════════════════════════════════════════

TEST(SCAIpcConvolution, CubeOfTheSameKernelIsTheStaticKernel)
{
    // Every pixel of the cube holds the same asymmetric kernel: same image as with the static
    // one, to float precision, border included (the image is not null there).
    const std::vector<double> kernel = {0.01, 0.02, 0.03,
                                        0.04, 0.73, 0.05,
                                        0.06, 0.02, 0.04};
    constexpr size_t W = 256, H = 192;
    const FileRemover f{writeCube<float>("ipc_cube_static.fits", W, H, 9, [&kernel](size_t, size_t, size_t k) { return kernel[k]; })};

    nisp::SCA fixed(W, H), dynamic(W, H);
    fixed  .initialize(squareCfg(kernel));
    dynamic.initialize(fileCfg(f.path));
    ASSERT_FALSE(fixed  .isIPCDynamic());
    ASSERT_TRUE (dynamic.isIPCDynamic());

    std::valarray<double> scene(0.0, W * H);
    for (size_t i = 0; i < scene.size(); ++i)
        scene[i] = static_cast<double>((7919 * i) % 1000) / 1000.0;   // [0, 1), no structure

    for (const Window& w : {Window{0, 0, W, H}, Window{17, 33, 120, 90}})
    {
        const std::valarray<double> in = crop(scene, W, w);
        std::valarray<double>       a = in, b = in;
        fixed  .applyIPC(a, w.c0, w.r0, w.ncols, w.nrows);
        dynamic.applyIPC(b, w.c0, w.r0, w.ncols, w.nrows);
        saveImages(in, a, w.ncols, "static_"  + fileLabel(w));
        saveImages(in, b, w.ncols, "dynamic_" + fileLabel(w));
        expectSameImage(b, a, w.ncols, 1e-6, name(w));
    }
}

TEST(SCAIpcConvolution, CubeFullFrameAndWindow)
{
    // Kernels that differ from pixel to pixel: each source gives an image of its own kernel.
    // Full-detector cube: 151 MB, removed at the end of the test.
    const FileRemover f{writeCube<float>("ipc_cube_2048.fits", kN, kN, 9, pixelWeight)};
    nisp::SCA sca(kN, kN);
    sca.initialize(fileCfg(f.path));
    ASSERT_TRUE(sca.isIPCDynamic());

    // ── Full frame ──
    const std::valarray<double> sources = pointSources(kN, kN);
    std::valarray<double>       frame   = sources;
    sca.applyIPC(frame, 0, 0, kN, kN);
    saveImages(sources, frame, kN, "frame");

    // Two sources, with the weights of their own kernels: (32, 32) and (96, 32).
    for (size_t x : {32, 96})
    {
        const double a = 0.01 * static_cast<double>(1 + (x + 3 * 32) % 7);
        EXPECT_NEAR(frame[32 * kN + x],     1.0 - 3.75 * a, 1e-6) << "x = " << x;
        EXPECT_NEAR(frame[32 * kN + x + 1], a,              1e-6) << "right, x = " << x;
        EXPECT_NEAR(frame[32 * kN + x - 1], 0.50 * a,       1e-6) << "left, x = "  << x;
        EXPECT_NEAR(frame[33 * kN + x],     2.00 * a,       1e-6) << "up, x = "    << x;
        EXPECT_NEAR(frame[31 * kN + x],     0.25 * a,       1e-6) << "down, x = "  << x;
    }

    expectSameImage(frame, expectedIpc(sources, kN, kN, 3, 3, pixelWeight), kN, 1e-6, "full frame");
    EXPECT_NEAR(frame.sum(), sources.sum(), 1e-3) << "charge";   // 32 × 32 sources

    // ── Windows: the kernels are those of the window pixels on the matrix ──
    for (const Window& w : {kInside, kOnSource})
    {
        const std::valarray<double> in = crop(sources, kN, w);
        std::valarray<double>       image = in;
        sca.applyIPC(image, w.c0, w.r0, w.ncols, w.nrows);
        saveImages(in, image, w.ncols, fileLabel(w));

        auto weight = [&w](size_t x, size_t y, size_t k) { return pixelWeight(w.c0 + x, w.r0 + y, k); };
        expectSameImage(image, expectedIpc(in, w.ncols, w.nrows, 3, 3, weight), w.ncols, 1e-6, name(w));

        // No source near the border of this window: the same as in the full frame.
        if (w.c0 == kInside.c0)
            expectSameImage(image, crop(frame, kN, w), w.ncols, 1e-12, name(w) + " and full frame");
    }
}
