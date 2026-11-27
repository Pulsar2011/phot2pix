//
// Created by gillard on 2026-11-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_pixelConvolution.cxx
//
// Unit tests of phot2pix::PixelConvolution, without detector nor file: the FITS cubes are made
// in memory and given by the loader.
// (a) initialize(): static and dynamic kernels, loader, crop, normalization, invalid kernels;
// (b) apply() with a static kernel;
// (c) apply() with a dynamic kernel;
// (d) setKernel(): kernels that are not normalized (channel crosstalk).
//
// gTest_nispIPC.cxx tests it within nisp::SCA, with FITS files and full-detector images.
//

#include <gtest/gtest.h>
#include <phot2pix/DetectorEffect.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <valarray>
#include <vector>

namespace
{
    using phot2pix::PixelConfig;
    using phot2pix::PixelConvolution;
    using CubePtr = std::shared_ptr<DSL::FITScube>;

    // ── Kernels ──────────────────────────────────────────────────────────────────

    // 9 different weights: a flip, a rotation or a transposition of the kernel is seen. Sum = 45.
    const std::vector<double> kDistinct = {1, 2, 3,
                                           4, 5, 6,
                                           7, 8, 9};

    std::vector<double> normalized(std::vector<double> weights)
    {
        double sum = 0.0;
        for (double w : weights) sum += w;
        for (double& w : weights) w /= sum;
        return weights;
    }

    // Weight k of the kernel of pixel (x, y), of nK weights: different from pixel to pixel and
    // from neighbour to neighbour. Sum = 1.
    double pixelKernel(size_t x, size_t y, size_t k, size_t nK)
    {
        const auto raw = [x, y, nK](size_t j) { return (j == nK / 2) ? 100.0 : 1.0 + static_cast<double>((3 * x + 5 * y + 7 * j) % 11); };
        double sum = 0.0;
        for (size_t j = 0; j < nK; ++j)
            sum += raw(j);
        return raw(k) / sum;
    }

    double pixelKernel9(size_t x, size_t y, size_t k) { return pixelKernel(x, y, k, 9); }

    // Cube of W × H × nK values of type T, in memory: value(x, y, k).
    template<typename T, typename F>
    CubePtr makeCube(size_t W, size_t H, size_t nK, F value)
    {
        auto cube = std::make_shared<DSL::FITSimg<T>>(std::vector<size_t>{W, H, nK});
        std::valarray<T>& data = *(cube->template GetData<T>());
        for (size_t k = 0; k < nK; ++k)
            for (size_t y = 0; y < H; ++y)
                for (size_t x = 0; x < W; ++x)
                    data[(k * H + y) * W + x] = static_cast<T>(value(x, y, k));
        return cube;
    }

    // ── Configurations ───────────────────────────────────────────────────────────

    PixelConfig jsonCfg(const nlohmann::json& ipc)
    {
        PixelConfig cfg;
        cfg.set("ipc", ipc);
        return cfg;
    }

    PixelConfig arrayCfg(const std::vector<double>& weights) { return jsonCfg(nlohmann::json(weights)); }

    PixelConfig objectCfg(const std::vector<double>& weights, int w, int h)
    {
        nlohmann::json ipc;
        ipc["kernel"] = weights;
        ipc["width"]  = w;
        ipc["height"] = h;
        return jsonCfg(ipc);
    }

    // The path is never opened: the loader gives the cube.
    PixelConfig fileCfg()
    {
        PixelConfig cfg;
        cfg.set("ipc", std::string("cube.fits"));
        return cfg;
    }

    // Loader of the kernels that are not in a file.
    CubePtr notCalled()
    {
        ADD_FAILURE() << "the FITS loader must not be called";
        return nullptr;
    }

    // ── Convolutions of a W × H matrix, which is the whole detector ──────────────

    void initialize(PixelConvolution& c, const PixelConfig& cfg, size_t W, size_t H, const PixelConvolution::CubeLoader& loader = notCalled)
    {
        c.initialize(cfg, phot2pix::Geometry(W, H), W, H, loader, "ipc");
    }

    PixelConvolution staticConvolution(const PixelConfig& cfg, size_t W, size_t H)
    {
        PixelConvolution c;
        initialize(c, cfg, W, H);
        return c;
    }

    PixelConvolution dynamicConvolution(const CubePtr& cube, size_t W, size_t H)
    {
        PixelConvolution c;
        initialize(c, fileCfg(), W, H, [cube]() { return cube; });
        return c;
    }

    // ── Checks ───────────────────────────────────────────────────────────────────

    void expectOff(const PixelConvolution& c)
    {
        EXPECT_FALSE(c.isEnabled());
        EXPECT_FALSE(c.isDynamic());
        EXPECT_FALSE(c.getKernelCube());
        ASSERT_EQ(c.getKernel().size(), 1u);
        EXPECT_DOUBLE_EQ(c.getKernel().data()[0], 1.0);
    }

    void expectStatic(const PixelConvolution& c, const std::vector<double>& weights, size_t w, size_t h)
    {
        EXPECT_TRUE (c.isEnabled());
        EXPECT_FALSE(c.isDynamic());
        EXPECT_FALSE(c.getKernelCube());
        ASSERT_EQ(c.getKernel().sizes().size(), 2u);
        EXPECT_EQ(c.getKernel().size(0), w) << "columns";
        EXPECT_EQ(c.getKernel().size(1), h) << "rows";
        ASSERT_EQ(c.getKernel().size(), weights.size());
        for (size_t k = 0; k < weights.size(); ++k)
            EXPECT_NEAR(c.getKernel().data()[k], weights[k], 1e-14) << "weight " << k;
    }

    template<typename F>
    void expectDynamic(const PixelConvolution& c, size_t W, size_t H, size_t nK, F weight)
    {
        EXPECT_TRUE(c.isEnabled());
        EXPECT_TRUE(c.isDynamic());

        const auto cube = c.getKernelCube();
        ASSERT_TRUE(cube);
        const std::valarray<float>* data = cube->GetData<float>();
        ASSERT_NE(data, nullptr) << "the cube is stored in float";
        ASSERT_EQ(data->size(), W * H * nK);
        for (size_t k = 0; k < nK; ++k)
            for (size_t y = 0; y < H; ++y)
                for (size_t x = 0; x < W; ++x)
                    ASSERT_NEAR((*data)[(k * H + y) * W + x], weight(x, y, k), 1e-6) << "x = " << x << ", y = " << y << ", layer " << k;
    }

    void expectSame(const std::valarray<double>& image, const std::valarray<double>& expected, size_t ncols, double tolerance)
    {
        ASSERT_EQ(image.size(), expected.size());
        for (size_t i = 0; i < image.size(); ++i)
            ASSERT_NEAR(image[i], expected[i], tolerance) << "x = " << i % ncols << ", y = " << i / ncols;
    }

    // ── Images ───────────────────────────────────────────────────────────────────

    // 1 in one pixel, 0 elsewhere.
    std::valarray<double> point(size_t W, size_t H, size_t x, size_t y)
    {
        std::valarray<double> image(0.0, W * H);
        image[y * W + x] = 1.0;
        return image;
    }

    // Every pixel between 0 and 1, without structure.
    std::valarray<double> scene(size_t W, size_t H)
    {
        std::valarray<double> image(0.0, W * H);
        for (size_t i = 0; i < image.size(); ++i)
            image[i] = static_cast<double>((7919 * i + 13) % 1000) / 1000.0;
        return image;
    }

    std::valarray<double> section(const std::valarray<double>& image, size_t W, size_t c0, size_t r0, size_t ncols, size_t nrows)
    {
        return image[std::gslice(r0 * W + c0, std::valarray<size_t>{nrows, ncols}, std::valarray<size_t>{W, 1})];
    }

    // Reference: each pixel distributes its charge with its own kernel of kw × kh weights,
    // weight(x, y, k). A pixel that misses neighbours is divided by the sum of its weights for the
    // neighbours it has.
    template<typename Weight>
    std::valarray<double> reference(const std::valarray<double>& image, size_t ncols, size_t nrows, size_t kw, size_t kh, Weight weight)
    {
        using idx = std::ptrdiff_t;
        const idx W = static_cast<idx>(ncols), H = static_cast<idx>(nrows), Kw = static_cast<idx>(kw), Kh = static_cast<idx>(kh);
        const auto inside = [W, H](idx x, idx y) { return x >= 0 && x < W && y >= 0 && y < H; };

        std::valarray<double> out(0.0, ncols * nrows), norm(0.0, ncols * nrows);
        for (idx y = 0; y < H; ++y)
            for (idx x = 0; x < W; ++x)
                for (idx k = 0; k < Kw * Kh; ++k)
                {
                    const idx    dx = k % Kw - Kw / 2, dy = k / Kw - Kh / 2;
                    const double w  = weight(static_cast<size_t>(x), static_cast<size_t>(y), static_cast<size_t>(k));
                    if (inside(x + dx, y + dy))
                        out[static_cast<size_t>((y + dy) * W + x + dx)] += image[static_cast<size_t>(y * W + x)] * w;
                    if (inside(x - dx, y - dy))
                        norm[static_cast<size_t>(y * W + x)] += w;
                }
        return out / norm;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (a) initialize()
// ═══════════════════════════════════════════════════════════════════════════════

TEST(PixelConvolution, OffByDefault)
{
    PixelConvolution c;
    expectOff(c);

    // No "ipc" key.
    initialize(c, PixelConfig{{"gain", 1.0}}, 5, 4);
    expectOff(c);

    // apply() does nothing, whatever the window and the size of the image.
    std::valarray<double> image = scene(5, 4);
    const std::valarray<double> before = image;
    EXPECT_NO_THROW(c.apply(image, 0, 0, 5, 4, 1));
    EXPECT_NO_THROW(c.apply(image, 100, 100, 50, 50, 1));
    expectSame(image, before, 5, 0.0);
}

TEST(PixelConvolution, SquareKernelFromAnArray)
{
    const PixelConvolution c = staticConvolution(arrayCfg(kDistinct), 8, 6);
    expectStatic(c, normalized(kDistinct), 3, 3);
}

TEST(PixelConvolution, SquareKernelOf5x5)
{
    std::vector<double> weights(25);
    for (size_t k = 0; k < weights.size(); ++k)
        weights[k] = static_cast<double>(k);   // a null weight is valid
    expectStatic(staticConvolution(arrayCfg(weights), 8, 6), normalized(weights), 5, 5);
}

TEST(PixelConvolution, RectangularKernelFromAnObject)
{
    const std::vector<double> weights = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    expectStatic(staticConvolution(objectCfg(weights, 5, 3), 8, 6), normalized(weights), 5, 3);
    expectStatic(staticConvolution(objectCfg(weights, 3, 5), 8, 6), normalized(weights), 3, 5);
    expectStatic(staticConvolution(objectCfg(weights, 15, 1), 8, 6), normalized(weights), 15, 1);
}

TEST(PixelConvolution, LoaderIsOnlyCalledForAFile)
{
    constexpr size_t W = 5, H = 4;
    size_t calls = 0;
    const PixelConvolution::CubeLoader loader = [&calls]() { ++calls; return makeCube<float>(W, H, 9, pixelKernel9); };

    PixelConvolution c;
    for (const PixelConfig& cfg : {PixelConfig{}, arrayCfg(kDistinct), objectCfg(kDistinct, 3, 3), jsonCfg(nlohmann::json::array())})
        initialize(c, cfg, W, H, loader);
    EXPECT_EQ(calls, 0u);

    // "ipc": "<path>", {"file", "hdu"} and "ipc.file".
    nlohmann::json file;
    file["file"] = "cube.fits";
    file["hdu"]  = 2;
    PixelConfig flat;
    flat.set("ipc.file", std::string("cube.fits"));

    size_t expected = 0;
    for (const PixelConfig& cfg : {fileCfg(), jsonCfg(file), flat})
    {
        PixelConvolution d;
        initialize(d, cfg, W, H, loader);
        EXPECT_EQ(calls, ++expected);
        expectDynamic(d, W, H, 9, pixelKernel9);
    }
}

TEST(PixelConvolution, CubeOfAnyTypeIsCopiedInFloat)
{
    constexpr size_t W = 5, H = 4;
    expectDynamic(dynamicConvolution(makeCube<float> (W, H, 9, pixelKernel9), W, H), W, H, 9, pixelKernel9);
    expectDynamic(dynamicConvolution(makeCube<double>(W, H, 9, pixelKernel9), W, H), W, H, 9, pixelKernel9);

    // Integers: the kernel that keeps everything in the central pixel.
    auto identity = [](size_t, size_t, size_t k) { return k == 4 ? 1.0 : 0.0; };
    expectDynamic(dynamicConvolution(makeCube<uint8_t>(W, H, 9, identity), W, H), W, H, 9, identity);

    // 5 × 5 kernel.
    auto kernel25 = [](size_t x, size_t y, size_t k) { return pixelKernel(x, y, k, 25); };
    expectDynamic(dynamicConvolution(makeCube<float>(W, H, 25, kernel25), W, H), W, H, 25, kernel25);
}

TEST(PixelConvolution, CubeOfTheLoaderIsNotModified)
{
    // Sums of 0.5: the stored cube is normalized, not the one of the loader.
    constexpr size_t W = 5, H = 4;
    const CubePtr source = makeCube<float>(W, H, 9, [](size_t x, size_t y, size_t k) { return 0.5 * pixelKernel9(x, y, k); });
    const std::valarray<float> before = *(source->GetData<float>());

    const PixelConvolution c = dynamicConvolution(source, W, H);
    expectDynamic(c, W, H, 9, pixelKernel9);

    EXPECT_NE(c.getKernelCube().get(), source.get());
    const std::valarray<float>& after = *(source->GetData<float>());
    ASSERT_EQ(after.size(), before.size());
    for (size_t i = 0; i < before.size(); ++i)
        ASSERT_EQ(after[i], before[i]);
}

TEST(PixelConvolution, DetectorSizedCubeIsCropped)
{
    // 16 × 12 detector, 8 × 6 matrix at (3, 2).
    phot2pix::Geometry geometry(16, 12);
    geometry.origin = std::pair<size_t, size_t>{3, 2};

    const CubePtr cube = makeCube<float>(16, 12, 9, pixelKernel9);
    PixelConvolution c;
    c.initialize(fileCfg(), geometry, 8, 6, [cube]() { return cube; }, "ipc");
    expectDynamic(c, 8, 6, 9, [](size_t x, size_t y, size_t k) { return pixelKernel9(3 + x, 2 + y, k); });

    // A matrix-sized cube is used as is.
    const CubePtr matrix = makeCube<float>(8, 6, 9, pixelKernel9);
    c.initialize(fileCfg(), geometry, 8, 6, [matrix]() { return matrix; } ,"ipc");
    expectDynamic(c, 8, 6, 9, pixelKernel9);
}

TEST(PixelConvolution, SumsUpToOneAreNormalizedWithoutWarning)
{
    constexpr size_t W = 5, H = 4;
    auto weight = [](size_t x, size_t y, size_t k) { return pixelKernel9(x, y, k) * ((x + y) % 2 ? 0.9 : 1.0); };

    testing::internal::CaptureStderr();
    const PixelConvolution c = dynamicConvolution(makeCube<float>(W, H, 9, weight), W, H);
    const std::string messages = testing::internal::GetCapturedStderr();

    EXPECT_EQ(messages.find("warning"), std::string::npos) << messages;
    expectDynamic(c, W, H, 9, pixelKernel9);
}

TEST(PixelConvolution, SumsAboveOneAreNormalizedWithAWarning)
{
    // 3 pixels of 20 above 1, the largest at 1.25.
    constexpr size_t W = 5, H = 4;
    auto scale = [](size_t x, size_t y)
    {
        if (x == 0 && y == 0) return 1.25;
        if (x == 4 && y == 3) return 1.10;
        if (x == 2 && y == 1) return 1.01;
        return 1.0;
    };
    auto weight = [scale](size_t x, size_t y, size_t k) { return pixelKernel9(x, y, k) * scale(x, y); };

    testing::internal::CaptureStderr();
    const PixelConvolution c = dynamicConvolution(makeCube<double>(W, H, 9, weight), W, H);
    const std::string messages = testing::internal::GetCapturedStderr();

    EXPECT_NE(messages.find("sum to more than 1 in 3 of 20 pixels (largest sum: 1.25)"), std::string::npos) << messages;
    expectDynamic(c, W, H, 9, pixelKernel9);
}

TEST(PixelConvolution, KernelsReplaceEachOther)
{
    constexpr size_t W = 5, H = 4;
    const CubePtr cube = makeCube<float>(W, H, 9, pixelKernel9);
    const PixelConvolution::CubeLoader loader = [cube]() { return cube; };

    PixelConvolution c;
    initialize(c, arrayCfg(kDistinct), W, H);
    expectStatic(c, normalized(kDistinct), 3, 3);

    initialize(c, fileCfg(), W, H, loader);
    expectDynamic(c, W, H, 9, pixelKernel9);

    initialize(c, PixelConfig{}, W, H);   // no "ipc": unchanged
    expectDynamic(c, W, H, 9, pixelKernel9);

    initialize(c, arrayCfg(kDistinct), W, H);
    expectStatic(c, normalized(kDistinct), 3, 3);

    initialize(c, jsonCfg(nlohmann::json::array()), W, H);
    expectOff(c);
}

TEST(PixelConvolution, InvalidStaticKernelThrowsAndKeepsThePreviousOne)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    PixelConvolution c = staticConvolution(arrayCfg(kDistinct), 8, 6);
    auto expectThrow = [&c](const PixelConfig& cfg, const std::string& what)
    {
        SCOPED_TRACE(what);
        EXPECT_THROW(initialize(c, cfg, 8, 6), std::invalid_argument);
        expectStatic(c, normalized(kDistinct), 3, 3);
    };

    expectThrow(arrayCfg({1, 2}),                        "2 weights");
    expectThrow(arrayCfg({1, 2, 3}),                     "3 weights");
    expectThrow(arrayCfg(std::vector<double>(4, 1.0)),   "2 × 2");
    expectThrow(arrayCfg(std::vector<double>(16, 1.0)),  "4 × 4");
    expectThrow(arrayCfg(std::vector<double>(9, 0.0)),   "null sum");
    expectThrow(arrayCfg({0, -1, 0, 1, 6, 1, 0, 1, 0}),  "negative weight");
    expectThrow(arrayCfg({0, nan, 0, 1, 6, 1, 0, 1, 0}), "NaN weight");
    expectThrow(arrayCfg({0, inf, 0, 1, 6, 1, 0, 1, 0}), "infinite weight");
    expectThrow(objectCfg(kDistinct, 5, 3),              "9 weights for 5 × 3");
    expectThrow(objectCfg(kDistinct, 9, 0),              "no row");
    expectThrow(objectCfg(std::vector<double>(6, 1.0), 2, 3), "even number of columns");
    expectThrow(objectCfg(std::vector<double>(6, 1.0), 3, 2), "even number of rows");
    expectThrow(objectCfg(kDistinct, -3, -3),            "negative sizes");
    for (const char* text : {R"([0, "a", 0, 1, 6, 1, 0, 1, 0])",
                             R"([[0, 1, 0], [1, 6, 1], [0, 1, 0]])",
                             R"({"kernel": [1, 2, 1], "width": 3})",
                             R"({"kernel": [1, 2, 1], "height": 1})",
                             R"({"kernel": [1, 2, 1], "width": "3", "height": 1})",
                             R"({"kernel": [1, 2, 1], "width": 1.5, "height": 2})",
                             R"({"kernel": [1, 2, 1], "width": 3, "height": 1, "file": "cube.fits"})",
                             R"({"kernel": [], "width": 1, "height": 1})",
                             R"({"kernel": 1, "width": 1, "height": 1})"})
        expectThrow(jsonCfg(nlohmann::json::parse(text)), text);
}

TEST(PixelConvolution, InvalidCubeThrowsAndKeepsThePreviousOne)
{
    constexpr size_t W = 5, H = 4;
    PixelConvolution c = dynamicConvolution(makeCube<float>(W, H, 9, pixelKernel9), W, H);

    auto expectThrow = [&c](const CubePtr& cube, const std::string& what)
    {
        SCOPED_TRACE(what);
        EXPECT_THROW(initialize(c, fileCfg(), W, H, [cube]() { return cube; }), std::invalid_argument);
        expectDynamic(c, W, H, 9, pixelKernel9);
    };

    expectThrow(nullptr, "no image: a number or a boolean");
    expectThrow(CubePtr(std::make_shared<DSL::FITSimg<float>>(std::vector<size_t>{W, H})), "2-D image");
    expectThrow(makeCube<float>(W, H, 8,  pixelKernel9), "8 layers");
    expectThrow(makeCube<float>(W, H, 5,  pixelKernel9), "5 layers");
    expectThrow(makeCube<float>(W, H, 4,  pixelKernel9), "4 layers");
    expectThrow(makeCube<float>(H, W, 9,  pixelKernel9), "transposed");
    expectThrow(makeCube<float>(W + 1, H, 9, pixelKernel9), "one column too many");

    auto replaced = [](double value, size_t kBad)   // in pixel (2, 3): layer kBad, or all of them (9)
    {
        return [value, kBad](size_t x, size_t y, size_t k) { return (x == 2 && y == 3 && (k == kBad || kBad == 9)) ? value : pixelKernel9(x, y, k); };
    };
    expectThrow(makeCube<float>(W, H, 9, replaced(0.0, 9)),  "null kernel");
    expectThrow(makeCube<float>(W, H, 9, replaced(-0.1, 9)), "negative sum");
    expectThrow(makeCube<float>(W, H, 9, replaced(std::numeric_limits<double>::quiet_NaN(), 0)), "NaN weight");
    expectThrow(makeCube<float>(W, H, 9, replaced(std::numeric_limits<double>::infinity(),  0)), "infinite weight");

    // An error of the loader itself is not hidden.
    EXPECT_THROW(initialize(c, fileCfg(), W, H, []() -> CubePtr { throw std::runtime_error("no such file"); }), std::runtime_error);
    expectDynamic(c, W, H, 9, pixelKernel9);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (b) apply() with a static kernel
// ═══════════════════════════════════════════════════════════════════════════════

TEST(PixelConvolution, StaticPointGivesTheKernel)
{
    // kernel[0] goes to the neighbour (x − 1, y − 1), kernel[8] to (x + 1, y + 1).
    const PixelConvolution c = staticConvolution(arrayCfg(kDistinct), 7, 5);
    std::valarray<double>  image = point(7, 5, 3, 2);
    c.apply(image, 0, 0, 7, 5, 1);

    for (size_t k = 0; k < 9; ++k)
        EXPECT_NEAR(image[(1 + k / 3) * 7 + 2 + k % 3], kDistinct[k] / 45.0, 1e-14) << "weight " << k;
    EXPECT_NEAR(image.sum(), 1.0, 1e-14) << "charge";
}

TEST(PixelConvolution, StaticRectangularPointGivesTheKernel)
{
    const std::vector<double> weights = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};   // sum = 120
    const PixelConvolution c = staticConvolution(objectCfg(weights, 5, 3), 9, 5);
    std::valarray<double>  image = point(9, 5, 4, 2);
    c.apply(image, 0, 0, 9, 5, 1);

    for (size_t k = 0; k < 15; ++k)
        EXPECT_NEAR(image[(1 + k / 5) * 9 + 2 + k % 5], weights[k] / 120.0, 1e-14) << "weight " << k;
    EXPECT_NEAR(image.sum(), 1.0, 1e-14) << "charge";
}

TEST(PixelConvolution, StaticByHandOnOneRow)
{
    // A pixel gives 0.2 to its left, keeps 0.5 and gives 0.3 to its right.
    const PixelConvolution c = staticConvolution(objectCfg({0.2, 0.5, 0.3}, 3, 1), 5, 1);
    std::valarray<double>  image = {0.0, 1.0, 0.0, 2.0, 0.0};
    c.apply(image, 0, 0, 5, 1, 1);

    EXPECT_NEAR(image[0], 0.2 / 0.7,             1e-14) << "no left neighbour: 0.3 is missing";
    EXPECT_NEAR(image[1], 0.5,                   1e-14);
    EXPECT_NEAR(image[2], 0.3 * 1.0 + 0.2 * 2.0, 1e-14) << "from both sides";
    EXPECT_NEAR(image[3], 1.0,                   1e-14);
    EXPECT_NEAR(image[4], 0.6 / 0.8,             1e-14) << "no right neighbour: 0.2 is missing";
}

TEST(PixelConvolution, StaticMatchesTheReference)
{
    constexpr size_t W = 13, H = 9;
    const std::valarray<double> in = scene(W, H);

    const std::vector<double> square = normalized(kDistinct);
    std::valarray<double> image = in;
    staticConvolution(arrayCfg(kDistinct), W, H).apply(image, 0, 0, W, H, 1);
    expectSame(image, reference(in, W, H, 3, 3, [&square](size_t, size_t, size_t k) { return square[k]; }), W, 1e-13);

    const std::vector<double> raw = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, tall = normalized(raw);
    image = in;
    staticConvolution(objectCfg(raw, 3, 5), W, H).apply(image, 0, 0, W, H, 1);
    expectSame(image, reference(in, W, H, 3, 5, [&tall](size_t, size_t, size_t k) { return tall[k]; }), W, 1e-13);
}

TEST(PixelConvolution, KernelOfOneWeightChangesNothing)
{
    const std::valarray<double> in = scene(7, 5);
    std::valarray<double> image = in;
    staticConvolution(arrayCfg({3.0}), 7, 5).apply(image, 0, 0, 7, 5, 1);
    expectSame(image, in, 7, 1e-14);
}

TEST(PixelConvolution, StaticDoesNotDependOnThePositionOfTheWindow)
{
    const PixelConvolution c = staticConvolution(arrayCfg(kDistinct), 20, 15);
    std::valarray<double> a = scene(7, 5), b = a;
    c.apply(a, 0, 0, 7, 5, 1);
    c.apply(b, 13, 10, 7, 5, 1);   // up to the last pixel of the matrix
    expectSame(b, a, 7, 0.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (c) apply() with a dynamic kernel
// ═══════════════════════════════════════════════════════════════════════════════

TEST(PixelConvolution, DynamicPointGivesTheKernelOfItsPixel)
{
    constexpr size_t W = 9, H = 7;
    const PixelConvolution c = dynamicConvolution(makeCube<float>(W, H, 9, pixelKernel9), W, H);

    for (size_t y = 2; y + 2 < H; ++y)   // the image of the kernel does not reach the border
        for (size_t x = 2; x + 2 < W; ++x)
        {
            std::valarray<double> image = point(W, H, x, y);
            c.apply(image, 0, 0, W, H, 1);
            for (size_t k = 0; k < 9; ++k)
                ASSERT_NEAR(image[(y - 1 + k / 3) * W + x - 1 + k % 3], pixelKernel9(x, y, k), 1e-6)
                    << "x = " << x << ", y = " << y << ", weight " << k;
            ASSERT_NEAR(image.sum(), 1.0, 1e-6) << "charge";
        }
}

TEST(PixelConvolution, DynamicEachNeighbourUsesItsOwnKernel)
{
    // The pixel (3, 2), empty, between two pixels that hold 1 and 10: it receives the "right"
    // weight of its left neighbour and the "left" weight of its right neighbour.
    constexpr size_t W = 7, H = 5;
    const PixelConvolution c = dynamicConvolution(makeCube<float>(W, H, 9, pixelKernel9), W, H);

    std::valarray<double> image = point(W, H, 2, 2) + 10.0 * point(W, H, 4, 2);
    c.apply(image, 0, 0, W, H, 1);
    EXPECT_NEAR(image[2 * W + 3], 1.0 * pixelKernel9(2, 2, 5) + 10.0 * pixelKernel9(4, 2, 3), 1e-6);
}

TEST(PixelConvolution, DynamicMatchesTheReference)
{
    constexpr size_t W = 13, H = 9;
    const std::valarray<double> in = scene(W, H);

    std::valarray<double> image = in;
    dynamicConvolution(makeCube<float>(W, H, 9, pixelKernel9), W, H).apply(image, 0, 0, W, H, 1);
    expectSame(image, reference(in, W, H, 3, 3, pixelKernel9), W, 1e-6);

    auto kernel25 = [](size_t x, size_t y, size_t k) { return pixelKernel(x, y, k, 25); };
    image = in;
    dynamicConvolution(makeCube<float>(W, H, 25, kernel25), W, H).apply(image, 0, 0, W, H, 1);
    expectSame(image, reference(in, W, H, 5, 5, kernel25), W, 1e-6);
}

TEST(PixelConvolution, DynamicWindowUsesTheKernelsOfItsPixels)
{
    constexpr size_t W = 20, H = 15, c0 = 6, r0 = 4, w = 9, h = 7;
    const PixelConvolution      c     = dynamicConvolution(makeCube<float>(W, H, 9, pixelKernel9), W, H);
    const std::valarray<double> frame = scene(W, H);
    const std::valarray<double> in    = section(frame, W, c0, r0, w, h);

    std::valarray<double> window = in;
    c.apply(window, c0, r0, w, h, 1);
    expectSame(window, reference(in, w, h, 3, 3, [](size_t x, size_t y, size_t k) { return pixelKernel9(c0 + x, r0 + y, k); }), w, 1e-6);

    // Away from the border of the window: the same as in the full frame.
    std::valarray<double> full = frame;
    c.apply(full, 0, 0, W, H, 1);
    expectSame(section(window, w, 1, 1, w - 2, h - 2), section(full, W, c0 + 1, r0 + 1, w - 2, h - 2), w - 2, 1e-12);
}

TEST(PixelConvolution, CubeOfTheSameKernelIsTheStaticKernel)
{
    constexpr size_t W = 13, H = 9;
    const std::vector<double> kernel = normalized(kDistinct);
    const CubePtr cube = makeCube<float>(W, H, 9, [&kernel](size_t, size_t, size_t k) { return kernel[k]; });

    std::valarray<double> a = scene(W, H), b = a;   // border included
    staticConvolution(arrayCfg(kDistinct), W, H).apply(a, 0, 0, W, H, 1);
    dynamicConvolution(cube, W, H)              .apply(b, 0, 0, W, H, 1);
    expectSame(b, a, W, 1e-6);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (b) and (c): properties of both kernels
// ═══════════════════════════════════════════════════════════════════════════════

namespace
{
    // The two convolutions of a W × H matrix: static, and dynamic with kernels that differ from
    // pixel to pixel.
    std::vector<PixelConvolution> bothConvolutions(size_t W, size_t H)
    {
        return {staticConvolution(arrayCfg(kDistinct), W, H),
                dynamicConvolution(makeCube<float>(W, H, 9, pixelKernel9), W, H)};
    }
}

TEST(PixelConvolution, IsLinear)
{
    constexpr size_t W = 13, H = 9;
    for (const PixelConvolution& c : bothConvolutions(W, H))
    {
        SCOPED_TRACE(c.isDynamic() ? "dynamic" : "static");
        std::valarray<double> a = scene(W, H), b = point(W, H, 0, 0) + 5.0 * point(W, H, 6, 4), both = 2.0 * a + 3.0 * b;
        c.apply(a,    0, 0, W, H, 1);
        c.apply(b,    0, 0, W, H, 1);
        c.apply(both, 0, 0, W, H, 1);
        expectSame(both, 2.0 * a + 3.0 * b, W, 1e-12);
    }
}

TEST(PixelConvolution, ResultDoesNotDependOnTheNumberOfThreads)
{
    // Large enough for 7 threads (parallelFor() needs 100 pixels per thread), whose chunks start
    // and end in the middle of a row.
    constexpr size_t W = 61, H = 47;
    for (const PixelConvolution& c : bothConvolutions(W, H))
    {
        SCOPED_TRACE(c.isDynamic() ? "dynamic" : "static");
        std::valarray<double> single = scene(W, H);
        c.apply(single, 0, 0, W, H, 1);

        for (size_t nThreads : {0, 2, 3, 7})
        {
            SCOPED_TRACE(std::to_string(nThreads) + " threads");
            std::valarray<double> image = scene(W, H);
            c.apply(image, 0, 0, W, H, nThreads);
            expectSame(image, single, W, 1e-14);
        }
    }
}

TEST(PixelConvolution, FlatImageStaysFlatWithASingleKernel)
{
    // The baseline: same kernel everywhere, border included, also in windows smaller than the kernel.
    constexpr size_t W = 12, H = 10;
    const std::vector<double> kernel = normalized(kDistinct);
    const CubePtr cube = makeCube<float>(W, H, 9, [&kernel](size_t, size_t, size_t k) { return kernel[k]; });

    for (const PixelConvolution& c : {staticConvolution(arrayCfg(kDistinct), W, H), dynamicConvolution(cube, W, H)})
    {
        SCOPED_TRACE(c.isDynamic() ? "dynamic" : "static");
        for (const auto& [w, h] : std::vector<std::pair<size_t, size_t>>{{W, H}, {1, 1}, {2, 2}, {1, 5}, {5, 1}, {2, 7}})
        {
            SCOPED_TRACE(std::to_string(w) + " × " + std::to_string(h) + " window");
            std::valarray<double> image(1000.0, w * h);
            c.apply(image, 0, 0, w, h, 1);
            expectSame(image, std::valarray<double>(1000.0, w * h), w, c.isDynamic() ? 1e-3 : 1e-9);   // float weights: 1e-6 of 1000
        }
    }
}

TEST(PixelConvolution, ValuesAfterTheWindowAreUnchanged)
{
    constexpr size_t W = 7, H = 5;
    for (const PixelConvolution& c : bothConvolutions(W, H))
    {
        SCOPED_TRACE(c.isDynamic() ? "dynamic" : "static");
        std::valarray<double> exact = scene(W, H), longer(-1.0, W * H + 10);
        longer[std::slice(0, W * H, 1)] = exact;

        c.apply(exact,  0, 0, W, H, 1);
        c.apply(longer, 0, 0, W, H, 1);
        ASSERT_EQ(longer.size(), W * H + 10);
        expectSame(std::valarray<double>(longer[std::slice(0, W * H, 1)]), exact, W, 0.0);
        expectSame(std::valarray<double>(longer[std::slice(W * H, 10, 1)]), std::valarray<double>(-1.0, 10), 10, 0.0);
    }
}

TEST(PixelConvolution, InvalidWindowThrowsAndKeepsTheImage)
{
    constexpr size_t W = 7, H = 5;
    for (const PixelConvolution& c : bothConvolutions(W, H))
    {
        SCOPED_TRACE(c.isDynamic() ? "dynamic" : "static");
        const std::valarray<double> in = scene(W, H);
        std::valarray<double> image = in;

        EXPECT_THROW(c.apply(image, 1, 0, W, H, 1),     std::out_of_range);
        EXPECT_THROW(c.apply(image, 0, 1, W, H, 1),     std::out_of_range);
        EXPECT_THROW(c.apply(image, W, 0, 1, 1, 1),     std::out_of_range);
        EXPECT_THROW(c.apply(image, 0, 0, W + 1, 1, 1), std::out_of_range);
        expectSame(image, in, W, 0.0);

        std::valarray<double> tooSmall(0.0, W * H - 1);
        EXPECT_THROW(c.apply(tooSmall, 0, 0, W, H, 1), std::invalid_argument);

        EXPECT_NO_THROW(c.apply(image, W - 1, H - 1, 1, 1, 1));   // last pixel
    }
}

TEST(PixelConvolution, EmptyWindowDoesNothing)
{
    constexpr size_t W = 7, H = 5;
    for (const PixelConvolution& c : bothConvolutions(W, H))
    {
        SCOPED_TRACE(c.isDynamic() ? "dynamic" : "static");
        const std::valarray<double> in = scene(W, H);
        std::valarray<double> image = in, empty;

        EXPECT_NO_THROW(c.apply(image, 0, 0, 0, 0, 1));
        EXPECT_NO_THROW(c.apply(image, 3, 2, 0, 3, 1));
        EXPECT_NO_THROW(c.apply(image, 3, 2, 4, 0, 1));
        EXPECT_NO_THROW(c.apply(empty, 0, 0, 0, 0, 1));
        expectSame(image, in, W, 0.0);
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  (d) setKernel()
// ═══════════════════════════════════════════════════════════════════════════════

namespace
{
    std::valarray<double> toArray(const std::vector<double>& v)
    {
        return std::valarray<double>(v.data(), v.size());
    }

    // Kernel of a crosstalk between channels of `columns` columns, all read in the same direction:
    // one row, 1 in the centre, coupling[d − 1] at ± d channels.
    std::valarray<double> crosstalkKernel(const std::vector<double>& coupling, size_t columns)
    {
        const size_t centre = coupling.size() * columns;
        std::valarray<double> weights(0.0, 2 * centre + 1);
        weights[centre] = 1.0;
        for (size_t d = 1; d <= coupling.size(); ++d)
            weights[centre - d * columns] = weights[centre + d * columns] = coupling[d - 1];
        return weights;
    }
}

TEST(PixelConvolution, SetKernelNormalizesByDefault)
{
    constexpr size_t W = 13, H = 9;
    PixelConvolution c;
    c.setKernel(toArray(kDistinct), 3, 3, W, H);
    expectStatic(c, normalized(kDistinct), 3, 3);

    // The same as the array of the configuration, border included.
    std::valarray<double> a = scene(W, H), b = a;
    c.apply(a, 0, 0, W, H, 1);
    staticConvolution(arrayCfg(kDistinct), W, H).apply(b, 0, 0, W, H, 1);
    expectSame(a, b, W, 0.0);
}

TEST(PixelConvolution, KernelNotNormalizedIsUsedAsItIs)
{
    PixelConvolution c;
    c.setKernel(toArray(kDistinct), 3, 3, 7, 5, false);
    expectStatic(c, kDistinct, 3, 3);

    std::valarray<double> image = point(7, 5, 3, 2);
    c.apply(image, 0, 0, 7, 5, 1);
    for (size_t k = 0; k < 9; ++k)
        EXPECT_NEAR(image[(1 + k / 3) * 7 + 2 + k % 3], kDistinct[k], 1e-13) << "weight " << k;
    EXPECT_NEAR(image.sum(), 45.0, 1e-12);
}

TEST(PixelConvolution, KernelNotNormalizedByHandOnOneRow)
{
    // Keeps 1, adds 0.2 to the left and 0.3 to the right; nothing is divided at the border.
    PixelConvolution c;
    c.setKernel({0.2, 1.0, 0.3}, 3, 1, 5, 1, false);

    std::valarray<double> image = {0.0, 1.0, 0.0, 2.0, 0.0};
    c.apply(image, 0, 0, 5, 1, 1);
    expectSame(image, {0.2, 1.0, 0.3 * 1.0 + 0.2 * 2.0, 2.0, 0.6}, 5, 1e-14);
}

TEST(PixelConvolution, CrosstalkGivesGhostsInTheOtherChannels)
{
    // 4 channels of 8 columns. The kernel (49 columns) is wider than the matrix (32).
    constexpr size_t W = 32, H = 5;
    PixelConvolution c;
    c.setKernel(crosstalkKernel({1.0e-2, 5.0e-3, 2.0e-3}, 8), 49, 1, W, H, false);

    // Source at position 3 of channel 1, row 2: ghosts at position 3 of channels 0, 2 and 3.
    std::valarray<double> image = 1000.0 * point(W, H, 8 + 3, 2);
    c.apply(image, 0, 0, W, H, 1);

    std::valarray<double> expected(0.0, W * H);
    expected[2 * W + 3]      = 10.0;     // distance 1
    expected[2 * W + 8 + 3]  = 1000.0;   // the source is unchanged
    expected[2 * W + 16 + 3] = 10.0;     // distance 1
    expected[2 * W + 24 + 3] = 5.0;      // distance 2
    expectSame(image, expected, W, 1e-12);

    // From the first channel: distances 1, 2 and 3.
    image = 1000.0 * point(W, H, 7, 0);
    c.apply(image, 0, 0, W, H, 1);
    expected = 0.0;
    expected[7]      = 1000.0;
    expected[8 + 7]  = 10.0;
    expected[16 + 7] = 5.0;
    expected[24 + 7] = 2.0;
    expectSame(image, expected, W, 1e-12);

    // A window of one channel has no ghost.
    const std::valarray<double> in = scene(8, H);
    std::valarray<double> window = in;
    c.apply(window, 8, 0, 8, H, 1);
    expectSame(window, in, 8, 0.0);
}

TEST(PixelConvolution, CrosstalkOfTheFullDetector)
{
    // 32 channels of 64 columns, a coupling that depends on the distance: each ghost is at its
    // place, with several threads as with one.
    constexpr size_t W = 2048, H = 16;
    std::vector<double> coupling(31);
    for (size_t d = 1; d <= coupling.size(); ++d)
        coupling[d - 1] = 1.0e-3 / static_cast<double>(d);

    PixelConvolution c;
    c.setKernel(crosstalkKernel(coupling, 64), 3969, 1, W, H, false);

    const size_t channel = 10, position = 17, row = 5;
    std::valarray<double> image = point(W, H, channel * 64 + position, row), threaded = image;
    c.apply(image,    0, 0, W, H, 1);
    c.apply(threaded, 0, 0, W, H, 8);

    std::valarray<double> expected(0.0, W * H);
    for (size_t ch = 0; ch < 32; ++ch)
    {
        const size_t d = (ch > channel) ? ch - channel : channel - ch;
        expected[row * W + ch * 64 + position] = (d == 0) ? 1.0 : coupling[d - 1];
    }
    expectSame(image,    expected, W, 1e-15);
    expectSame(threaded, expected, W, 1e-15);
}

TEST(PixelConvolution, InitializeAfterSetKernelIsNormalizedAgain)
{
    constexpr size_t W = 13, H = 9;
    PixelConvolution c;
    c.setKernel(toArray(kDistinct), 3, 3, W, H, false);

    // Flat image: unchanged by a normalized kernel only.
    const std::valarray<double> flat(1000.0, W * H);
    initialize(c, arrayCfg(kDistinct), W, H);
    std::valarray<double> image = flat;
    c.apply(image, 0, 0, W, H, 1);
    expectSame(image, flat, W, 1e-9);

    c.setKernel(toArray(kDistinct), 3, 3, W, H, false);
    const CubePtr cube = makeCube<float>(W, H, 9, pixelKernel9);
    initialize(c, fileCfg(), W, H, [cube]() { return cube; });
    expectDynamic(c, W, H, 9, pixelKernel9);

    c.setKernel(toArray(kDistinct), 3, 3, W, H, false);   // and back from the cube
    expectStatic(c, kDistinct, 3, 3);

    initialize(c, jsonCfg(nlohmann::json::array()), W, H);
    expectOff(c);
}

TEST(PixelConvolution, InvalidSetKernelThrowsAndKeepsThePreviousOne)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    PixelConvolution c = staticConvolution(arrayCfg(kDistinct), 8, 6);
    for (bool normalize : {true, false})
    {
        auto expectThrow = [&c, normalize](const std::vector<double>& weights, size_t w, size_t h, const std::string& what)
        {
            SCOPED_TRACE(what + (normalize ? ", normalized" : ", not normalized"));
            EXPECT_THROW(c.setKernel(toArray(weights), w, h, 8, 6, normalize), std::invalid_argument);
            expectStatic(c, normalized(kDistinct), 3, 3);
        };

        expectThrow({1, 1},            2, 1, "even number of columns");
        expectThrow({1, 1},            1, 2, "even number of rows");
        expectThrow({},                0, 0, "no weight");
        expectThrow({1, 2, 3},         3, 3, "3 weights for 3 × 3");
        expectThrow(kDistinct,         9, 3, "9 weights for 9 × 3");
        expectThrow({0.1, 1.0, -0.1},  3, 1, "negative weight");
        expectThrow({0.1, nan, 0.1},   3, 1, "NaN weight");
        expectThrow({0.1, inf, 0.1},   3, 1, "infinite weight");
        expectThrow({0.0, 0.0, 0.0},   3, 1, "null kernel");
    }
}
