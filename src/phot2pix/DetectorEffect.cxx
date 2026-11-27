//
// Created by gillard on 2026-11-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// DetectorEffect.cxx
//

#include "DetectorEffect.h"
#include "DetectorTools.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace phot2pix
{

#pragma region - Inter-pixel capacitance
    namespace
    {
        // Side of the square kernel of `count` weights. Throws unless it is an odd integer.
        size_t kernelSquareSide(size_t count)
        {
            const size_t side = static_cast<size_t>(std::lround(std::sqrt(static_cast<double>(count))));
            if (side * side != count || side % 2 == 0)
                throw std::invalid_argument("PixelConvolution: 'kernel': " + std::to_string(count) + " weights do not make a "
                                            "square kernel with an odd number of rows (1, 9, 25, ...).");
            return side;
        }

        // "width" or "height" of the kernel object j. Throws unless it is an odd integer.
        size_t kernelSize(const json& j, const char* name)
        {
            const double v = (j.contains(name) && j.at(name).is_number()) ? j.at(name).get<double>() : 0.0;
            if (!(v >= 1.0) || std::fmod(v, 2.0) != 1.0)
                throw std::invalid_argument(std::string("PixelConvolution: 'kernel': \"") + name + "\" must be an odd integer >= 1.");
            return static_cast<size_t>(v);
        }
    }

    void PixelConvolution::initialize(const PixelConfig& cfg, const Geometry& geometry, size_t width, size_t height,
                                           const CubeLoader& loadCube, const std::string& what)
    {
        if (!cfg.has(what) && !cfg.has(what+".file"))
            return;

        const std::string _what = "PixelConvolution: '" + what + "'";

        // The JSON value of "kernel"; null when it is a path (or a number, rejected below).
        const json* node = nullptr;
        if (cfg.has(what))
        {
            try                                { node = &cfg.getJson(what); }
            catch (const std::invalid_argument&) {}
        }

        // ── No Kernel: [] ──
        if (node && node->is_array() && node->empty())
        {
            m_kernel    = kernel::one({1, 1});
            m_cube      = nullptr;
            m_enabled   = false;
            m_normalize = true;
            return;
        }

        // ── Static kernel: [...] (square) or {"kernel": [...], "width": w, "height": h} ──
        if (node && (node->is_array() || (node->is_object() && node->contains("kernel"))))
        {
            const bool  isObject = node->is_object();
            const json& values   = isObject ? node->at("kernel") : *node;
            if (!values.is_array())
                throw std::invalid_argument(_what + ": \"kernel\" must be an array of numbers.");

            size_t w = 0, h = 0;
            if (isObject)
            {
                for (auto it = node->begin(); it != node->end(); ++it)
                    if (it.key() != "kernel" && it.key() != "width" && it.key() != "height")
                        throw std::invalid_argument(_what + ": unknown key '" + it.key() + "'.");

                w = kernelSize(*node, "width");
                h = kernelSize(*node, "height");
                if (w * h != values.size())
                    throw std::invalid_argument(_what + ": " + std::to_string(values.size()) + " weights for a " +
                                                std::to_string(w) + " × " + std::to_string(h) + " kernel.");
            }
            else
                w = h = kernelSquareSide(values.size());

            std::valarray<double> weights(0.0, values.size());
            for (size_t k = 0; k < values.size(); ++k)
            {
                if (!values[k].is_number())
                    throw std::invalid_argument(_what + ": weight " + std::to_string(k) + " is not a number.");
                weights[k] = values[k].get<double>();
            }

            try                                      { setKernel(weights, w, h, width, height); }
            catch (const std::invalid_argument& err) { throw std::invalid_argument(_what + ": " + err.what()); }
            return;
        }

        // ── Dynamic kernel: FITS cube of n + 1 layers (n neighbours and the central pixel) ──
        const std::shared_ptr<DSL::FITScube> image = loadCube();   // null: a number or a boolean
        if (!image)
            throw std::invalid_argument(_what + " must be an array of weights, {\"kernel\", \"width\", \"height\"}, "
                                        "or a FITS file (a path or {\"file\", \"hdu\"}).");
        if (image->GetDimension() != 3)
            throw std::invalid_argument(_what + ": the FITS kernel must be a cube (columns × rows × weights), not a " +
                                        std::to_string(image->GetDimension()) + "-D image.");

        const size_t nPix = width * height;
        const size_t mapW = image->Size(1), mapH = image->Size(2), nK = image->Size(3);
        kernelSquareSide(nK);   // throws if the layers do not make a square kernel

        auto cube = std::make_shared< DSL::FITSimg<float> >(std::vector<size_t>{width, height, nK});   // float: memory
        std::valarray<float>* out = cube->GetData<float>();
        if (!out)
            throw std::runtime_error(_what + ": cannot allocate the kernel cube.");

        // Section of this matrix, layer by layer, and sum of the layers of every pixel. The
        // middle layer is the fraction kept by the central pixel: it is summed like the others.
        std::valarray<double> sum(0.0, nPix);
        const bool known = withCubeData(*image, [&](const auto& data)
        {
            const size_t nMap = mapW * mapH;
            std::valarray<double> layer(0.0, nMap);
            for (size_t k = 0; k < nK; ++k)
            {
                std::copy_n(&data[k * nMap], nMap, std::begin(layer));                // T → double
                const std::valarray<double> weight = geometry.extract(layer, mapW, mapH, width, height, _what);
                std::copy(std::begin(weight), std::end(weight), &(*out)[k * nPix]);   // double → float
                sum += weight;                                                        // pixel by pixel, in double
            }
        });
        if (!known)
            throw std::invalid_argument(_what + ": the FITS kernel is not an integer or floating point cube.");

        // The layers of a pixel should not sum to more than 1 (1e-6: float rounding). The errors
        // of a measured kernel can make them do so: a warning, not an error.
        size_t nBad = 0, nAbove = 0;
        double largest = 0.0;
        for (size_t i = 0; i < nPix; ++i)
        {
            if (!(sum[i] > 0.0) || !std::isfinite(sum[i]))   // also NaN weights
                ++nBad;
            else if (sum[i] > 1.0 + 1e-6)
            {
                ++nAbove;
                largest = std::max(largest, sum[i]);
            }
        }

        if (nBad > 0)
            throw std::invalid_argument(_what + ": the sum of the layers must be finite and > 0: " + std::to_string(nBad) +
                                        " of " + std::to_string(nPix) + " pixels do not comply.");
        if (nAbove > 0)
            std::cerr << _what << ": warning: the layers sum to more than 1 in " << nAbove << " of " << nPix
                      << " pixels (largest sum: " << largest << "); the kernels are normalized to 1." << std::endl;

        // Every layer divided by the sums: the kernel of each pixel is normalized to 1.
        std::valarray<float> norm(0.0f, nPix);
        std::copy(std::begin(sum), std::end(sum), std::begin(norm));   // double → float

        for (size_t k = 0; k < nK; ++k)
            (*out)[std::slice(k * nPix, nPix, 1)] /= norm;

        m_kernel    = kernel::one({1, 1});   // not used
        m_cube      = cube;
        m_enabled   = true;
        m_normalize = true;
        m_width     = width;
        m_height    = height;
    }

    void PixelConvolution::setKernel(const std::valarray<double>& weights, size_t w, size_t h, size_t width, size_t height,
                                     bool normalize)
    {
        if (w % 2 == 0 || h % 2 == 0)
            throw std::invalid_argument("PixelConvolution::setKernel: the kernel must have odd numbers of columns and rows, not " +
                                        std::to_string(w) + " × " + std::to_string(h) + ".");
        if (weights.size() != w * h)
            throw std::invalid_argument("PixelConvolution::setKernel: " + std::to_string(weights.size()) + " weights for a " +
                                        std::to_string(w) + " × " + std::to_string(h) + " kernel.");

        for (size_t k = 0; k < weights.size(); ++k)
            if (!(weights[k] >= 0.0) || !std::isfinite(weights[k]))
                throw std::invalid_argument("PixelConvolution::setKernel: weight " + std::to_string(k) + " must be finite and >= 0.");

        const double sum = weights.sum();
        if (!(sum > 0.0) || !std::isfinite(sum))
            throw std::invalid_argument("PixelConvolution::setKernel: the sum of the weights must be finite and > 0.");

        std::valarray<double> kept = weights;
        if (normalize)
            kept /= sum;

        m_kernel    = kernel(kept, std::vector<size_t>{w, h});   // columns, rows
        m_cube      = nullptr;
        m_enabled   = true;
        m_normalize = normalize;
        m_width     = width;
        m_height    = height;
    }

    void PixelConvolution::apply(std::valarray<double>& image, size_t c0, size_t r0, size_t ncols, size_t nrows, size_t nThreads) const
    {
        if (!m_enabled || ncols == 0 || nrows == 0)
            return;

        if (c0 + ncols > m_width || r0 + nrows > m_height)
            throw std::out_of_range("PixelConvolution::apply: the window exceeds the matrix.");
        if (image.size() < ncols * nrows)
            throw std::invalid_argument("PixelConvolution::apply: the image is smaller than the window.");
        if (ncols == 0 || nrows == 0)
            return;

        if (m_cube)
            _applyDynamic(image, c0, r0, ncols, nrows, nThreads);
        else
            _applyStatic(image, ncols, nrows, nThreads);   // same kernel for every pixel: c0 and r0 do not matter
    }

    void PixelConvolution::_applyStatic(std::valarray<double>& image, size_t ncols, size_t nrows, size_t nThreads) const
    {
        const std::vector<size_t>& kSize = m_kernel.sizes();   // columns, rows: odd (initialize())

        // pixel_correlate() does not flip the kernel: the convolution is the correlation
        // with the kernel reversed along both axes.
        std::valarray<double> weights = m_kernel.data();
        std::reverse(std::begin(weights), std::end(weights));

        const size_t         imgSize[2] = {ncols, nrows};
        const std::ptrdiff_t center[2]  = {static_cast<std::ptrdiff_t>(kSize[0] / 2), static_cast<std::ptrdiff_t>(kSize[1] / 2)};

        std::valarray<double> out(0.0, ncols * nrows);   // the convolution cannot be done in place

        // Pixels outside the window are ignored with a normalized kernel, and count as 0 otherwise.
        parallelFor(ncols * nrows, nThreads, [&](size_t begin, size_t end)
        {
            for (size_t i = begin; i < end; ++i)
                out[i] = DST::Morphology::ops::pixel_correlate<double, double>(&image[0], nullptr, &weights[0], kSize.data(), 2, weights.size(),
                                                                               center, imgSize, 2, i, m_normalize);
        });

        std::copy(std::begin(out), std::end(out), std::begin(image));
    }

    void PixelConvolution::_applyDynamic(std::valarray<double>& image, size_t c0, size_t r0, size_t ncols, size_t nrows, size_t nThreads) const
    {
        const std::valarray<float>* cube = m_cube->GetData<float>();
        if (!cube)
            throw std::runtime_error("PixelConvolution::apply: the kernel cube is not a float cube.");

        using idx = std::ptrdiff_t;
        const idx nK   = static_cast<idx>(m_cube->Size(3));
        const idx side = static_cast<idx>(kernelSquareSide(m_cube->Size(3)));
        const idx half = side / 2;
        const idx nc   = static_cast<idx>(ncols), nr = static_cast<idx>(nrows);

        // Weights of layer k on row y of the window, from its first column.
        const size_t width = m_width, nPix = m_width * m_height;
        const auto weights = [&](idx k, idx y) { return &(*cube)[static_cast<size_t>(k) * nPix + (r0 + static_cast<size_t>(y)) * width + c0]; };

        std::valarray<double> out(0.0, ncols * nrows);   // the convolution cannot be done in place

        parallelFor(ncols * nrows, nThreads, [&](size_t begin, size_t end)
        {
            for (idx y = static_cast<idx>(begin) / nc; y * nc < static_cast<idx>(end); ++y)
            {
                // Columns [xa, xb) of this row that belong to this thread.
                const idx xa  = std::max<idx>(0,  static_cast<idx>(begin) - y * nc);
                const idx xb  = std::min<idx>(nc, static_cast<idx>(end)   - y * nc);
                double*   dst = &out[y * nc];

                // Pixel q gives the fraction K_k(q) of its charge to pixel q + (dx, dy): pixel x of
                // this row receives it from q = (x − dx, y − dy), with the kernel of q. Layer by
                // layer, to read the cube in the order it is stored.
                for (idx k = 0; k < nK; ++k)
                {
                    const idx dx = k % side - half, dy = k / side - half;
                    const idx yq = y - dy;
                    if (yq < 0 || yq >= nr)
                        continue;

                    const double* src = &image[yq * nc];
                    const float*  w   = weights(k, yq);
                    const idx     x1  = std::min(xb, nc + dx);
                    for (idx x = std::max(xa, dx); x < x1; ++x)
                        dst[x] += src[x - dx] * w[x - dx];
                }

                // Border: some sources are outside the window. Divide by the sum of the weights of
                // the others, like _applyStatic() with a normalized kernel.
                const bool inside = (y >= half && y + half < nr);
                for (idx x = xa; x < xb; ++x)
                {
                    if (inside && x >= half && x + half < nc)
                        continue;

                    double norm = 0.0;
                    for (idx k = 0; k < nK; ++k)
                    {
                        const idx xq = x - (k % side - half), yq = y - (k / side - half);
                        if (xq >= 0 && xq < nc && yq >= 0 && yq < nr)
                            norm += weights(k, y)[x];
                    }
                    if (norm > 0.0)
                        dst[x] /= norm;
                }
            }
        });

        std::copy(std::begin(out), std::end(out), std::begin(image));
    }
#pragma endregion

} // namespace phot2pix
