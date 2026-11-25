//
// Created by gillard on 2026-11-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// DetectorEffect.h
//
/// @file DetectorEffect.h
/// @brief Effects that couple the pixels of a detector: phot2pix::PixelConvolution.
///
/// An effect is a member of the detector that uses it (see nisp::SCA), which decides when to apply
/// it: it depends on the readout. Not included by phot2pix.h, because DSTmorpho.h brings
/// `using namespace DSL` at global scope.

#ifndef P2P_DETECTOR_EFFECT_H
#define P2P_DETECTOR_EFFECT_H

#include "Geometry.h"
#include "PixelConfig.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <valarray>

#include <DSTfits/FITSimg.h>
#include <DSTools/DSTmorpho.h>

namespace phot2pix
{

#pragma region - Pixel convolution
    /// @brief Convolution of an image by a kernel: each pixel gives a part of its signal to others.
    /// @details Used for the inter-pixel capacitance (IPC, read from the configuration by initialize())
    ///          and for the crosstalk between output channels (kernel given by the detector with
    ///          setKernel(), see nisp::SCA).
    ///
    ///          The kernel is the response to a signal in its central pixel: weights row by row
    ///          (columns fastest, first row y − 1), odd numbers of columns and rows.
    ///          - Static: the same for every pixel; normalized to a sum of 1, unless setKernel() says
    ///            otherwise.
    ///          - Dynamic: one kernel per pixel, from a FITS cube, normalized pixel by pixel.
    ///
    ///          A default-constructed convolution is disabled: apply() does nothing.
    class PixelConvolution
    {
    public:
        /// @brief Static kernel (columns, rows).
        using kernel = DST::Morphology::kernel<double>;
        /// @brief Returns the FITS image of "ipc"; null when the key does not name a file.
        using CubeLoader = std::function<std::shared_ptr<DSL::FITScube>()>;

        /// @brief Read the kernel from cfg "ipc". Nothing changes if the key is absent.
        /// @verbatim
        ///   "ipc": [k0, k1, ...]                                → static kernel, square (n × n weights)
        ///   "ipc": {"kernel": [...], "width": w, "height": h}   → static kernel, w columns × h rows
        ///   "ipc": "<path>" [+ "ipc.hdu": h]                    → dynamic kernel, FITS cube
        ///   "ipc": {"file": "<path>", "hdu": h}                 → dynamic kernel, FITS cube
        ///   "ipc": []                                           → no IPC
        /// @endverbatim
        /// @details - Static: the weights are >= 0.
        ///          - Dynamic: cube of n + 1 layers, the n neighbours and the central pixel of a square
        ///            kernel (9 layers for 3 × 3), layer k being weight k of every pixel. The middle
        ///            layer is the fraction kept by the central pixel. Matrix-sized or detector-sized
        ///            (see Geometry::extract()); the section of the matrix is kept, in float. The layers
        ///            of a pixel should not sum to more than 1; measurement errors can make them do so:
        ///            a warning is printed, not an error, and the kernel is normalized.
        /// @param cfg      Configuration.
        /// @param geometry Physical detector, and placement of the matrix on it.
        /// @param width    Matrix columns.
        /// @param height   Matrix rows.
        /// @param loadCube Loads the FITS image of "ipc", e.g. with DetectorSimulator::loadImage(); only
        ///                 called when "ipc" is not an array nor a {"kernel", …} object.
        /// @throws std::invalid_argument if the kernel is not valid; the IPC is then unchanged.
        void initialize(const PixelConfig& cfg, const Geometry& geometry, size_t width, size_t height,
                        const CubeLoader& loadCube);

        /// @brief Set a static kernel, e.g. built by the detector from its own configuration.
        /// @details The kernel may be wider than the matrix. A null weight costs little (a comparison per pixel,
        ///          see DST::Morphology::ops::pixel_correlate()), but a short kernel is still faster.
        /// @param weights   w × h weights >= 0, row by row (columns fastest, first row y − 1).
        /// @param w         Kernel columns, odd.
        /// @param h         Kernel rows, odd.
        /// @param width     Matrix columns.
        /// @param height    Matrix rows.
        /// @param normalize True (signal shared between pixels, as IPC): the weights are divided by their
        ///                  sum, and the border of a window is normalized, see apply(). False (signal added
        ///                  to other pixels, as crosstalk): the weights are used as they are, and what
        ///                  is outside the window counts as 0.
        /// @throws std::invalid_argument if the kernel is not valid; the convolution is then unchanged.
        void setKernel(const std::valarray<double>& weights, size_t w, size_t h, size_t width, size_t height,
                       bool normalize = true);

        /// @brief True if a kernel is configured.
        bool isEnabled() const { return m_enabled; }
        /// @brief True if the kernel depends on the pixel (FITS cube).
        bool isDynamic() const { return m_cube != nullptr; }
        /// @brief Static kernel, normalized to 1 unless setKernel() says otherwise; 1 × 1 without static kernel.
        const kernel& getKernel() const { return m_kernel; }
        /// @brief Dynamic kernel, width × height × (n + 1) float cube normalized pixel by pixel; null without
        ///        dynamic kernel.
        std::shared_ptr<const DSL::FITScube> getKernelCube() const { return m_cube; }

        /// @brief Convolve a window of the matrix, in place. Nothing is done without kernel, nor on an empty window.
        /// @details With a normalized kernel, charge is preserved, except at the border of the window: the
        ///          pixels that miss neighbours are divided by the sum of the weights of those they have.
        ///          With a kernel that is not normalized, missing neighbours count as 0.
        /// @param image    Window, ncols × nrows pixels, row-major (a longer array is left unchanged after them).
        /// @param c0       Matrix column of the first window pixel.
        /// @param r0       Matrix row of the first window pixel.
        /// @param ncols    Window columns.
        /// @param nrows    Window rows.
        /// @param nThreads Number of threads.
        /// @throws std::out_of_range if the window exceeds the matrix, std::invalid_argument if image is
        ///         smaller than the window.
        void apply(std::valarray<double>& image, size_t c0, size_t r0, size_t ncols, size_t nrows, size_t nThreads) const;

    private:
        kernel                         m_kernel = kernel::one({1, 1});   ///< Static kernel
        std::shared_ptr<DSL::FITScube> m_cube;                           ///< Dynamic kernel (float cube); null: static kernel
        bool                           m_enabled = false;                ///< A kernel is configured
        bool                           m_normalize = true;               ///< Normalized kernel (see setKernel())
        size_t                         m_width   = 0;                    ///< Matrix columns
        size_t                         m_height  = 0;                    ///< Matrix rows

        /// @brief apply() with the static kernel.
        void _applyStatic(std::valarray<double>& image, size_t ncols, size_t nrows, size_t nThreads) const;
        /// @brief apply() with the kernels of the cube.
        void _applyDynamic(std::valarray<double>& image, size_t c0, size_t r0, size_t ncols, size_t nrows, size_t nThreads) const;
    };
#pragma endregion

} // namespace phot2pix

#endif // P2P_DETECTOR_EFFECT_H
