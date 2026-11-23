//
// Created by gillard on 2026-06-05.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// DetectorSimulator.h
//
/// @file DetectorSimulator.h
/// @brief phot2pix::DetectorSimulator, the pixel array: geometry, per-pixel maps, threaded simulation.
#pragma once

#include "PixelSimulator.h"
#include "Geometry.h"

#include <algorithm>
#include <memory>
#include <optional>
#include <vector>
#include <span>
#include <ranges>
#include <version>      // __cpp_lib_mdspan: defined when the standard library has std::mdspan
#if defined(__cpp_lib_mdspan)
#include <mdspan>       // C++23; libstdc++ provides it from GCC 15 on (see StridedView below)
#endif
#include <type_traits>
#include <array>
#include <stdexcept>
#include <string>
#include <thread>       // std::thread, hardware_concurrency()
#include <valarray>
#include <DSTfits/FITSmanager.h>  // DSL::FITSmanager
#include <DSTfits/FITSimg.h>  // DSL::FITSmanager

namespace phot2pix
{
    /// @brief Write header cards (keyword, value, comment) into the header of a FITS image.
    /// @param image Image whose header receives the cards.
    /// @param cards Cards to write.
    inline void writeCards(DSL::FITScube& image, const std::vector<HeaderCard>& cards)
    {
        for (const auto& c : cards)
            std::visit([&image, &c](const auto& v) { image.HDU().ValueForKey(c.key, v, c.comment); }, c.value);
    }

#pragma region - StridedView: std::mdspan, or a fallback without <mdspan>
    // View of Rank dimensions on a flat array, with a stride per dimension: the columns and
    // windows of the pixel array (DetectorSimulator::Col(), window()). Built by makeStridedView().
    // This StridedView as been designed by CLAUDE OPUS 5.5 in order to make the code compatible with C++23 and C++20. It is a fallback for the standard libraries without <mdspan> (libstdc++ before GCC 15, used by clang on Linux as well): the part of std::mdspan<T, std::dextents<size_t, Rank>, std::layout_stride> listed below, under the same names. Code that keeps to these members compiles with either.
#if defined(__cpp_lib_mdspan)
    /// @brief View of Rank dimensions on a flat array, with a stride per dimension (std::mdspan).
    template<typename T, size_t Rank>
    using StridedView = std::mdspan<T, std::dextents<size_t, Rank>, std::layout_stride>;

    /// @brief Build a StridedView.
    /// @param data    First element of the view.
    /// @param extents Number of elements along each dimension.
    /// @param strides Distance between two elements along each dimension.
    template<typename T, size_t Rank>
    StridedView<T, Rank> makeStridedView(T* data, const std::array<size_t, Rank>& extents,
                                         const std::array<size_t, Rank>& strides)
    {
        using ext_t = std::dextents<size_t, Rank>;
        return StridedView<T, Rank>(data, std::layout_stride::mapping<ext_t>(ext_t(extents), strides));
    }
#else
    /// @brief View of Rank dimensions on a flat array, with a stride per dimension (fallback).
    /// @details Fallback for the standard libraries without &lt;mdspan&gt; (libstdc++ before GCC 15, used by
    /// clang on Linux as well): the part of std::mdspan<T, std::dextents<size_t, Rank>,
    /// std::layout_stride> listed below, under the same names. Code that keeps to these members
    /// compiles with either.
    template<typename T, size_t Rank>
    class StridedView
    {
    public:
        using element_type     = T;                      ///< Element type
        using value_type       = std::remove_cv_t<T>;    ///< Element type without const / volatile
        using index_type       = size_t;                 ///< Index type
        using size_type        = size_t;                 ///< Size type
        using data_handle_type = T*;                     ///< Pointer to the elements
        using reference        = T&;                     ///< Reference to an element

        /// @brief View on data (see makeStridedView()).
        StridedView(T* data, const std::array<size_t, Rank>& extents, const std::array<size_t, Rank>& strides)
            : m_data(data), m_extents(extents), m_strides(strides) {}

        /// @brief Number of dimensions.
        static constexpr size_t rank() noexcept { return Rank; }
        /// @brief Number of elements along dimension r.
        size_t extent(size_t r) const noexcept { return m_extents[r]; }
        /// @brief Distance between two elements along dimension r.
        size_t stride(size_t r) const noexcept { return m_strides[r]; }
        /// @brief First element of the view.
        T*     data_handle()    const noexcept { return m_data; }

        /// @brief Number of elements of the view.
        size_t size() const noexcept
        {
            size_t n = 1;
            for (size_t e : m_extents) n *= e;
            return n;
        }
        /// @brief True if the view has no element.
        bool empty() const noexcept { return size() == 0; }

        /// @brief view[i] (rank 1), view[row, col] (rank 2): one index per dimension, no bounds check.
        template<typename... I>
            requires (sizeof...(I) == Rank && (std::is_convertible_v<I, size_t> && ...))
        T& operator[](I... indices) const
        {
            const std::array<size_t, Rank> idx{ static_cast<size_t>(indices)... };
            size_t offset = 0;
            for (size_t r = 0; r < Rank; ++r)
                offset += idx[r] * m_strides[r];
            return m_data[offset];
        }

    private:
        T*                       m_data;      ///< First element
        std::array<size_t, Rank> m_extents;   ///< Elements along each dimension
        std::array<size_t, Rank> m_strides;   ///< Distance between two elements along each dimension
    };

    /// @brief Build a StridedView.
    /// @param data    First element of the view.
    /// @param extents Number of elements along each dimension.
    /// @param strides Distance between two elements along each dimension.
    template<typename T, size_t Rank>
    StridedView<T, Rank> makeStridedView(T* data, const std::array<size_t, Rank>& extents,
                                         const std::array<size_t, Rank>& strides)
    {
        return StridedView<T, Rank>(data, extents, strides);
    }
#endif
#pragma endregion

#pragma region - DetectorSimulator class definition
    /**
     * @class DetectorSimulator
     * @brief Detector: an array of pixels, its geometry, its per-pixel maps and its simulation.
     * @details Flat row-major pixel array with 2D access.
     *          Non-const pixel accessors are protected so only derived classes can
     *          modify individual pixels. Public API exposes const (read) access and
     *          batch configuration methods (gain map, dark map, ...) that operate on
     *          the whole detector coherently.
     *
     *          Parallelism is handled internally by simulate() via std::thread.
     *          Derived classes customise physics by overriding _engine(), not simulate().
     * @tparam PixelT Pixel model, derived from pixsim.
     */
    template<typename PixelT = pixsim>
    class DetectorSimulator
    {
        static_assert(std::is_base_of_v<pixsim, PixelT>, "PixelT must derive from pixsim");

#pragma region - private
// ── Private: thread dispatch ───────────────────────────────────────────────
    private:

        /// @brief The whole pixel array.
        auto view()       { return std::span<PixelT>      (m_pixels.data(), m_pixels.size()); }
        /// @brief The whole pixel array, read only.
        auto view() const { return std::span<const PixelT>(m_pixels.data(), m_pixels.size()); }

        /// @brief pixelresponse() with cfg, or without it if PixelT has no such overload.
        /// @details Falls back to the 2-arg form when PixelT does not expose the 3-arg overload
        ///          (e.g. because a derived class overrides only the 2-arg form, hiding the
        ///          base-class 3-arg via C++ name hiding).
        template<typename P>
        static double _dispatchPixelResponse(const P& px, double flux, double dt, const PixelConfig& cfg)
        {
            if constexpr (requires { px.pixelresponse(flux, dt, cfg); })
                return px.pixelresponse(flux, dt, cfg);
            else
                return px.pixelresponse(flux, dt);
        }

        /// @brief Check that a map has one value per pixel.
        /// @throws std::invalid_argument otherwise, naming the map.
        void _checkMapSize(size_t sz, const char* name) const
        {
            if (sz != m_pixels.size())
                throw std::invalid_argument(std::string(name) + " size does not match detector size.");
        }

        /// @brief Number of threads to use for n pixels.
        /// @details Below 100 * nCores the thread-spawn overhead exceeds the gain: 1.
        size_t _computeThreadCount(size_t n) const
        {
            const size_t maxT = std::thread::hardware_concurrency();
            return (n > 100 * maxT) ? std::min(maxT, n) : 1;
        }

        /// @brief Full-detector dispatch: partition the flat pixel index range across threads.
        void _dispatch_cpu(const std::valarray<double>& flux,
                           double                        dt,
                           std::valarray<double>&         frame,
                           size_t                         nThreads) const
        {
            const size_t n    = m_pixels.size();
            const size_t base = n / nThreads;
            const size_t rem  = n % nThreads;

            std::vector<std::thread> workers;
            workers.reserve(nThreads);

            for (size_t t = 0; t < nThreads; ++t)
            {
                const size_t start = t * base + std::min(t, rem);
                const size_t end   = start + base + (t < rem ? 1 : 0);
                workers.emplace_back([this, &flux, dt, start, end, &frame]()
                {
                    _engine(flux, dt, start, end, frame); // one virtual call per thread
                });
            }
            for (auto& th : workers) th.join();
        }

        /// @brief Windowed dispatch: partition the window rows across threads.
        void _dispatch_cpu_window(const std::valarray<double>& flux,
                                   double                        dt,
                                   size_t c0, size_t r0,
                                   size_t ncols, size_t nrows,
                                   std::valarray<double>&         frame,
                                   size_t                         nThreads) const
        {
            const size_t base = nrows / nThreads;
            const size_t rem  = nrows % nThreads;

            std::vector<std::thread> workers;
            workers.reserve(nThreads);

            for (size_t t = 0; t < nThreads; ++t)
            {
                const size_t row_start = t * base + std::min(t, rem);
                const size_t row_end   = row_start + base + (t < rem ? 1 : 0);
                workers.emplace_back([this, &flux, dt, c0, r0, row_start, row_end, ncols, &frame]()
                {
                    _engine(flux, dt, c0, r0, row_start, row_end, ncols, frame);
                });
            }
            for (auto& th : workers) th.join();
        }

        /// @brief Full-detector dispatch with cfg, captured by const-ref in the thread lambda.
        void _dispatch_cpu_cfg(const std::valarray<double>& flux,
                               double                        dt,
                               std::valarray<double>&         frame,
                               size_t                         nThreads,
                               const PixelConfig&             cfg) const
        {
            const size_t n    = m_pixels.size();
            const size_t base = n / nThreads;
            const size_t rem  = n % nThreads;

            std::vector<std::thread> workers;
            workers.reserve(nThreads);

            for (size_t t = 0; t < nThreads; ++t)
            {
                const size_t start = t * base + std::min(t, rem);
                const size_t end   = start + base + (t < rem ? 1 : 0);
                workers.emplace_back([this, &flux, dt, start, end, &frame, &cfg]()
                {
                    _engine(flux, dt, start, end, frame, cfg);
                });
            }
            for (auto& th : workers) th.join();
        }

        /// @brief Windowed dispatch with cfg, captured by const-ref in the thread lambda.
        void _dispatch_cpu_window_cfg(const std::valarray<double>& flux,
                                       double                        dt,
                                       size_t c0, size_t r0,
                                       size_t ncols, size_t nrows,
                                       std::valarray<double>&         frame,
                                       size_t                         nThreads,
                                       const PixelConfig&             cfg) const
        {
            const size_t base = nrows / nThreads;
            const size_t rem  = nrows % nThreads;

            std::vector<std::thread> workers;
            workers.reserve(nThreads);

            for (size_t t = 0; t < nThreads; ++t)
            {
                const size_t row_start = t * base + std::min(t, rem);
                const size_t row_end   = row_start + base + (t < rem ? 1 : 0);
                workers.emplace_back([this, &flux, dt, c0, r0, row_start, row_end, ncols, &frame, &cfg]()
                {
                    _engine(flux, dt, c0, r0, row_start, row_end, ncols, frame, cfg);
                });
            }
            for (auto& th : workers) th.join();
        }

        /// @brief FITS image of a simulated frame, in the pixel type of the detector.
        /// @details Factored out to avoid duplication between the standard and cfg-aware
        ///          simulate() overloads.
        /// @param frame  Values, row-major: width × height of them (the window size for a window).
        /// @param width  Image columns.
        /// @param height Image rows.
        /// @throws std::invalid_argument if the pixel type is not a FITS type.
        std::shared_ptr<DSL::FITScube> _convertFrame(const std::valarray<double>& frame,
                                                      size_t width, size_t height) const
        {
            std::shared_ptr<DSL::FITScube> img;
            switch (m_pixelType)
            {
                case SBYTE_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<int8_t>({width, height}));
                    *(img) += convert2int8(frame);   break;
                case BYTE_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<uint8_t>({width, height}));
                    *(img) += convert2uint8(frame);  break;
                case SHORT_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<int16_t>({width, height}));
                    *(img) += convert2int16(frame);  break;
                case USHORT_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<uint16_t>({width, height}));
                    *(img) += convert2uint16(frame); break;
                case LONG_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<int32_t>({width, height}));
                    *(img) += convert2int32(frame);  break;
                case ULONG_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<uint32_t>({width, height}));
                    *(img) += convert2uint32(frame); break;
                case LONGLONG_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<int64_t>({width, height}));
                    *(img) += convert2int64(frame);  break;
                case ULONGLONG_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<uint64_t>({width, height}));
                    *(img) += convert2uint64(frame); break;
                case FLOAT_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<float>({width, height}));
                    *(img) += convert2float(frame);  break;
                case DOUBLE_IMG:
                    img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<double>({width, height}));
                    *(img) += frame;                 break;
                default:
                    throw std::invalid_argument("Unrecognised FITS data type for output.");
            }
            return img;
        }
#pragma endregion

// ── Protected: data, pixel write access, and virtual physics hook ──────────
#pragma region - protected
    protected:
#pragma region -> members
        const size_t        m_width;  ///< Matrix columns
        const size_t        m_height; ///< Matrix rows
        std::vector<PixelT> m_pixels; ///< Flat row-major pixel array
        const int m_pixelType;        ///< FITS storage type for output cube
        size_t m_nThreads;            ///< Thread budget for simulate(); set via setThreadCount()
        std::vector<HeaderCard> m_inputCards;   ///< Simulation inputs as header cards (see getInputCards())

        /// Physical detector and placement of the matrix on it (see Geometry.h). Replaced by
        /// initialize(), never modified in place, so copies of a detector may share it.
        std::shared_ptr<const Geometry> m_geometry;

#pragma endregion
#pragma region -> member functions

        /// @brief Add an input card, or replace the card that has the same keyword.
        void setInputCard(const HeaderCard& card)
        {
            for (auto& c : m_inputCards)
                if (c.key == card.key) { c = card; return; }
            m_inputCards.push_back(card);
        }

        /// @brief True if there is an input card with this keyword.
        bool hasInputCard(const std::string& key) const
        {
            return std::any_of(m_inputCards.begin(), m_inputCards.end(),
                               [&key](const HeaderCard& c) { return c.key == key; });
        }

        /// @brief Remove the input cards whose keyword starts with prefix (before re-adding a list).
        void removeInputCards(const std::string& prefix)
        {
            std::erase_if(m_inputCards, [&prefix](const HeaderCard& c) { return c.key.rfind(prefix, 0) == 0; });
        }

        /// @brief Record the card of a per-pixel input.
        /// @details The value or map file cfg gives for `key`; a key that was never given gets
        ///          the value the pixels hold (the PixelT default).
        /// @param cfg     Configuration.
        /// @param key     Configuration key.
        /// @param fitsKey Keyword of the card.
        /// @param comment Comment of the card.
        /// @param get     Reads the parameter of a pixel.
        template<typename Getter>
        void _recordInput(const PixelConfig& cfg, const std::string& key, const std::string& fitsKey,
                          const std::string& comment, Getter get)
        {
            if (const auto value = cfg.card(key))
                setInputCard({fitsKey, *value, comment});
            else if (!hasInputCard(fitsKey) && !m_pixels.empty())
                setInputCard({fitsKey, static_cast<double>(get(m_pixels.front())), comment});
        }

        // Non-const pixel accessors — derived classes only.
        // Use these during initialisation (defect maps, special pixel regions, etc.).

        /// @brief Pixel of flat index index (row-major). No bounds check.
        PixelT& at(size_t index)           { return view()[index]; }
        /// @brief Pixel at column col and row row. No bounds check.
        PixelT& at(size_t col, size_t row) { return view()[col + row * m_width]; }

        /// @brief Pixels of one row.
        std::span<PixelT> Row(size_t row)
        {
            return { m_pixels.data() + row * m_width, m_width };
        }

        /// @brief Pixels of one column.
        StridedView<PixelT, 1> Col(size_t col)
        {
            return makeStridedView<PixelT, 1>(m_pixels.data() + col, { m_height }, { m_width });
        }

        /// @brief Pixels of a window, indexed [row, col] from its corner (c0, r0).
        StridedView<PixelT, 2> window(size_t c0, size_t r0, size_t ncols, size_t nrows)
        {
            return makeStridedView<PixelT, 2>(m_pixels.data() + r0 * m_width + c0,
                                              { nrows, ncols }, { m_width, 1 });
        }

        // ── Virtual physics hooks ────────────────────────────────────────────
        //
        // Derived classes override _engine() to add detector-level effects
        // (cross-talk, charge diffusion, persistence, column glow, …).
        //
        // Per-pixel physics (noise, dark current, QE, gain, …) belong in
        // pixsim::pixelresponse() — NOT here.
        //
        // The default implementation calls pixsim::pixelresponse() for each pixel.
        // Each thread calls _engine() exactly once with its assigned range,
        // so virtual dispatch overhead is O(nThreads), not O(nPixels).

        /// @brief Full-detector hook: the thread owns pixels [begin, end).
        /// @param flux  Photons/s of every pixel of the matrix, row-major.
        /// @param dt    Integration time (s).
        /// @param begin First pixel.
        /// @param end   Pixel after the last one.
        /// @param frame Output, one value per pixel of the matrix.
        virtual void _engine(const std::valarray<double>& flux,
                              double                        dt,
                              size_t                        begin,
                              size_t                        end,
                              std::valarray<double>&         frame) const
        {
            for (size_t i = begin; i < end; ++i)
                frame[i] = m_pixels[i].pixelresponse(flux[i], dt);
        }

        /// @brief Windowed hook: the thread owns the window rows [row_begin, row_end).
        /// @param flux      Photons/s of the window pixels (ncols * nrows elements, row-major).
        /// @param dt        Integration time (s).
        /// @param c0        Matrix column of the first window pixel.
        /// @param r0        Matrix row of the first window pixel.
        /// @param row_begin First window row.
        /// @param row_end   Window row after the last one.
        /// @param ncols     Window columns.
        /// @param frame     Output, window-sized like flux.
        virtual void _engine(const std::valarray<double>& flux,
                              double                        dt,
                              size_t                        c0,
                              size_t                        r0,
                              size_t                        row_begin,
                              size_t                        row_end,
                              size_t                        ncols,
                              std::valarray<double>&         frame) const
        {
            for (size_t r = row_begin; r < row_end; ++r)
            {
                const size_t det_base   = (r0 + r) * m_width + c0;
                const size_t frame_base = r * ncols;
                for (size_t c = 0; c < ncols; ++c)
                    frame[frame_base + c] = m_pixels[det_base + c].pixelresponse(flux[frame_base + c], dt);
            }
        }

        // cfg-aware engine hooks — pass PixelConfig to pixelresponse() for runtime
        // parameters (e.g. persistence dT).  Derived classes that need per-engine
        // cfg logic can override these; the defaults delegate per-pixel dispatch.

        /// @brief Full-detector hook with runtime parameters.
        /// @param flux  Photons/s of every pixel of the matrix, row-major.
        /// @param dt    Integration time (s).
        /// @param begin First pixel.
        /// @param end   Pixel after the last one.
        /// @param frame Output, one value per pixel of the matrix.
        /// @param cfg   Runtime parameters.
        virtual void _engine(const std::valarray<double>& flux,
                              double                        dt,
                              size_t                        begin,
                              size_t                        end,
                              std::valarray<double>&        frame,
                              const PixelConfig&            cfg) const
        {
            for (size_t i = begin; i < end; ++i)
                frame[i] = _dispatchPixelResponse(m_pixels[i], flux[i], dt, cfg);
        }

        /// @brief Windowed hook with runtime parameters.
        /// @param flux      Photons/s of the window pixels (ncols * nrows elements, row-major).
        /// @param dt        Integration time (s).
        /// @param c0        Matrix column of the first window pixel.
        /// @param r0        Matrix row of the first window pixel.
        /// @param row_begin First window row.
        /// @param row_end   Window row after the last one.
        /// @param ncols     Window columns.
        /// @param frame     Output, window-sized like flux.
        /// @param cfg       Runtime parameters.
        virtual void _engine(const std::valarray<double>& flux,
                              double                        dt,
                              size_t                        c0,
                              size_t                        r0,
                              size_t                        row_begin,
                              size_t                        row_end,
                              size_t                        ncols,
                              std::valarray<double>&         frame,
                              const PixelConfig&             cfg) const
        {
            for (size_t r = row_begin; r < row_end; ++r)
            {
                const size_t det_base   = (r0 + r) * m_width + c0;
                const size_t frame_base = r * ncols;
                for (size_t c = 0; c < ncols; ++c)
                    frame[frame_base + c] = _dispatchPixelResponse(m_pixels[det_base + c], flux[frame_base + c], dt, cfg);
            }
        }

#pragma endregion
#pragma region -> array format convertion
        // Format convertion functions for outputting the detector response to FITS files.
        // They are marked as virtual so that derived classes can override them to implement custom clipping or scaling.

        /// @brief Convert a valarray of doubles to a valarray of int8_t, with clipping to the range of int8_t.
        /// @param arr The input valarray of doubles to convert.
        /// @return A valarray of int8_t containing the converted values.
        virtual std::valarray<int8_t> convert2int8(const std::valarray<double>& arr) const
        {
            std::valarray<int8_t> out((int8_t) 0, arr.size());
            std::transform(std::begin(arr), std::end(arr), std::begin(out), [](double x) 
                {
                    if( (x+0.5) < static_cast<double>(std::numeric_limits<int8_t>::min()))
                        return std::numeric_limits<int8_t>::min();
                    else if( (x+0.5) > static_cast<double>(std::numeric_limits<int8_t>::max()))
                        return std::numeric_limits<int8_t>::max();
                    else
                        return static_cast<int8_t>(x+0.5);
                }
            );
                
            return out;
        }

        /// @brief Convert a valarray of doubles to a valarray of uint8_t, with clipping to the range of uint8_t.
        /// @param arr The input valarray of doubles to convert.
        /// @return A valarray of uint8_t containing the converted values.
        virtual std::valarray<uint8_t> convert2uint8(const std::valarray<double>& arr) const
        {
            std::valarray<uint8_t> out((uint8_t) 0, arr.size());
            std::transform(std::begin(arr), std::end(arr), std::begin(out), [](double x) 
                {
                    if( (x+0.5) < 0.0)
                        return (uint8_t) 0;
                    else if( (x+0.5) > static_cast<double>(std::numeric_limits<uint8_t>::max()))
                        return std::numeric_limits<uint8_t>::max();
                    else
                        return static_cast<uint8_t>(x+0.5);
                }
            );
            
            return out;
        }

        /// @brief Convert a valarray of doubles to a valarray of int16_t, with clipping to the range of int16_t.
        /// @param arr The input valarray of doubles to convert.
        /// @return A valarray of int16_t containing the converted values.
        virtual std::valarray<int16_t> convert2int16(const std::valarray<double>& arr) const
        {
            std::valarray<int16_t> out((int16_t) 0, arr.size());
            std::transform(std::begin(arr), std::end(arr), std::begin(out), [](double x) 
                {
                    if( (x+0.5) < static_cast<double>(std::numeric_limits<int16_t>::min()))
                        return std::numeric_limits<int16_t>::min();
                    else if( (x+0.5) > static_cast<double>(std::numeric_limits<int16_t>::max()))
                        return std::numeric_limits<int16_t>::max();
                    else
                        return static_cast<int16_t>(x+0.5);
                }
            );
                
            return out;
        }

        /// @brief Convert a valarray of doubles to a valarray of uint16_t, with clipping to the range of uint16_t.
        /// @param arr The input valarray of doubles to convert.
        /// @return A valarray of uint16_t containing the converted values.
        virtual std::valarray<uint16_t> convert2uint16(const std::valarray<double>& arr) const
        {
            std::valarray<uint16_t> out((uint16_t) 0, arr.size());
            std::transform(std::begin(arr), std::end(arr), std::begin(out), [](double x) 
                {
                    if( (x+0.5) < 0.0)
                        return (uint16_t) 0;
                    else if( (x+0.5) > static_cast<double>(std::numeric_limits<uint16_t>::max()))
                        return std::numeric_limits<uint16_t>::max();
                    else
                        return static_cast<uint16_t>(x+0.5);
                }
            );
            
            return out;
        }

        /// @brief Convert a valarray of doubles to a valarray of int32_t, with clipping to the range of int32_t.
        /// @param arr The input valarray of doubles to convert.
        /// @return A valarray of int32_t containing the converted values.
        virtual std::valarray<int32_t> convert2int32(const std::valarray<double>& arr) const
        {
            std::valarray<int32_t> out((int32_t) 0, arr.size());
            std::transform(std::begin(arr), std::end(arr), std::begin(out), [](double x) 
                {
                    if( (x+0.5) < static_cast<double>(std::numeric_limits<int32_t>::min()))
                        return std::numeric_limits<int32_t>::min();
                    else if( (x+0.5) > static_cast<double>(std::numeric_limits<int32_t>::max()))
                        return std::numeric_limits<int32_t>::max();
                    else
                        return static_cast<int32_t>(x+0.5);
                }
            );
                
            return out;
        }

        /// @brief Convert a valarray of doubles to a valarray of uint32_t, with clipping to the range of uint32_t.
        /// @param arr The input valarray of doubles to convert.
        /// @return A valarray of uint32_t containing the converted values.
        virtual std::valarray<uint32_t> convert2uint32(const std::valarray<double>& arr) const
        {
            std::valarray<uint32_t> out((uint32_t) 0, arr.size());
            std::transform(std::begin(arr), std::end(arr), std::begin(out), [](double x) 
                {
                    if( (x+0.5) < 0.0)
                        return (uint32_t) 0;
                    else if( (x+0.5) > static_cast<double>(std::numeric_limits<uint32_t>::max()))
                        return std::numeric_limits<uint32_t>::max();
                    else
                        return static_cast<uint32_t>(x+0.5);
                }
            );
            
            return out;
        }

        /// @brief Convert a valarray of doubles to a valarray of int64_t, with clipping to the range of int64_t.
        /// @param arr The input valarray of doubles to convert.
        /// @return A valarray of int64_t containing the converted values.
        virtual std::valarray<int64_t> convert2int64(const std::valarray<double>& arr) const
        {
            std::valarray<int64_t> out((int64_t) 0, arr.size());
            std::transform(std::begin(arr), std::end(arr), std::begin(out), [](double x) 
                {
                    if( (x+0.5) < static_cast<double>(std::numeric_limits<int64_t>::min()))
                        return std::numeric_limits<int64_t>::min();
                    else if( (x+0.5) > static_cast<double>(std::numeric_limits<int64_t>::max()))
                        return std::numeric_limits<int64_t>::max();
                    else
                        return static_cast<int64_t>(x+0.5);
                }
            );
                
            return out;
        }

        /// @brief Convert a valarray of doubles to a valarray of uint64_t, with clipping to the range of uint64_t.
        /// @param arr The input valarray of doubles to convert.
        /// @return A valarray of uint64_t containing the converted values.
        virtual std::valarray<uint64_t> convert2uint64(const std::valarray<double>& arr) const
        {
            std::valarray<uint64_t> out((uint64_t) 0, arr.size());
            std::transform(std::begin(arr), std::end(arr), std::begin(out), [](double x) 
                {
                    if( (x+0.5) < 0.0)
                        return (uint64_t) 0;
                    else if( (x+0.5) > static_cast<double>(std::numeric_limits<uint64_t>::max()))
                        return std::numeric_limits<uint64_t>::max();
                    else
                        return static_cast<uint64_t>(x+0.5);
                }
            );
            
            return out;
        }

        /// @brief Convert a valarray of doubles to a valarray of float, without clipping.
        /// @param arr The input valarray of doubles to convert.
        /// @return A valarray of float containing the converted values.
        virtual std::valarray<float> convert2float(const std::valarray<double>& arr) const
        {
            std::valarray<float> out(0.0f, arr.size());
            std::transform(std::begin(arr), std::end(arr), std::begin(out), [](double x) { return static_cast<float>(x); });
            return out;
        }

#pragma endregion
#pragma endregion

// ── Public API ─────────────────────────────────────────────────────────────
#pragma region - public
    public:

#pragma region - ctor/dtor
        /// @brief Matrix of width × height default pixels.
        /// @param width    Matrix columns, > 0.
        /// @param height   Matrix rows, > 0.
        /// @param dataType FITS pixel type of the images simulate() returns (DOUBLE_IMG, USHORT_IMG, …).
        /// @param geometry Physical detector the matrix belongs to. A derived Geometry adds its own
        ///                 layout (reference pixels, …) and configuration keys. Null → a detector of
        ///                 the matrix size.
        /// @throws std::invalid_argument if width or height is 0.
        DetectorSimulator(size_t width, size_t height, const int& dataType = DOUBLE_IMG,
                          std::shared_ptr<const Geometry> geometry = nullptr)
            : m_width(width), m_height(height), m_pixels(width * height), m_pixelType(dataType),
              m_nThreads(_computeThreadCount(width * height)), m_geometry(std::move(geometry))
        {
            if (width == 0 || height == 0)
                throw std::invalid_argument("Width and height must be positive integers.");
            if (!m_geometry)
                m_geometry = std::make_shared<Geometry>(width, height);
        }

        /// @brief Destructor.
        virtual ~DetectorSimulator() = default;
#pragma endregion

#pragma region - geometry
        /// @brief Matrix columns.
        size_t getWidth()  const { return m_width; }
        /// @brief Matrix rows.
        size_t getHeight() const { return m_height; }
        /// @brief Number of pixels.
        size_t size()      const { return m_pixels.size(); }

        /// @brief Physical detector and placement of the matrix on it (see Geometry.h).
        const Geometry& getGeometry() const { return *m_geometry; }

        /// @brief Position (x0, y0) of the matrix on the physical detector.
        /// @details "detectorOrigin" if given, otherwise (0, 0). Matrix pixel (x, y) is physical
        ///          pixel (x0 + x, y0 + y).
        std::pair<size_t, size_t> getOrigin() const { return m_geometry->place(m_width, m_height); }

        /// @brief Per-pixel values of the matrix from a 2-D image.
        /// @details A matrix-sized image is used as is, a detector-sized one is cropped at getOrigin().
        /// @param image Image.
        /// @param what  Name of the image in the error messages.
        /// @return One value per pixel, row-major.
        /// @throws std::invalid_argument for any other size, or an image without data.
        std::valarray<double> toMatrix(const DSL::FITScube& image, const std::string& what) const
        {
            std::shared_ptr<DSL::FITScube> doubleImage;
            if (image.GetEqBitPerPixel() != DOUBLE_IMG || image.GetBitPerPixel() != DOUBLE_IMG)
                doubleImage = image.cast<double>();

            const DSL::FITScube* img = doubleImage ? doubleImage.get() : &image;
            const std::valarray<double>* data = img->GetData<double>();
            if (!data)
                throw std::invalid_argument(what + ": no pixel data.");

            return m_geometry->extract(*data, img->Size(1), img->Size(2), m_width, m_height, what);
        }

        /// @brief Thread budget used by simulate().
        size_t getThreadCount() const { return m_nThreads; }

        /// @brief Set the thread budget used by simulate().
        /// @param n Number of threads; 0, or more than hardware_concurrency(), restores the default
        ///          heuristic.
        void setThreadCount(size_t n)
        {
            const size_t hw = std::thread::hardware_concurrency();
            m_nThreads = (n == 0 || n > hw) ? _computeThreadCount(m_pixels.size())
                                             : std::max(size_t(1), n);
        }
#pragma endregion

#pragma region - read-only pixel access
        /// @brief Pixel of flat index index (row-major), read only. No bounds check.
        const PixelT& at(size_t index)           const { return view()[index]; }
        /// @brief Pixel at column col and row row, read only. No bounds check.
        const PixelT& at(size_t col, size_t row) const { return view()[col + row * m_width]; }

        /// @brief Pixels of one row, read only.
        std::span<const PixelT> Row(size_t row) const
        {
            return { m_pixels.data() + row * m_width, m_width };
        }

        /// @brief Pixels of one column, read only.
        StridedView<const PixelT, 1> Col(size_t col) const
        {
            return makeStridedView<const PixelT, 1>(m_pixels.data() + col, { m_height }, { m_width });
        }

        /// @brief Pixels of a window, read only, indexed [row, col] from its corner (c0, r0).
        StridedView<const PixelT, 2> window(size_t c0, size_t r0, size_t ncols, size_t nrows) const
        {
            return makeStridedView<const PixelT, 2>(m_pixels.data() + r0 * m_width + c0,
                                                    { nrows, ncols }, { m_width, 1 });
        }
#pragma endregion

#pragma region - initialisation (public)
        /// @brief Initialise every pixel from a single prototype.
        void fill(const PixelT& proto)
        {
            std::fill(m_pixels.begin(), m_pixels.end(), proto);
        }

        /// @brief Set the geometry and the pixel parameters from a configuration.
        /// @details Known keys: "gain", "dark", "readNoise", "QE", "KTCnoise", "baseline", each read
        ///          by loadMap():
        ///          - a number is applied uniformly to all pixels;
        ///          - a FITS image is a per-pixel map, "<path>" or {"file": "<path>", "hdu": h}, h an
        ///            HDU index or an EXTNAME (default: 1, the primary HDU); loadImage() lists every
        ///            accepted form. A map is either matrix-sized or detector-sized (see Geometry.h).
        ///
        ///          The geometry keys ("detectorSize", "detectorOrigin", and those of a derived
        ///          Geometry) are read first, so the maps are cropped at the new origin. The inputs
        ///          are recorded as header cards (see getInputCards()).
        ///          Derived classes override initialize() to handle additional keys; call
        ///          DetectorSimulator<PixelT>::initialize(cfg) to apply the standard keys first.
        /// @param cfg Configuration; absent keys leave the current value unchanged.
        /// @throws std::invalid_argument if a value or a map is not valid, or if the matrix does not
        ///         fit on the detector.
        virtual void initialize(const PixelConfig& cfg)
        {
            // Configure a copy and commit it only if the matrix fits on the detector.
            {
                std::shared_ptr<Geometry> geometry = m_geometry->clone();
                geometry->configure(cfg);
                geometry->place(m_width, m_height);   // throws if the matrix does not fit
                m_geometry = std::move(geometry);
            }

            // Each setter is guarded by if constexpr so that PixelT types that
            // hide pixsim's interface (e.g. via private inheritance) compile cleanly:
            // the call is not emitted if the method is not accessible on PixelT.
            if constexpr (requires(PixelT& px, double v) { px.setGain(v); })
                _applyMap(cfg, "gain",      [](PixelT& px, double v) { px.setGain(v); });
            if constexpr (requires(PixelT& px, double v) { px.setDarkCurrent(v); })
                _applyMap(cfg, "dark",      [](PixelT& px, double v) { px.setDarkCurrent(v); });
            if constexpr (requires(PixelT& px, double v) { px.setReadNoise(v); })
                _applyMap(cfg, "readNoise", [](PixelT& px, double v) { px.setReadNoise(v); });
            if constexpr (requires(PixelT& px, double v) { px.setQuantumEfficiency(v); })
                _applyMap(cfg, "QE",        [](PixelT& px, double v) { px.setQuantumEfficiency(v); });
            if constexpr (requires(PixelT& px, double v) { px.setKTCnoise(v); })
                _applyMap(cfg, "KTCnoise",  [](PixelT& px, double v) { px.setKTCnoise(v); });
            if constexpr (requires(PixelT& px, double v) { px.setBaseline(v); })
                _applyMap(cfg, "baseline",  [](PixelT& px, double v) { px.setBaseline(v); });

            // The same inputs as header cards (see getInputCards()).
            if constexpr (requires(const PixelT& px) { px.getGain(); })
                _recordInput(cfg, "gain",      "PIXGAIN", "[e-/ADU] Pixel gain",        [](const PixelT& px) { return px.getGain(); });
            if constexpr (requires(const PixelT& px) { px.getDarkCurrent(); })
                _recordInput(cfg, "dark",      "PIXDARK", "[e-/s] Pixel dark current",  [](const PixelT& px) { return px.getDarkCurrent(); });
            if constexpr (requires(const PixelT& px) { px.getReadNoise(); })
                _recordInput(cfg, "readNoise", "PIXRN",   "[e-] Pixel read noise",      [](const PixelT& px) { return px.getReadNoise(); });
            if constexpr (requires(const PixelT& px) { px.getQuantumEfficiency(); })
                _recordInput(cfg, "QE",        "PIXQE",   "Pixel quantum efficiency",   [](const PixelT& px) { return px.getQuantumEfficiency(); });
            if constexpr (requires(const PixelT& px) { px.getKTCnoise(); })
                _recordInput(cfg, "KTCnoise",  "PIXKTC",  "[e-] Pixel kTC (reset) noise", [](const PixelT& px) { return px.getKTCnoise(); });
            if constexpr (requires(const PixelT& px) { px.getBaseline(); })
                _recordInput(cfg, "baseline",  "PIXBASE", "[ADU] Pixel baseline",       [](const PixelT& px) { return px.getBaseline(); });
        }

        /// @brief initialize() from a JSON file path.
        virtual void initialize(const std::string& jsonPath) { initialize(PixelConfig(jsonPath)); }
#pragma endregion

#pragma region - per-pixel inputs from a configuration
        /// @brief The FITS image a configuration key names.
        /// @details h is an HDU index (1 = primary) or an EXTNAME:
        /// @verbatim
        ///   "key": "<path>"                         → primary HDU
        ///   "key": "<path>", "key.hdu": h           → HDU h
        ///   "key": {"file": "<path>"}               → primary HDU
        ///   "key": {"file": "<path>", "hdu": h}     → HDU h
        ///   "key.file": "<path>" [, "key.hdu": h]   → flat notation, as above
        /// @endverbatim
        /// A string starting with '{' is read as a JSON-encoded object.
        /// @param cfg Configuration.
        /// @param key Key.
        /// @return The image; nullptr when the key is absent or holds a number.
        /// @throws std::invalid_argument for a malformed entry or an HDU without an image.
        std::shared_ptr<DSL::FITScube> loadImage(const PixelConfig& cfg, const std::string& key) const
        {
            const std::optional<json> source = _fitsSource(cfg, key);
            return source ? _openImage(key, *source) : nullptr;
        }

        /// @brief Per-pixel values of a configuration key.
        /// @details A number is applied to every pixel, a FITS image (see loadImage()) is
        ///          matrix-sized or detector-sized (see toMatrix()).
        /// @param cfg Configuration.
        /// @param key Key.
        /// @return One value per pixel, row-major; empty when the key is absent.
        /// @throws std::invalid_argument for a malformed entry or a map of another size.
        std::valarray<double> loadMap(const PixelConfig& cfg, const std::string& key) const
        {
            if (const std::optional<json> source = _fitsSource(cfg, key))
                return toMatrix(*_openImage(key, *source),
                                "FITS map for '" + key + "' ('" + source->at("file").get<std::string>() + "')");
            if (!cfg.has(key))
                return {};
            return std::valarray<double>(cfg.getAsDouble(key), m_pixels.size());
        }
#pragma endregion

#pragma region - header cards
        /// @brief Simulation inputs recorded by initialize(), as FITS header cards.
        /// @details One per per-pixel parameter (PIXGAIN, PIXDARK, PIXRN, PIXQE, PIXKTC, PIXBASE),
        ///          holding the value it was given or the file name of its map ("dark.fits[1]"),
        ///          plus the cards a derived simulator adds with setInputCard(). Write them into an
        ///          output with writeInputCards().
        const std::vector<HeaderCard>& getInputCards() const { return m_inputCards; }

        /// @brief Write the input cards into the header of image.
        void writeInputCards(DSL::FITScube& image) const { writeCards(image, m_inputCards); }
#pragma endregion

#pragma region - named map setters (protected — use initialize() from public API)
    protected:
        /// @brief Set the gain of every pixel (e-/ADU). @throws std::invalid_argument on a wrong size or value.
        void setGainMap(std::span<const double> gains)
        {
            _checkMapSize(gains.size(), "Gain map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setGain(gains[i]);
        }

        /// @brief Set the dark current of every pixel (e-/s). @throws std::invalid_argument on a wrong size or value.
        void setDarkCurrentMap(std::span<const double> darks)
        {
            _checkMapSize(darks.size(), "Dark current map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setDarkCurrent(darks[i]);
        }

        /// @brief Set the read noise of every pixel (e- RMS). @throws std::invalid_argument on a wrong size or value.
        void setReadNoiseMap(std::span<const double> readNoises)
        {
            _checkMapSize(readNoises.size(), "Read noise map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setReadNoise(readNoises[i]);
        }

        /// @brief Set the quantum efficiency of every pixel. @throws std::invalid_argument on a wrong size or value.
        void setQEMap(std::span<const double> qes)
        {
            _checkMapSize(qes.size(), "QE map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setQuantumEfficiency(qes[i]);
        }

        /// @brief Set the baseline of every pixel (ADU). @throws std::invalid_argument on a wrong size or value.
        void setBaselineMap(std::span<const double> baselines)
        {
            _checkMapSize(baselines.size(), "Baseline map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setBaseline(baselines[i]);
        }

        /// @brief Set the kTC noise of every pixel (e- RMS). @throws std::invalid_argument on a wrong size or value.
        void setKTCnoiseMap(std::span<const double> ktcs)
        {
            _checkMapSize(ktcs.size(), "KTC noise map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setKTCnoise(ktcs[i]);
        }

    private:
        /// @brief The FITS file and HDU a configuration key names, as {"file": path[, "hdu": h]}.
        /// @return nullopt when the key is absent or holds a number. See loadImage() for the forms.
        static std::optional<json> _fitsSource(const PixelConfig& cfg, const std::string& key)
        {
            // A path, with the HDU, if given, in the "key.hdu" sibling.
            const auto fromPath = [&cfg, &key](const std::string& path)
            {
                json source{{"file", path}};
                if (cfg.has(key + ".hdu"))
                {
                    try                                  { source["hdu"] = static_cast<int>(cfg.getAsDouble(key + ".hdu")); }
                    catch (const std::invalid_argument&) { source["hdu"] = cfg.get<std::string>(key + ".hdu"); }
                }
                return source;
            };

            if (!cfg.has(key))
            {
                if (!cfg.has(key + ".file"))
                    return std::nullopt;
                return std::make_optional(fromPath(cfg.get<std::string>(key + ".file")));
            }

            // The JSON value is checked first: getAsDouble() would not reject an object, since
            // a json converts to double (and throws json::type_error doing so).
            json source;
            try
            {
                source = cfg.getJson(key);                        // {"file": ..., "hdu": ...}
            }
            catch (const std::invalid_argument&)                  // a path, or a number
            {
                std::string s;
                try                                    { s = cfg.get<std::string>(key); }
                catch (const std::bad_variant_access&) { return std::nullopt; }
                if (s.empty() || s.front() != '{')
                    return std::make_optional(fromPath(s));
                source = json::parse(s);                          // JSON-encoded object
            }

            if (source.is_number())
                return std::nullopt;
            if (!source.is_object() || !source.contains("file") || !source["file"].is_string())
                throw std::invalid_argument("'" + key + "': a FITS image is given as {\"file\": \"<path>\", "
                                            "\"hdu\": <index or EXTNAME>}, not " + source.dump() + ".");
            return std::make_optional(std::move(source));
        }

        /// @brief The image of a {"file": path[, "hdu": h]} source; the primary HDU (1) when h is absent.
        /// @details The image keeps its data once the file is closed.
        static std::shared_ptr<DSL::FITScube> _openImage(const std::string& key, const json& source)
        {
            const std::string path = source.at("file").get<std::string>();
            const json        hdu  = source.contains("hdu") ? source.at("hdu") : json(1);

            DSL::FITSmanager mgr(path, true);   // read only
            std::shared_ptr<DSL::FITScube> image;
            if (hdu.is_string())
                image = mgr.GetImage(hdu.get<std::string>());
            else if (hdu.is_number_integer() && hdu.get<long long>() >= 1)
                image = mgr.GetImageAtIndex(hdu.get<int>());
            else
                throw std::invalid_argument("'" + key + "': \"hdu\" is an HDU index (1 = primary) or an EXTNAME, not " +
                                            hdu.dump() + ".");

            if (!image)
                throw std::invalid_argument("'" + key + "': no image in HDU " + hdu.dump() + " of '" + path + "'.");
            return image;
        }

        /// @brief Apply setter to every pixel with the value of a key (see loadMap()); nothing if it is absent.
        template<typename Setter>
        void _applyMap(const PixelConfig& cfg,
                       const std::string& key,
                       Setter              setter)
        {
            const std::valarray<double> map = loadMap(cfg, key);
            if (map.size() == 0) return;
            for (size_t i = 0; i < m_pixels.size(); ++i)
                setter(m_pixels[i], map[i]);
        }

    public:
#pragma endregion

#pragma region - simulation
        // simulate() is NOT virtual. It owns thread dispatch and calls _engine().
        // To customise physics, override _engine() in a derived class.

        /// @brief Full matrix, uniform flux.
        /// @param flux         Photons/s/pixel.
        /// @param exposureTime Exposure time (s).
        /// @return Image of getWidth() × getHeight() pixels.
        std::shared_ptr<DSL::FITScube> simulate(double flux, double exposureTime) const
        {
            return simulate(std::valarray<double>(flux, m_pixels.size()), exposureTime);
        }

        /// @brief Full matrix, flux of a double-precision FITS image.
        /// @param flux         Photons/s/pixel, matrix-sized.
        /// @param exposureTime Exposure time (s).
        /// @return Image of getWidth() × getHeight() pixels.
        std::shared_ptr<DSL::FITScube> simulate(const std::shared_ptr<DSL::FITScube>& flux, double exposureTime) const
        {
            std::valarray<double> fluxArray(0.0, m_pixels.size());
            fluxArray += *flux->GetData<double>();
            return simulate(fluxArray, exposureTime);
        }

        /// @brief Full matrix, per-pixel flux.
        /// @details Dispatches to _engine() across std::thread workers when the detector is large.
        /// @param flux         Photons/s of every pixel, row-major.
        /// @param exposureTime Exposure time (s).
        /// @return Image of getWidth() × getHeight() pixels.
        std::shared_ptr<DSL::FITScube> simulate(const std::valarray<double>& flux, double exposureTime) const
        {
            const size_t n        = m_pixels.size();
            const size_t nThreads = m_nThreads;

            std::valarray<double> frame(0.0, n);

            if (nThreads > 1)
                _dispatch_cpu(flux, exposureTime, frame, nThreads);
            else
                _engine(flux, exposureTime, 0, n, frame);

            return _convertFrame(frame, m_width, m_height);
        }

        /// @brief Window: simulate only a sub-region of the matrix.
        /// @details Use when a variable/moving source is localised in a small region: run the full
        ///          detector at coarse cadence, then re-run the affected window at finer time steps
        ///          and composite the results.
        /// @param flux         Photons/s of the window pixels (ncols × nrows, row-major).
        /// @param exposureTime Exposure time (s).
        /// @param c0           Matrix column of the first window pixel.
        /// @param r0           Matrix row of the first window pixel.
        /// @param ncols        Window columns.
        /// @param nrows        Window rows.
        /// @return Image of ncols × nrows pixels.
        /// @throws std::out_of_range if the window exceeds the matrix, std::invalid_argument if flux
        ///         does not have the size of the window.
        std::shared_ptr<DSL::FITScube> simulate(const std::valarray<double>& flux,   // ncols*nrows elements, row-major
                  const double&                        exposureTime,
                  const size_t&                        c0,
                  const size_t&                        r0,
                  const size_t&                        ncols,
                  const size_t&                        nrows) const
        {

            if (r0 + nrows > m_height || c0 + ncols > m_width)
                throw std::out_of_range("Simulation window exceeds detector boundaries.");
            if (flux.size() != nrows * ncols)
                throw std::invalid_argument("Flux size does not match window size (nrows*ncols).");

            const size_t nThreads = std::min(m_nThreads, _computeThreadCount(nrows * ncols));
            std::valarray<double> frame(0.0, nrows * ncols);

            if (nThreads > 1)
                _dispatch_cpu_window(flux, exposureTime, c0, r0, ncols, nrows, frame, nThreads);
            else
                _engine(flux, exposureTime, c0, r0, 0, nrows, ncols, frame);

            return _convertFrame(frame, ncols, nrows);
        }

        /// @brief Full matrix, with runtime parameters.
        /// @details cfg (e.g. persistence dT) is passed down to pixelresponse() via the cfg-aware
        ///          _engine() hooks.
        /// @param flux         Photons/s of every pixel, row-major.
        /// @param exposureTime Exposure time (s).
        /// @param cfg          Runtime parameters.
        /// @return Image of getWidth() × getHeight() pixels.
        std::shared_ptr<DSL::FITScube> simulate(const std::valarray<double>& flux,
                                                 double                        exposureTime,
                                                 const PixelConfig&            cfg) const
        {
            const size_t n        = m_pixels.size();
            const size_t nThreads = m_nThreads;
            std::valarray<double> frame(0.0, n);

            if (nThreads > 1)
                _dispatch_cpu_cfg(flux, exposureTime, frame, nThreads, cfg);
            else
                _engine(flux, exposureTime, 0, n, frame, cfg);

            return _convertFrame(frame, m_width, m_height);
        }

        /// @brief Window, with runtime parameters.
        /// @param flux         Photons/s of the window pixels (ncols × nrows, row-major).
        /// @param exposureTime Exposure time (s).
        /// @param c0           Matrix column of the first window pixel.
        /// @param r0           Matrix row of the first window pixel.
        /// @param ncols        Window columns.
        /// @param nrows        Window rows.
        /// @param cfg          Runtime parameters.
        /// @return Image of ncols × nrows pixels.
        /// @throws std::out_of_range if the window exceeds the matrix, std::invalid_argument if flux
        ///         does not have the size of the window.
        std::shared_ptr<DSL::FITScube> simulate(const std::valarray<double>& flux,
                                                 double                        exposureTime,
                                                 const size_t&                 c0,
                                                 const size_t&                 r0,
                                                 const size_t&                 ncols,
                                                 const size_t&                 nrows,
                                                 const PixelConfig&            cfg) const
        {
            if (r0 + nrows > m_height || c0 + ncols > m_width)
                throw std::out_of_range("Simulation window exceeds detector boundaries.");
            if (flux.size() != nrows * ncols)
                throw std::invalid_argument("Flux size does not match window size (nrows*ncols).");

            const size_t nThreads = std::min(m_nThreads, _computeThreadCount(nrows * ncols));
            std::valarray<double> frame(0.0, nrows * ncols);

            if (nThreads > 1)
                _dispatch_cpu_window_cfg(flux, exposureTime, c0, r0, ncols, nrows, frame, nThreads, cfg);
            else
                _engine(flux, exposureTime, c0, r0, 0, nrows, ncols, frame, cfg);

            return _convertFrame(frame, ncols, nrows);
        }
#pragma endregion

#pragma endregion
    };
#pragma endregion
} // namespace phot2pix
