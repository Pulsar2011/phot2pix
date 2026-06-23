//
// Created by gillard on 2026-06-05.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// DetectorSimulator.h
//
#pragma once

#include "PixelSimulator.h"

#include <vector>
#include <span>
#include <ranges>
#include <mdspan>
#include <array>
#include <stdexcept>
#include <string>
#include <thread>       // std::thread, hardware_concurrency()
#include <valarray>
#include <DSTfits/FITSmanager.h>  // DSL::FITSmanager
#include <DSTfits/FITSimg.h>  // DSL::FITSmanager

namespace phot2pix
{
#pragma region - DetectorSimulator class definition
    template<typename PixelT = pixsim>
    class DetectorSimulator
    {
        /*
        @class DetectorSimulator
        @brief Default detector simulator class to simulate the detector response to photons.
        @details Flat row-major pixel array with 2D access.
                 Non-const pixel accessors are protected so only derived classes can
                 modify individual pixels. Public API exposes const (read) access and
                 batch configuration methods (gain map, dark map, ...) that operate on
                 the whole detector coherently.

                 Parallelism is handled internally by simulate() via std::thread.
                 Derived classes customise physics by overriding _engine(), not simulate().
        */

        static_assert(std::is_base_of_v<pixsim, PixelT>, "PixelT must derive from pixsim");

#pragma region – private
// ── Private: thread dispatch ───────────────────────────────────────────────
    private:

        auto view()       { return std::span<PixelT>      (m_pixels.data(), m_pixels.size()); }
        auto view() const { return std::span<const PixelT>(m_pixels.data(), m_pixels.size()); }

        void _checkMapSize(size_t sz, const char* name) const
        {
            if (sz != m_pixels.size())
                throw std::invalid_argument(std::string(name) + " size does not match detector size.");
        }

        // Returns the number of threads to use for n pixels.
        // Below 100 * nCores the thread-spawn overhead exceeds the gain.
        size_t _computeThreadCount(size_t n) const
        {
            const size_t maxT = std::thread::hardware_concurrency();
            return (n > 100 * maxT) ? std::min(maxT, n) : 1;
        }

        // Full-detector dispatch: partition flat pixel index range across threads.
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

        // Windowed dispatch: partition window rows across threads.
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
#pragma endregion

// ── Protected: data, pixel write access, and virtual physics hook ──────────
#pragma region - protected
    protected:
#pragma region -> members
        const size_t        m_width;
        const size_t        m_height;
        std::vector<PixelT> m_pixels; ///< Flat row-major pixel array
        const int m_pixelType;        ///< FITS storage type for output cube
        
#pragma endregion
#pragma region -> member functions

        // Non-const pixel accessors — derived classes only.
        // Use these during initialisation (defect maps, special pixel regions, etc.).
        PixelT& at(size_t index)           { return view()[index]; }
        PixelT& at(size_t col, size_t row) { return view()[col + row * m_width]; }

        std::span<PixelT> Row(size_t row)
        {
            return { m_pixels.data() + row * m_width, m_width };
        }

        auto Col(size_t col)
        {
            using ext_t = std::dextents<size_t, 1>;
            std::array<size_t, 1> strides = { m_width };
            return std::mdspan<PixelT, ext_t, std::layout_stride>(
                m_pixels.data() + col,
                std::layout_stride::mapping<ext_t>{ ext_t{m_height}, strides }
            );
        }

        auto window(size_t c0, size_t r0, size_t ncols, size_t nrows)
        {
            using ext_t = std::dextents<size_t, 2>;
            std::array<size_t, 2> strides = { m_width, 1 };
            return std::mdspan<PixelT, ext_t, std::layout_stride>(
                m_pixels.data() + r0 * m_width + c0,
                std::layout_stride::mapping<ext_t>{ ext_t{nrows, ncols}, strides }
            );
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

        // Full-detector hook: thread owns pixels [begin, end).
        virtual void _engine(const std::valarray<double>& flux,
                              double                        dt,
                              size_t                        begin,
                              size_t                        end,
                              std::valarray<double>&         frame) const
        {
            for (size_t i = begin; i < end; ++i)
                frame[i] = m_pixels[i].pixelresponse(flux[i], dt);
        }

        // Windowed hook: thread owns window rows [row_begin, row_end).
        // flux and frame are window-sized (ncols * nrows elements, row-major).
        // c0, r0 are the detector coordinates of the window top-left corner.
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

#pragma endregion
#pragma region -> array format convertion
        // Format convertion functions for outputting the detector response to FITS files.
        // They are marked as virtual so that derived classes can override them to implement custom clipping or scaling.

        //@brief Convert a valarray of doubles to a valarray of int8_t, with clipping to the range of int8_t.
        //@param arr The input valarray of doubles to convert.
        //@return A valarray of int8_t containing the converted values.
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

        //@brief Convert a valarray of doubles to a valarray of uint8_t, with clipping to the range of uint8_t.
        //@param arr The input valarray of doubles to convert.
        //@return A valarray of uint8_t containing the converted values.
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

        //@brief Convert a valarray of doubles to a valarray of int16_t, with clipping to the range of int16_t.
        //@param arr The input valarray of doubles to convert.
        //@return A valarray of int16_t containing the converted values.
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

        //@brief Convert a valarray of doubles to a valarray of uint16_t, with clipping to the range of uint16_t.
        //@param arr The input valarray of doubles to convert.
        //@return A valarray of uint16_t containing the converted values.
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

        //@brief Convert a valarray of doubles to a valarray of int32_t, with clipping to the range of int32_t.
        //@param arr The input valarray of doubles to convert.
        //@return A valarray of int32_t containing the converted values.
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

        //@brief Convert a valarray of doubles to a valarray of uint32_t, with clipping to the range of uint32_t.
        //@param arr The input valarray of doubles to convert.
        //@return A valarray of uint32_t containing the converted values.
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

        //@brief Convert a valarray of doubles to a valarray of int64_t, with clipping to the range of int64_t.
        //@param arr The input valarray of doubles to convert.
        //@return A valarray of int64_t containing the converted values.
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

        //@brief Convert a valarray of doubles to a valarray of uint64_t, with clipping to the range of uint64_t.
        //@param arr The input valarray of doubles to convert.
        //@return A valarray of uint64_t containing the converted values.
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

        //@brief Convert a valarray of doubles to a valarray of float, with clipping to the range of float.
        //@param arr The input valarray of doubles to convert.
        //@return A valarray of float containing the converted values.
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
        DetectorSimulator(size_t width, size_t height, const int& dataType = DOUBLE_IMG)
            : m_width(width), m_height(height), m_pixels(width * height), m_pixelType(dataType)
        {
            if (width == 0 || height == 0)
                throw std::invalid_argument("Width and height must be positive integers.");
        }

        virtual ~DetectorSimulator() = default;
#pragma endregion

#pragma region - geometry
        size_t getWidth()  const { return m_width; }
        size_t getHeight() const { return m_height; }
        size_t size()      const { return m_pixels.size(); }
#pragma endregion

#pragma region - read-only pixel access
        const PixelT& at(size_t index)           const { return view()[index]; }
        const PixelT& at(size_t col, size_t row) const { return view()[col + row * m_width]; }

        std::span<const PixelT> Row(size_t row) const
        {
            return { m_pixels.data() + row * m_width, m_width };
        }

        auto Col(size_t col) const
        {
            using ext_t = std::dextents<size_t, 1>;
            std::array<size_t, 1> strides = { m_width };
            return std::mdspan<const PixelT, ext_t, std::layout_stride>(m_pixels.data() + col,
                                                                        std::layout_stride::mapping<ext_t>{ ext_t{m_height}, strides });
        }

        auto window(size_t c0, size_t r0, size_t ncols, size_t nrows) const
        {
            using ext_t = std::dextents<size_t, 2>;
            std::array<size_t, 2> strides = { m_width, 1 };
            return std::mdspan<const PixelT, ext_t, std::layout_stride>(m_pixels.data() + r0 * m_width + c0,
                                                                        std::layout_stride::mapping<ext_t>{ ext_t{nrows, ncols}, strides });
        }
#pragma endregion

#pragma region - initialisation (public)
        // Initialise every pixel from a single prototype.
        void fill(const PixelT& proto)
        {
            std::fill(m_pixels.begin(), m_pixels.end(), proto);
        }

        // Primary pixel-property initialisation from a typed config.
        // Each known key is interpreted as:
        //   numeric value  → applied uniformly to all pixels
        //   string value   → path to a FITS file containing a per-pixel map;
        //                    an optional "<key>.hdu" sibling specifies the HDU
        //                    as an integer index or an EXTNAME string (default: 0).
        //
        // Known keys: "gain", "dark", "readNoise", "QE", "KTCnoise", "baseline"
        // Unknown keys are forwarded to _applyCustomConfig() for derived classes.
        void initialize(const PixelConfig& cfg)
        {
            _applyMap(cfg, "gain",      [](PixelT& px, double v) { px.setGain(v); });
            _applyMap(cfg, "dark",      [](PixelT& px, double v) { px.setDarkCurrent(v); });
            _applyMap(cfg, "readNoise", [](PixelT& px, double v) { px.setReadNoise(v); });
            _applyMap(cfg, "QE",        [](PixelT& px, double v) { px.setQuantumEfficiency(v); });
            _applyMap(cfg, "KTCnoise",  [](PixelT& px, double v) { px.setKTCnoise(v); });
            _applyMap(cfg, "baseline",  [](PixelT& px, double v) { px.setBaseline(v); });
            _applyCustomConfig(cfg);
        }

        // Convenience: load config from a JSON file path.
        void initialize(const std::string& jsonPath) { initialize(PixelConfig(jsonPath)); }
#pragma endregion

#pragma region - named map setters (protected — use initialize() from public API)
    protected:
        void setGainMap(std::span<const double> gains)
        {
            _checkMapSize(gains.size(), "Gain map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setGain(gains[i]);
        }

        void setDarkCurrentMap(std::span<const double> darks)
        {
            _checkMapSize(darks.size(), "Dark current map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setDarkCurrent(darks[i]);
        }

        void setReadNoiseMap(std::span<const double> readNoises)
        {
            _checkMapSize(readNoises.size(), "Read noise map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setReadNoise(readNoises[i]);
        }

        void setQEMap(std::span<const double> qes)
        {
            _checkMapSize(qes.size(), "QE map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setQuantumEfficiency(qes[i]);
        }

        void setBaselineMap(std::span<const double> baselines)
        {
            _checkMapSize(baselines.size(), "Baseline map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setBaseline(baselines[i]);
        }

        void setKTCnoiseMap(std::span<const double> ktcs)
        {
            _checkMapSize(ktcs.size(), "KTC noise map");
            for (size_t i = 0; i < m_pixels.size(); ++i)
                m_pixels[i].setKTCnoise(ktcs[i]);
        }

        // Hook for derived DetectorSimulator classes: handle config keys that
        // correspond to properties of a derived pixsim subclass.
        // Called at the end of initialize() after all known keys are processed.
        virtual void _applyCustomConfig(const PixelConfig&) {}

    private:
        // Resolve a config key to a per-pixel valarray:
        //   - empty valarray if key is absent (caller skips the property)
        //   - uniform valarray if the value is numeric (via PixelConfig::getAsDouble)
        //   - per-pixel valarray if the value is a FITS path string
        std::valarray<double> _resolveMap(const PixelConfig& cfg,
                                           const std::string& key) const
        {
            if (!cfg.has(key)) return {};

            // Scalar: PixelConfig owns the any-numeric-type → double conversion
            try
            {
                return std::valarray<double>(cfg.getAsDouble(key), m_pixels.size());
            }
            catch (const std::invalid_argument&) {}   // not numeric → try FITS path

            // FITS map
            const std::string path = cfg.get<std::string>(key);

            DSL::FITSmanager mgr(path);
            std::shared_ptr<DSL::FITScube> cube;

            if (cfg.has(key + ".hdu"))
            {
                // HDU spec: integer index or EXTNAME string
                try
                {
                    cube = mgr.GetImageAtIndex(
                        static_cast<int>(cfg.getAsDouble(key + ".hdu")));
                }
                catch (const std::invalid_argument&)
                {
                    cube = mgr.GetImage(cfg.get<std::string>(key + ".hdu"));
                }
            }
            else
            {
                cube = mgr.GetImageAtIndex(0);
            }

            if (cube->Nelements() != m_pixels.size())
                throw std::invalid_argument(
                    "FITS map for '" + key + "': " +
                    std::to_string(cube->Nelements()) +
                    " pixels, detector has " +
                    std::to_string(m_pixels.size()) + ".");

            // cast<double>() converts any FITS storage type to double
            auto doubleCube = cube->cast<double>();
            return *doubleCube->GetData<double>();
        }

        template<typename Setter>
        void _applyMap(const PixelConfig& cfg,
                       const std::string& key,
                       Setter              setter)
        {
            auto map = _resolveMap(cfg, key);
            if (map.size() == 0) return;
            for (size_t i = 0; i < m_pixels.size(); ++i)
                setter(m_pixels[i], map[i]);
        }

    public:
#pragma endregion

#pragma region - simulation
        // simulate() is NOT virtual. It owns thread dispatch and calls _engine().
        // To customise physics, override _engine() in a derived class.

        // Convenience: uniform flux over the full detector.
        std::shared_ptr<DSL::FITScube> simulate(double flux, double exposureTime) const
        {
            return simulate(std::valarray<double>(flux, m_pixels.size()), exposureTime);
        }

        // Convenience: per-pixel flux from a FITScube.
        std::shared_ptr<DSL::FITScube> simulate(const std::shared_ptr<DSL::FITScube>& flux, double exposureTime) const
        {
            std::valarray<double> fluxArray(0.0, m_pixels.size());
            fluxArray += *flux->GetData<double>();
            return simulate(fluxArray, exposureTime);
        }

        // Core: per-pixel flux, full detector.
        // Dispatches to _engine() across std::thread workers when detector is large.
        std::shared_ptr<DSL::FITScube> simulate(const std::valarray<double>& flux, double exposureTime) const
        {
            const size_t n        = m_pixels.size();
            const size_t nThreads = _computeThreadCount(n);

            std::valarray<double> frame(0.0, n);

            if (nThreads > 1)
                _dispatch_cpu(flux, exposureTime, frame, nThreads);
            else
                _engine(flux, exposureTime, 0, n, frame);

            std::shared_ptr< DSL::FITScube > img;

            switch (m_pixelType)
            {
                case SBYTE_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<int8_t>({m_width, m_height}));
                    *(img) += convert2int8(frame);
                    break;
                
                case BYTE_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<uint8_t>({m_width, m_height}));
                    *(img) += convert2uint8(frame);
                    break;

                case SHORT_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<int16_t>({m_width, m_height}));
                    *(img) += convert2int16(frame);
                    break;

                case USHORT_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<uint16_t>({m_width, m_height}));
                    *(img) += convert2uint16(frame);
                    break;

                case LONG_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<int32_t>({m_width, m_height}));
                    *(img) += convert2int32(frame);
                    break;
                
                case ULONG_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<uint32_t>({m_width, m_height}));
                    *(img) += convert2uint32(frame);
                    break;

                case LONGLONG_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<int64_t>({m_width, m_height}));
                    *(img) += convert2int64(frame);
                    break;

                case ULONGLONG_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<uint64_t>({m_width, m_height}));
                    *(img) += convert2uint64(frame);
                    break;

                case FLOAT_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<float>({m_width, m_height}));
                    *(img) += convert2float(frame);
                    break;

                case DOUBLE_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<double>({m_width, m_height}));
                    *(img) += frame;
                    break;
                
                default:
                    throw std::invalid_argument("Unrecognised FITS data type for output.");
            }

            return img;
        }

        // Windowed: simulate only a detector sub-region.
        // Returns a FITScube of size ncols × nrows.
        // Use when a variable/moving source is localised in a small region:
        // run the full detector at coarse cadence, then re-run the affected
        // window at finer time steps and composite the results.
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

            const size_t nThreads = _computeThreadCount(nrows * ncols);
            std::valarray<double> frame(0.0, nrows * ncols);

            if (nThreads > 1)
                _dispatch_cpu_window(flux, exposureTime, c0, r0, ncols, nrows, frame, nThreads);
            else
                _engine(flux, exposureTime, c0, r0, 0, nrows, ncols, frame);

            std::shared_ptr< DSL::FITScube > img;

            switch (m_pixelType)
            {
                case SBYTE_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<int8_t>({m_width, m_height}));
                    *(img) += convert2int8(frame);
                    break;
                
                case BYTE_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<uint8_t>({m_width, m_height}));
                    *(img) += convert2uint8(frame);
                    break;

                case SHORT_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<int16_t>({m_width, m_height}));
                    *(img) += convert2int16(frame);
                    break;

                case USHORT_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<uint16_t>({m_width, m_height}));
                    *(img) += convert2uint16(frame);
                    break;

                case LONG_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<int32_t>({m_width, m_height}));
                    *(img) += convert2int32(frame);
                    break;
                
                case ULONG_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<uint32_t>({m_width, m_height}));
                    *(img) += convert2uint32(frame);
                    break;

                case LONGLONG_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<int64_t>({m_width, m_height}));
                    *(img) += convert2int64(frame);
                    break;

                case ULONGLONG_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<uint64_t>({m_width, m_height}));
                    *(img) += convert2uint64(frame);
                    break;

                case FLOAT_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<float>({m_width, m_height}));
                    *(img) += convert2float(frame);
                    break;

                case DOUBLE_IMG:
                    img = std::shared_ptr< DSL::FITScube >(new DSL::FITSimg<double>({m_width, m_height}));
                    *(img) += frame;
                    break;
                
                default:
                    throw std::invalid_argument("Unrecognised FITS data type for output.");
            }

            return img;
        }
#pragma endregion

#pragma endregion
    };
#pragma endregion
} // namespace phot2pix
