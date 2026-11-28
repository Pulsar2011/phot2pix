//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispSCA.cxx
//

#include "nispSCA.h"
#include "nispPixel.h"

#include <phot2pix/DetectorTools.h>
#include <phot2pix/CosmicDeposit.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <future>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <valarray>

namespace nisp
{
// ── MACC parameters ───────────────────────────────────────────────────────────
#pragma region - macc
    namespace
    {
        // Validate that v is a non-negative integer and return it as size_t.
        // Goes through double so it accepts whatever integer type PixelConfig
        // stored (JSON gives uint64_t, which is not size_t on macOS).
        size_t toCount(double v, const char* name)
        {
            if (!(v >= 0.0) || v != std::floor(v))
                throw std::invalid_argument(std::string("macc::fromConfig: '") + name + "' must be a non-negative integer.");
            return static_cast<size_t>(v);
        }
    }

    bool macc::inConfig(const phot2pix::PixelConfig& cfg)
    {
        return cfg.has("macc") || cfg.has("ng") || cfg.has("nr") || cfg.has("nd");
    }

    macc macc::fromConfig(const phot2pix::PixelConfig& cfg)
    {
        const bool hasArray = cfg.has("macc");
        const bool hasKeys  = cfg.has("ng") || cfg.has("nr") || cfg.has("nd");

        if (hasArray && hasKeys)
            throw std::invalid_argument("macc::fromConfig: give either 'macc'=[ng,nr,nd] or 'ng','nr','nd', not both.");
        if (!hasArray && !hasKeys)
            throw std::invalid_argument("macc::fromConfig: no 'macc'=[ng,nr,nd] or 'ng','nr','nd' in configuration.");

        const char* names[3] = {"ng", "nr", "nd"};
        size_t      v[3]     = {0, 0, 0};

        if (hasArray)
        {
            const phot2pix::json& j = cfg.getJson("macc");
            if (!j.is_array() || j.size() != 3)
                throw std::invalid_argument("macc::fromConfig: 'macc' must be an array [ng, nr, nd].");

            for (size_t i = 0; i < 3; ++i)
            {
                if (!j[i].is_number())
                    throw std::invalid_argument(std::string("macc::fromConfig: 'macc' element '") + names[i] + "' is not a number.");
                v[i] = toCount(j[i].get<double>(), names[i]);
            }
        }
        else
        {
            for (size_t i = 0; i < 3; ++i)
            {
                if (!cfg.has(names[i]))
                    throw std::invalid_argument(std::string("macc::fromConfig: missing key '") + names[i] + "'.");
                v[i] = toCount(cfg.getAsDouble(names[i]), names[i]);
            }
        }

        if (v[0] == 0 || v[1] == 0)
            throw std::invalid_argument("macc::fromConfig: ng and nr must be >= 1.");

        return macc{v[0], v[1], v[2]};   // {nGroups, nFrames, nDrop}
    }
#pragma endregion

// ── Readout mode and exposure parameters ──────────────────────────────────────
#pragma region - readout
    mode modeFromString(const std::string& name)
    {
        std::string key = name;
        std::transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

        if      (key == "MACC") return mode::MACC;
        else if (key == "UTR")  return mode::UTR;
        else if (key == "CDS")  return mode::CDS;
        else if (key == "DPU")  return mode::DPU;
        else if (key == "FLUX") return mode::FLUX;

        throw std::invalid_argument("unknown exposureMode '" + name + "' (expected MACC, UTR, CDS, DPU or FLUX).");
    }

    readout readout::fromConfig(const phot2pix::PixelConfig& cfg)
    {
        const bool hasTime = cfg.has("exposureTime");
        const bool hasMacc = macc::inConfig(cfg);

        if (hasTime && hasMacc)
            throw std::invalid_argument("readout::fromConfig: give either 'exposureTime' or the MACC parameters "
                                        "('macc'=[ng,nr,nd] or 'ng','nr','nd'), not both.");
        if (!hasTime && !hasMacc)
            throw std::invalid_argument("readout::fromConfig: no 'exposureTime' or MACC parameters "
                                        "('macc'=[ng,nr,nd] or 'ng','nr','nd').");

        readout r;

        if (cfg.has("exposureMode"))
        {
            std::string name;
            try                                    { name = cfg.get<std::string>("exposureMode"); }
            catch (const std::bad_variant_access&) { throw std::invalid_argument("readout::fromConfig: 'exposureMode' must be "
                                                                                 "MACC, UTR, CDS, DPU or FLUX."); }
            r.exposureMode = modeFromString(name);
        }

        if (hasMacc)
            r.maccParams = macc::fromConfig(cfg);
        else
        {
            // getAsDouble() throws std::invalid_argument for a string and json::type_error for
            // a JSON array or object.
            double t = -1.0;
            try                          { t = cfg.getAsDouble("exposureTime"); }
            catch (const std::exception&) {}
            if (!(t > 0.0))
                throw std::invalid_argument("readout::fromConfig: 'exposureTime' must be a number of seconds > 0.");
            r.exposureTime = t;
        }

        return r;
    }
#pragma endregion

    namespace persistence
    {

// ── model ─────────────────────────────────────────────────────────────────────
#pragma region - model
        void model::insertParameter(const double& factor, const double& decayTime)
        {
            const size_t idx = _factor.size();
            _factor[idx]    = factor;
            _decayTime[idx] = decayTime;
        }

        double model::pixelresponse(const double& flux, const double& exposureTime) const
        {
            phot2pix::PixelConfig cfg;
            cfg.set("dT", dT);
            return pixelresponse(flux, exposureTime, cfg);
        }

        // Charge released during a step that starts timeDelay after the end of the
        // exposure that filled the traps and lasts exposureTime. Trap k releases
        // E·f_k·exp(−t/τ_k) per second, t being the time since that exposure, so over
        // the step [Δ, Δ + T]:
        //   ∫ f e^(−t/τ) dt = f·τ·e^(−Δ/τ)·(1 − e^(−T/τ))
        // Consecutive steps add up to the release over their total duration, and the
        // whole release (T → ∞ from Δ = 0) is f·τ of the stored charge.
        double model::pixelresponse(const double& flux, const double& exposureTime, const phot2pix::PixelConfig& cfg) const
        {
            const double timeDelay = cfg.has("dT") ? cfg.getAsDouble("dT") : dT;
            double signal = 0.0;
            for (const auto& [idx, f] : _factor)
            {
                const double tau = _decayTime.at(idx);
                signal += f * tau * std::exp(-timeDelay / tau) * (1.0 - std::exp(-exposureTime / tau));
            }
            return signal * flux / getGain();   // e- → ADU (default gain 1: e-)
        }
#pragma endregion
// ── pixelArray ────────────────────────────────────────────────────────────────

#pragma region - pixelArray
        void pixelArray::incrementTime(const double& dT)
        {
            for (size_t i = 0; i < size(); ++i)
                at(i).incrementTime(dT);
        }

        void pixelArray::setGain(const std::valarray<double>& gain)
        {
            setGainMap(std::span<const double>(std::begin(gain), gain.size()));   // throws on a wrong size or value
        }

        void pixelArray::initialize(const phot2pix::PixelConfig& cfg)
        {
            if (!cfg.has("persistence"))
                return;

            const auto& persistJson = cfg.getJson("persistence");

            // Two supported formats:
            //   array:  [{"f": ..., "t": ...}, ...]
            //   object: {"traps": [...], "dT": <seconds, optional>}
            // "dT" starts the clock of the pixel models, for a standalone pixelArray. Within
            // the SCA, the scheduler passes its own clock and applies "dT" as the delay
            // between two exposures (scheduler::insertCharges).
            double initialDT = 0.0;
            if (persistJson.is_object() && persistJson.contains("dT"))
                initialDT = persistJson.at("dT").get<double>();

            if(initialDT > std::numeric_limits<double>::epsilon())
                incrementTime(initialDT);

            const auto& traps = (persistJson.is_object() && persistJson.contains("traps"))
                ? persistJson.at("traps") : persistJson;

            for (size_t i = 0; i < size(); ++i)
            {
                at(i).clearTraps();
                for (const auto& entry : traps)
                    at(i).insertParameter(entry.at("f").get<double>(), entry.at("t").get<double>());
            }
        }

        void pixelArray::_engine(const std::valarray<double>& flux,
                     double dt, size_t begin, size_t end,
                     std::valarray<double>& frame,
                     const phot2pix::PixelConfig& cfg) const
        {
            for (size_t i = begin; i < end; ++i)
                frame[i] = at(i).pixelresponse(flux[i], dt, cfg);
        }

        void pixelArray::_engine(const std::valarray<double>& flux,
                     double dt, size_t c0, size_t r0,
                     size_t row_begin, size_t row_end, size_t ncols,
                     std::valarray<double>& frame,
                     const phot2pix::PixelConfig& cfg) const
        {
            for (size_t r = row_begin; r < row_end; ++r)
            {
                const size_t det_base   = (r0 + r) * getWidth() + c0;
                const size_t frame_base = r * ncols;
                for (size_t c = 0; c < ncols; ++c)
                    frame[frame_base + c] = at(det_base + c).pixelresponse(flux[frame_base + c], dt, cfg);
            }
        }

#pragma endregion

// ── scheduler ─────────────────────────────────────────────────────────────────
#pragma region - scheduler
        scheduler::scheduler(const phot2pix::PixelConfig& cfg)
            : m_pixelModel(nullptr), m_fluxHistory(), m_cfg(cfg)
        {
            if (!cfg.has("persistence"))
                return;

            const auto& p = cfg.getJson("persistence");
            if (p.is_object() && p.contains("dT"))
            {
                const auto& d = p.at("dT");
                if (!d.is_number() || !(d.get<double>() >= 0.0))
                    throw std::invalid_argument("persistence::scheduler: 'dT' must be a number >= 0 (seconds).");
                m_delay = d.get<double>();
            }
        }

        void scheduler::setGain(const std::valarray<double>& gain)
        {
            if (m_pixelModel)
                m_pixelModel->setGain(gain.size() > 0 ? gain : std::valarray<double>(1.0, m_pixelModel->size()));
            m_gain = gain;
        }

        void scheduler::insertCharges(const std::shared_ptr<DSL::FITScube>& flux, const double& t0)
        {
            if (!m_pixelModel)
            {
                // Built aside, kept only when complete: a bad gain leaves the scheduler unchanged.
                auto pixels = std::make_unique<pixelArray>(flux->Size(1), flux->Size(2));
                pixels->initialize(m_cfg);
                if (m_gain.size() > 0)
                    pixels->setGain(m_gain);   // released charge in ADU
                if (m_threadCount > 0)
                    pixels->setThreadCount(m_threadCount);
                m_pixelModel = std::move(pixels);
            }
            else if (flux->Size(1) != m_pixelModel->getWidth() ||
                     flux->Size(2) != m_pixelModel->getHeight())
            {
                throw std::invalid_argument(
                    "scheduler::insertCharges: flux dimensions do not match model.");
            }

            // The SCA simulates exposures back to back: the time between the end of this
            // exposure and the start of the next one (dithers, overheads…) is m_delay,
            // during which every stored exposure releases charge that nobody reads.
            for (auto& record : m_fluxHistory)
                record.t0 += m_delay;

            m_fluxHistory.push_back({t0 + m_delay, flux});
        }

        void scheduler::removeCharges(const double& t0)
        {
            auto it = std::find_if(m_fluxHistory.begin(), m_fluxHistory.end(),
                    [&t0](const ExposureRecord& r)
                            {
                                const double scale = std::abs(r.t0) + std::abs(t0) + 1.0;
                                return std::abs(r.t0 - t0) <= std::numeric_limits<double>::epsilon() * scale;
                            }
                        );

            if (it != m_fluxHistory.end())
                m_fluxHistory.erase(it);
        }

        std::shared_ptr<DSL::FITScube> scheduler::computePersistenceSignal(const double& dT, const double& exposureTime)
        {
            if (m_fluxHistory.empty())
                return nullptr;

            std::shared_ptr<DSL::FITScube> result;

            for (auto& record : m_fluxHistory)
            {
                phot2pix::PixelConfig cfg;
                record.t0 += dT;
                cfg.set("dT", record.t0);

                std::valarray<double> fluxArray(0.0, m_pixelModel->size());
                fluxArray += *(record.flux->GetData<double>());

                std::shared_ptr<DSL::FITScube> signal = m_pixelModel->simulate(fluxArray, exposureTime, cfg);

                record.t0 += exposureTime;

                if (!result)
                    result = std::move(signal);
                else
                    *(result->GetData<double>()) += *(signal->GetData<double>());
            }

            return result;
        }

        std::shared_ptr<DSL::FITScube> scheduler::computePersistenceSignal(
            const double& dT, const double& exposureTime,
            const size_t& c0, const size_t& r0,
            const size_t& ncols, const size_t& nrows)
        {
            if (m_fluxHistory.empty())
                return nullptr;

            std::shared_ptr<DSL::FITScube> result;

            const size_t stride = m_pixelModel->getWidth();

            for (auto& record : m_fluxHistory)
            {
                // Same clock as the full-frame overload: the gap dT is part of the elapsed time.
                phot2pix::PixelConfig cfg;
                record.t0 += dT;
                cfg.set("dT", record.t0);

                // Extract the sub-window from the stored full-frame flux.
                const auto& fullFlux = *(record.flux->GetData<double>());
                std::valarray<double> windowFlux(0.0, ncols * nrows);
                for (size_t r = 0; r < nrows; ++r)
                    windowFlux[std::slice(r * ncols, ncols, 1)] =
                        fullFlux[std::slice((r0 + r) * stride + c0, ncols, 1)];

                std::shared_ptr<DSL::FITScube> signal = m_pixelModel->simulate(windowFlux, exposureTime, c0, r0, ncols, nrows, cfg);

                record.t0 += exposureTime;

                if (!result)
                    result = std::move(signal);
                else
                    *(result->GetData<double>()) += *(signal->GetData<double>());
            }

            return result;
        }

    } // namespace persistence
#pragma endregion

// ── SCA ───────────────────────────────────────────────────────────────────────
#pragma region - SCA

    namespace
    {
        // Same result as FITScube::cast<uint16_t>(): round to nearest, clamp to [0, 65535], NaN → 0.
        inline uint16_t toUint16(double x)
        {
            if (!(x >= 0.5))  return 0;        // also catches NaN
            if (x >= 65534.5) return 65535;
            return static_cast<uint16_t>(x + 0.5);
        }

        // Short text form of a number, for string header values ("0.05", "1e-05").
        std::string number(double v)
        {
            std::ostringstream ss;
            ss << std::setprecision(6) << v;
            return ss.str();
        }

        // FITS section "[x1:x2,y1:y2]" (1-based, inclusive) of [x0, x0+w) × [y0, y0+h).
        std::string section(size_t x0, size_t y0, size_t w, size_t h)
        {
            return "[" + std::to_string(x0 + 1) + ":" + std::to_string(x0 + w) + "," +
                         std::to_string(y0 + 1) + ":" + std::to_string(y0 + h) + "]";
        }

        // Median of v: middle value (odd size) or mean of the two middle values (even
        // size); NaN if v is empty. v is taken by copy since nth_element reorders it.
        template <typename T>
        double median(std::valarray<T> v)
        {
            const size_t n = v.size();
            if (n == 0)
                return std::numeric_limits<double>::quiet_NaN();

            T* first = std::begin(v);
            T* mid   = first + n / 2;
            std::nth_element(first, mid, std::end(v));

            if (n % 2)
                return static_cast<double>(*mid);
            return 0.5 * (static_cast<double>(*mid) + static_cast<double>(*std::max_element(first, mid)));
        }

        // Shared with the detector effects (DetectorTools.h).
        using phot2pix::parallelFor;
        using phot2pix::withCubeData;

        // v − ref stored in v's type: as is for floating point; for integer types rounded
        // to the nearest value (not truncated) and clamped to the type's range. The largest
        // value of an unsigned type (65535 for uint16) marks a saturated read: left unchanged.
        template <typename T>
        T subtractReference(T v, double ref)
        {
            if constexpr (std::is_floating_point_v<T>)
                return static_cast<T>(static_cast<double>(v) - ref);
            else
            {
                constexpr T lo = std::numeric_limits<T>::lowest();
                constexpr T hi = std::numeric_limits<T>::max();

                if constexpr (std::is_unsigned_v<T>)
                {
                    if (v == hi)
                        return v;
                }

                const double x = std::round(static_cast<double>(v) - ref);
                if (!(x > static_cast<double>(lo))) return lo;   // also NaN
                if (x >= static_cast<double>(hi))   return hi;
                return static_cast<T>(x);
            }
        }

        // On / off switch read from cfg: true / false, or "ON" / "OFF" ("TRUE" / "FALSE"
        // accepted too), in any case. Throws std::invalid_argument for anything else.
        bool toSwitch(const phot2pix::PixelConfig& cfg, const std::string& key)
        {
            try { return cfg.get<bool>(key); }
            catch (const std::bad_variant_access&) {}

            std::string s;
            try { s = cfg.get<std::string>(key); }
            catch (const std::bad_variant_access&) {}

            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });

            if (s == "ON"  || s == "TRUE")  return true;
            if (s == "OFF" || s == "FALSE") return false;
            throw std::invalid_argument("SCA::initialize: '" + key + "' must be true / false or \"ON\" / \"OFF\".");
        }
    }

    // One pass over memory: accumulate the charge and, for a read, write the measured
    // value straight into the output layer. A std::slice / std::gslice version would need
    // temporaries for the double → uint16 conversion and run on a single thread.
    // The readout noise map is in electrons: it is converted with each pixel's gain,
    // like the read noise in pixel::read(). The RTN map is in ADU, a voltage offset: added as is.
    // With IPC or crosstalk, a pixel needs the charge of others: all the pixels are accumulated
    // first, and read in a second pass.
    void SCA::_accumulate(const std::valarray<double>& charge, std::valarray<double>& cumul,
                          uint16_t* layer, const std::valarray<double>* oneOverF,
                          const std::valarray<double>* rtn) const
    {

        // Measured value of pixel i for the given signal (ADU).
        const auto readPixel = [this, oneOverF, rtn](size_t i, double signal)
        {
            double value = m_pixels[i].read(signal, m_nl);
            if (oneOverF)
                value += (*oneOverF)[i] / m_pixels[i].getGain();   // e- → ADU
            if (rtn)
            {
                value += (*rtn)[i];                                // ADU
                value = (value < 0.0) ? 0.0 : value;               // clamp to [0, 65535]
            }
            return toUint16(value);
        };

        const bool coupled = layer && _pixelsAreCoupled();

        parallelFor(cumul.size(), getThreadCount(), [&](size_t begin, size_t end)
        {
            if (!layer || coupled)          // dropped frames, or read below: charge only
                for (size_t i = begin; i < end; ++i)
                    cumul[i] += charge[i];
            else
                for (size_t i = begin; i < end; ++i)
                {
                    cumul[i] += charge[i];
                    layer[i]  = readPixel(i, cumul[i]);
                }
        });

        if (!coupled)
            return;

        // IPC and crosstalk act on a copy, at each read: the charge in the pixels is not affected.
        std::valarray<double> signal = cumul;
        _couplePixels(signal, 0, 0, getWidth(), getHeight());

        parallelFor(signal.size(), getThreadCount(), [&](size_t begin, size_t end)
        {
            for (size_t i = begin; i < end; ++i)
                layer[i] = readPixel(i, signal[i]);
        });
    }

    void SCA::_couplePixels(std::valarray<double>& signal, size_t c0, size_t r0, size_t ncols, size_t nrows) const
    {
        applyIPC(signal, c0, r0, ncols, nrows);
        applyCrosstalk(signal, c0, r0, ncols, nrows);
    }

    void SCA::applyCrosstalk(std::valarray<double>& signal, size_t c0, size_t r0, size_t ncols, size_t nrows) const
    {
        if (!m_crosstalk.isEnabled())
            return;

        if (c0 + ncols > getWidth() || r0 + nrows > getHeight())
            throw std::out_of_range("SCA::applyCrosstalk: the window exceeds the matrix.");
        if (signal.size() < ncols * nrows)
            throw std::invalid_argument("SCA::applyCrosstalk: the image is smaller than the window.");

        // Signal above the baseline output by every pixel (ADU), then the ghosts it gives. The
        // baseline is removed: it would add an offset that depends on the channel. Full well and
        // non-linearity as in pixel::read(): a saturated pixel cannot give more than its full well.
        const size_t width = getWidth();
        std::valarray<double> ghosts(0.0, ncols * nrows);
        parallelFor(ncols * nrows, getThreadCount(), [&](size_t begin, size_t end)
        {
            for (size_t i = begin; i < end; ++i)
            {
                const nisp::pixel& px = m_pixels[(r0 + i / ncols) * width + c0 + i % ncols];
                const double       Q  = px.storedCharge((signal[i] - px.getBaseline()) * px.getGain());   // ADU → e-
                ghosts[i] = px.nlResponse(Q, m_nl) / px.getGain();                                        // e- → ADU
            }
        });

        m_crosstalk.apply(ghosts, c0, r0, ncols, nrows, getThreadCount());

        parallelFor(ncols * nrows, getThreadCount(), [&](size_t begin, size_t end)
        {
            for (size_t i = begin; i < end; ++i)
                signal[i] += ghosts[i];
        });
    }

    void SCA::_read(std::valarray<double>& charge, double tStart) const
    {
        // The full frame is the window covering the whole matrix.
        _read(charge, 0, 0, getWidth(), getHeight(), tStart);
    }

    void SCA::_read(std::valarray<double>& charge, size_t c0, size_t r0, size_t ncols, size_t nrows, double tStart) const
    {
        _couplePixels(charge, c0, r0, ncols, nrows);

        // The windowed generate() still computes the full channel time series for each
        // overlapping channel (AR(1) correlation from frame start), but only returns
        // the ncols × nrows window portion.
        std::future<std::valarray<double>> oneOverF;
        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
            oneOverF = std::async(std::launch::async,
                [this, c0, r0, ncols, nrows, tStart]() { return m_oneOverFNoise.generate(getWidth(), getHeight(), c0, r0, ncols, nrows, tStart); });

        // Split by pixel, not by row: parallelFor() would run a window (or a frame) of few
        // rows on a single thread. i: pixel of the window, row-major.
        const size_t width = getWidth();
        parallelFor(ncols * nrows, getThreadCount(),
            [this, &charge, c0, r0, ncols, width](size_t begin, size_t end)
            {
                for (size_t i = begin; i < end; ++i)
                    charge[i] = m_pixels[(r0 + i / ncols) * width + c0 + i % ncols].read(charge[i], m_nl);
            });

        if (oneOverF.valid())
        {
            const std::valarray<double> noiseMap = oneOverF.get();   // ncols × nrows, electrons
            parallelFor(ncols * nrows, getThreadCount(),
                [this, &charge, &noiseMap, c0, r0, ncols, width](size_t begin, size_t end)
                {
                    for (size_t i = begin; i < end; ++i)
                        charge[i] += noiseMap[i] / m_pixels[(r0 + i / ncols) * width + c0 + i % ncols].getGain();   // e- → ADU
                });
        }
    }



    void SCA::_checkResponse(const nonlinearity::functor& f) const
    {
        if (!f)
            return;

        constexpr int nSteps = 32;   // charges tested per pixel, from 0 to the full well
        std::atomic<size_t> nBad{0};

        parallelFor(size(), getThreadCount(),
            [this, &f, &nBad](size_t begin, size_t end)
            {
                for (size_t i = begin; i < end; ++i)
                {
                    const double fullWell = m_pixels[i].getFullWell();
                    double previous = 0.0;
                    bool   good     = std::isfinite(fullWell);
                    for (int k = 1; good && k <= nSteps; ++k)
                    {
                        const double S = m_pixels[i].nlResponse(fullWell * k / nSteps, f);
                        good     = std::isfinite(S) && S > previous;
                        previous = S;
                    }
                    if (!good)
                        ++nBad;
                }
            });

        if (nBad > 0)
            throw std::invalid_argument("SCA: the non-linearity model needs a finite \"fullWell\" and must increase up to it: "
                                        + std::to_string(nBad.load()) + " of " + std::to_string(size()) + " pixels do not comply.");
    }

    void SCA::initialize(const phot2pix::PixelConfig& cfg)
    {
        // Also configures the geometry (nisp::Geometry: "detectorSize", "detectorOrigin",
        // "refBorder", "colPerChannel") before loading the maps, and throws if the matrix does not fit.
        phot2pix::DetectorSimulator<nisp::pixel>::initialize(cfg);

        // Placement of this matrix on the physical detector ((0, 0) by default), and output channels.
        const auto [x0, y0] = getOrigin();
        m_oneOverFNoise.setGeometry(_geometry());
        m_oneOverFNoise.setColPerChannel(_geometry().colPerChannel);

        // DPU flux-estimation parameters, written as GAIN / RN in the output header:
        // explicit "dpu_gain" (e-/ADU) / "dpu_noise" (e-) if given, otherwise the mean
        // pixel gain / read noise (uniform value or FITS map), so the header matches
        // the simulated detector.
        if (cfg.has("dpu_gain"))
            dpuGain = cfg.getAsDouble("dpu_gain");
        else
        {
            double sum = 0.0;
            for (size_t i = 0; i < size(); ++i)
                sum += at(i).getGain();
            dpuGain = sum / static_cast<double>(size());
        }

        if (cfg.has("dpu_noise"))
            dpuNoise = cfg.getAsDouble("dpu_noise");
        else
        {
            double sum = 0.0;
            for (size_t i = 0; i < size(); ++i)
                sum += at(i).getReadNoise();
            dpuNoise = sum / static_cast<double>(size());
        }

        if (cfg.has("offset"))
            dpuOffset = cfg.getAsDouble("offset");
        else if (cfg.has("dpu_offset"))
            dpuOffset = cfg.getAsDouble("dpu_offset");

        if (cfg.has("scale"))
            dpuScale = cfg.getAsDouble("scale");
        else if (cfg.has("dpu_scale"))
            dpuScale = cfg.getAsDouble("dpu_scale");

        if (cfg.has("fr_time"))
            frameTime = cfg.getAsDouble("fr_time");
        else if (cfg.has("frame_time"))
            frameTime = cfg.getAsDouble("frame_time");

        // The readout noise needs the frame time to place each read in time (pedestal drift).
        m_oneOverFNoise.setFrameTime(frameTime);

        if (cfg.has("exposureMode"))
        {
            try                                     { m_mode = modeFromString(cfg.get<std::string>("exposureMode")); }
            catch (const std::invalid_argument& err) { throw std::invalid_argument(std::string("SCA::initialize: ") + err.what()); }
        }

        // Reference pixel correction of the DPU mode (Eq. 3, mode ALL): on by default.
        if (cfg.has("refPixCorr"))
            m_refPixCorr = toSwitch(cfg, "refPixCorr");

        if(cfg.has("persistence"))
        {
            m_enable_persistence |= true;
            if (!persistenceModelling)
                persistenceModelling = std::make_shared<persistence::scheduler>(cfg);
        }

        // The persistence signal is added to ADU frames: the traps of each pixel get its gain.
        _setPersistenceGain();

        // Readout noise added at each read: per-channel 1/f ("oneOverF"), 1/f common to
        // every channel ("bias_voltage_noise") and channel offsets ("pedestal_drift").
        if(cfg.has("oneOverF") || cfg.has("bias_voltage_noise") || cfg.has("pedestal_drift"))
        {
            m_enable_oneOverF |= true;
            m_oneOverFNoise.initialize(cfg);
        }

        // Fail here rather than at the first exposure if the channel layout does not
        // match the detector (colPerChannel, nChannels).
        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
            m_oneOverFNoise.checkGeometry(getWidth(), getHeight());

        // Random telegraph noise of the source followers ("rtn"): the RTN pixels, their time
        // constants and the sign of their offset are drawn here, once.
        if (cfg.has("rtn"))
            m_rtn.initialize(cfg, getWidth(), getHeight());

        // Reference pixels: the physical pixels within refBorder (4) of the detector edge
        // (nisp::Geometry). A 2048 × 2048 matrix has 4 reference pixels on each side; a
        // smaller one has those of the part of the detector it covers (at the default
        // origin (0, 0): 4 columns on the left and 4 rows at the bottom; a 2040 × 2040
        // one at [4, 4]: none). They have no photodiode: QE = 0. Their dark current
        // comes from "dark" (value or map) unless "refDark" (e-/s) is given.
        if (_geometry().refBorder > 0)
        {
            const nisp::Geometry& geometry = _geometry();
            std::optional<double> refDark;
            if (cfg.has("refDark"))
                refDark = cfg.getAsDouble("refDark");

            const size_t W = getWidth(), H = getHeight();
            for (size_t row = 0; row < H; ++row)
                for (size_t col = 0; col < W; ++col)
                    if (geometry.isReference(x0 + col, y0 + row))
                    {
                        nisp::pixel& px = m_pixels[row * W + col];
                        px.setQuantumEfficiency(0.0);
                        if (refDark)
                            px.setDarkCurrent(*refDark);   // throws if < 0
                    }
        }

        // Full well "fullWell" (e-), and non-linearity: model "nlModel" with the parameters
        // "nl0", "nl1", … of each pixel. Each number may be a FITS map. No full well and a
        // linear response by default.
        {
            const std::valarray<double> fullWell = loadMap(cfg, "fullWell");
            for (size_t i = 0; i < fullWell.size(); ++i)
                m_pixels[i].setFullWell(static_cast<float>(fullWell[i]));   // throws if not > 0

            for (size_t k = 0; k < nonlinearity::nParam; ++k)
            {
                const std::valarray<double> param = loadMap(cfg, "nl" + std::to_string(k));
                for (size_t i = 0; i < param.size(); ++i)
                    m_pixels[i].setNlParam(k, static_cast<float>(param[i]));
            }

            if (cfg.has("nlModel"))
                setNonLinearity(nonlinearity::fromString(cfg.get<std::string>("nlModel")));
            else
                _checkResponse(m_nl);   // a model already set, with the new values
        }

        // Inter-pixel capacitance kernel "ipc": none by default.
        m_ipc.initialize(cfg, getGeometry(), getWidth(), getHeight(), [this, &cfg]() { return loadImage(cfg, "ipc"); }, "ipc");

        // Crosstalk between the output channels "crosstalk": none by default.
        _configureCrosstalk(cfg);

        _recordInputs(cfg);

        m_cosmicDeposit.configure(cfg);

        // Cosmic ray inputs in the output header (CRFLUX, CRACC, ...): only with a "cosmic" object.
        if (m_cosmicDeposit.isConfigured())
            for (const auto& card : m_cosmicDeposit.inputCards())
                setInputCard(card);
    }

    void SCA::_configureCrosstalk(const phot2pix::PixelConfig& cfg)
    {
        if (!cfg.has("crosstalk"))
        {
            // "colPerChannel" may have changed: same couplings, with the new channels.
            if (m_crosstalk.isEnabled() && m_crosstalkColumns != _geometry().colPerChannel)
                _setCrosstalk(m_crosstalkCouplings);
            return;
        }

        const std::string what = "SCA::initialize: 'crosstalk'";

        const phot2pix::json* node = nullptr;
        try                                { node = &cfg.getJson("crosstalk"); }
        catch (const std::invalid_argument&) {}   // not a JSON value
        if (!node || !(node->is_array() || (node->is_object() && node->contains("kernel"))))
            throw std::invalid_argument(what + " must be an array of couplings, or {\"kernel\": [...]}.");

        // ── No crosstalk: [] ──
        if (node->is_array() && node->empty())
        {
            m_crosstalk        = phot2pix::PixelConvolution();
            m_crosstalkColumns = 0;
            m_crosstalkCouplings.clear();
            return;
        }

        if (node->is_object())
        {
            for (auto it = node->begin(); it != node->end(); ++it)
            {
                if (it.key() == "nChannels" || it.key() == "colPerChannel")
                    throw std::invalid_argument(what + ": '" + it.key() + "' is not a key of crosstalk: "
                                                "set it at the top level of the configuration.");
                if (it.key() != "kernel")
                    throw std::invalid_argument(what + ": unknown key '" + it.key() + "'.");
            }
        }

        // ── Couplings, one per distance between channels ──
        const phot2pix::json& values = node->is_object() ? node->at("kernel") : *node;
        if (!values.is_array() || values.empty())
            throw std::invalid_argument(what + ": the kernel must be an array of 1 to " + std::to_string(_geometry().nChannels()) +
                                        " couplings (at most one per channel).");

        std::vector<double> couplings;
        for (size_t d = 0; d < values.size(); ++d)
        {
            const double k = values[d].is_number() ? values[d].get<double>() : -1.0;
            if (!(k >= 0.0) || !std::isfinite(k))
                throw std::invalid_argument(what + ": coupling " + std::to_string(d) + " must be a finite number >= 0.");
            if (d == 0 && k != 1.0)
                throw std::invalid_argument(what + ": coupling 0, the fraction of its signal that a channel keeps, must be 1 (no loss).");
            couplings.push_back(k);
        }

        _setCrosstalk(couplings);
    }

    void SCA::_setCrosstalk(const std::vector<double>& couplings)
    {
        const std::string what    = "SCA::initialize: 'crosstalk'";
        const size_t      columns = _geometry().colPerChannel, channels = _geometry().nChannels();
        if (channels * columns != _geometry().width)
            throw std::invalid_argument(what + ": " + std::to_string(channels) + " channels of " + std::to_string(columns) +
                                        " columns do not cover the " + std::to_string(_geometry().width) + " columns of the detector.");
        if (couplings.size() > channels)
            throw std::invalid_argument(what + ": " + std::to_string(couplings.size()) + " couplings for " +
                                        std::to_string(channels) + " channels.");

        // Kernel of the ghosts: one row, coupling d at ± d channels from the centre, 0 elsewhere
        // (and in the centre: the source keeps its signal, the ghosts are added to it). The
        // couplings that are not given, at larger distances, are 0.
        const size_t centre = (couplings.size() - 1) * columns;
        std::valarray<double> weights(0.0, 2 * centre + 1);
        for (size_t d = 1; d < couplings.size(); ++d)
        {
            weights[centre - d * columns] = couplings[d];
            weights[centre + d * columns] = couplings[d];
        }

        if (weights.max() > 0.0)
            m_crosstalk.setKernel(weights, weights.size(), 1, getWidth(), getHeight(), false);   // ghosts only: not normalized
        else
            m_crosstalk = phot2pix::PixelConvolution();   // every coupling is 0: no crosstalk
        m_crosstalkColumns   = columns;
        m_crosstalkCouplings = couplings;
    }

    void SCA::_recordInputs(const phot2pix::PixelConfig& cfg)
    {
        // ── Detector geometry and reference pixels ──
        const nisp::Geometry& geometry = _geometry();
        setInputCard({"DETSIZE", section(0, 0, geometry.width, geometry.height), "Physical detector size"});
        setInputCard({"REFBORD", geometry.refBorder, "[pix] Reference pixels on each side of the detector"});
        if (cfg.has("refDark"))
            setInputCard({"REFDARK", cfg.getAsDouble("refDark"), "[e-/s] Dark current of the reference pixels"});

        // ── Full well and non-linearity: the keys that are given ──
        if (const auto value = cfg.card("fullWell"))
            setInputCard({"PIXFW", *value, "[e-] Pixel full well"});
        if (cfg.has("nlModel"))
            setInputCard({"NLMODEL", cfg.get<std::string>("nlModel"), "Non-linearity model"});
        for (size_t k = 0; k < nonlinearity::nParam; ++k)
            if (const auto value = cfg.card("nl" + std::to_string(k)))
                setInputCard({"PIXNL" + std::to_string(k), *value, "Parameter of the non-linearity model"});

        // ── Channel crosstalk: the couplings that are not 0 ──
        removeInputCards("XTALK");
        if (m_crosstalk.isEnabled())
        {
            const std::valarray<double>& weights = m_crosstalk.getKernel().data();
            const size_t centre      = weights.size() / 2;
            const size_t maxDistance = centre / m_crosstalkColumns;   // in channels
            setInputCard({"XTALKCOL", m_crosstalkColumns,                  "[pix] Crosstalk: columns per output channel"});
            setInputCard({"XTALKNCH", geometry.width / m_crosstalkColumns, "Crosstalk: output channels"});
            for (size_t d = 1; d <= maxDistance; ++d)
                if (const double k = weights[centre + d * m_crosstalkColumns]; k > 0.0)
                    setInputCard({"XTALK" + std::to_string(d), k, "Crosstalk: coupling at this distance in channels"});
        }

        // ── Readout noise: the components that are on ──
        for (const char* prefix : {"CHANWID", "OOF", "BIAS", "PED"})
            removeInputCards(prefix);

        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
        {
            const OneOverFNoise& n = m_oneOverFNoise;
            setInputCard({"CHANWID", n.getColPerChannel(), "[pix] Columns per output channel"});

            if (n.getAmplitude() > 0.0)
            {
                setInputCard({"OOFAMP",   n.getAmplitude(), "[e-] 1/f noise of each channel: amplitude"});
                setInputCard({"OOFSLOPE", n.getSlope(),     "1/f noise of each channel: PSD slope"});
                setInputCard({"OOFPOLES", n.getNPoles(),    "1/f noise of each channel: AR(1) poles"});
            }

            if (n.getBiasVoltageNoise() > 0.0)
            {
                setInputCard({"BIASAMP",  n.getBiasVoltageNoise(), "[e-] Bias voltage noise (all channels): amplitude"});
                setInputCard({"BIASSLOP", n.getBiasSlope(),        "Bias voltage noise: PSD slope"});
                setInputCard({"BIASPOLE", n.getBiasNPoles(),       "Bias voltage noise: AR(1) poles"});
            }

            if (n.getPedestalSigma() > 0.0 || n.hasPedestalDrift())
            {
                setInputCard({"PEDSIGMA", n.getPedestalSigma(), "[e-] Pedestal noise per channel and frame"});
                setInputCard({"PEDT0",    n.getPedestalT0(),    "[s] Time origin of the pedestal drift"});

                if (!n.getDriftPoly().empty())
                {
                    std::string poly;
                    for (double a : n.getDriftPoly())
                        poly += (poly.empty() ? "" : ", ") + number(a);
                    setInputCard({"PEDPOLY", poly, "[e-, e-/s, ...] Pedestal drift: sum a_k t^k"});
                }

                size_t k = 0;
                for (const auto& e : n.getDriftExp())
                    setInputCard({"PEDEXP" + std::to_string(++k),
                                  "A=" + number(e.amplitude) + " tau=" + number(e.tau),
                                  "[e-, s] Pedestal drift: A exp(-t/tau)"});

                k = 0;
                for (const auto& s : n.getDriftSin())
                    setInputCard({"PEDSIN" + std::to_string(++k),
                                  "B=" + number(s.amplitude) + " P=" + number(s.period) + " phi=" + number(s.phase),
                                  "[e-, s, rad] Pedestal drift: B sin(2 pi t/P + phi)"});
            }
        }

        // ── Random telegraph noise ──
        removeInputCards("RTN");
        if (m_rtn.isEnabled())
        {
            setInputCard({"RTNFRAC", m_rtn.getFraction(),         "RTN: fraction of the pixels"});
            setInputCard({"RTNNPIX", m_rtn.getPixels().size(),    "[pix] RTN: pixels drawn"});
            setInputCard({"RTNTMIN", m_rtn.getTauMin(),           "[s] RTN: shortest time constant"});
            setInputCard({"RTNTMAX", m_rtn.getTauMax(),           "[s] RTN: longest time constant"});
            setInputCard({"RTNAMP",  m_rtn.getAmplitude(),        "[ADU] RTN: largest amplitude, uniform from 1 ADU"});
            setInputCard({"RTNPOS",  m_rtn.getPositiveFraction(), "RTN: fraction of positive offsets"});
            if (const auto seed = m_rtn.getSeed())
                setInputCard({"RTNSEED", static_cast<size_t>(*seed), "RTN: seed of the draw"});
        }

        // ── Persistence: [{"f", "t"}, …] or {"traps": [...], "dT": s} ──
        if (cfg.has("persistence"))
        {
            removeInputCards("PERTR");
            removeInputCards("PERDT");

            const auto& p     = cfg.getJson("persistence");
            const auto& traps = (p.is_object() && p.contains("traps")) ? p.at("traps") : p;

            if (p.is_object() && p.contains("dT") && p.at("dT").is_number())
                setInputCard({"PERDT", p.at("dT").get<double>(), "[s] Persistence: delay between two exposures"});

            size_t k = 0;
            if (traps.is_array())
                for (const auto& trap : traps)
                    if (trap.is_object())
                        setInputCard({"PERTR" + std::to_string(++k),
                                      "f=" + number(trap.value("f", 0.0)) + " tau=" + number(trap.value("t", 0.0)),
                                      "[1/s, s] Persistence trap: release rate, decay time"});
        }
    }

    void SCA::_writeHeader(DSL::FITScube& image, size_t c0, size_t r0, size_t ncols, size_t nrows) const
    {
        auto& hdu = image.HDU();

        hdu.ValueForKey("OFFSET",  dpuOffset, "[ADU] DPU offset: DN = value x SFACTOR + OFFSET");
        hdu.ValueForKey("SFACTOR", dpuScale,  "DPU scale factor: DN = value x SFACTOR + OFFSET");
        hdu.ValueForKey("GAIN",    dpuGain,   "[e-/ADU] Mean gain used by the DPU");
        hdu.ValueForKey("RN",      dpuNoise,  "[e-] Mean read noise used by the DPU");
        hdu.ValueForKey("FRTIME",  frameTime, "[s] NISP HAWAII-2RG frame time");
        hdu.ValueForKey("PHOT2PIX", std::string("v")+std::to_string(PHOT2PIX_VERSION_MAJOR)
                                    +std::string(".")+std::to_string(PHOT2PIX_VERSION_MINOR)
                                    +std::string(".")+std::to_string(PHOT2PIX_VERSION_PATCH)
                                    +std::string("-")+PHOT2PIX_VERSION, DSL::key_type::fChar, "phot2pix version");

        // Where this image sits on the physical detector (NOAO DETSEC convention).
        try
        {
            const auto [x0, y0] = getOrigin();
            hdu.ValueForKey("DETSEC", section(x0 + c0, y0 + r0, ncols, nrows), "Detector section of this image");
        }
        catch (const std::invalid_argument&) {}   // matrix larger than the default detector (initialize() not called)

        writeInputCards(image);
    }

    std::shared_ptr<DSL::FITScube> SCA::_toUint16(const std::shared_ptr<DSL::FITScube>& doubleFrame) const
    {
        const std::valarray<double>& data = *(doubleFrame->GetData<double>());
        std::valarray<uint16_t> out((uint16_t) 0, data.size());
        std::transform(std::begin(data), std::end(data), std::begin(out),
            [this](double x)
            {
                const double v = x * dpuScale + dpuOffset + 0.5;
                if (v < 0.0)         return (uint16_t) 0;
                if (v >= 65535.0)    return (uint16_t) 65535;
                return static_cast<uint16_t>(v);
            }
        );

        auto img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<uint16_t>({doubleFrame->Size(1), doubleFrame->Size(2)}));
        *(img) += out;
        return img;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose(
        const std::valarray<double>& fullFlux, double exposureTime,
        size_t c0, size_t r0, size_t ncols, size_t nrows)
    {

        std::shared_ptr<DSL::FITScube> frame;

        m_cosmicCharge.resize(0);   // cosmic hits of this exposure only (persistence)
        _simulate(fullFlux, exposureTime, frame, c0, r0, ncols, nrows);
        // Single read: read noise + readout noise, the read being the last frame of the exposure.
        _read(*(frame->GetData<double>()), c0, r0, ncols, nrows, std::max(0.0, exposureTime - frameTime));
        // Random telegraph noise of the read (ADU), levels drawn for this read alone.
        if (m_rtn.isEnabled())
            m_rtn.beginStream().addTo(*(frame->GetData<double>()), c0, r0, ncols, nrows);

        // Register full-frame charge (electrons) for future persistence computation.
        if (persistenceModelling && m_enable_persistence)
        {
            auto chargeImg = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<double>({getWidth(), getHeight()}));

            std::valarray<double> chargeMap(0.0, size());
            for (size_t i = 0; i < size(); ++i)
                chargeMap[i] = fullFlux[i] * at(i).getQuantumEfficiency() * exposureTime;
            if (m_cosmicCharge.size() == chargeMap.size())
                chargeMap += m_cosmicCharge;   // + charge of the cosmic hits (e-)
            
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        std::shared_ptr<DSL::FITScube> result = _toUint16(frame);
        result->HDU().ValueForKey("EXPTIME", exposureTime, "[s] Exposure time, reset to read");
        _writeHeader(*result, c0, r0, ncols, nrows);
        return result;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose(const std::valarray<double>& fullFlux, const nisp::macc& maccParams)
    {
        switch (m_mode)
        {
            case mode::FLUX:
                return expose_flux(fullFlux, (maccParams.nGroups*(maccParams.nFrames+maccParams.nDrop)-maccParams.nDrop) * frameTime);
                break;

            case mode::CDS:
                return expose_cds(fullFlux, (maccParams.nGroups*(maccParams.nFrames+maccParams.nDrop)-maccParams.nDrop) * frameTime);

            case mode::DPU:
                return expose_dpu(fullFlux, maccParams);
                break;
            
            default:
                return expose_macc(fullFlux, maccParams);
                break;
        }
    }

    std::shared_ptr<DSL::FITScube> SCA::expose(const std::valarray<double>& fullFlux, const double& exposureTime)
    {
        switch (m_mode)
        {
            case mode::MACC:
            case mode::UTR:
            case mode::DPU:
                {
                    macc maccParams = {static_cast<size_t>(exposureTime / frameTime), 1, 0};
                    if(m_mode == mode::DPU)
                        return expose_dpu(fullFlux, maccParams);
                    else    
                        return expose_macc(fullFlux, maccParams);
                }
                break;

            case mode::CDS:
                return expose_cds(fullFlux, exposureTime);

            default:
                return expose_flux(fullFlux, exposureTime);
                break;
        }
    }

    std::shared_ptr<DSL::FITScube> SCA::expose(const std::valarray<double>& fullFlux, const nisp::readout& r)
    {
        // The mode of r is for this exposure only.
        const nisp::mode scaMode = m_mode;
        if (r.exposureMode)
            m_mode = *r.exposureMode;

        std::shared_ptr<DSL::FITScube> result;
        try
        {
            result = r.maccParams ? expose(fullFlux, *r.maccParams) : expose(fullFlux, r.exposureTime);
        }
        catch (...)
        {
            m_mode = scaMode;
            throw;
        }

        m_mode = scaMode;
        return result;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose_macc(const std::valarray<double>& fullFlux, const nisp::macc& maccParams)
    {
        const size_t nG      = maccParams.nGroups;
        const size_t nF      = maccParams.nFrames;
        const size_t nD      = maccParams.nDrop;
        const size_t nPix    = size();
        const size_t nLayers = nG * nF;

        if(nLayers == 0)
            throw std::invalid_argument("SCA::expose: MACC needs nGroups >= 1 and nFrames >= 1.");

        double exposureTime = 0;
        double integrationTime = frameTime;

        // Reset level (baseline + kTC) belongs to the reset: drawn in the very first
        // frame of the ramp only (pixel::resetLevel()), switched off for every later increment.
        const phot2pix::PixelConfig firstRead;
        phot2pix::PixelConfig cfg;
        cfg.set("baseline", std::string("off"));

        // Output ramp allocated once (uint16, W × H × nLayers). Each read is written
        // straight into its own layer: no AddLayer() regrowth of the whole cube,
        // no double-precision copy of the ramp and no final cast<uint16_t>().
        std::vector<size_t> axes = {getWidth(), getHeight()};
        if(nLayers > 1)
            axes.push_back(nLayers);

        auto out = std::make_shared<DSL::FITSimg<uint16_t>>(axes);
        std::valarray<uint16_t>* outData = out->GetData<uint16_t>();
        if(!outData)
            throw std::runtime_error("SCA::expose: cannot allocate the output ramp.");

        // Charge accumulated in the pixels (ADU): signal + dark + shot noise + reset level.
        // Read noise and 1/f noise are added to each read's layer only, never to cumul.
        std::valarray<double>          cumul(0.0, nPix);
        m_cosmicCharge.resize(0);   // cosmic hits of this exposure only (persistence)
        std::shared_ptr<DSL::FITScube> frame    = nullptr;
        std::shared_ptr<DSL::FITScube> skipped  = nullptr;
        size_t                         layer    = 0;

        // One continuous 1/f noise stream for the whole ramp: every frame (read or
        // dropped) consumes the next chunk of a single per-channel AR(1) process
        // instead of restarting the amplifier state at zero for each read.
        std::optional<nisp::OneOverFNoise::Stream> oneOverFStream;
        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
            oneOverFStream = m_oneOverFNoise.beginStream(getWidth(), getHeight());

        // Random telegraph noise: the levels of the RTN pixels over the whole ramp, moved forward at
        // every frame (read or dropped), and their offsets (ADU) written in rtnMap before each read.
        std::optional<nisp::RTN::Stream> rtnStream;
        std::valarray<double>            rtnMap;
        if (m_rtn.isEnabled())
        {
            rtnStream = m_rtn.beginStream();
            rtnMap.resize(nPix, 0.0);
        }

        for (size_t g = 0; g < nG; ++g)
        {
            for(size_t f = 0; f < nF; ++f)
            {
                integrationTime = frameTime;
                exposureTime   += integrationTime;

                // This read's 1/f chunk, generated concurrently with the pixel charge.
                std::future<std::valarray<double>> oneOverFFuture;
                if (oneOverFStream)
                    oneOverFFuture = std::async(std::launch::async,
                        [&stream = *oneOverFStream]() { return stream.next(); });
                
                _simulate(fullFlux, integrationTime, frame, ((layer == 0) ? firstRead : cfg));
                

                if(!frame)
                    throw std::runtime_error("SCA::expose: failed to simulate frame ["+std::to_string(g)+","+std::to_string(f)+"]");

                std::valarray<double> oneOverF;
                if (oneOverFFuture.valid())
                    oneOverF = oneOverFFuture.get();

                // RTN levels at this read: one frame after the previous frame, read or dropped.
                if (rtnStream)
                {
                    if (layer > 0)
                        rtnStream->advance(frameTime);
                    rtnStream->writeTo(rtnMap);
                }

                // cumul += this frame's charge, then read it: layer = pixel read (read noise) + 1/f + RTN.
                _accumulate(*(frame->GetData<double>()), cumul, &(*outData)[layer * nPix],
                            oneOverF.size() ? &oneOverF : nullptr, m_rtn.isEnabled() ? &rtnMap : nullptr);
                frame.reset();
                ++layer;
            }

            if(nD <= 0)
                continue;

            integrationTime = frameTime * nD;
            exposureTime += integrationTime;

            if(g+1 < nG)
            {
                // Dropped frames are read and discarded: their nD 1/f chunks are consumed
                // so the amplifier process stays continuous in time, but never stored.
                std::future<void> oneOverFSkip;
                if (oneOverFStream)
                    oneOverFSkip = std::async(std::launch::async,
                        [&stream = *oneOverFStream, nD]() { for (size_t k = 0; k < nD; ++k) stream.next(); });

                _simulate(fullFlux, integrationTime, skipped, cfg);
                if(!skipped)
                    throw std::runtime_error("SCA::expose: failed to simulate skipped frame ["+std::to_string(g)+","+std::to_string(nF-1)+",0]");

                _accumulate(*(skipped->GetData<double>()), cumul, nullptr, nullptr, nullptr);   // dropped frames: charge only, no layer
                skipped.reset();

                // The RTN levels go on changing during the dropped frames.
                if (rtnStream)
                    rtnStream->advance(static_cast<double>(nD) * frameTime);

                if (oneOverFSkip.valid())
                    oneOverFSkip.get();
            }
        }

        cumul.resize(0);

        // Register accumulated charge (electrons) for future persistence computation.
        if (persistenceModelling && m_enable_persistence)
        {
            auto chargeImg = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<double>({getWidth(), getHeight()}));
            std::valarray<double> chargeMap(0.0, size());
            for (size_t i = 0; i < size(); ++i)
                chargeMap[i] = fullFlux[i] * at(i).getQuantumEfficiency() * integrationTime;
            if (m_cosmicCharge.size() == chargeMap.size())
                chargeMap += m_cosmicCharge;   // + charge of the cosmic hits (e-)
            
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        std::shared_ptr<DSL::FITScube> result = out;
        result->HDU().ValueForKey("INTTIME", static_cast<double>(nG*(nF+nD)-nD)*frameTime, "[s] Integration time, reset to last read");
        result->HDU().ValueForKey("EXPTIME", exposureTime, "[s] Exposure time");
        result->HDU().ValueForKey("MODE", std::string("MACC"), "Simulation mode");
        result->HDU().ValueForKey("NG", nG, "Number of groups");
        result->HDU().ValueForKey("NR", nF, "Number of reads/group");
        result->HDU().ValueForKey("ND", nD, "Number of drops/group");
        if(m_enable_persistence)
            result->HDU().ValueForKey("PERSIST", std::string("ON"), "Persistence modelling enabled");

        _writeHeader(*result);

        return result;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose_cds(const std::valarray<double>& fullFlux, const double& exposureTime)
    {
        if(exposureTime < 2*frameTime)
            throw std::invalid_argument("SCA::expose_cds: exposure time must be at least 2*"+std::to_string(frameTime)+" seconds for CDS mode.");

        // Small margin so that rounding (75.9999… instead of 76) cannot drop a frame.
        const double nFrames          = std::floor(exposureTime / frameTime + 1e-9);
        const double eff_exposureTime = nFrames * frameTime;            // reset → read N
        const double cdsTime          = eff_exposureTime - frameTime;   // read 1 → read N

        std::shared_ptr<DSL::FITScube> f0 = nullptr;   // read 1
        std::shared_ptr<DSL::FITScube> f1 = nullptr;   // read N

        m_cosmicCharge.resize(0);   // cosmic hits of this exposure only (persistence)

        // Charge at read 1: reset level (baseline + kTC) + one frame of signal.
        _simulate(fullFlux, frameTime, f0);
        if(!f0)
            throw std::runtime_error("SCA::expose_cds: failed to simulate frame [f0]");

        // Charge collected between read 1 and read N (baseline off: no new reset).
        phot2pix::PixelConfig cfg;
        cfg.set("baseline", std::string("off"));

        _simulate(fullFlux, cdsTime, f1, cfg);
        if(!f1)
            throw std::runtime_error("SCA::expose_cds: failed to simulate frame [f1]");

        std::valarray<double>& read1 = *(f0->GetData<double>());
        std::valarray<double>& readN = *(f1->GetData<double>());

        readN += read1;   // charge at read N = charge at read 1 + charge collected since

        // Each read adds its own, independent read noise (pixel::read) and readout noise.
        // Read 1 is frame 0, read N is frame N−1: it starts cdsTime = (N−1)·t_f after the reset.
        // The RTN levels (ADU, added to each read) are those of one process, cdsTime apart.
        std::optional<nisp::RTN::Stream> rtnStream;
        if (m_rtn.isEnabled())
            rtnStream = m_rtn.beginStream();

        _read(read1, 0.0);
        if (rtnStream)
            rtnStream->addTo(read1);

        _read(readN, cdsTime);
        if (rtnStream)
        {
            rtnStream->advance(cdsTime);
            rtnStream->addTo(readN);
        }

        // CDS = read N − read 1: reset level and kTC cancel, the read noise of both reads remains.
        *(f1->GetData<double>()) -= *(f0->GetData<double>());

        f0.reset();

        // Register accumulated charge (electrons) for future persistence computation.
        if (persistenceModelling && m_enable_persistence)
        {
            auto chargeImg = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<double>({getWidth(), getHeight()}));
            std::valarray<double> chargeMap(0.0, size());
            for (size_t i = 0; i < size(); ++i)
                chargeMap[i] = fullFlux[i] * at(i).getQuantumEfficiency() * eff_exposureTime;
            if (m_cosmicCharge.size() == chargeMap.size())
                chargeMap += m_cosmicCharge;   // + charge of the cosmic hits (e-)
            
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        // DPU scale / offset applied like the FLUX mode, so negative CDS values (noise on
        // dark pixels) survive when an offset is set, instead of being clipped to 0.
        std::shared_ptr<DSL::FITScube> result = _toUint16(f1);
        result->HDU().ValueForKey("INTTIME", cdsTime, "[s] Integration time, read 1 to read N");
        result->HDU().ValueForKey("EXPTIME", eff_exposureTime, "[s] Exposure time, reset to read N");
        result->HDU().ValueForKey("MODE", std::string("CDS"), "Simulation mode");
        if(m_enable_persistence)
            result->HDU().ValueForKey("PERSIST", std::string("ON"), "Persistence modelling enabled");

        _writeHeader(*result);

        f1.reset();
        return result;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose_flux(const std::valarray<double>& flux, const double& exposureTime)
    {
        std::shared_ptr<DSL::FITScube> frame;

        double effectiveExposureTime = std::floor(exposureTime/frameTime)*frameTime;

        m_cosmicCharge.resize(0);   // cosmic hits of this exposure only (persistence)
        _simulate(flux, effectiveExposureTime, frame);
        // Single read: read noise + readout noise, the read being the last frame of the exposure.
        _read(*(frame->GetData<double>()), std::max(0.0, effectiveExposureTime - frameTime));
        // Random telegraph noise of the read (ADU), levels drawn for this read alone.
        if (m_rtn.isEnabled())
            m_rtn.beginStream().addTo(*(frame->GetData<double>()));

        // Register accumulated charge (electrons) for future persistence computation.
        if (persistenceModelling && m_enable_persistence)
        {
            auto chargeImg = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<double>({getWidth(), getHeight()}));
            std::valarray<double> chargeMap(0.0, size());
            for (size_t i = 0; i < size(); ++i)
                chargeMap[i] = flux[i] * at(i).getQuantumEfficiency() * effectiveExposureTime;
            if (m_cosmicCharge.size() == chargeMap.size())
                chargeMap += m_cosmicCharge;   // + charge of the cosmic hits (e-)
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        std::shared_ptr<DSL::FITScube> result = _toUint16(frame);
        result->HDU().ValueForKey("INTTIME", effectiveExposureTime, "[s] Integration time, reset to read");
        result->HDU().ValueForKey("EXPTIME", effectiveExposureTime, "[s] Exposure time");
        result->HDU().ValueForKey("MODE", std::string("FLUX"), "Simulation mode");
        if(m_enable_persistence)
            result->HDU().ValueForKey("PERSIST", std::string("ON"), "Persistence modelling enabled");

        _writeHeader(*result);

        frame.reset();
        return result;
    }

    std::valarray<double> SCA::getReferencePixelCorrection(const std::shared_ptr<DSL::FITScube>& ramp, const refPixCorrection& mode,
                                                           const nisp::Geometry& geometry)
    {
        if(!ramp)
            throw std::invalid_argument("SCA::getReferencePixelCorrection: ramp is null.");

        if(mode == refPixCorrection::NONE)
            return {};

        // Detector of W × H pixels, border of b reference pixels, nCh channels of w columns.
        const size_t W = geometry.width, H = geometry.height, b = geometry.refBorder;
        const size_t w = geometry.colPerChannel, nCh = geometry.nChannels();

        if(ramp->Size(1) != W || ramp->Size(2) != H)
        {
            std::cerr << "SCA::getReferencePixelCorrection: warning: reference pixel correction is only applied to the full "
                      << W << "x" << H << " detector." << std::endl;
            return {};
        }

        // Every channel needs columns outside the left / right border (both in a single channel).
        if(b == 0 || nCh * w != W || w <= (nCh == 1 ? 2 * b : b) || H <= 2 * b)
        {
            std::cerr << "SCA::getReferencePixelCorrection: warning: reference pixel correction is not applied with a border of "
                      << b << " reference pixels and channels of " << w << " columns." << std::endl;
            return {};
        }

        const size_t nPix     = W * H;
        const size_t nLines   = H - 2 * b;   // lines that are corrected: [b, H − b)
        const size_t nReads   = (ramp->GetDimension() >= 3) ? ramp->Size(3) : 1;
        const size_t nThreads = std::thread::hardware_concurrency();

        // c[(y * nCh + ch) * nReads + k]: the reads of one line and channel are contiguous.
        // The b first and last lines keep 0 (not corrected).
        std::valarray<double> c(0.0, H * nCh * nReads);

        const bool known = withCubeData(*ramp, [&](const auto& data)
        {
            using T = typename std::decay_t<decltype(data)>::value_type;

            // T^(ch), B^(ch) of every read k, at [k * nCh + ch]: median of the b top / bottom
            // reference lines of the channel, without the left / right reference columns
            // (NISP, channels 0 and 31: 60 columns).
            std::valarray<double> topRef(0.0, nCh * nReads);
            std::valarray<double> bottomRef(0.0, nCh * nReads);

            parallelFor(nCh * nReads, nThreads, [&](size_t begin, size_t end)
            {
                for(size_t j = begin; j < end; ++j)
                {
                    const size_t k  = j / nCh;
                    const size_t ch = j % nCh;
                    const size_t L0 = k * nPix;   // first pixel of read k

                    size_t size = w;
                    size_t x0   = ch * w;
                    if(ch == 0)
                    {
                        x0   += b;
                        size -= b;
                    }
                    if(ch + 1 == nCh)
                        size -= b;

                    bottomRef[j] = median(std::valarray<T>(data[std::gslice(L0 + x0,               {size, b}, {1, W})]));
                    topRef[j]    = median(std::valarray<T>(data[std::gslice(L0 + x0 + (H - b) * W, {size, b}, {1, W})]));
                }
            });

            // For every read k and line y in [b, H − b): L_y, R_y, median of the b left / right
            // reference columns over lines y−4 … y+4 (those of the detector: all of them if
            // b >= 4), then c_y^(ch) of the nCh channels (Eq. 3).
            parallelFor(nLines * nReads, nThreads, [&](size_t begin, size_t end)
            {
                for(size_t j = begin; j < end; ++j)
                {
                    const size_t k  = j / nLines;
                    const size_t y  = b + j % nLines;
                    const size_t L0 = k * nPix;

                    const size_t first = (y > 4) ? y - 4 : 0;
                    const size_t n     = std::min(y + 4, H - 1) - first + 1;

                    const double leftRef  = median(std::valarray<T>(data[std::gslice(L0           + first * W, {b, n}, {1, W})]));
                    const double rightRef = median(std::valarray<T>(data[std::gslice(L0 + (W - b) + first * W, {b, n}, {1, W})]));

                    for(size_t ch = 0; ch < nCh; ++ch)
                    {
                        const double top    = topRef[k * nCh + ch];
                        const double bottom = bottomRef[k * nCh + ch];

                        // mode is a bitmask: each pair of bits picks its own reference,
                        // so every combination (TOP|LEFT, BOTTOM|RIGHTnLEFT, …) is handled.
                        double ref = 0.0;

                        switch(mode & refPixCorrection::TOPnBOTTOM)
                        {
                            case refPixCorrection::TOP:        ref += top;                        break;
                            case refPixCorrection::BOTTOM:     ref += bottom;                     break;
                            case refPixCorrection::TOPnBOTTOM: ref += (top + bottom) / 2.0;       break;
                            default:                                                              break;
                        }

                        switch(mode & refPixCorrection::RIGHTnLEFT)
                        {
                            case refPixCorrection::LEFT:       ref += leftRef;                    break;
                            case refPixCorrection::RIGHT:      ref += rightRef;                   break;
                            case refPixCorrection::RIGHTnLEFT: ref += (rightRef + leftRef) / 2.0; break;
                            default:                                                              break;
                        }

                        c[(y * nCh + ch) * nReads + k] = ref;
                    }
                }
            });
        });

        if(!known)
            throw std::invalid_argument("SCA::getReferencePixelCorrection: ramp is not an integer or floating point cube.");

        return c;
    }

    void SCA::referencePixelsCorrection(std::shared_ptr<DSL::FITScube>& ramp, const refPixCorrection& mode, const nisp::Geometry& geometry)
    {
        const std::valarray<double> c = getReferencePixelCorrection(ramp, mode, geometry);
        if(c.size() == 0)
            return;   // NONE, or the correction does not apply

        const size_t W        = geometry.width, H = geometry.height, b = geometry.refBorder;
        const size_t w        = geometry.colPerChannel, nCh = geometry.nChannels();
        const size_t nPix     = W * H;
        const size_t nLines   = H - 2 * b;
        const size_t nReads   = c.size() / (H * nCh);
        const size_t nThreads = std::thread::hardware_concurrency();

        withCubeData(*ramp, [&](auto& data)
        {
            using T = typename std::decay_t<decltype(data)>::value_type;

            // One task per read k and line y in [b, H − b); the other lines are not corrected.
            parallelFor(nLines * nReads, nThreads, [&](size_t begin, size_t end)
            {
                for(size_t j = begin; j < end; ++j)
                {
                    const size_t k    = j / nLines;
                    const size_t y    = b + j % nLines;
                    T*           line = &data[k * nPix + y * W];

                    for(size_t ch = 0; ch < nCh; ++ch)
                    {
                        // Unsigned cubes cannot hold v − c (≈ −baseline with ALL): they get
                        // v − (c_k − c_0), the change of the correction since the first read.
                        const double* cy  = &c[(y * nCh + ch) * nReads];
                        const double  ref = std::is_unsigned_v<T> ? cy[k] - cy[0] : cy[k];

                        for(size_t x = ch * w; x < ch * w + w; ++x)
                            line[x] = subtractReference(line[x], ref);
                    }
                }
            });
        });
    }
    

    std::shared_ptr<DSL::FITScube> SCA::expose_dpu(const std::valarray<double>& flux, const nisp::macc& maccParams)
    {
        const size_t nG = maccParams.nGroups;
        const size_t nR = maccParams.nFrames;
        const size_t nD = maccParams.nDrop;
        if(nG < 2)
            throw std::invalid_argument("SCA::expose_dpu: the DPU slope estimator needs at least 2 groups.");

        std::shared_ptr<DSL::FITScube> MACC = expose_macc(flux, maccParams);

        const std::valarray<uint16_t>* rampData = MACC->GetData<uint16_t>();
        if(!rampData)
            throw std::runtime_error("SCA::expose_dpu: the MACC ramp is not a uint16 cube.");

        // Reference pixel correction (Eq. 3) of every read, refCorr[(y * nCh + ch) * nG·nR + k],
        // with the output channels of the geometry.
        // It is subtracted in double below, while each pixel's reads are copied for the fit:
        // the uint16 ramp is not modified, nothing wraps, and the saturation test sees the
        // raw reads. Empty when "refPixCorr" is off, or for a matrix that is not the full
        // detector: no correction.
        const size_t columns = _geometry().colPerChannel, nCh = _geometry().nChannels();
        const std::valarray<double> refCorr = m_refPixCorr ? getReferencePixelCorrection(MACC, refPixCorrection::ALL, _geometry())
                                                           : std::valarray<double>();

        const size_t W    = MACC->Size(1);
        const size_t nPix = MACC->Size(1) * MACC->Size(2);

        std::shared_ptr<DSL::FITScube> result = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<uint16_t>({MACC->Size(1), MACC->Size(2),3}));
        *(result) *=0;
        std::valarray<uint16_t>& outData = *(result->GetData<uint16_t>());

        // nisp::Ramp works on doubles, in the units of the data: the read noise is
        // converted from e- to ADU like the ramp.
        const double sigmaADU = dpuNoise / dpuGain;

        setThreadCount(0);

        // One pixel per iteration, pixels split over the threads. Each thread owns two
        // double buffers reused for all its pixels (no allocation per pixel):
        //   reads  : the pixel's uint16 reads minus their reference pixel correction, in double
        //   groups : their group averages, filled by Ramp::coadd
        parallelFor(nPix, getThreadCount(), [&](size_t begin, size_t end)
        {
            nisp::rawLines reads(0.0, nG * nR);
            nisp::rawLines groups(0.0, nG);

            for(size_t i = begin; i < end; ++i)
            {
                // This pixel's corrections, one per read: line i / W, channel (i % W) / columns.
                const double* c = (refCorr.size() > 0) ? &refCorr[((i / W) * nCh + (i % W) / columns) * nG * nR] : nullptr;

                // Copy the reads. As on the flight DPU, a saturated read (65535) anywhere in the
                // ramp flags the pixel: no slope is fitted, the slope plane is 65535 (no scale,
                // no offset) and the other planes stay 0.
                bool saturated = false;
                for(size_t k = 0; k < nG * nR && !saturated; ++k)
                {
                    const uint16_t v = (*rampData)[k * nPix + i];
                    saturated = (v >= 65535);
                    reads[k] = static_cast<double>(v) - (c ? c[k] : 0.0);
                }

                if(saturated)
                {
                    outData[i] = 65535;
                    continue;
                }

                // The MACC cube holds the reads only (dropped frames are not stored), so
                // the reads are co-added with no drops; nD is then set for the estimator.
                const nisp::Ramp raw(&reads, std::slice(0, nG * nR, 1), sigmaADU, dpuGain);
                nisp::Ramp       grouped = raw.coadd(groups, nG, nR, 0);
                grouped.drops_per_group(static_cast<int16_t>(nD));

                double slope = 0., chi2 = 0., err = 0.;
                grouped.getSlope(slope, chi2, err);

                outData[i] = toUint16(slope * dpuScale + dpuOffset);   // rounded, clamped to [0, 65535], NaN → 0
                outData[i+nPix] = toUint16(chi2 * static_cast<double>(nG - 2));   // reduced chi2 → chi2 (2 groups: NaN → 0)
                outData[i+2*nPix] = toUint16(err * 100. * dpuScale + dpuOffset);            // x100 with same conversion as the slope
            }
        });

        // The slope is the signal per group interval, (nr + nd)·t_f; the integrated signal
        // it stands for spans the (ng − 1) intervals from the first to the last group
        // (87.2 s for MACC(4,16,4), 549.6 s for MACC(15,16,11)). INTTIME is that span, so
        // that slope × (ng − 1) / INTTIME is the signal rate in ADU/s.
        result->HDU().ValueForKey("INTTIME", static_cast<double>((maccParams.nGroups - 1) * (maccParams.nFrames + maccParams.nDrop)) * frameTime, "[s] Integration time, first to last group");
        result->HDU().ValueForKey("EXPTIME", static_cast<double>(maccParams.nGroups*(maccParams.nFrames+maccParams.nDrop)) * frameTime, "[s] Exposure time");
        result->HDU().ValueForKey("MODE", std::string("DPU"), "Simulation mode");
        result->HDU().ValueForKey("NG", maccParams.nGroups, "Number of groups");
        result->HDU().ValueForKey("NR", maccParams.nFrames, "Number of reads/group");
        result->HDU().ValueForKey("ND", maccParams.nDrop, "Number of drops/group");
        result->HDU().ValueForKey("REFCORR", std::string(refCorr.size() > 0 ? "ON" : "OFF"), "Reference pixel correction applied (Eq. 3, ALL)");
        if(m_enable_persistence)
            result->HDU().ValueForKey("PERSIST", std::string("ON"), "Persistence modelling enabled");

        _writeHeader(*result);

        return result;
    }

    void SCA::_simulate(const std::valarray<double>& flux, const double& exposureTime,
                        std::shared_ptr<DSL::FITScube>& frame, phot2pix::PixelConfig cfg)
    {
        // Full frame = the window covering the whole matrix.
        _simulate(flux, exposureTime, frame, 0, 0, getWidth(), getHeight(), cfg);
    }

    void SCA::_simulate(const std::valarray<double>& flux,
                        const double& exposureTime,
                        std::shared_ptr<DSL::FITScube>& frame,
                        size_t c0, size_t r0, size_t ncols, size_t nrows, phot2pix::PixelConfig cfg)
    {
        std::shared_ptr<DSL::FITScube>      persistenceSignal;

        // Here we draw the cosmic ray event. This gives us which pixel are affected and how much energy is deposited in each pixel.
        // Once done the cosmic ray signal is added to the measured signal before the persistence signal is added.
        // The draw runs in parallel with the simulation; get() waits for it just before the hits are added.
        auto cosmicFuture = std::async(std::launch::async,
            [this, exposureTime]() -> std::vector<phot2pix::Deposit>
            {
                std::vector<phot2pix::Deposit> hits;

                if(!m_cosmicDeposit.isConfigured())
                    return hits;

                m_cosmicDeposit.generate(exposureTime, _geometry(), hits);
                return hits;
            });

        // Extract window flux first — pure data copy, no compute cost.
        std::valarray<double> windowFlux(0.0, ncols * nrows);
        for (size_t r = 0; r < nrows; ++r)
            windowFlux[std::slice(r * ncols, ncols, 1)] = flux[std::slice((r0 + r) * getWidth() + c0, ncols, 1)];

        if (persistenceModelling && persistenceModelling->getHistory().size() > 0)
        {
            const size_t hw   = std::thread::hardware_concurrency();
            const size_t half = std::max(size_t(1), hw / 2);
            setThreadCount(half);
            persistenceModelling->setThreadCount((hw > half + 1) ? hw - half - 1 : 1);

            auto persFuture = std::async(std::launch::async,
                [this, exposureTime, c0, r0, ncols, nrows]() -> std::shared_ptr<DSL::FITScube>
                {
                    return persistenceModelling->computePersistenceSignal(
                        0.0, exposureTime, c0, r0, ncols, nrows);
                });

            frame = phot2pix::DetectorSimulator<nisp::pixel>::simulate(windowFlux, exposureTime, c0, r0, ncols, nrows, cfg);
            
            _addCosmicHits(*(frame->GetData<double>()), cosmicFuture.get(), c0, r0, ncols, nrows);
            
            persistenceSignal = persFuture.get();

            setThreadCount(0);
            persistenceModelling->setThreadCount(0);
        }
        else
        {
            frame = phot2pix::DetectorSimulator<nisp::pixel>::simulate(windowFlux, exposureTime, c0, r0, ncols, nrows, cfg);
            
            _addCosmicHits(*(frame->GetData<double>()), cosmicFuture.get(), c0, r0, ncols, nrows);
        }

        // Charge released by the past exposures, in ADU (gains set by _setPersistenceGain()).
        if (persistenceSignal && m_enable_persistence)
            *(frame->GetData<double>()) += *(persistenceSignal->GetData<double>());
    }

    void SCA::_setPersistenceGain()
    {
        if (!persistenceModelling)
            return;

        std::valarray<double> gain(size());
        for (size_t i = 0; i < size(); ++i)
            gain[i] = at(i).getGain();
        persistenceModelling->setGain(gain);
    }

    void SCA::_addCosmicHits(std::valarray<double>& frame, const std::vector<phot2pix::Deposit>& hits,
                             size_t c0, size_t r0, size_t ncols, size_t nrows)
    {
        if (hits.empty())
            return;

        // Charge of the hits kept for persistence (e-), only when persistence is modelled.
        const bool keep = persistenceModelling && m_enable_persistence;
        if (keep && m_cosmicCharge.size() != size())
            m_cosmicCharge.resize(size(), 0.0);

        const size_t W      = _geometry().width;   // hit.idx = Y * W + X on the physical detector
        const auto [x0, y0] = getOrigin();         // matrix placement on the detector

        for (const auto& hit : hits)
        {
            const size_t X = hit.idx % W, Y = hit.idx / W;
            if (X < x0 || X >= x0 + getWidth() || Y < y0 || Y >= y0 + getHeight())
                continue;                                            // outside this matrix

            const size_t col = X - x0, row = Y - y0;
            const size_t i   = row * getWidth() + col;               // matrix pixel, reference pixels included

            if (keep)
                m_cosmicCharge[i] += hit.electrons;                  // the whole matrix sees the hit

            if (col < c0 || col >= c0 + ncols || row < r0 || row >= r0 + nrows)
                continue;                                            // outside this window

            frame[(row - r0) * ncols + (col - c0)] += hit.electrons / m_pixels[i].getGain();   // e- → ADU
        }
    }

#pragma endregion
} // namespace nisp
