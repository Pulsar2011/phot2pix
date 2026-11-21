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

#include <algorithm>
#include <cctype>
#include <cmath>
#include <future>
#include <limits>
#include <optional>
#include <stdexcept>
#include <thread>
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

        double model::pixelresponse(const double& flux, const double& exposureTime, const phot2pix::PixelConfig& cfg) const
        {
            const double timeDelay = cfg.has("dT") ? cfg.getAsDouble("dT") : dT;
            double signal = 0.0;
            for (const auto& [idx, f] : _factor)
            {
                const double tau = _decayTime.at(idx);
                signal += f * std::exp(-timeDelay / tau) * tau * (std::exp(exposureTime / tau) - 1.0);
            }
            return signal * flux;
        }
#pragma endregion
// ── pixelArray ────────────────────────────────────────────────────────────────

#pragma region - pixelArray
        void pixelArray::incrementTime(const double& dT)
        {
            for (size_t i = 0; i < size(); ++i)
                at(i).incrementTime(dT);
        }

        void pixelArray::initialize(const phot2pix::PixelConfig& cfg)
        {
            if (!cfg.has("persistence"))
                return;

            const auto& persistJson = cfg.getJson("persistence");

            // Two supported formats:
            //   array:  [{"f": ..., "t": ...}, ...]
            //   object: {"traps": [...], "dT": <seconds, optional>}
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
        void scheduler::insertCharges(const std::shared_ptr<DSL::FITScube>& flux, const double& t0)
        {
            if (!m_pixelModel)
            {
                m_pixelModel = std::make_unique<pixelArray>(flux->Size(1), flux->Size(2));
                m_pixelModel->initialize(m_cfg);
                if (m_threadCount > 0)
                    m_pixelModel->setThreadCount(m_threadCount);
            }
            else if (flux->Size(1) != m_pixelModel->getWidth() ||
                     flux->Size(2) != m_pixelModel->getHeight())
            {
                throw std::invalid_argument(
                    "scheduler::insertCharges: flux dimensions do not match model.");
            }

            m_fluxHistory.push_back({t0, flux});
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
                phot2pix::PixelConfig cfg;
                cfg.set("dT", record.t0 + dT);

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

        // Run work(begin, end) over [0, n), split into contiguous chunks on nThreads
        // std::threads (so the compiler can vectorise each loop); a single call when the
        // job is too small to be worth threading.
        template <typename Work>
        void parallelFor(size_t n, size_t nThreads, const Work& work)
        {
            if (nThreads <= 1 || n < 100 * nThreads)
            {
                work(size_t(0), n);
                return;
            }

            const size_t chunk = (n + nThreads - 1) / nThreads;
            std::vector<std::thread> workers;
            workers.reserve(nThreads);
            for (size_t b = 0; b < n; b += chunk)
                workers.emplace_back(work, b, std::min(n, b + chunk));
            for (auto& w : workers)
                w.join();
        }
    }

    // One pass over memory: accumulate the charge and, for a read, write the measured
    // value straight into the output layer. A std::slice / std::gslice version would need
    // temporaries for the double → uint16 conversion and run on a single thread.
    // The readout noise map is in electrons: it is converted with each pixel's gain,
    // like the read noise in pixel::read().
    void SCA::_accumulate(const std::valarray<double>& charge, std::valarray<double>& cumul,
                          uint16_t* layer, const std::valarray<double>* oneOverF) const
    {
        parallelFor(cumul.size(), getThreadCount(),
            [this, &charge, &cumul, layer, oneOverF](size_t begin, size_t end)
            {
                if (!layer)                 // dropped frames: charge only
                    for (size_t i = begin; i < end; ++i)
                        cumul[i] += charge[i];
                else if (!oneOverF)
                    for (size_t i = begin; i < end; ++i)
                    {
                        cumul[i] += charge[i];
                        layer[i]  = toUint16(m_pixels[i].read(cumul[i]));
                    }
                else
                    for (size_t i = begin; i < end; ++i)
                    {
                        cumul[i] += charge[i];
                        layer[i]  = toUint16(m_pixels[i].read(cumul[i]) + (*oneOverF)[i] / m_pixels[i].getGain());   // e- → ADU
                    }
            });
    }

    void SCA::_read(std::valarray<double>& charge, double tStart) const
    {
        // Readout noise (1/f, bias voltage, pedestal) generated concurrently with the pixel read noise.
        std::future<std::valarray<double>> oneOverF;
        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
            oneOverF = std::async(std::launch::async,
                [this, tStart]() { return m_oneOverFNoise.generate(getWidth(), getHeight(), tStart); });

        parallelFor(size(), getThreadCount(),
            [this, &charge](size_t begin, size_t end)
            {
                for (size_t i = begin; i < end; ++i)
                    charge[i] = m_pixels[i].read(charge[i]);
            });

        if (oneOverF.valid())
        {
            const std::valarray<double> noiseMap = oneOverF.get();   // W × H, electrons
            parallelFor(size(), getThreadCount(),
                [this, &charge, &noiseMap](size_t begin, size_t end)
                {
                    for (size_t i = begin; i < end; ++i)
                        charge[i] += noiseMap[i] / m_pixels[i].getGain();   // e- → ADU
                });
        }
    }

    void SCA::_read(std::valarray<double>& charge, size_t c0, size_t r0, size_t ncols, size_t nrows, double tStart) const
    {
        // The windowed generate() still computes the full channel time series for each
        // overlapping channel (AR(1) correlation from frame start), but only returns
        // the ncols × nrows window portion.
        std::future<std::valarray<double>> oneOverF;
        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
            oneOverF = std::async(std::launch::async,
                [this, c0, r0, ncols, nrows, tStart]() { return m_oneOverFNoise.generate(getWidth(), getHeight(), c0, r0, ncols, nrows, tStart); });

        const size_t width = getWidth();
        parallelFor(nrows, getThreadCount(),
            [this, &charge, c0, r0, ncols, width](size_t rowBegin, size_t rowEnd)
            {
                for (size_t r = rowBegin; r < rowEnd; ++r)
                    for (size_t c = 0; c < ncols; ++c)
                        charge[r * ncols + c] = m_pixels[(r0 + r) * width + c0 + c].read(charge[r * ncols + c]);
            });

        if (oneOverF.valid())
        {
            const std::valarray<double> noiseMap = oneOverF.get();   // ncols × nrows, electrons
            parallelFor(nrows, getThreadCount(),
                [this, &charge, &noiseMap, c0, r0, ncols, width](size_t rowBegin, size_t rowEnd)
                {
                    for (size_t r = rowBegin; r < rowEnd; ++r)
                        for (size_t c = 0; c < ncols; ++c)
                            charge[r * ncols + c] += noiseMap[r * ncols + c] / m_pixels[(r0 + r) * width + c0 + c].getGain();   // e- → ADU
                });
        }
    }

    void SCA::initialize(const phot2pix::PixelConfig& cfg)
    {
        phot2pix::DetectorSimulator<nisp::pixel>::initialize(cfg);

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
            const std::string modeName = cfg.get<std::string>("exposureMode");
            std::string key = modeName;
            std::transform(key.begin(), key.end(), key.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

            if      (key == "MACC") m_mode = nisp::mode::MACC;
            else if (key == "UTR")  m_mode = nisp::mode::UTR;
            else if (key == "CDS")  m_mode = nisp::mode::CDS;
            else if (key == "DPU")  m_mode = nisp::mode::DPU;
            else if (key == "FLUX") m_mode = nisp::mode::FLUX;
            else
                throw std::invalid_argument("SCA::initialize: unknown exposureMode '" + modeName +
                                            "' (expected MACC, UTR, CDS, DPU or FLUX).");
        }

        if(cfg.has("persistence"))
        {
            m_enable_persistence |= true;
            if (!persistenceModelling)
                persistenceModelling = std::make_shared<persistence::scheduler>(cfg);
        }

        // Readout noise added at each read: per-channel 1/f ("oneOverF"), 1/f common to
        // every channel ("bias_voltage_noise") and channel offsets ("pedestal_drift").
        if(cfg.has("oneOverF") || cfg.has("bias_voltage_noise") || cfg.has("pedestal_drift"))
        {
            m_enable_oneOverF |= true;
            m_oneOverFNoise.initialize(cfg);
        }
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

        _simulate(fullFlux, exposureTime, frame, c0, r0, ncols, nrows);
        // Single read: read noise + readout noise, the read being the last frame of the exposure.
        _read(*(frame->GetData<double>()), c0, r0, ncols, nrows, std::max(0.0, exposureTime - frameTime));

        // Register full-frame charge (electrons) for future persistence computation.
        if (persistenceModelling && m_enable_persistence)
        {
            auto chargeImg = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<double>({getWidth(), getHeight()}));

            std::valarray<double> chargeMap(0.0, size());
            for (size_t i = 0; i < size(); ++i)
                chargeMap[i] = fullFlux[i] * at(i).getQuantumEfficiency() * exposureTime;
            
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        return _toUint16(frame);
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
        std::shared_ptr<DSL::FITScube> frame    = nullptr;
        std::shared_ptr<DSL::FITScube> skipped  = nullptr;
        size_t                         layer    = 0;

        // One continuous 1/f noise stream for the whole ramp: every frame (read or
        // dropped) consumes the next chunk of a single per-channel AR(1) process
        // instead of restarting the amplifier state at zero for each read.
        std::optional<nisp::OneOverFNoise::Stream> oneOverFStream;
        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
            oneOverFStream = m_oneOverFNoise.beginStream(getWidth(), getHeight());

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

                // cumul += this frame's charge, then read it: layer = pixel read (read noise) + 1/f.
                _accumulate(*(frame->GetData<double>()), cumul, &(*outData)[layer * nPix],
                            oneOverF.size() ? &oneOverF : nullptr);
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

                _accumulate(*(skipped->GetData<double>()), cumul, nullptr, nullptr);   // dropped frames: charge only, no layer
                skipped.reset();

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
            
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        std::shared_ptr<DSL::FITScube> result = out;
        result->HDU().ValueForKey("INTTIME", static_cast<double>(nG*(nF+nD)-nD)*frameTime, "Integration time in seconds");
        result->HDU().ValueForKey("EXPTIME", exposureTime, "Exposure time in seconds");
        result->HDU().ValueForKey("MODE", std::string("MACC"), "Simulation mode");
        result->HDU().ValueForKey("NG", nG, "Number of groups");
        result->HDU().ValueForKey("NR", nF, "Number of reads/group");
        result->HDU().ValueForKey("ND", nD, "Number of drops/group");
        if(m_enable_persistence)
            result->HDU().ValueForKey("PERSIST", std::string("ON"), "Persistence modelling enabled");

        result->HDU().ValueForKey("OFFSET", dpuOffset, "DPU offset in ADU");
        result->HDU().ValueForKey("SFACTOR", dpuScale, "DPU scale factor");
        result->HDU().ValueForKey("GAIN", dpuGain, "Mean DPU gain");
        result->HDU().ValueForKey("RN", dpuNoise, "Mean DPU read noise");

        result->HDU().ValueForKey("FRTIME", frameTime, "NISP HAWAII-2RG frame time in seconds");
        result->HDU().ValueForKey("PHOT2PIX", std::string("v")+std::to_string(PHOT2PIX_VERSION_MAJOR)
                                             +std::string(".")+std::to_string(PHOT2PIX_VERSION_MINOR)
                                             +std::string(".")+std::to_string(PHOT2PIX_VERSION_PATCH)
                                             +std::string("-")+PHOT2PIX_VERSION, DSL::key_type::fChar,  "phot2pix version");
        

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
        _read(read1, 0.0);
        _read(readN, cdsTime);

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
            
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        // DPU scale / offset applied like the FLUX mode, so negative CDS values (noise on
        // dark pixels) survive when an offset is set, instead of being clipped to 0.
        std::shared_ptr<DSL::FITScube> result = _toUint16(f1);
        result->HDU().ValueForKey("INTTIME", cdsTime, "Integration time in seconds");   // read 1 → read N
        result->HDU().ValueForKey("EXPTIME", eff_exposureTime, "Exposure time in seconds");
        result->HDU().ValueForKey("MODE", std::string("CDS"), "Simulation mode");
        if(m_enable_persistence)
            result->HDU().ValueForKey("PERSIST", std::string("ON"), "Persistence modelling enabled");

        result->HDU().ValueForKey("OFFSET", dpuOffset, "DPU offset in ADU");
        result->HDU().ValueForKey("SFACTOR", dpuScale, "DPU scale factor");
        result->HDU().ValueForKey("GAIN", dpuGain, "Mean DPU gain");
        result->HDU().ValueForKey("RN", dpuNoise, "Mean DPU read noise");

        result->HDU().ValueForKey("FRTIME", frameTime, "NISP HAWAII-2RG frame time in seconds");
        result->HDU().ValueForKey("POT2PIX",  std::string("v")+std::to_string(PHOT2PIX_VERSION_MAJOR)
                                             +std::string(".")+std::to_string(PHOT2PIX_VERSION_MINOR)
                                             +std::string(".")+std::to_string(PHOT2PIX_VERSION_PATCH)
                                             +std::string("-")+PHOT2PIX_VERSION, DSL::key_type::fChar, "phot2pix version");

        f1.reset();
        return result;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose_flux(const std::valarray<double>& flux, const double& exposureTime)
    {
        std::shared_ptr<DSL::FITScube> frame;

        double effectiveExposureTime = std::floor(exposureTime/frameTime)*frameTime;

        _simulate(flux, effectiveExposureTime, frame);
        // Single read: read noise + readout noise, the read being the last frame of the exposure.
        _read(*(frame->GetData<double>()), std::max(0.0, effectiveExposureTime - frameTime));

        // Register accumulated charge (electrons) for future persistence computation.
        if (persistenceModelling && m_enable_persistence)
        {
            auto chargeImg = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<double>({getWidth(), getHeight()}));
            std::valarray<double> chargeMap(0.0, size());
            for (size_t i = 0; i < size(); ++i)
                chargeMap[i] = flux[i] * at(i).getQuantumEfficiency() * effectiveExposureTime;
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        std::shared_ptr<DSL::FITScube> result = _toUint16(frame);
        result->HDU().ValueForKey("INTTIME", effectiveExposureTime, "Integration time in seconds");
        result->HDU().ValueForKey("EXPTIME", effectiveExposureTime, "Exposure time in seconds");
        result->HDU().ValueForKey("MODE", std::string("FLUX"), "Simulation mode");
        if(m_enable_persistence)
            result->HDU().ValueForKey("PERSIST", std::string("ON"), "Persistence modelling enabled");

        result->HDU().ValueForKey("OFFSET", dpuOffset, "DPU offset in ADU");
        result->HDU().ValueForKey("SFACTOR", dpuScale, "DPU scale factor");
        result->HDU().ValueForKey("GAIN", dpuGain, "Mean DPU gain");
        result->HDU().ValueForKey("RN", dpuNoise, "Mean DPU read noise");

        result->HDU().ValueForKey("FRTIME", frameTime, "NISP HAWAII-2RG frame time in seconds");
        result->HDU().ValueForKey("POT2PIX",  std::string("v")+std::to_string(PHOT2PIX_VERSION_MAJOR)
                                             +std::string(".")+std::to_string(PHOT2PIX_VERSION_MINOR)
                                             +std::string(".")+std::to_string(PHOT2PIX_VERSION_PATCH)
                                             +std::string("-")+PHOT2PIX_VERSION, DSL::key_type::fChar, "phot2pix version");

        frame.reset();
        return result;
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

        const size_t nPix = MACC->Size(1) * MACC->Size(2);

        std::shared_ptr<DSL::FITScube> result = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<uint16_t>({MACC->Size(1), MACC->Size(2)}));
        *(result) *=0;
        std::valarray<uint16_t>& outData = *(result->GetData<uint16_t>());

        // nisp::Ramp works on doubles, in the units of the data: the read noise is
        // converted from e- to ADU like the ramp.
        const double sigmaADU = dpuNoise / dpuGain;

        setThreadCount(0);

        // One pixel per iteration, pixels split over the threads. Each thread owns two
        // double buffers reused for all its pixels (no allocation per pixel):
        //   reads  : the pixel's uint16 reads, converted to double for nisp::Ramp
        //   groups : their group averages, filled by Ramp::coadd
        parallelFor(nPix, getThreadCount(), [&](size_t begin, size_t end)
        {
            nisp::rawLines reads(0.0, nG * nR);
            nisp::rawLines groups(0.0, nG);

            for(size_t i = begin; i < end; ++i)
            {
                // Copy the reads group by group, up to the first saturated group: a
                // saturated read makes the rest of the ramp unusable, only the groups
                // before it are fitted.
                size_t nGood = 0;
                for(size_t g = 0; g < nG; ++g)
                {
                    bool saturated = false;
                    for(size_t r = 0; r < nR; ++r)
                    {
                        const uint16_t v = (*rampData)[(g * nR + r) * nPix + i];
                        saturated |= (v >= 65535);
                        reads[g * nR + r] = static_cast<double>(v);
                    }
                    if(saturated)
                        break;
                    ++nGood;
                }

                if(nGood < 2)
                {
                    outData[i] = 65535;   // saturated before two groups: no slope
                    continue;
                }

                // The MACC cube holds the reads only (dropped frames are not stored), so
                // the reads are co-added with no drops; nD is then set for the estimator.
                const nisp::Ramp raw(&reads, std::slice(0, nGood * nR, 1), sigmaADU, dpuGain);
                nisp::Ramp       grouped = raw.coadd(groups, nGood, nR, 0);
                grouped.drops_per_group(static_cast<int16_t>(nD));

                double slope = 0., chi2 = 0., err = 0.;
                grouped.getSlope(slope, chi2, err);

                outData[i] = toUint16(slope * dpuScale + dpuOffset);   // rounded, clamped to [0, 65535], NaN → 0
            }
        });

        result->HDU().ValueForKey("INTTIME", static_cast<double>(maccParams.nGroups*(maccParams.nFrames+maccParams.nDrop) - maccParams.nDrop) * frameTime, "Integration time in seconds");
        result->HDU().ValueForKey("EXPTIME", static_cast<double>(maccParams.nGroups*(maccParams.nFrames+maccParams.nDrop)) * frameTime, "Exposure time in seconds");
        result->HDU().ValueForKey("MODE", std::string("DPU"), "Simulation mode");
        result->HDU().ValueForKey("NG", maccParams.nGroups, "Number of groups");
        result->HDU().ValueForKey("NR", maccParams.nFrames, "Number of reads/group");
        result->HDU().ValueForKey("ND", maccParams.nDrop, "Number of drops/group");
        if(m_enable_persistence)
            result->HDU().ValueForKey("PERSIST", std::string("ON"), "Persistence modelling enabled");

        result->HDU().ValueForKey("OFFSET", dpuOffset, "DPU offset in ADU");
        result->HDU().ValueForKey("SFACTOR", dpuScale, "DPU scale factor");
        result->HDU().ValueForKey("GAIN", dpuGain, "Mean DPU gain");
        result->HDU().ValueForKey("RN", dpuNoise, "Mean DPU read noise");

        result->HDU().ValueForKey("FRTIME", frameTime, "NISP HAWAII-2RG frame time in seconds");
        result->HDU().ValueForKey("POT2PIX",  std::string("v")+std::to_string(PHOT2PIX_VERSION_MAJOR)
                                             +std::string(".")+std::to_string(PHOT2PIX_VERSION_MINOR)
                                             +std::string(".")+std::to_string(PHOT2PIX_VERSION_PATCH)
                                             +std::string("-")+PHOT2PIX_VERSION, DSL::key_type::fChar, "phot2pix version");

        return result;
    }



    void SCA::_simulate(const std::valarray<double>& flux, const double& exposureTime,
                        std::shared_ptr<DSL::FITScube>& frame, phot2pix::PixelConfig cfg)
    {
        std::shared_ptr<DSL::FITScube>       persistenceSignal;

        if (persistenceModelling && persistenceModelling->getHistory().size() > 0)
        {
            const size_t hw   = std::thread::hardware_concurrency();
            const size_t half = std::max(size_t(1), hw / 2);
            setThreadCount(half);
            persistenceModelling->setThreadCount(hw - half);

            auto persFuture = std::async(std::launch::async,
                [this, exposureTime]() -> std::shared_ptr<DSL::FITScube>
                {
                    return persistenceModelling->computePersistenceSignal(0.0, exposureTime);
                });

            frame = phot2pix::DetectorSimulator<nisp::pixel>::simulate(flux, exposureTime, cfg);
            persistenceSignal = persFuture.get();

            setThreadCount(0);
            persistenceModelling->setThreadCount(0);
        }
        else
        {
            frame = phot2pix::DetectorSimulator<nisp::pixel>::simulate(flux, exposureTime, cfg);
        }

        if (persistenceSignal && m_enable_persistence)
        {
            *(frame->GetData<double>()) += *(persistenceSignal->GetData<double>());
            persistenceSignal.reset();
        }
    }

    void SCA::_simulate(const std::valarray<double>& flux,
                        const double& exposureTime,
                        std::shared_ptr<DSL::FITScube>& frame,
                        size_t c0, size_t r0, size_t ncols, size_t nrows, phot2pix::PixelConfig cfg)
    {
        std::shared_ptr<DSL::FITScube>      persistenceSignal;

        // Extract window flux first — pure data copy, no compute cost.
        std::valarray<double> windowFlux(0.0, ncols * nrows);
        for (size_t r = 0; r < nrows; ++r)
            windowFlux[std::slice(r * ncols, ncols, 1)] = flux[std::slice((r0 + r) * getWidth() + c0, ncols, 1)];

        if (persistenceModelling && persistenceModelling->getHistory().size() > 0)
        {
            const size_t hw   = std::thread::hardware_concurrency();
            const size_t half = std::max(size_t(1), hw / 2);
            setThreadCount(half);
            persistenceModelling->setThreadCount(hw - half);

            auto persFuture = std::async(std::launch::async,
                [this, exposureTime, c0, r0, ncols, nrows]() -> std::shared_ptr<DSL::FITScube>
                {
                    return persistenceModelling->computePersistenceSignal(
                        0.0, exposureTime, c0, r0, ncols, nrows);
                });

            frame = phot2pix::DetectorSimulator<nisp::pixel>::simulate(windowFlux, exposureTime, c0, r0, ncols, nrows, cfg);
            persistenceSignal = persFuture.get();

            setThreadCount(0);
            persistenceModelling->setThreadCount(0);
        }
        else
        {
            frame = phot2pix::DetectorSimulator<nisp::pixel>::simulate(windowFlux, exposureTime, c0, r0, ncols, nrows, cfg);
        }

        if (persistenceSignal && m_enable_persistence)
            *(frame->GetData<double>()) += *(persistenceSignal->GetData<double>());
    }

#pragma endregion
} // namespace nisp
