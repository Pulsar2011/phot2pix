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
#include <cmath>
#include <future>
#include <limits>
#include <optional>
#include <stdexcept>
#include <thread>
#include <valarray>

namespace nisp
{
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

    void SCA::initialize(const phot2pix::PixelConfig& cfg)
    {
        phot2pix::DetectorSimulator<nisp::pixel>::initialize(cfg);

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

        if(cfg.has("persistence"))
        {
            m_enable_persistence |= true;
            if (!persistenceModelling)
                persistenceModelling = std::make_shared<persistence::scheduler>(cfg);
        }

        if(cfg.has("oneOverF"))
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
        size_t nG = maccParams.nGroups;
        size_t nF = maccParams.nFrames;
        size_t nD = maccParams.nDrop;

        std::shared_ptr<DSL::FITScube> ramp        = nullptr;
        std::shared_ptr<DSL::FITScube> skipped     = nullptr;
        std::shared_ptr<DSL::FITScube> frame       = nullptr;
        std::shared_ptr<DSL::FITScube> cummul      = nullptr;

        double exposureTime = 0;
        double integrationTime = frameTime;
        resetSimulationState();

        phot2pix::PixelConfig cfg;
        cfg.set("baseline", std::string("off"));

        // One continuous 1/f noise stream for the whole ramp: every read
        // (kept or dropped) consumes the next chunk of a single per-channel
        // AR(1) process instead of each _simulate() call independently
        // restarting the amplifier state at zero.
        std::optional<nisp::OneOverFNoise::Stream> oneOverFStream;
        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
            oneOverFStream = m_oneOverFNoise.beginStream(getWidth(), getHeight());

        for (size_t g = 0; g < nG; ++g)
        {
            for(size_t f = 0; f < nF; ++f)
            {
                integrationTime = frameTime;
                exposureTime   += integrationTime;

                _simulate(fullFlux, integrationTime, frame, ((g < 1)? phot2pix::PixelConfig() : cfg),
                          oneOverFStream ? &(*oneOverFStream) : nullptr);
                if(!frame)
                    throw std::runtime_error("SCA::expose: failed to simulate frame ["+std::to_string(g)+","+std::to_string(f)+"]");

                if(!ramp)
                {
                    ramp = frame;
                    cummul = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<double>({frame->Size(1), frame->Size(2)}));
                    *(cummul) += *(frame->GetData<double>());
                }
                else
                {
                    *(frame->GetData<double>()) += *(cummul->GetData<double>());
                    ramp->AddLayer(*(frame));

                    std::swap(cummul,frame);
                    frame.reset();
                }
            }

            if(nD <= 0)
                continue;

            integrationTime = frameTime * nD;
            exposureTime += integrationTime;

            if(g+1 < nG)
            {
                _simulate(fullFlux, integrationTime, skipped, cfg,
                          oneOverFStream ? &(*oneOverFStream) : nullptr);
                if(!skipped)
                    throw std::runtime_error("SCA::expose: failed to simulate skipped frame ["+std::to_string(g)+","+std::to_string(nF-1)+",0]");
                
                if(!cummul)
                {
                    cummul = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<double>({getWidth(), getHeight()}));
                    *(cummul->GetData<double>()) *=0.;
                }
                
                *(cummul) += *(skipped->GetData<double>());
                skipped.reset();
            }
        }

        cummul.reset();
        skipped.reset();
        frame.reset();

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

        std::shared_ptr<DSL::FITScube> result = ramp->cast<uint16_t>();
        result->HDU().ValueForKey("INTTIME", integrationTime, "Integration time in seconds");
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
        result->HDU().ValueForKey("POT2PIX", PHOT2PIX_VERSION, "phot2pix version");
        

        ramp.reset();
        return result;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose_cds(const std::valarray<double>& fullFlux, const double& exposureTime)
    {
        if(exposureTime < 2*frameTime)
            throw std::invalid_argument("SCA::expose_cds: exposure time must be at least 2*"+std::to_string(frameTime)+" seconds for CDS mode.");

        resetSimulationState();

        double nFrames = std::floor(exposureTime / frameTime);
        double eff_exposureTime = nFrames * frameTime;
        
        std::shared_ptr<DSL::FITScube> f0 = nullptr;
        std::shared_ptr<DSL::FITScube> f1 = nullptr;
        
        double this_exposureTime = frameTime;

        
        _simulate(fullFlux, this_exposureTime, f0);
        if(!f0)
            throw std::runtime_error("SCA::expose_cds: failed to simulate frame [f0]");

        this_exposureTime = eff_exposureTime - frameTime;

        phot2pix::PixelConfig cfg;
        cfg.set("baseline", std::string("off"));
        
        _simulate(fullFlux, this_exposureTime, f1, cfg);
        if(!f1)
            throw std::runtime_error("SCA::expose_cds: failed to simulate frame [f1]");
        
        *(f0->GetData<double>()) -= *(f1->GetData<double>());

        f1.reset();

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

        std::shared_ptr<DSL::FITScube> result = f0->cast<uint16_t>();
        result->HDU().ValueForKey("INTTIME", eff_exposureTime, "Integration time in seconds");
        result->HDU().ValueForKey("EXPTIME", eff_exposureTime, "Exposure time in seconds");
        result->HDU().ValueForKey("MODE", std::string("CDS"), "Simulation mode");
        if(m_enable_persistence)
            result->HDU().ValueForKey("PERSIST", std::string("ON"), "Persistence modelling enabled");

        result->HDU().ValueForKey("OFFSET", dpuOffset, "DPU offset in ADU");
        result->HDU().ValueForKey("SFACTOR", dpuScale, "DPU scale factor");
        result->HDU().ValueForKey("GAIN", dpuGain, "Mean DPU gain");
        result->HDU().ValueForKey("RN", dpuNoise, "Mean DPU read noise");

        result->HDU().ValueForKey("FRTIME", frameTime, "NISP HAWAII-2RG frame time in seconds");
        result->HDU().ValueForKey("POT2PIX", PHOT2PIX_VERSION, "phot2pix version");

        f0.reset();
        return result;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose_flux(const std::valarray<double>& flux, const double& exposureTime)
    {
        std::shared_ptr<DSL::FITScube> frame;

        double effectiveExposureTime = std::floor(exposureTime/frameTime)*frameTime;
        
        resetSimulationState();
        _simulate(flux, effectiveExposureTime, frame);

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
        result->HDU().ValueForKey("POT2PIX", PHOT2PIX_VERSION, "phot2pix version");

        frame.reset();
        return result;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose_dpu(const std::valarray<double>& flux, const nisp::macc& maccParams)
    {
        std::shared_ptr<DSL::FITScube> MACC = expose_macc(flux, maccParams);
        if(MACC->Size(3) < 2)
            throw std::runtime_error("SCA::expose_dpu: insufficient number of groups for DPU flux estimation.");

        std::shared_ptr<DSL::FITScube> result = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<uint16_t>({MACC->Size(1), MACC->Size(2)}));
        *(result) *=0;

        setThreadCount(0);

        for(size_t i = 0; i < MACC->Size(1) * MACC->Size(2); ++i)
        {
            double slope=0.;
            double chi2 =0.;
            double err  =0.;

            if((*(MACC->GetData<double>()))[i*MACC->Size(1)*MACC->Size(2)] >= 65535) // Reset the first group to zero for slope estimation
            {
                (*(result->GetData<uint16_t>()))[i] = 65535;
                continue;
            }

            nisp::Ramp dpu = nisp::Ramp(MACC->GetData<double>(), std::slice(i, MACC->Size(3), MACC->Size(1)*MACC->Size(2)),dpuNoise, dpuGain);
            dpu.getSlope(slope, chi2, err);

            uint16_t dpuValue = (slope >= 65535) ? 65535 : ((slope <= 0) ? 0 : static_cast<uint16_t>(slope * dpuScale + dpuOffset +0.5));
            (*(result->GetData<uint16_t>()))[i] = dpuValue;
        }

        result->HDU().ValueForKey("INTTIME", maccParams.nFrames * frameTime, "Integration time in seconds");
        result->HDU().ValueForKey("EXPTIME", maccParams.nFrames * frameTime, "Exposure time in seconds");
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
        result->HDU().ValueForKey("POT2PIX", PHOT2PIX_VERSION, "phot2pix version");

        return result;
    }



    void SCA::_simulate(const std::valarray<double>& flux, const double& exposureTime,
                        std::shared_ptr<DSL::FITScube>& frame, phot2pix::PixelConfig cfg,
                        nisp::OneOverFNoise::Stream* oneOverFStream)
    {
        std::shared_ptr<DSL::FITScube>       persistenceSignal;
        std::future<std::valarray<double>>   oneOverFFuture;

        // Launch 1/f noise generation before pixel simulation so both run concurrently.
        // generate()/Stream::next() are themselves parallel over H2RG output channels.
        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
        {
            if (oneOverFStream)
                oneOverFFuture = std::async(std::launch::async,
                    [oneOverFStream]() { return oneOverFStream->next(); });
            else
                oneOverFFuture = std::async(std::launch::async,
                    [this]() { return m_oneOverFNoise.generate(getWidth(), getHeight()); });
        }

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

        if (oneOverFFuture.valid())
            *(frame->GetData<double>()) += oneOverFFuture.get();
    }

    void SCA::_simulate(const std::valarray<double>& flux,
                        const double& exposureTime,
                        std::shared_ptr<DSL::FITScube>& frame,
                        size_t c0, size_t r0, size_t ncols, size_t nrows, phot2pix::PixelConfig cfg)
    {
        std::shared_ptr<DSL::FITScube>      persistenceSignal;
        std::future<std::valarray<double>>  oneOverFFuture;

        // Extract window flux first — pure data copy, no compute cost.
        std::valarray<double> windowFlux(0.0, ncols * nrows);
        for (size_t r = 0; r < nrows; ++r)
            windowFlux[std::slice(r * ncols, ncols, 1)] = flux[std::slice((r0 + r) * getWidth() + c0, ncols, 1)];

        // Launch 1/f noise concurrently with pixel simulation.
        // The windowed generate() still computes the full channel time series for each
        // overlapping channel (preserving AR(1) correlation from frame start), but only
        // returns the ncols × nrows window portion.
        if (m_enable_oneOverF && m_oneOverFNoise.isEnabled())
            oneOverFFuture = std::async(std::launch::async,
                [this, c0, r0, ncols, nrows]() {
                    return m_oneOverFNoise.generate(getWidth(), getHeight(), c0, r0, ncols, nrows);
                });

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

        // frame->GetData() has W×H elements; window data occupies the first ncols*nrows.
        // Use a slice to add only to the window portion.
        if (oneOverFFuture.valid())
        {
            auto noiseMap = oneOverFFuture.get();   // ncols × nrows
            (*frame->GetData<double>())[std::slice(0, noiseMap.size(), 1)] += noiseMap;
        }
    }

    void SCA::resetSimulationState()
    {
        for (size_t i = 0; i < size(); ++i)
            phot2pix::DetectorSimulator<nisp::pixel>::at(i).resetState();
    }

#pragma endregion
} // namespace nisp
