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
#include <limits>
#include <stdexcept>
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

    std::shared_ptr<DSL::FITScube> SCA::expose( const std::valarray<double>& flux, double exposureTime, double dT)
    {
        // 1. Persistence contribution from past exposures (electrons)
        std::shared_ptr<DSL::FITScube> persistenceSignal;
        if (persistenceModelling)
            persistenceSignal = persistenceModelling->computePersistenceSignal(dT, exposureTime);

        // 2. Pixel response: noisy electrons in double precision
        std::shared_ptr<DSL::FITScube> frame =
            phot2pix::DetectorSimulator<nisp::pixel>::simulate(flux, exposureTime);

        // 3. Add persistence electrons before ADU conversion
        if (persistenceSignal)
            *(frame->GetData<double>()) += *(persistenceSignal->GetData<double>());

        // 4. Register accumulated charge (electrons) for future persistence computation
        if (persistenceModelling)
        {
            auto chargeImg = std::shared_ptr<DSL::FITScube>( new DSL::FITSimg<double>({getWidth(), getHeight()}) );
            std::valarray<double> chargeMap(0.0, size());
            for (size_t i = 0; i < size(); ++i)
                chargeMap[i] = flux[i] * at(i).getQuantumEfficiency() * exposureTime;
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        std::shared_ptr<DSL::FITScube> result = _toUint16(frame);
        
        frame.reset(); // free memory
        return result;
    }

    std::shared_ptr<DSL::FITScube> SCA::expose( const std::valarray<double>& fullFlux, double exposureTime, double dT,
                                                size_t c0, size_t r0, size_t ncols, size_t nrows)
    {
        // 1. Windowed persistence signal
        std::shared_ptr<DSL::FITScube> persistenceSignal;
        if (persistenceModelling)
            persistenceSignal = persistenceModelling->computePersistenceSignal(dT, exposureTime, c0, r0, ncols, nrows);

        // 2. Extract window flux and simulate pixel response (electrons)
        std::valarray<double> windowFlux(0.0, ncols * nrows);
        for (size_t r = 0; r < nrows; ++r)
            windowFlux[std::slice(r * ncols, ncols, 1)] =
                fullFlux[std::slice((r0 + r) * getWidth() + c0, ncols, 1)];

        std::shared_ptr<DSL::FITScube> frame = phot2pix::DetectorSimulator<nisp::pixel>::simulate(windowFlux, exposureTime, c0, r0, ncols, nrows);

        // 3. Add windowed persistence electrons before ADU conversion
        if (persistenceSignal)
            *(frame->GetData<double>()) += *(persistenceSignal->GetData<double>());

        // 4. Register full-frame charge (electrons) for future persistence computation
        if (persistenceModelling)
        {
            auto chargeImg = std::shared_ptr<DSL::FITScube>( new DSL::FITSimg<double>({getWidth(), getHeight()}));
            std::valarray<double> chargeMap(0.0, size());
            for (size_t i = 0; i < size(); ++i)
                chargeMap[i] = fullFlux[i] * at(i).getQuantumEfficiency() * exposureTime;
            *(chargeImg) += chargeMap;
            persistenceModelling->insertCharges(chargeImg, 0.0);
        }

        // 5. Apply DPU scale/offset and saturate to uint16
        return _toUint16(frame);
    }

#pragma endregion
} // namespace nisp
