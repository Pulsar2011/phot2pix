//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispSCA.h
//

#ifndef _NISP_SCA_H_
#define _NISP_SCA_H_

#include "nispPixel.h"
#include "nispOneOverFNoise.h"
#include <vector>
#include <map>
#include <memory>
#include <algorithm>

#include <phot2pix/phot2pix.h>

namespace nisp
{
    enum class mode
    {
        CDS,   ///< Correlated Double Sampling
        MACC,  ///< Multiple Accumulated Correlated Double Sampling
        UTR,   ///< Sample Up The Ramp
        DPU,   ///< Flux estimation emulating the onboard DPU algorithm
        FLUX   ///< Fast flux estimation (last frame only)
    };

    struct macc
    {
        size_t nGroups;
        size_t nFrames;  ///< Number of frames to accumulate
        size_t nDrop;    ///< Number of frames to skip between accumulations

        macc UTR(const size_t& nFrames) { return {1, nFrames, 0}; }
        macc MACC(const size_t& nGroups, const size_t& nFrames, const size_t& nDrop) { return {nGroups, nFrames, nDrop}; }
        macc CDS(const size_t nDrops) { return {2, 1, nDrops}; }
    };

// ── Persistence model ─────────────────────────────────────────────────────────
#pragma region - persistence model
    namespace persistence
    {

#pragma region - pixel model
        class model : private phot2pix::pixsim
        {
            protected:
                typedef std::map<size_t, double> parameter;
                double dT = 0.0; ///< Elapsed time since last exposure (seconds), used as fallback when model is used outside of a scheduler context (e.g. standalone pixelArray).

            private:
                parameter _factor;    ///< index → trap fraction
                parameter _decayTime; ///< index → trap decay time (seconds)

            public:
                model() : phot2pix::pixsim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0) {}
                model(const parameter& factor, const parameter& decayTime): phot2pix::pixsim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0),_factor(factor), _decayTime(decayTime) {}
                model(const double& capturedFraction, const parameter& factor, const parameter& decayTime): phot2pix::pixsim(0.0, 0.0, 1.0, capturedFraction, 0.0, 0.0),_factor(factor), _decayTime(decayTime) {}
                model(const model& other): phot2pix::pixsim(other), dT(other.dT),_factor(other._factor), _decayTime(other._decayTime) {}

                virtual ~model() = default;

                void insertParameter(const double& factor, const double& decayTime);
                void clearTraps() { _factor.clear(); _decayTime.clear(); }

                inline const parameter& getFactor()    const { return _factor; }
                inline const parameter& getDecayTime() const { return _decayTime; }

                inline const double& getCapturedFraction() const { return getQuantumEfficiency(); }
                inline void setCapturedFraction(const double& fR) { setQuantumEfficiency(fR); }

                inline double getTimeSinceLastExposure() const { return dT; }
                void incrementTime(const double& bydT) { dT += bydT; }

                // Standard path: uses this->dT as the elapsed time.
                double pixelresponse(const double& flux, const double& exposureTime) const override;

                // cfg-aware path: reads "dT" from cfg if present, falls back to this->dT.
                // The scheduler always provides cfg with "dT" set; the fallback supports
                // standalone pixelArray use without a scheduler.
                double pixelresponse(const double& flux, const double& exposureTime, const phot2pix::PixelConfig& cfg) const override;
        };
#pragma endregion

#pragma region - pixel array
        class pixelArray : public phot2pix::DetectorSimulator<nisp::persistence::model>
        {
            public:
                pixelArray() : phot2pix::DetectorSimulator<nisp::persistence::model>(2048, 2048) {}
                pixelArray(size_t w, size_t h)
                    : phot2pix::DetectorSimulator<nisp::persistence::model>(w, h) {}
                virtual ~pixelArray() = default;

                void incrementTime(const double& dT);
                void initialize(const phot2pix::PixelConfig& cfg) override;

            protected:
                void _engine(const std::valarray<double>& flux,
                             double dt, size_t begin, size_t end,
                             std::valarray<double>& frame,
                             const phot2pix::PixelConfig& cfg) const override;

                void _engine(const std::valarray<double>& flux,
                             double dt, size_t c0, size_t r0,
                             size_t row_begin, size_t row_end, size_t ncols,
                             std::valarray<double>& frame,
                             const phot2pix::PixelConfig& cfg) const override;
        };
#pragma endregion

#pragma region - scheduler
        // Scheduler holds one shared pixelArray (trap parameters) plus the flux
        // history of past exposures.  Memory cost: O(N * W*H * sizeof(double))
        class scheduler
        {
            private:
                struct ExposureRecord
                {
                    double                         t0;   ///< elapsed time since this exposure (incremented after each compute)
                    std::shared_ptr<DSL::FITScube> flux; ///< stored charge image
                };

                std::unique_ptr<pixelArray>  m_pixelModel;
                std::vector<ExposureRecord>  m_fluxHistory;
                const phot2pix::PixelConfig& m_cfg;
                size_t                       m_threadCount{0}; ///< 0 = use pixelArray default

            public:
                explicit scheduler(const phot2pix::PixelConfig& cfg)
                    : m_pixelModel(nullptr), m_fluxHistory(), m_cfg(cfg) {}
                virtual ~scheduler() = default;

                // Limit the thread budget used by the internal pixelArray.
                // Mirrors DetectorSimulator::setThreadCount(); pass 0 to restore the default.
                // Has no effect before the first insertCharges() call (pixelArray not yet created).
                void setThreadCount(size_t n)
                {
                    m_threadCount = n;
                    if (m_pixelModel) m_pixelModel->setThreadCount(n);
                }

                // Store a new exposure. t0 is the elapsed time since this exposure
                // (0 if it just finished, positive if there was a gap before insertion).
                void insertCharges(const std::shared_ptr<DSL::FITScube>& flux, const double& t0);

                // Remove the entry whose current elapsed time matches t0.
                void removeCharges(const double& t0);

                const std::vector<ExposureRecord>& getHistory() const { return m_fluxHistory; }

                // Compute the summed persistence signal for a new exposure that starts
                // dT seconds after the last inserted exposure and lasts exposureTime seconds.
                std::shared_ptr<DSL::FITScube> computePersistenceSignal(const double& dT, const double& exposureTime);
                std::shared_ptr<DSL::FITScube> computePersistenceSignal(const double& dT, const double& exposureTime, const size_t& c0, const size_t& r0, const size_t& ncols, const size_t& nrows);

        };
#pragma endregion

    } // namespace persistence
#pragma endregion

// ── NISP SCA ──────────────────────────────────────────────────────────────────
#pragma region - NISP SCA
    // SCA uses private inheritance so that DetectorSimulator::simulate() is not
    // directly reachable by callers. The only public simulation entry point is
    // expose(), which handles persistence, then converts electrons to uint16 ADU.
    class SCA : private phot2pix::DetectorSimulator<nisp::pixel>
    {
    private:
        double dpuOffset = 0.0;
        double dpuScale  = 1.0;
        double dpuGain   = 2.0;  // Mean gain in ADU/electron for the DPU flux estimation algorithm
        double dpuNoise  = 10.0; // Mean read noise in electrons (RMS) for the DPU flux estimation algorithm
        double frameTime = 1.45408; ///< NISP HAWAII-2RG frame time (seconds)
        bool   m_enable_persistence = false;
        bool   m_enable_oneOverF = true;
        nisp::mode m_mode = nisp::mode::MACC;

        std::shared_ptr<persistence::scheduler> persistenceModelling;
        nisp::OneOverFNoise m_oneOverFNoise;

        // Apply DPU scale/offset and saturate to uint16.
        // Output dimensions are read from doubleFrame->Size(1/2).
        std::shared_ptr<DSL::FITScube> _toUint16(
            const std::shared_ptr<DSL::FITScube>& doubleFrame) const;

        void _engine(const std::valarray<double>& flux,
                     double                        exposureTime,
                     size_t                        begin,
                     size_t                        end,
                     std::valarray<double>&         frame) const override
        {
            for (size_t i = begin; i < end; ++i)
                frame[i] = at(i).pixelresponse(flux[i], exposureTime);
        }

        // oneOverFStream: when non-null (and 1/f noise enabled), consumes the
        // next chunk of a continuous per-ramp AR(1) process instead of
        // independently generating a fresh one. Passed by expose_macc() so
        // every read in a ramp shares one continuous noise stream rather than
        // each read restarting the amplifier state at zero.
        void _simulate(const std::valarray<double>& flux,
                          const double& exposureTime,
                          std::shared_ptr<DSL::FITScube>& frame, phot2pix::PixelConfig cfg = phot2pix::PixelConfig(),
                          nisp::OneOverFNoise::Stream* oneOverFStream = nullptr);

        void _simulate(const std::valarray<double>& flux,
                          const double& exposureTime, 
                          std::shared_ptr<DSL::FITScube>& frame,
                          size_t c0, size_t r0, size_t ncols, size_t nrows, phot2pix::PixelConfig cfg = phot2pix::PixelConfig());

        std::shared_ptr<DSL::FITScube> expose_macc(const std::valarray<double>& flux,
                          const nisp::macc& maccParams);

        std::shared_ptr<DSL::FITScube> expose_cds(const std::valarray<double>& flux,
                          const double& exposureTime);

        std::shared_ptr<DSL::FITScube> expose_flux(const std::valarray<double>& flux,
                          const double& exposureTime);

        std::shared_ptr<DSL::FITScube> expose_dpu(const std::valarray<double>& flux,
                           const nisp::macc& maccParams);

    public:
        // ── Selectively re-exported base API ─────────────────────────────────
        using phot2pix::DetectorSimulator<nisp::pixel>::getWidth;
        using phot2pix::DetectorSimulator<nisp::pixel>::getHeight;
        using phot2pix::DetectorSimulator<nisp::pixel>::size;
        using phot2pix::DetectorSimulator<nisp::pixel>::fill;
        using phot2pix::DetectorSimulator<nisp::pixel>::initialize; // exposes string overload

        // Const-only pixel read access; non-const write access is intentionally hidden.
        const nisp::pixel& at(size_t i) const
            { return phot2pix::DetectorSimulator<nisp::pixel>::at(i); }
        const nisp::pixel& at(size_t col, size_t row) const
            { return phot2pix::DetectorSimulator<nisp::pixel>::at(col, row); }

        // ── Constructors ─────────────────────────────────────────────────────
        // DOUBLE_IMG: simulate() produces electron counts in double precision;
        // expose() applies DPU scale/offset and converts to uint16 ADU.
        SCA()                   : phot2pix::DetectorSimulator<nisp::pixel>(2048, 2048, DOUBLE_IMG) {}
        SCA(size_t w, size_t h) : phot2pix::DetectorSimulator<nisp::pixel>(w, h, DOUBLE_IMG) {}
        virtual ~SCA() = default;

        // ── Configuration ─────────────────────────────────────────────────────
        void initialize(const phot2pix::PixelConfig& cfg) override;

        void insertPersistenceModel(const std::shared_ptr<persistence::scheduler>& model)
        { persistenceModelling = model; }

        // ── Exposure ─────────────────────────────────────────────────────────
        // Windowed exposure: only the [c0, r0, ncols, nrows] sub-region is simulated.
        // fullFlux must cover the entire detector and is used for scheduler bookkeeping.
        // Returns a uint16 ADU FITS image of size ncols × nrows.
        std::shared_ptr<DSL::FITScube> expose(
            const std::valarray<double>& fullFlux,
            double exposureTime,
            size_t c0, size_t r0, size_t ncols, size_t nrows);

        std::shared_ptr<DSL::FITScube> expose(
            const std::valarray<double>& fullFlux,
            const nisp::macc& maccParams);

        std::shared_ptr<DSL::FITScube> expose(
            const std::valarray<double>& fullFlux,
            const double& exposureTime);
        
        void resetSimulationState();
    };

#pragma endregion

} // namespace nisp

#endif // _NISP_SCA_H_
