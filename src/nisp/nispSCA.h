//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispSCA.h
//
/// @file nispSCA.h
/// @brief nisp::SCA, the NISP detector and its readout modes; its exposure parameters and its persistence model.

#ifndef _NISP_SCA_H_
#define _NISP_SCA_H_

#include "nispPixel.h"
#include "nispOneOverFNoise.h"
#include "nispRTN.h"
#include <vector>
#include <map>
#include <memory>
#include <algorithm>
#include <optional>
#include <string>

#include <phot2pix/phot2pix.h>
#include <phot2pix/DetectorEffect.h>
#include <phot2pix/CosmicDeposit.h>

/// @namespace nisp
/// @brief The Euclid NISP HAWAII-2RG detectors, built on the phot2pix core.
namespace nisp
{
    /// @brief Readout mode: what SCA::expose() computes.
    enum class mode
    {
        CDS,   ///< Correlated Double Sampling
        MACC,  ///< Multiple Accumulated Correlated Double Sampling
        UTR,   ///< Sample Up The Ramp
        DPU,   ///< Flux estimation emulating the onboard DPU algorithm
        FLUX   ///< Fast flux estimation (last frame only)
    };

    /// @brief MACC parameters: groups of frames that are read, separated by dropped frames.
    struct macc
    {
        size_t nGroups;  ///< Number of groups
        size_t nFrames;  ///< Number of frames to accumulate
        size_t nDrop;    ///< Number of frames to skip between accumulations

        /// @brief One group of nFrames reads.
        macc UTR(const size_t& nFrames) { return {1, nFrames, 0}; }
        /// @brief nGroups groups of nFrames reads, separated by nDrop dropped frames.
        macc MACC(const size_t& nGroups, const size_t& nFrames, const size_t& nDrop) { return {nGroups, nFrames, nDrop}; }
        /// @brief Two reads separated by nDrops dropped frames.
        macc CDS(const size_t nDrops) { return {2, 1, nDrops}; }

        /// @brief True if cfg holds MACC parameters, as "macc"=[ng,nr,nd] or as "ng","nr","nd" keys.
        static bool inConfig(const phot2pix::PixelConfig& cfg);

        /// @brief Read MACC parameters from cfg.
        /// @details Given either as "macc"=[ng,nr,nd] or as "ng"=ng, "nr"=nr, "nd"=nd
        ///          (ng → nGroups, nr → nFrames, nd → nDrop).
        /// @param cfg Configuration.
        /// @throws std::invalid_argument if absent, if both forms are given, or if a
        ///         value is not a non-negative integer (ng and nr must be >= 1).
        static macc fromConfig(const phot2pix::PixelConfig& cfg);
    };

    /// @brief Readout mode named by an "exposureMode" value, case-insensitive: "MACC", "UTR", "CDS",
    ///        "DPU" or "FLUX".
    /// @throws std::invalid_argument for any other name.
    mode modeFromString(const std::string& name);

    /// @brief Parameters of one exposure (see SCA::expose(fullFlux, readout)).
    /// @details Its timing, given either as an exposure time or as MACC parameters, and the readout
    ///          mode it uses, if not the SCA's.
    struct readout
    {
        std::optional<nisp::mode> exposureMode;        ///< "exposureMode"; empty → the SCA's mode
        std::optional<macc>       maccParams;          ///< "macc"=[ng,nr,nd] or "ng","nr","nd"
        double                    exposureTime = 0.0;  ///< "exposureTime" (s), used when maccParams is empty

        /// @brief Read the parameters of one exposure from cfg.
        /// @details Exactly one timing, "exposureTime" (s, > 0) or MACC parameters (see
        ///          macc::fromConfig()), and an optional "exposureMode" (see modeFromString()).
        /// @param cfg Configuration.
        /// @throws std::invalid_argument if both timings or neither are given, or a value is invalid.
        static readout fromConfig(const phot2pix::PixelConfig& cfg);
    };

// ── Persistence model ─────────────────────────────────────────────────────────
#pragma region - persistence model
    /// @brief Persistence: charge trapped during an exposure and released during the following ones.
    namespace persistence
    {

#pragma region - pixel model
        /// @brief Charge released by the traps of one pixel.
        /// @details Trap k releases a stored charge E at the rate E·f_k·exp(−t/τ_k), t being the time
        ///          since the end of the exposure that filled it.
        class model : private phot2pix::pixsim
        {
            protected:
                typedef std::map<size_t, double> parameter;   ///< Trap index → value
                double dT = 0.0; ///< Elapsed time since last exposure (seconds), used as fallback when model is used outside of a scheduler context (e.g. standalone pixelArray).

            private:
                parameter _factor;    ///< index → release rate f_k per stored electron (1/s); f_k·τ_k is the fraction released in total
                parameter _decayTime; ///< index → trap decay time τ_k (seconds)

            public:
                /// @brief No trap.
                model() : phot2pix::pixsim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0) {}
                /// @brief Traps of release rates factor (1/s per stored electron) and decay times decayTime (s).
                model(const parameter& factor, const parameter& decayTime): phot2pix::pixsim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0),_factor(factor), _decayTime(decayTime) {}
                /// @brief Same, with a captured fraction in [0, 1] (stored, not used by pixelresponse()).
                model(const double& capturedFraction, const parameter& factor, const parameter& decayTime): phot2pix::pixsim(0.0, 0.0, 1.0, capturedFraction, 0.0, 0.0),_factor(factor), _decayTime(decayTime) {}
                /// @brief Copy constructor.
                model(const model& other): phot2pix::pixsim(other), dT(other.dT),_factor(other._factor), _decayTime(other._decayTime) {}

                /// @brief Destructor.
                virtual ~model() = default;

                /// @brief Add a trap.
                /// @param factor    Release rate per stored electron (1/s).
                /// @param decayTime Decay time (s).
                void insertParameter(const double& factor, const double& decayTime);
                /// @brief Remove every trap.
                void clearTraps() { _factor.clear(); _decayTime.clear(); }

                /// @brief Release rates of the traps (1/s per stored electron).
                inline const parameter& getFactor()    const { return _factor; }
                /// @brief Decay times of the traps (s).
                inline const parameter& getDecayTime() const { return _decayTime; }

                /// @brief Captured fraction, in [0, 1].
                inline const double& getCapturedFraction() const { return getQuantumEfficiency(); }
                /// @brief Set the captured fraction. @throws std::invalid_argument if it is not in [0, 1].
                inline void setCapturedFraction(const double& fR) { setQuantumEfficiency(fR); }

                /// @brief Gain of the pixel the traps belong to (e-/ADU), 1 by default: pixelresponse()
                ///        returns the released charge divided by it.
                using phot2pix::pixsim::getGain;
                /// @brief Set the gain (e-/ADU). @throws std::invalid_argument if it is not valid.
                using phot2pix::pixsim::setGain;

                /// @brief Elapsed time since the last exposure (s).
                inline double getTimeSinceLastExposure() const { return dT; }
                /// @brief Add bydT seconds to the elapsed time since the last exposure.
                void incrementTime(const double& bydT) { dT += bydT; }

                /// @brief Released charge, with the elapsed time of the pixel (getTimeSinceLastExposure()).
                /// @param flux         Charge stored during the past exposure (e-).
                /// @param exposureTime Duration of the step (s).
                /// @return Charge released during the step (ADU: e- / getGain(); e- with the default gain 1).
                double pixelresponse(const double& flux, const double& exposureTime) const override;

                /// @brief Released charge, with the elapsed time "dT" of cfg.
                /// @details Falls back to the elapsed time of the pixel if cfg has no "dT". The scheduler
                ///          always provides it; the fallback supports standalone pixelArray use without
                ///          a scheduler.
                /// @param flux         Charge stored during the past exposure (e-).
                /// @param exposureTime Duration of the step (s).
                /// @param cfg          Runtime parameters: "dT", time from the end of the past exposure to
                ///                     the start of the step (s).
                /// @return Charge released during the step (ADU: e- / getGain(); e- with the default gain 1).
                double pixelresponse(const double& flux, const double& exposureTime, const phot2pix::PixelConfig& cfg) const override;
        };
#pragma endregion

#pragma region - pixel array
        /// @brief Array of persistence pixels.
        class pixelArray : public phot2pix::DetectorSimulator<nisp::persistence::model>
        {
            public:
                /// @brief 2048 × 2048 pixels.
                pixelArray() : phot2pix::DetectorSimulator<nisp::persistence::model>(2048, 2048) {}
                /// @brief w × h pixels.
                pixelArray(size_t w, size_t h)
                    : phot2pix::DetectorSimulator<nisp::persistence::model>(w, h) {}
                /// @brief Destructor.
                virtual ~pixelArray() = default;

                /// @brief Add dT seconds to the elapsed time of every pixel.
                void incrementTime(const double& dT);
                /// @brief Set the traps of every pixel from "persistence".
                /// @details Either [{"f", "t"}, …] or {"traps": [{"f", "t"}, …], "dT": s}. "dT" starts the
                ///          clock of the pixels, for a standalone array. Nothing without the key.
                /// @param cfg Configuration.
                void initialize(const phot2pix::PixelConfig& cfg) override;

                /// @brief Set the gain of every pixel (e-/ADU): the released charge is then in ADU.
                /// @param gain One value per pixel, row-major.
                /// @throws std::invalid_argument if gain has not one value per pixel, or a value is not valid.
                void setGain(const std::valarray<double>& gain);

            protected:
                /// @brief Full-array hook: charge released by the pixels [begin, end).
                void _engine(const std::valarray<double>& flux,
                             double dt, size_t begin, size_t end,
                             std::valarray<double>& frame,
                             const phot2pix::PixelConfig& cfg) const override;

                /// @brief Windowed hook: charge released by the window rows [row_begin, row_end).
                void _engine(const std::valarray<double>& flux,
                             double dt, size_t c0, size_t r0,
                             size_t row_begin, size_t row_end, size_t ncols,
                             std::valarray<double>& frame,
                             const phot2pix::PixelConfig& cfg) const override;

                /// @brief No charge source: the array models the charge stored by past exposures, not a
                ///        detector, so it never receives cosmic-ray hits (see DetectorSimulator::expose()).
                std::shared_ptr<phot2pix::ChargeDeposit> _makeChargeDeposit() const override { return nullptr; }
        };
#pragma endregion

#pragma region - scheduler
        /// @brief History of the past exposures, and the persistence they leave in the next ones.
        /// @details Holds one shared pixelArray (trap parameters) plus the flux history of past
        ///          exposures. Memory cost: O(N * W*H * sizeof(double)).
        /// @warning It keeps a reference to the configuration it is built from, which must outlive it.
        class scheduler
        {
            private:
                /// @brief One past exposure.
                struct ExposureRecord
                {
                    double                         t0;   ///< elapsed time since this exposure (incremented after each compute)
                    std::shared_ptr<DSL::FITScube> flux; ///< stored charge image
                };

                std::unique_ptr<pixelArray>  m_pixelModel;     ///< Trap parameters, created by the first insertCharges()
                std::vector<ExposureRecord>  m_fluxHistory;    ///< Past exposures
                const phot2pix::PixelConfig& m_cfg;            ///< Configuration (not owned)
                size_t                       m_threadCount{0}; ///< 0 = use pixelArray default
                double                       m_delay{0.0};     ///< "dT": time between the end of an exposure and the start of the next (s)
                std::valarray<double>        m_gain;           ///< Gain of every pixel (e-/ADU), given to the pixelArray; empty: 1

            public:
                /// @brief Empty history.
                /// @details Reads the delay between exposures from the "persistence" object form,
                ///          {"traps": [...], "dT": seconds} (0 when absent).
                /// @param cfg Configuration, kept by reference: it must outlive the scheduler.
                /// @throws std::invalid_argument if "dT" is not a number >= 0.
                explicit scheduler(const phot2pix::PixelConfig& cfg);
                /// @brief Destructor.
                virtual ~scheduler() = default;

                /// @brief Delay between two exposures ("dT"), seconds.
                double getDelay() const { return m_delay; }

                /// @brief Limit the thread budget used by the internal pixelArray.
                /// @details Mirrors DetectorSimulator::setThreadCount(). Before the first insertCharges()
                ///          call the value is kept, and applied when the pixelArray is created.
                /// @param n Number of threads; 0 restores the default.
                void setThreadCount(size_t n)
                {
                    m_threadCount = n;
                    if (m_pixelModel) m_pixelModel->setThreadCount(n);
                }

                /// @brief Gain of every pixel (e-/ADU): computePersistenceSignal() then returns ADU.
                /// @details Before the first insertCharges() call the gains are kept, and applied when the
                ///          pixelArray is created.
                /// @param gain One value per pixel, row-major; empty: 1 (released charge in e-).
                /// @throws std::invalid_argument if the pixelArray exists and gain has not one value per pixel
                ///         (nothing changed), or a value is not valid.
                void setGain(const std::valarray<double>& gain);

                /// @brief Store a new exposure.
                /// @details The next exposure starts getDelay() later: every stored exposure, this one
                ///          included, ages by that delay.
                /// @param flux Charge stored in every pixel (e-, double precision).
                /// @param t0   Elapsed time since this exposure (s): 0 if it just finished, positive if
                ///             there was a gap before insertion.
                /// @throws std::invalid_argument if its size is not that of the previous exposures.
                void insertCharges(const std::shared_ptr<DSL::FITScube>& flux, const double& t0);

                /// @brief Remove the entry whose current elapsed time matches t0 (s).
                void removeCharges(const double& t0);

                /// @brief The past exposures.
                const std::vector<ExposureRecord>& getHistory() const { return m_fluxHistory; }

                /// @brief Persistence signal of a step, summed over the past exposures.
                /// @details Advances the clock of every past exposure by dT + exposureTime.
                /// @param dT           Time from the end of the previous step to the start of this one (s).
                /// @param exposureTime Duration of the step (s).
                /// @return Released charge of every pixel (ADU with setGain(), e- otherwise); nullptr without past exposure.
                std::shared_ptr<DSL::FITScube> computePersistenceSignal(const double& dT, const double& exposureTime);
                /// @brief Same for the window [c0, c0+ncols) × [r0, r0+nrows).
                /// @return Released charge of the window pixels (ADU with setGain(), e- otherwise); nullptr without past exposure.
                std::shared_ptr<DSL::FITScube> computePersistenceSignal(const double& dT, const double& exposureTime, const size_t& c0, const size_t& r0, const size_t& ncols, const size_t& nrows);

        };
#pragma endregion

    } // namespace persistence
#pragma endregion

// ── Cosmic rays ───────────────────────────────────────────────────────────────
#pragma region - cosmic rays
    /// @brief Cosmic ray deposits in a NISP H2RG: the phot2pix model (phot2pix/CosmicDeposit.h) with a
    ///        Hg(1-x)Cd(x)Te sensitive layer, x = cdFraction.
    /// @details The material is not a configuration key: it is set here, for every NISP configuration.
    ///          The other parameters and the spectrum come from the "cosmic" object.
    class CosmicDeposit : public phot2pix::CosmicRay::CosmicDeposit
    {
    public:
        /// @brief Cd fraction x of the sensitive layer.
        static constexpr double cdFraction = 0.3;

        /// @brief Not configured (see phot2pix::CosmicRay::CosmicDeposit), HgCdTe(cdFraction) material.
        CosmicDeposit() { m_params.material = phot2pix::CosmicRay::target::material::HgCdTe(cdFraction); }

        /// @brief Copy of this source with its dynamic type.
        std::shared_ptr<phot2pix::ChargeDeposit> clone() const override { return std::make_shared<CosmicDeposit>(*this); }
    };
#pragma endregion

// ── NISP SCA ──────────────────────────────────────────────────────────────────
#pragma region - NISP SCA
    /// @brief NISP detector (sensor chip assembly): its pixels, inter-pixel capacitance, crosstalk, readout
    ///        noise, persistence and readout modes.
    /// @details Private inheritance: DetectorSimulator::simulate() is not directly reachable by
    ///          callers. The only public simulation entry point is expose(), which handles
    ///          persistence, reads the charge in the readout mode of the SCA and returns uint16 ADU.
    class SCA : private phot2pix::DetectorSimulator<nisp::pixel>
    {
    private:
        double dpuOffset = 0.0;  ///< "offset" (ADU): DN = value × scale + offset (not applied to MACC / UTR layers)
        double dpuScale  = 1.0;  ///< "scale": DN = value × scale + offset (not applied to MACC / UTR layers)
        double dpuGain   = 2.0;  ///< Mean gain in e-/ADU for the DPU flux estimation algorithm; set by initialize() ("dpu_gain" or mean pixel gain)
        double dpuNoise  = 10.0; ///< Mean read noise in electrons (RMS) for the DPU flux estimation algorithm; set by initialize() ("dpu_noise" or mean pixel read noise)
        double frameTime = 1.45408; ///< NISP HAWAII-2RG frame time (seconds)
        bool   m_enable_persistence = false;   ///< Persistence, on when "persistence" is configured
        bool   m_enable_oneOverF = true;       ///< Readout noise, when m_oneOverFNoise has a component
        bool   m_refPixCorr = true;       ///< DPU reference pixel correction ("refPixCorr": true/false or "ON"/"OFF"; on by default)

        nisp::mode m_mode = nisp::mode::MACC;   ///< Readout mode ("exposureMode")
        nonlinearity::functor m_nl;       ///< Non-linearity model of every pixel ("nlModel" or setNonLinearity()); empty: linear


        std::shared_ptr<persistence::scheduler> persistenceModelling;   ///< Past exposures; null without persistence
        nisp::OneOverFNoise m_oneOverFNoise;                            ///< Readout noise of the output channels
        nisp::RTN           m_rtn;                                      ///< Random telegraph noise of the pixels ("rtn"); none by default

        /// @brief Inter-pixel capacitance ("ipc", see PixelConvolution::initialize()).
        /// @details Applied before the read noise: from detector characterisation it is not possible to
        ///          decorrelate IPC from the read noise, which is measured with the IPC effect.
        phot2pix::PixelConvolution m_ipc;

        /// @brief Crosstalk between the output channels ("crosstalk", see _configureCrosstalk()).
        /// @details Gives the ghosts alone: static kernel, not normalized, null in its centre.
        phot2pix::PixelConvolution m_crosstalk;
        size_t                     m_crosstalkColumns = 0;   ///< Crosstalk: columns per output channel of the kernel
        std::vector<double>        m_crosstalkCouplings;     ///< Crosstalk: coupling at each distance in channels (1 at 0)

        nisp::CosmicDeposit m_cosmicDeposit;   ///< Cosmic ray deposit model ("cosmic"), HgCdTe sensitive layer; empty: no deposit

        /// @brief No charge source in DetectorSimulator: the SCA draws its cosmic rays itself
        ///        (m_cosmicDeposit, HgCdTe) and adds them at each step of its readout (_addCosmicHits()),
        ///        so the base initialize() leaves "cosmic" to SCA::initialize().
        std::shared_ptr<phot2pix::ChargeDeposit> _makeChargeDeposit() const override { return nullptr; }

        /// @brief Charge of the cosmic hits on the matrix during the current exposure (e-), for persistence.
        /// @details Filled by _addCosmicHits() when persistence is modelled, cleared at the start of each
        ///          exposure and added to the charge registered with insertCharges(); empty: no hit.
        std::valarray<double> m_cosmicCharge;

        /// @brief Physical detector and placement of this matrix on it ("detectorSize",
        ///        "detectorOrigin", "refBorder"), held by the base class.
        /// @details The constructors install a nisp::Geometry and initialize() replaces it with a
        ///          configured clone(), so the dynamic type is always nisp::Geometry.
        const nisp::Geometry& _geometry() const
            { return static_cast<const nisp::Geometry&>(getGeometry()); }

        // ── Output header ────────────────────────────────────────────────────

        /// @brief Record the NISP inputs as header cards (setInputCard()).
        /// @details Detector geometry (DETSIZE, REFBORD), reference pixels (REFDARK), full well and
        ///          non-linearity (PIXFW, NLMODEL, PIXNLn), readout noise (CHANWID, OOF*, BIAS*, PED*),
        ///          crosstalk (XTALKCOL, XTALKNCH, XTALKn)
        ///          and persistence (PERDT, PERTRn). The pixel model is recorded by the base, see
        ///          DetectorSimulator::getInputCards(). Called at the end of initialize().
        /// @param cfg Configuration.
        void _recordInputs(const phot2pix::PixelConfig& cfg);

        /// @brief Write the keywords common to every output.
        /// @details DPU parameters (OFFSET, SFACTOR, GAIN, RN), FRTIME, version, DETSEC (the detector
        ///          section of the image section [c0, c0+ncols) × [r0, r0+nrows) of the matrix), then
        ///          every input card.
        void _writeHeader(DSL::FITScube& image, size_t c0, size_t r0, size_t ncols, size_t nrows) const;
        /// @brief Same for an image of the whole matrix.
        void _writeHeader(DSL::FITScube& image) const { _writeHeader(image, 0, 0, getWidth(), getHeight()); }

        /// @brief Apply DPU scale/offset and saturate to uint16.
        /// @details Output dimensions are read from doubleFrame->Size(1/2).
        std::shared_ptr<DSL::FITScube> _toUint16(
            const std::shared_ptr<DSL::FITScube>& doubleFrame) const;

        // ── Charge engines ───────────────────────────────────────────────────
        // Every frame the SCA simulates is charge only, built from the pixel stages
        // resetLevel() + integrate(); read noise and 1/f noise are added per read by
        // _read() / _accumulate() through pixel::read(). All four DetectorSimulator
        // hooks are overridden: the default ones call pixelresponse(), which is a
        // complete single read (read noise included).

        /// @brief Whether a step starts with a reset, from cfg "baseline".
        /// @details Absent or "on" → resetLevel() + integrate() (first frame after a reset), "off" →
        ///          integrate() only (later ramp increments, the reset level is already in the
        ///          accumulated charge). Decoded once per thread, not per pixel.
        static bool _withReset(const phot2pix::PixelConfig& cfg)
        {
            return !(cfg.has("baseline") && cfg.get<std::string>("baseline") == "off");
        }

        /// @brief Full-array hook: reset level and charge of the pixels [begin, end) (ADU).
        void _engine(const std::valarray<double>& flux,
                     double                        exposureTime,
                     size_t                        begin,
                     size_t                        end,
                     std::valarray<double>&         frame) const override
        {
            for (size_t i = begin; i < end; ++i)
                frame[i] = m_pixels[i].resetLevel() + m_pixels[i].integrate(flux[i], exposureTime);
        }

        /// @brief Full-array hook: charge of the pixels [begin, end) (ADU), with the reset level unless
        ///        cfg says otherwise (see _withReset()).
        void _engine(const std::valarray<double>& flux,
                     double                        exposureTime,
                     size_t                        begin,
                     size_t                        end,
                     std::valarray<double>&         frame,
                     const phot2pix::PixelConfig&  cfg) const override
        {
            const bool withReset = _withReset(cfg);
            for (size_t i = begin; i < end; ++i)
                frame[i] = (withReset ? m_pixels[i].resetLevel() : 0.0) + m_pixels[i].integrate(flux[i], exposureTime, cfg);
        }

        /// @brief Windowed hook: reset level and charge of the window rows [row_begin, row_end) (ADU).
        /// @details flux and frame are window-sized (ncols × nrows, row-major).
        void _engine(const std::valarray<double>& flux,
                     double                        exposureTime,
                     size_t                        c0,
                     size_t                        r0,
                     size_t                        row_begin,
                     size_t                        row_end,
                     size_t                        ncols,
                     std::valarray<double>&         frame) const override
        {
            for (size_t r = row_begin; r < row_end; ++r)
                for (size_t c = 0; c < ncols; ++c)
                {
                    const nisp::pixel& px = m_pixels[(r0 + r) * m_width + c0 + c];
                    frame[r * ncols + c]  = px.resetLevel() + px.integrate(flux[r * ncols + c], exposureTime);
                }
        }

        /// @brief Windowed hook, with the reset level unless cfg says otherwise (see _withReset()).
        void _engine(const std::valarray<double>& flux,
                     double                        exposureTime,
                     size_t                        c0,
                     size_t                        r0,
                     size_t                        row_begin,
                     size_t                        row_end,
                     size_t                        ncols,
                     std::valarray<double>&         frame,
                     const phot2pix::PixelConfig&  cfg) const override
        {
            const bool withReset = _withReset(cfg);
            for (size_t r = row_begin; r < row_end; ++r)
                for (size_t c = 0; c < ncols; ++c)
                {
                    const nisp::pixel& px = m_pixels[(r0 + r) * m_width + c0 + c];
                    frame[r * ncols + c]  = (withReset ? px.resetLevel() : 0.0) + px.integrate(flux[r * ncols + c], exposureTime, cfg);
                }
        }

        /// @brief Charge collected over exposureTime (ADU, double): pixel charge + cosmic hits (_addCosmicHits()) + persistence.
        /// @details The pixel charge is resetLevel() + integrate(), see the engines above. No read
        ///          noise and no 1/f noise: those belong to a read and are added by _read() /
        ///          _accumulate(), so they never accumulate along a ramp.
        /// @param flux         Photons/s of every pixel.
        /// @param exposureTime Duration of the step (s).
        /// @param frame        Output, the charge map.
        /// @param cfg          Runtime parameters (see _withReset()).
        void _simulate(const std::valarray<double>& flux,
                          const double& exposureTime,
                          std::shared_ptr<DSL::FITScube>& frame, phot2pix::PixelConfig cfg = phot2pix::PixelConfig());

        /// @brief Same for the window [c0, c0+ncols) × [r0, r0+nrows); flux covers the whole matrix.
        void _simulate(const std::valarray<double>& flux,
                          const double& exposureTime,
                          std::shared_ptr<DSL::FITScube>& frame,
                          size_t c0, size_t r0, size_t ncols, size_t nrows, phot2pix::PixelConfig cfg = phot2pix::PixelConfig());

        /// @brief Add the cosmic hits of a step (CosmicDeposit::generate()) to the charge of a window.
        /// @details hit.idx is a physical pixel (Y * width + X). The hits inside the window
        ///          [c0, c0+ncols) × [r0, r0+nrows) are added to frame, divided by the gain of their pixel
        ///          (e- → ADU), reference pixels included. With persistence, every hit on the matrix is
        ///          also added to m_cosmicCharge (e-), which the expose functions register with the charge
        ///          of the flux: the hits leave persistence in the next exposures.
        /// @param frame Charge of the window (ADU), ncols × nrows, row-major.
        /// @param hits  Deposits of the step.
        void _addCosmicHits(std::valarray<double>& frame, const std::vector<phot2pix::Deposit>& hits,
                            size_t c0, size_t r0, size_t ncols, size_t nrows);

        /// @brief Give the gains of the pixels to the persistence scheduler, so that the persistence
        ///        signal added to the frames is in ADU. Nothing without a scheduler.
        void _setPersistenceGain();

        /// @brief Read a full-frame charge map in place (ADU).
        /// @details The pixels are first coupled (_couplePixels(): IPC, crosstalk). Every pixel is then read
        ///          (pixel::read(): full well, non-linearity, read noise), and one independent
        ///          readout-noise realisation is added when enabled (1/f, bias voltage noise, pedestal
        ///          drift; generated in electrons and divided by each pixel's gain).
        /// @param charge Accumulated charge (ADU), replaced by the read values.
        /// @param tStart Time since the reset at which this frame's read starts (s), used by the
        ///               pedestal drift.
        void _read(std::valarray<double>& charge, double tStart = 0.0) const;

        /// @brief Same for a window: the first ncols × nrows values of charge (row-major) are the
        ///        matrix region [c0, c0+ncols) × [r0, r0+nrows).
        void _read(std::valarray<double>& charge, size_t c0, size_t r0, size_t ncols, size_t nrows, double tStart = 0.0) const;

        /// @brief Ramp step: cumul += charge (charge only), and a read when layer is not null.
        /// @details The read is written as layer[i] = uint16(pixel_i.read(cumul[i], m_nl) +
        ///          oneOverF[i] / gain_i + rtn[i]): the read noise, readout noise and random telegraph
        ///          noise go to that layer only, never into cumul. With IPC or crosstalk
        ///          (_couplePixels()), a coupled copy of cumul is read instead of cumul, which is never
        ///          coupled itself. Threaded over getThreadCount().
        /// @param charge   Charge of the step (ADU).
        /// @param cumul    Accumulated charge (ADU), updated.
        /// @param layer    Layer the read is written to; null for dropped frames.
        /// @param oneOverF Readout-noise map in electrons, from OneOverFNoise::Stream; null when the
        ///                 readout noise is disabled.
        /// @param rtn      Random telegraph noise map in ADU, from RTN::Stream::writeTo(), not divided by
        ///                 the gain; null without RTN.
        void _accumulate(const std::valarray<double>& charge, std::valarray<double>& cumul,
                         uint16_t* layer, const std::valarray<double>* oneOverF,
                         const std::valarray<double>* rtn) const;


        /// @brief Ramp step: cumul += charge (charge only), and a read when layer is not null.
        /// @details The read is written as layer[i] = uint16(pixel_i.read(cumul[i], m_nl) +
        ///          oneOverF[i] / gain_i + rtn[i]): the read noise, readout noise and random telegraph
        ///          noise go to that layer only, never into cumul. With IPC or crosstalk
        ///          (_couplePixels()), a coupled copy of cumul is read instead of cumul, which is never
        ///          coupled itself. Threaded over getThreadCount().
        /// @param charge   Charge of the step (ADU).
        /// @param cumul    Accumulated charge (ADU), updated.
        /// @param layer    Layer the read is written to; null for dropped frames.
        /// @param oneOverF Readout-noise map in electrons, from OneOverFNoise::Stream; null when the
        ///                 readout noise is disabled.
        void _accumulate(const std::valarray<double>& charge, std::valarray<double>& cumul,
                         uint16_t* layer, const std::valarray<double>* oneOverF) const
                        {
                            _accumulate(charge, cumul, layer, oneOverF, nullptr);
                        }

        /// @brief True if the read of a pixel depends on other pixels (IPC or crosstalk).
        bool _pixelsAreCoupled() const { return m_ipc.isEnabled() || m_crosstalk.isEnabled(); }

        /// @brief Effects that couple the pixels during a read, in place: applyIPC(), then applyCrosstalk().
        /// @param signal Accumulated charge (ADU) of the window, replaced by the signal to read.
        /// @param c0     Matrix column of the first window pixel.
        /// @param r0     Matrix row of the first window pixel.
        /// @param ncols  Window columns.
        /// @param nrows  Window rows.
        void _couplePixels(std::valarray<double>& signal, size_t c0, size_t r0, size_t ncols, size_t nrows) const;

        /// @brief Read the channel crosstalk from cfg "crosstalk". Nothing changes if the key is absent.
        /// @verbatim
        ///   "crosstalk": [k0, k1, ...]
        ///   "crosstalk": {"kernel": [k0, k1, ...]}   → the same
        ///   "crosstalk": []                          → no crosstalk
        /// @endverbatim
        /// @details kd is the fraction of the signal of a pixel that is added to the pixel read at the same
        ///          time in a channel at a distance of d channels: same row, d × colPerChannel columns
        ///          away on each side (every channel is read in the same direction). The values are >= 0,
        ///          at most one per channel: the missing ones are 0, so [1, 0.5] only couples a channel with
        ///          its two neighbours. k0, the fraction of its own signal that a channel keeps, is 1 (no
        ///          loss). The convolution kernel is built here: 0 between the channels.
        ///          The channels are those of the instrument: top-level "colPerChannel" / "nChannels", see
        ///          nisp::Geometry.
        /// @param cfg Configuration.
        /// @throws std::invalid_argument if the crosstalk is not valid; it is then unchanged.
        void _configureCrosstalk(const phot2pix::PixelConfig& cfg);

        /// @brief Build the crosstalk kernel from the couplings and the channels of the geometry.
        /// @param couplings Coupling at each distance in channels, 1 at distance 0.
        /// @throws std::invalid_argument if they are more than the channels, or if the channels do not cover
        ///         the detector; the crosstalk is then unchanged.
        void _setCrosstalk(const std::vector<double>& couplings);

        /// @brief Check a non-linearity model against the full well and the parameters of every pixel.
        /// @details Nothing is checked if f is empty (linear).
        /// @param f Model.
        /// @throws std::invalid_argument unless every pixel has a finite full well, up to which f, with
        ///         the parameters of that pixel, is finite and increasing.
        void _checkResponse(const nonlinearity::functor& f) const;

        /// @brief MACC and UTR: the ramp, one layer per read.
        /// @return uint16 cube of W × H × (nGroups · nFrames) reads, without scale nor offset.
        /// @throws std::invalid_argument if nGroups or nFrames is 0.
        std::shared_ptr<DSL::FITScube> expose_macc(const std::valarray<double>& flux,
                          const nisp::macc& maccParams);

        /// @brief CDS: last read − first read, over exposureTime (s) rounded down to whole frames.
        /// @return uint16 image, with scale and offset.
        /// @throws std::invalid_argument if exposureTime is less than 2 frames.
        std::shared_ptr<DSL::FITScube> expose_cds(const std::valarray<double>& flux,
                          const double& exposureTime);

        /// @brief FLUX: a single read after exposureTime (s) rounded down to whole frames.
        /// @return uint16 image, with scale and offset.
        std::shared_ptr<DSL::FITScube> expose_flux(const std::valarray<double>& flux,
                          const double& exposureTime);

        /// @brief DPU: ramp, reference pixel correction, then slope of every pixel (Ramp::getSlope()).
        /// @details As on the flight DPU, a pixel with a saturated read (65535 in the raw ramp) is not
        ///          fitted. The planes of the output are:
        ///          - 0: slope per group interval, × scale + offset; 65535 for a saturated pixel (no scale, no offset);
        ///          - 1: chi2 (not reduced), 0 with 2 groups or for a saturated pixel;
        ///          - 2: squared error of the slope, × 100 x scale + offset; 0 for a saturated pixel.
        /// @return uint16 cube of W × H × 3.
        /// @throws std::invalid_argument with less than 2 groups.
        std::shared_ptr<DSL::FITScube> expose_dpu(const std::valarray<double>& flux,
                           const nisp::macc& maccParams);

    public:

        /// @brief Terms of the reference pixel correction (bitmask), see getReferencePixelCorrection().
        enum refPixCorrection
        {
            NONE=0,                      ///< No correction
            TOP=0x1,                     ///< Top reference lines
            BOTTOM=0x2,                  ///< Bottom reference lines
            LEFT=0x4,                    ///< Left reference columns
            RIGHT=0x8,                   ///< Right reference columns
            RIGHTnLEFT=LEFT|RIGHT,       ///< Mean of the left and right reference columns
            TOPnBOTTOM=TOP|BOTTOM,       ///< Mean of the top and bottom reference lines
            ALL=TOP|BOTTOM|LEFT|RIGHT,   ///< Every term: Eq. 3 of arXiv:2507.11326
        };

        /// @brief Reference pixel correction, Eq. 3 of arXiv:2507.11326 (from Kubik et al. 2014).
        /// @details For line y in [b, H − b) of output channel ch, with the detector of H lines, its
        ///          border of b reference pixels and its nCh channels (NISP: 2048, 4, 32 of 64 columns):
        ///            c_y^(ch) = ½ (T^(ch) + B^(ch)) + ½ (L_y + R_y)
        ///
        ///          T, B: median of the b top / bottom reference lines of the channel (without the
        ///          left / right reference columns); L_y, R_y: median of the b left / right reference
        ///          columns over lines y−4 … y+4 (those of the detector).
        /// @param ramp     Cube of reads (any integer or floating point type), of the full detector.
        /// @param mode     Terms (bitmask: TOP, BOTTOM or TOPnBOTTOM, plus LEFT, RIGHT or RIGHTnLEFT;
        ///                 ALL = Eq. 3).
        /// @param geometry Detector size, reference border and channels. Default: NISP.
        /// @return c of every read k of the cube, in its units: c[(y * nCh + ch) * nReads + k], the reads
        ///         of one line and channel being contiguous. The reference lines hold 0 (not
        ///         corrected). Empty if mode is NONE, or, with a warning, if it does not apply: image that
        ///         is not the full detector, no reference border, channels that do not cover the detector
        ///         or have no column outside the border.
        static std::valarray<double> getReferencePixelCorrection(const std::shared_ptr<DSL::FITScube>& ramp, const refPixCorrection& mode = refPixCorrection::ALL,
                                                                 const nisp::Geometry& geometry = nisp::Geometry());

        /// @brief Apply getReferencePixelCorrection() to a cube in place, e.g. to a MACC ramp from expose().
        /// @details Floating point and signed cubes: v − c (Eq. 3), rounded and clamped for
        ///          integer types. Unsigned cubes cannot hold the negative values Eq. 3 gives on raw
        ///          reads (c ≈ 2 × baseline with ALL): they get v − (c_k − c_0), c_0 being the correction
        ///          of the first read. Differences between reads, and so slopes, are those of Eq. 3; the
        ///          level stays the raw one (a single-read unsigned image is therefore left unchanged).
        ///          Saturated reads of unsigned cubes (the type's largest value) are left unchanged.
        /// @param image    Cube of reads, corrected in place.
        /// @param mode     Terms of the correction.
        /// @param geometry Detector size, reference border and channels. Default: NISP.
        static void referencePixelsCorrection(std::shared_ptr<DSL::FITScube>& image, const refPixCorrection& mode = refPixCorrection::ALL,
                                              const nisp::Geometry& geometry = nisp::Geometry());

        // ── Selectively re-exported base API ─────────────────────────────────
        using phot2pix::DetectorSimulator<nisp::pixel>::getWidth;
        using phot2pix::DetectorSimulator<nisp::pixel>::getHeight;
        using phot2pix::DetectorSimulator<nisp::pixel>::size;
        using phot2pix::DetectorSimulator<nisp::pixel>::fill;
        using phot2pix::DetectorSimulator<nisp::pixel>::initialize; // exposes string overload
        using phot2pix::DetectorSimulator<nisp::pixel>::getInputCards; // inputs written in the output header
        using phot2pix::DetectorSimulator<nisp::pixel>::getGeometry;   // physical detector (a nisp::Geometry)
        using phot2pix::DetectorSimulator<nisp::pixel>::getOrigin;     // (x0, y0) of this matrix on the detector
        using phot2pix::DetectorSimulator<nisp::pixel>::toMatrix;      // matrix- or detector-sized image → matrix values
        using phot2pix::DetectorSimulator<nisp::pixel>::loadImage;     // FITS image named by a config key ({"file", "hdu"}, …)
        using phot2pix::DetectorSimulator<nisp::pixel>::loadMap;       // config key (number or FITS image) → matrix values

        // Const-only pixel read access; non-const write access is intentionally hidden.

        /// @brief Pixel of flat index i (row-major), read only.
        const nisp::pixel& at(size_t i) const
            { return phot2pix::DetectorSimulator<nisp::pixel>::at(i); }
        /// @brief Pixel at column col and row row, read only.
        const nisp::pixel& at(size_t col, size_t row) const
            { return phot2pix::DetectorSimulator<nisp::pixel>::at(col, row); }

        // ── Constructors ─────────────────────────────────────────────────────
        // DOUBLE_IMG: simulate() produces ADU in double precision;
        // expose() applies DPU scale/offset and converts to uint16 ADU.
        // The matrix is placed on a 2048 × 2048 HAWAII-2RG with 4 reference pixels on each
        // side (nisp::Geometry), at detector pixel (0, 0) unless "detectorOrigin" is given.

        /// @brief Full detector, 2048 × 2048 pixels.
        SCA()                   : phot2pix::DetectorSimulator<nisp::pixel>(2048, 2048, DOUBLE_IMG, std::make_shared<nisp::Geometry>()) {}
        /// @brief Matrix of w × h pixels, part of a 2048 × 2048 detector.
        SCA(size_t w, size_t h) : phot2pix::DetectorSimulator<nisp::pixel>(w, h, DOUBLE_IMG, std::make_shared<nisp::Geometry>()) {}
        /// @brief Destructor.
        virtual ~SCA() = default;

        // ── Configuration ─────────────────────────────────────────────────────

        /// @brief Configure the detector: geometry and output channels, pixels, readout, persistence, readout
        ///        noise, random telegraph noise, reference pixels, full well and non-linearity, inter-pixel
        ///        capacitance, crosstalk.
        /// @details The keys are described in the wiki (NISP example, configuration reference). With
        ///          "persistence", cfg must outlive the SCA: the scheduler keeps a reference to it.
        /// @param cfg Configuration; absent keys leave the current value unchanged.
        /// @throws std::invalid_argument if a value is not valid.
        void initialize(const phot2pix::PixelConfig& cfg) override;

        /// @brief Replace the persistence scheduler, e.g. to share the history of another SCA.
        /// @details Only used when persistence is enabled ("persistence" in the configuration). The
        ///          scheduler takes the gains of this SCA (_setPersistenceGain()).
        void insertPersistenceModel(const std::shared_ptr<persistence::scheduler>& model)
        { persistenceModelling = model; _setPersistenceGain(); }

        /// @brief Set the non-linearity model applied at each read.
        /// @details Call it after initialize(), which sets the full wells and the parameters of the
        ///          pixels: the model is checked against them (see _checkResponse()) and replaces the
        ///          one of "nlModel".
        /// @param f Model; an empty functor is linear.
        /// @throws std::invalid_argument if f is not valid for a pixel; the model is then unchanged.
        void setNonLinearity(const nonlinearity::functor& f)
        {
            _checkResponse(f);
            m_nl = f;
            setInputCard({"NLMODEL", std::string(f ? "custom" : "linear"), "Non-linearity model"});
        }

        // ── Random telegraph noise ────────────────────────────────────────────

        /// @brief Random telegraph noise of the pixels ("rtn"): the RTN pixels drawn by initialize().
        const nisp::RTN& getRTN() const { return m_rtn; }

        // ── Inter-pixel capacitance ───────────────────────────────────────────

        /// @brief True if an IPC kernel is configured ("ipc").
        bool isIPCEnabled() const { return m_ipc.isEnabled(); }
        /// @brief True if the IPC kernel depends on the pixel (FITS cube).
        bool isIPCDynamic() const { return m_ipc.isDynamic(); }
        /// @brief Static IPC kernel (columns, rows), normalized to 1; 1 × 1 without static kernel.
        const phot2pix::PixelConvolution::kernel& getIPCKernel() const { return m_ipc.getKernel(); }
        /// @brief Dynamic IPC kernel, W × H × (n + 1) float cube normalized pixel by pixel; null without dynamic kernel.
        std::shared_ptr<const DSL::FITScube> getIPCKernelCube() const { return m_ipc.getKernelCube(); }

        /// @brief Apply the IPC to a window of the matrix, in place. Nothing is done without IPC.
        /// @param image Window, ncols × nrows pixels, row-major.
        /// @param c0    Matrix column of the first window pixel.
        /// @param r0    Matrix row of the first window pixel.
        /// @param ncols Window columns.
        /// @param nrows Window rows.
        /// @throws std::out_of_range if the window exceeds the matrix, std::invalid_argument if image
        ///         is smaller than the window.
        void applyIPC(std::valarray<double>& image, size_t c0, size_t r0, size_t ncols, size_t nrows) const
        {
            m_ipc.apply(image, c0, r0, ncols, nrows, getThreadCount());
        }

        // ── Channel crosstalk ─────────────────────────────────────────────────

        /// @brief True if a crosstalk between the output channels is configured ("crosstalk").
        bool isCrosstalkEnabled() const { return m_crosstalk.isEnabled(); }
        /// @brief Crosstalk kernel: one row, the coupling of the channels at their distance in columns, 0
        ///        elsewhere; 1 × 1 without crosstalk.
        const phot2pix::PixelConvolution::kernel& getCrosstalkKernel() const { return m_crosstalk.getKernel(); }

        /// @brief Add the crosstalk ghosts to a window of the matrix, in place. Nothing is done without crosstalk.
        /// @details The ghosts are those of the signal above the baseline that the source pixel outputs,
        ///          limited by its full well and its non-linearity (as in pixel::read()). The source itself
        ///          is unchanged. A window has no ghost of the pixels outside of it.
        /// @param signal Window (ADU), ncols × nrows pixels, row-major.
        /// @param c0     Matrix column of the first window pixel.
        /// @param r0     Matrix row of the first window pixel.
        /// @param ncols  Window columns.
        /// @param nrows  Window rows.
        /// @throws std::out_of_range if the window exceeds the matrix, std::invalid_argument if signal
        ///         is smaller than the window.
        void applyCrosstalk(std::valarray<double>& signal, size_t c0, size_t r0, size_t ncols, size_t nrows) const;

        // ── Exposure ─────────────────────────────────────────────────────────

        /// @brief Windowed exposure: a single read of a sub-region, whatever the readout mode.
        /// @details The inter-pixel capacitance and the crosstalk are those of the pixels of the window: what
        ///          is outside of it gives nothing.
        /// @param fullFlux     Photons/s of every pixel of the matrix, row-major: the whole of it is
        ///                     stored for the persistence of the next exposures.
        /// @param exposureTime Exposure time (s), not rounded to whole frames.
        /// @param c0           Matrix column of the first window pixel.
        /// @param r0           Matrix row of the first window pixel.
        /// @param ncols        Window columns.
        /// @param nrows        Window rows.
        /// @return uint16 image of ncols × nrows pixels (ADU, with scale and offset).
        /// @pre The window lies within the matrix.
        std::shared_ptr<DSL::FITScube> expose(
            const std::valarray<double>& fullFlux,
            double exposureTime,
            size_t c0, size_t r0, size_t ncols, size_t nrows);

        /// @brief Exposure with MACC parameters, in the readout mode of the SCA.
        /// @details MACC, UTR and DPU use the pattern; CDS and FLUX use its duration,
        ///          (nGroups·(nFrames + nDrop) − nDrop) frames.
        /// @param fullFlux   Photons/s of every pixel of the matrix, row-major.
        /// @param maccParams MACC parameters.
        /// @return uint16: cube of reads (MACC, UTR), image (CDS, FLUX), or slope, chi2 and error
        ///         planes (DPU).
        std::shared_ptr<DSL::FITScube> expose(
            const std::valarray<double>& fullFlux,
            const nisp::macc& maccParams);

        /// @brief Exposure of a given duration, in the readout mode of the SCA.
        /// @details MACC, UTR and DPU use a ramp of exposureTime / frame time groups of one read.
        /// @param fullFlux     Photons/s of every pixel of the matrix, row-major.
        /// @param exposureTime Exposure time (s).
        /// @return uint16, as expose(fullFlux, maccParams).
        std::shared_ptr<DSL::FITScube> expose(
            const std::valarray<double>& fullFlux,
            const double& exposureTime);

        /// @brief Exposure with its own parameters (see readout::fromConfig()).
        /// @details expose(fullFlux, macc) or expose(fullFlux, exposureTime), in the readout mode of r
        ///          if it gives one. That mode is used for this exposure only: the SCA's mode is
        ///          restored afterwards.
        /// @param fullFlux Photons/s of every pixel of the matrix, row-major.
        /// @param r        Timing and readout mode.
        /// @return uint16, as expose(fullFlux, maccParams).
        std::shared_ptr<DSL::FITScube> expose(
            const std::valarray<double>& fullFlux,
            const nisp::readout& r);
    };

#pragma endregion

} // namespace nisp

#endif // _NISP_SCA_H_
