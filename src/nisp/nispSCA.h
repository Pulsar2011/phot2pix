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
#include <optional>
#include <string>

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

        /// True if cfg holds MACC parameters, as "macc"=[ng,nr,nd] or as "ng","nr","nd" keys.
        static bool inConfig(const phot2pix::PixelConfig& cfg);

        /// Read MACC parameters from cfg, given either as "macc"=[ng,nr,nd] or as
        /// "ng"=ng, "nr"=nr, "nd"=nd (ng → nGroups, nr → nFrames, nd → nDrop).
        /// Throws std::invalid_argument if absent, if both forms are given, or if a
        /// value is not a non-negative integer (ng and nr must be >= 1).
        static macc fromConfig(const phot2pix::PixelConfig& cfg);
    };

    /// Readout mode named by an "exposureMode" value, case-insensitive: "MACC", "UTR", "CDS",
    /// "DPU" or "FLUX". Throws std::invalid_argument for any other name.
    mode modeFromString(const std::string& name);

    /// Parameters of one exposure (see SCA::expose(fullFlux, readout)): its timing, given either
    /// as an exposure time or as MACC parameters, and the readout mode it uses, if not the SCA's.
    struct readout
    {
        std::optional<nisp::mode> exposureMode;        ///< "exposureMode"; empty → the SCA's mode
        std::optional<macc>       maccParams;          ///< "macc"=[ng,nr,nd] or "ng","nr","nd"
        double                    exposureTime = 0.0;  ///< "exposureTime" (s), used when maccParams is empty

        /// Read from cfg: exactly one timing, "exposureTime" (s, > 0) or MACC parameters (see
        /// macc::fromConfig()), and an optional "exposureMode" (see modeFromString()).
        /// Throws std::invalid_argument if both timings or neither are given, or a value is invalid.
        static readout fromConfig(const phot2pix::PixelConfig& cfg);
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
                parameter _factor;    ///< index → release rate f_k per stored electron (1/s); f_k·τ_k is the fraction released in total
                parameter _decayTime; ///< index → trap decay time τ_k (seconds)

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
                double                       m_delay{0.0};     ///< "dT": time between the end of an exposure and the start of the next (s)

            public:
                // Reads the delay between exposures from the "persistence" object form,
                // {"traps": [...], "dT": seconds} (0 when absent). Throws
                // std::invalid_argument if "dT" is not a number >= 0.
                explicit scheduler(const phot2pix::PixelConfig& cfg);
                virtual ~scheduler() = default;

                // Delay between two exposures ("dT"), seconds.
                double getDelay() const { return m_delay; }

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
                // The next exposure starts getDelay() later: every stored exposure,
                // this one included, ages by that delay.
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
        double dpuGain   = 2.0;  // Mean gain in e-/ADU for the DPU flux estimation algorithm; set by initialize() ("dpu_gain" or mean pixel gain)
        double dpuNoise  = 10.0; // Mean read noise in electrons (RMS) for the DPU flux estimation algorithm; set by initialize() ("dpu_noise" or mean pixel read noise)
        double frameTime = 1.45408; ///< NISP HAWAII-2RG frame time (seconds)
        bool   m_enable_persistence = false;
        bool   m_enable_oneOverF = true;
        bool   m_refPixCorr = true;       ///< DPU reference pixel correction ("refPixCorr": true/false or "ON"/"OFF"; on by default)
        nisp::mode m_mode = nisp::mode::MACC;

        std::shared_ptr<persistence::scheduler> persistenceModelling;
        nisp::OneOverFNoise m_oneOverFNoise;

        // Physical detector and placement of this matrix on it ("detectorSize",
        // "detectorOrigin", "refBorder"), held by the base class. The constructors install
        // a nisp::Geometry and initialize() replaces it with a configured clone(), so the
        // dynamic type is always nisp::Geometry.
        const nisp::Geometry& _geometry() const
            { return static_cast<const nisp::Geometry&>(getGeometry()); }

        // ── Output header ────────────────────────────────────────────────────
        // Record the NISP inputs as header cards (setInputCard(), see
        // DetectorSimulator::getInputCards(); the pixel model is recorded by the base):
        // detector geometry (DETSIZE, REFBORD), reference pixels (REFDARK), readout
        // noise (CHANWID, OOF*, BIAS*, PED*) and persistence (PERDT, PERTRn).
        // Called at the end of initialize().
        void _recordInputs(const phot2pix::PixelConfig& cfg);

        // Write the keywords common to every output: DPU parameters (OFFSET, SFACTOR,
        // GAIN, RN), FRTIME, version, DETSEC (the detector section of the image section
        // [c0, c0+ncols) × [r0, r0+nrows) of the matrix), then every input card.
        void _writeHeader(DSL::FITScube& image, size_t c0, size_t r0, size_t ncols, size_t nrows) const;
        void _writeHeader(DSL::FITScube& image) const { _writeHeader(image, 0, 0, getWidth(), getHeight()); }

        // Apply DPU scale/offset and saturate to uint16.
        // Output dimensions are read from doubleFrame->Size(1/2).
        std::shared_ptr<DSL::FITScube> _toUint16(
            const std::shared_ptr<DSL::FITScube>& doubleFrame) const;

        // ── Charge engines ───────────────────────────────────────────────────
        // Every frame the SCA simulates is charge only, built from the pixel stages
        // resetLevel() + integrate(); read noise and 1/f noise are added per read by
        // _read() / _accumulate() through pixel::read(). All four DetectorSimulator
        // hooks are overridden: the default ones call pixelresponse(), which is a
        // complete single read (read noise included).
        //
        // cfg "baseline": absent or "on" → resetLevel() + integrate() (first frame after
        // a reset), "off" → integrate() only (later ramp increments, the reset level is
        // already in the accumulated charge). Decoded once per thread, not per pixel.
        static bool _withReset(const phot2pix::PixelConfig& cfg)
        {
            return !(cfg.has("baseline") && cfg.get<std::string>("baseline") == "off");
        }

        void _engine(const std::valarray<double>& flux,
                     double                        exposureTime,
                     size_t                        begin,
                     size_t                        end,
                     std::valarray<double>&         frame) const override
        {
            for (size_t i = begin; i < end; ++i)
                frame[i] = m_pixels[i].resetLevel() + m_pixels[i].integrate(flux[i], exposureTime);
        }

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

        // Windowed hooks: flux and frame are window-sized (ncols × nrows, row-major).
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

        // Charge collected over exposureTime (ADU, double): pixel charge (resetLevel() +
        // integrate(), see the engines above) + persistence. No read noise and no 1/f noise:
        // those belong to a read and are added by _read() / _accumulate(), so they never
        // accumulate along a ramp.
        void _simulate(const std::valarray<double>& flux,
                          const double& exposureTime,
                          std::shared_ptr<DSL::FITScube>& frame, phot2pix::PixelConfig cfg = phot2pix::PixelConfig());

        void _simulate(const std::valarray<double>& flux,
                          const double& exposureTime,
                          std::shared_ptr<DSL::FITScube>& frame,
                          size_t c0, size_t r0, size_t ncols, size_t nrows, phot2pix::PixelConfig cfg = phot2pix::PixelConfig());

        // Read a full-frame charge map in place (ADU): every pixel adds its own read
        // noise (pixel::read), then one independent readout-noise realisation when
        // enabled (1/f, bias voltage noise, pedestal drift; generated in electrons and
        // divided by each pixel's gain). tStart is the time since the reset at which
        // this frame's read starts (s), used by the pedestal drift.
        void _read(std::valarray<double>& charge, double tStart = 0.0) const;

        // Same for a window: the first ncols × nrows values of charge (row-major) are the
        // detector region [c0, c0+ncols) × [r0, r0+nrows).
        void _read(std::valarray<double>& charge, size_t c0, size_t r0, size_t ncols, size_t nrows, double tStart = 0.0) const;

        // Ramp step: cumul += charge (charge only). When layer is not null this is also a
        // read, written as layer[i] = uint16(pixel_i.read(cumul[i]) + oneOverF[i] / gain_i):
        // the read noise and readout noise go to that layer only, never into cumul.
        // oneOverF is the readout-noise map in electrons, from OneOverFNoise::Stream; it
        // may be null (readout noise disabled). Threaded over getThreadCount().
        void _accumulate(const std::valarray<double>& charge, std::valarray<double>& cumul,
                         uint16_t* layer, const std::valarray<double>* oneOverF) const;

        std::shared_ptr<DSL::FITScube> expose_macc(const std::valarray<double>& flux,
                          const nisp::macc& maccParams);

        std::shared_ptr<DSL::FITScube> expose_cds(const std::valarray<double>& flux,
                          const double& exposureTime);

        std::shared_ptr<DSL::FITScube> expose_flux(const std::valarray<double>& flux,
                          const double& exposureTime);

        std::shared_ptr<DSL::FITScube> expose_dpu(const std::valarray<double>& flux,
                           const nisp::macc& maccParams);

    public:

        enum refPixCorrection
        {
            NONE=0,
            TOP=0x1,
            BOTTOM=0x2,
            LEFT=0x4,
            RIGHT=0x8,
            RIGHTnLEFT=LEFT|RIGHT,
            TOPnBOTTOM=TOP|BOTTOM,
            ALL=TOP|BOTTOM|LEFT|RIGHT,
        };

        // Reference pixel correction, Eq. 3 of arXiv:2507.11326 (from Kubik et al. 2014),
        // for line y in [4, 2044) of output channel ch (64 columns):
        //   c_y^(ch) = ½ (T^(ch) + B^(ch)) + ½ (L_y + R_y)
        // T, B: median of the 4 top / bottom reference lines of the channel (without the
        // left / right reference columns); L_y, R_y: median of the 4 left / right reference
        // columns over lines y−4 … y+4. mode picks the terms (bitmask: TOP, BOTTOM or
        // TOPnBOTTOM, plus LEFT, RIGHT or RIGHTnLEFT; ALL = Eq. 3).
        //
        // Returns c of every read k of the cube (any integer or floating point type), in its
        // units: c[(y * 32 + ch) * nReads + k], the reads of one line and channel being
        // contiguous. Lines 0-3 and 2044-2047 hold 0 (not corrected). Empty if mode is NONE
        // or the image is smaller than 2048 × 2048.
        static std::valarray<double> getReferencePixelCorrection(const std::shared_ptr<DSL::FITScube>& ramp, const refPixCorrection& mode = refPixCorrection::ALL);

        // Applies getReferencePixelCorrection() to the cube in place, e.g. to a MACC ramp from
        // expose(). Floating point and signed cubes: v − c (Eq. 3), rounded and clamped for
        // integer types. Unsigned cubes cannot hold the negative values Eq. 3 gives on raw
        // reads (c ≈ 2 × baseline with ALL): they get v − (c_k − c_0), c_0 being the correction
        // of the first read. Differences between reads, and so slopes, are those of Eq. 3; the
        // level stays the raw one (a single-read unsigned image is therefore left unchanged).
        // Saturated reads of unsigned cubes (the type's largest value) are left unchanged.
        static void referencePixelsCorrection(std::shared_ptr<DSL::FITScube>& image, const refPixCorrection& mode = refPixCorrection::ALL);

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
        const nisp::pixel& at(size_t i) const
            { return phot2pix::DetectorSimulator<nisp::pixel>::at(i); }
        const nisp::pixel& at(size_t col, size_t row) const
            { return phot2pix::DetectorSimulator<nisp::pixel>::at(col, row); }

        // ── Constructors ─────────────────────────────────────────────────────
        // DOUBLE_IMG: simulate() produces electron counts in double precision;
        // expose() applies DPU scale/offset and converts to uint16 ADU.
        // The matrix is placed on a 2048 × 2048 HAWAII-2RG with 4 reference pixels on each
        // side (nisp::Geometry), at detector pixel (0, 0) unless "detectorOrigin" is given.
        SCA()                   : phot2pix::DetectorSimulator<nisp::pixel>(2048, 2048, DOUBLE_IMG, std::make_shared<nisp::Geometry>()) {}
        SCA(size_t w, size_t h) : phot2pix::DetectorSimulator<nisp::pixel>(w, h, DOUBLE_IMG, std::make_shared<nisp::Geometry>()) {}
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

        // Exposure with its own parameters (see readout::fromConfig()): expose(fullFlux, macc)
        // or expose(fullFlux, exposureTime), in the readout mode of r if it gives one. That mode
        // is used for this exposure only: the SCA's mode is restored afterwards.
        std::shared_ptr<DSL::FITScube> expose(
            const std::valarray<double>& fullFlux,
            const nisp::readout& r);
    };

#pragma endregion

} // namespace nisp

#endif // _NISP_SCA_H_
