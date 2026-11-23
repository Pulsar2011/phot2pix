//
// Created by gillard on 2026-06-26.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispOneOverFNoise.h
//
/// @file nispOneOverFNoise.h
/// @brief Readout noise of the NISP HAWAII-2RG output channels, added at each read (nisp::OneOverFNoise).
///
/// Three components, each enabled by its own configuration key:
///
/// - oneOverF: 1/f (flicker / pink) noise of each channel's own amplifier.
///   Independent between channels.
/// - bias_voltage_noise: 1/f noise of the bias voltages shared by every channel: ONE
///   time series, seen by all channels at the same pixel clock.
///   This is what the left/right reference columns track.
/// - pedestal_drift: offset of each channel: a mean drift f(t), common to all
///   channels and given as an equation, plus white noise drawn
///   per channel and per frame. This is what the top/bottom
///   reference rows track.
///
/// All amplitudes are in electrons: the generated maps are in electrons and
/// nisp::SCA converts them to ADU with each pixel's gain, as it does for read noise.
///
/// @par Physical detector
/// The noise is defined on the physical detector (nisp::Geometry, 2048 × 2048 by
/// default), not on the simulated matrix. The detector is read through
/// width / colPerChannel output channels (32 × 64 columns for NISP); channel c covers
/// the physical columns [c·w, (c+1)·w), w = colPerChannel. A simulated matrix is
/// placed on the detector (at (0, 0) by default), and each of its pixels takes the noise
/// of the physical pixel it covers. The channel boundaries, the readout timing and the
/// spectrum of the noise are therefore the same whatever the size of the matrix: a
/// smaller matrix is only a smaller part of the same detector.
///
/// @par Time model
/// All channels are read in parallel, in the same direction (increasing column
/// within a row, then increasing row), with no row or frame overhead. Within one
/// frame, channel c reads physical pixel (X, Y) at pixel-clock sample
/// @verbatim
///     t = Y·w + (X − c·w)            (one frame = L = height · w samples)
/// @endverbatim
/// and frame k after the reset (k = 0, 1, …, read or dropped) starts at k·L. In
/// seconds, row Y of frame k is read at  k·t_f + (Y + 0.5)/height · t_f.
///
/// @par 1/f noise
/// PSD approximation by a sum of nPoles AR(1) processes with corner
/// frequencies log-spaced in [1/(2L), 0.5] cycles per sample — no FFT dependency.
/// Each pole carries a variance ∝ f_k^(1 − slope), so that equal log steps give
/// @verbatim
///   S(f) ≈ A² / f^slope          (slope 1 = pink, 2 = Brownian / random walk)
/// @endverbatim
/// Only the samples inside the matrix (or window) are generated: the samples of the
/// rows and columns outside it are skipped exactly, an AR(1) process being advanced
/// over n samples in one step,
/// @verbatim
///   x ← ρⁿ x + σ √((1 − ρ²ⁿ)/(1 − ρ²)) ε,
/// @endverbatim
/// so the cost scales with the matrix area. Every process starts from its stationary
/// state: the amplifier noise does not reset with the pixels.
///
/// The bias voltage series is indexed by the pixel clock (t above), so at a given
/// sample every channel gets the same value. The pedestal drift f(t) is evaluated once
/// per row, at the row's read time.

#ifndef _NISP_ONE_OVER_F_NOISE_H_
#define _NISP_ONE_OVER_F_NOISE_H_

#include "nispGeometry.h"

#include <phot2pix/PixelConfig.h>
#include <cstddef>
#include <valarray>
#include <vector>

namespace nisp
{

// ── OneOverFNoise ──────────────────────────────────────────────────────────────
#pragma region - OneOverFNoise

    /// @brief Readout noise of the output channels: 1/f noise of each channel, bias voltage noise and
    ///        pedestal drift (see nispOneOverFNoise.h).
    class OneOverFNoise
    {
    public:
        /// @brief One exponential term of the pedestal drift: amplitude (e⁻) · exp(−t/tau), tau in s.
        struct DriftExp { double amplitude; double tau; };
        /// @brief One periodic term of the pedestal drift: amplitude (e⁻) · sin(2π t/period + phase),
        ///        period in s, phase in rad.
        struct DriftSin { double amplitude; double period; double phase; };

        /// @brief Default-constructed instance is disabled (every component off).
        OneOverFNoise() = default;

        /// @brief Per-channel 1/f noise only.
        /// @details bias_voltage_noise and pedestal_drift are configured through initialize().
        /// @param amplitude     Per-channel RMS in electrons (typical NISP: 5—10 e⁻).
        /// @param slope         PSD exponent α (1.0 = pure 1/f, 2.0 = Brownian / random walk).
        /// @param nChannels     Number of H2RG output channels; 0 = detector width / colPerChannel.
        ///                      When given, it must equal detector width / colPerChannel.
        /// @param colPerChannel Channel width in pixels (64 for NISP); must divide the detector width.
        /// @param nPoles        AR(1) poles for PSD approximation (about 3 per decade: 16 covers the
        ///                      5 decades of a frame).
        explicit OneOverFNoise(double amplitude,
                               double slope        = 1.0,
                               size_t nChannels    = 0,
                               size_t colPerChannel = 64,
                               size_t nPoles       = 16);

        /// @brief Configure from a PixelConfig. Each key is optional; absent keys leave
        ///        their component unchanged.
        /// @details
        /// @verbatim
        /// Per-channel 1/f noise, scalar shorthand (amplitude only) or object:
        ///   "oneOverF": 7.0
        ///   "oneOverF": {
        ///     "amplitude"    : 7.0,
        ///     "slope"        : 1.0,
        ///     "nChannels"    : 0,
        ///     "colPerChannel": 64,
        ///     "nPoles"       : 16
        ///   }
        ///
        /// Bias voltage noise (1/f common to all channels), scalar shorthand or object:
        ///   "bias_voltage_noise": 3.0
        ///   "bias_voltage_noise": { "amplitude": 3.0, "slope": 1.0, "nPoles": 16 }
        ///
        /// Pedestal drift, scalar shorthand (white noise only) or object:
        ///   "pedestal_drift": 4.0
        ///   "pedestal_drift": {
        ///     "sigma": 4.0,                  // RMS around the mean, per channel per frame (e⁻)
        ///     "t0"   : 0.0,                  // time added to the time since reset (s)
        ///     "drift": {                     // mean offset f(t) in e⁻, every term optional
        ///       "poly": [a0, a1, a2],                                    // Σ a_k t^k
        ///       "exp" : [{"amplitude": A, "tau": τ}],                    // Σ A exp(−t/τ)
        ///       "sin" : [{"amplitude": B, "period": P, "phase": φ}]      // Σ B sin(2π t/P + φ), φ in rad
        ///     }
        ///   }
        /// @endverbatim
        /// "exp" and "sin" also accept a single object instead of an array.
        ///
        /// The detector geometry is not read here: nisp::SCA passes it with setGeometry().
        /// @param cfg Configuration.
        /// @throws std::invalid_argument on a malformed or out-of-range value.
        void initialize(const phot2pix::PixelConfig& cfg);

        /// @brief Set the frame time t_f (default: NISP 1.45408 s). nisp::SCA sets it from its own "frame_time".
        /// @param frameTime Frame time (s), > 0.
        /// @throws std::invalid_argument if not > 0.
        void setFrameTime(double frameTime);

        /// @brief Set the physical detector and the placement of the matrix (default: 2048 × 2048,
        ///        origin (0, 0)). nisp::SCA sets it from its own configuration.
        void setGeometry(const Geometry& geometry) { m_geometry = geometry; }

        /// @brief Check that a width × height matrix fits on the detector and that the channel
        ///        layout matches it (colPerChannel divides the detector width, nChannels agrees).
        /// @throws std::invalid_argument otherwise; nisp::SCA calls it from initialize().
        void checkGeometry(size_t width, size_t height) const;

        // ── Accessors ──────────────────────────────────────────────────────────

        /// @brief True when at least one component produces noise.
        bool   isEnabled()     const
        {
            return m_amplitude > 0.0 || m_biasAmplitude > 0.0 || m_pedestalSigma > 0.0 || hasPedestalDrift();
        }
        double getAmplitude()  const { return m_amplitude; }          ///< Per-channel 1/f noise: RMS (e⁻)
        double getSlope()      const { return m_slope; }              ///< Per-channel 1/f noise: PSD exponent
        size_t getNChannels()  const { return m_nChannels; }          ///< Number of channels; 0 = detector width / colPerChannel
        size_t getColPerChannel() const { return m_colPerChannel; }   ///< Columns per channel
        size_t getNPoles()     const { return m_nPoles; }             ///< Per-channel 1/f noise: AR(1) poles
        double getFrameTime()  const { return m_frameTime; }          ///< Frame time (s)
        const Geometry& getGeometry() const { return m_geometry; }    ///< Physical detector and placement of the matrix

        double getBiasVoltageNoise() const { return m_biasAmplitude; }   ///< Bias voltage noise: RMS (e⁻)
        double getBiasSlope()        const { return m_biasSlope; }       ///< Bias voltage noise: PSD exponent
        size_t getBiasNPoles()       const { return m_biasPoles; }       ///< Bias voltage noise: AR(1) poles

        double getPedestalSigma()    const { return m_pedestalSigma; }   ///< Pedestal: RMS per channel and frame (e⁻)
        double getPedestalT0()       const { return m_pedestalT0; }      ///< Pedestal: time added to the time since reset (s)
        /// @brief True if a drift equation is configured.
        bool   hasPedestalDrift()    const { return !m_driftPoly.empty() || !m_driftExp.empty() || !m_driftSin.empty(); }
        const std::vector<double>&   getDriftPoly() const { return m_driftPoly; }   ///< [a0, a1, …], e⁻ / e⁻ s⁻¹ / …
        const std::vector<DriftExp>& getDriftExp()  const { return m_driftExp; }    ///< Exponential terms of the drift
        const std::vector<DriftSin>& getDriftSin()  const { return m_driftSin; }    ///< Periodic terms of the drift

        /// @brief Mean pedestal offset f(t0 + t); 0 when no drift equation is configured.
        /// @param t Time since the reset (s).
        /// @return Offset (e⁻).
        double pedestalMean(double t) const;

        // ── Noise generation ───────────────────────────────────────────────────

        /// @brief Noise map of a width × height matrix: one read, independent of any other call.
        /// @details The matrix is placed on the detector by the geometry. Channels are generated in
        ///          parallel (one std::async task per channel).
        /// @param width  Matrix columns.
        /// @param height Matrix rows.
        /// @param tStart Time since the reset at which this frame's read starts (s); it only matters
        ///               for the pedestal drift.
        /// @return width × height values, row-major (e⁻).
        /// @throws std::invalid_argument if the matrix does not fit on the detector or
        ///         colPerChannel does not divide the detector width.
        std::valarray<double> generate(size_t width, size_t height, double tStart = 0.0) const;

        /// @brief Noise map of the window [c0, c0+ncols) × [r0, r0+nrows) of a width × height matrix.
        /// @details Only the samples of the window are generated.
        /// @return ncols × nrows values, row-major (e⁻).
        /// @throws std::out_of_range if the window exceeds the matrix.
        std::valarray<double> generate(size_t width, size_t height,
                                       size_t c0, size_t r0,
                                       size_t ncols, size_t nrows,
                                       double tStart = 0.0) const;

    private:
        double m_amplitude    = 0.0;       ///< Per-channel 1/f noise: RMS (e⁻)
        double m_slope        = 1.0;       ///< Per-channel 1/f noise: PSD exponent
        size_t m_nChannels    = 0;         ///< 0 = detector width / m_colPerChannel
        size_t m_colPerChannel = 64;       ///< Columns per channel
        size_t m_nPoles       = 16;        ///< Per-channel 1/f noise: AR(1) poles
        double m_frameTime    = 1.45408;   ///< NISP HAWAII-2RG frame time (s)
        Geometry m_geometry;               ///< physical detector and placement of the matrix

        // bias_voltage_noise: 1/f series common to every channel
        double m_biasAmplitude = 0.0;      ///< Bias voltage noise: RMS (e⁻)
        double m_biasSlope     = 1.0;      ///< Bias voltage noise: PSD exponent
        size_t m_biasPoles     = 16;       ///< Bias voltage noise: AR(1) poles

        // pedestal_drift: f(t0 + t) + sigma · ε per channel per frame
        double                m_pedestalSigma = 0.0;   ///< RMS per channel and frame (e⁻)
        double                m_pedestalT0    = 0.0;   ///< Time added to the time since reset (s)
        std::vector<double>   m_driftPoly;             ///< Polynomial terms of the drift
        std::vector<DriftExp> m_driftExp;              ///< Exponential terms of the drift
        std::vector<DriftSin> m_driftSin;              ///< Periodic terms of the drift

        /// @brief A rectangle of physical pixels [px0, px0+ncols) × [py0, py0+nrows): the
        ///        matrix, or a window of it. The output map is ncols × nrows, row-major.
        struct Region { size_t px0, py0, ncols, nrows; };

        /// @brief The part of one output channel inside a region.
        struct Segment {
            size_t channel;   ///< physical channel index
            size_t outCol;    ///< first column, in region (output map) coordinates
            size_t width;     ///< number of columns
            size_t offset;    ///< pixel-clock index of the first column within a row
        };

        /// @brief AR(1) descriptor for one pole of the PSD approximation.
        struct Pole {
            double rho;   ///< AR(1) coefficient: rho = exp(−2π f_k)
            double sigma; ///< Input white-noise σ, normalized so Σ Var_k = 1
        };

        /// @brief Exact advance of one pole over n samples: x ← a·x + b·ε.
        struct Jump { double a; double b; };

        /// @brief One AR(1) bank following a series through the frames: `width` samples per
        ///        row over the `nrows` rows of a region, the other samples being skipped.
        struct Process {
            std::vector<Pole>   poles;       ///< empty: component off
            std::vector<double> state;       ///< current value of each pole
            std::vector<Jump>   rowJump;     ///< from a row's last sample to the next row's first
            std::vector<Jump>   frameJump;   ///< from the region's last row to its first row in the next frame
            size_t              rowGap   = 0;   ///< samples skipped between two rows
            size_t              frameGap = 0;   ///< samples skipped between two frames
            bool enabled() const { return !poles.empty(); }   ///< True if the component is on
        };

        /// @brief Channels of the detector overlapping `region`, from left to right.
        /// @throws std::invalid_argument if colPerChannel does not divide the detector width, if
        ///         nChannels does not match, or if the region is outside the detector.
        std::vector<Segment> _segments(const Region& region) const;

        /// @brief Pixel-clock span [offset, offset + width) covered by the segments: the part of
        ///        each row the bias voltage series must provide.
        static std::pair<size_t, size_t> _clockSpan(const std::vector<Segment>& segments);

        /// @brief Process of `width` samples per row over `nrows` rows, with the spectrum
        ///        (slope, nPoles) of a physical frame, started from its stationary state.
        /// @return The process, disabled (no poles) when amplitude <= 0.
        Process _process(double amplitude, double slope, size_t nPoles, size_t width, size_t nrows) const;

        /// @brief One per-channel 1/f process per segment.
        std::vector<Process> _channelProcesses(const std::vector<Segment>& segments, size_t nrows) const;

        /// @brief Build nPoles AR(1) descriptors for a series of `channelLength` samples and
        ///        PSD exponent `slope`. Corner frequencies log-spaced in [1/(2N), 0.5].
        static std::vector<Pole> _buildPoles(size_t channelLength, double slope, size_t nPoles);

        /// @brief Coefficients advancing every pole over n samples.
        static std::vector<Jump> _jumps(const std::vector<Pole>& poles, size_t n);

        /// @brief Advance a state by one jump (one draw per pole).
        static void _advance(std::vector<double>& state, const std::vector<Jump>& jump);

        /// @brief Fill out[0..length) with consecutive samples scaled by `amplitude`,
        ///        advancing `state` in place.
        /// @details Uses a thread_local mt19937_64 — safe to call from concurrent threads.
        static void _fillChannel(const std::vector<Pole>& poles, std::vector<double>& state,
                                 size_t length, double amplitude, double* out);

        /// @brief One frame of a process: nrows × width samples into out, skipping the samples
        ///        between rows, and between frames unless firstFrame.
        static void _run(Process& p, size_t nrows, size_t width, double amplitude,
                         bool firstFrame, double* out);

        /// @brief Pedestal mean of the rows [py0, py0 + nrows) of the frame read from tStart
        ///        (electrons, one value per row). Empty when no drift equation is configured.
        std::vector<double> _pedestalRows(size_t py0, size_t nrows, double tStart) const;

        /// @brief One frame of readout noise over `region`, written into out (ncols × nrows,
        ///        electrons).
        /// @details For every segment (one std::async task each):
        ///          per-channel 1/f (channels[i])
        ///          + bias voltage noise (bias, one series for all channels, pixel-clock span
        ///            [clockOffset, clockOffset + clockWidth), generated in its own task)
        ///          + pedestal: drift of the row + one N(0, sigma) draw for the channel.
        ///          The processes are advanced in place.
        void _frame(const Region& region, const std::vector<Segment>& segments,
                    std::vector<Process>& channels, Process& bias,
                    size_t clockOffset, size_t clockWidth,
                    bool firstFrame, double tStart, std::valarray<double>& out) const;

        /// @brief One independent read of a region (fresh, stationary processes).
        std::valarray<double> _read(const Region& region, double tStart) const;

    public:
        // ── Continuous multi-frame stream ───────────────────────────────────────

        /// @brief ONE continuous readout-noise process spanning several consecutive frame reads
        ///        (e.g. every read within a single ramp exposure).
        /// @details Each next() call advances every AR(1) state (per-channel
        ///          and bias voltage) by one physical frame instead of restarting —
        ///          matching how a real amplifier's 1/f noise persists across
        ///          non-destructive reads within one ramp. It also counts the frames, so
        ///          the pedestal drift is evaluated at the right time: call next() once
        ///          per frame, read or dropped.
        ///
        ///          Only call beginStream() when isEnabled() is true.
        class Stream
        {
        public:
            /// @brief Move constructor.
            Stream(Stream&&) noexcept = default;
            /// @brief Move assignment.
            Stream& operator=(Stream&&) noexcept = default;
            /// @brief Not copyable: a stream is one process.
            Stream(const Stream&) = delete;
            /// @brief Not copyable: a stream is one process.
            Stream& operator=(const Stream&) = delete;

            /// @brief Next frame (width × height, row-major, electrons), continuing the AR(1)
            ///        states from the previous call.
            std::valarray<double> next();

        private:
            friend class OneOverFNoise;
            /// @brief Stream of a width × height matrix, see beginStream().
            Stream(const OneOverFNoise& owner, size_t width, size_t height);

            const OneOverFNoise*  m_owner;             ///< noise model (not owned)
            Region                m_region;            ///< the matrix, in physical coordinates
            std::vector<Segment>  m_segments;          ///< channels of the matrix
            std::vector<Process>  m_channels;          ///< per-channel 1/f, one per segment
            Process               m_bias;              ///< bias voltage noise
            size_t                m_clockOffset = 0;   ///< first pixel-clock sample of a row
            size_t                m_clockWidth  = 0;   ///< pixel-clock samples of a row
            size_t                m_frame       = 0;   ///< frames consumed since the reset
        };

        /// @brief Begin a new continuous noise stream for a width × height matrix (e.g. one
        ///        ramp), starting at the reset (frame 0).
        /// @details The stream keeps a pointer to this object, which must outlive it.
        Stream beginStream(size_t width, size_t height) const;
    };

#pragma endregion

} // namespace nisp

#endif // _NISP_ONE_OVER_F_NOISE_H_
