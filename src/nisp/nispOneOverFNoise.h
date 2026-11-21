//
// Created by gillard on 2026-06-26.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispOneOverFNoise.h
//
// Readout noise of the NISP HAWAII-2RG output channels, added at each read.
// Three components, each enabled by its own configuration key:
//
//   oneOverF            1/f (flicker / pink) noise of each channel's own amplifier.
//                       Independent between channels.
//   bias_voltage_noise  1/f noise of the bias voltages shared by every channel: ONE
//                       time series, seen by all channels at the same pixel clock.
//                       This is what the left/right reference columns track.
//   pedestal_drift      offset of each channel: a mean drift f(t), common to all
//                       channels and given as an equation, plus white noise drawn
//                       per channel and per frame. This is what the top/bottom
//                       reference rows track.
//
// All amplitudes are in electrons: the generated maps are in electrons and
// nisp::SCA converts them to ADU with each pixel's gain, as it does for read noise.
//
// Physical detector
// ─────────────────
// The noise is defined on the physical detector (nisp::Geometry, 2048 × 2048 by
// default), not on the simulated matrix. The detector is read through
// width / colPerChannel output channels (32 × 64 columns for NISP); channel c covers
// the physical columns [c·w, (c+1)·w), w = colPerChannel. A simulated matrix is
// placed on the detector (at (0, 0) by default), and each of its pixels takes the noise
// of the physical pixel it covers. The channel boundaries, the readout timing and the
// spectrum of the noise are therefore the same whatever the size of the matrix: a
// smaller matrix is only a smaller part of the same detector.
//
// Time model
// ──────────
// All channels are read in parallel, in the same direction (increasing column
// within a row, then increasing row), with no row or frame overhead. Within one
// frame, channel c reads physical pixel (X, Y) at pixel-clock sample
//
//     t = Y·w + (X − c·w)            (one frame = L = height · w samples)
//
// and frame k after the reset (k = 0, 1, …, read or dropped) starts at k·L. In
// seconds, row Y of frame k is read at  k·t_f + (Y + 0.5)/height · t_f.
//
// 1/f noise: PSD approximation by a sum of nPoles AR(1) processes with corner
// frequencies log-spaced in [1/(2L), 0.5] cycles per sample — no FFT dependency.
// Each pole carries a variance ∝ f_k^(1 − slope), so that equal log steps give
//
//   S(f) ≈ A² / f^slope          (slope 1 = pink, 2 = Brownian / random walk)
//
// Only the samples inside the matrix (or window) are generated: the samples of the
// rows and columns outside it are skipped exactly, an AR(1) process being advanced
// over n samples in one step,
//
//   x ← ρⁿ x + σ √((1 − ρ²ⁿ)/(1 − ρ²)) ε,
//
// so the cost scales with the matrix area. Every process starts from its stationary
// state: the amplifier noise does not reset with the pixels.
//
// The bias voltage series is indexed by the pixel clock (t above), so at a given
// sample every channel gets the same value. The pedestal drift f(t) is evaluated once
// per row, at the row's read time.

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

    class OneOverFNoise
    {
    public:
        /// One exponential term of the pedestal drift: amplitude · exp(−t/tau).
        struct DriftExp { double amplitude; double tau; };
        /// One periodic term of the pedestal drift: amplitude · sin(2π t/period + phase).
        struct DriftSin { double amplitude; double period; double phase; };

        // Default-constructed instance is disabled (every component off).
        OneOverFNoise() = default;

        // amplitude    : per-channel RMS in electrons (typical NISP: 5—10 e⁻)
        // slope        : PSD exponent α  (1.0 = pure 1/f,  2.0 = Brownian / random walk)
        // nChannels    : number of H2RG output channels; 0 = detector width / colPerChannel.
        //                When given, it must equal detector width / colPerChannel.
        // colPerChannel: channel width in pixels (64 for NISP); must divide the detector width
        // nPoles       : AR(1) poles for PSD approximation (about 3 per decade: 16 covers the 5 decades of a frame)
        // Only the per-channel 1/f noise is set here; bias_voltage_noise and
        // pedestal_drift are configured through initialize().
        explicit OneOverFNoise(double amplitude,
                               double slope        = 1.0,
                               size_t nChannels    = 0,
                               size_t colPerChannel = 64,
                               size_t nPoles       = 16);

        // Configure from a PixelConfig. Each key is optional; absent keys leave
        // their component unchanged.
        //
        // Per-channel 1/f noise, scalar shorthand (amplitude only) or object:
        //   "oneOverF": 7.0
        //   "oneOverF": {
        //     "amplitude"    : 7.0,
        //     "slope"        : 1.0,
        //     "nChannels"    : 0,
        //     "colPerChannel": 64,
        //     "nPoles"       : 16
        //   }
        //
        // Bias voltage noise (1/f common to all channels), scalar shorthand or object:
        //   "bias_voltage_noise": 3.0
        //   "bias_voltage_noise": { "amplitude": 3.0, "slope": 1.0, "nPoles": 16 }
        //
        // Pedestal drift, scalar shorthand (white noise only) or object:
        //   "pedestal_drift": 4.0
        //   "pedestal_drift": {
        //     "sigma": 4.0,                  // RMS around the mean, per channel per frame (e⁻)
        //     "t0"   : 0.0,                  // time added to the time since reset (s)
        //     "drift": {                     // mean offset f(t) in e⁻, every term optional
        //       "poly": [a0, a1, a2],                                    // Σ a_k t^k
        //       "exp" : [{"amplitude": A, "tau": τ}],                    // Σ A exp(−t/τ)
        //       "sin" : [{"amplitude": B, "period": P, "phase": φ}]      // Σ B sin(2π t/P + φ), φ in rad
        //     }
        //   }
        // "exp" and "sin" also accept a single object instead of an array.
        //
        // The detector geometry is not read here: nisp::SCA passes it with setGeometry().
        // Throws std::invalid_argument on a malformed or out-of-range value.
        void initialize(const phot2pix::PixelConfig& cfg);

        // Frame time t_f in seconds (default: NISP 1.45408 s). nisp::SCA sets it from
        // its own "frame_time". Throws std::invalid_argument if not > 0.
        void setFrameTime(double frameTime);

        // Physical detector and placement of the matrix (default: 2048 × 2048, origin (0, 0)).
        // nisp::SCA sets it from its own configuration.
        void setGeometry(const Geometry& geometry) { m_geometry = geometry; }

        // Check that a width × height matrix fits on the detector and that the channel
        // layout matches it (colPerChannel divides the detector width, nChannels agrees).
        // Throws std::invalid_argument otherwise; nisp::SCA calls it from initialize().
        void checkGeometry(size_t width, size_t height) const;

        // ── Accessors ──────────────────────────────────────────────────────────
        // True when at least one component produces noise.
        bool   isEnabled()     const
        {
            return m_amplitude > 0.0 || m_biasAmplitude > 0.0 || m_pedestalSigma > 0.0 || hasPedestalDrift();
        }
        double getAmplitude()  const { return m_amplitude; }
        double getSlope()      const { return m_slope; }
        size_t getNChannels()  const { return m_nChannels; }    // 0 = detector width / colPerChannel
        size_t getColPerChannel() const { return m_colPerChannel; }
        size_t getNPoles()     const { return m_nPoles; }
        double getFrameTime()  const { return m_frameTime; }
        const Geometry& getGeometry() const { return m_geometry; }

        double getBiasVoltageNoise() const { return m_biasAmplitude; }
        double getBiasSlope()        const { return m_biasSlope; }
        size_t getBiasNPoles()       const { return m_biasPoles; }

        double getPedestalSigma()    const { return m_pedestalSigma; }
        double getPedestalT0()       const { return m_pedestalT0; }
        bool   hasPedestalDrift()    const { return !m_driftPoly.empty() || !m_driftExp.empty() || !m_driftSin.empty(); }
        const std::vector<double>&   getDriftPoly() const { return m_driftPoly; }   ///< [a0, a1, …], e⁻ / e⁻ s⁻¹ / …
        const std::vector<DriftExp>& getDriftExp()  const { return m_driftExp; }
        const std::vector<DriftSin>& getDriftSin()  const { return m_driftSin; }

        // Mean pedestal offset f(t0 + t) in electrons, t = time since the reset (s).
        // 0 when no drift equation is configured.
        double pedestalMean(double t) const;

        // ── Noise generation ───────────────────────────────────────────────────

        // Noise map of a width × height matrix (row-major, electrons), placed on the
        // detector by the geometry: one read, independent of any other call.
        // tStart is the time since the reset at which this frame's read starts
        // (seconds); it only matters for the pedestal drift.
        // Channels are generated in parallel (one std::async task per channel).
        // Throws std::invalid_argument if the matrix does not fit on the detector or
        // colPerChannel does not divide the detector width.
        std::valarray<double> generate(size_t width, size_t height, double tStart = 0.0) const;

        // Windowed noise map (ncols × nrows, row-major, electrons) of the window
        // [c0, c0+ncols) × [r0, r0+nrows) of a width × height matrix. Only the
        // samples of the window are generated. Throws std::out_of_range if the window
        // exceeds the matrix.
        std::valarray<double> generate(size_t width, size_t height,
                                       size_t c0, size_t r0,
                                       size_t ncols, size_t nrows,
                                       double tStart = 0.0) const;

    private:
        double m_amplitude    = 0.0;
        double m_slope        = 1.0;
        size_t m_nChannels    = 0;    // 0 = detector width / m_colPerChannel
        size_t m_colPerChannel = 64;
        size_t m_nPoles       = 16;
        double m_frameTime    = 1.45408;   // NISP HAWAII-2RG frame time (s)
        Geometry m_geometry;               // physical detector and placement of the matrix

        // bias_voltage_noise: 1/f series common to every channel
        double m_biasAmplitude = 0.0;
        double m_biasSlope     = 1.0;
        size_t m_biasPoles     = 16;

        // pedestal_drift: f(t0 + t) + sigma · ε per channel per frame
        double                m_pedestalSigma = 0.0;
        double                m_pedestalT0    = 0.0;
        std::vector<double>   m_driftPoly;
        std::vector<DriftExp> m_driftExp;
        std::vector<DriftSin> m_driftSin;

        // A rectangle of physical pixels [px0, px0+ncols) × [py0, py0+nrows): the
        // matrix, or a window of it. The output map is ncols × nrows, row-major.
        struct Region { size_t px0, py0, ncols, nrows; };

        // The part of one output channel inside a region.
        struct Segment {
            size_t channel;   // physical channel index
            size_t outCol;    // first column, in region (output map) coordinates
            size_t width;     // number of columns
            size_t offset;    // pixel-clock index of the first column within a row
        };

        // AR(1) descriptor for one pole of the PSD approximation.
        struct Pole {
            double rho;   // AR(1) coefficient: rho = exp(−2π f_k)
            double sigma; // Input white-noise σ, normalized so Σ Var_k = 1
        };

        // Exact advance of one pole over n samples: x ← a·x + b·ε.
        struct Jump { double a; double b; };

        // One AR(1) bank following a series through the frames: `width` samples per
        // row over the `nrows` rows of a region, the other samples being skipped.
        struct Process {
            std::vector<Pole>   poles;       // empty: component off
            std::vector<double> state;
            std::vector<Jump>   rowJump;     // from a row's last sample to the next row's first
            std::vector<Jump>   frameJump;   // from the region's last row to its first row in the next frame
            size_t              rowGap   = 0;   // samples skipped between two rows
            size_t              frameGap = 0;   // samples skipped between two frames
            bool enabled() const { return !poles.empty(); }
        };

        // Channels of the detector overlapping `region`, from left to right. Throws
        // std::invalid_argument if colPerChannel does not divide the detector width, if
        // nChannels does not match, or if the region is outside the detector.
        std::vector<Segment> _segments(const Region& region) const;

        // Pixel-clock span [offset, offset + width) covered by the segments: the part of
        // each row the bias voltage series must provide.
        static std::pair<size_t, size_t> _clockSpan(const std::vector<Segment>& segments);

        // Process of `width` samples per row over `nrows` rows, with the spectrum
        // (slope, nPoles) of a physical frame, started from its stationary state.
        // Disabled (no poles) when amplitude <= 0.
        Process _process(double amplitude, double slope, size_t nPoles, size_t width, size_t nrows) const;

        // One per-channel 1/f process per segment.
        std::vector<Process> _channelProcesses(const std::vector<Segment>& segments, size_t nrows) const;

        // Build nPoles AR(1) descriptors for a series of `channelLength` samples and
        // PSD exponent `slope`. Corner frequencies log-spaced in [1/(2N), 0.5].
        static std::vector<Pole> _buildPoles(size_t channelLength, double slope, size_t nPoles);

        // Coefficients advancing every pole over n samples.
        static std::vector<Jump> _jumps(const std::vector<Pole>& poles, size_t n);

        // Advance a state by one jump (one draw per pole).
        static void _advance(std::vector<double>& state, const std::vector<Jump>& jump);

        // Fill out[0..length) with consecutive samples scaled by `amplitude`,
        // advancing `state` in place.
        // Uses a thread_local mt19937_64 — safe to call from concurrent threads.
        static void _fillChannel(const std::vector<Pole>& poles, std::vector<double>& state,
                                 size_t length, double amplitude, double* out);

        // One frame of a process: nrows × width samples into out, skipping the samples
        // between rows, and between frames unless firstFrame.
        static void _run(Process& p, size_t nrows, size_t width, double amplitude,
                         bool firstFrame, double* out);

        // Pedestal mean of the rows [py0, py0 + nrows) of the frame read from tStart
        // (electrons, one value per row). Empty when no drift equation is configured.
        std::vector<double> _pedestalRows(size_t py0, size_t nrows, double tStart) const;

        // One frame of readout noise over `region`, written into out (ncols × nrows,
        // electrons). For every segment (one std::async task each):
        //   per-channel 1/f (channels[i])
        //   + bias voltage noise (bias, one series for all channels, pixel-clock span
        //     [clockOffset, clockOffset + clockWidth), generated in its own task)
        //   + pedestal: drift of the row + one N(0, sigma) draw for the channel.
        // The processes are advanced in place.
        void _frame(const Region& region, const std::vector<Segment>& segments,
                    std::vector<Process>& channels, Process& bias,
                    size_t clockOffset, size_t clockWidth,
                    bool firstFrame, double tStart, std::valarray<double>& out) const;

        // One independent read of a region (fresh, stationary processes).
        std::valarray<double> _read(const Region& region, double tStart) const;

    public:
        // ── Continuous multi-frame stream ───────────────────────────────────────
        //
        // Represents ONE continuous readout-noise process spanning several
        // consecutive frame reads (e.g. every read within a single ramp
        // exposure). Each next() call advances every AR(1) state (per-channel
        // and bias voltage) by one physical frame instead of restarting —
        // matching how a real amplifier's 1/f noise persists across
        // non-destructive reads within one ramp. It also counts the frames, so
        // the pedestal drift is evaluated at the right time: call next() once
        // per frame, read or dropped.
        //
        // Only call beginStream() when isEnabled() is true.
        class Stream
        {
        public:
            Stream(Stream&&) noexcept = default;
            Stream& operator=(Stream&&) noexcept = default;
            Stream(const Stream&) = delete;
            Stream& operator=(const Stream&) = delete;

            // Next frame (width × height, row-major, electrons), continuing the AR(1)
            // states from the previous call.
            std::valarray<double> next();

        private:
            friend class OneOverFNoise;
            Stream(const OneOverFNoise& owner, size_t width, size_t height);

            const OneOverFNoise*  m_owner;
            Region                m_region;        // the matrix, in physical coordinates
            std::vector<Segment>  m_segments;
            std::vector<Process>  m_channels;      // per-channel 1/f, one per segment
            Process               m_bias;          // bias voltage noise
            size_t                m_clockOffset = 0;
            size_t                m_clockWidth  = 0;
            size_t                m_frame       = 0;   // frames consumed since the reset
        };

        // Begin a new continuous noise stream for a width × height matrix (e.g. one
        // ramp), starting at the reset (frame 0).
        Stream beginStream(size_t width, size_t height) const;
    };

#pragma endregion

} // namespace nisp

#endif // _NISP_ONE_OVER_F_NOISE_H_
