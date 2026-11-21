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
// The H2RG reads through nChannels parallel output amplifiers, each covering
// a vertical stripe of colPerChannel columns.  Within a stripe, pixels are
// read row by row, so the amplifier's 1/f noise drifts along the readout
// time axis and maps back as horizontal banding on the 2D frame.
//
// PSD approximation: sum of nPoles AR(1) processes with log-spaced corner
// frequencies — no FFT dependency, O(N × nPoles) per channel.
//
//   S(f) ≈ Σ_k  σ_k² / |1 − ρ_k exp(−2πif)|²  ≈  A² / f^slope
//
// Channel layout
// ─────────────
// The number of channels is derived at generate() time from the detector
// width and colPerChannel (default 64 for NISP).  Two special cases apply:
//
//   width == nCh × colPerChannel       uniform layout (e.g. 2048-wide full frame)
//   width == nCh × colPerChannel − 8   reference-pixel layout (e.g. 2040-wide
//                                       science frame): the first and last channels
//                                       are 4 columns narrower because the reference
//                                       pixel columns (4 on each side) are stripped.
//
// Any other width falls back to a uniform layout with nCh = round(width/colPerChannel).
//
// Time model
// ──────────
// All channels are read in parallel, in the same direction (increasing column
// within a row, then increasing row). Frame k after the reset (k = 0, 1, …, read or
// dropped) is read during [k·t_f, (k+1)·t_f], and row y of that frame at
//
//     t = k·t_f + (y + 0.5)/H · t_f        (t_f = frame time, H = detector height)
//
// The bias voltage series is indexed by the pixel clock within a row, so column x
// of channel c is sample (x − x_c + time offset): the time offset is 4 for the
// left edge channel of a 2040-wide frame (its 4 reference columns were read
// first) and 0 otherwise.
//
// The pedestal drift f(t) is evaluated once per row, at the row's read time.

#ifndef _NISP_ONE_OVER_F_NOISE_H_
#define _NISP_ONE_OVER_F_NOISE_H_

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
        // nChannels    : number of H2RG output channels; 0 = auto-derive from width/colPerChannel
        // colPerChannel: nominal channel width in pixels (64 for NISP 2048-wide SCA)
        // nPoles       : AR(1) poles for PSD approximation (16 → < 0.5 dB residual over 4 decades)
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
        // Throws std::invalid_argument on a malformed or out-of-range value.
        void initialize(const phot2pix::PixelConfig& cfg);

        // Frame time t_f in seconds (default: NISP 1.45408 s). nisp::SCA sets it from
        // its own "frame_time". Throws std::invalid_argument if not > 0.
        void setFrameTime(double frameTime);

        // ── Accessors ──────────────────────────────────────────────────────────
        // True when at least one component produces noise.
        bool   isEnabled()     const
        {
            return m_amplitude > 0.0 || m_biasAmplitude > 0.0 || m_pedestalSigma > 0.0 || hasPedestalDrift();
        }
        double getAmplitude()  const { return m_amplitude; }
        double getSlope()      const { return m_slope; }
        size_t getNChannels()  const { return m_nChannels; }    // 0 = auto at generate() time
        size_t getColPerChannel() const { return m_colPerChannel; }
        size_t getNPoles()     const { return m_nPoles; }
        double getFrameTime()  const { return m_frameTime; }

        double getBiasVoltageNoise() const { return m_biasAmplitude; }
        double getBiasSlope()        const { return m_biasSlope; }
        size_t getBiasNPoles()       const { return m_biasPoles; }

        double getPedestalSigma()    const { return m_pedestalSigma; }
        double getPedestalT0()       const { return m_pedestalT0; }
        bool   hasPedestalDrift()    const { return !m_driftPoly.empty() || !m_driftExp.empty() || !m_driftSin.empty(); }

        // Mean pedestal offset f(t0 + t) in electrons, t = time since the reset (s).
        // 0 when no drift equation is configured.
        double pedestalMean(double t) const;

        // ── Noise generation ───────────────────────────────────────────────────

        // Full-frame noise map (width × height, row-major, electrons): one read,
        // independent of any other call (every AR(1) state starts at zero).
        // tStart is the time since the reset at which this frame's read starts
        // (seconds); it only matters for the pedestal drift.
        // Channels are generated in parallel (one std::async task per channel).
        // width must satisfy the channel-layout rules described in the header comment.
        std::valarray<double> generate(size_t width, size_t height, double tStart = 0.0) const;

        // Windowed noise map (ncols × nrows, row-major, electrons).
        // Only channels that overlap [c0, r0, c0+ncols, r0+nrows) are generated.
        // Each overlapping channel's full time series is computed so that the
        // AR(1) state (and thus the spatial correlation) is correct relative to
        // the start of the frame readout.  Channels still run in parallel.
        std::valarray<double> generate(size_t width, size_t height,
                                       size_t c0, size_t r0,
                                       size_t ncols, size_t nrows,
                                       double tStart = 0.0) const;

    private:
        double m_amplitude    = 0.0;
        double m_slope        = 1.0;
        size_t m_nChannels    = 0;    // 0 = auto-derive from width / m_colPerChannel
        size_t m_colPerChannel = 64;
        size_t m_nPoles       = 16;
        double m_frameTime    = 1.45408;   // NISP HAWAII-2RG frame time (s)

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

        // Descriptor for one output channel stripe.
        struct ChannelStripe {
            size_t col_start;   // first detector column of this channel
            size_t col_width;   // number of columns in this channel
            size_t time_offset; // pixel-clock index of col_start within a row (see "Time model")
        };

        // Build the channel layout for a detector of given width.
        // Detects the reference-pixel case automatically.
        std::vector<ChannelStripe> _buildLayout(size_t width) const;

        // Pixel-clock samples per row: length of one row of the bias voltage series.
        static size_t _rowSamples(const std::vector<ChannelStripe>& layout);

        // AR(1) descriptor for one pole of the PSD approximation.
        struct Pole {
            double rho;   // AR(1) coefficient: rho = exp(−2π f_k)
            double sigma; // Input white-noise σ, normalized so Σ Var_k = 1
        };

        // Build nPoles AR(1) descriptors for a series of `channelLength` samples and
        // PSD exponent `slope`. Corner frequencies log-spaced in [1/(2N), 0.5].
        static std::vector<Pole> _buildPoles(size_t channelLength, double slope, size_t nPoles);

        // Fill out[0..length) with one 1/f noise realization scaled by `amplitude`,
        // advancing `state` in place. Pass a fresh zero-initialized vector for
        // an independent draw (generate()), or a vector that persists across
        // calls to continue the AR(1) recursion instead of restarting it
        // (Stream::next()).
        // Uses a thread_local mt19937_64 — safe to call from concurrent threads.
        static void _fillChannel(const std::vector<Pole>& poles, std::vector<double>& state,
                                 size_t length, double amplitude, double* out);

        // Poles and zero (cold-start) states of the per-channel 1/f noise for the
        // channels listed in `channels` (vectors indexed by channel, sized to the
        // layout; empty entries when the per-channel 1/f noise is off).
        void _initChannels(const std::vector<ChannelStripe>& layout,
                           const std::vector<size_t>& channels, size_t height,
                           std::vector<std::vector<Pole>>& poles,
                           std::vector<std::vector<double>>& states) const;

        // Poles and zero state of the bias voltage series (empty when it is off).
        void _initCommon(const std::vector<ChannelStripe>& layout, size_t height,
                         std::vector<Pole>& poles, std::vector<double>& state) const;

        // Pedestal mean of every row of the frame read from tStart (electrons, one
        // value per row). Empty when no drift equation is configured.
        std::vector<double> _pedestalRows(size_t height, double tStart) const;

        // One frame of readout noise, shared by generate() and Stream::next().
        // For every channel listed in `channels` (one std::async task each), builds
        // the channel's time series (height × col_width samples, electrons):
        //   per-channel 1/f (chanPoles[ch], chanStates[ch])
        //   + bias voltage noise (commonPoles, commonState; one series for all channels)
        //   + pedestal: drift of the row + one N(0, sigma) draw for the channel
        // then calls scatter(stripe, series) to place it in the output map.
        // The states are advanced in place. scatter is called concurrently, once per
        // channel, and must only write that channel's columns.
        template <typename Scatter>
        void _frame(const std::vector<ChannelStripe>& layout,
                    const std::vector<size_t>& channels,
                    const std::vector<std::vector<Pole>>& chanPoles,
                    std::vector<std::vector<double>>& chanStates,
                    const std::vector<Pole>& commonPoles,
                    std::vector<double>& commonState,
                    size_t height, double tStart,
                    const Scatter& scatter) const;

    public:
        // ── Continuous multi-frame stream ───────────────────────────────────────
        //
        // Represents ONE continuous readout-noise process spanning several
        // consecutive frame reads (e.g. every read within a single ramp
        // exposure). Each next() call advances every AR(1) state (per-channel
        // and bias voltage) by one frame's worth of samples instead of
        // restarting from zero — matching how a real amplifier's 1/f noise
        // persists across non-destructive reads within one ramp, rather than
        // resetting on every read as generate() effectively does when called
        // once per read. It also counts the frames, so the pedestal drift is
        // evaluated at the right time: call next() once per frame, read or dropped.
        //
        // Only call beginStream() when isEnabled() is true.
        class Stream
        {
        public:
            Stream(Stream&&) noexcept = default;
            Stream& operator=(Stream&&) noexcept = default;
            Stream(const Stream&) = delete;
            Stream& operator=(const Stream&) = delete;

            // Next full-frame chunk (width × height, row-major, electrons),
            // continuing the AR(1) states from the previous call.
            std::valarray<double> next();

        private:
            friend class OneOverFNoise;
            Stream(const OneOverFNoise& owner, size_t width, size_t height);

            const OneOverFNoise*              m_owner;
            size_t                            m_width;
            size_t                            m_height;
            size_t                            m_frame = 0;   // frames consumed since the reset
            std::vector<ChannelStripe>        m_layout;
            std::vector<size_t>               m_channels;    // 0 … nCh−1
            std::vector<std::vector<Pole>>    m_poles;       // one pole set per channel, built once
            std::vector<std::vector<double>>  m_state;       // persistent AR(1) state per channel, carried across next() calls
            std::vector<Pole>                 m_commonPoles; // bias voltage noise poles
            std::vector<double>               m_commonState; // bias voltage noise AR(1) state
        };

        // Begin a new continuous noise stream for a multi-read exposure (e.g. one
        // ramp), starting at the reset (frame 0). Poles are built once per channel;
        // each Stream::next() call advances the persistent AR(1) states rather
        // than restarting them.
        Stream beginStream(size_t width, size_t height) const;
    };

#pragma endregion

} // namespace nisp

#endif // _NISP_ONE_OVER_F_NOISE_H_
