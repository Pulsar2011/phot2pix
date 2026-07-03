//
// Created by gillard on 2026-06-26.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispOneOverFNoise.h
//
// 1/f (flicker / pink) noise model for the NISP HAWAII-2RG.
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
        // Default-constructed instance is disabled (amplitude == 0).
        OneOverFNoise() = default;

        // amplitude    : per-channel RMS in electrons (typical NISP: 5—10 e⁻)
        // slope        : PSD exponent α  (1.0 = pure 1/f,  2.0 = Brownian / random walk)
        // nChannels    : number of H2RG output channels; 0 = auto-derive from width/colPerChannel
        // colPerChannel: nominal channel width in pixels (64 for NISP 2048-wide SCA)
        // nPoles       : AR(1) poles for PSD approximation (16 → < 0.5 dB residual over 4 decades)
        explicit OneOverFNoise(double amplitude,
                               double slope        = 1.0,
                               size_t nChannels    = 0,
                               size_t colPerChannel = 64,
                               size_t nPoles       = 16);

        // Configure from a PixelConfig entry keyed "oneOverF".
        //
        // Scalar shorthand (amplitude only, all other parameters keep their defaults):
        //   "oneOverF": 7.0
        //
        // Full object form:
        //   "oneOverF": {
        //     "amplitude"    : 7.0,
        //     "slope"        : 1.0,
        //     "nChannels"    : 0,
        //     "colPerChannel": 64,
        //     "nPoles"       : 16
        //   }
        void initialize(const phot2pix::PixelConfig& cfg);

        // ── Accessors ──────────────────────────────────────────────────────────
        bool   isEnabled()     const { return m_amplitude > 0.0; }
        double getAmplitude()  const { return m_amplitude; }
        double getSlope()      const { return m_slope; }
        size_t getNChannels()  const { return m_nChannels; }    // 0 = auto at generate() time
        size_t getColPerChannel() const { return m_colPerChannel; }
        size_t getNPoles()     const { return m_nPoles; }

        // ── Noise generation ───────────────────────────────────────────────────

        // Full-frame noise map (width × height, row-major, electrons).
        // Channels are generated in parallel (one std::async task per channel).
        // width must satisfy the channel-layout rules described in the header comment.
        std::valarray<double> generate(size_t width, size_t height) const;

        // Windowed noise map (ncols × nrows, row-major, electrons).
        // Only channels that overlap [c0, r0, c0+ncols, r0+nrows) are generated.
        // Each overlapping channel's full time series is computed so that the
        // AR(1) state (and thus the spatial correlation) is correct relative to
        // the start of the frame readout.  Channels still run in parallel.
        std::valarray<double> generate(size_t width, size_t height,
                                       size_t c0, size_t r0,
                                       size_t ncols, size_t nrows) const;

    private:
        double m_amplitude    = 0.0;
        double m_slope        = 1.0;
        size_t m_nChannels    = 0;    // 0 = auto-derive from width / m_colPerChannel
        size_t m_colPerChannel = 64;
        size_t m_nPoles       = 16;

        // Descriptor for one output channel stripe.
        struct ChannelStripe {
            size_t col_start; // first detector column of this channel
            size_t col_width; // number of columns in this channel
        };

        // Build the channel layout for a detector of given width.
        // Detects the reference-pixel case automatically.
        std::vector<ChannelStripe> _buildLayout(size_t width) const;

        // AR(1) descriptor for one pole of the PSD approximation.
        struct Pole {
            double rho;   // AR(1) coefficient: rho = exp(−2π f_k)
            double sigma; // Input white-noise σ, normalized so Σ Var_k = 1
        };

        // Build nPoles AR(1) descriptors for a channel of `channelLength` samples.
        // Corner frequencies log-spaced in [1/(2N), 0.5].
        std::vector<Pole> _buildPoles(size_t channelLength) const;

        // Fill out[0..length) with one channel's 1/f noise realization,
        // advancing `state` in place. Pass a fresh zero-initialized vector for
        // an independent draw (generate()), or a vector that persists across
        // calls to continue the AR(1) recursion instead of restarting it
        // (Stream::next()). Output is scaled by m_amplitude.
        // Uses thread_local mt19937_64 — safe to call from concurrent threads.
        void _fillChannel(const std::vector<Pole>& poles, std::vector<double>& state,
                          size_t length, double* out) const;

    public:
        // ── Continuous multi-frame stream ───────────────────────────────────────
        //
        // Represents ONE continuous per-channel AR(1) process spanning several
        // consecutive frame reads (e.g. every read within a single ramp
        // exposure). Each next() call advances the AR(1) state by one frame's
        // worth of samples instead of restarting from zero — matching how a
        // real amplifier's 1/f noise persists across non-destructive reads
        // within one ramp, rather than resetting on every read as generate()
        // effectively does when called once per read.
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
            // continuing the AR(1) state from the previous call.
            std::valarray<double> next();

        private:
            friend class OneOverFNoise;
            Stream(const OneOverFNoise& owner, size_t width, size_t height);

            const OneOverFNoise*              m_owner;
            size_t                            m_width;
            size_t                            m_height;
            std::vector<ChannelStripe>        m_layout;
            std::vector<std::vector<Pole>>    m_poles;   // one pole set per channel, built once
            std::vector<std::vector<double>>  m_state;   // persistent AR(1) state per channel, carried across next() calls
        };

        // Begin a new continuous 1/f noise stream for a multi-read exposure
        // (e.g. one ramp). Poles are built once per channel; each
        // Stream::next() call advances the shared per-channel AR(1) state
        // rather than restarting it.
        Stream beginStream(size_t width, size_t height) const;
    };

#pragma endregion

} // namespace nisp

#endif // _NISP_ONE_OVER_F_NOISE_H_
