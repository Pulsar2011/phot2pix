//
// Created by gillard on 2026-11-27.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispRTN.h
//
/// @file nispRTN.h
/// @brief Random telegraph noise of the NISP HAWAII-2RG pixels.
///
/// A trap in the source follower of a pixel captures and releases charges or fluctuation 
/// of the threshold voltage: the output of the pixel jumps therefore between two levels.
/// It translate into a voltage offset of the output (ADU) at  read time and leaves no persistence.
/// 
///
/// @par Model of one RTN pixel
/// Two levels: rest (offset 0) and shifted (offset s·A, A > 0 in ADU, s = −1 or +1). A is the
/// amplitude of this pixel. τ_0 and τ_1 are the mean times spent in the rest and shifted levels.
/// Chao CY, Yeh SF, Wu MH, et al. Sensors (Basel). 2019 Dec 10;19(24):5447.
///
/// @par Time evolution
/// The level follows a two-state Markov process. Over a time dt, with
/// @verbatim
///   λ = 1/τ_0 + 1/τ_1,   p = τ_1 / (τ_0 + τ_1),   e = exp(−λ·dt)
///   P(shifted | was shifted)  = p + (1 − p)·e
///   P(shifted | was at rest)  = p·(1 − e)
/// @endverbatim
/// p is the fraction of the time spent in the shifted
/// level: the state of a pixel at the reset is drawn from it, the trap does not reset with the pixel.
///
/// @par Drawing the pixels
/// initialize() draws, once, round(fraction · W·H) distinct pixels of the matrix, uniformly
/// (reference pixels included). Each gets τ_0 and τ_1, independently uniform in
/// [tauMin, tauMax], its amplitude A, uniform in [1 ADU, amplitude], and s = +1 with the
/// probability positiveFraction.
/// nisp::SCA then begins one Stream per exposure, advanced at every frame, read or dropped.

#ifndef _NISP_RTN_H_
#define _NISP_RTN_H_

#include <phot2pix/PixelConfig.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <random>
#include <valarray>
#include <vector>

namespace nisp
{

// ── RTN ───────────────────────────────────────────────────────────────────────
#pragma region - RTN

    /// @brief Random telegraph noise of the pixels: which pixels, their time constants and the sign
    ///        of their offset (see nispRTN.h).
    class RTN
    {
    public:
        /// @brief One pixel with RTN.
        struct Pixel
        {
            size_t idx;    ///< Matrix pixel, row-major: row · width + column
            double tau0;   ///< Mean time in the rest level (s)
            double tau1;        ///< Mean time in the shifted level (s)
            double amplitude;   ///< Amplitude A of the offset s·A (ADU), in [1, "amplitude"]
            int    sign;        ///< Sign s of the offset s·A: −1 (down) or +1 (up)

            /// @brief Fraction of the time spent in the shifted level, τ_1 / (τ_0 + τ_1).
            double shiftedProbability() const { return tau1 / (tau0 + tau1); }
        };

        /// @brief No RTN.
        RTN() = default;

        /// @brief Configure from "rtn" and draw the RTN pixels of a width × height matrix.
        /// @details
        /// @verbatim
        ///   "rtn": {
        ///     "fraction"        : 1e-3,    // fraction of the matrix pixels with RTN, in [0, 1]
        ///     "tauMin"          : 0.5,     // s, > 0
        ///     "tauMax"          : 100.0,   // s, >= tauMin
        ///     "amplitude"       : 30.0,    // ADU, >= 1: A of each RTN pixel uniform in [1, amplitude]
        ///     "positiveFraction": 0.4,     // optional (0.4): fraction of the RTN pixels with +A
        ///     "seed"            : 1        // optional: integer >= 0, reproducible draw
        ///   }
        /// @endverbatim
        /// "fraction", "tauMin", "tauMax" and "amplitude" are required; unknown keys are rejected.
        /// Every "rtn" replaces the whole previous configuration and draws the pixels again.
        /// Nothing changes without "rtn". fraction 0: no RTN.
        /// @param cfg    Configuration.
        /// @param width  Matrix columns.
        /// @param height Matrix rows.
        /// @throws std::invalid_argument if a value is missing, malformed or out of range; nothing
        ///         is then changed.
        void initialize(const phot2pix::PixelConfig& cfg, size_t width, size_t height);

        // ── Accessors ──────────────────────────────────────────────────────────

        /// @brief True when some pixels have RTN, with an amplitude > 0.
        bool isEnabled() const { return m_amplitude > 0.0 && !m_pixels.empty(); }

        double getFraction()         const { return m_fraction; }           ///< Fraction of the pixels with RTN
        double getTauMin()           const { return m_tauMin; }             ///< Shortest time constant (s)
        double getTauMax()           const { return m_tauMax; }             ///< Longest time constant (s)
        double getAmplitude()        const { return m_amplitude; }          ///< Largest amplitude A of a pixel (ADU)
        double getPositiveFraction() const { return m_positiveFraction; }   ///< Fraction of the RTN pixels with +A
        std::optional<uint64_t> getSeed() const { return m_seed; }          ///< Seed of the draw; empty: none
        size_t getWidth()            const { return m_width; }              ///< Matrix columns
        size_t getHeight()           const { return m_height; }             ///< Matrix rows

        /// @brief The RTN pixels, in increasing idx.
        const std::vector<Pixel>& getPixels() const { return m_pixels; }

        // ── Stream ─────────────────────────────────────────────────────────────

        /// @brief Levels of the RTN pixels during one exposure.
        /// @details Begun at the reset, in the stationary state; advance() moves it in time, call it
        ///          for every frame, read or dropped. Only begin it when isEnabled() is true.
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

            /// @brief Move the levels dt seconds forward (see nispRTN.h).
            /// @throws std::invalid_argument if dt is not >= 0.
            void advance(double dt);

            /// @brief Add the current offset s·A (ADU) of the shifted pixels to a matrix-sized map.
            /// @param values width × height values, row-major.
            /// @throws std::invalid_argument if values is smaller than the matrix.
            void addTo(std::valarray<double>& values) const;

            /// @brief Same for the window [c0, c0+ncols) × [r0, r0+nrows) of the matrix.
            /// @param values ncols × nrows values, row-major.
            /// @throws std::out_of_range if the window exceeds the matrix, std::invalid_argument if
            ///         values is smaller than the window.
            void addTo(std::valarray<double>& values, size_t c0, size_t r0, size_t ncols, size_t nrows) const;

            /// @brief Write the current offset of every RTN pixel (ADU, 0 at rest) in a matrix-sized
            ///        map; the other pixels are not touched.
            /// @throws std::invalid_argument if map does not have one value per pixel of the matrix.
            void writeTo(std::valarray<double>& map) const;

            /// @brief Level of each RTN pixel, in the order of getPixels(): 1 shifted, 0 at rest.
            const std::vector<uint8_t>& states() const { return m_shifted; }

        private:
            friend class RTN;
            /// @brief Stream of owner, see beginStream().
            explicit Stream(const RTN& owner);

            const RTN*           m_owner;     ///< RTN model (not owned)
            std::mt19937_64      m_rng;       ///< Generator of this stream
            std::vector<uint8_t> m_shifted;   ///< Level of each RTN pixel
        };

        /// @brief Begin the stream of one exposure, at the reset.
        /// @details The stream keeps a pointer to this object, which must outlive it.
        Stream beginStream() const { return Stream(*this); }

    private:
        /// @brief Generator seeded by seed, or by the random device without seed.
        static std::mt19937_64 makeRng(const std::optional<uint64_t>& seed);

        double                  m_fraction         = 0.0;   ///< "fraction"
        double                  m_tauMin           = 0.0;   ///< "tauMin" (s)
        double                  m_tauMax           = 0.0;   ///< "tauMax" (s)
        double                  m_amplitude        = 0.0;   ///< "amplitude": largest amplitude (ADU)
        double                  m_positiveFraction = 0.4;   ///< "positiveFraction"
        std::optional<uint64_t> m_seed;                     ///< "seed"; empty: none
        size_t                  m_width            = 0;     ///< Matrix columns
        size_t                  m_height           = 0;     ///< Matrix rows
        std::vector<Pixel>      m_pixels;                   ///< The RTN pixels, in increasing idx
    };

#pragma endregion

} // namespace nisp

#endif // _NISP_RTN_H_
