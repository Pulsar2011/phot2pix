//
// Created by gillard on 2026-06-26.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispOneOverFNoise.cxx
//

#include "nispOneOverFNoise.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <future>
#include <numbers>
#include <random>
#include <stdexcept>
#include <thread>

namespace nisp
{

// ── OneOverFNoise ──────────────────────────────────────────────────────────────
#pragma region - OneOverFNoise

    OneOverFNoise::OneOverFNoise(double amplitude, double slope,
                                 size_t nChannels, size_t colPerChannel, size_t nPoles)
        : m_amplitude(amplitude), m_slope(slope),
          m_nChannels(nChannels), m_colPerChannel(colPerChannel), m_nPoles(nPoles)
    {
        if (m_colPerChannel < 1)
            throw std::invalid_argument("OneOverFNoise: colPerChannel must be >= 1.");
        if (m_nPoles < 1)
            throw std::invalid_argument("OneOverFNoise: nPoles must be >= 1.");
    }

    void OneOverFNoise::initialize(const phot2pix::PixelConfig& cfg)
    {
        if (!cfg.has("oneOverF"))
            return;

        // Scalar shorthand: "oneOverF": 7.0 → amplitude only, keep defaults
        try
        {
            m_amplitude = cfg.getAsDouble("oneOverF");
            return;
        }
        catch (const std::invalid_argument&) {}      // not a scalar — fall through to object form
        catch (const nlohmann::json::type_error&) {} // stored as json object — fall through to object form

        const auto& j = cfg.getJson("oneOverF");

        if (j.contains("amplitude"))     m_amplitude    = j.at("amplitude").get<double>();
        if (j.contains("slope"))         m_slope        = j.at("slope").get<double>();
        if (j.contains("nChannels"))     m_nChannels    = j.at("nChannels").get<size_t>();
        if (j.contains("colPerChannel")) m_colPerChannel = j.at("colPerChannel").get<size_t>();
        if (j.contains("nPoles"))        m_nPoles       = j.at("nPoles").get<size_t>();

        if (m_colPerChannel < 1)
            throw std::invalid_argument("OneOverFNoise: colPerChannel must be >= 1.");
        if (m_nPoles < 1)
            throw std::invalid_argument("OneOverFNoise: nPoles must be >= 1.");
    }

// ── _buildLayout ──────────────────────────────────────────────────────────────
//
// Rules (in priority order):
//
//  1. nCh = m_nChannels if set explicitly (> 0), otherwise
//     nCh = round(width / m_colPerChannel)
//
//  2. If nCh × m_colPerChannel == width          → uniform layout
//     If nCh × m_colPerChannel == width + 8      → reference-pixel layout:
//                                                   first and last channels are
//                                                   4 columns narrower (ref cols stripped)
//     Otherwise                                  → uniform layout with exact division
//
// Reference-pixel layout for NISP 2040-wide science frame (32 × 64 − 8 = 2040):
//   channel 0    : cols [0,          59]   (60 cols, 4 ref cols stripped on left)
//   channels 1-30: cols [60+(k-1)×64, ...]  (64 cols each)
//   channel 31   : cols [1980,       2039]  (60 cols, 4 ref cols stripped on right)

    std::vector<OneOverFNoise::ChannelStripe> OneOverFNoise::_buildLayout(size_t width) const
    {
        // Determine channel count
        const size_t nCh = (m_nChannels > 0)
                           ? m_nChannels
                           : static_cast<size_t>(
                                 std::lround(static_cast<double>(width) /
                                             static_cast<double>(m_colPerChannel)));

        if (nCh == 0)
            throw std::invalid_argument(
                "OneOverFNoise::_buildLayout: derived nChannels == 0 for width " +
                std::to_string(width) + " and colPerChannel " + std::to_string(m_colPerChannel) + ".");

        const size_t nominal_total = nCh * m_colPerChannel;
        const bool   has_ref_cols  = (nominal_total == width + 8);
        const bool   is_uniform    = (nominal_total == width);

        std::vector<ChannelStripe> layout(nCh);

        if (has_ref_cols)
        {
            // 4 reference columns stripped from each side.
            const size_t inner_w = m_colPerChannel;
            const size_t edge_w  = m_colPerChannel - 4;

            layout[0] = {0, edge_w};
            size_t col = edge_w;
            for (size_t ch = 1; ch + 1 < nCh; ++ch)
            {
                layout[ch] = {col, inner_w};
                col += inner_w;
            }
            layout[nCh - 1] = {col, edge_w};
        }
        else
        {
            // Uniform layout (also used as fallback for non-standard widths)
            size_t col_w = is_uniform ? m_colPerChannel : (width / nCh);
            size_t rem   = is_uniform ? 0              : (width % nCh);
            size_t col   = 0;
            for (size_t ch = 0; ch < nCh; ++ch)
            {
                const size_t w = col_w + (ch < rem ? 1 : 0);
                layout[ch] = {col, w};
                col += w;
            }
        }

        return layout;
    }

// ── _buildPoles ───────────────────────────────────────────────────────────────
//
// Pole k has corner frequency f_k log-spaced in [f_min, 0.5]:
//   f_min = 1 / (2 × N)   one cycle per full channel readout
//   f_k   = exp( log(f_min) + k × Δ(log f) )
//   ρ_k   = exp( −2π f_k )
//
// Weight rule: to make Σ S_k(f) ≈ 1/f^slope, each pole's spectral energy at
// its corner frequency must scale as f_k^(−slope).  Because the pole's
// Lorentzian contributes S_k(f_k) ∝ σ_k²/f_k and Δ(log f) is constant:
//   σ_k ∝ f_k^((1−slope)/2)
//
// Poles are rescaled so Σ_k Var_k = Σ_k σ_k²/(1−ρ_k²) = 1 (unit variance).
// The caller multiplies by m_amplitude to obtain electrons.

    std::vector<OneOverFNoise::Pole> OneOverFNoise::_buildPoles(size_t channelLength) const
    {
        constexpr double twoPi = 2.0 * std::numbers::pi;

        const double f_min    = 1.0 / (2.0 * static_cast<double>(channelLength));
        const double lf_min   = std::log(f_min);
        const double lf_max   = std::log(0.5);
        const double log_step = (m_nPoles > 1)
                                ? (lf_max - lf_min) / static_cast<double>(m_nPoles - 1)
                                : 0.0;

        std::vector<Pole> poles(m_nPoles);
        double total_var = 0.0;

        for (size_t k = 0; k < m_nPoles; ++k)
        {
            const double fk        = std::exp(lf_min + static_cast<double>(k) * log_step);
            const double rho       = std::exp(-twoPi * fk);
            const double raw_sigma = std::pow(fk, 0.5 * (1.0 - m_slope));

            poles[k].rho   = rho;
            poles[k].sigma = raw_sigma;
            total_var += (raw_sigma * raw_sigma) / (1.0 - rho * rho);
        }

        const double scale = (total_var > 0.0) ? 1.0 / std::sqrt(total_var) : 1.0;
        for (auto& p : poles)
            p.sigma *= scale;

        return poles;
    }

// ── _fillChannel ──────────────────────────────────────────────────────────────
//
// Advances all nPoles AR(1) poles sample by sample, mutating `state` in place.
// Callers own `state`'s lifetime: generate() passes a fresh zero-initialized
// vector (cold start every call); Stream::next() passes a persistent vector
// so the recursion continues across calls instead of restarting.
// thread_local RNG so concurrent channels on different threads are independent.

    void OneOverFNoise::_fillChannel(const std::vector<Pole>& poles, std::vector<double>& state,
                                     size_t length, double* out) const
    {
        static thread_local std::mt19937_64 rng(
            static_cast<uint64_t>(
                std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        static thread_local std::normal_distribution<double> randn(0.0, 1.0);

        const size_t nP = poles.size();

        for (size_t t = 0; t < length; ++t)
        {
            double sample = 0.0;
            for (size_t k = 0; k < nP; ++k)
            {
                state[k] = poles[k].rho * state[k] + poles[k].sigma * randn(rng);
                sample  += state[k];
            }
            out[t] = sample * m_amplitude;
        }
    }

// ── generate (full frame) ─────────────────────────────────────────────────────
//
// One std::async task per output channel.  Each task:
//   1. Builds AR(1) poles for its channel length.
//   2. Fills a local buffer with the noise time series.
//   3. Scatters the buffer back onto the correct columns of the output map.
//
// Tasks write to non-overlapping column ranges → no data race.

    std::valarray<double> OneOverFNoise::generate(size_t width, size_t height) const
    {
        std::valarray<double> noise(0.0, width * height);
        if (!isEnabled())
            return noise;

        const auto layout = _buildLayout(width);
        const size_t nCh  = layout.size();

        std::vector<std::future<void>> futures;
        futures.reserve(nCh);

        for (size_t ch = 0; ch < nCh; ++ch)
        {
            const ChannelStripe st = layout[ch];   // copy by value for capture

            futures.push_back(std::async(std::launch::async,
                [this, st, &noise, width, height]()
                {
                    const size_t channel_len = height * st.col_width;
                    const auto   poles       = _buildPoles(channel_len);
                    std::vector<double> state(poles.size(), 0.0);

                    std::vector<double> buf(channel_len);
                    _fillChannel(poles, state, channel_len, buf.data());

                    // Scatter time series back to 2D map (row-major readout order)
                    for (size_t row = 0; row < height; ++row)
                    {
                        const double* src = buf.data() + row * st.col_width;
                        double*       dst = &noise[row * width + st.col_start];
                        for (size_t s = 0; s < st.col_width; ++s)
                            dst[s] = src[s];
                    }
                }));
        }

        for (auto& f : futures)
            f.get();   // re-throws any exception from the channel task

        return noise;
    }

// ── generate (windowed) ───────────────────────────────────────────────────────
//
// Only channels that overlap [c0, r0, c0+ncols, r0+nrows) are generated.
// Each overlapping channel still produces its full time series so that the
// AR(1) correlation from the start of the readout is preserved.
// Output is ncols × nrows, row-major.

    std::valarray<double> OneOverFNoise::generate(size_t width, size_t height,
                                                   size_t c0, size_t r0,
                                                   size_t ncols, size_t nrows) const
    {
        std::valarray<double> noise(0.0, ncols * nrows);
        if (!isEnabled())
            return noise;

        if (c0 + ncols > width || r0 + nrows > height)
            throw std::out_of_range(
                "OneOverFNoise::generate: window [" +
                std::to_string(c0) + "+" + std::to_string(ncols) + ", " +
                std::to_string(r0) + "+" + std::to_string(nrows) +
                ") exceeds detector (" +
                std::to_string(width) + "×" + std::to_string(height) + ").");

        const auto layout = _buildLayout(width);

        // Collect channels that overlap the requested column range
        std::vector<size_t> overlap;
        for (size_t ch = 0; ch < layout.size(); ++ch)
        {
            const auto& st = layout[ch];
            if (st.col_start + st.col_width > c0 && st.col_start < c0 + ncols)
                overlap.push_back(ch);
        }

        std::vector<std::future<void>> futures;
        futures.reserve(overlap.size());

        for (size_t ch : overlap)
        {
            const ChannelStripe st = layout[ch];

            futures.push_back(std::async(std::launch::async,
                [this, st, &noise, width, height, c0, r0, ncols, nrows]()
                {
                    const size_t channel_len = height * st.col_width;
                    const auto   poles       = _buildPoles(channel_len);
                    std::vector<double> state(poles.size(), 0.0);

                    std::vector<double> buf(channel_len);
                    _fillChannel(poles, state, channel_len, buf.data());

                    // Intersection of this channel's stripe with the window
                    const size_t win_c0 = std::max(st.col_start, c0);
                    const size_t win_c1 = std::min(st.col_start + st.col_width, c0 + ncols);

                    for (size_t r = 0; r < nrows; ++r)
                    {
                        // Pointer to the row in the full channel time series
                        const double* src_row = buf.data() + (r0 + r) * st.col_width;

                        for (size_t c = win_c0; c < win_c1; ++c)
                        {
                            // Local column index within the channel stripe
                            const size_t sample_idx = c - st.col_start;
                            // Row-major index in the output window
                            const size_t map_idx    = r * ncols + (c - c0);
                            noise[map_idx] = src_row[sample_idx];
                        }
                    }
                }));
        }

        for (auto& f : futures)
            f.get();

        return noise;
    }

// ── Stream ────────────────────────────────────────────────────────────────────
//
// Builds the channel layout and one AR(1) pole set per channel once, and owns
// a persistent per-channel state vector that next() advances in place —
// unlike generate(), which restarts every channel's AR(1) state at zero on
// every call.

    OneOverFNoise::Stream::Stream(const OneOverFNoise& owner, size_t width, size_t height)
        : m_owner(&owner), m_width(width), m_height(height)
    {
        m_layout = owner._buildLayout(width);
        m_poles.reserve(m_layout.size());
        m_state.reserve(m_layout.size());

        for (const auto& st : m_layout)
        {
            const size_t channel_len = height * st.col_width;
            m_poles.push_back(owner._buildPoles(channel_len));
            m_state.emplace_back(m_poles.back().size(), 0.0);   // cold start
        }
    }

    OneOverFNoise::Stream OneOverFNoise::beginStream(size_t width, size_t height) const
    {
        return Stream(*this, width, height);
    }

    std::valarray<double> OneOverFNoise::Stream::next()
    {
        std::valarray<double> noise(0.0, m_width * m_height);

        std::vector<std::future<void>> futures;
        futures.reserve(m_layout.size());

        for (size_t ch = 0; ch < m_layout.size(); ++ch)
        {
            futures.push_back(std::async(std::launch::async,
                [this, ch, &noise]()
                {
                    const ChannelStripe&  st    = m_layout[ch];
                    const auto&           poles = m_poles[ch];
                    std::vector<double>&  state = m_state[ch];   // persists across next() calls

                    const size_t channel_len = m_height * st.col_width;
                    std::vector<double> buf(channel_len);
                    m_owner->_fillChannel(poles, state, channel_len, buf.data());

                    for (size_t row = 0; row < m_height; ++row)
                    {
                        const double* src = buf.data() + row * st.col_width;
                        double*       dst = &noise[row * m_width + st.col_start];
                        for (size_t s = 0; s < st.col_width; ++s)
                            dst[s] = src[s];
                    }
                }));
        }

        for (auto& f : futures)
            f.get();   // re-throws any exception from the channel task

        return noise;
    }

#pragma endregion

} // namespace nisp
