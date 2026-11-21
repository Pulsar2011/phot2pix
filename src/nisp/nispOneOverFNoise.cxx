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
#include <cstdint>
#include <functional>
#include <future>
#include <initializer_list>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>

namespace nisp
{

    namespace
    {
        // ── Random numbers ────────────────────────────────────────────────────
        //
        // One generator per thread, seeded from the clock and the thread id, so that
        // channel tasks running concurrently draw independent sequences.
        std::mt19937_64& threadRng()
        {
            static thread_local std::mt19937_64 rng = []()
            {
                const uint64_t clock = static_cast<uint64_t>(
                    std::chrono::high_resolution_clock::now().time_since_epoch().count());
                const uint64_t id    = static_cast<uint64_t>(
                    std::hash<std::thread::id>{}(std::this_thread::get_id()));
                std::seed_seq seq{static_cast<uint32_t>(clock), static_cast<uint32_t>(clock >> 32),
                                  static_cast<uint32_t>(id),    static_cast<uint32_t>(id >> 32)};
                return std::mt19937_64(seq);
            }();
            return rng;
        }

        std::normal_distribution<double>& threadNormal()
        {
            static thread_local std::normal_distribution<double> randn(0.0, 1.0);
            return randn;
        }

        double gaussian()
        {
            return threadNormal()(threadRng());
        }

        // ── Configuration helpers ─────────────────────────────────────────────

        // Value of `key` when it holds a plain number; nullopt when it holds a JSON
        // object or array (or anything else that is not a number).
        std::optional<double> scalarValue(const phot2pix::PixelConfig& cfg, const std::string& key)
        {
            try { return cfg.getAsDouble(key); }
            catch (const std::invalid_argument&) {}       // not numeric
            catch (const nlohmann::json::type_error&) {}  // stored as a JSON object/array
            return std::nullopt;
        }

        // Reject keys of the JSON object j that are not in `allowed` (catches typos).
        void checkKeys(const nlohmann::json& j, std::initializer_list<const char*> allowed,
                       const std::string& context)
        {
            for (auto it = j.begin(); it != j.end(); ++it)
                if (std::none_of(allowed.begin(), allowed.end(),
                                 [&it](const char* k) { return it.key() == k; }))
                    throw std::invalid_argument(context + ": unknown key '" + it.key() + "'.");
        }

        // Number `name` of the JSON object j, or `fallback` when absent.
        double jsonNumber(const nlohmann::json& j, const char* name, double fallback,
                          const std::string& context)
        {
            if (!j.contains(name))
                return fallback;
            const auto& v = j.at(name);
            if (!v.is_number())
                throw std::invalid_argument(context + ": '" + name + "' must be a number.");
            return v.get<double>();
        }

        // Integer >= 1 `name` of the JSON object j, or `fallback` when absent.
        size_t jsonCount(const nlohmann::json& j, const char* name, size_t fallback,
                         const std::string& context)
        {
            const double v = jsonNumber(j, name, static_cast<double>(fallback), context);
            if (!(v >= 1.0) || v != std::floor(v))
                throw std::invalid_argument(context + ": '" + name + "' must be an integer >= 1.");
            return static_cast<size_t>(v);
        }

        // The terms listed under `name` in the drift object: one object or an array of objects.
        std::vector<nlohmann::json> jsonTerms(const nlohmann::json& drift, const char* name,
                                              const std::string& context)
        {
            std::vector<nlohmann::json> out;
            if (!drift.contains(name))
                return out;

            const auto& v = drift.at(name);
            if (v.is_object())
                out.push_back(v);
            else if (v.is_array())
                for (const auto& term : v)
                {
                    if (!term.is_object())
                        throw std::invalid_argument(context + ": every term must be an object.");
                    out.push_back(term);
                }
            else
                throw std::invalid_argument(context + " must be an object or an array of objects.");

            return out;
        }

        // Scatter for a full-frame map: one channel's series (row-major within its
        // stripe) is copied onto its columns of the width × height map.
        auto fullFrameScatter(std::valarray<double>& noise, size_t width, size_t height)
        {
            return [&noise, width, height](const auto& st, const std::vector<double>& series)
            {
                for (size_t row = 0; row < height; ++row)
                {
                    const double* src = series.data() + row * st.col_width;
                    double*       dst = &noise[row * width + st.col_start];
                    for (size_t s = 0; s < st.col_width; ++s)
                        dst[s] = src[s];
                }
            };
        }
    }

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
        // ── oneOverF: per-channel 1/f noise ───────────────────────────────────
        if (cfg.has("oneOverF"))
        {
            // Scalar shorthand: "oneOverF": 7.0 → amplitude only, keep the other parameters.
            if (const auto amplitude = scalarValue(cfg, "oneOverF"))
                m_amplitude = *amplitude;
            else
            {
                const auto& j = cfg.getJson("oneOverF");

                if (j.contains("amplitude"))     m_amplitude    = j.at("amplitude").get<double>();
                if (j.contains("slope"))         m_slope        = j.at("slope").get<double>();
                if (j.contains("nChannels"))     m_nChannels    = j.at("nChannels").get<size_t>();
                if (j.contains("colPerChannel")) m_colPerChannel = j.at("colPerChannel").get<size_t>();
                if (j.contains("nPoles"))        m_nPoles       = j.at("nPoles").get<size_t>();
            }

            if (m_amplitude < 0.0)
                throw std::invalid_argument("OneOverFNoise: oneOverF amplitude must be >= 0.");
            if (m_colPerChannel < 1)
                throw std::invalid_argument("OneOverFNoise: colPerChannel must be >= 1.");
            if (m_nPoles < 1)
                throw std::invalid_argument("OneOverFNoise: nPoles must be >= 1.");
        }

        // ── bias_voltage_noise: 1/f noise common to every channel ─────────────
        if (cfg.has("bias_voltage_noise"))
        {
            const std::string ctx = "OneOverFNoise: bias_voltage_noise";

            if (const auto amplitude = scalarValue(cfg, "bias_voltage_noise"))
                m_biasAmplitude = *amplitude;   // scalar shorthand: amplitude only
            else
            {
                const auto& j = cfg.getJson("bias_voltage_noise");
                if (!j.is_object())
                    throw std::invalid_argument(ctx + " must be a number or an object.");
                checkKeys(j, {"amplitude", "slope", "nPoles"}, ctx);

                m_biasAmplitude = jsonNumber(j, "amplitude", m_biasAmplitude, ctx);
                m_biasSlope     = jsonNumber(j, "slope",     m_biasSlope,     ctx);
                m_biasPoles     = jsonCount (j, "nPoles",    m_biasPoles,     ctx);
            }

            if (m_biasAmplitude < 0.0)
                throw std::invalid_argument(ctx + ": amplitude must be >= 0.");
        }

        // ── pedestal_drift: f(t) + sigma · ε per channel and per frame ────────
        if (cfg.has("pedestal_drift"))
        {
            const std::string ctx = "OneOverFNoise: pedestal_drift";

            // Every form replaces the whole pedestal configuration.
            m_pedestalSigma = 0.0;
            m_pedestalT0    = 0.0;
            m_driftPoly.clear();
            m_driftExp.clear();
            m_driftSin.clear();

            if (const auto sigma = scalarValue(cfg, "pedestal_drift"))
                m_pedestalSigma = *sigma;   // scalar shorthand: white noise only, no drift
            else
            {
                const auto& j = cfg.getJson("pedestal_drift");
                if (!j.is_object())
                    throw std::invalid_argument(ctx + " must be a number or an object.");
                checkKeys(j, {"sigma", "t0", "drift"}, ctx);

                m_pedestalSigma = jsonNumber(j, "sigma", 0.0, ctx);
                m_pedestalT0    = jsonNumber(j, "t0",    0.0, ctx);

                if (j.contains("drift"))
                {
                    const auto&       d    = j.at("drift");
                    const std::string dctx = ctx + ": 'drift'";
                    if (!d.is_object())
                        throw std::invalid_argument(dctx + " must be an object with 'poly', 'exp' and/or 'sin'.");
                    checkKeys(d, {"poly", "exp", "sin"}, dctx);

                    if (d.contains("poly"))
                    {
                        const auto& p = d.at("poly");
                        if (!p.is_array())
                            throw std::invalid_argument(dctx + ": 'poly' must be an array [a0, a1, ...].");
                        for (const auto& a : p)
                        {
                            if (!a.is_number())
                                throw std::invalid_argument(dctx + ": every 'poly' coefficient must be a number.");
                            m_driftPoly.push_back(a.get<double>());
                        }
                    }

                    const std::string ectx = dctx + ".exp";
                    for (const auto& term : jsonTerms(d, "exp", ectx))
                    {
                        checkKeys(term, {"amplitude", "tau"}, ectx);
                        const double tau = jsonNumber(term, "tau", 0.0, ectx);
                        if (!(tau > 0.0))
                            throw std::invalid_argument(ectx + ": 'tau' must be given and > 0.");
                        m_driftExp.push_back({jsonNumber(term, "amplitude", 0.0, ectx), tau});
                    }

                    const std::string sctx = dctx + ".sin";
                    for (const auto& term : jsonTerms(d, "sin", sctx))
                    {
                        checkKeys(term, {"amplitude", "period", "phase"}, sctx);
                        const double period = jsonNumber(term, "period", 0.0, sctx);
                        if (!(period > 0.0))
                            throw std::invalid_argument(sctx + ": 'period' must be given and > 0.");
                        m_driftSin.push_back({jsonNumber(term, "amplitude", 0.0, sctx), period,
                                              jsonNumber(term, "phase", 0.0, sctx)});
                    }
                }
            }

            if (m_pedestalSigma < 0.0)
                throw std::invalid_argument(ctx + ": sigma must be >= 0.");
        }
    }

    void OneOverFNoise::setFrameTime(double frameTime)
    {
        if (!(frameTime > 0.0))
            throw std::invalid_argument("OneOverFNoise::setFrameTime: frame time must be > 0 (got " +
                                        std::to_string(frameTime) + ").");
        m_frameTime = frameTime;
    }

// ── Pedestal drift ────────────────────────────────────────────────────────────
//
//   f(t) = Σ_k a_k t^k + Σ_j A_j exp(−t/τ_j) + Σ_m B_m sin(2π t/P_m + φ_m),   t = t0 + time since reset

    double OneOverFNoise::pedestalMean(double t) const
    {
        constexpr double twoPi = 2.0 * std::numbers::pi;
        const double     tt    = m_pedestalT0 + t;

        double f = 0.0;
        for (size_t k = m_driftPoly.size(); k-- > 0; )   // Horner
            f = f * tt + m_driftPoly[k];
        for (const auto& e : m_driftExp)
            f += e.amplitude * std::exp(-tt / e.tau);
        for (const auto& s : m_driftSin)
            f += s.amplitude * std::sin(twoPi * tt / s.period + s.phase);

        return f;
    }

    std::vector<double> OneOverFNoise::_pedestalRows(size_t height, double tStart) const
    {
        std::vector<double> rows;
        if (!hasPedestalDrift())
            return rows;

        // Row y of the frame is read at tStart + (y + 0.5)/H · t_f.
        rows.resize(height);
        const double h = static_cast<double>(height);
        for (size_t y = 0; y < height; ++y)
            rows[y] = pedestalMean(tStart + (static_cast<double>(y) + 0.5) / h * m_frameTime);

        return rows;
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
//   channel 0    : cols [0,          59]   (60 cols, 4 ref cols stripped on left, time offset 4)
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

            layout[0] = {0, edge_w, inner_w - edge_w};   // its reference columns were read first
            size_t col = edge_w;
            for (size_t ch = 1; ch + 1 < nCh; ++ch)
            {
                layout[ch] = {col, inner_w, 0};
                col += inner_w;
            }
            layout[nCh - 1] = {col, edge_w, 0};          // its reference columns are read last
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
                layout[ch] = {col, w, 0};
                col += w;
            }
        }

        return layout;
    }

    size_t OneOverFNoise::_rowSamples(const std::vector<ChannelStripe>& layout)
    {
        size_t n = 0;
        for (const auto& st : layout)
            n = std::max(n, st.time_offset + st.col_width);
        return n;
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
// The caller multiplies by the component's amplitude to obtain electrons.

    std::vector<OneOverFNoise::Pole> OneOverFNoise::_buildPoles(size_t channelLength, double slope, size_t nPoles)
    {
        constexpr double twoPi = 2.0 * std::numbers::pi;

        const double f_min    = 1.0 / (2.0 * static_cast<double>(channelLength));
        const double lf_min   = std::log(f_min);
        const double lf_max   = std::log(0.5);
        const double log_step = (nPoles > 1)
                                ? (lf_max - lf_min) / static_cast<double>(nPoles - 1)
                                : 0.0;

        std::vector<Pole> poles(nPoles);
        double total_var = 0.0;

        for (size_t k = 0; k < nPoles; ++k)
        {
            const double fk        = std::exp(lf_min + static_cast<double>(k) * log_step);
            const double rho       = std::exp(-twoPi * fk);
            const double raw_sigma = std::pow(fk, 0.5 * (1.0 - slope));

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
                                     size_t length, double amplitude, double* out)
    {
        // Thread-local generator and distribution looked up once, not once per sample.
        std::mt19937_64&                  rng   = threadRng();
        std::normal_distribution<double>& randn = threadNormal();

        const size_t nP = poles.size();

        for (size_t t = 0; t < length; ++t)
        {
            double sample = 0.0;
            for (size_t k = 0; k < nP; ++k)
            {
                state[k] = poles[k].rho * state[k] + poles[k].sigma * randn(rng);
                sample  += state[k];
            }
            out[t] = sample * amplitude;
        }
    }

// ── AR(1) states ──────────────────────────────────────────────────────────────
//
// Poles and zero states for the channels listed in `channels` (per-channel 1/f)
// and for the bias voltage series. A component that is off gets empty vectors:
// _frame() never advances it.

    void OneOverFNoise::_initChannels(const std::vector<ChannelStripe>& layout,
                                      const std::vector<size_t>& channels, size_t height,
                                      std::vector<std::vector<Pole>>& poles,
                                      std::vector<std::vector<double>>& states) const
    {
        poles.assign(layout.size(), std::vector<Pole>());
        states.assign(layout.size(), std::vector<double>());
        if (m_amplitude <= 0.0)
            return;

        for (size_t ch : channels)
        {
            poles[ch] = _buildPoles(height * layout[ch].col_width, m_slope, m_nPoles);
            states[ch].assign(poles[ch].size(), 0.0);   // cold start
        }
    }

    void OneOverFNoise::_initCommon(const std::vector<ChannelStripe>& layout, size_t height,
                                    std::vector<Pole>& poles, std::vector<double>& state) const
    {
        poles.clear();
        state.clear();
        if (m_biasAmplitude <= 0.0)
            return;

        poles = _buildPoles(height * _rowSamples(layout), m_biasSlope, m_biasPoles);
        state.assign(poles.size(), 0.0);   // cold start
    }

// ── _frame ────────────────────────────────────────────────────────────────────
//
// One frame of readout noise. The bias voltage series is generated in its own
// task, concurrently with the channel tasks; every channel task waits for it
// through its own copy of the shared_future, then adds it at the channel's
// pixel-clock index. Tasks write to non-overlapping column ranges → no data race.

    template <typename Scatter>
    void OneOverFNoise::_frame(const std::vector<ChannelStripe>& layout,
                               const std::vector<size_t>& channels,
                               const std::vector<std::vector<Pole>>& chanPoles,
                               std::vector<std::vector<double>>& chanStates,
                               const std::vector<Pole>& commonPoles,
                               std::vector<double>& commonState,
                               size_t height, double tStart,
                               const Scatter& scatter) const
    {
        const size_t rowSamples = _rowSamples(layout);

        // Bias voltage noise: one series for every channel.
        std::shared_future<std::vector<double>> bias;
        if (m_biasAmplitude > 0.0)
            bias = std::async(std::launch::async,
                [this, &commonPoles, &commonState, height, rowSamples]()
                {
                    std::vector<double> series(height * rowSamples);
                    _fillChannel(commonPoles, commonState, series.size(), m_biasAmplitude, series.data());
                    return series;
                }).share();

        // Pedestal drift of each row (same for every channel); empty without a drift equation.
        const std::vector<double> rowDrift = _pedestalRows(height, tStart);

        std::vector<std::future<void>> futures;
        futures.reserve(channels.size());

        for (size_t ch : channels)
        {
            futures.push_back(std::async(std::launch::async,
                [this, ch, bias, &layout, &chanPoles, &chanStates, &rowDrift, &scatter, height, rowSamples]()
                {
                    const ChannelStripe& st  = layout[ch];
                    const size_t         len = height * st.col_width;

                    // 1/f noise of the channel's own amplifier.
                    std::vector<double> series(len, 0.0);
                    if (m_amplitude > 0.0)
                        _fillChannel(chanPoles[ch], chanStates[ch], len, m_amplitude, series.data());

                    // Pedestal noise: one draw for the whole channel, for this frame.
                    const double offset = (m_pedestalSigma > 0.0) ? m_pedestalSigma * gaussian() : 0.0;

                    const std::vector<double>* common = bias.valid() ? &bias.get() : nullptr;

                    if (common || offset != 0.0 || !rowDrift.empty())
                        for (size_t row = 0; row < height; ++row)
                        {
                            double*      s     = series.data() + row * st.col_width;
                            const double level = offset + (rowDrift.empty() ? 0.0 : rowDrift[row]);

                            if (common)
                            {
                                const double* b = common->data() + row * rowSamples + st.time_offset;
                                for (size_t k = 0; k < st.col_width; ++k)
                                    s[k] += b[k] + level;
                            }
                            else
                                for (size_t k = 0; k < st.col_width; ++k)
                                    s[k] += level;
                        }

                    scatter(st, series);
                }));
        }

        for (auto& f : futures)
            f.get();   // re-throws any exception from a channel task (or from the bias task)
    }

// ── generate (full frame) ─────────────────────────────────────────────────────
//
// One independent read: every AR(1) state starts from zero.

    std::valarray<double> OneOverFNoise::generate(size_t width, size_t height, double tStart) const
    {
        std::valarray<double> noise(0.0, width * height);
        if (!isEnabled())
            return noise;

        const auto layout = _buildLayout(width);

        std::vector<size_t> channels(layout.size());
        for (size_t ch = 0; ch < channels.size(); ++ch)
            channels[ch] = ch;

        std::vector<std::vector<Pole>>   chanPoles;
        std::vector<std::vector<double>> chanStates;
        std::vector<Pole>                commonPoles;
        std::vector<double>              commonState;
        _initChannels(layout, channels, height, chanPoles, chanStates);
        _initCommon(layout, height, commonPoles, commonState);

        _frame(layout, channels, chanPoles, chanStates, commonPoles, commonState,
               height, tStart, fullFrameScatter(noise, width, height));

        return noise;
    }

// ── generate (windowed) ───────────────────────────────────────────────────────
//
// Only channels that overlap [c0, r0, c0+ncols, r0+nrows) are generated.
// Each overlapping channel still produces its full time series so that the
// AR(1) correlation from the start of the readout is preserved; the bias
// voltage series is always generated in full.
// Output is ncols × nrows, row-major.

    std::valarray<double> OneOverFNoise::generate(size_t width, size_t height,
                                                   size_t c0, size_t r0,
                                                   size_t ncols, size_t nrows,
                                                   double tStart) const
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

        std::vector<std::vector<Pole>>   chanPoles;
        std::vector<std::vector<double>> chanStates;
        std::vector<Pole>                commonPoles;
        std::vector<double>              commonState;
        _initChannels(layout, overlap, height, chanPoles, chanStates);
        _initCommon(layout, height, commonPoles, commonState);

        _frame(layout, overlap, chanPoles, chanStates, commonPoles, commonState, height, tStart,
            [&noise, c0, r0, ncols, nrows](const ChannelStripe& st, const std::vector<double>& series)
            {
                // Intersection of this channel's stripe with the window
                const size_t win_c0 = std::max(st.col_start, c0);
                const size_t win_c1 = std::min(st.col_start + st.col_width, c0 + ncols);

                for (size_t r = 0; r < nrows; ++r)
                {
                    // Row r0 + r of the channel's full time series
                    const double* src_row = series.data() + (r0 + r) * st.col_width;

                    for (size_t c = win_c0; c < win_c1; ++c)
                        noise[r * ncols + (c - c0)] = src_row[c - st.col_start];
                }
            });

        return noise;
    }

// ── Stream ────────────────────────────────────────────────────────────────────
//
// Builds the channel layout and the AR(1) poles once, and owns persistent
// states (one per channel, one for the bias voltage series) that next()
// advances in place — unlike generate(), which restarts every state at zero
// on every call. It also counts the frames, so frame k is read from k·t_f.

    OneOverFNoise::Stream::Stream(const OneOverFNoise& owner, size_t width, size_t height)
        : m_owner(&owner), m_width(width), m_height(height)
    {
        m_layout = owner._buildLayout(width);

        m_channels.resize(m_layout.size());
        for (size_t ch = 0; ch < m_channels.size(); ++ch)
            m_channels[ch] = ch;

        owner._initChannels(m_layout, m_channels, height, m_poles, m_state);
        owner._initCommon(m_layout, height, m_commonPoles, m_commonState);
    }

    OneOverFNoise::Stream OneOverFNoise::beginStream(size_t width, size_t height) const
    {
        return Stream(*this, width, height);
    }

    std::valarray<double> OneOverFNoise::Stream::next()
    {
        std::valarray<double> noise(0.0, m_width * m_height);

        const double tStart = static_cast<double>(m_frame) * m_owner->m_frameTime;
        m_owner->_frame(m_layout, m_channels, m_poles, m_state, m_commonPoles, m_commonState,
                        m_height, tStart, fullFrameScatter(noise, m_width, m_height));
        ++m_frame;

        return noise;
    }

#pragma endregion

} // namespace nisp
