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

    void OneOverFNoise::checkGeometry(size_t width, size_t height) const
    {
        const auto [x0, y0] = m_geometry.place(width, height);
        _segments(Region{x0, y0, width, height});
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

    std::vector<double> OneOverFNoise::_pedestalRows(size_t py0, size_t nrows, double tStart) const
    {
        std::vector<double> rows;
        if (!hasPedestalDrift())
            return rows;

        // Physical row Y of the frame is read at tStart + (Y + 0.5)/height · t_f.
        rows.resize(nrows);
        const double h = static_cast<double>(m_geometry.height);
        for (size_t r = 0; r < nrows; ++r)
            rows[r] = pedestalMean(tStart + (static_cast<double>(py0 + r) + 0.5) / h * m_frameTime);

        return rows;
    }

// ── Channel layout ────────────────────────────────────────────────────────────
//
// Channel c covers the physical columns [c·w, (c+1)·w), w = colPerChannel. A region
// overlaps the channels px0/w … (px0+ncols−1)/w; the first and the last may be
// partial. For the NISP 2048-wide detector and a matrix at x0:
//   2048 wide (x0 = 0) : 32 full channels
//   2040 wide (x0 = 4) : channels 0 and 31 lose their 4 reference columns (60 columns;
//                        channel 0 starts at pixel-clock offset 4)
//   2044 wide (x0 = 2) : channels 0 and 31 keep 2 reference columns (62 columns)
//   256 wide  (x0 = 896): channels 14 … 17

    std::vector<OneOverFNoise::Segment> OneOverFNoise::_segments(const Region& region) const
    {
        const size_t w    = m_colPerChannel;
        const size_t detW = m_geometry.width;

        if (w < 1 || detW % w != 0)
            throw std::invalid_argument("OneOverFNoise: colPerChannel (" + std::to_string(w) +
                                        ") must divide the detector width (" + std::to_string(detW) + ").");
        if (m_nChannels > 0 && m_nChannels != detW / w)
            throw std::invalid_argument("OneOverFNoise: nChannels (" + std::to_string(m_nChannels) +
                                        ") must equal detector width / colPerChannel (" +
                                        std::to_string(detW / w) + ").");
        if (region.px0 + region.ncols > detW || region.py0 + region.nrows > m_geometry.height)
            throw std::invalid_argument("OneOverFNoise: region outside the " + std::to_string(detW) + " × " +
                                        std::to_string(m_geometry.height) + " detector.");

        std::vector<Segment> segments;
        if (region.ncols == 0 || region.nrows == 0)
            return segments;

        const size_t x1 = region.px0 + region.ncols;
        for (size_t c = region.px0 / w; c * w < x1; ++c)
        {
            const size_t a = std::max(region.px0, c * w);
            const size_t b = std::min(x1, (c + 1) * w);
            segments.push_back({c, a - region.px0, b - a, a - c * w});
        }

        return segments;
    }

    std::pair<size_t, size_t> OneOverFNoise::_clockSpan(const std::vector<Segment>& segments)
    {
        if (segments.empty())
            return {0, 0};

        size_t lo = segments.front().offset, hi = 0;
        for (const auto& s : segments)
        {
            lo = std::min(lo, s.offset);
            hi = std::max(hi, s.offset + s.width);
        }
        return {lo, hi - lo};
    }

// ── _buildPoles ───────────────────────────────────────────────────────────────
//
// Pole k has corner frequency f_k log-spaced in [f_min, 0.5]:
//   f_min = 1 / (2 × N)   one cycle per two frames (N = samples per channel per frame)
//   f_k   = exp( log(f_min) + k × Δ(log f) )
//   ρ_k   = exp( −2π f_k )
//
// Weight rule: a 1/f^slope spectrum holds a variance ∝ f^(1−slope) per log step of
// frequency. The poles are log-spaced, so pole k is given the variance
//   Var_k ∝ f_k^(1−slope)        (slope 1: the same variance for every pole)
// and the AR(1) innovation that produces it, σ_k = √(Var_k · (1 − ρ_k²)).
//
// Poles are rescaled so Σ_k Var_k = 1 (unit variance).
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
            const double fk  = std::exp(lf_min + static_cast<double>(k) * log_step);
            const double rho = std::exp(-twoPi * fk);
            const double var = std::pow(fk, 1.0 - slope);

            poles[k].rho   = rho;
            poles[k].sigma = std::sqrt(var * (1.0 - rho * rho));
            total_var += var;
        }

        const double scale = (total_var > 0.0) ? 1.0 / std::sqrt(total_var) : 1.0;
        for (auto& p : poles)
            p.sigma *= scale;

        return poles;
    }

// ── AR(1) processes ───────────────────────────────────────────────────────────
//
// A process follows one series (one channel's 1/f, or the bias voltage series)
// through the frames, generating only `width` samples per row over `nrows` rows:
//
//   row r of the region  : width samples
//   between two rows     : rowGap   = w − width samples skipped
//   between two frames   : frameGap = (height − nrows)·w + rowGap samples skipped
//
// A skip of n samples is exact: x ← ρⁿ x + σ √((1 − ρ²ⁿ)/(1 − ρ²)) ε.
// The process starts from its stationary state, x_k ~ N(0, σ_k²/(1 − ρ_k²)).

    OneOverFNoise::Process OneOverFNoise::_process(double amplitude, double slope, size_t nPoles,
                                                   size_t width, size_t nrows) const
    {
        Process p;
        if (!(amplitude > 0.0) || width == 0 || nrows == 0)
            return p;

        const size_t w = m_colPerChannel;
        p.poles     = _buildPoles(m_geometry.height * w, slope, nPoles);   // spectrum of a physical frame
        p.rowGap    = w - width;
        p.frameGap  = (m_geometry.height - nrows) * w + p.rowGap;
        p.rowJump   = _jumps(p.poles, p.rowGap);
        p.frameJump = _jumps(p.poles, p.frameGap);

        p.state.resize(p.poles.size());
        for (size_t k = 0; k < p.poles.size(); ++k)
        {
            const Pole& q = p.poles[k];
            p.state[k] = std::sqrt(q.sigma * q.sigma / (1.0 - q.rho * q.rho)) * gaussian();   // stationary start
        }

        return p;
    }

    std::vector<OneOverFNoise::Process> OneOverFNoise::_channelProcesses(const std::vector<Segment>& segments,
                                                                         size_t nrows) const
    {
        std::vector<Process> processes;
        processes.reserve(segments.size());
        for (const auto& s : segments)
            processes.push_back(_process(m_amplitude, m_slope, m_nPoles, s.width, nrows));
        return processes;
    }

    std::vector<OneOverFNoise::Jump> OneOverFNoise::_jumps(const std::vector<Pole>& poles, size_t n)
    {
        std::vector<Jump> jumps(poles.size(), Jump{1.0, 0.0});
        if (n == 0)
            return jumps;

        const double dn = static_cast<double>(n);
        for (size_t k = 0; k < poles.size(); ++k)
        {
            const double rho = poles[k].rho;
            const double an  = std::pow(rho, dn);
            jumps[k] = {an, poles[k].sigma * std::sqrt((1.0 - an * an) / (1.0 - rho * rho))};
        }
        return jumps;
    }

    void OneOverFNoise::_advance(std::vector<double>& state, const std::vector<Jump>& jump)
    {
        std::mt19937_64&                  rng   = threadRng();
        std::normal_distribution<double>& randn = threadNormal();

        for (size_t k = 0; k < state.size(); ++k)
            state[k] = jump[k].a * state[k] + jump[k].b * randn(rng);
    }

// ── _fillChannel ──────────────────────────────────────────────────────────────
//
// Advances all nPoles AR(1) poles sample by sample, mutating `state` in place.
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

    void OneOverFNoise::_run(Process& p, size_t nrows, size_t width, double amplitude,
                             bool firstFrame, double* out)
    {
        for (size_t row = 0; row < nrows; ++row)
        {
            if (row > 0)
            {
                if (p.rowGap > 0)
                    _advance(p.state, p.rowJump);
            }
            else if (!firstFrame && p.frameGap > 0)
                _advance(p.state, p.frameJump);

            _fillChannel(p.poles, p.state, width, amplitude, out + row * width);
        }
    }

// ── _frame ────────────────────────────────────────────────────────────────────
//
// One frame of readout noise. The bias voltage series is generated in its own
// task, concurrently with the channel tasks; every channel task waits for it
// through its own copy of the shared_future, then adds it at the channel's
// pixel-clock offsets. Tasks write to non-overlapping column ranges → no data race.

    void OneOverFNoise::_frame(const Region& region, const std::vector<Segment>& segments,
                               std::vector<Process>& channels, Process& bias,
                               size_t clockOffset, size_t clockWidth,
                               bool firstFrame, double tStart, std::valarray<double>& out) const
    {
        const size_t nrows = region.nrows;
        const size_t ncols = region.ncols;

        // Bias voltage noise: one series for every channel.
        std::shared_future<std::vector<double>> biasSeries;
        if (bias.enabled())
            biasSeries = std::async(std::launch::async,
                [this, &bias, nrows, clockWidth, firstFrame]()
                {
                    std::vector<double> series(nrows * clockWidth);
                    _run(bias, nrows, clockWidth, m_biasAmplitude, firstFrame, series.data());
                    return series;
                }).share();

        // Pedestal drift of each row (same for every channel); empty without a drift equation.
        const std::vector<double> rowDrift = _pedestalRows(region.py0, nrows, tStart);

        std::vector<std::future<void>> futures;
        futures.reserve(segments.size());

        for (size_t i = 0; i < segments.size(); ++i)
        {
            futures.push_back(std::async(std::launch::async,
                [this, i, biasSeries, &segments, &channels, &rowDrift, &out,
                 nrows, ncols, clockOffset, clockWidth, firstFrame]()
                {
                    const Segment& sg = segments[i];

                    // 1/f noise of the channel's own amplifier.
                    std::vector<double> series(nrows * sg.width, 0.0);
                    if (channels[i].enabled())
                        _run(channels[i], nrows, sg.width, m_amplitude, firstFrame, series.data());

                    // Pedestal noise: one draw for the whole channel, for this frame.
                    const double offset = (m_pedestalSigma > 0.0) ? m_pedestalSigma * gaussian() : 0.0;

                    const std::vector<double>* common = biasSeries.valid() ? &biasSeries.get() : nullptr;

                    for (size_t r = 0; r < nrows; ++r)
                    {
                        const double  level = offset + (rowDrift.empty() ? 0.0 : rowDrift[r]);
                        const double* s     = series.data() + r * sg.width;
                        double*       dst   = &out[r * ncols + sg.outCol];

                        if (common)
                        {
                            const double* b = common->data() + r * clockWidth + (sg.offset - clockOffset);
                            for (size_t k = 0; k < sg.width; ++k)
                                dst[k] = s[k] + b[k] + level;
                        }
                        else
                            for (size_t k = 0; k < sg.width; ++k)
                                dst[k] = s[k] + level;
                    }
                }));
        }

        for (auto& f : futures)
            f.get();   // re-throws any exception from a channel task (or from the bias task)
    }

    std::valarray<double> OneOverFNoise::_read(const Region& region, double tStart) const
    {
        std::valarray<double> noise(0.0, region.ncols * region.nrows);

        const auto segments = _segments(region);
        if (segments.empty())
            return noise;

        // Fresh processes, from their stationary state: this read is independent of any other call.
        std::vector<Process> channels = _channelProcesses(segments, region.nrows);
        const auto [clockOffset, clockWidth] = _clockSpan(segments);
        Process bias = _process(m_biasAmplitude, m_biasSlope, m_biasPoles, clockWidth, region.nrows);

        _frame(region, segments, channels, bias, clockOffset, clockWidth, true, tStart, noise);
        return noise;
    }

// ── generate ──────────────────────────────────────────────────────────────────

    std::valarray<double> OneOverFNoise::generate(size_t width, size_t height, double tStart) const
    {
        if (!isEnabled())
            return std::valarray<double>(0.0, width * height);

        const auto [x0, y0] = m_geometry.place(width, height);
        return _read(Region{x0, y0, width, height}, tStart);
    }

    std::valarray<double> OneOverFNoise::generate(size_t width, size_t height,
                                                   size_t c0, size_t r0,
                                                   size_t ncols, size_t nrows,
                                                   double tStart) const
    {
        if (!isEnabled())
            return std::valarray<double>(0.0, ncols * nrows);

        if (c0 + ncols > width || r0 + nrows > height)
            throw std::out_of_range(
                "OneOverFNoise::generate: window [" +
                std::to_string(c0) + "+" + std::to_string(ncols) + ", " +
                std::to_string(r0) + "+" + std::to_string(nrows) +
                ") exceeds the matrix (" +
                std::to_string(width) + "×" + std::to_string(height) + ").");

        const auto [x0, y0] = m_geometry.place(width, height);
        return _read(Region{x0 + c0, y0 + r0, ncols, nrows}, tStart);
    }

// ── Stream ────────────────────────────────────────────────────────────────────
//
// Builds the segments and the processes once; next() advances them by one
// physical frame (skipping the rows outside the matrix) instead of restarting.
// It also counts the frames, so frame k is read from k·t_f.

    OneOverFNoise::Stream::Stream(const OneOverFNoise& owner, size_t width, size_t height)
        : m_owner(&owner), m_region{0, 0, width, height}
    {
        const auto [x0, y0] = owner.m_geometry.place(width, height);
        m_region.px0 = x0;
        m_region.py0 = y0;

        m_segments = owner._segments(m_region);
        m_channels = owner._channelProcesses(m_segments, height);

        const auto span = OneOverFNoise::_clockSpan(m_segments);
        m_clockOffset = span.first;
        m_clockWidth  = span.second;
        m_bias = owner._process(owner.m_biasAmplitude, owner.m_biasSlope, owner.m_biasPoles, m_clockWidth, height);
    }

    OneOverFNoise::Stream OneOverFNoise::beginStream(size_t width, size_t height) const
    {
        return Stream(*this, width, height);
    }

    std::valarray<double> OneOverFNoise::Stream::next()
    {
        std::valarray<double> noise(0.0, m_region.ncols * m_region.nrows);

        if (!m_segments.empty())
        {
            const double tStart = static_cast<double>(m_frame) * m_owner->m_frameTime;
            m_owner->_frame(m_region, m_segments, m_channels, m_bias, m_clockOffset, m_clockWidth,
                            m_frame == 0, tStart, noise);
        }
        ++m_frame;

        return noise;
    }

#pragma endregion

} // namespace nisp
