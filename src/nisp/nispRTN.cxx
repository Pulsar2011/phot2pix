//
// Created by gillard on 2026-11-27.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispRTN.cxx
//

#include "nispRTN.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace nisp
{

    namespace
    {
        const std::string kContext = "RTN: 'rtn'";

        constexpr double kAmplitudeMin = 1.0;   // Smallest amplitude of an RTN pixel (ADU)

        // Finite number `name` of the JSON object j; fallback when absent and not required.
        double jsonNumber(const nlohmann::json& j, const char* name, bool required, double fallback)
        {
            if (!j.contains(name))
            {
                if (required)
                    throw std::invalid_argument(kContext + ": '" + name + "' is required.");
                return fallback;
            }

            const auto& v = j.at(name);
            if (!v.is_number() || !std::isfinite(v.get<double>()))
                throw std::invalid_argument(kContext + ": '" + name + "' must be a finite number.");
            return v.get<double>();
        }
    }

// ── RTN ───────────────────────────────────────────────────────────────────────
#pragma region - RTN

    // Generator seeded by seed, or by the random device without seed.
    std::mt19937_64 RTN::makeRng(const std::optional<uint64_t>& seed)
    {
        if (seed)
            return std::mt19937_64(*seed);

        std::random_device rd;
        std::seed_seq      seq{rd(), rd(), rd(), rd()};
        return std::mt19937_64(seq);
    }

    void RTN::initialize(const phot2pix::PixelConfig& cfg, size_t width, size_t height)
    {
        if (!cfg.has("rtn"))
            return;

        const nlohmann::json* node = nullptr;
        try{
            node = &cfg.getJson("rtn");
        }
        catch (const std::invalid_argument&)
        {}   // not a JSON value

        if (!node || !node->is_object())
            throw std::invalid_argument(kContext + " must be an object {\"fraction\", \"tauMin\", \"tauMax\", \"amplitude\"}.");
        const nlohmann::json& j = *node;

        for (auto it = j.begin(); it != j.end(); ++it)
            if (it.key() != "fraction" && it.key() != "tauMin" && it.key() != "tauMax" &&
                it.key() != "amplitude" && it.key() != "positiveFraction" && it.key() != "seed")
                throw std::invalid_argument(kContext + ": unknown key '" + it.key() + "'.");

        // ── Values, all checked before anything is changed ──
        const double fraction         = jsonNumber(j, "fraction",         true,  0.0);
        const double tauMin           = jsonNumber(j, "tauMin",           true,  0.0);
        const double tauMax           = jsonNumber(j, "tauMax",           true,  0.0);
        const double amplitude        = jsonNumber(j, "amplitude",        true,  0.0);
        const double positiveFraction = jsonNumber(j, "positiveFraction", false, 0.4);

        if (!(fraction >= 0.0 && fraction <= 1.0))
            throw std::invalid_argument(kContext + ": 'fraction' must be in [0, 1].");
        if (!(tauMin > 0.0))
            throw std::invalid_argument(kContext + ": 'tauMin' must be > 0.");
        if (!(tauMax >= tauMin))
            throw std::invalid_argument(kContext + ": 'tauMax' must be >= 'tauMin'.");
        if (!(amplitude >= kAmplitudeMin))
            throw std::invalid_argument(kContext + ": 'amplitude' must be >= 1 ADU.");
        if (!(positiveFraction >= 0.0 && positiveFraction <= 1.0))
            throw std::invalid_argument(kContext + ": 'positiveFraction' must be in [0, 1].");

        std::optional<uint64_t> seed;
        if (j.contains("seed"))
        {
            const auto& v = j.at("seed");

            if (!v.is_number_unsigned() && !(v.is_number_integer() && v.get<int64_t>() >= 0))
                throw std::invalid_argument(kContext + ": 'seed' must be an integer >= 0.");

            seed = v.get<uint64_t>();
        }

        // ── Pixels: n distinct indices, uniform (Floyd's algorithm) ──
        const size_t nPix = width * height;
        const size_t n    = std::min(nPix, static_cast<size_t>(std::llround(fraction * static_cast<double>(nPix))));

        std::mt19937_64   rng = RTN::makeRng(seed);
        std::vector<bool> chosen(nPix, false);
        for (size_t k = nPix - n; k < nPix; ++k)
        {
            const size_t t = std::uniform_int_distribution<size_t>(0, k)(rng);
            chosen[chosen[t] ? k : t] = true;   // k itself is never chosen before this step
        }

        // ── Time constants, amplitude and sign of every RTN pixel, in increasing idx ──
        std::uniform_real_distribution<double> tau(tauMin, tauMax);
        std::uniform_real_distribution<double> amp(kAmplitudeMin, amplitude);
        std::bernoulli_distribution            positive(positiveFraction);

        std::vector<Pixel> pixels;
        pixels.reserve(n);
        for (size_t i = 0; i < nPix; ++i)
            if (chosen[i])
            {
                const double tau0 = (tauMax > tauMin) ? tau(rng) : tauMin;
                const double tau1 = (tauMax > tauMin) ? tau(rng) : tauMin;
                const double A    = (amplitude > kAmplitudeMin) ? amp(rng) : kAmplitudeMin;
                pixels.push_back({i, tau0, tau1, A, positive(rng) ? +1 : -1});
            }

        m_fraction         = fraction;
        m_tauMin           = tauMin;
        m_tauMax           = tauMax;
        m_amplitude        = amplitude;
        m_positiveFraction = positiveFraction;
        m_seed             = seed;
        m_width            = width;
        m_height           = height;
        m_pixels           = std::move(pixels);
    }

#pragma endregion

// ── Stream ────────────────────────────────────────────────────────────────────
#pragma region - Stream

    RTN::Stream::Stream(const RTN& owner)
        : m_owner(&owner), m_rng(makeRng(std::nullopt)), m_shifted(owner.m_pixels.size(), 0)
    {
        // Stationary state at the reset: shifted with the probability τ_1 / (τ_0 + τ_1).
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        for (size_t k = 0; k < m_shifted.size(); ++k)
            m_shifted[k] = (u01(m_rng) < owner.m_pixels[k].shiftedProbability()) ? 1 : 0;
    }

    void RTN::Stream::advance(double dt)
    {
        if (!(dt >= 0.0))
            throw std::invalid_argument("RTN::Stream::advance: dt must be >= 0 (got " + std::to_string(dt) + ").");
        if (dt == 0.0)
            return;

        std::uniform_real_distribution<double> u01(0.0, 1.0);
        const std::vector<Pixel>& pixels = m_owner->m_pixels;
        for (size_t k = 0; k < m_shifted.size(); ++k)
        {
            const Pixel& px     = pixels[k];
            const double lambda = 1.0 / px.tau0 + 1.0 / px.tau1;
            const double p      = px.shiftedProbability();
            const double q      = -std::expm1(-lambda * dt);   // 1 − e, accurate when λ·dt is small

            // P(shifted at t + dt): 1 − (1 − p)·(1 − e) if shifted at t, p·(1 − e) if at rest.
            const double pShifted = m_shifted[k] ? 1.0 - (1.0 - p) * q : p * q;
            m_shifted[k] = (u01(m_rng) < pShifted) ? 1 : 0;
        }
    }

    void RTN::Stream::addTo(std::valarray<double>& values) const
    {
        addTo(values, 0, 0, m_owner->m_width, m_owner->m_height);
    }

    void RTN::Stream::addTo(std::valarray<double>& values, size_t c0, size_t r0, size_t ncols, size_t nrows) const
    {
        const size_t width = m_owner->m_width;
        if (c0 + ncols > width || r0 + nrows > m_owner->m_height)
            throw std::out_of_range("RTN::Stream::addTo: the window exceeds the matrix.");
        if (values.size() < ncols * nrows)
            throw std::invalid_argument("RTN::Stream::addTo: the map is smaller than the window.");

        const std::vector<Pixel>& pixels = m_owner->m_pixels;
        for (size_t k = 0; k < m_shifted.size(); ++k)
        {
            if (!m_shifted[k])
                continue;

            const size_t col = pixels[k].idx % width, row = pixels[k].idx / width;
            if (col >= c0 && col < c0 + ncols && row >= r0 && row < r0 + nrows)
                values[(row - r0) * ncols + (col - c0)] += pixels[k].sign * pixels[k].amplitude;
        }
    }

    void RTN::Stream::writeTo(std::valarray<double>& map) const
    {
        if (map.size() != m_owner->m_width * m_owner->m_height)
            throw std::invalid_argument("RTN::Stream::writeTo: the map must have one value per pixel of the matrix.");

        const std::vector<Pixel>& pixels = m_owner->m_pixels;
        for (size_t k = 0; k < m_shifted.size(); ++k)
            map[pixels[k].idx] = m_shifted[k] ? pixels[k].sign * pixels[k].amplitude : 0.0;
    }

#pragma endregion

} // namespace nisp
