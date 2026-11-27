//
// Created by gillard on 2026-11-27.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// ChargeDeposit.h
//
/// @file ChargeDeposit.h
/// @brief Sources of charge deposited over several pixels (phot2pix::Deposit, phot2pix::ChargeDeposit).
///
/// A source draws, for one step of an exposure, the charge it deposits on the pixels of the physical
/// detector: cosmic rays and radioactivity (CosmicRay::CosmicDeposit, see CosmicDeposit.h), or any
/// source a detector family writes. DetectorSimulator::expose() adds the deposits of its source to the
/// simulated image. An instrument with a readout of several steps owns its source instead and adds the
/// deposits at each step (nisp::SCA).
///
/// This header depends only on Geometry.h and PixelConfig.h, so that DetectorSimulator.h, and therefore
/// phot2pix.h, can include it. DetectorEffect.h, which includes DSTmorpho.h, includes it as well.

#ifndef P2P_CHARGE_DEPOSIT_H
#define P2P_CHARGE_DEPOSIT_H

#include "Geometry.h"
#include "PixelConfig.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

namespace phot2pix
{
#pragma region - charge deposition

    /// @brief Charge deposited in one pixel of the physical detector.
    struct Deposit
    {
        size_t idx;         ///< Physical pixel, row-major: Y * geometry.width + X
        double electrons;   ///< Charge (e-), > 0
    };

    class ChargeDeposit
    {
        /**
        * @brief source of charge deposition over pixel for one step of dureation dt.
        * @details Abstract class use to draw pixels and charge deposition. Pixel are drawn from the whole detector, and charge is drawn from a distribution.
        */
#pragma region -- private
        private:
            uint64_t                      m_seed = 0;   ///< Seed; 0: none
            mutable std::atomic<uint64_t> m_calls{0};    ///< Calls of makeRng() so far
#pragma endregion
#pragma region -- public
        public:
#pragma region -- ctor/dtor
            ChargeDeposit() = default;
            ChargeDeposit(const ChargeDeposit& other):m_seed(other.m_seed),m_calls(other.m_calls.load()) {};

            ChargeDeposit& operator=(const ChargeDeposit& other)
            {
                m_seed = other.m_seed;
                m_calls.store(other.m_calls.load());
                return *this;
            };

            virtual ~ChargeDeposit() = default;
            virtual std::shared_ptr<ChargeDeposit> clone() const = 0;
#pragma endregion

#pragma region -- methods
            virtual void configure(const PixelConfig& cfg) = 0;

            /// @brief True when the source can draw, so that generate() may be called.
            /// @details Default: true. A source that needs a configuration first (CosmicDeposit) overrides it.
            virtual bool isConfigured() const { return true; }

            /// @brief Header cards that describe the inputs of the source; DetectorSimulator records them.
            /// @details Default: none.
            virtual std::vector<HeaderCard> inputCards() const { return {}; }

            /// @brief Append the deposits of a step of dt seconds on the physical detector.
            /// @param dt       Duration of the step (s).
            /// @param geometry Physical detector: its pixels are the only ones that receive charge.
            /// @param out      Deposits are appended to it, repeated pixels included.
            /// @throws std::logic_error if the source is not configured.
            virtual void generate(double dt, const Geometry& geometry, std::vector<Deposit>& out) const = 0;
#pragma endregion

#pragma region -- modifiers
            /// @brief Seed of the random streams; 0 (default): not reproducible.
            /// @details With a seed, the n-th call of generate() always draws the same stream, whatever
            ///          the threads: a sequence of frames is reproducible, as long as the calls are
            ///          made in the same order.
            void setSeed(uint64_t seed) { m_seed = seed; m_calls.store(0); }
#pragma endregion

#pragma region -- accessors
            /// @brief Seed of the random streams.
            uint64_t seed() const { return m_seed; }
#pragma endregion

#pragma region - protected
        protected:

            std::mt19937_64 makeRng() const
            {
                if (m_seed == 0)
                {
                    std::random_device rd;
                    std::seed_seq seq{rd(), rd(), rd(), rd()};
                    return std::mt19937_64(seq);
                }

                const uint64_t call = m_calls.fetch_add(1, std::memory_order_relaxed);
                std::seed_seq seq{static_cast<uint32_t>(m_seed), static_cast<uint32_t>(m_seed >> 32),
                                  static_cast<uint32_t>(call),   static_cast<uint32_t>(call >> 32)};
                return std::mt19937_64(seq);
            }

            static void push(std::vector<Deposit>& out, const Geometry& geometry, const size_t& X, const size_t& Y, const double& e)
            {
                if (!(e > 0.0) || X < 0 || Y < 0 ||
                    X >= geometry.width || Y >= geometry.height)
                    return;

                out.push_back({Y * geometry.width + X, e});
            }

#pragma endregion

#pragma endregion
    };

    /// @brief Default charge source of a detector: an unconfigured CosmicRay::CosmicDeposit (silicon).
    /// @details Used by DetectorSimulator::_makeChargeDeposit(). It is defined in CosmicDeposit.cxx, so that
    ///          this header does not need CosmicDeposit.h.
    std::shared_ptr<ChargeDeposit> defaultChargeDeposit();
#pragma endregion
} // namespace phot2pix

#endif // P2P_CHARGE_DEPOSIT_H
