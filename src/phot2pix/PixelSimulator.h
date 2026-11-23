//
// Created by gillard on 2026-06-05.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// PixelSimulator.h
//
/// @file PixelSimulator.h
/// @brief phot2pix::pixsim, the three stages of a pixel: reset, integration and read.


#ifndef P2P_PIXEL_SIMULATOR_H
#define P2P_PIXEL_SIMULATOR_H

#include "PixelConfig.h"
#include <DSTmath/DSTstatistic.h>

namespace phot2pix
{
    /// @namespace phot2pix
    /// @brief Photon to pixel conversion: the instrument-agnostic core of the simulator.

#pragma region - pixel class definition
    /**
     * @class pixsim
     * @brief Default pixel: converts the photons it receives into ADU.
     * @details The response is split into three stages, each one a virtual method:
     *          resetLevel(), integrate() and read(). Derived classes override the stage whose
     *          physics they change. pixelresponse() chains the three stages for a single read.
     *          A pixel holds its parameters only, never the accumulated charge: the caller keeps it.
     */
    class pixsim
    {
        private:

#pragma region - protected member
        protected:
            /// @brief Check the parameters, then call _userValidation().
            /// @details Call it from the constructors of derived classes, after their own members are
            ///          initialized, so that the derived vtable is active and _userValidation()
            ///          dispatches correctly.
            /// @throws std::invalid_argument if a parameter is out of range: dark current, read noise,
            ///         kTC noise or baseline < 0, gain <= 0, quantum efficiency outside [0, 1].
            void _validateState();
            double m_darkCurrent; ///< Dark current in electrons per second
            double m_readNoise; ///< Read noise in electrons
            double m_gain; ///< Gain in electrons per ADU
            double m_quantumEfficiency; ///< Quantum efficiency (0.0 to 1.0)
            double m_KTCnoise; ///< KTC noise in electrons
            double m_baseline; ///< Baseline in ADU

            /// @brief Validation of the parameters added by a derived class. Does nothing by default.
            /// @details Called by _validateState(), so by every setter and every copy: keep it cheap.
            virtual void _userValidation() {};

#pragma endregion
#pragma region - public member
        public:

#pragma region - ctor/dtor
            /// @brief Ideal pixel: gain 1 e-/ADU, quantum efficiency 1, no noise, no dark current, baseline 0.
            pixsim();

            /// @brief Pixel with its six parameters.
            /// @param darkCurrent Dark current (e-/s), >= 0.
            /// @param readNoise   Read noise (e- RMS per read), >= 0.
            /// @param gain        Gain (e-/ADU), > 0.
            /// @param QE          Quantum efficiency, in [0, 1].
            /// @param m_KTCnoise  Reset (kTC) noise (e- RMS per reset), >= 0.
            /// @param baseline    Baseline (ADU), >= 0.
            /// @throws std::invalid_argument if a parameter is out of range.
            explicit pixsim(const double& darkCurrent, const double& readNoise, const double& gain, const double& QE, const double& m_KTCnoise, const double& baseline);

            /// @brief Copy constructor.
            pixsim(const pixsim& other);
            /// @brief Copy assignment.
            pixsim& operator=(const pixsim& other);

            /// @brief Destructor.
            virtual ~pixsim();
#pragma endregion

#pragma region modifier
            // Every setter throws std::invalid_argument if the value is out of range.

            /// @brief Set the dark current (e-/s), >= 0. @throws std::invalid_argument otherwise.
            void setDarkCurrent(const double& darkCurrent) { m_darkCurrent = darkCurrent; _validateState();}
            /// @brief Set the read noise (e- RMS per read), >= 0. @throws std::invalid_argument otherwise.
            void setReadNoise(const double& readNoise)     { m_readNoise = readNoise; _validateState();}
            /// @brief Set the gain (e-/ADU), > 0. @throws std::invalid_argument otherwise.
            void setGain(const double& gain)               { m_gain = gain; _validateState(); }
            /// @brief Set the quantum efficiency, in [0, 1]. @throws std::invalid_argument otherwise.
            void setQuantumEfficiency(const double& QE)    { m_quantumEfficiency = QE; _validateState();}
            /// @brief Set the reset (kTC) noise (e- RMS per reset), >= 0. @throws std::invalid_argument otherwise.
            void setKTCnoise(const double& KTCnoise)       { m_KTCnoise = KTCnoise; _validateState();}
            /// @brief Set the baseline (ADU), >= 0. @throws std::invalid_argument otherwise.
            void setBaseline(const double& baseline)       { m_baseline = baseline; _validateState();}
#pragma endregion

#pragma region accessor
            /// @brief Dark current (e-/s).
            const double& getDarkCurrent() const { return m_darkCurrent;}
            /// @brief Read noise (e- RMS per read).
            const double& getReadNoise() const { return m_readNoise;}
            /// @brief Gain (e-/ADU).
            const double& getGain() const { return m_gain;}
            /// @brief Quantum efficiency, in [0, 1].
            const double& getQuantumEfficiency() const { return m_quantumEfficiency;}
            /// @brief Reset (kTC) noise (e- RMS per reset).
            const double& getKTCnoise() const { return m_KTCnoise;}
            /// @brief Baseline (ADU).
            const double& getBaseline() const { return m_baseline;}
#pragma endregion

#pragma region - pixel stages
            // The pixel response is split into the three physical stages of a read-out.
            // Each one can be overridden on its own; pixelresponse() chains them for a
            // single-read exposure, and a multi-read simulator (e.g. a MACC / UTR ramp)
            // calls them directly:
            //   read 1 : charge  = resetLevel() + integrate(...)   value = read(charge)
            //   read k : charge += integrate(...)                  value = read(charge)

            /// @brief Reset level: baseline + kTC noise draw.
            /// @details Drawn once per reset, common to every read until the next reset.
            /// @return Reset level (ADU).
            virtual double resetLevel() const;

            /// @brief Charge collected over exposureTime: signal + dark current with their shot noise.
            /// @details This is what accumulates in the pixel between reads: no reset level, no read
            ///          noise. The shot noise is the Gaussian approximation of Poisson.
            /// @param flux         Photons/s received by the pixel.
            /// @param exposureTime Integration time (s).
            /// @return Collected charge (ADU).
            virtual double integrate(const double& flux, const double& exposureTime) const;

            /// @brief integrate() with runtime parameters. The base class ignores cfg.
            /// @param flux         Photons/s received by the pixel.
            /// @param exposureTime Integration time (s).
            /// @param cfg          Runtime parameters, for derived classes.
            /// @return Collected charge (ADU).
            virtual double integrate(const double& flux, const double& exposureTime, const PixelConfig& cfg) const
            {
                return integrate(flux, exposureTime);
            }

            /// @brief One read of the accumulated charge: charge + N(0, readNoise)/gain.
            /// @details Called once per read, with an independent draw each time. The result is a
            ///          measurement and is never fed back into the charge.
            /// @param charge Accumulated charge, reset level included (ADU).
            /// @return Measured value (ADU).
            virtual double read(const double& charge) const;
#pragma endregion

#pragma region - pixel response
            /// @brief Single-read exposure: read(resetLevel() + integrate(flux, exposureTime)).
            /// @details By default a purely linear response with shot noise, reset level (baseline +
            ///          kTC) and read noise. Override the stages above to change one physical effect,
            ///          or pixelresponse() itself for a completely different single-read model.
            /// @param flux         Photons/s received by the pixel.
            /// @param exposureTime Exposure time (s).
            /// @return Read value (ADU).
            virtual double pixelresponse(const double& flux, const double& exposureTime) const;

            /// @brief pixelresponse() with runtime parameters (e.g. "dT" for persistence).
            /// @details The base class ignores cfg and delegates to the two-argument form, so a model
            ///          that only overrides pixelresponse(flux, exposureTime) is used by both.
            /// @param flux         Photons/s received by the pixel.
            /// @param exposureTime Exposure time (s).
            /// @param cfg          Runtime parameters, for derived classes.
            /// @return Read value (ADU).
            virtual double pixelresponse(const double& flux, const double& exposureTime, const PixelConfig& cfg) const
            {
                return pixelresponse(flux, exposureTime);
            }
#pragma endregion
#pragma endregion
    };
#pragma endregion
    
} // namespace phot2pix



#endif