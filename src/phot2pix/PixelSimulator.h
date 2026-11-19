//
// Created by gillard on 2026-06-05.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// PixelSimulator.h
//


#ifndef P2P_PIXEL_SIMULATOR_H
#define P2P_PIXEL_SIMULATOR_H

#include "PixelConfig.h"
#include <DSTmath/DSTstatistic.h>

namespace phot2pix
{
    //@namespace phot2pix
    //@brief phot2pix namespace for the photon to pixel conversion

#pragma region - pixel class definition
    class pixsim
    {
        /*
        @class pixel
        @brief Default pixel class to simulate the pixel response to photons
        @details This class is default pixel class behaviour which used to simulate pixel response to photons. It aim to be inherited by daughter class with specific overriding of the pixel responses.
        */

        private:

#pragma region - protected member
        protected:
            // @brief Validate base + derived parameters; call from derived ctors
            //        after their own members are initialized so that the derived
            //        vtable is active and _userValidation() dispatches correctly.
            void _validateState();
            double m_darkCurrent; ///< Dark current in electrons per second
            double m_readNoise; ///< Read noise in electrons
            double m_gain; ///< Gain in electrons per ADU
            double m_quantumEfficiency; ///< Quantum efficiency (0.0 to 1.0)
            double m_KTCnoise; ///< KTC noise in electrons
            double m_baseline; ///< Baseline in ADU

            // @brief User-defined validation for derived classes
            virtual void _userValidation() {};

#pragma endregion
#pragma region - public member
        public:

#pragma region - ctor/dtor
            pixsim();
            explicit pixsim(const double& darkCurrent, const double& readNoise, const double& gain, const double& QE, const double& m_KTCnoise, const double& baseline);

            pixsim(const pixsim& other);
            pixsim& operator=(const pixsim& other);

            virtual ~pixsim();
#pragma endregion

#pragma region modifier
            void setDarkCurrent(const double& darkCurrent) { m_darkCurrent = darkCurrent; _validateState();}
            void setReadNoise(const double& readNoise)     { m_readNoise = readNoise; _validateState();}
            void setGain(const double& gain)               { m_gain = gain; _validateState(); }
            void setQuantumEfficiency(const double& QE)    { m_quantumEfficiency = QE; _validateState();}
            void setKTCnoise(const double& KTCnoise)       { m_KTCnoise = KTCnoise; _validateState();}
            void setBaseline(const double& baseline)       { m_baseline = baseline; _validateState();}
#pragma endregion

#pragma region accessor
            const double& getDarkCurrent() const { return m_darkCurrent;}
            const double& getReadNoise() const { return m_readNoise;}
            const double& getGain() const { return m_gain;}
            const double& getQuantumEfficiency() const { return m_quantumEfficiency;}
            const double& getKTCnoise() const { return m_KTCnoise;}
            const double& getBaseline() const { return m_baseline;}
#pragma endregion

#pragma region - pixel stages
            // The pixel response is split into the three physical stages of a read-out.
            // Each one can be overridden on its own; pixelresponse() chains them for a
            // single-read exposure, and a multi-read simulator (e.g. a MACC / UTR ramp)
            // calls them directly:
            //   read 1 : charge  = resetLevel() + integrate(...)   value = read(charge)
            //   read k : charge += integrate(...)                  value = read(charge)

            // @brief Reset level, in ADU: baseline + kTC noise draw. Set once per reset,
            //        common to every read until the next reset.
            virtual double resetLevel() const;

            // @brief Charge collected over exposureTime, in ADU: signal + dark current with
            //        their shot noise. This is what accumulates in the pixel between reads:
            //        no reset level, no read noise.
            virtual double integrate(const double& flux, const double& exposureTime) const;

            // @brief Runtime-parameter overload of integrate(). The base class ignores cfg.
            virtual double integrate(const double& flux, const double& exposureTime, const PixelConfig& cfg) const
            {
                return integrate(flux, exposureTime);
            }

            // @brief One read of the charge accumulated in the pixel, in ADU: charge + one
            //        independent read-noise draw N(0, readNoise)/gain. Called once per read;
            //        its result is a measurement and is never fed back into the charge.
            virtual double read(const double& charge) const;
#pragma endregion

#pragma region - pixel response
            // @details Simulate a single-read exposure: read(resetLevel() + integrate(flux, exposureTime)).
            //          Default behaviour is a purely linear response with shot noise, reset
            //          level (baseline + kTC) and read noise. Override the stages above to change
            //          one physical effect, or override pixelresponse() itself to provide a
            //          completely different single-read model.
            virtual double pixelresponse(const double& flux, const double& exposureTime) const;

            // Optional runtime-parameter overload. Derived classes can override this to
            // extract extra parameters (e.g. "dT" for persistence) from cfg.
            // The base class ignores cfg and delegates to the two-argument form, so a
            // model that only overrides pixelresponse(flux, exposureTime) is used by both.
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