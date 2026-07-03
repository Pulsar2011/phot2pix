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

#pragma region - pixel response
            // @bdetails Simulate the pixel response to incoming photons. In default behaviour this function emulate purely linear response with shot noise and read noise. This function need to be overridden by derived classes to implement specific pixel response behaviour.
            virtual double pixelresponse(const double& flux, const double& exposureTime) const;

            // Optional runtime-parameter overload. Derived classes can override this to
            // extract extra parameters (e.g. "dT" for persistence) from cfg.
            // The base class ignores cfg and delegates to the two-argument form.
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