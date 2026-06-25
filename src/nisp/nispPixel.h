//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispPixel.h
//

#ifndef _NISP_HAWAII_PIXEL_H_
#define _NISP_HAWAII_PIXEL_H_

#include "nispPixel.h"
#include <PixelSimulator.h>

namespace nisp
{
    class pixel : public phot2pix::pixsim
    {
    
        protected:
            void _userValidation() override
            { }
        
        public:
            pixel() : phot2pix::pixsim() { _userValidation(); }
            pixel(const double& darkCurrent, const double& readNoise, const double& gain, const double& QE, const double& KTCnoise, const double& baseline) : phot2pix::pixsim(darkCurrent, readNoise, gain, QE, KTCnoise, baseline) { _userValidation(); }

            pixel(const pixel& other) : phot2pix::pixsim(other) {};

            virtual ~pixel() = default;


    };
}

#endif