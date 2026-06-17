//
// Created by gillard on 2024-06-05.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// PixelSimulator.h
//


#ifndef P2P_PIXEL_SIMULATOR_H
#define P2P_PIXEL_SIMULATOR_H

#include "PixelConfig.h"

namespace phot2pix
{
    //@namespace phot2pix
    //@brief phot2pix namespace for the photon to pixel conversion


    class pixel
    {
        /*
        @class pixel
        @brief pixel class to simulate the pixel response to photons
        @details This class is used to simulate the pixel response to photons. It is used in the phot2pix project to simulate single pixel behaviour.
        */

        protected:
            double m_darkCurrent; ///< Dark current in electrons per second
            double m_readNoise; ///< Read noise in electrons
            double m_gain; ///< Gain in electrons per ADU

            virtual double _pixelresponse(const double& flux, const double& exposureTime);

        public:

            pixel();
            explicit pixel(const PixelConfig& config);
            explicit pixel(const std::string& configFile);  ///< loads JSON from path

            pixel(const pixel& other);
            pixel& operator=(const pixel& other);

            virtual ~pixel();
    };
    
} // namespace phot2pix



#endif