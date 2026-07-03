//
// Created by gillard on 2024-06-05.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// PixelSimulator.h
//

#include "PixelSimulator.h"

#include<cstring>

namespace phot2pix
{

    void pixsim::_validateState()
    {
        if (m_darkCurrent < 0.0)
            throw std::invalid_argument("Dark current must be non-negative. D="+std::to_string(m_darkCurrent));
        
        if (m_readNoise < 0.0)
            throw std::invalid_argument("Read noise must be non-negative. RN="+std::to_string(m_readNoise));
        
        if (m_gain <= std::numeric_limits<double>::min())
            throw std::invalid_argument("Gain must be positive. G="+std::to_string(m_gain));
        
        if (m_quantumEfficiency < 0.0 || m_quantumEfficiency > 1.0 + std::numeric_limits<double>::epsilon())
            throw std::invalid_argument("Quantum efficiency must be between 0 and 1. QE="+std::to_string(m_quantumEfficiency));
        
        m_quantumEfficiency = std::clamp(m_quantumEfficiency, 0.0, 1.0);
        
        
        if (m_baseline < 0.0)
            throw std::invalid_argument("Baseline must be positive number. B="+std::to_string(m_baseline));
        
        if (m_KTCnoise < 0.0)
            throw std::invalid_argument("KTC noise must be non-negative. KTC="+std::to_string(m_KTCnoise));

        _userValidation();
    }


    pixsim::pixsim() : m_darkCurrent(0.0), m_readNoise(0.0), m_gain(1.0),m_quantumEfficiency(1.0), m_KTCnoise(0.0), m_baseline(0.0)
    {
        _validateState();
    }

    pixsim::pixsim(const double& darkCurrent, const double& readNoise, const double& gain, const double& QE, const double& KTCnoise, const double& baseline) : m_darkCurrent(darkCurrent), m_readNoise(readNoise), m_gain(gain), m_quantumEfficiency(QE), m_KTCnoise(KTCnoise), m_baseline(baseline)
    {
        _validateState();
    }

    pixsim::pixsim(const pixsim& other) : m_darkCurrent(other.m_darkCurrent), m_readNoise(other.m_readNoise), m_gain(other.m_gain), m_quantumEfficiency(other.m_quantumEfficiency), m_KTCnoise(other.m_KTCnoise), m_baseline(other.m_baseline)
    {
        _validateState();
    }

    pixsim& pixsim::operator=(const pixsim& other)
    {
        if (this != &other)
        {
            m_darkCurrent = other.m_darkCurrent;
            m_readNoise = other.m_readNoise;
            m_gain = other.m_gain;
            m_quantumEfficiency = other.m_quantumEfficiency;
            m_KTCnoise = other.m_KTCnoise;
            m_baseline = other.m_baseline;

            _validateState();
        }
        return *this;
    }

    pixsim::~pixsim() {}

    double pixsim::pixelresponse(const double& flux, const double& exposureTime) const
    {
        //Randome noise distribution:
        // Simulate the pixsim response to incoming flux over the exposure time.
        // This is a placeholder implementation; actual simulation logic would go here.
        double signal = flux * exposureTime * m_quantumEfficiency ; // Simple linear response expressed in electrons

        signal += m_darkCurrent * exposureTime; // Add dark current contribution
      
        double var = sqrt(signal + m_readNoise * m_readNoise); // Total noise (Poisson + read noise)

        return ((var > 0.0) ? DST::Math::sample_normal_distribution(signal, var) / m_gain : 0.0) + // Convert to ADU using gain
               ((m_KTCnoise/m_gain > 0.0) ? DST::Math::sample_normal_distribution(m_baseline, m_KTCnoise/m_gain) : m_baseline); // Add KTC noise contribution
    }
} // namespace phot2pix