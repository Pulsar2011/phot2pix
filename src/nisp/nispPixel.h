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
#include <phot2pix/phot2pix.h>

namespace nisp
{
    class pixel : public phot2pix::pixsim
    {
    
        protected:
            void _userValidation() override
            { }

            bool m_computeBaseline;
        
        public:
            pixel() : phot2pix::pixsim(), m_computeBaseline(true) { _userValidation(); }
            pixel(const double& darkCurrent, const double& readNoise, const double& gain, const double& QE, const double& KTCnoise, const double& baseline) : phot2pix::pixsim(darkCurrent, readNoise, gain, QE, KTCnoise, baseline), m_computeBaseline(true) { _userValidation(); }

            pixel(const pixel& other) : phot2pix::pixsim(other), m_computeBaseline(other.m_computeBaseline) {};

            virtual ~pixel() = default;

            double pixelresponse(const double& flux, const double& exposureTime) const override;
            double pixelresponse(const double& flux, const double& exposureTime, const phot2pix::PixelConfig& cfg) const override;

            void resetState() { m_computeBaseline |= true; }

    };


#pragma region - Ramp processing

    typedef std::valarray<double> rawLines;

    class Ramp
    {
    /**
     * \class Ramp DetectorData "gacen/DetectorData.h"
     * \details Helper class to process UTR or MACC
     **/
#pragma region -- private member
    private:
        rawLines* const fdata;   ///< Pointer to the raw pixel data array (not owned; shared with the parent datacube).
        const std::slice fstr;   ///< Slice into \c fdata selecting the frames that belong to this pixel's ramp.
        double fgain;            ///< Detector gain in e⁻/ADU. Set to 1 when the data is already in electrons.
        double fsigma;           ///< Read noise in the same units as the pixel data (e⁻ or ADU).
        size_t nreads;           ///< Number of non-destructive reads per coadded group (MACC parameter).
        size_t ndrops;           ///< Number of frames dropped between consecutive groups (MACC parameter).
        
    public:
#pragma endregion
#pragma region --  ctor/dtor
        /**
         @brief Default constructor
         @param cube pointer to a 3D datacude that contains ramps for 1 or more pixels
         @param str std::slice that allow to select data (frames) related to the pixel to process.
         **/
        Ramp(rawLines* const cube, const std::slice& str): fdata(cube),fstr(str),fgain(1.),fsigma(0),nreads(0),ndrops(0){};
        
        /**
         @brief Default constructor
         @param cube pointer to a 3D datacude that contains ramps for 1 or more pixels
         @param str std::slice that allow to select data (frames) related to the pixel to process.
         @param sigma Pixel read noise. Shall be given in the same units as the cube (either e- or ADU).
         **/
        Ramp(rawLines* const cube, const std::slice& str, double sigma): fdata(cube),fstr(str),fgain(1.),fsigma(sigma),nreads(0),ndrops(0){};
        
        /**
         @brief Default constructor
         @param cube pointer to a 3D datacude that contains ramps for 1 or more pixels
         @param str std::slice that allow to select data (frames) related to the pixel to process.
         @param sigma Pixel read noise. Shall be given in the same units as the cube (either electron or ADU).
         @param gain Pixel gain in e-/ADU. When the ramp is in electron, the gain shall be 1.
         **/
        Ramp(rawLines* const cube, const std::slice& str, double sigma, double gain): fdata(cube),fstr(str),fgain(gain),fsigma(sigma),nreads(0),ndrops(0){};
        
        /**
         @brief Copy constructor
         **/
        Ramp(const Ramp& r): fdata(r.fdata),fstr(r.fstr),fgain(r.fgain),fsigma(r.fsigma),nreads(r.nreads),ndrops(r.ndrops){};
        ~Ramp(){};

#pragma endregion
#pragma region --  modifier

        /// Set the number of reads per coadded group. Must match the MACC pattern used in coadd() / getSlope().
        inline void reads_per_group(const int16_t &n) {nreads = static_cast<size_t>(n);}
        /// Set the number of dropped frames per group. Must match the MACC pattern used in coadd() / getSlope().
        inline void drops_per_group(const int16_t &n) {ndrops = static_cast<size_t>(n);}

        inline const size_t& get_reads_per_group() const { return nreads;} ///< Return the current reads-per-group setting.
        inline const size_t& get_drops_per_group() const { return ndrops;} ///< Return the current drops-per-group setting.
        inline const std::slice& get_slice      () const { return fstr;}   ///< Return the slice that addresses this pixel's frames in the data array.
        
#pragma endregion

#pragma region --   accessor
        /// Return a copy of the raw ramp data (the frames selected by \c fstr from the parent array).
        inline const rawLines getRamp() const {return (rawLines) (*fdata)[fstr];}
        
#pragma endregion
#pragma region --   methods
        
        /**
         @brief Compute the co-added (MACC) version of this UTR ramp.
         @details Produces a new Ramp whose \p ngr groups each contain the average of the
                  first \p nreads frames within that group, with \p ndrops frames discarded
                  between consecutive groups. The input ramp is assumed to be in UTR order
                  (all frames recorded, no prior averaging). The returned Ramp is heap-allocated
                  and its reads_per_group / drops_per_group fields are pre-set so that
                  getSlope() can be called directly on it.
         @param ngr    Number of output groups.
         @param nreads Number of reads to average per group.
         @param ndrops Number of frames to discard between groups.
         @return Reference to a newly heap-allocated Ramp representing the MACC output.
         @note The caller is responsible for managing the lifetime of the returned object.
         **/
        Ramp& coadd(const size_t& ngr, const size_t& nreads, const size_t& ndrops) const;
        
        
        /**
         @brief Compute CDS from ramp
         @param cds Upon return average CDS value.
         @param noise Upon return CDS noise.
         @note CDS and CDS noise returned are of the same unit as the pixel content
         */
        void getCDS(double &cds, double &noise) const;
        
        /**
         @brief Compute EUCLID Flux and Chi2 from COADDED group.
         @details Compute the EUCLID Flux, Chi2 quality factor and associated error on the flux based on Kubic et al, Publications of the Astronomical Society of the Pacific, 128:104504, 2016. The flux is estimated from equation (20). The quality factor, as the reduce chi2, is estimated from equation (19) and flux error form equation (12).
         @param flux Upon return flux value in pix_unit/Goup, where pix_unit is the pixel unit, either e- when fe is 1 or ADU if fe is defferent to 1.
         @param chi2 Upon return reduced chi2 value.
         @param err Upon return flux error estimation value.
         @note If the ramp is not already in the form of coadded groups, number of reads, number of drops and number of groups have to be tunned to this ramp otherwise the flux estimator will be biased.
         */
        void getSlope(double &flux, double &chi2, double& err) const;
        
        /**
         @brief Compute Flux from a fowler
         @param flux Upon return slope value
         @note  flux returned is of the same unit as the pixel content
         */
        void getFowler(double &flux) const;
#pragma endregion
    };
#pragma endregion

}

#endif