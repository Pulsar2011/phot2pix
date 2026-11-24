//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispPixel.h
//
/// @file nispPixel.h
/// @brief nisp::pixel (the H2RG pixel), its non-linearity models, and nisp::Ramp (slope estimators).

#ifndef _NISP_HAWAII_PIXEL_H_
#define _NISP_HAWAII_PIXEL_H_

#include "nispPixel.h"
#include <phot2pix/phot2pix.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>

namespace nisp
{
#pragma region - Non-linearity models
    /// @brief Models of the non-linear response of a pixel.
    namespace nonlinearity
    {
        constexpr size_t nParam = 3;   ///< Parameters held by each pixel ("nl0", "nl1", "nl2"); 0 when not given

        /// @brief Response of a pixel: measured signal (e-) for a stored charge Q (e-), 0 < Q <= full well.
        /// @details It must tend to Q at low charge (otherwise it changes the gain) and increase up to the
        ///          full well. Called from several threads: it must not throw nor modify what it captures.
        using functor = std::function<double(std::span<const float> param, double Q)>;

        /// @brief Polynomial response Q·(1 + c2·Q + c3·Q²).
        /// @param param {c2 (1/e-), c3 (1/e-²)}; the third parameter is not used.
        /// @param Q     Stored charge (e-).
        /// @return Measured signal (e-).
        double polynomial(std::span<const float> param, double Q);

        /// @brief Model named by an "nlModel" value ("polynomial").
        /// @throws std::invalid_argument for any other name.
        functor fromString(const std::string& name);
    }
#pragma endregion

    /// @brief NISP H2RG pixel: the phot2pix::pixsim stages, with a full well and a non-linear response.
    /// @details
    ///   - resetLevel(): baseline + kTC, once per reset;
    ///   - integrate(): signal + dark with shot noise, accumulates between reads;
    ///   - read(charge, nl): full well and non-linearity of the accumulated charge, then
    ///     pixsim::read(): + read noise, once per read (never accumulated).
    ///
    /// nisp::SCA chains them for every readout mode (MACC, UTR, CDS, FLUX).
    /// Override them here to add NISP-specific pixel physics.
    class pixel : public phot2pix::pixsim
    {
        private:
            float                                   m_fullWell = std::numeric_limits<float>::infinity(); ///< Full well (e-); infinity: none
            std::array<float, nonlinearity::nParam> m_nlParam{};                                         ///< Parameters of the non-linearity model

        protected:
            /// @brief Nothing to validate here: setFullWell() and SCA check the added parameters.
            void _userValidation() override
            { }

        public:
            /// @brief Ideal pixel, without full well.
            pixel() : phot2pix::pixsim() { _userValidation(); }
            /// @brief Pixel with the six phot2pix::pixsim parameters, without full well.
            /// @param darkCurrent Dark current (e-/s), >= 0.
            /// @param readNoise   Read noise (e- RMS per read), >= 0.
            /// @param gain        Gain (e-/ADU), > 0.
            /// @param QE          Quantum efficiency, in [0, 1].
            /// @param KTCnoise    Reset (kTC) noise (e- RMS per reset), >= 0.
            /// @param baseline    Baseline (ADU), >= 0.
            /// @throws std::invalid_argument if a parameter is out of range.
            pixel(const double& darkCurrent, const double& readNoise, const double& gain, const double& QE, const double& KTCnoise, const double& baseline) : phot2pix::pixsim(darkCurrent, readNoise, gain, QE, KTCnoise, baseline) { _userValidation(); }

            /// @brief Copy constructor: full well and non-linearity parameters included.
            pixel(const pixel& other) = default;

            /// @brief Destructor.
            virtual ~pixel() = default;

            /// @brief Set the full well.
            /// @param fullWell Full well (e-), > 0; infinity for none.
            /// @throws std::invalid_argument if it is not > 0.
            void setFullWell(float fullWell)
            {
                if (!(fullWell > 0.0f))   // also rejects NaN
                    throw std::invalid_argument("Full well must be positive. FW=" + std::to_string(fullWell));
                m_fullWell = fullWell;
            }
            /// @brief Set parameter k of the non-linearity model. @throws std::out_of_range if k >= nonlinearity::nParam.
            void setNlParam(size_t k, float value) { m_nlParam.at(k) = value; }

            /// @brief Full well (e-); infinity: none.
            float getFullWell() const { return m_fullWell; }
            /// @brief Parameter k of the non-linearity model. @throws std::out_of_range if k >= nonlinearity::nParam.
            float getNlParam(size_t k) const { return m_nlParam.at(k); }

            /// @brief Charge the pixel holds (e-) when Q electrons were generated.
            double storedCharge(double Q) const { return std::min(Q, static_cast<double>(m_fullWell)); }

            /// @brief Measured signal (e-) for a stored charge Q (e-).
            /// @details Linear if nl is empty, and for Q <= 0 (kTC noise), which a model does not have
            ///          to support.
            /// @param Q  Stored charge (e-).
            /// @param nl Non-linearity model.
            double nlResponse(double Q, const nonlinearity::functor& nl) const
            {
                return (nl && Q > 0.0) ? nl(m_nlParam, Q) : Q;
            }

            using phot2pix::pixsim::read;   ///< read(charge): read noise only, no full well, linear

            /// @brief One read of the accumulated charge: full well, then non-linearity, then read noise.
            /// @details The baseline is an electronic offset: it does not fill the well. The charge
            ///          itself is not modified. Without full well, this is read(charge).
            /// @param charge Total accumulated charge (ADU, linear, reset level included).
            /// @param nl     Non-linearity model; empty for a linear response.
            /// @return Measured value (ADU).
            double read(const double& charge, const nonlinearity::functor& nl) const
            {
                if (std::isinf(m_fullWell))
                    return read(charge);

                const double Q = storedCharge((charge - m_baseline) * m_gain);   // ADU → e-
                return read(m_baseline + nlResponse(Q, nl) / m_gain);            // e- → ADU, + read noise
            }
    };


#pragma region - Ramp processing

    /// @brief Reads of one or several pixels, as a flat array.
    typedef std::valarray<double> rawLines;

    /**
     * @class Ramp
     * @brief The reads of one pixel: co-addition into groups, and CDS, Fowler and slope estimators.
     * @details A view on the reads of one pixel (UTR or MACC), selected in a larger array by a
     *          std::slice. It does not own the data.
     **/
    class Ramp
    {
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
        /// @brief Destructor; the data are not freed.
        ~Ramp(){};

#pragma endregion
#pragma region --  modifier

        /// @brief Set the number of reads per coadded group. Must match the MACC pattern used in coadd() / getSlope().
        inline void reads_per_group(const int16_t &n) {nreads = static_cast<size_t>(n);}
        /// @brief Set the number of dropped frames per group. Must match the MACC pattern used in coadd() / getSlope().
        inline void drops_per_group(const int16_t &n) {ndrops = static_cast<size_t>(n);}

        inline const size_t& get_reads_per_group() const { return nreads;} ///< Return the current reads-per-group setting.
        inline const size_t& get_drops_per_group() const { return ndrops;} ///< Return the current drops-per-group setting.
        inline const std::slice& get_slice      () const { return fstr;}   ///< Return the slice that addresses this pixel's frames in the data array.
        
#pragma endregion

#pragma region --   accessor
        /// @brief Copy of the raw ramp data (the frames selected by \c fstr from the parent array).
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
               Deleting the Ramp does not free its group array (Ramp does not own its data):
               for a per-pixel loop use the buffer overload below.
         **/
        Ramp& coadd(const size_t& ngr, const size_t& nreads, const size_t& ndrops) const;

        /**
         @brief Same co-addition as coadd(ngr, nreads, ndrops), without any allocation.
         @details The \p ngr group averages are written into \p groups (resized only if it
                  is smaller than \p ngr) and the returned Ramp points to it, with
                  reads_per_group / drops_per_group pre-set. Meant for per-pixel loops:
                  the caller reuses one buffer for every pixel.
         @param groups Output buffer for the group averages; must outlive the returned Ramp.
         @param ngr    Number of output groups.
         @param nreads Number of reads to average per group.
         @param ndrops Number of frames to discard between groups.
         @return Ramp over \p groups.
         **/
        Ramp coadd(rawLines& groups, const size_t& ngr, const size_t& nreads, const size_t& ndrops) const;
        
        
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
         @param flux Upon return flux value in pix_unit/Group, where pix_unit is the pixel unit, either e- when the gain is 1 or ADU if it is different from 1.
         @param chi2 Upon return reduced chi2 value: chi2 / (number of groups - 2), not defined for 2 groups.
         @param err2 Upon return squarred flux error estimation value, in the same unit as flux at the power of 2.
         @throws std::runtime_error if the ramp has less than 2 groups.
         @note If the ramp is not already in the form of coadded groups, number of reads, number of drops and number of groups have to be tunned to this ramp otherwise the flux estimator will be biased.
         */
        void getSlope(double &flux, double &chi2, double& err2) const;
        
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