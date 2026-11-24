//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispPixel.cxx
//

#include "nispPixel.h"

namespace nisp
{
    // nisp::pixel uses phot2pix::pixsim::resetLevel() / integrate() as they are; read(charge, nl)
    // applies the full well and the non-linearity before phot2pix::pixsim::read().

#pragma region - Non-linearity models
    // To add a model: write its function here, declare it in nispPixel.h and name it in fromString().
    namespace nonlinearity
    {

        double polynomial(std::span<const float> param, double Q)
        {
            // NOTE one din't use the DST::Math::polynomial::polyev() because it is not constexpr and we want to be able to use it in a constexpr context.
            // It also not obtimized for large number of pixel and will significantly slow down the simulation.
            return Q * (1.0 + Q * (param[0] + Q * param[1]));
        }

        functor fromString(const std::string& name)
        {
            if (name == "polynomial") return polynomial;

            throw std::invalid_argument("unknown nlModel '" + name + "' (expected polynomial).");
        }
    }
#pragma endregion

#pragma region - Ramp class implementation
#pragma region -- methods
    // Correlated Double Sampling estimator. Two branches:
    //
    //  Standard CDS (ndrops ≤ 1, UTR or single-drop):
    //    Computes (frame[n+1] − frame[n]) for every consecutive pair across the whole
    //    ramp, then averages all differences for cds and computes their std-dev for noise.
    //
    //  MACC CDS (nreads > 1):
    //    Within each group of nreads frames, computes (frame[k+1] − frame[k]) for
    //    k = 0 … nreads−2, collects all differences across all groups, then derives
    //    mean and std-dev. This avoids the group-boundary discontinuities that would
    //    bias a naïve all-frames CDS in a MACC ramp.
    void Ramp::getCDS(double& cds, double& noise) const
    {
        if(fstr.size() < 2)
            throw std::runtime_error("[Ramp::getCDS] Incorrect dimension ["+std::to_string(__LINE__)+"]");
        
        if(ndrops <= 1)
        {
            nisp::rawLines _cds1 = (*fdata)[std::slice(fstr.start()              , fstr.size()-1, fstr.stride())];
            nisp::rawLines _cds2 = (*fdata)[std::slice(fstr.start()+fstr.stride(), fstr.size()-1, fstr.stride())];
        
            nisp::rawLines _CDS   = (_cds2-_cds1);
            nisp::rawLines _CDSSQ = _CDS*_CDS;
        
            cds = _CDS.sum() / static_cast<double>(_CDS.size());
            _CDSSQ -= cds*cds;
        
            noise = sqrt(_CDSSQ.sum() / static_cast<double>(_CDS.size() - 1));
        
            _cds1 .resize(0);
            _cds2 .resize(0);
            _CDS  .resize(0);
            _CDSSQ.resize(0);
        }
        else if(nreads > 1)
        {
            std::vector<double> _CDS;

            for(size_t iS = 0; iS < fstr.size()/nreads; iS++)
            {
                nisp::rawLines _cds1 = (*fdata)[std::slice(fstr.start()+iS              , nreads-1, fstr.stride())];
                nisp::rawLines _cds2 = (*fdata)[std::slice(fstr.start()+iS+fstr.stride(), nreads-1, fstr.stride())];
        
                nisp::rawLines this_cds   = (_cds2-_cds1);
                for(size_t ik=0; ik<this_cds.size(); ik++)
                    _CDS.push_back(this_cds[ik]);
                
                this_cds.resize(0);
                _cds1 .resize(0);
                _cds2 .resize(0);
            }
                
            nisp::rawLines _CDS_(_CDS.data(),_CDS.size()); 
            nisp::rawLines _CDSSQ = _CDS_*_CDS_;

            cds = _CDS_.sum() / static_cast<double>(_CDS_.size());
            _CDSSQ -= cds*cds;
        
            noise = sqrt(_CDSSQ.sum() / static_cast<double>(_CDS.size() - 1));
        
            _CDS_ .resize(0);
            _CDSSQ.resize(0);
        }
    }
    
    // Build the MACC co-added ramp from a UTR input.
    // For each output group gr = 0 … ngr−1, the starting frame index within fdata is
    //   f0 = fstr.start() + gr * (nreads + ndrops) * fstr.stride()
    // The nreads frames starting at f0 are summed and divided by nreads to produce
    // the group value. The ndrops frames that follow are skipped (charge still
    // accumulated in the detector, but not read, so they do not appear in the slice).
    // The result is stored in a freshly heap-allocated rawLines array of size ngr.
    Ramp& Ramp::coadd(const size_t& ngr, const size_t& nreads, const size_t& ndrops) const
    {
        nisp::rawLines* COADD = new nisp::rawLines(0.0, static_cast<size_t>(ngr));
        return *(new Ramp(coadd(*COADD, ngr, nreads, ndrops)));
    }

    // Same co-addition into a caller-owned buffer: no allocation (unless groups is too
    // small), no temporary slice arrays. Used once per pixel by SCA::expose_dpu().
    Ramp Ramp::coadd(rawLines& groups, const size_t& ngr, const size_t& nreads, const size_t& ndrops) const
    {
        if (groups.size() < ngr)
            groups.resize(ngr);

        const size_t step = fstr.stride();
        for (size_t gr = 0; gr < ngr; ++gr)
        {
            const size_t f0 = fstr.start() + gr * (nreads + ndrops) * step;
            double sum = 0.0;
            for (size_t r = 0; r < nreads; ++r)
                sum += (*fdata)[f0 + r * step];
            groups[gr] = sum / static_cast<double>(nreads);
        }

        Ramp macc(&groups, std::slice(0, ngr, 1), fsigma, fgain);
        macc.reads_per_group(static_cast<int16_t>(nreads));
        macc.drops_per_group(static_cast<int16_t>(ndrops));
        return macc;
    }
    
    // Kubik et al. 2016 (PASP 128:104504) slope estimator for a MACC ramp.
    //
    // Notation follows the paper (ng = number of groups, nf = reads/group, nd = drops/group):
    //   alpha = (1 − nf²) / (3·nf·(nf+nd))          [eq. A2, correlation correction]
    //   beta  = 2·σ²·fe / (nf·(alpha+1))            [effective noise per CDS difference]
    //   gamma = 2·fe / (alpha+1)                    [flux-to-CDS conversion factor]
    //   fe    = gain                                [conversion gain, e⁻/ADU]
    //
    // The CDS differences d[i] = group[i+1] − group[i] (i = 0 … ng−2) are formed.
    // A pseudo-flux estimate is derived from the mean of (d+β)² (eq. 20), then
    // corrected for the non-linear bias to obtain the final flux (eq. 20, full form).
    // The reduced χ² (eq. 19) uses a Fowler (first-last) estimate as reference.
    // The flux error (eq. 12) is propagated from the CDS scatter.
    void Ramp::getSlope(double& flux, double& chi2, double &err2) const
    {
        //std::cout<<"read noise :"<<fsigma<<std::endl;
        //std::cout<<"gain       :"<<fgain<<std::endl;
        
        //std::cout<<"read noise :"<<ron<<" ADU"<<std::endl;
        
        double ng = static_cast<double>( fstr.size() );
        double nf = static_cast<double>( nreads );
        double nd = static_cast<double>( ndrops );
        double fe = fgain;

        double alpha = (1. - nf*nf) / (3. * nf * (nf + nd));
        double beta  = 2. * (fsigma*fsigma) * fe / (nf * (alpha + 1.));
        double gamma = (2. * fe) / (alpha + 1.);
        
        //std::cout<<"alpha :"<<alpha<<std::endl;
        //std::cout<<"beta  :"<<beta<<std::endl;
        //std::cout<<"gamma :"<<gamma<<std::endl;

        if(fstr.size() < 2)
            throw std::runtime_error("[Ramp::getSlope] At least 2 groups are needed ["+std::to_string(__LINE__)+"]");

        // Mean of (d + beta)² over the CDS differences d[k] = group[k+1] − group[k],
        // accumulated in one pass: no temporary arrays (called once per pixel).
        const size_t s0   = fstr.start();
        const size_t st   = fstr.stride();
        const double nCds = static_cast<double>(fstr.size() - 1);

        double ffit = 0.0;
        for(size_t k = 0; k + 1 < fstr.size(); ++k)
        {
            const double d = (*fdata)[s0 + (k + 1) * st] - (*fdata)[s0 + k * st] + beta;
            ffit += d * d;
        }
        ffit /= nCds;

        double lf;
        this->getFowler(lf);

        double pseudo_flux = sqrt(ffit) - beta;

        flux = 1./gamma * (sqrt(1. + gamma*gamma*ffit ) - 1) - beta;
    //        flux /= (nf + nd);

        chi2 = gamma * (nCds*pseudo_flux - lf);
        chi2 /= ( ng - 2. ); // <- Get Normalized chi2

        err2 = 2*(pseudo_flux+beta)*(pseudo_flux+beta)/nCds;
        err2 /= (1.+gamma*(pseudo_flux+beta));
    }

    // Fowler sampling: simplest possible slope estimate — last frame minus first frame.
    // Used as a reference value in getSlope() to anchor the reduced χ² calculation.
    // Returns the result in the same units as the pixel data (e⁻ or ADU depending on fgain).
    void Ramp::getFowler(double &flux) const
    {
        flux = (*fdata)[fstr.start() + (fstr.size() - 1) * fstr.stride()] - (*fdata)[fstr.start()];
    }

#pragma endregion
#pragma endregion
}
