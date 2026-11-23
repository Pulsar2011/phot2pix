//
// Created by gillard on 2026-11-22.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispScenerySequencer.h
//
/// @file nispScenerySequencer.h
/// @brief Sequence of NISP exposures (nisp::ScenerySequencer).
///
/// See phot2pix::ScenerySequencer for the "flux" array. An entry gives its timing either as an
/// exposure time or as MACC parameters, and may give its own readout mode:
/// @verbatim
///   "flux": [ {"file": "<path>", "hdu": h, "exposureMode": "MACC", "macc": [4, 16, 4], "output": "<path>"},
///             {"file": "<path>", "hdu": h, "exposureMode": "FLUX", "exposureTime": 87.2}, ... ]
///
///   "exposureTime" or "macc" = [ng, nr, nd] (or "ng", "nr", "nd")   exactly one (see nisp::readout)
///   "exposureMode"   MACC, UTR, CDS, DPU or FLUX, for this exposure only; default: the SCA's
///                    mode (the top-level "exposureMode")
/// @endverbatim
/// run(sca) simulates every exposure with SCA::expose() on the same SCA, so the persistence of
/// an exposure reaches the following ones.

#ifndef _NISP_SCENERY_SEQUENCER_H_
#define _NISP_SCENERY_SEQUENCER_H_

#include "nispSCA.h"

#include <phot2pix/phot2pix.h>

namespace nisp
{
    /// @brief Sequence of NISP exposures, each with its own timing and readout mode
    ///        (see nispScenerySequencer.h).
    class ScenerySequencer : public phot2pix::ScenerySequencer
    {
    public:
        /// @brief Empty sequence: call initialize().
        ScenerySequencer() = default;
        /// @brief Destructor.
        virtual ~ScenerySequencer() = default;

        /// @brief Simulate the exposures in array order.
        /// @details sca.expose(flux, readout::fromConfig(entry)) for each, saved when it has an "output".
        /// @param sca Detector, already initialized.
        /// @throws std::logic_error if initialize() has not been called.
        void run(SCA& sca) const;

    protected:
        /// @brief Check a NISP entry: exactly one of "exposureTime" or the MACC parameters, an
        ///        optional "exposureMode" (readout::fromConfig()).
        /// @param e Entry.
        /// @throws std::invalid_argument otherwise.
        void _checkExposure(const Exposure& e) const override;
    };
}

#endif // _NISP_SCENERY_SEQUENCER_H_
