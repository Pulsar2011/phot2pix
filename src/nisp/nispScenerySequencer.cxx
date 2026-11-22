//
// Created by gillard on 2026-11-22.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispScenerySequencer.cxx
//

#include "nispScenerySequencer.h"

#include <valarray>

namespace nisp
{
    void ScenerySequencer::run(SCA& sca) const
    {
        _checkInitialized();
        for (const Exposure& e : getExposures())
        {
            const std::valarray<double> flux = sca.loadMap(e.cfg, "flux");   // later: + scene generators
            _save(e, *sca.expose(flux, readout::fromConfig(e.cfg)));
        }
    }

    void ScenerySequencer::_checkExposure(const Exposure& e) const
    {
        (void) readout::fromConfig(e.cfg);
    }
}
