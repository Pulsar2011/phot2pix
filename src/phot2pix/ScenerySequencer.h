//
// Created by gillard on 2026-06-22.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// ScenerySequencer.h — NOT YET IMPLEMENTED
//
// Planned responsibility:
//   Owns the time loop and dynamic sky model for a sequence of exposures.
//   Calls DetectorSimulator::simulate() (full or windowed) at the appropriate
//   cadence per source type:
//
//     Static scene    — one simulate() call per exposure at full-detector resolution
//
//     Variable source — repeated windowed simulate() calls at fine sub-exposure
//                       time steps covering only the affected detector region;
//                       results are composited into the accumulated full frame
//
//     Moving object   — windowed simulate() whose (r0, c0) position is updated
//                       at each time step following the object's trajectory
//
//   By restricting high-cadence re-simulation to the affected window, computation
//   cost scales with source area rather than full-detector area.
//
// TODO: design Sky / SceneObject model before implementing.
//       Key questions:
//         - How is the PSF represented (analytical / numerical)?
//         - How is object motion parameterised (linear trail, ephemeris)?
//         - How are variable sources described (light curve / callable)?
//         - What is the output format (FITScube with one layer per readout)?
//
#pragma once

namespace phot2pix
{
    class ScenerySequencer;  // placeholder — implementation pending
}
