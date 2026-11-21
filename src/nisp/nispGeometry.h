//
// Created by gillard on 2026-11-21.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// nispGeometry.h
//
// Physical layout of the NISP HAWAII-2RG (see phot2pix::Geometry for the placement of
// the simulated matrix and the mapping of the per-pixel inputs).
//
// The physical detector is 2048 × 2048 pixels, framed by refBorder reference pixels on
// each side (4). It is read through output channels of colPerChannel columns (64,
// configured with the 1/f noise, see nispOneOverFNoise.h).
//
// The matrix starts at the first detector pixel unless "detectorOrigin" is given: a
// 2048 × 2048 matrix is the full detector, a 2040 × 2040 one at [4, 4] its active area,
// and a 2044-wide one at x0 = 2 keeps 2 reference columns on each side.
//
// Everything that depends on the detector (reference pixels, output channels, readout
// timing) is derived from physical coordinates, so it does not depend on the size of
// the matrix: a smaller matrix is simply a smaller part of the same detector.
//
// Configuration (all optional; absent keys leave the current value unchanged):
//   "detectorSize"  : [2048, 2048]   physical width and height
//   "detectorOrigin": [x0, y0]       placement of the matrix (default: [0, 0])
//   "refBorder"     : 4              reference pixels on each side (0: none)

#ifndef _NISP_GEOMETRY_H_
#define _NISP_GEOMETRY_H_

#include <phot2pix/Geometry.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>

namespace nisp
{

#pragma region - Geometry

    class Geometry : public phot2pix::Geometry
    {
    public:
        size_t refBorder = 4;   ///< reference pixels on each side of the physical detector

        Geometry() : phot2pix::Geometry(2048, 2048) {}

        std::shared_ptr<phot2pix::Geometry> clone() const override { return std::make_shared<Geometry>(*this); }

        /// Read "detectorSize", "detectorOrigin" (phot2pix::Geometry) and "refBorder" from cfg.
        /// Throws std::invalid_argument on a malformed value; the geometry is then unchanged.
        void configure(const phot2pix::PixelConfig& cfg) override
        {
            Geometry g = *this;
            g.phot2pix::Geometry::configure(cfg);

            if (cfg.has("refBorder"))
            {
                const double b = cfg.getAsDouble("refBorder");
                if (!(b >= 0.0) || b != std::floor(b))
                    throw std::invalid_argument("nisp::Geometry: 'refBorder' must be a non-negative integer.");
                g.refBorder = static_cast<size_t>(b);
            }

            if (2 * g.refBorder > std::min(g.width, g.height))
                throw std::invalid_argument("nisp::Geometry: 'refBorder' (" + std::to_string(g.refBorder) +
                                            ") does not fit on a " + std::to_string(g.width) + " × " +
                                            std::to_string(g.height) + " detector.");

            *this = g;
        }

        /// True if physical pixel (px, py) is within refBorder of the detector edge.
        bool isReference(size_t px, size_t py) const override
        {
            return px < refBorder || px + refBorder >= width ||
                   py < refBorder || py + refBorder >= height;
        }
    };

#pragma endregion

} // namespace nisp

#endif // _NISP_GEOMETRY_H_
