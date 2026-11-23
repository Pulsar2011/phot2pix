//
// Created by gillard on 2026-11-21.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// Geometry.h
//
/// @file Geometry.h
/// @brief Physical detector, and placement of the simulated matrix on it (phot2pix::Geometry).
///
/// The physical detector is width × height pixels. A simulated matrix of W × H pixels
/// is a part of that detector, placed at (x0, y0): matrix pixel (x, y) is physical
/// pixel (x0 + x, y0 + y). Unless "detectorOrigin" is given, (x0, y0) = (0, 0): the
/// first pixel of the matrix is the first pixel of the detector.
///
/// The geometry maps the per-pixel inputs onto the matrix (extract()): a map is either
/// matrix-sized (W × H, used as is) or detector-sized (width × height, the section
/// [x0, x0+W) × [y0, y0+H) is used). A full-detector calibration map can then be used
/// to simulate any part of the detector.
///
/// By default (DetectorSimulator without an explicit geometry) the detector is the
/// matrix. A detector family derives from Geometry to add its own layout (reference
/// pixels, …) by overriding configure(), isReference() and clone(): see nisp::Geometry.
///
/// Configuration (all optional; absent keys leave the current value unchanged):
/// @verbatim
///   "detectorSize"  : [width, height]   physical width and height
///   "detectorOrigin": [x0, y0]          placement of the matrix (default: [0, 0])
/// @endverbatim

#ifndef P2P_GEOMETRY_H
#define P2P_GEOMETRY_H

#include "PixelConfig.h"

#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <valarray>

namespace phot2pix
{

#pragma region - Geometry

    /// @brief Physical detector and placement of the simulated matrix on it (see Geometry.h).
    class Geometry
    {
    public:
        size_t width;    ///< physical detector columns
        size_t height;   ///< physical detector rows
        std::optional<std::pair<size_t, size_t>> origin;   ///< explicit placement of the matrix; (0, 0) when empty

        /// @brief Physical detector of w × h pixels.
        /// @param w Columns, > 0.
        /// @param h Rows, > 0.
        /// @throws std::invalid_argument if w or h is 0.
        Geometry(size_t w, size_t h) : width(w), height(h)
        {
            if (w == 0 || h == 0)
                throw std::invalid_argument("phot2pix::Geometry: the detector size must be > 0.");
        }

        /// @brief Destructor.
        virtual ~Geometry() = default;

        /// @brief Copy of this geometry with its dynamic type (configure a copy, then commit it).
        virtual std::shared_ptr<Geometry> clone() const { return std::make_shared<Geometry>(*this); }

        /// @brief Read "detectorSize" and "detectorOrigin" from cfg (see Geometry.h).
        /// @param cfg Configuration; absent keys leave the current value unchanged.
        /// @throws std::invalid_argument on a malformed value; the geometry is then unchanged.
        virtual void configure(const PixelConfig& cfg)
        {
            Geometry g = *this;

            if (cfg.has("detectorSize"))
            {
                const auto [w, h] = _pair(cfg, "detectorSize");
                if (w == 0 || h == 0)
                    throw std::invalid_argument("phot2pix::Geometry: 'detectorSize' must be > 0.");
                g.width  = w;
                g.height = h;
            }

            if (cfg.has("detectorOrigin"))
                g.origin = _pair(cfg, "detectorOrigin");

            width  = g.width;
            height = g.height;
            origin = g.origin;
        }

        /// @brief True if physical pixel (px, py) is a reference pixel (no photodiode). None by default.
        virtual bool isReference(size_t /*px*/, size_t /*py*/) const { return false; }

        /// @brief Origin (x0, y0) of a W × H matrix on the detector: "detectorOrigin" if given, otherwise (0, 0).
        /// @param W Matrix columns.
        /// @param H Matrix rows.
        /// @throws std::invalid_argument if the matrix does not fit.
        std::pair<size_t, size_t> place(size_t W, size_t H) const
        {
            const auto [x0, y0] = origin.value_or(std::pair<size_t, size_t>{0, 0});
            if (x0 + W > width || y0 + H > height)
                throw std::invalid_argument("phot2pix::Geometry: a " + std::to_string(W) + " × " + std::to_string(H) +
                                            " matrix at detectorOrigin [" + std::to_string(x0) + ", " + std::to_string(y0) +
                                            "] does not fit on the " + std::to_string(width) + " × " +
                                            std::to_string(height) + " detector (see \"detectorSize\").");
            return {x0, y0};
        }

        /// @brief Values of the W × H matrix from a mapW × mapH map (row-major).
        /// @details A W × H map is used as is (matrix coordinates); a width × height map gives its
        ///          section [x0, x0+W) × [y0, y0+H) (detector coordinates).
        /// @param map  Values of the map, row-major.
        /// @param mapW Map columns.
        /// @param mapH Map rows.
        /// @param W    Matrix columns.
        /// @param H    Matrix rows.
        /// @param what Name of the map in the error messages.
        /// @return W × H values, row-major.
        /// @throws std::invalid_argument if the map has any other size.
        std::valarray<double> extract(const std::valarray<double>& map, size_t mapW, size_t mapH,
                                      size_t W, size_t H, const std::string& what) const
        {
            if (map.size() != mapW * mapH)
                throw std::invalid_argument(what + ": " + std::to_string(map.size()) + " values for a " +
                                            std::to_string(mapW) + " × " + std::to_string(mapH) + " image.");

            if (mapW == W && mapH == H)
                return map;

            if (mapW == width && mapH == height)
            {
                const auto [x0, y0] = place(W, H);
                return map[std::gslice(y0 * mapW + x0, std::valarray<size_t>{H, W}, std::valarray<size_t>{mapW, 1})];
            }

            throw std::invalid_argument(what + ": a " + std::to_string(mapW) + " × " + std::to_string(mapH) +
                                        " map matches neither the " + std::to_string(W) + " × " + std::to_string(H) +
                                        " matrix nor the " + std::to_string(width) + " × " + std::to_string(height) +
                                        " detector.");
        }

    protected:
        /// @brief Value of a key holding [a, b], two non-negative integers.
        /// @throws std::invalid_argument if the value has another form.
        static std::pair<size_t, size_t> _pair(const PixelConfig& cfg, const char* key)
        {
            const std::string what = std::string("phot2pix::Geometry: '") + key + "'";

            const auto& j = cfg.getJson(key);   // throws std::invalid_argument if not an array/object
            if (!j.is_array() || j.size() != 2 || !j[0].is_number() || !j[1].is_number())
                throw std::invalid_argument(what + " must be an array of two numbers.");

            const double a = j[0].get<double>(), b = j[1].get<double>();
            if (!(a >= 0.0) || !(b >= 0.0) || a != std::floor(a) || b != std::floor(b))
                throw std::invalid_argument(what + " must hold non-negative integers.");

            return {static_cast<size_t>(a), static_cast<size_t>(b)};
        }
    };

#pragma endregion

} // namespace phot2pix

#endif // P2P_GEOMETRY_H
