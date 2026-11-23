//
// Created by gillard on 2026-06-22.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// ScenerySequencer.h
//
/// @file ScenerySequencer.h
/// @brief A sequence of exposures read from the configuration, simulated in order on one detector
///        (phot2pix::ScenerySequencer).
///
/// @verbatim
///   "flux": [ {"file": "<path>", "hdu": <index or EXTNAME>, "exposureTime": <s>, "output": "<path>"}, ... ]
///
///   "file"          FITS image of the flux (photons/s/pixel), matrix-sized or detector-sized
///                   (see DetectorSimulator::toMatrix()).
///   "hdu"           HDU index (1 = primary) or EXTNAME of "file". Default: 1.
///   "flux"          Uniform flux (photons/s/pixel, >= 0) applied to every pixel, instead of
///                   "file" (e.g. 0 for a dark exposure). Exactly one of "file" or "flux".
///   "exposureTime"  Exposure time (s), > 0. Required by this class; a detector family may
///                   accept another timing instead (NISP: "macc", see nisp::ScenerySequencer).
///   "output"        File the exposure is written to (replaced if it exists). Without it the
///                   exposure is simulated but not saved: it still counts for the effects that
///                   remember the past exposures (persistence).
/// @endverbatim
///
/// Each exposure is saved as soon as it is simulated, so only one is in memory at a time, with
/// the cards FLUX ("file.fits[hdu]" or the uniform value), SEQIDX (its number, from 1) and SEQNUM
/// (number of exposures in the sequence).
///
/// A detector family derives from ScenerySequencer to accept its own keys (_checkExposure())
/// and adds a run() for its own detector, which calls that detector's exposure function
/// (NISP: run(SCA&) calls SCA::expose()).
//
// Later:
//   The scene will be built by generators owned by this class (sky, sky background,
//   transient sources, moving objects), whose flux is added where run() loads "flux":
//
//     Static scene    — one exposure at full-detector resolution
//
//     Variable source — repeated windowed exposures at fine sub-exposure time steps
//                       covering only the affected detector region; results are
//                       composited into the accumulated full frame
//
//     Moving object   — windowed exposures whose (r0, c0) position is updated at each
//                       time step following the object's trajectory
//
//   By restricting high-cadence re-simulation to the affected window, computation cost
//   scales with source area rather than full-detector area. This needs the time window of
//   each exposure (its duration and the delay between exposures).
//
#pragma once

#include "PixelConfig.h"
#include "DetectorSimulator.h"

#include <cstddef>
#include <memory>
#include <string>
#include <valarray>
#include <vector>

#include <DSTfits/FITSimg.h>

namespace phot2pix
{
    /// @brief Sequence of exposures of the "flux" array, simulated in order on one detector
    ///        (see ScenerySequencer.h).
    class ScenerySequencer
    {
    public:
        /// @brief One entry of the "flux" array.
        struct Exposure
        {
            size_t      index;   ///< position in the "flux" array, from 0
            PixelConfig cfg;     ///< the entry's keys, "file" / "hdu" replaced by "flux" = {"file", "hdu"}
                                 ///< (or "flux" = the uniform value): the forms
                                 ///< DetectorSimulator::loadMap(cfg, "flux") and cfg.card("flux") read
            std::string output;  ///< the entry's "output"; empty → not saved
        };

        /// @brief Empty sequence: call initialize().
        ScenerySequencer() = default;
        /// @brief Destructor.
        virtual ~ScenerySequencer() = default;

        /// @brief True if cfg holds a sequence of exposures: "flux" is a JSON array.
        static bool inConfig(const PixelConfig& cfg);

        /// @brief Read and check every entry of the "flux" array (see ScenerySequencer.h), then
        ///        _checkExposure() for each.
        /// @param cfg Configuration.
        /// @throws std::invalid_argument, naming the entry, if one is malformed; the previous
        ///         sequence is then kept.
        virtual void initialize(const PixelConfig& cfg);

        /// @brief Number of exposures.
        size_t size() const { return m_exposures.size(); }
        /// @brief The exposures, in array order.
        const std::vector<Exposure>& getExposures() const { return m_exposures; }

        /// @brief Set the directory the outputs are written to.
        /// @details A relative "output" is taken relative to it, an absolute one is kept. Created
        ///          when the first output is written. Default: none (the current directory).
        /// @param directory Directory; empty for the current one.
        void setOutputDirectory(const std::string& directory) { m_outputDirectory = directory; }
        /// @brief Directory the outputs are written to; empty for the current one.
        const std::string& getOutputDirectory() const { return m_outputDirectory; }

        /// @brief File an exposure is written to: its "output" in the output directory.
        /// @return Path; empty when the exposure is not saved.
        std::string outputPath(const Exposure& e) const;

        /// @brief Simulate the exposures in array order.
        /// @details simulate(flux, "exposureTime") for each, saved when it has an "output".
        ///          Not virtual: DetectorSimulator is a template. A detector family adds its own run()
        ///          overload, which hides this one (nisp::ScenerySequencer::run(SCA&)).
        /// @tparam PixelT Pixel model of the detector.
        /// @param detector Detector, already initialized.
        /// @throws std::logic_error if initialize() has not been called.
        template<typename PixelT>
        void run(const DetectorSimulator<PixelT>& detector) const
        {
            _checkInitialized();
            for (const Exposure& e : m_exposures)
            {
                const std::valarray<double> flux = detector.loadMap(e.cfg, "flux");   // later: + scene generators
                _save(e, *detector.simulate(flux, e.cfg.getAsDouble("exposureTime")));
            }
        }

    protected:
        /// @brief Check the detector-specific keys of one entry; called by initialize().
        /// @details This class requires "exposureTime".
        /// @param e Entry.
        /// @throws std::exception (any) if the entry cannot be simulated.
        virtual void _checkExposure(const Exposure& e) const;

        /// @brief Write the FLUX, SEQIDX and SEQNUM cards into image and, when e.output is set,
        ///        write image to outputPath(e).
        /// @details An existing file is replaced; its directory is created if missing.
        /// @param e     Exposure.
        /// @param image Simulated image of that exposure.
        void _save(const Exposure& e, DSL::FITScube& image) const;

        /// @brief Check that there is a sequence.
        /// @throws std::logic_error if initialize() has not been called.
        void _checkInitialized() const;

    private:
        std::vector<Exposure> m_exposures;         ///< The exposures, in array order
        std::string           m_outputDirectory;   ///< see setOutputDirectory(); empty → current directory
    };
}
