//
// Created by gillard on 2026-11-22.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// ScenerySequencer.cxx
//

#include "ScenerySequencer.h"

#include <filesystem>
#include <set>
#include <stdexcept>
#include <utility>

namespace phot2pix
{
    bool ScenerySequencer::inConfig(const PixelConfig& cfg)
    {
        if (!cfg.has("flux"))
            return false;
        try                                  { return cfg.getJson("flux").is_array(); }
        catch (const std::invalid_argument&) { return false; }   // a number or a path
    }

    void ScenerySequencer::initialize(const PixelConfig& cfg)
    {
        if (!inConfig(cfg))
            throw std::invalid_argument("ScenerySequencer: \"flux\" must be an array of exposures "
                                        "[{\"file\": ..., \"hdu\": ..., \"exposureTime\": ..., \"output\": ...}, ...].");

        const json& entries = cfg.getJson("flux");
        if (entries.empty())
            throw std::invalid_argument("ScenerySequencer: the \"flux\" array is empty.");

        // Parsed into a copy, committed only if every entry is valid.
        std::vector<Exposure> exposures;
        exposures.reserve(entries.size());
        std::set<std::string> outputs;   // normalised paths, to catch two entries writing the same file

        for (size_t i = 0; i < entries.size(); ++i)
        {
            const json&       entry = entries[i];
            const std::string where = "ScenerySequencer: 'flux'[" + std::to_string(i) + "]";

            if (!entry.is_object())
                throw std::invalid_argument(where + ": an exposure is an object {\"file\", \"hdu\", \"exposureTime\", "
                                            "\"output\"}, not " + entry.dump() + ".");

            // Flux: a FITS image {"file": path[, "hdu": h]} or a uniform value, the forms
            // DetectorSimulator::loadMap() reads.
            const bool hasFile    = entry.contains("file");
            const bool hasUniform = entry.contains("flux");
            if (hasFile && hasUniform)
                throw std::invalid_argument(where + ": give the flux either as a FITS image (\"file\") or as a uniform "
                                            "\"flux\" (ph/s/pixel), not both.");
            if (!hasFile && !hasUniform)
                throw std::invalid_argument(where + ": no flux, give a FITS image (\"file\") or a uniform \"flux\" (ph/s/pixel).");

            json source;
            if (hasUniform)
            {
                source = entry["flux"];
                if (!source.is_number() || !(source.get<double>() >= 0.0))
                    throw std::invalid_argument(where + ": \"flux\" must be a number of photons/s/pixel >= 0, not " +
                                                source.dump() + ".");
                if (entry.contains("hdu"))
                    throw std::invalid_argument(where + ": \"hdu\" goes with a FITS image (\"file\"), not with a uniform \"flux\".");
            }
            else
            {
                if (!entry["file"].is_string() || entry["file"].get<std::string>().empty())
                    throw std::invalid_argument(where + ": \"file\" must name the FITS image of the flux.");

                source = json{{"file", entry["file"]}};
                if (entry.contains("hdu"))
                {
                    const json& hdu   = entry["hdu"];
                    const bool  index = hdu.is_number_integer() && hdu.get<long long>() >= 1;
                    const bool  name  = hdu.is_string() && !hdu.get<std::string>().empty();
                    if (!index && !name)
                        throw std::invalid_argument(where + ": \"hdu\" is an HDU index (1 = primary) or an EXTNAME, not " +
                                                    hdu.dump() + ".");
                    source["hdu"] = hdu;
                }
            }

            if (entry.contains("exposureTime"))
            {
                const json& t = entry["exposureTime"];
                if (!t.is_number() || !(t.get<double>() > 0.0))
                    throw std::invalid_argument(where + ": \"exposureTime\" must be a number of seconds > 0, not " +
                                                t.dump() + ".");
            }

            std::string output;
            if (entry.contains("output"))
            {
                const json& o = entry["output"];
                if (!o.is_string() || o.get<std::string>().empty())
                    throw std::invalid_argument(where + ": \"output\" must be a file name, not " + o.dump() + ".");
                output = o.get<std::string>();
                if (!outputs.insert(std::filesystem::path(output).lexically_normal().string()).second)
                    throw std::invalid_argument(where + ": \"output\" '" + output + "' is the output of another exposure.");
            }

            json keys = entry;
            keys.erase("file");
            keys.erase("hdu");
            keys.erase("output");
            keys["flux"] = source;

            Exposure exposure{i, PixelConfig(keys), output};
            try                              { _checkExposure(exposure); }
            catch (const std::exception& err) { throw std::invalid_argument(where + ": " + err.what()); }

            exposures.push_back(std::move(exposure));
        }

        m_exposures = std::move(exposures);
    }

    void ScenerySequencer::_checkExposure(const Exposure& e) const
    {
        if (!e.cfg.has("exposureTime"))
            throw std::invalid_argument("no \"exposureTime\" (s).");
    }

    std::string ScenerySequencer::outputPath(const Exposure& e) const
    {
        if (e.output.empty() || m_outputDirectory.empty())
            return e.output;
        return (std::filesystem::path(m_outputDirectory) / e.output).string();   // an absolute output is kept
    }

    void ScenerySequencer::_save(const Exposure& e, DSL::FITScube& image) const
    {
        writeCards(image, {
            {"FLUX",   *e.cfg.card("flux"),    "[ph/s/pix] Input flux (value or file)"},
            {"SEQIDX", e.index + 1,            "Exposure number in the sequence (from 1)"},
            {"SEQNUM", m_exposures.size(),     "Number of exposures in the sequence"}
        });

        const std::string path = outputPath(e);
        if (path.empty())
            return;

        // CFITSIO does not create directories.
        const std::filesystem::path directory = std::filesystem::path(path).parent_path();
        if (!directory.empty())
            std::filesystem::create_directories(directory);

        image.Write(path, true);
    }

    void ScenerySequencer::_checkInitialized() const
    {
        if (m_exposures.empty())
            throw std::logic_error("ScenerySequencer::run: no exposure, call initialize() first.");
    }
}
