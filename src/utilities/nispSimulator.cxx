//
// Created by GILLARD William on 2026-11-19.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
//

#include <nisp/nispSCA.h>
#include <nisp/nispScenerySequencer.h>
#include <phot2pix/phot2pix.h>
#include <DSTCore/OptionManager.h>

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <chrono>
#include <functional>
#include <numeric>
#include <valarray>
#include <variant>
#include <vector>

using namespace DSL;

const std::string DSF_need = "v1.5.0";
const std::string ThisVersion = "v0.1.0";

class InputParameters : public DOption
{
    public:
    
    InputParameters():DOption()
    {
        AddDescription("Generate NISP exposure simulation for single NISP detectors.");
        AddParam("JSON_Parameters","JSON configuration file.");

        
        DOption::prg_name += "nispSimulator -- "+ThisVersion;

        if(DSF::gGIT::version() < DSF_need)
            throw std::runtime_error("[ERROR] UserParameters::UserParameters: Incompatible DSF version, please update your DSF installation to "+DSF_need+" at minimum.");
    }

    virtual ~InputParameters(){};

    bool UserParameters(std::vector<std::string>::const_iterator& it, const std::vector<std::string>& arg)
    {
        return true;
    }
};

// ── Benchmark helpers ─────────────────────────────────────────────────────────
namespace
{
    using Clock = std::chrono::steady_clock;

    const std::string kWarnTag  = "\033[33;40m[WARNING]\033[0m ";
    const std::string kInfoTag  = "\033[32m[INFO]\033[0m ";
    const std::string kDebugTag = "\033[35m[DEBUG]\033[0m ";

    double secondsSince(const Clock::time_point& t0)
    {
        return std::chrono::duration<double>(Clock::now() - t0).count();
    }

    // Print the timing of each step and the statistics of the sca.expose() runs.
    // tWrite < 0 means the output was not written (no -o).
    void printBenchmark(const std::string& label, size_t width, size_t height, size_t nLayers,
                        double tSetup, double tFlux, std::vector<double> tExpose, double tWrite)
    {
        const size_t n    = tExpose.size();
        const double mean = std::accumulate(tExpose.begin(), tExpose.end(), 0.0) / n;
        double var = 0.0;
        for (double t : tExpose) var += (t - mean) * (t - mean);
        const double stddev = (n > 1) ? std::sqrt(var / (n - 1)) : 0.0;

        std::sort(tExpose.begin(), tExpose.end());
        const double median = (n % 2) ? tExpose[n / 2] : 0.5 * (tExpose[n / 2 - 1] + tExpose[n / 2]);

        const double nPix = static_cast<double>(width) * height;

        std::cout << std::fixed << std::setprecision(3)
                  << "\n── nispSimulator benchmark ─────────────────────────────\n"
                  << "  detector     : " << width << " x " << height << "\n"
                  << "  readout      : " << label << "  → " << nLayers << " output layer(s)\n"
                  << "  runs         : " << n << "\n"
                  << "  setup        : " << tSetup << " s   (SCA construction + initialize)\n"
                  << "  flux load    : " << tFlux  << " s\n"
                  << "  sca.expose   : mean " << mean << " s ± " << stddev
                  << "   (min " << tExpose.front() << ", median " << median << ", max " << tExpose.back() << ")\n"
                  << "  per layer    : " << mean / nLayers << " s   ("
                  << std::setprecision(1) << nPix * nLayers / mean / 1e6 << " Mpix/s)\n"
                  << std::setprecision(3);
        if (tWrite < 0.0)
            std::cout << "  write        : skipped (no -o output file)\n";
        else
            std::cout << "  write        : " << tWrite << " s\n";
        std::cout << "────────────────────────────────────────────────────────\n";
    }
}

int main(int argc, char **argv)
{
    try
    {
        InputParameters fOPT = InputParameters();
        if (!fOPT.InitParameters(argc, argv))
            return EXIT_FAILURE;

        const auto tStart = Clock::now();

        // Debug level (-d x):
        //   0  warnings only (+ benchmark report if "benchmark" > 0)
        //   1  [INFO]  run summary: detector, readout, flux source, effects, output, timings
        //   2  [DEBUG] + pixel/SCA parameters, flux statistics, per-step and per-run timings
        //   3+ [DEBUG] + full configuration dump, DSTfits verbose output
        const int16_t dbg = fOPT.DebugLevel();
        switch(dbg)
        {
            case 0:  DSL::verbose = DSL::verboseLevel::VERBOSE_NONE;   break;
            case 1:  DSL::verbose = DSL::verboseLevel::VERBOSE_BASIC;  break;
            case 2:  DSL::verbose = DSL::verboseLevel::VERBOSE_DETAIL; break;
            default: DSL::verbose = DSL::verboseLevel::VERBOSE_ALL;    break;
        }

        const std::string outputFile = fOPT.GetOutput();
        const std::string cfgFile = fOPT.GetFiles()[0];
        phot2pix::PixelConfig cfg(cfgFile);

        // A "flux" array is a sequence of exposures, each with its own timing and "output"
        // (see phot2pix::ScenerySequencer); otherwise a single exposure written to -o.
        const bool sequence = phot2pix::ScenerySequencer::inConfig(cfg);

        // Warn early, before a possibly long simulation.
        if(!sequence && outputFile.empty())
            std::cerr << kWarnTag << "No output file given (-o/--output): the simulated exposure will NOT be written." << std::endl;

        // JSON view of the configuration, used only for printing values.
        const phot2pix::json cfgJson = cfg.toJson();
        auto show = [&cfgJson](const std::string& key) -> std::string
        { return cfgJson.contains(key) ? cfgJson.at(key).dump() : std::string("not set (default)"); };

        if(dbg >= 1)
            std::cout << kInfoTag << "nispSimulator " << ThisVersion << ", configuration '" << cfgFile << "'" << std::endl;
        if(dbg >= 3)
            std::cout << kDebugTag << "full configuration:\n" << cfgJson.dump(4) << std::endl;

        // JSON integers are stored as uint64_t, which is not size_t on macOS:
        // read through getAsDouble() instead of get<size_t>().
        auto getCount = [&cfg](const std::string& key, size_t def) -> size_t
        { return cfg.has(key) ? static_cast<size_t>(cfg.getAsDouble(key)) : def; };

        size_t width  = getCount("width",  2048);
        size_t height = getCount("height", 2048);

        // "benchmark": N → run sca.expose() N times and print timings; absent or 0 → single silent run.
        const size_t nRuns = getCount("benchmark", 0);

        auto tStep = Clock::now();
        nisp::SCA sca(width, height);
        sca.initialize(cfg);
        const double tSetup = secondsSince(tStep);

        // Position of the matrix on the physical detector ("detectorOrigin", (0, 0) by default),
        // parsed by the SCA.
        const auto [x0, y0] = sca.getOrigin();

        if(dbg >= 1)
            std::cout << kInfoTag << "detector       : " << width << " x " << height
                      << " pixels at detectorOrigin (x0, y0) = (" << x0 << ", " << y0 << ")"
                      << (cfg.has("detectorOrigin") ? "" : " (default)") << std::endl;

        if(dbg >= 1)
        {
            bool oneOverF = false;
            if(cfgJson.contains("oneOverF"))
            {
                const auto& o = cfgJson.at("oneOverF");
                oneOverF = o.is_number() ? (o.get<double>() > 0.0)
                                         : (o.is_object() && o.value("amplitude", 0.0) > 0.0);
            }
            std::cout << kInfoTag << "exposure mode  : " << show("exposureMode") << std::endl
                      << kInfoTag << "persistence    : " << (cfg.has("persistence") ? "on" : "off") << std::endl
                      << kInfoTag << "1/f noise      : " << (oneOverF ? "on" : "off") << std::endl;
        }
        if(dbg >= 2)
        {
            for(const char* key : {"gain", "dark", "readNoise", "QE", "KTCnoise", "baseline"})
                std::cout << kDebugTag << std::left << std::setw(15) << key << ": " << show(key) << std::endl;
            for(const char* key : {"frame_time", "fr_time", "offset", "dpu_offset", "scale", "dpu_scale", "oneOverF", "bias_voltage_noise", "pedestal_drift", "persistence",
                                     "detectorSize", "detectorOrigin", "refBorder", "refDark"})
                if(cfg.has(key))
                    std::cout << kDebugTag << std::left << std::setw(15) << key << ": " << show(key) << std::endl;
            std::cout << kDebugTag << "SCA setup      : " << tSetup << " s" << std::endl;
        }

        // ── Sequence of exposures ("flux" array) ─────────────────────────────
        if(sequence)
        {
            if(cfg.has("exposureTime") || nisp::macc::inConfig(cfg))
                std::cerr << kWarnTag << "\"flux\" is an array: the top-level exposureTime / MACC parameters are ignored, "
                                         "each exposure gives its own." << std::endl;
            if(nRuns > 0)
                std::cerr << kWarnTag << "\"benchmark\" is ignored for a sequence of exposures." << std::endl;

            nisp::ScenerySequencer sequencer;
            sequencer.initialize(cfg);   // checks every exposure before the first one is simulated

            // -o is the directory of the outputs: each exposure's "output" is written into it.
            if(!outputFile.empty())
                sequencer.setOutputDirectory(outputFile);

            const auto& exposures = sequencer.getExposures();
            if(std::none_of(exposures.begin(), exposures.end(), [](const auto& e) { return !e.output.empty(); }))
                std::cerr << kWarnTag << "No exposure of the \"flux\" array has an \"output\": nothing will be written." << std::endl;

            if(dbg >= 1)
            {
                std::cout << kInfoTag << "sequence       : " << exposures.size() << " exposure(s)" << std::endl
                          << kInfoTag << "output dir     : " << (outputFile.empty() ? std::string("current directory") : outputFile) << std::endl;
                for(const auto& e : exposures)
                {
                    const nisp::readout r = nisp::readout::fromConfig(e.cfg);
                    const std::string timing = r.maccParams
                        ? "MACC(" + std::to_string(r.maccParams->nGroups) + "," + std::to_string(r.maccParams->nFrames) +
                          "," + std::to_string(r.maccParams->nDrop) + ")"
                        : "exposureTime = " + std::to_string(r.exposureTime) + " s";
                    const std::string mode = e.cfg.has("exposureMode") ? e.cfg.get<std::string>("exposureMode")
                                                                       : show("exposureMode");
                    const phot2pix::CardValue fluxCard = *e.cfg.card("flux");   // "file.fits[hdu]" or a uniform value
                    const std::string fluxText = std::holds_alternative<std::string>(fluxCard)
                        ? std::get<std::string>(fluxCard)
                        : "uniform " + std::to_string(e.cfg.getAsDouble("flux")) + " ph/s/pixel";
                    std::cout << kInfoTag << "  " << e.index + 1 << "/" << exposures.size()
                              << "  flux " << fluxText
                              << ", mode " << mode << ", " << timing
                              << " → " << (e.output.empty() ? std::string("not saved") : sequencer.outputPath(e)) << std::endl;
                }
                std::cout << kInfoTag << "simulating..." << std::endl;
            }

            tStep = Clock::now();
            sequencer.run(sca);

            if(dbg >= 1)
                std::cout << kInfoTag << "sequence       : " << secondsSince(tStep) << " s" << std::endl
                          << kInfoTag << "total          : " << secondsSince(tStart) << " s" << std::endl;
            return EXIT_SUCCESS;
        }

        // Flux (photons/s/pixel), read by the SCA as a per-pixel parameter (see
        // DetectorSimulator::loadImage() for every form):
        //   absent                                    → 0 everywhere (dark exposure)
        //   number                                    → uniform flux
        //   "<path>" or {"file": "<path>", "hdu": h}  → FITS image in HDU h, an index or an EXTNAME
        //                                               (default: 1, the primary HDU): matrix-sized (used
        //                                               as is) or detector-sized (cropped to
        //                                               [x0, x0+width) × [y0, y0+height)), as the per-pixel maps
        tStep = Clock::now();
        std::valarray<double> fluxArray(0.0, width*height);
        std::string fluxSource = "none → dark exposure (flux = 0)";
        std::pair<size_t, size_t> fluxCorner{0, 0};   // corner of the flux image section used
        if(const std::shared_ptr<DSL::FITScube> inputFlux = sca.loadImage(cfg, "flux"))
        {
            const std::string fluxName = std::get<std::string>(*cfg.card("flux"));   // "file.fits[hdu]"
            fluxArray = sca.toMatrix(*inputFlux, "Flux image '" + fluxName + "'");

            const bool matrixSized = inputFlux->Size(1) == width && inputFlux->Size(2) == height;
            if(!matrixSized)
                fluxCorner = {x0, y0};
            fluxSource = "FITS image '" + fluxName + "' (" + std::to_string(inputFlux->Size(1)) + " x " +
                         std::to_string(inputFlux->Size(2)) +
                         (matrixSized ? ", matrix-sized)" : ", detector-sized, cropped at the matrix origin)");
        }
        else if(cfg.has("flux"))
        {
            fluxArray  = sca.loadMap(cfg, "flux");
            fluxSource = "uniform " + show("flux") + " ph/s/pixel";
        }
        const double tFlux = secondsSince(tStep);

        if(dbg >= 1)
            std::cout << kInfoTag << "flux           : " << fluxSource << std::endl;
        if(dbg >= 2)
            std::cout << kDebugTag << "flux stats     : min " << fluxArray.min()
                      << ", mean " << fluxArray.sum() / fluxArray.size()
                      << ", max " << fluxArray.max() << " ph/s/pixel  (load " << tFlux << " s)" << std::endl;

        if(cfg.has("exposureTime") && nisp::macc::inConfig(cfg))
            std::cerr << kWarnTag << "Both 'exposureTime' and MACC parameters are set: using exposureTime, MACC parameters ignored." << std::endl;

        // Pick the exposure call once, so the benchmark loop times only sca.expose().
        std::function<std::shared_ptr<DSL::FITScube>()> exposeFn;
        std::string label;
        if(cfg.has("exposureTime"))
        {
            const double exposureTime = cfg.getAsDouble("exposureTime");
            label    = "exposureTime = " + std::to_string(exposureTime) + " s";
            exposeFn = [&sca, &fluxArray, exposureTime]() { return sca.expose(fluxArray, exposureTime); };
        }
        else if(nisp::macc::inConfig(cfg))
        {
            const nisp::macc maccParams = nisp::macc::fromConfig(cfg);
            label    = "MACC(" + std::to_string(maccParams.nGroups) + "," + std::to_string(maccParams.nFrames) +
                       "," + std::to_string(maccParams.nDrop) + ")";
            exposeFn = [&sca, &fluxArray, maccParams]() { return sca.expose(fluxArray, maccParams); };
        }
        else
            throw std::invalid_argument("No exposure time or MACC parameters provided in configuration.");

        if(dbg >= 1)
            std::cout << kInfoTag << "readout        : " << label << std::endl
                      << kInfoTag << "output         : " << (outputFile.empty() ? "none (not written)" : outputFile) << std::endl
                      << kInfoTag << "simulating..." << std::endl;

        const size_t                   nExpose = std::max<size_t>(1, nRuns);
        std::vector<double>            tExpose;
        std::shared_ptr<DSL::FITScube> result;
        for(size_t i = 0; i < nExpose; ++i)
        {
            result.reset();   // free the previous ramp outside the timed region
            const auto t0 = Clock::now();
            result = exposeFn();
            tExpose.push_back(secondsSince(t0));

            if(dbg >= 2 && nExpose > 1)
                std::cout << kDebugTag << "sca.expose run " << i + 1 << "/" << nExpose << ": " << tExpose.back() << " s" << std::endl;
        }

        // Input flux in the header, next to the detector inputs written by the SCA.
        {
            const phot2pix::CardValue flux = cfg.card("flux").value_or(phot2pix::CardValue(0.0));
            std::vector<phot2pix::HeaderCard> cards{{"FLUX", flux, "[ph/s/pix] Input flux (value or file)"}};
            if(std::holds_alternative<std::string>(flux))
            {
                const auto [fx0, fy0] = fluxCorner;
                cards.push_back({"FLUXSEC",
                                 "[" + std::to_string(fx0 + 1) + ":" + std::to_string(fx0 + width) + "," +
                                       std::to_string(fy0 + 1) + ":" + std::to_string(fy0 + height) + "]",
                                 "Section of the flux image used"});
            }
            phot2pix::writeCards(*result, cards);
        }

        const size_t nLayers = std::max<size_t>(1, result->Size(3));
        if(dbg >= 1)
            std::cout << kInfoTag << "sca.expose     : " << tExpose.back() << " s, "
                      << result->Size(1) << " x " << result->Size(2) << " x " << nLayers << " output" << std::endl;

        double tWrite = -1.0;
        if(!outputFile.empty())
        {
            tStep = Clock::now();
            result->Write(outputFile, true);
            tWrite = secondsSince(tStep);

            if(dbg >= 1)
                std::cout << kInfoTag << "written        : " << outputFile << " (" << tWrite << " s)" << std::endl;
        }

        if(nRuns > 0)
            printBenchmark(label, width, height, nLayers, tSetup, tFlux, tExpose, tWrite);

        if(dbg >= 1)
            std::cout << kInfoTag << "total          : " << secondsSince(tStart) << " s" << std::endl;

    }
    catch (std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        throw;
    }
    catch (...)
    {
        return EXIT_FAILURE;
    }

    EXIT_SUCCESS;
}