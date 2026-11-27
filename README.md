![Project logo](Phot2Pix.png)

# phot2pix

**Photon-to-Pixel CCD/CMOS sensor simulator library**

*Centre de Physique des Particules de Marseille (CPPM)*
*Licence: [PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0)*

---

## What is phot2pix?

`phot2pix` is a C++23 library that simulates the conversion of photon flux into detector pixel values. It provides:

- A **generic pixel model** (`pixsim`) implementing shot noise, dark current, read noise, KTC noise, gain, and quantum efficiency.
- A **generic detector** (`DetectorSimulator<PixelT>`) that owns a 2-D array of pixels, dispatches simulation across CPU threads, and writes output to FITS images via the `DSTfits` library.
- A **typed configuration bag** (`PixelConfig`) that can be loaded from JSON files, inline brace-init lists, or FITS-based per-pixel maps.
- An **extension model** based on C++ inheritance and virtual hooks, allowing instrument teams to add detector-specific physics (persistence, crosstalk, charge diffusion, custom readout modes) without modifying the core library.

The library is instrument-agnostic by design. The `src/nisp/` subdirectory provides a worked example for the **Euclid/NISP HAWAII-2RG** infrared detector array.

---

## Dependencies

| Library | Version | Role | Notes |
|---|---|---|---|
| [`DSTmath`]:https://gitlab.in2p3.fr/deepskyproject/deepskymath | v0.2.2 | RNG and math primitives | Required; part of the DeepSky toolkit |
| [`DSTfits`]:https://gitlab.in2p3.fr/deepskyproject/dstfits |  v1.5.1 | FITS I/O (read/write) | Required; wraps `cfitsio` + `wcslib` + `Minuit2` — provided as submodules |
| [`DSTcore`]:https://gitlab.in2p3.fr/deepskyproject/deepskysys | v2.0.0 | Command Line Tools Helpers | Required to build utilities |
| [`DSTutils`](https://gitlab.in2p3.fr/deepskyproject/deepskyutilities) | v2.0.2 | *Optional* | Required for `kernelConvolution`. If not found, CMake fetches only the code needed for the convolution, instead of cloning the whole repository, which itself depends on [`DSTdisplay`](https://gitlab.in2p3.fr/deepskyproject/deepskydisplay). |
| `nlohmann/json` | | JSON parsing for `PixelConfig` | Fetched automatically if not found |
| GoogleTest | | Unit test runner | Git submodule; only needed when `-DBUILD_TESTS=ON` |

C++ standard: **C++23**. CMake ≥ 3.24.

## Doxygen website:
 https://phot2pix-b41428.pages.in2p3.fr

---

## Getting the source

```bash
git clone --recurse-submodules https://gitlab.in2p3.fr/gillard/phot2pix.git
cd phot2pix
```

The `--recurse-submodules` flag fetches the bundled GoogleTest submodule used by the test suite.

---

## Building

`phot2pix` is built like the other DeepSky projects. The DeepSky libraries it needs (available on the [DeepSkyProject]:https://gitlab.in2p3.fr/deepskyproject git repositories) are only looked for in one install tree, `DeepSkyROOT`: the environment variable of that name, or the directory above the project when it is not set. CMake checks their version:

| Library | Repository | Minimum version |
|---|---|---|
| `DSTmath` | deepskymath | 0.2.0 |
| `DSTfits` (and `cfitsio`, `wcslib`, `Minuit2`, installed with it) | deepskyfits | 1.5.0 |
| `DSTCore` | deepskysys | 2.0.1 — optional: without it, `nispSimulator` is not built |
| `DSTools` headers | deepskyutilities | 2.2.0 — when no `DSTools` is installed, its three headers are downloaded instead |

The minimum versions are the `<library>_MIN_VERSION` variables at the top of the *DeepSky libraries* section of `CMakeLists.txt`.

```bash
cmake --preset default                  # configure: Release build in build/
cmake --build build -j                  # build
ctest --preset default                  # run the test suite
cmake --build build --target install    # install into DeepSkyROOT
```

Other presets: `debug` (Debug build in `build-debug/`) and `asan` (Debug build with AddressSanitizer in `build-asan/`). Without the presets, `cmake -B build` gives the same Release build. In-source builds are refused; `cmake -P cmake/CleanInSourceBuild.cmake` removes what a former one left in the project directory.

The version of phot2pix is the git tag `vX.Y.Z` (`git describe`). A clone without a tag builds unversioned libraries and installs no CMake package.

### CMake options

| Option | Default | Description |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release` enables LTO; `Debug` builds with `-O0 -g` |
| `WITH_ASAN` | `OFF` | Build with AddressSanitizer (the `asan` preset) |
| `BUILD_TESTS` | `ON` | Build and register the GoogleTest suite |
| `CMAKE_INSTALL_PREFIX` | `DeepSkyROOT` | Where headers (`include/`), libraries (`lib/`) and `nispSimulator` (`exe/`) are installed |
| `DSTools_REF` | `main` | Branch, tag or commit of deepskyutilities the `DSTools` headers are downloaded from |
| `CI_CIWORKER` | `OFF` | Build for the GitLab CI runner |

### Using phot2pix in your own CMake project

After installation, `find_package(phot2pix)` is available:

```cmake
find_package(phot2pix REQUIRED)
target_link_libraries(myTarget PRIVATE phot2pix::phot2pix)
```

---

## Core library — `src/phot2pix/`

### `PixelConfig` — typed key-value configuration

`PixelConfig` is a JSON-compatible property bag. Values are stored as `std::variant` over all common numeric types plus `std::string` and raw `nlohmann::json` (for arrays and objects).

```cpp
#include <phot2pix/PixelConfig.h>

// Inline construction
phot2pix::PixelConfig cfg = {
    {"dark",      0.001},   // dark current in e-/s
    {"readNoise", 5.0  },   // read noise in e-
    {"gain",      1.5  },   // gain in e-/ADU
    {"QE",        0.85 },
    {"KTCnoise",  20.0 },
    {"baseline",  1000.0}
};

// From a JSON file
phot2pix::PixelConfig cfg2("detector.json");

double rn = cfg.get<double>("readNoise");          // throws if absent
double g  = cfg.get<double>("gain", 1.0);          // with a default
double d  = cfg.getAsDouble("dark");               // type-erased numeric cast
```

Per-pixel maps can be provided as FITS paths inside the config:

```json
{
  "gain":      "gain_map.fits",
  "gain.hdu":  "GAIN",
  "readNoise": 5.0,
  "dark":      "dark_map.fits"
}
```

### `pixsim` — generic pixel model

`phot2pix::pixsim` models a single pixel's response to an incoming photon flux. The default implementation is a linear detector with Poisson shot noise, Gaussian read noise, dark current, KTC noise, gain, and quantum efficiency.

```cpp
#include <phot2pix/PixelSimulator.h>

phot2pix::pixsim px(
    /*darkCurrent*/ 0.001,
    /*readNoise  */ 5.0,
    /*gain       */ 1.5,
    /*QE         */ 0.85,
    /*KTCnoise   */ 20.0,
    /*baseline   */ 1000.0
);

double adu = px.pixelresponse(/*flux [ph/s]*/ 1000.0, /*exposureTime [s]*/ 100.0);
```

The virtual method `pixelresponse()` is the single customisation point for per-pixel physics. A second overload accepting a `PixelConfig` allows runtime parameters (e.g. elapsed time since last exposure) to be passed through the simulation chain without altering the pixel object's state.

### `DetectorSimulator<PixelT>` — generic detector array

`DetectorSimulator<PixelT>` owns a flat row-major array of `PixelT` pixels and orchestrates the simulation across CPU threads. It is parameterised on any class derived from `pixsim`.

```cpp
#include <phot2pix/DetectorSimulator.h>

// A 4096×4096 detector of default pixels, output as 64-bit floating-point FITS
phot2pix::DetectorSimulator<phot2pix::pixsim> det(4096, 4096, DOUBLE_IMG);

// Uniform initialisation from a config
phot2pix::PixelConfig cfg = {{"gain", 1.0}, {"readNoise", 5.0}, {"dark", 0.001},
                              {"QE", 0.9}, {"KTCnoise", 0.0}, {"baseline", 0.0}};
det.initialize(cfg);

// Simulate: uniform flux of 500 ph/s for 100 s → FITS image
std::shared_ptr<DSL::FITScube> frame = det.simulate(500.0, 100.0);

// Per-pixel flux from a FITS image
std::shared_ptr<DSL::FITScube> skyFlux = /* load from file */;
std::shared_ptr<DSL::FITScube> frame2  = det.simulate(skyFlux, 100.0);

// Windowed simulation: re-simulate only a 64×64 region at (col=100, row=200)
std::valarray<double> windowFlux(500.0, 64 * 64);
std::shared_ptr<DSL::FITScube> patch = det.simulate(windowFlux, 100.0, 100, 200, 64, 64);
```

**Thread dispatch** is automatic: `simulate()` partitions the pixel range across `std::thread::hardware_concurrency()` workers when the detector is large enough that parallelism pays off. Each thread calls `_engine()` once — virtual dispatch overhead is O(nThreads), not O(nPixels).  `_engine()` is the second customisation point for detector-level physics which allows for the implementation of crosstalk, charge diffusion, persistence, custom readout modes ...

**Initialisation from config**: `initialize(cfg)` understands the standard keys `gain`, `dark`, `readNoise`, `QE`, `KTCnoise`, `baseline`. Each key can be either a scalar (applied uniformly) or a FITS image holding a per-pixel map: `"dark": "dark.fits"` (primary HDU), `"dark": {"file": "dark.fits", "hdu": h}`, or `"dark": "dark.fits", "dark.hdu": h`, where `h` is an HDU index (1 = primary) or an EXTNAME. `loadMap(cfg, key)` and `loadImage(cfg, key)` read any other key in the same forms (e.g. a flux map).

**Geometry**: the matrix is a part of a physical detector (`phot2pix::Geometry`, `Geometry.h`). `detectorSize` (`[width, height]`, default: the matrix size) and `detectorOrigin` (`[x0, y0]`, default `[0, 0]`) place it: matrix pixel (x, y) is detector pixel (x0 + x, y0 + y). A FITS map is either matrix-sized (used as is) or detector-sized (the section `[x0, x0+W) × [y0, y0+H)` is used), so a full-detector calibration map can drive the simulation of any part of the detector. `getOrigin()` and `toMatrix(image)` give the placement and apply the same rule to any image (e.g. a flux map). A detector family derives from `Geometry` (overriding `configure()`, `isReference()` and `clone()`) to add its own layout, and passes it to the constructor: `DetectorSimulator<PixelT>(W, H, DOUBLE_IMG, std::make_shared<MyGeometry>())`. `nisp::SCA` does this with `nisp::Geometry` (2048 × 2048 HAWAII-2RG, 4 reference pixels on each side, `refBorder`).

### `ScenerySequencer` — sequence of exposures

`ScenerySequencer` reads a sequence of exposures from a `"flux"` array and simulates them in order on one detector:

```json
"flux": [ {"file": "scene1.fits", "hdu": "DET13.SCI", "exposureTime": 100.0, "output": "exp1.fits"},
          {"file": "scene2.fits", "exposureTime": 50.0} ]
```

`file` is the flux image (matrix- or detector-sized, see Geometry), `hdu` an HDU index (1 = primary, the default) or an EXTNAME; `"flux": <ph/s/pixel>` gives a uniform flux instead of `file` (e.g. `0` for a dark); `exposureTime` the exposure time in seconds. An exposure with an `output` is written to that file as soon as it is simulated, with the cards `FLUX`, `SEQIDX` and `SEQNUM`; one without `output` is simulated but not saved (it still counts for persistence). `setOutputDirectory(dir)` puts the relative outputs into `dir` (created if missing); `nispSimulator` uses its `-o` option as that directory. Every entry is checked by `initialize(cfg)` before the first exposure.

```cpp
phot2pix::ScenerySequencer seq;
seq.initialize(cfg);
seq.run(det);            // det: a DetectorSimulator<PixelT>, simulate(flux, exposureTime) per exposure
```

A detector family derives from it to accept its own keys (`_checkExposure()`) and adds a `run()` for its detector. `nisp::ScenerySequencer` accepts `"macc": [ng, nr, nd]` (or `ng`, `nr`, `nd`) instead of `exposureTime`, and an `exposureMode` per exposure; its `run(sca)` calls `SCA::expose()`. Sky generation, variable sources and moving objects (windowed re-simulation) are planned.

---

## Extending phot2pix — building a custom pixel simulator

The design separates two concerns:

- **Per-pixel physics** → subclass `pixsim`, override `pixelresponse()`.
- **Detector-level physics** (crosstalk, charge diffusion, persistence, custom readout) → subclass `DetectorSimulator<YourPixel>`, override `_engine()` and/or `initialize()`.

### Step 1: custom pixel

```cpp
#include <phot2pix/PixelSimulator.h>

namespace myinstrument {

class MyPixel : public phot2pix::pixsim
{
    double m_cteFactor = 1.0;   // charge transfer efficiency

protected:
    void _userValidation() override
    {
        if (m_cteFactor < 0.0 || m_cteFactor > 1.0)
            throw std::invalid_argument("CTE factor must be in [0, 1].");
    }

public:
    MyPixel() : phot2pix::pixsim() {}
    MyPixel(double dark, double rn, double gain, double QE,
            double ktc, double baseline, double cte)
        : phot2pix::pixsim(dark, rn, gain, QE, ktc, baseline), m_cteFactor(cte)
    {
        _userValidation();
    }

    void setCTE(double cte) { m_cteFactor = cte; _userValidation(); }

    // Override the response to incorporate CTE loss
    double pixelresponse(const double& flux, const double& expTime) const override
    {
        double baseSignal = phot2pix::pixsim::pixelresponse(flux, expTime);
        return baseSignal * m_cteFactor;
    }
};

} // namespace myinstrument
```

### Step 2: custom detector

```cpp
#include <phot2pix/DetectorSimulator.h>
#include "MyPixel.h"

namespace myinstrument {

class MyDetector : public phot2pix::DetectorSimulator<MyPixel>
{
public:
    MyDetector(size_t w, size_t h)
        : phot2pix::DetectorSimulator<MyPixel>(w, h, USHORT_IMG) {}

    // Extend initialize() to handle custom config keys
    void initialize(const phot2pix::PixelConfig& cfg) override
    {
        // Apply standard keys first (gain, dark, readNoise, QE, KTCnoise, baseline)
        phot2pix::DetectorSimulator<MyPixel>::initialize(cfg);

        // Apply custom key
        if (cfg.has("CTE"))
        {
            const double cte = cfg.getAsDouble("CTE");
            for (size_t i = 0; i < size(); ++i)
                at(i).setCTE(cte);
        }
    }

protected:
    // Override _engine() to add detector-level effects (e.g. column glow)
    void _engine(const std::valarray<double>& flux,
                 double dt, size_t begin, size_t end,
                 std::valarray<double>& frame) const override
    {
        for (size_t i = begin; i < end; ++i)
            frame[i] = at(i).pixelresponse(flux[i], dt);

        // Add column glow to the last 10 columns
        for (size_t i = begin; i < end; ++i)
        {
            const size_t col = i % m_width;
            if (col >= m_width - 10)
                frame[i] += 2.0; // 2 e- glow
        }
    }
};

} // namespace myinstrument
```

---

## Example: Euclid/NISP (`src/nisp/`)

The `src/nisp/` directory is an in-progress implementation of a pixel simulator for the Euclid **NISP** instrument, which uses **HAWAII-2RG** HgCdTe infrared detector arrays (SCAs).

> **Status:** the persistence model and SCA class are implemented; the full readout mode chain (CDS, MACC, UTR, DPU) is declared but not yet filled in.

### Class hierarchy

```
phot2pix::pixsim
  └── nisp::pixel                    — HAWAII-2RG pixel (currently thin wrapper)

phot2pix::pixsim
  └── nisp::persistence::model       — persistence trap pixel (private inheritance)

phot2pix::DetectorSimulator<nisp::persistence::model>
  └── nisp::persistence::pixelArray  — full-array persistence simulator

phot2pix::DetectorSimulator<nisp::pixel>
  └── nisp::SCA                      — full NISP SCA (2048×2048, uint16 output)
```

### `nisp::pixel`

A thin subclass of `pixsim` reserving room for HAWAII-2RG-specific pixel behaviour:

```cpp
// src/nisp/nispPixel.h
nisp::pixel px(
    /*dark*/     0.0007,
    /*readNoise*/12.0,
    /*gain*/     2.0,
    /*QE*/       0.85,
    /*KTCnoise*/ 0.0,
    /*baseline*/ 0.0
);
```

### `nisp::SCA`

A 2048×2048 detector producing `uint16` FITS output. `initialize()` recognises the standard pixel keys plus NISP-specific ones (`offset`/`dpu_offset`, `scale`/`dpu_scale`, `fr_time`/`frame_time`).

```cpp
#include <nisp/nispSCA.h>

nisp::SCA sca;   // 2048×2048

phot2pix::PixelConfig cfg = {
    {"gain",       2.0   },
    {"readNoise",  12.0  },
    {"dark",       0.0007},
    {"QE",         0.85  },
    {"KTCnoise",   0.0   },
    {"baseline",   0.0   },
    {"dpu_offset", 0.0   },
    {"dpu_scale",  1.0   },
    {"frame_time", 1.45408}
};
sca.initialize(cfg);

std::shared_ptr<DSL::FITScube> frame = sca.simulate(/*flux ph/s*/ 200.0, /*exp s*/ 100.0);
```

### `nisp::persistence` — trap-based persistence model

The persistence sub-namespace implements a multi-trap exponential decay model:

```
signal(dT, expTime) = flux × Σ_k  f_k × e^{−dT/τ_k} × τ_k × (e^{expTime/τ_k} − 1)
```

where `f_k` is the trapped fraction and `τ_k` is the decay time for trap species `k`.

The **`scheduler`** accumulates the charge history of past exposures and computes the total persistence contribution for a new exposure:

```cpp
phot2pix::PixelConfig persistCfg = {
    /* persistence trap parameters as a JSON array */
};
// The JSON value must contain a "persistence" key with trap list:
// { "persistence": [{"f": 0.02, "t": 100.0}, {"f": 0.005, "t": 1000.0}] }

auto sched = std::make_shared<nisp::persistence::scheduler>(persistCfg);

// After each science exposure, register its charge image:
sched->insertCharges(scienceFrame, /*t0=*/ 0.0);

// Before the next exposure, compute the persistence signal:
double gapSeconds    = 30.0;
double nextExpTime   = 100.0;
auto persistSignal   = sched->computePersistenceSignal(gapSeconds, nextExpTime);
// Add persistSignal to the simulated science frame.
```

The `nisp::SCA` can host a `persistence::scheduler` via `insertPersistenceModel()`.

### Initialising persistence from JSON config

```json
{
  "gain":      2.0,
  "readNoise": 12.0,
  "dark":      0.0007,
  "QE":        0.85,
  "KTCnoise":  0.0,
  "baseline":  0.0,
  "persistence": {
    "dT": 0.0,
    "traps": [
      {"f": 0.020, "t":  100.0},
      {"f": 0.005, "t": 1000.0}
    ]
  }
}
```

---

## Running the tests

```bash
cmake --build build -j
ctest --preset default
```

The tests run from the build directory: the files they write (`test_data/`) stay there.

Test sources live in `src/test/`:

| File | Covers |
|---|---|
| `gTest_pixelconfig.cxx` | `PixelConfig` construction, accessors, JSON round-trip |
| `gTest_pixels.cxx` | `pixsim` response, setters, copy semantics |
| `gTest_detector.cxx` | `DetectorSimulator` geometry, initialisation, simulate |
| `gTest_persistence.cxx` | `nisp::persistence` model and scheduler |
| `gTest_scenerySequencer.cxx` | `ScenerySequencer`: the `flux` array, run and outputs |
| `gTest_nispScenerySequencer.cxx` | `nisp::ScenerySequencer`: MACC / exposureMode per exposure |

---

## Project structure

```
phot2pix/
├── CMakeLists.txt
├── CMakePresets.json              ← default, debug and asan builds
├── src/
│   ├── phot2pix/                  ← generic library (instrument-agnostic)
│   │   ├── PixelConfig.h          ← typed key-value config bag
│   │   ├── PixelSimulator.h/.cxx  ← pixsim base class
│   │   ├── Geometry.h             ← physical detector, matrix placement, map cropping
│   │   ├── DetectorSimulator.h    ← DetectorSimulator<PixelT> template
│   │   └── ScenerySequencer.h/.cxx ← sequence of exposures ("flux" array)
│   ├── nisp/                      ← Euclid/NISP instrument example
│   │   ├── nispPixel.h/.cxx       ← nisp::pixel
│   │   ├── nispGeometry.h         ← nisp::Geometry (H2RG, reference pixels)
│   │   ├── nispSCA.h/.cxx         ← nisp::SCA + persistence model
│   │   └── nispScenerySequencer.h/.cxx ← nisp::ScenerySequencer (MACC, exposureMode per exposure)
│   └── test/                      ← GoogleTest suite
├── googletest/                    ← git submodule
└── cmake/
    ├── phot2pixConfig.cmake.in    ← CMake package (find_package(phot2pix))
    ├── GetGitRevisionDescription.cmake(.in) ← version from the git tag
    └── CleanInSourceBuild.cmake   ← removes what a former in-source build left
```

---

## Authors

William Gillard — CPPM, Aix-Marseille Université / CNRS

*The Phot2Pix logo has been generated by OpenAI's ChatGPT Image Generation tool (ImageGen).*
