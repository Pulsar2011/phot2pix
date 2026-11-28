# Changelog

Notable changes of phot2pix, most recent first. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). No version has been tagged yet: the changes that are already in the repository are listed by date.

## 2026-11-28

### Added

- **Random Telegraph Noise (RTN)** for the nispSCA model.  

## 2026-11-27

### Added

- **Cosmic-ray and radioactivity hits** are simulated by `phot2pix::CosmicRay::CosmicDeposit` (`CosmicDeposit.h`), configured with the `"cosmic"` key. For each step of an exposure, it draws Poisson hits from the spectrum of each species at the detector (CSV or FITS table; energies per nucleon for a nucleus, per particle otherwise), with an angular model (`isotropic`, `surface`, `cosn`, or a tabulated zenith distribution), a sensor tilted by `tilt_deg`, and the faces each species reaches (`front`, `back`, `both`). A particle that stops deposits all its kinetic energy, a thick crossing a Gaussian loss with the Bohr variance, and a thin crossing a Landau loss, from the Bethe-Bloch (heavy particles) and ICRU 37 (e±) stopping powers. The charge follows the CSDA loss along the track, Bragg peak included, and spreads laterally with a core and an optional `halo` per species, and an optional δ-ray halo (`deltaRays`). Tracks longer than `maxTrackPixels` are not drawn (`CRTRUNC`). With a `seed`, the hits are reproducible whatever the threads. The model is described in the new wiki page *CosmicDeposit*.
- `CosmicDeposit` also gives the header cards of its inputs (`CRFLUX`, `CRACC`, …, `CRSEED`), the table of the hits of a step (`generate(…, &events)`, `toColumns()`), rates and probabilities without simulating (`eventRate()`, `expectedEvents()`, `probability()`, `geometricFactor()`, `truncatedFraction()`, `deltaRayShare()`), the sensor materials of `target::material` (silicon by default, `compound()`, `HgCdTe(x)`), and the free functions `incidence()`, `effectiveArea()`, `geometricFactor()` and `etendue()` for the geometry of a flat sensor.
- `nisp::SCA` simulates cosmic rays with `nisp::CosmicDeposit`, whose sensitive layer is Hg₀.₇Cd₀.₃Te. The hits are added to the charge of each step of every readout mode, reference pixels included, so that they accumulate along the ramp; they also leave persistence in the following exposures, and their inputs are recorded in the header.
- `phot2pix::ChargeDeposit` (`ChargeDeposit.h`, included by `phot2pix.h`) is the base of the sources of charge over several pixels: `configure()`, `generate()`, `clone()`, `isConfigured()`, `inputCards()`, and random streams that a seed makes reproducible (`setSeed()`).
- `DetectorSimulator::expose()`, with the same overloads as `simulate()`, returns the pixel response plus the deposits of the detector's charge source, drawn in a thread of their own while the pixels are simulated and divided by the gain of each pixel. `initialize()` creates the source from the `"cosmic"` key with the virtual function `_makeChargeDeposit()`, which returns a silicon `CosmicDeposit` by default, or `nullptr` for no source (`nisp::SCA` and `nisp::persistence::pixelArray`). `getChargeDeposit()` gives read access to the source, and `_addDeposits()` adds a list of deposits to a frame. `simulate()` is unchanged and never adds deposits.
- Reference spectra in `ressource/`: galactic cosmic rays at solar minimum and maximum, solar particle events and an iron spectrum, in free space and behind 10 g/cm² of aluminium (`cosmic_reference_spectra.csv`); alphas and betas of the U-238, Th-232 and K-40 natural radioactivity (`natural_radioactivity_spectra.csv`); and a grazing angular table (`grazing_angular_table.csv`).
- Test configurations with cosmic rays, `nispSimulator_{CDS,DPU,FLUX,MACC}_CR_*.json`, which exercise the angular models, the faces and the spectra. Their spectrum paths are relative to the root of the repository.
- Tests: `gTest_chargeDeposit`, `gTest_cosmicDeposit` and `gTest_detectorExpose`, and new tests of the DPU saturation, of the gain of the persistence signal, and of `ScenerySequencer::run()` with a charge source.
- tag Version number v0.1.0,
- **Inter-pixel capacitance (IPC)** of the NISP detectors, applied at each read in every readout mode. The `ipc` key gives one kernel (array, or `kernel` / `width` / `height`) or one kernel per pixel (FITS cube).
- **Crosstalk** between the output channels of the NISP detectors, applied at each read in every readout mode. The `crosstalk` key gives the couplings by distance between channels; they are recorded in the header (`XTALKCOL`, `XTALKNCH`, `XTALKn`). The ghosts are limited by the full well of their source.
- `colPerChannel` and `nChannels`, top-level configuration keys that set the output channels of the instrument. They are held by `nisp::Geometry` and used by the readout noise, the crosstalk, the reference pixel correction and the DPU. *They are no more part of `oneOverF`.* The default is 32 channels of 64 columns each, as in the NISP H2RG.
- `phot2pix::PixelConvolution` (`DetectorEffect.h`): convolution of an image by one kernel, or by one kernel per pixel. To be included separately from `phot2pix.h`. Its `initialize()` reads the kernel of the configuration key it is given (`ipc` for the IPC of NISP), and its error messages name that key.
- `DetectorTools.h`: two function templates shared by the detectors and their effects. They were private to `nispSCA.cxx` and are unchanged.
  - `phot2pix::parallelFor(n, nThreads, work)` splits the range [0, n) into contiguous chunks and calls `work(begin, end)` on each, one thread per chunk. It returns when all of them are done. With `nThreads` ≤ 1 or fewer than 100 elements per thread, `work(0, n)` is called once, in the calling thread. Used for the charge accumulation, the reads, the reference pixel correction and the DPU of `nisp::SCA`, and now by `PixelConvolution`.
  - `phot2pix::withCubeData(cube, f)` calls `f` with the data of a `DSL::FITScube`, as a `std::valarray` of the type the cube is stored in (8- to 64-bit integers, `float` or `double`), without a copy. `f` must accept any of these types, for example a lambda with an `auto&` parameter. Returns `false`, without calling `f`, for any other type. Used to read and correct ramps in the reference pixel correction, and now to read the kernel files of `PixelConvolution`.
- `phot2pix::parallelFor()` and `withCubeData()` (`DetectorTools.h`), shared by the detectors and their effects.
- `nisp::SCA`: `applyIPC()` and `applyCrosstalk()`, which act on an image in double precision, and accessors to the kernels. Both are applyed in `nisp::SCA::_read` before readnoise and 1/f is applyed to the matrix.
- CMake: when `DSTools/DSTmorpho.h` ([`deepskyutilities`](https://gitlab.in2p3.fr/deepskyproject/deepskyutils)) is not installed, it is downloaded from `deepskyutilities` with the two headers it includes. `DSTools_REF` chooses the branch, tag or commit (default: `main`).
- Tests: `gTest_pixelConvolution`, `gTest_nispIPC` and `gTest_nispCrosstalk`. The environment variable `NISP_IPC_OUTPUT` gives a directory in which the last two test case save their images. If this environement variable is not specified, the generated image aren't saved.
- This changelog.

### Changed

- `ScenerySequencer::run(DetectorSimulator&)` calls `expose()` instead of `simulate()`, so the exposures of a sequence on a plain detector receive the hits of its charge source.
- The DPU of `nisp::SCA` flags a pixel with a saturated read (65535) in any group, as the flight DPU does: its slope plane is 65535, with no scale and no offset, and its χ² and error planes are 0. It used to fit the groups before the first saturated one.
- **CMake, as in the other DeepSky projects:** presets (`default`, `debug`, `asan`, in `CMakePresets.json`). The DeepSky libraries are found with `find_package()` in `DeepSkyROOT` only, at a minimum version: `DSTmath` 0.2.0, `DSTfits` 1.5.1, `DSTCore` 2.0.1 (optional), `DSTools` 2.2.0 (headers only; downloaded when not installed, error when an older one is). The installed `phot2pixConfig.cmake` finds `DSTmath` and `DSTfits` again at these versions. Version from the git tag through `git describe`; no tag, no versioned library and no CMake package. `BUILD_TESTING` is now `BUILD_TESTS`; AddressSanitizer is `WITH_ASAN` (the `asan` preset) instead of every Debug build. Install prefix by default `DeepSkyROOT`, `nispSimulator` installed into `exe/`, build products in `build/bin` and `build/lib`. The `distclean` target is gone.
- **Breaking:** `colPerChannel` and `nChannels` are no longer accepted inside `oneOverF`. Move them to the top level of the configuration. `initialize()` throws if they are still there.
- The DSTools headers (`DSTmorpho.h`, `DSTmorph_ops.h`, `DST_version.h`) are a new dependency of the core. The whole `libDSTools` from [`deepskyutilities`](https://gitlab.in2p3.fr/deepskyproject/deepskyutils) is not required and isn't linked.
- `SCA::getReferencePixelCorrection()` and `SCA::referencePixelsCorrection()` take an optional `nisp::Geometry` (default: NISP). The correction uses its size, reference border and channels, and is only applied to the full detector.
- The test configurations `nispSimulator_CDS`, `nispSimulator_DPU` and `nispSimulator_FLUX` have IPC and crosstalk on. `nispSimulator_FLUX` simulates the full 2048 × 2048 detector instead of a 30 × 50 matrix, with a uniform flux of 10 photons/s/pixel instead of an external FITS image, and `nispSimulator_DPU` uses MACC(4, 16, 4).

### Fixed

- The persistence signal of `nisp::SCA` was added in electrons to frames in ADU. The persistence scheduler now takes the gain of each pixel (`scheduler::setGain()`, `pixelArray::setGain()`), so the released charge is in ADU.
- Out-of-bounds read in the DPU for a matrix larger than 2048 × 2048 with the reference pixel correction on.
- Windowed reads with few rows ran on a single thread: the work is now split by pixel. The full-frame read uses the same code.

## 2026-11-24

### Added

- Photon transfer curve tests: gain and read noise recovered from CDS and DPU exposures.

### Changed

- `nisp::Ramp::getSlope()` returns the squared error of the slope instead of the error. The third plane of the DPU output holds it, × 100 × `scale` + `offset`.

## 2026-11-23

### Added

- Full well and non-linear response of the NISP pixel.
- Doxygen documentation, logo and code of conduct.

## 2026-11-22

### Added

- Sequences of exposures (`ScenerySequencer`).
- Publication of the wiki from the repository.

### Changed

- The code compiles with standard libraries that do not provide `<mdspan>`.

## 2026-11-21

### Added

- Physical detector and placement of the simulated matrix (`phot2pix::Geometry`, `nisp::Geometry`): `detectorSize`, `detectorOrigin`, `refBorder`.
- Reference pixels (`refDark`) and their correction, Eq. 3 of arXiv:2507.11326 (`refPixCorr`).
- Simulation inputs recorded as header cards.

### Changed

- **Breaking:** `x0` and `y0` are replaced by `"detectorOrigin": [x0, y0]`.
- The readout noise is defined on the physical channels.

### Fixed

- `POT2PIX` header keyword, now `PHOT2PIX`.

## 2026-11-20

### Changed

- Revision of the architecture.

### Added

- Wiki.

## 2026-06-15 to 2026-07-03

### Added

- First versions: single pixel, full array, NISP H2RG pixel and detector simulators.
