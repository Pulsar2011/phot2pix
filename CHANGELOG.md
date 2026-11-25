# Changelog

Notable changes of phot2pix, most recent first. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). No version has been tagged yet: the changes that are already in the repository are listed by date.

## [Unreleased]

### Added

- tag Version number v0.1.0,
- **Inter-pixel capacitance (IPC)** of the NISP detectors, applied at each read in every readout mode. The `ipc` key gives one kernel (array, or `kernel` / `width` / `height`) or one kernel per pixel (FITS cube).
- **Crosstalk** between the output channels of the NISP detectors, applied at each read in every readout mode. The `crosstalk` key gives the couplings by distance between channels; they are recorded in the header (`XTALKCOL`, `XTALKNCH`, `XTALKn`). The ghosts are limited by the full well of their source.
- `colPerChannel` and `nChannels`, top-level configuration keys that set the output channels of the instrument. They are held by `nisp::Geometry` and used by the readout noise, the crosstalk, the reference pixel correction and the DPU. *They are no more part of `oneOverF`.* The default is 32 channels of 64 columns each, as in the NISP H2RG.
- `phot2pix::PixelConvolution` (`DetectorEffect.h`): convolution of an image by one kernel, or by one kernel per pixel. To be included separately from `phot2pix.h`.
- `DetectorTools.h`: two function templates shared by the detectors and their effects. They were private to `nispSCA.cxx` and are unchanged.
  - `phot2pix::parallelFor(n, nThreads, work)` splits the range [0, n) into contiguous chunks and calls `work(begin, end)` on each, one thread per chunk. It returns when all of them are done. With `nThreads` ≤ 1 or fewer than 100 elements per thread, `work(0, n)` is called once, in the calling thread. Used for the charge accumulation, the reads, the reference pixel correction and the DPU of `nisp::SCA`, and now by `PixelConvolution`.
  - `phot2pix::withCubeData(cube, f)` calls `f` with the data of a `DSL::FITScube`, as a `std::valarray` of the type the cube is stored in (8- to 64-bit integers, `float` or `double`), without a copy. `f` must accept any of these types, for example a lambda with an `auto&` parameter. Returns `false`, without calling `f`, for any other type. Used to read and correct ramps in the reference pixel correction, and now to read the kernel files of `PixelConvolution`.
- `phot2pix::parallelFor()` and `withCubeData()` (`DetectorTools.h`), shared by the detectors and their effects.
- `nisp::SCA`: `applyIPC()` and `applyCrosstalk()`, which act on an image in double precision, and accessors to the kernels. Both are applyed in `nisp::SCA::_read` before readnoise and 1/f is applyed to the matrix.
- CMake: when `DSTools/DSTmorpho.h` ([`deepskyutilities`](https://gitlab.in2p3.fr/deepskyproject/deepskyutils)) is not installed, it is downloaded from `deepskyutilities` with the two headers it includes. `DSTools_REF` chooses the branch, tag or commit (default: `main`).
- Tests: `gTest_pixelConvolution`, `gTest_nispIPC` and `gTest_nispCrosstalk`. The environment variable `NISP_IPC_OUTPUT` gives a directory in which the last two test case save their images. If this environement variable is not specified, the generated image aren't saved.
- This changelog.

### Changed

- **Breaking:** `colPerChannel` and `nChannels` are no longer accepted inside `oneOverF`. Move them to the top level of the configuration. `initialize()` throws if they are still there.
- The DSTools headers (`DSTmorpho.h`, `DSTmorph_ops.h`, `DST_version.h`) are a new dependency of the core. The whole `libDSTools` from [`deepskyutilities`](https://gitlab.in2p3.fr/deepskyproject/deepskyutils) is not required and isn't linked.
- `SCA::getReferencePixelCorrection()` and `SCA::referencePixelsCorrection()` take an optional `nisp::Geometry` (default: NISP). The correction uses its size, reference border and channels, and is only applied to the full detector.
- The test configurations `nispSimulator_CDS`, `nispSimulator_DPU` and `nispSimulator_FLUX` have IPC and crosstalk on. `nispSimulator_FLUX` simulates the full 2048 × 2048 detector instead of a 30 × 50 matrix.
- Wiki and Doxygen documentation: IPC, crosstalk, output channels, and where an effect that couples pixels belongs.

### Fixed

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
