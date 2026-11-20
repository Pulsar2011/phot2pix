# NISP example: the Euclid HAWAII-2RG detectors

`src/nisp/` is a complete instrument built on the phot2pix core: a simulator of the infrared detectors of Euclid's **NISP** (Near-Infrared Spectrometer and Photometer). It shows how an instrument is assembled from the generic building blocks described in [Philosophy](Philosophy).

---

## 1. The instrument in a few lines

- NISP has 16 **Teledyne HAWAII-2RG** (H2RG) HgCdTe arrays. Each SCA (sensor chip assembly) is 2048 × 2048 pixels of 18 µm, including a 4-pixel border of reference pixels (2040 × 2040 sensitive).
- An SCA is read **non-destructively**, one full frame every **1.45408 s**, through **32 output channels** of 64 columns each.
- Reads are taken in **MACC(ng, nr, nd)** mode: *ng* groups of *nr* consecutive reads, separated by *nd* frames that are read but dropped.
  - Photometry uses MACC(4,16,4) (87.2 s between the first and last group).
  - Spectroscopy uses MACC(15,16,11) (549.6 s).
- The onboard **DPU** averages each group and fits a slope to the group averages (Kubik et al. 2016, PASP 128, 104504). Only that slope image is sent to the ground.

---

## 2. Map of the code

| File | Content |
|---|---|
| `nispPixel.h/.cxx` | `nisp::pixel`: the H2RG pixel. `nisp::Ramp`: the DPU slope estimator (group co-adding, Kubik slope, Fowler, CDS) |
| `nispSCA.h/.cxx` | `nisp::SCA`: the detector and its readout modes. `nisp::macc`: MACC parameters. `nisp::persistence`: the persistence model (`model`, `pixelArray`, `scheduler`) |
| `nispOneOverFNoise.h/.cxx` | `nisp::OneOverFNoise`: correlated 1/f noise of the output channels |
| `../utilities/nispSimulator.cxx` | Command-line tool driven by a JSON file |
| `../utilities/*.json` | A default configuration with averaged NISP values, and one test configuration per readout mode |

```mermaid
classDiagram
    pixsim <|-- pixel : nisp
    pixsim <|-- model : nisp persistence
    DetectorSimulator~pixel~ <|-- SCA : private
    DetectorSimulator~model~ <|-- pixelArray
    SCA o-- scheduler : persistence
    scheduler o-- pixelArray
    SCA o-- OneOverFNoise : 1/f noise
    SCA ..> Ramp : DPU slope
```

How the layers map onto the [philosophy](Philosophy):

- **Pixel:** `nisp::pixel` uses the three `pixsim` stages unchanged (reset level, integration, read). NISP-specific pixel physics would go there, by overriding a stage.
- **Detector:** `nisp::SCA` inherits `DetectorSimulator<nisp::pixel>` **privately**. Its simulation engines produce **charge only** (`resetLevel() + integrate()`), never a complete read.
- **Instrument:** `nisp::SCA` sequences the reads of each readout mode, adds the read noise (through `pixel::read()`) and the 1/f noise at each read, and runs the DPU processing. Its only public simulation entry point is `expose()`.

---

## 3. Configuration

A `nisp::SCA` is configured with one `PixelConfig`, typically a JSON file.

| Key | Meaning | Default |
|---|---|---|
| `gain`, `dark`, `readNoise`, `QE`, `KTCnoise`, `baseline` | Pixel parameters (units in [Philosophy §3](Philosophy#3-units)). Each can be a number or a FITS map, with optional `"<key>.hdu"` | `pixsim` defaults |
| `exposureMode` | `MACC`, `UTR`, `CDS`, `FLUX` or `DPU` (case-insensitive); anything else throws | `MACC` |
| `macc` | `[ng, nr, nd]`; or the three keys `ng`, `nr`, `nd`; not both | — |
| `frame_time` (or `fr_time`) | Frame time, s | 1.45408 |
| `offset` (or `dpu_offset`) | ADU added before conversion to `uint16` | 0 |
| `scale` (or `dpu_scale`) | Factor applied before conversion to `uint16` | 1 |
| `dpu_gain`, `dpu_noise` | Gain (e⁻/ADU) and read noise (e⁻) assumed by the DPU slope estimator, written as `GAIN` and `RN` in the header | Mean pixel gain and read noise |
| `oneOverF` | 1/f noise: a number (amplitude), or `{"amplitude", "slope", "nChannels", "colPerChannel", "nPoles"}` | Off |
| `persistence` | Trap list `[{"f": …, "t": …}, …]`, or `{"traps": [...], "dT": seconds}` | Off |

`dpu_gain` and `dpu_noise` exist to simulate a DPU whose calibration differs from the real detector. Without them, the DPU uses the detector's own mean values.

Example, averaged NISP values (Euclid Collaboration 2025, arXiv:2507.11326), from `src/utilities/nispSimulator_default.json`:

```json
{
    "width": 2048, "height": 2048,
    "flux": 0.0,
    "exposureMode": "MACC",
    "macc": [4, 16, 4],
    "frame_time": 1.45408,
    "gain": 1.92, "dark": 0.01, "readNoise": 13.0,
    "QE": 0.94, "KTCnoise": 48.0, "baseline": 10000.0
}
```

> When the configuration contains `persistence`, keep the `PixelConfig` object alive as long as the `SCA`: the persistence scheduler keeps a reference to it.

---

## 4. Readout modes

### Entry points

```cpp
nisp::SCA sca(2048, 2048);
sca.initialize(phot2pix::PixelConfig("config.json"));
std::valarray<double> flux(10.0, sca.size());            // photons / s / pixel

auto out  = sca.expose(flux, nisp::macc{4, 16, 4});      // readout pattern given explicitly
auto out2 = sca.expose(flux, 100.0);                     // exposure time in seconds
auto win  = sca.expose(flux, 100.0, c0, r0, nc, nr);     // single read of a window
```

The mode set by `exposureMode` decides what `expose()` computes:

| Mode | `expose(flux, macc)` | `expose(flux, t)` | Output (`uint16`) |
|---|---|---|---|
| **MACC** | Full ramp, groups and drops | UTR ramp of ⌊t / t_frame⌋ reads | W × H × (ng·nr) cube of reads |
| **UTR** | Same as MACC | UTR ramp of ⌊t / t_frame⌋ reads | W × H × N cube of reads |
| **CDS** | Last read − first read over the MACC duration | Last − first over *t* | W × H |
| **FLUX** | One read after the MACC duration | One read after *t* | W × H |
| **DPU** | Ramp, then onboard slope fitting | Slope fitted on a UTR ramp | W × H slope map |

The MACC duration is (ng·(nr+nd) − nd)·t_frame, i.e. 76 frames for MACC(4,16,4).

### How a ramp is built (MACC, UTR, and the input of DPU)

```
frame : 1 ······ 16 │ 17 ·· 20 │ 21 ······ 36 │ 37 ·· 40 │ ··· │ 61 ······ 76
        group 0     │  drops   │ group 1      │  drops   │     │ group 3
read  : R  R ··· R  │  · · · · │ R  R ··· R   │  · · · · │     │ R  R ··· R     → 64 layers
```

1. One array `cumul` holds the **charge** of every pixel, in ADU.
2. For each frame that is read:
   - `cumul += resetLevel() + integrate()` on the first frame, `cumul += integrate()` afterwards;
   - then the layer is written as `uint16( pixel.read(cumul) + 1/f )`.
3. For dropped frames: `cumul += integrate()` over nd frames; nothing is read or stored.
4. Read noise and 1/f noise are added to the layer **only**, never to `cumul`. Every read has RN², and two consecutive reads differ by 2·RN², as on the real detector.

The output cube is allocated once at its final size. Each ramp step (accumulate, read, convert) is a single threaded pass over memory.

### CDS

Correlated double sampling, computed literally:

```
f0 = read( charge at read 1 )            charge = reset level + 1 frame of signal
f1 = read( charge at read N )            charge = charge at read 1 + (N−1) frames
CDS = f1 − f0
```

The reset level and kTC noise cancel. The read noise and 1/f noise of the **two** reads remain (variance 2·RN²). The result goes through `offset`/`scale`, so negative values (noise on dark pixels) survive when an offset is set.

### FLUX

A single read after the whole exposure: reset level + integration over *t* + one read. This is the fastest mode and is useful for quick image simulations.

### DPU: the onboard slope

The MACC ramp is simulated first, then for each pixel (threaded, with no allocation per pixel):

1. **Saturation:** the reads are copied group by group up to the first group containing a saturated read (65535). Only the groups before it are used. With fewer than 2 good groups, the pixel is set to 65535.
2. **Co-adding:** each group is averaged (`Ramp::coadd`).
3. **Slope:** the Kubik et al. 2016 estimator (`Ramp::getSlope`) is applied, with σ = `dpu_noise` / `dpu_gain` (in ADU).
4. **Output:** `uint16(slope × scale + offset)`, in **ADU per group interval** ((nr+nd) frames).

To convert to a flux: flux [e⁻/s] = slope × gain / ((nr+nd)·t_frame).

### Header keywords

| Keyword | Content |
|---|---|
| `MODE` | Readout mode |
| `NG`, `NR`, `ND` | MACC parameters (MACC, DPU) |
| `INTTIME`, `EXPTIME` | Integration and exposure times, s |
| `FRTIME` | Frame time, s |
| `GAIN`, `RN` | DPU gain (e⁻/ADU) and read noise (e⁻) |
| `OFFSET`, `SFACTOR` | DPU offset and scale |
| `PERSIST` | `ON` when persistence is modelled |
| `PHOT2PIX` / `POT2PIX` | Library version |

---

## 5. Noise and detector effects

| Effect | Where | Model |
|---|---|---|
| Shot noise (signal, dark) | `pixsim::integrate()` | Gaussian approximation of Poisson; accumulates |
| Reset (kTC) noise, baseline | `pixsim::resetLevel()` | Drawn once per ramp, at the first frame |
| Read noise | `pixsim::read()` | Independent at each read |
| 1/f noise | `nisp::OneOverFNoise` | Per output channel, independent at each read, continuous in time over a ramp |
| Persistence | `nisp::persistence` | Charge released by traps filled during previous exposures |
| Digitisation | `nisp::SCA` | Rounding and clamping to `uint16`, after `offset`/`scale` |

### 1/f noise

The H2RG is read through parallel output amplifiers, each covering a vertical stripe of `colPerChannel` columns (64 by default, so 32 channels on 2048 columns). Within a stripe, pixels are read row by row, so the amplifier's low-frequency noise appears as **horizontal banding**.

The noise is generated per channel as a sum of `nPoles` AR(1) processes with log-spaced corner frequencies. Together they approximate a power spectrum A²/f^slope without any FFT. Channels are generated in parallel.

Within a ramp, a single continuous stream is used: each read, and each dropped frame, consumes the next chunk of the same process, so the amplifier noise is correlated in time as on the real detector. A 2040-wide science frame (reference pixels removed) is recognised automatically.

### Persistence

Charge trapped during an exposure is released during the following ones. For each trap species *k* (trapped fraction `f`, decay time `t` = τ_k), the signal released into an exposure of duration *T* that starts *dT* after a past exposure is

```
signal = charge × Σ_k  f_k · e^(−dT/τ_k) · τ_k · (e^(T/τ_k) − 1)
```

The `persistence::scheduler` stores the charge map of every exposure. During each later frame it computes the released signal from all of them, on half of the available threads, while the pixel charge is simulated on the other half, and adds it to the charge. The optional `"dT"` sets an initial delay.

---

## 6. The `nispSimulator` tool

`src/utilities/nispSimulator.cxx` simulates one SCA from a JSON configuration (built when the DSTCore library is available; installed in `bin/`).

```bash
nispSimulator -o exposure.fits -d 1 my_config.json
```

| Option / key | Effect |
|---|---|
| `-o <file>` | Output FITS file. Without it, the simulation runs but nothing is written (a warning is printed) |
| `-d <level>` | 0: warnings only; 1: run summary; 2: parameters, flux statistics, timings; 3+: full configuration dump |
| `"flux"` | Absent: dark exposure (flux 0). A number: uniform flux in ph/s/pixel. A string: FITS image, cropped to `[x0, x0+width) × [y0, y0+height)` |
| `"width"`, `"height"`, `"x0"`, `"y0"` | Detector size and position of the window in the flux image |
| `"exposureTime"` or `"macc"` | Readout duration or pattern; `exposureTime` wins if both are given (with a warning) |
| `"benchmark": N` | Runs the exposure N times and prints mean ± σ, min/median/max, time per layer and Mpix/s |

Ready-made configurations:

- `src/utilities/nispSimulator_default.json`: averaged NISP values.
- `src/utilities/testcases/`: one configuration per readout mode, with every detector effect off and 10 ph/s.

---

## 7. Validation

`src/test/gTest_nispReadoutModes.cxx` checks every readout mode on a full 2048 × 2048 detector, with MACC(4,16,4), 10 ph/s, gain 1, QE 1, and every other effect off. Each frame then adds S = 10 × 1.45408 = 14.5408 ADU, with Poisson variance S. Means and variances are checked over the 4.2 million pixels within 5σ.

| Mode | Expected | Flux recovered from |
|---|---|---|
| MACC | Layer 16g+f = (20g+f+1)·S | (G₃ − G₀) / (3·20·t_frame) |
| UTR | Read k = (k+1)·S, 76 reads | Least-squares slope |
| CDS | 75·S = 1090.56 ADU | CDS / (75·t_frame) |
| FLUX | 76·S = 1105.10 ADU | Value / (76·t_frame) |
| DPU | ≈ 290.8 ADU per group, scatter ≈ 9.4 ADU | Value / (20·t_frame) |

`src/test/gTest_nispSCA.cxx` checks the noise budget:

- Every MACC read carries RN² and read-to-read differences carry 2·RN² (read noise does not accumulate).
- A CDS carries 2·RN² while the reset level and kTC noise cancel.
- 1/f noise produces column-correlated banding.
- Persistence raises the signal of the next exposure.

`src/test/gTest_nispPixel.cxx` checks the pixel stages and the `Ramp` estimator.

---

## 8. Known limitations

- With MACC parameters, the UTR mode produces the same output as MACC, and the MACC output always writes `MODE = MACC`.
- The gain convention inside `Ramp::getSlope()` has only been validated at gain 1.
- The 1/f amplitude is documented in electrons but added in ADU.
- The charge stored for persistence after a MACC ramp uses the last integration interval instead of the total exposure.
- Saturation is applied at `uint16` conversion only; there is no full-well model on the accumulated charge yet.
- The version keyword is `PHOT2PIX` in MACC output and `POT2PIX` elsewhere.

<!-- This page is generated from wiki/NISP-example.md in the phot2pix repository. Edit it there: changes made in the web interface are overwritten on the next push. -->
