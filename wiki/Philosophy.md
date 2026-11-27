# Philosophy of phot2pix

This page explains how phot2pix is organised and **why**. It is not an API reference; the headers in `src/phot2pix/` hold the details and the API also accessible [here](https://phot2pix-b41428.pages.in2p3.fr) .

---

## 1. Phot2pix : a layered detector model

A detector simulation mixes very different effects: the physics of one pixel, the geometry of an array, the readout sequence of a given instrument, its configuration, and the observing sequence. phot2pix gives each of them its own layer, so that one can be changed without touching the others.

| Layer | Class | Responsible for | Knows nothing about |
|---|---|---|---|
| **Configuration** | `PixelConfig` | The model parameters, read from a JSON files. It holds the physical parameters to configure a simulation run | Physics |
| **Pixel** | `pixsim` and subclasses | The physics of **one** single pixel: reset, charge collection, read. This class holds the physical response of the pixel to illumination. | Neighbouring pixels, time sequences, readout patterns |
| **Detector** | `DetectorSimulator<PixelT>` | The pixel array, per-pixel maps, multi-threading, FITS output, effects that couple pixels | Readout patterns |
| **Instrument** | e.g. `nisp::SCA` | The readout modes, electronics noise, onboard processing, output headers | Scenery |
| **Scene** | `ScenerySequencer` | The sequence of exposures (`"flux"` array, see [Sequences of exposures](ScenerySequencer)); later the sky, variable and moving sources | — |

The core (`src/phot2pix/`) contains the first three layers and is **instrument-agnostic**. An instrument lives in its own directory, as `src/nisp/` does, and builds on the core by inheritance.

<span style="color:#FF0000" bf>**At the current time,** `ScenerySequencer` **is a placeholder for the scene layer. It is not yet implemented, but it will be in the future.**</span>

---

## 2. The pixel: reset, integrate, read

### The physics

Current model have been initially desinged for CMOS sensor array. `phot2pix::pixsim` implements the three stages of a CMOS pixel:

1. **Reset:** the pixel is set to a reference level. That level is not exactly the same at every reset ($kTC$ or reset noise), and it stays fixed until the next reset.
2. **Integration:** photons and dark current deposit charge, which **accumulates**. Its fluctuation is Poisson (shot) noise, so it accumulates too.
3. **Read:** the electronics measure the charge. Each measurement adds its own electronic noise (read noise), and **the measurement does not change the charge**. Reading the same pixel twice gives two independent noise draws on the same charge.

### The model

`phot2pix::pixsim` implements exactly these three stages, one virtual method each:

| Stage | Method | Returns (ADU) | Default physics |
|---|---|---|---|
| Reset | `resetLevel()` | Reset level | baseline + N(0, kTC)/gain |
| Integration | `integrate(flux, t)` | Charge collected during *t* | N(S, √S)/gain, with S = (flux·QE + dark)·t electrons |
| Read | `read(charge)` | Measured value of a given charge | charge + N(0, RN)/gain |

A **single-read exposure** chains the three: `pixelresponse(flux, t) = read(resetLevel() + integrate(flux, t))`. For an input signal $S$ on a pixel of baseline $B$, the mean count in the pixel is $S/g + B$ and its variance $(S + RN^2 + kTC^2)/g^2$.

A **multi-read sequence** (a ramp or UTR) calls the stages directly:

```
read 1 :  charge  = resetLevel() + integrate(flux, t_frame)    value_1 = read(charge)
read k :  charge += integrate(flux, t_frame)                   value_k = read(charge)
```

`phot2pix::pixsim` class can easily be used for CCD pixel too by tuning the parameters for a CCD, and by **overriding the** `read()` **method** to implement the CCD readout model (which is different from CMOS).

### Why pixels are simulated on their own

phot2pix starts from the assumption that **pixels are independent**: the response of a pixel depends only on the flux it receives, on its own parameters and on the time it integrates. It does not depend on its neighbours or on what happened before. Everything that breaks this assumption is added *on top* of the pixel simulation, by the layer that old the physical model to correlate pixels.

**Array effects** couple pixels together. Examples are inter-pixel capacitance, crosstalk, the 1/f noise shared by all the pixels of one output channel, and column or amplifier glow. These effects need to know the geometry of the array and the way it is read out, so they belong to the detector or instrument layer, not to the pixel. The calculation itself is be generic: `phot2pix::PixelConvolution` function is used to convolves charge deposited in the pixel area with a dedicated kernel. The instrument model is then used to decide both the kernel and at which point of the readout the convolution is applyed. Default `phot2pix::DetectorSimulator` class **do not apply any convolution**.

**History effects** depend on time. Examples are the charge accumulated along a ramp, the persistence released by traps filled during previous exposures, and the readout sequence itself. These effects need a memory of the past, which is added on top of the pixel response. In the NISP example, the readout loop keeps the accumulated charge and the persistence scheduler keeps the past exposures.

This separation is what makes the simulator understandable:

- **Each layer can be validated on its own.** The pixel is checked against its analytic mean and variance. Each array or history effect is then checked as a difference from the independent-pixel result.
- **Each effect lives in one place**.
- **Independent pixels can be simulated in parallel** with no synchronisation, and that is where most of the computing time goes. This is verry usefull in testing the same pixel model is used in a multi-threaded ramp simulation millions of pixels.

The three pixel stages follow from the same reasoning applied to time. A ramp is a history: the charge accumulates from one read to the next, but the noise of a read does not. Because the read is separate from the charge, the readout loop can accumulate `integrate()` and call `read()` at each read. Read noise then stays at the readout level instead of growing along the ramp.

### Statistics

- Shot noise is drawn from a **Gaussian approximation of Poisson** (mean S, variance S), which is accurate for a few tens of electrons and above.
- Random numbers come from a `thread_local` generator (`DST::Math`), so simulations can run on many threads without locks.

---

## 3. Units

| Quantity | Unit | Config key |
|---|---|---|
| Incident flux | photons / s / pixel | (input array) |
| Quantum efficiency | 0 -- 1 | `QE` |
| Dark current | e⁻ / s | `dark` |
| Read noise | e⁻ (RMS, per read) | `readNoise` |
| Reset (kTC) noise | e⁻ (RMS, per reset) | `KTCnoise` |
| Gain | **e⁻ / ADU** | `gain` |
| Baseline | ADU | `baseline` |
| Stage outputs | ADU (double precision) | — |

Every stage returns ADU in double precision. Conversion to the final integer type (e.g. `uint16`) happens once, at output, in the detector or instrument layer.

**Inheritance can add their own parameters**. For example, the NISP SCA adds `fullWell` and `nonLinearity` to the pixel, and `exposureMode`, `frame_time` and `oneOverF` to the instrument.

---

## 4. Configuration-driven

All parameters enter through a `PixelConfig`, a typed key/value bag that can be built:

- **from a JSON file:** `PixelConfig cfg("detector.json");`
- **inline:** `PixelConfig cfg{{"gain", 1.9}, {"readNoise", 13.0}};`
- **with per-pixel maps:** any standard key can be a FITS file path instead of a number, with an optional `"<key>.hdu"` to select the HDU. For example it is possible to provide readnoise map as `"readNoise": <file>` coupled with `"readNoise.hdu": <idx>||<EXTNAME>` or `"readNoise": {"file":<file>, "hdu":<idx>||<EXTNAME>}`.

`DetectorSimulator::initialize(cfg)` understands the standard pixel keys (`gain`, `dark`, `readNoise`, `QE`, `KTCnoise`, `baseline`). Derived classes add their own keys by overriding `initialize()` and calling the base version first.

Two rules:

- **Configuration is parsed in the library, not in the tools.** A command-line tool, a test and a Python binding using the same config must get the same detector.
- **Read numbers with `getAsDouble()`.** JSON integers are stored as `uint64_t`, which is not the same type as `size_t` on every platform and `size_t` is integer size element of the array and time series; `getAsDouble()` accepts any numeric type and is converted afterward as `size_t` when needed.

---

## 5. Extending phot2pix: pick the right hook

| You want to… | Override | Example |
|---|---|---|
| Change **one** physical effect of a pixel | The stage it belongs to: `integrate()`, `resetLevel()` or `read()` | Cosmic-ray hits or charge trapping → `integrate()`; non-Gaussian read noise, full well and non-linearity → `read()` |
| Replace the whole single-read model | `pixelresponse()` (virtual) | A completely different detector technology |
| Add an effect that depends on the position of the pixel | `DetectorSimulator::_engine()` | Column or amplifier glow |
| Add an effect that couples pixels | No override: a member of the detector (`PixelConvolution`, in `DetectorEffect.h`), applied by the readout to the whole image. | Inter-pixel capacitance and crosstalk, in `nisp::SCA` |
| Add configuration keys | `initialize(cfg)`, calling the base first | NISP: `exposureMode`, `frame_time`, `oneOverF`… |
| Record the inputs in the output header | `setInputCard()` in `initialize()` (`PixelConfig::card(key)` gives a key's value or map file name), then `writeInputCards()` on the output. The base already records the six pixel parameters (`PIXGAIN`, `PIXDARK`, …) | NISP: readout noise, persistence, reference pixels, full well and non-linearity |
| Add a readout pattern | An instrument class that sequences the stages | NISP: MACC, UTR, CDS, DPU |
| Simulate a sequence of exposures on your instrument | A `ScenerySequencer` subclass: `_checkExposure()` for the keys of one exposure, and a `run()` that calls your detector's exposure function ([details](ScenerySequencer#4-building-your-own-sequencer)) | NISP: `macc` and `exposureMode` per exposure, `run(SCA&)` |

### Example: a pixel hit by cosmic rays

Cosmic-ray hits deposit charge during integration, so they are a charge effect and belong in `integrate()`. The reset level and the read noise are inherited unchanged.

```cpp
#include <phot2pix/PixelSimulator.h>
#include <random>

class CosmicPixel : public phot2pix::pixsim
{
    double m_rate   = 1.0e-4;   // hits / s / pixel
    double m_charge = 2000.0;   // e- per hit

public:
    using phot2pix::pixsim::pixsim;      // inherit the constructors
    using phot2pix::pixsim::integrate;   // keep the cfg overload visible (C++ name hiding)

    double integrate(const double& flux, const double& t) const override
    {
        thread_local std::mt19937_64 rng{std::random_device{}()};
        std::poisson_distribution<int> hits(m_rate * t);
        return phot2pix::pixsim::integrate(flux, t) + hits(rng) * m_charge / getGain();   // ADU
    }
};
```

`CosmicPixel::pixelresponse()` now gives a single read with cosmic rays, and any ramp simulator that uses the stages gets cosmic rays that accumulate correctly from read to read.

Things to keep in mind:

- **Override a stage, not `pixelresponse()`, when you can.** Multi-read simulators (such as the NISP readout modes) call the stages directly and never call `pixelresponse()`. A model that only overrides `pixelresponse()` works for single-read exposures (`DetectorSimulator::simulate()`) but not in ramps.
- **`integrate()` only sees the charge of one interval.** Effects that depend on the *total* accumulated charge (non-linearity, full-well saturation) need the accumulated charge. The instrument's readout loop accumulates it and passes it to `read()`, so these effects belong there: `nisp::pixel::read(charge, nl)` limits the charge to the full well and applies the non-linear response, while the accumulated charge itself stays linear ([NISP example §3.1](NISP-example#31-one-pixel)) - An other possibility would have been to change pixel response model with non linear function and make sure pixel is tracking its accumulated charge along integration time. This would have been more complex and less generic, so it was not chosen.
- **Pixels with state** Some physical effects, such as charge trapping in persistence models, depend on the history of individual pixels. To retain memory of its past behaviour across integration intervals and exposures, the associated state (e.g. trap occupancy) must be stored at the pixel level. It should not be stored in temporary variables inside integrate(), nor in shared global variables or thread-local storage because `phot2pix::DetectorSimulator` divides the detector into disjoint pixel ranges, each processed by a single thread. Consequently, each pixel is updated by only one thread at a time, avoiding the need for thread synchronisation when accessing pixel-specific state.

---

## 6. Performance principles

- **Parallelism is built into the detector.** `simulate()` splits the pixel range across `std::thread::hardware_concurrency()` threads once the array has more than 100 pixels per thread. Each thread makes **one** virtual call to `_engine()`, so the cost of virtual dispatch scales with the number of threads, not pixels.
- **No allocation in per-pixel loops.** Buffers are allocated once and reused; output cubes are allocated at their final size instead of being grown layer by layer.
- **Fused passes over memory.** For example, a NISP ramp step adds the charge, applies the read and converts to `uint16` in a single threaded loop.

---

## 7. Testing philosophy

phot2pix come with serrry of unitary test that validate the modelling (including validation of the NISP example). However, noise simulator as phot2pix cannot be tested value by value. The phot2pix unit tests are **statistical** and may occasionally fail due to random fluctuations, particularly when the random number generator is not seeded. The tests are designed to be **robust**, with **sufficiently conservative tolerances** to minimise such false failures. However, occasional failures due to statistical fluctuations cannot be entirely ruled out, even when the code is correct.

Tests use [GoogleTest](https://github.com/google/googletest) and run with `ctest`. Files written by tests for inspection go to `./build/` and are replaced on every run. [GoogleTest](https://github.com/google/googletest) is provided as gitsubmodule.

---

*The Phot2Pix logo has been generated by OpenAI's ChatGPT Image Generation tool (ImageGen).*

*This wiki page has been generated by AI (Claude Opus 5.5) by reading the code itself.*

<!-- This page is generated from wiki/Philosophy.md in the phot2pix repository. Edit it there: changes made in the web interface are overwritten on the next push. -->
