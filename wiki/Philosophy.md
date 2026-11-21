# Philosophy of phot2pix

This page explains how phot2pix is organised and **why**. It is not an API reference; the headers in `src/phot2pix/` hold the details.

---

## 1. Separate what changes for different reasons

A detector simulation mixes very different things: the physics of one pixel, the geometry of an array, the readout sequence of a given instrument, its configuration, and the observing sequence. phot2pix gives each of them its own layer, so that one can be changed without touching the others.

| Layer | Class | Responsible for | Knows nothing about |
|---|---|---|---|
| **Configuration** | `PixelConfig` | Parameters, from JSON files, inline lists or per-pixel FITS maps | Physics |
| **Pixel** | `pixsim` and subclasses | The physics of **one** pixel: reset, charge collection, read | Neighbouring pixels, time sequences, readout patterns |
| **Detector** | `DetectorSimulator<PixelT>` | The pixel array, per-pixel maps, multi-threading, FITS output, effects that couple pixels | Readout patterns |
| **Instrument** | e.g. `nisp::SCA` | The readout modes, electronics noise, onboard processing, output headers | — |
| **Scene** *(planned)* | `ScenerySequencer` | The time loop over exposures, variable and moving sources | — |

The core (`src/phot2pix/`) contains the first three layers and is **instrument-agnostic**. An instrument lives in its own directory, as `src/nisp/` does, and builds on the core by inheritance.

---

## 2. The pixel: reset, integrate, read

### The physics

A modern image sensor pixel (an infrared HgCdTe pixel, a CMOS pixel) is read **non-destructively**:

1. **Reset:** the pixel is set to a reference level. That level is not exactly the same at every reset (kTC or reset noise), and it stays fixed until the next reset.
2. **Integration:** photons and dark current deposit charge, which **accumulates**. Its fluctuation is Poisson (shot) noise, so it accumulates too.
3. **Read:** the electronics measure the charge. Each measurement adds its own electronic noise (read noise, 1/f noise), and **the measurement does not change the charge**. Reading the same pixel twice gives two independent noise draws on the same charge.

### The model

`phot2pix::pixsim` implements exactly these three stages, one virtual method each:

| Stage | Method | Returns (ADU) | Default physics |
|---|---|---|---|
| Reset | `resetLevel()` | Reset level | baseline + N(0, kTC)/gain |
| Integration | `integrate(flux, t)` | Charge collected during *t* | N(S, √S)/gain, with S = (flux·QE + dark)·t electrons |
| Read | `read(charge)` | Measured value of a given charge | charge + N(0, RN)/gain |

A **single-read exposure** chains the three: `pixelresponse(flux, t) = read(resetLevel() + integrate(flux, t))`. Its mean is S/g + B and its variance (S + RN² + kTC²)/g².

A **multi-read sequence** (a ramp) calls the stages directly:

```
read 1 :  charge  = resetLevel() + integrate(flux, t_frame)    value_1 = read(charge)
read k :  charge += integrate(flux, t_frame)                   value_k = read(charge)
```

### Why pixels are simulated on their own

phot2pix starts from the assumption that **pixels are independent**: the response of a pixel depends only on the flux it receives, on its own parameters and on the time it integrates. It does not depend on its neighbours or on what happened before. Everything that breaks this assumption is added *on top* of the pixel simulation, by the layer that has the information it needs.

**Array effects** couple pixels together. Examples are inter-pixel capacitance, crosstalk, the 1/f noise shared by all the pixels of one output channel, and column or amplifier glow. These effects need to know the geometry of the array and the way it is read out, so they belong to the detector or instrument layer, not to the pixel.

**History effects** depend on time. Examples are the charge accumulated along a ramp, the persistence released by traps filled during previous exposures, and the readout sequence itself. These effects need a memory of the past, which is added on top of the pixel response. In the NISP example, the readout loop keeps the accumulated charge and the persistence scheduler keeps the past exposures.

This separation is what makes the simulator understandable:

- **Each layer can be validated on its own.** The pixel is checked against its analytic mean and variance. Each array or history effect is then checked as a difference from the independent-pixel result.
- **Each effect lives in one place**, so reading or changing it does not require understanding the rest.
- **Independent pixels can be simulated in parallel** with no synchronisation, and that is where most of the computing time goes.

The three pixel stages follow from the same reasoning applied to time. A ramp is a history: the charge accumulates from one read to the next, but the noise of a read does not. Because the read is separate from the charge, the readout loop can accumulate `integrate()` and call `read()` at each read. Read noise then stays at the same level on every read instead of growing along the ramp.

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

---

## 4. Configuration-driven

All parameters enter through a `PixelConfig`, a typed key/value bag that can be built:

- **from a JSON file:** `PixelConfig cfg("detector.json");`
- **inline:** `PixelConfig cfg{{"gain", 1.9}, {"readNoise", 13.0}};`
- **with per-pixel maps:** any standard key can be a FITS file path instead of a number, with an optional `"<key>.hdu"` to select the HDU.

`DetectorSimulator::initialize(cfg)` understands the standard pixel keys (`gain`, `dark`, `readNoise`, `QE`, `KTCnoise`, `baseline`). Derived classes add their own keys by overriding `initialize()` and calling the base version first.

Two rules:

- **Configuration is parsed in the library, not in the tools.** A command-line tool, a test and a Python binding using the same config must get the same detector.
- **Read numbers with `getAsDouble()`.** JSON integers are stored as `uint64_t`, which is not the same type as `size_t` on every platform and `size_t` is integer size element of the array and time series; `getAsDouble()` accepts any numeric type and is converted afterward as `size_t` when needed.

---

## 5. Extending phot2pix: pick the right hook

| You want to… | Override | Example |
|---|---|---|
| Change **one** physical effect of a pixel | The stage it belongs to: `integrate()`, `resetLevel()` or `read()` | Cosmic-ray hits or charge trapping → `integrate()`; non-Gaussian read noise → `read()` |
| Replace the whole single-read model | `pixelresponse()` (virtual) | A completely different detector technology |
| Add an effect that couples pixels | `DetectorSimulator::_engine()` | Inter-pixel capacitance, crosstalk, column glow |
| Add configuration keys | `initialize(cfg)`, calling the base first | NISP: `exposureMode`, `frame_time`, `oneOverF`… |
| Record the inputs in the output header | `setInputCard()` in `initialize()` (`PixelConfig::card(key)` gives a key's value or map file name), then `writeInputCards()` on the output. The base already records the six pixel parameters (`PIXGAIN`, `PIXDARK`, …) | NISP: readout noise, persistence, reference pixels |
| Add a readout pattern | An instrument class that sequences the stages | NISP: MACC, UTR, CDS, DPU |

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
- **`integrate()` only sees the charge of one interval.** Effects that depend on the *total* accumulated charge (non-linearity, full-well saturation) need the accumulated charge; that is the job of the layer that accumulates it, the instrument's readout loop.
- **Pixels with state** (for example trap occupancy in a persistence model) must keep that state per pixel. `DetectorSimulator` splits the pixels into disjoint ranges per thread, so each pixel is touched by one thread only.

---

## 6. Performance principles

- **Parallelism is built into the detector.** `simulate()` splits the pixel range across `std::thread::hardware_concurrency()` threads once the array has more than 100 pixels per thread. Each thread makes **one** virtual call to `_engine()`, so the cost of virtual dispatch scales with the number of threads, not pixels.
- **No allocation in per-pixel loops.** Buffers are allocated once and reused; output cubes are allocated at their final size instead of being grown layer by layer.
- **Fused passes over memory.** For example, a NISP ramp step adds the charge, applies the read and converts to `uint16` in a single threaded loop.
- **Measure before optimising.** The `nispSimulator` tool has a `"benchmark"` mode that times the simulation, per run and per layer.

---

## 7. Testing philosophy

A simulator of noise cannot be tested value by value. phot2pix tests are **statistical**:

- **One effect at a time:** every other effect is switched off (gain 1, QE 1, no read noise, no kTC…) so that the expected mean and variance are known exactly.
- **Many samples:** means and variances are measured over thousands of draws or millions of pixels and compared with the expectation within **5σ**, plus a small margin for integer rounding.
- **Each stage on its own, then the composition:** for example, `integrate()` must contain no read noise; a ramp must show RN² on every read and 2·RN² on every read-to-read difference; a CDS must cancel the reset level.
- **Full scale:** the NISP readout modes are validated on the 2040 × 2040 active area of the detector against analytic expectations; see [NISP example → Validation](NISP-example#9-validation).

Tests use GoogleTest and run with `ctest`. Files written by tests for inspection go to `./build/` and are replaced on every run.

<!-- This page is generated from wiki/Philosophy.md in the phot2pix repository. Edit it there: changes made in the web interface are overwritten on the next push. -->
