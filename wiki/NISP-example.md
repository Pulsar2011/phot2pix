# NISP example: the Euclid/NISP HAWAII-2RG detectors

`src/nisp/` is a modelling of the Near Infrared Detector of the **NISP** (Near-Infrared Spectrometer and Photometer) instrument which is built on the `phot2pix` core. It shows how an realistic detector is assembled from the generic building blocks described in [Philosophy](Philosophy).

---

## 1. The instrument in a few lines

- NISP has 16 **Teledyne HAWAII-2RG** (H2RG) HgCdTe arrays. Each SCA (sensor chip assembly) is 2048 × 2048 pixels of 18 µm, including a 4-pixel border of reference pixels (2040 × 2040 sensitive). Reference pixels have no photodiode but are read like the others, so they record the drifts of the readout electronics.
- An SCA is read **non-destructively**, one full frame every **1.45408 s**, through **32 output channels** of 64 columns each, read in parallel. One channel reads 2048 × 64 pixels per frame, one pixel every 11.1 µs.
- Reads are taken in **MACC(ng, nr, nd)** mode: *ng* groups of *nr* consecutive reads, separated by *nd* dropped frames.
  - Photometry uses MACC(4,16,4) (87.2 s between the first and last group).
  - Spectroscopy uses MACC(15,16,11) (549.6 s).
- The onboard **DPU** groups read into groups to compute a per group average before it fits a slope to the group averages ([Kubik et al. 2016, PASP 128, 104504](https://iopscience-iop-org.lama.univ-amu.fr/article/10.1088/1538-3873/128/968/104504)). Only that slope image is sent to the ground.

While NISP instrument use fixed exposure parameter this model allows to simulate any exposure time and readout mode, with the same physics as the real instrument. The model includes:
- **per-pixel calibration maps** of gain, dark, read noise, quantum efficiency, kTC noise and baseline; Configuration can be a number or a FITS map of the matrix. When number is provided the same value is used for all pixels.
- **full well and non-linear response** of the pixels. Non linear response is given as list of polynomial coeficients.
- **inter-pixel capacitance** (IPC) and **crosstalk** between channels, applied at each read. **IPC** can be a 3D cube from a fits file. In that case, the convolution is performed dynamically with the kernel associated to each pixel. **Crosstalk** is a list of couplings between channels as a vector representing the coupling depending on the distance separating two channels. The crosstalk is applied to the signal after the IPC convolution.
- **correlated readout noise** of the output channels, including 1/f noise per channel, bias voltage noise common to all channels and pedestal drift per channel. The readout noise is generated as a continuous stream over the ramp.
- **persistence** of the pixels, with a configurable number of traps and decay times. The persistence is added to the charge at each step of the ramp.
- **reference pixel correction** and **on board NISP DPU processing** as well as various acquisition mode not used by NISP instrument which are native to the H2RG detectors.

---

## 2. Architecture

### 2.1 Map of the code

| File | Content |
|---|---|
| `../phot2pix/PixelConfig.h` | `phot2pix::PixelConfig`: the typed key/value configuration, read from JSON. `HeaderCard` and `PixelConfig::card(key)`: an input written as a FITS header card |
| `../phot2pix/Geometry.h` | `phot2pix::Geometry`: the physical detector, the placement of the simulated matrix on it, and the mapping of detector-sized images onto the matrix |
| `../phot2pix/PixelSimulator.h/.cxx` | `phot2pix::pixsim`: the three stages of a pixel (reset, integration, read) |
| `../phot2pix/DetectorSimulator.h` | `phot2pix::DetectorSimulator<PixelT>`: the pixel array, its geometry, per-pixel FITS maps, threaded simulation, and the input header cards |
| `../phot2pix/DetectorEffect.h/.cxx` | `phot2pix::PixelConvolution`: the convolution of an image by a kernel, the same for every pixel or one per pixel. Used for the inter-pixel capacitance and for the crosstalk between channels. This convolution isn't applyed onto `phot2pix.h` and kernel and when to apply this convolution need to be defined in derived model. In this example, `nisp::SCA` apply IPC (Inter Pixel Capacitance) and crosstalk convolution during the readout. |
| `../phot2pix/DetectorTools.h` | `phot2pix::parallelFor()` and `withCubeData()`: threaded loop and typed access to a FITS cube, shared by the detectors and their effects |
| `../phot2pix/ScenerySequencer.h/.cxx` | `phot2pix::ScenerySequencer`: a sequence of exposures on one detector ([Sequences of exposures](ScenerySequencer)) |
| `nispPixel.h/.cxx` | `nisp::pixel`: the H2RG pixel, with its full well and non-linear response. `nisp::nonlinearity`: the models of that response. `nisp::Ramp`: the DPU slope estimator (group co-adding, Kubik slope, Fowler, CDS) |
| `nispSCA.h/.cxx` | `nisp::SCA`: the detector, its readout modes, the inter-pixel capacitance and crosstalk of each read, the reference pixel correction and its output header. `nisp::macc`: MACC parameters. `nisp::readout`: the parameters of one exposure. `nisp::persistence`: the persistence model (`model`, `pixelArray`, `scheduler`) |
| `nispScenerySequencer.h/.cxx` | `nisp::ScenerySequencer`: a sequence of NISP exposures, each with its own MACC parameters or exposure time and readout mode ([details](ScenerySequencer#3-the-nisp-override)) |
| `nispGeometry.h` | `nisp::Geometry`: the H2RG layout (2048 × 2048, reference border, output channels), derived from `phot2pix::Geometry` |
| `nispOneOverFNoise.h/.cxx` | `nisp::OneOverFNoise`: readout noise of the output channels (1/f noise per channel, bias voltage noise common to all channels, pedestal drift) and its continuous `Stream` over a ramp |
| `../utilities/nispSimulator.cxx` | Command-line tool driven by a JSON file |
| `../utilities/*.json` | A default configuration with averaged NISP values, and one test configuration per readout mode |

### 2.2 Classes

Class diagram of the NISP model, with the main classes and their relations. The arrows show the direction of the dependency, not the direction of the data flow. The `DetectorSimulator` base class is templated on the pixel type, and `SCA` inherits it **privately** to hide the pixel type from the outside.

```mermaid
classDiagram
    pixsim <|-- pixel : nisp
    pixsim <|-- model : persistence, private
    Geometry <|-- NispGeometry : nisp, reference border, channels
    DetectorSimulator~pixel~ <|-- SCA : private
    DetectorSimulator~model~ <|-- pixelArray
    DetectorSimulator o-- Geometry : placement, maps
    DetectorSimulator ..> PixelConfig : initialize(cfg)
    DetectorSimulator o-- HeaderCard : input cards
    SCA ..> NispGeometry : constructed with
    SCA o-- OneOverFNoise : readout noise
    OneOverFNoise o-- NispGeometry : copy
    OneOverFNoise ..> Stream : beginStream()
    SCA o-- PixelConvolution : IPC, crosstalk
    SCA o-- functor : non-linearity model
    pixel ..> functor : read(charge, nl)
    SCA o-- scheduler : persistence
    scheduler o-- pixelArray
    SCA ..> Ramp : DPU slope
```

`NispGeometry` stands for `nisp::Geometry`, `functor` for `nisp::nonlinearity::functor`.

How the layers map onto the [philosophy](Philosophy):

- **Pixel:** `nisp::pixel` uses the reset and integration stages of `pixsim` unchanged, and adds the full well and the non-linear response to the read: `read(charge, nl)`. Other NISP-specific pixel physics would go there too, by overriding a stage.
- **Detector:** `nisp::SCA` inherits `DetectorSimulator<nisp::pixel>` **privately**, gives it a `nisp::Geometry`, and adds its own methods on top of it to model the behaviour specific to the NISP H2RG.
  - Used as is from `DetectorSimulator`: the geometry, the per-pixel maps (`loadMap()`, `loadImage()`), the header cards that describe them, and the multi-threaded dispatch `simulate()`. The members that `nisp::SCA` keeps public are re-exported with `using`.
  - Overridden: `initialize()`, which calls the base version first and then reads the NISP keys, and the `_engine()` hooks that `simulate()` dispatches to. The base hooks return a complete read of each pixel (`pixelresponse()`); those of `nisp::SCA` return its **charge only** (`resetLevel() + integrate()`), which the readout modes then read with `pixel::read()`.
  - On top of that, `nisp::SCA` adds its own methods to model the detector behaviour proper to the NISP H2RG: it marks the reference pixels, sequences the reads of each readout mode, applies the inter-pixel capacitance and the crosstalk between channels, adds the read noise (through `pixel::read()`) and the readout noise of the output channels at each read, runs the DPU processing (reference pixel correction, then slope) and writes the output header.
  - Its only public simulation entry point is `expose()`. The reference pixel correction is also public, as two static functions ([§5](#reference-pixel-correction-of-a-ramp)), and so are the inter-pixel capacitance and the crosstalk of an image ([§5](#inter-pixel-capacitance-and-crosstalk-of-an-image)).

### 2.3 The pipeline

```mermaid
flowchart TD
    JSON["JSON configuration"] --> CFG["phot2pix::PixelConfig"]

    subgraph INIT ["SCA::initialize(cfg)"]
        direction TB
        I1["① geometry, in DetectorSimulator::initialize<br/>detectorSize, detectorOrigin, refBorder, colPerChannel, nChannels<br/>place the W × H matrix on the detector (throws if it does not fit)"]
        I2["② pixel model, in DetectorSimulator::initialize<br/>gain, dark, readNoise, QE, KTCnoise, baseline<br/>numbers, or FITS maps of matrix or detector size → PIX* cards"]
        I3["③ readout and DPU settings<br/>exposureMode, frame_time, offset, scale, dpu_gain, dpu_noise"]
        I4["④ persistence scheduler<br/>persistence traps"]
        I5["⑤ readout noise, OneOverFNoise<br/>oneOverF, bias_voltage_noise, pedestal_drift<br/>channel layout checked against the detector"]
        I6["⑥ reference pixels<br/>QE = 0, refDark"]
        I7["⑦ full well and non-linearity<br/>fullWell, nlModel, nl0, nl1, nl2: numbers or FITS maps<br/>model checked to increase up to the full well of every pixel"]
        I8["⑧ coupling of the pixels, PixelConvolution<br/>ipc: kernel, or FITS cube of kernels<br/>crosstalk: couplings between the channels of the geometry"]
        I9["⑨ input cards<br/>DETSIZE, REFBORD, REFDARK, full well, non-linearity, crosstalk, readout noise, persistence"]
        I1 --> I2 --> I3 --> I4 --> I5 --> I6 --> I7 --> I8 --> I9
    end

    CFG --> INIT
    FLUX["flux: number, or FITS image of matrix or<br/>detector size (nispSimulator, SCA::toMatrix)"] --> EXP
    INIT --> EXP["SCA::expose(flux, macc) or expose(flux, T)"]

    EXP --> MODE{"exposureMode"}
    MODE -->|"MACC, UTR"| RAMP["ramp: one charge step per frame,<br/>reads kept per group<br/>IPC, crosstalk, read noise, readout noise per read"]
    MODE -->|"DPU"| RAMP2["ramp, reference pixel correction,<br/>then Kubik slope per pixel"]
    MODE -->|"CDS"| CDS["read N − read 1"]
    MODE -->|"FLUX"| ONE["single read"]

    RAMP --> OUT["uint16 FITS image or cube<br/>+ header: mode, times, DPU, geometry, inputs"]
    RAMP2 --> OUT
    CDS --> OUT
    ONE --> OUT
    OUT --> TOOL["nispSimulator adds FLUX, FLUXSEC and writes the file"]
```

1. **Configuration.** Every parameter lives in one `PixelConfig`, normally read from a JSON file ([§4](#4-configuration-reference)). The library reads all the detector keys; `nispSimulator` reads only the matrix size, the flux and the readout pattern.
2. **Initialisation** (`SCA::initialize`, in the order of the diagram). The base class configures the geometry first, so that the FITS maps are cropped at the right place, then applies the pixel model. The SCA then reads its own settings, creates the persistence scheduler, configures the readout noise (and fails if the channel layout does not fit the detector), sets the reference pixels, then the full well and the non-linearity, then the inter-pixel capacitance and the crosstalk, and records the inputs as header cards.
3. **Exposure** (`SCA::expose`). The charge of every pixel is simulated step by step and read according to `exposureMode` ([§3.3](#33-the-algorithm), [§5](#5-readout-modes)).
4. **Output.** A `uint16` image or cube, whose header holds the readout description and every input of the model ([§6](#6-output-header)).

---

## 3. The model in equations

This section writes down what the code computes. Every ε is an independent draw from a standard normal N(0, 1): draws are independent between pixels, between stages and between reads.

### 3.1 One pixel

| Symbol | Meaning | Unit | Config key |
|---|---|---|---|
| $`\phi`$ | Incident flux | ph s⁻¹ | input array |
| $`\eta`$ | Quantum efficiency | — | `QE` |
| $`D`$ | Dark current | e⁻ s⁻¹ | `dark` |
| $`g`$ | Gain | e⁻ ADU⁻¹ | `gain` |
| $`B`$ | Baseline | ADU | `baseline` |
| $`\sigma_{\mathrm{kTC}}`$ | Reset (kTC) noise | e⁻ | `KTCnoise` |
| $`\sigma_{\mathrm{RN}}`$ | Read noise | e⁻ | `readNoise` |
| $`Q_{\mathrm{FW}}`$ | Full well | e⁻ | `fullWell` |
| $`c_2, c_3`$ | Non-linearity (`polynomial`) | (e⁻)⁻¹, (e⁻)⁻² | `nl0`, `nl1` |
| $`t_f`$ | Frame time | s | `frame_time` |

**The three stages** (`phot2pix::pixsim`). Every stage returns ADU:

```math
\begin{aligned}
\text{reset:}\quad & R = B + \frac{\sigma_{\mathrm{kTC}}}{g}\,\varepsilon_R && \texttt{resetLevel()}\\[6pt]
\text{integration over } \Delta t\text{:}\quad & I(\Delta t) = \frac{S + \sqrt{S}\,\varepsilon_I}{g}, \qquad S = (\eta\,\phi + D)\,\Delta t && \texttt{integrate}(\phi, \Delta t)\\[6pt]
\text{read of a charge } Q\text{:}\quad & M(Q) = Q + \frac{\sigma_{\mathrm{RN}}}{g}\,\varepsilon_M && \texttt{read}(Q)
\end{aligned}
```

$S$ is the mean number of electrons collected, photo-electrons and dark current together. Its shot noise $\sqrt{S}$ is the Gaussian approximation of Poisson. $I(\Delta t) = 0$ only when $S = 0$, that is with no light **and** no dark current ($\eta\,\phi = 0$ and $D = 0$), or for $\Delta t = 0$: a dark exposure ($\phi = 0$) still integrates the dark current, with its shot noise.

A read does not measure the mean signal $S$ of one integration, but the charge $Q$ the pixel actually holds when it is read, in ADU: the reset level $R$ plus every integration $I(\Delta t_j)$ since the reset, each with its own noise draw (both defined below). For read $k$ after a reset:

```math
Q_k = R + \sum_{j=1}^{k} I(\Delta t_j)
```
A single read after one integration of duration $t$ is $Q = R + I(t)$. $Q$ accumulates from one read to the next; the read noise of $`M(Q)`$ does not. In the NISP readout, the persistence is also added to $Q$ ([§4.8](#48-persistence)), and the full well and the non-linearity are applied to it at each read (below).

**Full well and non-linearity** (`nisp::pixel::read(charge, nl)`, used by every readout mode). The accumulated charge *Q* stays linear; each read converts it to electrons above the baseline, limits it to the full well, and applies the non-linear response *P* before the read noise:

```math
q = \min\big(g\,(Q - B),\ Q_{\mathrm{FW}}\big),
\qquad
P(q) = q\,\big(1 + c_2\,q + c_3\,q^2\big),
\qquad
M(Q) = B + \frac{P(q)}{g} + \frac{\sigma_{\mathrm{RN}}}{g}\,\varepsilon_M
```

- *P* is the `polynomial` model. A model is any function of the charge and of the parameters of the pixel ([§5](#non-linearity-model)). P(q) = q for q ≤ 0 (kTC noise can put the charge below the baseline), and without a model.
- The kTC noise and the persistence are part of the charge, so they share the well. The baseline is an electronic offset: it does not fill it.
- The read noise and the readout noise are added after the response: they are not compressed. The shot noise is: the variance of a read is (P′(q)²·q + σ_RN²)/g² below the full well, and σ_RN²/g² above.

> **⚠ Warning: the non-linearity needs a full well.** Without a finite full well, `nisp::pixel::read(charge, nl)` is the linear read M(Q) above: both the full well **and** the non-linear response are skipped, whatever `nl` is. By default there is no full well (`fullWell` not given), so:
> - through `nisp::SCA`, a non-linearity model (`nlModel`, or `setNonLinearity()`) is refused unless every pixel has a finite `fullWell`: `initialize()` or `setNonLinearity()` throws `std::invalid_argument`;
> - `nl0`, `nl1`, … given without `nlModel` give a linear response, with no error;
> - a code that calls `nisp::pixel::read(charge, nl)` itself, outside `nisp::SCA`, gets a linear read with no error when the full well of the pixel was not set (`setFullWell()`).

The statistics below are those of the linear pixel, far from the full well.

**A single read** (`pixelresponse()`, the FLUX mode) chains the three stages:

```math
y = M\big(R + I(T)\big) = B + \frac{S + \sqrt{S}\,\varepsilon_I + \sigma_{\mathrm{kTC}}\,\varepsilon_R + \sigma_{\mathrm{RN}}\,\varepsilon_M}{g},
\qquad
\langle y\rangle = B + \frac{S}{g},
\qquad
\mathrm{Var}(y) = \frac{S + \sigma_{\mathrm{kTC}}^2 + \sigma_{\mathrm{RN}}^2}{g^2}
```

**A ramp** reads the pixel once per frame without resetting it. With $`S_f = (\eta\,\phi + D)\,t_f`$ the mean charge collected in one frame, the charge after frame *n* and the value of read *n* are

```math
Q_1 = R + I_1, \qquad Q_n = Q_{n-1} + I_n, \qquad y_n = M(Q_n) = Q_n + \frac{\sigma_{\mathrm{RN}}}{g}\,\varepsilon_n
```

where each $`I_n`$ is an independent draw of $`I(t_f)`$. The charge carries the memory of the ramp and the read noise does not. The statistics of the ramp follow:

```math
\langle y_n\rangle = B + n\,\frac{S_f}{g},
\qquad
\mathrm{Cov}(y_m, y_n) = \frac{\sigma_{\mathrm{kTC}}^2 + \min(m,n)\,S_f + \delta_{mn}\,\sigma_{\mathrm{RN}}^2}{g^2}
```

so that

```math
\mathrm{Var}(y_n) = \frac{\sigma_{\mathrm{kTC}}^2 + n\,S_f + \sigma_{\mathrm{RN}}^2}{g^2},
\qquad
\mathrm{Var}(y_{n+1} - y_n) = \frac{S_f + 2\,\sigma_{\mathrm{RN}}^2}{g^2},
\qquad
\mathrm{Var}(y_N - y_1) = \frac{(N-1)\,S_f + 2\,\sigma_{\mathrm{RN}}^2}{g^2}
```

The last relation is the CDS: the reset level *R* and its kTC noise cancel. These relations are the ones checked by `gTest_nispSCA.cxx` ([§9](#9-validation)).

### 3.2 The full array

Pixel *i* of the W × H matrix sits at column *x* and row *y*, with i = y·W + x. Every parameter of §3.1 can be a per-pixel FITS map: $`\eta_i, D_i, g_i, B_i, \sigma_{\mathrm{kTC},i}, \sigma_{\mathrm{RN},i}, Q_{\mathrm{FW},i}, c_{2,i}, c_{3,i}`$. Pixels remain independent, except through the terms that the SCA adds on top of the pixel model:

- the **inter-pixel capacitance** (IPC) and the **crosstalk** between output channels, which couple the pixels at each read and leave the charge unchanged. Both are a convolution of the image by a kernel (`phot2pix::PixelConvolution`): for the IPC, the kernel $`K`$ of the `ipc` key, the same for every pixel or one per pixel; for the crosstalk, a one-row kernel holding the couplings $`\kappa_d`$ at ± *d* channel widths, applied to the signal above the baseline and added to the image ([details](#inter-pixel-capacitance-and-crosstalk));
- the **readout noise** $`\mathcal{R}`$, in electrons, added at each read: the 1/f noise of each output channel, modelled as a sum of $`n_p`$ AR(1) processes whose spectrum approximates PSD ∝ 1/ν^a, with amplitude *A* and slope *a* ([details](#1f-noise-of-channel-c)), the bias voltage noise shared by all channels and the pedestal drift of each channel;
- the **persistence** *P*, which depends on the previous exposures and is added to the charge.

A quantiser then converts each value to `uint16`:

```math
\mathcal{Q}(x) = \min\big(65535,\ \max(0,\ \lfloor x + \tfrac{1}{2} \rfloor)\big),
\qquad
\mathcal{Q}_{so}(x) = \mathcal{Q}(s\,x + o)
```

with *s* = `scale` and *o* = `offset`.

#### Steps of a ramp

An exposure is simulated as a sequence of charge steps $m = 1 … M$ of duration $\Delta t_m$. In MACC($n_g$, $n_r$, $n_d$), each group contributes $n_r$ exposure step of $t_f$ integration time, each followed by a read. Between two groups there is one exposure step of $n_d·t_f$ integration time without any read (the dropped frames). Read *r* of group *g* is stored in layer $\ell = g\times n_r + r$ and taken at frame

```math
n_\ell = g\,(n_r + n_d) + r + 1
```

that is, n_ℓ·t_f after the reset. UTR is MACC(N, 1, 0).

#### Charge

```math
Q_i[m] = R_i + \sum_{m'=1}^{m} \Big( I_i(\Delta t_{m'}) + P_i[m'] \Big)
```

$`R_i`$ is drawn at the first step only; the later steps are run with `baseline = "off"`.

#### Inter-pixel capacitance and crosstalk

Both act on the signal at each read, in every readout mode, and never on the charge *Q* that the pixels hold. They are off by default.

The **inter-pixel capacitance** shares the signal of a pixel with its neighbours. $`K_j(\delta)`$ is the fraction of the signal of pixel *j* that appears at the position j + δ, with $`\sum_\delta K_j(\delta) = 1`$. The kernel is either the same for every pixel or one per pixel:

```math
V_i = \sum_\delta K_{i-\delta}(\delta)\; Q_{i-\delta}
```

Each neighbour gives its signal with **its own** kernel. A pixel at the border of the image misses neighbours: its sum is divided by the sum of the weights it has, so that the level of a uniform baseline is kept.

The **crosstalk** adds, to the pixels read at the same time in the other channels, a fraction $`\kappa_d`$ of the signal above the baseline that a pixel outputs. These pixels are on the same row, a multiple of the channel width *w* away, since every channel is read in the same direction. *d* is the distance between the channels:

```math
\tilde Q_i = V_i + \sum_{d \ge 1} \kappa_d\,\big(s_{i - d\,w} + s_{i + d\,w}\big),
\qquad
s_j = \frac{P_j\big(\min(g_j\,(V_j - B_j),\ Q_{\mathrm{FW},j})\big)}{g_j}
```

The source keeps its own signal (no loss), and a saturated pixel cannot give more than its full well. $`P_j`$ is the non-linear response of §3.1, not the persistence. Positions outside the image give nothing.

The read is then that of the coupled signal: **everywhere below, $`M_i(Q_i)`$ stands for $`M_i(\tilde Q_i)`$.**

#### MACC and UTR: one layer of the output cube

```math
Y_i[\ell] = \mathcal{Q}\Big( M_i\big(Q_i[m_\ell]\big) + \frac{\mathcal{R}_i[n_\ell]}{g_i} \Big)
```

where $`M_i`$ is the read of §3.1 (full well, non-linearity, read noise) and $`m_\ell`$ is the step that ends at frame $`n_\ell`$. MACC and UTR layers are rounded and clamped but **not** scaled or offset.

#### Physical detector and placement of the matrix

The readout noise, the crosstalk and the reference pixels are defined on the **physical detector** (`nisp::Geometry`), not on the simulated matrix: $`W_d \times H_d`$ = 2048 × 2048 pixels (`detectorSize`), a border of b = 4 reference pixels (`refBorder`), and $`W_d / w`$ = 32 channels (`nChannels`) of w = 64 columns (`colPerChannel`), channel *c* covering the physical columns [c·w, (c+1)·w). The W × H matrix is a part of that detector, placed at (x₀, y₀) (`detectorOrigin`, [0, 0] by default: the first matrix pixel is the first detector pixel):

```math
X_i = x_0 + x, \qquad Y_i = y_0 + y
```

Pixel *i* belongs to channel $`c(i) = \lfloor X_i / w \rfloor`$, and is a reference pixel when $`X_i`$ or $`Y_i`$ is within *b* of the detector edge. A 2048 × 2048 matrix is the full detector, a 2040 × 2040 one at [4, 4] its active area, a 2044-wide one at x₀ = 2 keeps 2 reference columns on each side, and a 256-wide one at x₀ = 896 covers the channels 14 to 17. At the default origin, a small matrix has the 4 reference columns on its left and the 4 reference rows at its bottom. A matrix that does not fit on the detector throws at `initialize()`. Since everything is derived from physical coordinates, a smaller matrix has the same channels, timing and noise spectrum as the corresponding part of the full detector: it is only cheaper.

The geometry also maps the per-pixel inputs onto the matrix (`phot2pix::Geometry`): a FITS map (`gain`, `dark`, `readNoise`, `QE`, `KTCnoise`, `baseline`) is either W × H, used as is, or $`W_d \times H_d`$, in which case the section [x₀, x₀ + W) × [y₀, y₀ + H) is used. Any other size throws at `initialize()`. A full-detector calibration map can therefore be used to simulate any part of the detector.

#### Readout noise

All channels are read in parallel, row by row, with no row or frame overhead. Within a frame, pixel *i* is read at pixel-clock sample $`t_i = Y_i\,w + (X_i - c(i)\,w)`$ of its channel, one frame being $`L = H_d\,w`$ samples (one sample every $`t_f / L`$ = 11.1 µs). The readout noise of pixel *i* in frame *n* is the sum of three components, in electrons:

```math
\mathcal{R}_i[n] =
\underbrace{\mathcal{N}_{c(i)}\big[(n - 1)\,L + t_i\big]}_{\texttt{oneOverF}}
+ \underbrace{\mathcal{C}\big[(n - 1)\,L + t_i\big]}_{\texttt{bias\_voltage\_noise}}
+ \underbrace{f\big(t_0 + t_{n,Y_i}\big) + \sigma_P\,\varepsilon_{c(i),n}}_{\texttt{pedestal\_drift}},
\qquad
t_{n,Y} = \Big(n - 1 + \frac{Y + \frac12}{H_d}\Big)\,t_f
```

$`t_{n,Y}`$ is the time since the reset at which physical row *Y* of frame *n* is read. Every component is off by default.

#### 1/f noise of channel *c*

The noise of channel *c* is a single time series over the whole ramp, a sum of n_p AR(1) processes:

```math
\mathcal{N}_c[t] = A \sum_{p=0}^{n_p-1} z_{c,p}[t],
\qquad
z_{c,p}[t] = \rho_p\, z_{c,p}[t-1] + \sigma_p\,\varepsilon_{c,p,t},
\qquad
z_{c,p}[0] \sim \mathcal{N}\Big(0,\ \frac{\sigma_p^2}{1 - \rho_p^2}\Big)
```

```math
\nu_p = \nu_{\min} \left(\frac{0.5}{\nu_{\min}}\right)^{p/(n_p-1)},
\quad
\nu_{\min} = \frac{1}{2 L},
\qquad
\rho_p = e^{-2\pi\nu_p},
\qquad
V_p = \frac{\sigma_p^2}{1 - \rho_p^2} \propto \nu_p^{\,1-a},
\quad
\sum_p V_p = 1
```

The corner frequencies ν_p are in cycles per pixel-clock sample and log-spaced (ν_0 = ν_min when n_p = 1). *A* is `amplitude` and *a* is `slope`. Each pole carries the variance $`V_p`$ that a spectrum ∝ 1/ν^a holds in its log step of frequency, so together the poles approximate that spectrum, with a stationary variance A². The processes start from their stationary state, since the amplifier does not reset with the pixels. The AR(1) processes of different channels are driven by independent draws, so channels are uncorrelated. [§4.7](#47-readout-noise) explains how to choose *A*, *a* and n_p.

Only the samples of the matrix are generated: the samples of the rows and columns outside it are skipped exactly, since an AR(1) process advances over *n* samples in one step,

```math
z[t + n] = \rho^n\, z[t] + \sigma \sqrt{\frac{1 - \rho^{2n}}{1 - \rho^2}}\;\varepsilon
```

so the cost scales with W × H.

#### Bias voltage noise

$`\mathcal{C}`$ is one series built like $`\mathcal{N}_c`$ (amplitude, slope and poles of `bias_voltage_noise`), but **shared by every channel**: it is indexed by the same pixel-clock sample $`t_i`$ in every channel. The pixels read at the same instant in all channels get the same value: this is the horizontal banding across the full width, repeating every 64 columns, that the left/right reference columns measure.

#### Pedestal drift

Each channel is offset by the same amount over all its pixels: a mean drift *f*, common to all channels and evaluated at the read time of each row, plus $`\sigma_P\,\varepsilon_{c,n}`$, one independent draw per channel and per frame (read or dropped). The top/bottom reference rows of a channel measure it.

```math
f(t) = \sum_k a_k\,t^k + \sum_j A_j\,e^{-t/\tau_j} + \sum_m B_m \sin\!\Big(\frac{2\pi t}{P_m} + \varphi_m\Big)
```

$`t_0`$ places the exposure in time (for example after a power-on or a mode change). A constant offset belongs to the `baseline` map instead: it is drawn at the reset and cancels in every difference, whereas the pedestal drift changes from one read to the next.

#### Persistence

Each past exposure *j* leaves a stored charge $`E_i^{(j)} = \eta_i\,\phi_i^{(j)}\,T^{(j)}`$, in electrons, where $`T^{(j)}`$ is the integration time registered for that exposure (see [§10](#10-known-limitations) for MACC). Trap *k* releases it at the rate $`E\,f_k\,e^{-t/\tau_k}`$, *t* being the time since the end of exposure *j*. The charge released during step *m*, which starts $`\Delta_j[m]`$ after the end of exposure *j* and lasts $`\Delta t_m`$, is the integral of that rate over the step:

```math
P_i[m] = \sum_j E_i^{(j)} \sum_k f_k\,\tau_k\, e^{-\Delta_j[m]/\tau_k} \left( 1 - e^{-\Delta t_m/\tau_k} \right),
\qquad
\Delta_j[m] = \sum_{m' < m} \Delta t_{m'} + n_j\,\delta
```

$`f_k`$ is therefore a release rate (s⁻¹) per stored electron, and $`f_k\,\tau_k`$ the fraction of the stored charge that trap *k* releases in total ([§4.8](#48-persistence)); consecutive steps add up to the release over their total duration. The elapsed time $`\Delta_j[m]`$ counts the simulated steps since exposure *j*, plus the delay δ = `dT` between two exposures, once for each of the $`n_j`$ exposures that have ended since *j* (exposure *j* included). $`P_i[m]`$ is added to the charge as it is, without conversion by $`g_i`$. $`E_i^{(j)}`$ is not limited by the full well ([§10](#10-known-limitations)).

#### FLUX and CDS

Both modes read the array once or twice. Each read draws a new, independent readout-noise realisation whose AR(1) processes start from their stationary state ($`\mathcal{R}`$ and $`\mathcal{R}'`$ below). The pedestal drift is evaluated at the frame of the read: the last frame for FLUX, frames 1 and N for CDS.

```math
\text{FLUX:}\qquad Y_i = \mathcal{Q}_{so}\Big( M_i\big(R_i + I_i(T') + P_i[1]\big) + \frac{\mathcal{R}_i[N]}{g_i} \Big),
\qquad T' = \lfloor T/t_f \rfloor\, t_f = N\,t_f
```

```math
\text{CDS:}\qquad
\begin{aligned}
Q_i^{(1)} &= R_i + I_i(t_f) + P_i[1], \qquad
Q_i^{(N)} = Q_i^{(1)} + I_i\big((N-1)\,t_f\big) + P_i[2], \qquad N = \lfloor T/t_f \rfloor\\[4pt]
Y_i &= \mathcal{Q}_{so}\Big( \Big[M_i(Q_i^{(N)}) + \frac{\mathcal{R}'_i[N]}{g_i}\Big] - \Big[M_i(Q_i^{(1)}) + \frac{\mathcal{R}_i[1]}{g_i}\Big] \Big)
\end{aligned}
```

In CDS, the pedestal drift leaves $`f(t_0 + t_{N,Y}) - f(t_0 + t_{1,Y})`$.

$`M_i`$ is the read of §3.1 with the parameters of pixel *i*. Each call draws a new ε. The windowed single read `expose(flux, t, c0, r0, nc, nr)` uses the FLUX expression on the window, with *T* not rounded to whole frames.

#### Reference pixel correction

It works on the full detector, one read at a time. The numbers below are those of NISP: 2048 × 2048 pixels, a border of 4 reference pixels, 32 channels of 64 columns. With other values of `detectorSize`, `refBorder` and `colPerChannel`, the same formula is used with them; only the sliding window of ± 4 lines is fixed (it stops at the edge of the detector).

Write $`Y[\ell](x, y)`$ for read ℓ at column *x* and line *y*. For output channel ch, which covers the columns [64·ch, 64·ch + 64), and for a line y ∈ [4, 2044), the correction is Eq. 3 of Euclid Collaboration (arXiv:2507.11326), from Kubik et al. (2014):

```math
c_y^{(\mathrm{ch})}[\ell] = \tfrac12\big(T^{(\mathrm{ch})}[\ell] + B^{(\mathrm{ch})}[\ell]\big) + \tfrac12\big(L_y[\ell] + R_y[\ell]\big)
```

Each term is a median over reference pixels of the same read:

```math
\begin{aligned}
T^{(\mathrm{ch})}[\ell] &= \operatorname{med}\big\{\,Y[\ell](x, y') : x \in X_{\mathrm{ch}},\ 2044 \le y' < 2048\,\big\}, &
B^{(\mathrm{ch})}[\ell] &= \operatorname{med}\big\{\,Y[\ell](x, y') : x \in X_{\mathrm{ch}},\ 0 \le y' < 4\,\big\},\\
L_y[\ell] &= \operatorname{med}\big\{\,Y[\ell](x, y') : 0 \le x < 4,\ |y' - y| \le 4\,\big\}, &
R_y[\ell] &= \operatorname{med}\big\{\,Y[\ell](x, y') : 2044 \le x < 2048,\ |y' - y| \le 4\,\big\}
\end{aligned}
```

- **Columns $`X_{\mathrm{ch}}`$:** the columns of the channel without the left and right reference columns. That is [4, 64) for channel 0, [1984, 2044) for channel 31, and [64·ch, 64·ch + 64) for the others.
- **Pixels per median:** T and B use 4 × 64 pixels (4 × 60 for channels 0 and 31). L and R use 4 columns × 9 lines.
- **Corrected read:** $`Y[\ell](x, y) - c_y^{(\mathrm{ch})}[\ell]`$, applied to the 64 columns of the channel, including the left and right reference columns. The top and bottom reference lines are not corrected ($`c = 0`$).

The `mode` bitmask selects the terms. `TOP` or `BOTTOM` alone replace the first term by T or B, and `TOPnBOTTOM` keeps their mean. `LEFT`, `RIGHT` and `RIGHTnLEFT` do the same for the second term. `ALL` is Eq. 3.

#### DPU: the onboard slope, per pixel

The DPU starts from the MACC layers $`Y_i[\ell]`$. Let $`n_{\mathrm{good}}`$ be the number of groups before the first group that contains a saturated read (65535 in the raw reads $`Y_i[\ell]`$). The reads are corrected for the reference pixels (`ALL`, previous paragraph) in double precision. With y(i) and ch(i) the line and channel of pixel *i*:

```math
Y'_i[\ell] = Y_i[\ell] - c_{y(i)}^{(\mathrm{ch}(i))}[\ell]
```

On a matrix that is not the full detector, $`Y'_i = Y_i`$. Then

```math
G_g = \frac{1}{n_r} \sum_{r=0}^{n_r-1} Y'_i[g\,n_r + r],
\qquad
d_g = G_{g+1} - G_g,
\qquad
g < n_{\mathrm{good}}
```

and the Kubik et al. (2016) estimator, with $`\sigma_{\mathrm{DPU}}`$ = `dpu_noise` (e⁻) and $`g_{\mathrm{DPU}}`$ = `dpu_gain` (e⁻ ADU⁻¹, the $`f_e`$ of the paper), is

```math
\alpha = \frac{1 - n_r^2}{3\,n_r\,(n_r + n_d)},
\qquad
\sigma = \frac{\sigma_{\mathrm{DPU}}}{g_{\mathrm{DPU}}},
\qquad
\beta = \frac{2\,\sigma^2\,g_{\mathrm{DPU}}}{n_r\,(1 + \alpha)},
\qquad
\gamma = \frac{2\,g_{\mathrm{DPU}}}{1 + \alpha}
```

```math
\bar F = \frac{1}{n_{\mathrm{good}} - 1} \sum_{g=0}^{n_{\mathrm{good}}-2} (d_g + \beta)^2,
\qquad
\hat F = \frac{\sqrt{1 + \gamma^2 \bar F} - 1}{\gamma} - \beta,
\qquad
Y_i =
\begin{cases}
\mathcal{Q}_{so}(\hat F) & n_{\mathrm{good}} \ge 2\\
65535 & \text{otherwise}
\end{cases}
```

$`\hat F`$ is in ADU per group interval, (n_r + n_d)·t_f. For an unsaturated linear pixel, $`\langle d_g \rangle = (n_r + n_d)\,S_{f,i}/g_i`$. The estimator also gives the χ² of the fit (not reduced, n_good − 2 degrees of freedom) and the square of the error of the slope estimate, in the square unit of the slope (ADU² per group interval²):

```math
\hat F_x = \sqrt{\bar F} - \beta,
\qquad
\chi^2 = \gamma\,\Big[(n_{\mathrm{good}} - 1)\,\hat F_x - \big(G_{n_{\mathrm{good}}-1} - G_0\big)\Big],
\qquad
\sigma_{\hat F}^2 = \frac{2\,(\hat F_x + \beta)^2}{(n_{\mathrm{good}} - 1)\,\big(1 + \gamma\,(\hat F_x + \beta)\big)}
```

The output has three planes: $`\mathcal{Q}_{so}(\hat F)`$, $`\mathcal{Q}(\chi^2)`$ and $`\mathcal{Q}_{so}(100\,\sigma_{\hat F}^2)`$. With $`n_{\mathrm{good}} = 2`$ there is no degree of freedom and the χ² is 0; with $`n_{\mathrm{good}} < 2`$ the last two planes are 0.

### 3.3 The algorithm

#### One pixel

The pixel provides three stages and holds no state. The charge *Q* is kept by the caller, the SCA readout loop, which decides when to reset, how long each step lasts and when to read.

<!-- Diagram labels use Mermaid's KaTeX syntax $$…$$ (Mermaid ≥ 10.9). -->

```mermaid
flowchart TD
    RST["① resetLevel()<br/>$$R = B + \frac{\sigma_{\mathrm{kTC}}}{g}\,\varepsilon_R$$"]
    INT["② integrate(φ, Δt)<br/>$$S = (\eta\,\phi + D)\,\Delta t$$<br/>$$I = \frac{S + \sqrt{S}\,\varepsilon_I}{g}$$"]
    Q[("charge $$Q$$ in ADU<br/>kept by the SCA")]
    RD["③ read(Q, nl)<br/>full well, non-linearity, read noise<br/>$$y = B + \frac{P(q)}{g} + \frac{\sigma_{\mathrm{RN}}}{g}\,\varepsilon_M$$"]
    Y["read value $$y_n$$"]

    RST -->|"first step only: $$Q = R$$"| Q
    INT -->|"every step: $$Q \leftarrow Q + I$$"| Q
    Q -->|"read steps only"| RD
    RD --> Y
    RD -.->|"$$Q$$ is not changed by the read"| Q
```

- **Single read (FLUX):** reset, one step of *T*, one read.
- **Ramp (MACC, UTR):** reset, then one step per frame, with a read after each frame that is not dropped. The n_d dropped frames are a single step of n_d·t_f.

#### Full array: one exposure

`nisp::SCA::expose()` dispatches on `exposureMode`. In the ramp loop, each step runs three tasks concurrently and then makes one threaded pass over memory (with IPC or crosstalk, a second one for the reads).

```mermaid
flowchart TD
    START["expose(flux, macc) or expose(flux, T)"] --> MODE{"exposureMode"}

    MODE -->|"FLUX"| FL["one charge step over $$T'$$, + persistence<br/>one read: $$M_i(Q_i) + \mathcal{R}_i/g_i$$<br/>output $$\mathcal{Q}_{so}$$, $$W \times H$$"]
    MODE -->|"CDS"| CD["charge at read 1, then at read N, + persistence<br/>two independent reads, difference<br/>output $$\mathcal{Q}_{so}$$, $$W \times H$$"]
    MODE -->|"MACC, UTR, DPU"| INIT["$$\mathrm{cumul}_i = 0$$ in ADU<br/>open the readout-noise stream: AR(1) states per channel and for the bias voltage<br/>allocate the uint16 cube $$W \times H \times n_g n_r$$"]

    INIT --> STEP

    subgraph STEP ["Read step: frame n_ℓ, Δt = t_f"]
        direction TB
        subgraph PAR ["Concurrent tasks"]
            direction LR
            A1["readout noise<br/>stream.next(), one task per channel + one for the bias voltage<br/>$$\mathcal{R}_i[n_\ell]$$"]
            A2["pixel charge<br/>DetectorSimulator::simulate<br/>$$R_i\ (\text{first step only}) + I_i(t_f)$$<br/>pixel ranges split over threads"]
            A3["persistence scheduler<br/>$$P_i[m]$$ over the past exposures $$j$$<br/>other half of the threads"]
        end
        SUM["$$\mathrm{frame}_i = \mathrm{charge}_i + P_i[m]$$"]
        ACC["_accumulate, threaded<br/>$$\mathrm{cumul}_i \leftarrow \mathrm{cumul}_i + \mathrm{frame}_i$$<br/>IPC and crosstalk, if any, on a copy of cumul<br/>$$Y_i[\ell] = \mathcal{Q}\big(M_i(\mathrm{cumul}_i) + \mathcal{R}_i/g_i\big)$$"]
        PAR --> SUM --> ACC
    end

    ACC --> LASTR{"last read<br/>of the group?"}
    LASTR -->|"no"| STEP
    LASTR -->|"yes"| LASTG{"last group?"}
    LASTG -->|"no"| DROP["dropped frames, if $$n_d \ge 1$$<br/>one charge step of $$n_d\,t_f$$, + persistence<br/>skip $$n_d$$ frames of the readout-noise stream<br/>$$\mathrm{cumul}_i \leftarrow \mathrm{cumul}_i + \mathrm{frame}_i$$, no layer written"]
    DROP --> STEP
    LASTG -->|"yes"| REG["store the charge map $$E_i = \eta_i\,\phi_i\,T$$<br/>in the persistence history"]

    REG --> DPUQ{"DPU?"}
    DPUQ -->|"no: MACC, UTR"| CUBE["uint16 cube<br/>$$W \times H \times n_g n_r$$"]
    DPUQ -->|"yes"| RPC["reference pixel correction, threaded<br/>$$c_y^{(\mathrm{ch})}[\ell]$$ of every read: medians of T, B, L, R<br/>the uint16 cube is not modified"]
    RPC --> DPU["per pixel, threaded<br/>cut at the first saturated group, on the raw reads: $$n_{\mathrm{good}}$$<br/>subtract the correction in double: $$Y' = Y - c$$<br/>co-add the groups $$G_g$$<br/>Kubik slope $$\hat F$$, $$\chi^2$$, squared error<br/>output $$\mathcal{Q}_{so}(\hat F)$$, $$\mathcal{Q}(\chi^2)$$, $$\mathcal{Q}_{so}(100 \sigma_{\hat F}^2)$$"]
    DPU --> MAP["uint16 slope, chi2 and squared error planes<br/>$$W \times H \times 3$$"]
```

Two invariants hold throughout the loop:

- `cumul` contains **charge only**: R_i, the integrated signal and dark, and the persistence. The IPC, the crosstalk, the read noise and the readout noise go to the layer that is written and are never fed back: the read noise does not accumulate along the ramp, and each read applies the IPC and the crosstalk once, to the whole accumulated signal.
- Pixels are independent within a step, so the charge pass needs no synchronisation. With IPC or crosstalk, a read needs the charge of other pixels: all the pixels are accumulated first, then read in a second pass. The other coupling between pixels is the readout-noise stream, generated per channel (plus the bias voltage series shared by all channels). The only coupling in time, apart from `cumul`, is the persistence history.

---

## 4. Configuration reference

### 4.1 Conventions

- One JSON object, read into a `PixelConfig`. Keys are case-sensitive; numbers can be written as integers or reals.
- **Units:** charges and noises in **electrons** (e⁻), levels and offsets in **ADU**, times in **seconds**, flux in **photons s⁻¹ pixel⁻¹** (before QE), sizes and positions in **pixels**. Every noise given in electrons is converted to ADU with the gain of the pixel it is added to.
- A FITS image given as input (a per-pixel map, the flux) is either the size of the matrix, used as is, or the size of the detector, cropped at `detectorOrigin`.
- A key that is absent keeps its default. The objects of `bias_voltage_noise`, `pedestal_drift`, `crosstalk` and `ipc` (kernel form) reject unknown keys, so a typo throws instead of being ignored.
- Every input is copied into the output header ([§6](#6-output-header)).

### 4.2 Matrix, flux and readout pattern (`nispSimulator`)

These keys are read by the tool, not by the library; a program using `nisp::SCA` directly passes the same quantities to the constructor and to `expose()`.

| Key | Unit | Default | Meaning |
|---|---|---|---|
| `width`, `height` | pix | 2048, 2048 | Size W × H of the simulated matrix, placed on the physical detector by `detectorOrigin` ([§4.3](#43-physical-detector-and-reference-pixels)). Smaller matrices are faster and keep the physics of the part of the detector they cover |
| `flux` | ph s⁻¹ pix⁻¹ | 0 (dark exposure) | Photons reaching each pixel, before QE. A number: uniform flux. A FITS image, in the forms of a per-pixel map ([§4.4](#44-pixel-model)): `"flux": "sky.fits"` (primary HDU), `"flux": {"file": "sky.fits", "hdu": 2}` or `{"file": "sky.fits", "hdu": "DET22.ZOD"}` (HDU index, 1 = primary, or EXTNAME), or `"flux": "sky.fits", "flux.hdu": 2`. The image is either W × H or detector-sized (then cropped at `detectorOrigin`; `FLUXSEC` in the header gives the section). An array: a sequence of exposures, each with its own flux, timing and output ([Sequences of exposures](ScenerySequencer#1-a-sequence-in-the-configuration)) |
| `macc` | — | — | Readout pattern `[ng, nr, nd]` ([§4.5](#45-readout-modes-and-timing)), or the three keys `ng`, `nr`, `nd` (not both) |
| `exposureTime` | s | — | Duration instead of a pattern: a UTR ramp of ⌊T/t_f⌋ reads in MACC, UTR and DPU; the time between reset and read in CDS and FLUX. Wins over `macc` if both are given (with a warning) |
| `benchmark` | — | 0 | Run the exposure N times and print timing statistics |

### 4.3 Physical detector and reference pixels

| Key | Unit | Default | Meaning and use |
|---|---|---|---|
| `detectorSize` | `[pix, pix]` | `[2048, 2048]` | Width and height of the physical detector. Defines the reference border, the output channels (its width must be a multiple of `colPerChannel`), the frame length of the readout noise and the size of a detector-sized map. Change it only for a detector other than the H2RG |
| `detectorOrigin` | `[pix, pix]` | `[0, 0]` | Physical column and row (0-based) of matrix pixel (0, 0). At the default, a small matrix starts in the corner of the detector and keeps its 4 reference columns on the left and 4 reference rows at the bottom. `[4, 4]` with 2040 × 2040 is the active area; ((2048 − W)/2, (2048 − H)/2) centres any matrix. The matrix must fit on the detector, otherwise `initialize()` throws |
| `colPerChannel` | pix | 64 | Columns read by one output channel. It fixes the channel boundaries and the readout order, and is used by the readout noise (frame length L = H_d · `colPerChannel` samples), the crosstalk, the reference pixel correction and the DPU. Must divide the detector width |
| `nChannels` | — | 32 | Number of output channels. `colPerChannel` × `nChannels` is the detector width, so one of the two keys is enough; `initialize()` throws if both are given and do not match. These two keys are the only ones that set the channels: they are refused inside `oneOverF` and `crosstalk` |
| `refBorder` | pix | 4 | Reference pixels on each side of the physical detector; 0 disables them. Every matrix pixel whose physical position is within the border gets QE = 0 |
| `refDark` | e⁻ s⁻¹ | value of `dark` | Dark current of the reference pixels. They have no photodiode, so their residual signal comes from the readout circuit; measure it on the border of dark ramps |

### 4.4 Pixel model

Each key is either a number, applied to every pixel, or a per-pixel FITS map, matrix-sized or detector-sized ([§4.1](#41-conventions)): `"dark": "dark.fits"` (primary HDU), `"dark": "dark.fits", "dark.hdu": 2` or `"DARK"` (index or EXTNAME), `"dark": {"file": "dark.fits", "hdu": "DARK"}`, or `"dark.file": "dark.fits"`.

| Key | Symbol | Unit | Default | NISP | Physical meaning | How it is used |
|---|---|---|---|---|---|---|
| `gain` | g | e⁻ ADU⁻¹ | 1 | 1.92 | Conversion gain of the whole chain (pixel, preamplifier, ADC): electrons per ADU | Divides every charge and every noise in electrons; > 0 |
| `dark` | D | e⁻ s⁻¹ | 0 | 0.01 | Thermally generated charge per second | Integrated with the signal, with its shot noise: S = (ηφ + D)Δt |
| `readNoise` | σ_RN | e⁻ RMS | 0 | 13 | White noise of one read (source follower, ADC) | Added independently at each read, never accumulated: a CDS carries √2·σ_RN |
| `QE` | η | 0 – 1 | 1 | 0.94 | Fraction of the incident photons converted to electrons | Multiplies the flux; the reference pixels get 0 |
| `KTCnoise` | σ_kTC | e⁻ RMS | 0 | 48 | Reset noise: random offset of the reset level | Drawn once per reset, shared by every read of the ramp, so it cancels in any difference of reads |
| `baseline` | B | ADU | 0 | 10 000 | Mean reset level (bias pedestal) of the pixel | Level of a read with no charge; constant over the ramp; ≥ 0. A constant offset per channel belongs here, as a map |
| `fullWell` | Q_FW | e⁻ | none | 130 000 | Largest charge the pixel can hold (arXiv:2507.11326) | Limits the charge at each read; > 0 for every pixel. With the NISP gain and baseline, the ADC saturates first (65 535 ADU) |
| `nlModel` | P | — | none (linear) | `polynomial` | Non-linear response of the pixel | A name, the same model for every pixel. Needs `fullWell`; must increase up to it for every pixel, otherwise `initialize()` throws. Only `polynomial` is built in ([§5](#non-linearity-model)) |
| `nl0`, `nl1`, `nl2` | c₂, c₃ for `polynomial` | (e⁻)⁻¹, (e⁻)⁻² | 0 | −6.25 × 10⁻⁷, 0 | Parameters of the model. From a ramp fit a₀ + a₁t + a₂t² + a₃t³ in ADU: c₂ = a₂/(a₁²·g), c₃ = a₃/(a₁³·g²) | `polynomial` uses `nl0` and `nl1`. −6.25 × 10⁻⁷ is a deviation of 5 % at 80 ke⁻, the limit for 95 % of the pixels (arXiv:2507.11326): a pessimistic pixel |

### 4.5 Readout modes and timing

| Key | Unit | Default | Meaning and use |
|---|---|---|---|
| `exposureMode` | — | `MACC` | `MACC`, `UTR`, `CDS`, `FLUX` or `DPU` (case-insensitive): what `expose()` computes ([§5](#5-readout-modes)) |
| `macc` | — | — | `[ng, nr, nd]`: *ng* groups (≥ 1) of *nr* reads (≥ 1), separated by *nd* dropped frames. One group interval is (nr + nd)·t_f; the first-to-last group span is (ng − 1)(nr + nd)·t_f (87.2 s for MACC(4,16,4)); reset to last read is (ng(nr + nd) − nd)·t_f (76 frames) |
| `frame_time` (or `fr_time`) | s | 1.45408 | Time to read the whole detector once, all channels in parallel. It is the time step between two reads (signal per frame S_f = (ηφ + D)·t_f) and the time base of the readout noise: one pixel-clock sample lasts t_f / (H_d · `colPerChannel`) = 11.1 µs |

### 4.6 Digitisation and DPU

| Key | Unit | Default | Meaning and use |
|---|---|---|---|
| `offset` (or `dpu_offset`) | ADU | 0 | Added before the `uint16` conversion: DN = value × `scale` + `offset`. Keeps negative values (noise on the CDS of a dark pixel) above 0. Applied in CDS, FLUX, DPU (slope and squared error × 100, not the χ²) and windowed reads, not to MACC/UTR layers |
| `scale` (or `dpu_scale`) | — | 1 | Multiplies the value before the conversion. For DPU, NISP uses SFACTOR = ng − 1, so DN − OFFSET is the signal integrated over `INTTIME` |
| `dpu_gain` | e⁻ ADU⁻¹ | mean pixel `gain` | Gain the DPU slope estimator assumes (header `GAIN`). Set it to simulate a DPU calibrated differently from the detector |
| `dpu_noise` | e⁻ | mean pixel `readNoise` | Read noise per read the estimator assumes (header `RN`). It does not include the readout noise: when the readout noise is large, set it to the total noise of one read |
| `refPixCorr` | — | `ON` | Reference pixel correction of the DPU (Eq. 3, mode `ALL`, [§3.2](#reference-pixel-correction)): `true`/`false`, or `"ON"`/`"OFF"` (any case). Applied only when the matrix is the full detector (header `REFCORR`) |

### 4.7 Readout noise

The three components are added at each read, to the reference pixels as to the active ones, and converted to ADU with the pixel gain. They are defined on the physical channels ([§3.2](#physical-detector-and-placement-of-the-matrix)), set by `colPerChannel` or `nChannels` ([§4.3](#43-physical-detector-and-reference-pixels)), so they do not depend on the size of the matrix. All are off by default.

#### `oneOverF`: 1/f noise of each output channel

```json
"oneOverF": { "amplitude": 15.0, "slope": 1.0, "nPoles": 16 }
```

The 1/f (flicker) noise of the preamplifier of each channel. The channels are independent. A number, `"oneOverF": 15.0`, sets the amplitude only.

| Key | Unit | Default | Meaning and use |
|---|---|---|---|
| `amplitude` | e⁻ RMS | 0 (off) | Total RMS of the channel noise in one read, over the whole band: what one pixel sees, in electrons. Over a ramp, the noise is continuous in time |
| `slope` | — | 1 | Spectral index *a*, PSD ∝ 1/ν^a. It decides how the variance is spread between pixel-scale noise and slow drifts (table below). 1 is pink noise |
| `nPoles` | — | 16 | Number of AR(1) processes that approximate the spectrum over its 5.1 decades. It does not change the total variance; fewer poles give a lumpier spectrum, 1 pole a single Lorentzian whatever the slope. About 3 per decade is enough: 16. The cost grows linearly with it (one random draw per pole and per pixel) |

The spectrum spans the pixel-clock frequencies from ν_min = 1/(2L) to 0.5 cycles per sample. With the NISP timing, each scale of the image corresponds to a frequency:

| ν (cycles / sample) | Frequency | Period | In the image |
|---|---|---|---|
| 0.5 | 45 kHz | 2 pixels | Pixel-to-pixel noise along a row |
| 1/64 | 1.4 kHz | 1 row (0.71 ms) | Bands one row tall, one channel wide |
| 1/6400 | 14 Hz | 100 rows (71 ms) | Broad bands |
| 1/(2L) | 0.34 Hz | 2 frames (2.9 s) | Level of the whole channel changing from read to read |

and `slope` sets the share of the variance in each scale (16 poles):

| `slope` | Pixel scale (< 1 row) | Bands (1 – 100 rows) | Slower (channel level) | Looks like |
|---|---|---|---|---|
| 0 | ≈ 98 % | ≈ 2 % | — | Extra white noise per pixel, no banding |
| 1 | ≈ 30 % | ≈ 40 % | ≈ 30 % | Pink noise: horizontal banding at all scales |
| 2 | — | ≈ 10 % | ≈ 90 % | Brownian: the level of each 64-column stripe changes between reads, little banding |

For visible horizontal banding, use `slope` 1 with an amplitude comparable to the read noise.

#### `bias_voltage_noise`: 1/f noise common to all channels

```json
"bias_voltage_noise": { "amplitude": 10.0, "slope": 1.0, "nPoles": 16 }
```

The noise of the bias voltages shared by every channel: **one** series, seen by all channels at the same pixel clock. In the image, full-width bands that repeat every 64 columns; the left and right reference columns measure them row by row. A number sets the amplitude only.

| Key | Unit | Default | Meaning and use |
|---|---|---|---|
| `amplitude` | e⁻ RMS | 0 (off) | Total RMS of the common series in one read |
| `slope` | — | 1 | Spectral index, as for `oneOverF` |
| `nPoles` | — | 16 | Number of AR(1) processes, as for `oneOverF` |

If the 1/f noise of a channel was measured as a whole, split it between the two components: `oneOverF.amplitude`² = total² − `bias_voltage_noise.amplitude`².

#### `pedestal_drift`: offset of each channel

```json
"pedestal_drift": {
    "sigma": 4.0,
    "t0":    0.0,
    "drift": {
        "poly": [0.0, 0.05],
        "exp":  [{"amplitude": 20.0, "tau": 30.0}],
        "sin":  [{"amplitude": 2.0, "period": 300.0, "phase": 0.0}]
    }
}
```

An offset shared by all the pixels of a channel, changing from read to read (preamplifier and ADC offsets drifting with temperature). The top and bottom reference rows measure it. A number, `"pedestal_drift": 4.0`, sets `sigma` only; every form replaces the whole previous pedestal configuration.

| Key | Unit | Default | Meaning and use |
|---|---|---|---|
| `sigma` | e⁻ RMS | 0 | Random part: one draw per channel and per frame (read or dropped), independent between channels. 64-column stripes whose level changes from read to read |
| `t0` | s | 0 | Time added to the time since the reset in f(t): places the exposure in a longer sequence (after a power-on, a mode change…) |
| `drift.poly` | e⁻, e⁻ s⁻¹, e⁻ s⁻² … | — | Coefficients [a₀, a₁, …] of Σ a_k t^k: a linear or polynomial drift |
| `drift.exp` | e⁻, s | — | Terms A·exp(−t/τ), `{"amplitude": A, "tau": τ}`: a thermal settling. One object or an array |
| `drift.sin` | e⁻, s, rad | — | Terms B·sin(2πt/P + φ), `{"amplitude": B, "period": P, "phase": φ}`: a periodic drift, e.g. thermal regulation. One object or an array |

The drift f(t) is the same for all channels and is evaluated at the read time of each row, t = `t0` + (n − 1 + (Y + ½)/H_d)·t_f. A drift within one frame therefore appears as a gradient along the columns. A constant f is pointless: it belongs to `baseline`.

### 4.8 Persistence

```json
"persistence": [{"f": 1e-6, "t": 80.0}, {"f": 5e-6, "t": 8.0}]
"persistence": {"traps": [{"f": 1e-6, "t": 80.0}], "dT": 0.0}
```

Charge trapped during an exposure and released during the following ones. After each exposure, the SCA stores the charge map E = η·φ·T (e⁻); every later step receives the charge released during that step by every trap species *k* ([§3.2](#persistence)).

| Key | Unit | Default | Meaning |
|---|---|---|---|
| `f` | s⁻¹ | — | Initial release rate of trap *k* per stored electron. f·τ is the fraction of the stored charge that the trap releases in total: f = 1e-6 s⁻¹ with τ = 80 s releases 8 × 10⁻⁵ of the charge |
| `t` | s | — | Decay time τ_k of the release |
| `dT` | s | 0 | Delay between the end of an exposure and the start of the next one (dithers, overheads). The SCA simulates exposures back to back; at the end of each exposure, every stored exposure ages by `dT`, so the charge released during the delay is lost, as on the detector. Must be ≥ 0 |

Keep the `PixelConfig` object alive as long as the `SCA` when it contains `persistence`: the persistence scheduler keeps a reference to it.

### 4.9 Inter-pixel capacitance and crosstalk

Both couple the pixels at each read ([§3.2](#inter-pixel-capacitance-and-crosstalk)) and are off by default. The values below are examples, not NISP measurements.

#### `ipc`: inter-pixel capacitance

```json
"ipc": [0.0, 0.01, 0.0,  0.01, 0.96, 0.01,  0.0, 0.01, 0.0]
"ipc": {"kernel": [0.01, 0.97, 0.02], "width": 3, "height": 1}
"ipc": "ipc.fits"
"ipc": {"file": "ipc.fits", "hdu": "IPC"}
"ipc": []
```

| Form | Kernel |
|---|---|
| Array | The same for every pixel, square: n × n weights |
| `kernel`, `width`, `height` | The same for every pixel, `width` columns × `height` rows |
| FITS cube (path, or `file` and `hdu`, or `"ipc"` + `"ipc.hdu"`) | One kernel per pixel |
| Empty array | No IPC |

- **Weights:** row by row, columns first, starting with the row below the pixel (y − 1). A kernel is the response to a signal in its central pixel: the first weight goes to the neighbour (x − 1, y − 1), the last one to (x + 1, y + 1). It has odd numbers of columns and rows.
- **Normalization:** the weights are divided by their sum, so a kernel can be given in any unit. A kernel given in the configuration must have weights ≥ 0.
- **FITS cube:** one layer per weight, in the same order, for a square kernel: 9 layers for 3 × 3, 25 for 5 × 5. The middle layer is the fraction that the pixel keeps. The cube is matrix-sized or detector-sized (cropped at `detectorOrigin`) and is normalized pixel by pixel. The layers of a pixel should not sum to more than 1: measurement errors can make them do so, and a warning then gives the number of pixels concerned. A pixel whose sum is not finite and positive is an error.
- **Memory:** the cube is kept in single precision, 189 MB for 2048 × 2048 pixels and 9 layers (values and FITS mask).

#### `crosstalk`: crosstalk between the output channels

```json
"crosstalk": [1, 1e-3, 5e-4]
"crosstalk": {"kernel": [1, 1e-3, 5e-4]}
"crosstalk": []
```

The index of a value is a **distance between channels**. A signal above the baseline in one channel adds a ghost to the pixels read at the same time in the others: here, 10⁻³ of it in the two neighbouring channels and 5 × 10⁻⁴ in the next two. With channels of 64 columns, the ghosts of column *x* are in columns x ± 64 and x ± 128 of the same row.

| Value | Meaning |
|---|---|
| First (distance 0) | 1: the source keeps its signal, the ghosts are added (no loss). Any other value throws |
| Next ones | Fraction of the signal added at that distance, on both sides; ≥ 0 |

- From 1 to `nChannels` values. The distances that are not given have no ghost, so `[1]` is no crosstalk.
- The channels are those of the instrument (`colPerChannel`, `nChannels`, [§4.3](#43-physical-detector-and-reference-pixels)). If they change without a new `crosstalk` key, the couplings are kept and applied to the new channels.
- There is no FITS form: the couplings are the same for every pixel.

### 4.10 A complete example

The NISP averages of `src/utilities/nispSimulator_default.json` (Euclid Collaboration 2025, arXiv:2507.11326), with the readout noise set to visible levels:

```json
{
    "width": 2048, "height": 2048,
    "detectorOrigin": [0, 0],
    "flux": 10.0,

    "exposureMode": "DPU",
    "macc": [4, 16, 4],
    "frame_time": 1.45408,
    "offset": 1024.0,
    "scale": 3.0,

    "gain": 1.92, "dark": 0.01, "readNoise": 13.0,
    "QE": 0.94, "KTCnoise": 48.0, "baseline": 10000.0,
    "fullWell": 130000.0, "nlModel": "polynomial", "nl0": -6.25e-7, "nl1": 0.0,
    "refBorder": 4, "colPerChannel": 64,

    "oneOverF":           { "amplitude": 15.0, "slope": 1.0, "nPoles": 16 },
    "bias_voltage_noise": { "amplitude": 10.0, "slope": 1.0, "nPoles": 16 },
    "pedestal_drift":     { "sigma": 8.0, "t0": 0.0, "drift": { "poly": [0.0, 0.05] } }
}
```

- The full detector: 4 reference pixels on each side, QE 0, dark current from `dark`.
- MACC(4,16,4) ramp fitted by the DPU. The output is DN = slope × 3 + 1024, so (DN − 1024) is the signal integrated over 87.2 s: 10 ph/s × 0.94 × 87.2 s / 1.92 ≈ 427 ADU.
- Full well of 130 ke⁻ and quadratic non-linearity, −5 % at 80 ke⁻ (negligible at this flux).
- Per read: 13 e⁻ read noise, then pink 1/f noise of 15 e⁻ per channel, 10 e⁻ of common bands, and channel offsets of 8 e⁻ RMS drifting by 0.05 e⁻/s.

---

## 5. Readout modes

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
| **DPU** | Ramp, then onboard slope fitting | Slope fitted on a UTR ramp | W × H × 3 cube: slope, χ², squared error of the slope (× 100) |

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
   - with IPC or crosstalk, they are applied to a copy of `cumul`, which is then the one that is read;
   - then the layer is written as `uint16( pixel.read(cumul, nl) + readout noise / gain )`, the read applying the full well, the non-linearity and the read noise. `cumul` itself stays linear and uncoupled.
3. For dropped frames: `cumul += integrate()` over nd frames; nothing is read or stored, but the readout-noise stream advances by nd frames.
4. Read noise and readout noise are added to the layer **only**, never to `cumul`. Every read has RN², and two consecutive reads differ by 2·RN², as on the real detector.

The output cube is allocated once at its final size. Each ramp step (accumulate, read, convert) is a single threaded pass over memory, or two with IPC or crosstalk, plus the convolutions.

### CDS

Correlated double sampling, computed literally:

```
f0 = read( charge at read 1 )            charge = reset level + 1 frame of signal
f1 = read( charge at read N )            charge = charge at read 1 + (N−1) frames
CDS = f1 − f0
```

The reset level and kTC noise cancel. The read noise and readout noise of the **two** reads remain (variance 2·RN²), as does the change of the pedestal drift between them. The result goes through `offset`/`scale`, so negative values (noise on dark pixels) survive when an offset is set.

### FLUX

A single read after the whole exposure: reset level + integration over *t* + one read. This is the fastest mode and is useful for quick image simulations.

### DPU: the onboard slope

The MACC ramp is simulated first. When `refPixCorr` is on (the default), `SCA::getReferencePixelCorrection(ramp)` then computes the reference pixel correction $`c_y^{(\mathrm{ch})}`$ of every read, with mode `ALL` ([§3.2](#reference-pixel-correction)). It uses the size, the border and the channels of the geometry. When the matrix is not the full detector there is no correction, and a warning is printed. The header keyword `REFCORR` records whether the correction was applied. Then, for each pixel (threaded, with no allocation per pixel):

1. **Saturation:** the reads are copied group by group, up to the first group that contains a saturated read (65535). Only the groups before it are used. With fewer than 2 good groups, the pixel is set to 65535. As on the flight DPU, 65535 is the only test: a pixel whose full well is reached below the ADC limit is not flagged, and its plateau is fitted as signal.
2. **Reference pixel correction** (when on): while the reads are copied, the pixel's correction is subtracted in double. The `uint16` ramp itself is never modified, so nothing wraps around below 0, and the saturation test sees the raw reads.
3. **Co-adding:** each group is averaged (`Ramp::coadd`).
4. **Slope:** the Kubik et al. 2016 estimator (`Ramp::getSlope`) is applied, with σ = `dpu_noise` / `dpu_gain` (in ADU). It uses only the differences between groups, so the level left by the correction does not matter.
5. **Output**, one plane each:
   - plane 0, `uint16(slope × scale + offset)`, the slope being in **ADU per group interval** ((nr+nd) frames);
   - plane 1, `uint16(χ²)`, not reduced: n_good − 2 degrees of freedom, 0 with 2 good groups;
   - plane 2, `uint16(error² × 100 × scale + offset)`, the squared error (variance) of the slope, in **(ADU per group interval)²**. The factor 100 limits the precision lost when rounding to `uint16`: with `scale` = 1, the step is 0.01 ADU².

The estimator is linear: with a non-linear response, the slope is that of the chord between the first and the last group, and the χ² grows with the curvature.

To convert to a flux: flux [e⁻/s] = slope × gain / ((nr+nd)·t_frame) = slope × (ng − 1) × gain / `INTTIME`, with slope = (DN − `OFFSET`) / `SFACTOR`. With `scale` = ng − 1, DN − `OFFSET` is directly the signal integrated over `INTTIME`.

### Non-linearity model

The response *P* ([§3.1](#31-one-pixel)) is a `nisp::nonlinearity::functor`, one for the SCA, called at each read with the parameters of the pixel (`nl0`, `nl1`, `nl2`) and its charge in electrons. It is chosen with `nlModel`, or given from C++ without changing the library:

```cpp
sca.initialize(cfg);                                  // sets fullWell, nl0, nl1, nl2

const std::vector<double> lut = readMyTable();        // measured response, captured by value
sca.setNonLinearity([lut](std::span<const float> p, double Q)
{
    return interpolate(lut, Q) * (1.0 + p[0] * Q);    // measured signal (e-)
});

sca.setNonLinearity(nisp::nonlinearity::functor{});   // empty: linear
```

- The function must tend to Q at low charge (otherwise it changes the gain) and increase up to the full well. It is checked on 32 charges for every pixel when it is set, and rejected with an exception otherwise; the previous model is then kept.
- It is called from several threads: it must not throw nor modify what it captures.
- To add a built-in model, write its function in `nispPixel.cxx`, declare it in `nispPixel.h` and name it in `nonlinearity::fromString()`.

### Reference pixel correction of a ramp

The DPU's correction is also available on its own, as two static functions of `nisp::SCA`. Use them, for example, on a MACC or UTR ramp returned by `expose()`:

```cpp
auto ramp = sca.expose(flux, nisp::macc{4, 16, 4});        // exposureMode MACC: uint16 cube

nisp::SCA::referencePixelsCorrection(ramp);                 // Eq. 3 (ALL), in place

// Or only the correction, without changing the cube:
std::valarray<double> c = nisp::SCA::getReferencePixelCorrection(ramp, nisp::SCA::TOPnBOTTOM);
const size_t nReads = c.size() / (2048 * 32);
double c_yk = c[(y * 32 + ch) * nReads + k];                // line y, channel ch, read k
```

| Function | Result |
|---|---|
| `getReferencePixelCorrection(ramp, mode, geometry)` | $`c_y^{(\mathrm{ch})}`$ of every read *k*, in the units of the cube, at `c[(y * nCh + ch) * nReads + k]`, nCh being the number of channels (the reads of one line and channel are contiguous). The top and bottom reference lines hold 0. Empty when `mode` is `NONE`, or, with a warning, when the correction does not apply: the image is not the full detector, there is no reference border, or the channels do not cover the detector or have no column outside the border |
| `referencePixelsCorrection(ramp, mode, geometry)` | Subtracts that correction from the cube, in place |

`geometry` is a `nisp::Geometry`: detector size, reference border and channels. By default it is the NISP one (2048 × 2048, 4, 32 channels of 64 columns), as in the example above. For another instrument, pass the one of the detector: `static_cast<const nisp::Geometry&>(sca.getGeometry())`.

`mode` is a `nisp::SCA::refPixCorrection` bitmask, `ALL` by default ([§3.2](#reference-pixel-correction)):

| Mode | Correction |
|---|---|
| `NONE` | None |
| `TOP`, `BOTTOM`, `TOPnBOTTOM` | T, B or ½(T + B) of the channel |
| `LEFT`, `RIGHT`, `RIGHTnLEFT` | L_y, R_y or ½(L_y + R_y) of the line |
| `ALL` | ½(T + B) + ½(L_y + R_y), Eq. 3 |

A top/bottom mode combined with a left/right one adds the two terms. `|` gives an `int`, so cast the result: `static_cast<nisp::SCA::refPixCorrection>(nisp::SCA::TOP | nisp::SCA::LEFT)`.

What the in-place correction stores depends on the type of the cube:

- **Floating point and signed integer:** v − c. For integer types the result is rounded to the nearest value (not truncated) and clamped to the type's range.
- **Unsigned, for example the `uint16` MACC ramp:** on raw reads, v − c is negative (c ≈ 2 × `baseline` with `ALL`), and an unsigned type cannot hold it. The function stores v − (c_k − c_0) instead, where c_0 is the correction of the first read. The differences between reads, and so the slopes, are those of Eq. 3, and the level stays the raw one. A single-read unsigned image is therefore left unchanged. Reads at the type's largest value (65535 for `uint16`) are saturated and are left unchanged too.

Both functions accept 2D (one read) and 3D cubes of any integer or floating-point type, and run on all hardware threads.

### Inter-pixel capacitance and crosstalk of an image

`expose()` applies both at each read. They can also be applied to any image, in double precision and in place, for example to look at their effect without noise:

```cpp
std::valarray<double> image(0.0, sca.size());           // ADU, row-major
image[1000 * sca.getWidth() + 500] = 1000.0;            // one source, baseline 0

sca.applyIPC(image, 0, 0, sca.getWidth(), sca.getHeight());
sca.applyCrosstalk(image, 0, 0, sca.getWidth(), sca.getHeight());
```

| Function | Result |
|---|---|
| `applyIPC(image, c0, r0, nc, nr)` | IPC of the window [c0, c0 + nc) × [r0, r0 + nr) of the matrix. Nothing without `ipc` |
| `applyCrosstalk(image, c0, r0, nc, nr)` | Adds the ghosts of the window. The image includes the `baseline` of the pixels, which gives no ghost. Nothing without `crosstalk` |
| `isIPCEnabled()`, `isIPCDynamic()`, `getIPCKernel()`, `getIPCKernelCube()` | The IPC that is configured: on or off, one kernel per pixel or not, the kernel, the cube of kernels |
| `isCrosstalkEnabled()`, `getCrosstalkKernel()` | The crosstalk that is configured, and its convolution kernel: one row, the couplings at ± *d* channels, 0 elsewhere and in the centre (it gives the ghosts alone) |

Both throw `std::out_of_range` when the window exceeds the matrix, and `std::invalid_argument` when the image is smaller than the window, if the effect is on.

---

## 6. Output header

Every comment starts with the unit in square brackets (FITS convention), e.g. `GAIN = 1.92 / [e-/ADU] Mean gain used by the DPU`.

**Readout and DPU**

| Keyword | Content |
|---|---|
| `MODE` | Readout mode |
| `NG`, `NR`, `ND` | MACC parameters (MACC, DPU) |
| `REFCORR` | `ON` when the DPU applied the reference pixel correction, `OFF` when `refPixCorr` is off or the matrix is not the full detector (DPU) |
| `INTTIME`, `EXPTIME` | Integration and exposure times, s. For DPU, `INTTIME` = (ng − 1)(nr + nd)·t_frame, the span from the first to the last group (87.2 s for MACC(4,16,4)) |
| `FRTIME` | Frame time, s |
| `GAIN`, `RN` | Gain (e⁻/ADU) and read noise (e⁻) used by the DPU |
| `OFFSET`, `SFACTOR` | DN = value × `SFACTOR` + `OFFSET` (ADU). Does not apply to the χ² plane of the DPU, nor to MACC/UTR layers |
| `PERSIST` | `ON` when persistence is modelled |
| `PHOT2PIX` | Library version |

**Simulation inputs** (written only for the components that are on)

| Keyword | Content |
|---|---|
| `PIXGAIN`, `PIXDARK`, `PIXRN`, `PIXQE`, `PIXKTC`, `PIXBASE` | Pixel model as given: the value (e⁻/ADU, e⁻/s, e⁻, –, e⁻, ADU), or the file name of the per-pixel map with its HDU, e.g. `'dark.fits[DARK]'`. Recorded by `DetectorSimulator`, so every simulator gets them |
| `DETSIZE`, `DETSEC` | Physical detector size, and the section of it covered by this image (`'[1:2048,1:2048]'`, 1-based) |
| `REFBORD`, `REFDARK` | Reference border (pix) and reference-pixel dark current (e⁻/s, when `refDark` is given) |
| `PIXFW`, `NLMODEL`, `PIXNL0` – `PIXNL2` | Full well (e⁻), non-linearity model and its parameters, for the keys that are given (value or map file). `NLMODEL` is `custom` after `setNonLinearity()`, `linear` if the functor is empty |
| `XTALKCOL`, `XTALKNCH`, `XTALKn` | Crosstalk: columns per output channel, number of channels, and coupling at a distance of *n* channels (those that are not 0) |
| `CHANWID` | Columns per output channel (readout noise) |
| `OOFAMP`, `OOFSLOPE`, `OOFPOLES` | `oneOverF` |
| `BIASAMP`, `BIASSLOP`, `BIASPOLE` | `bias_voltage_noise` |
| `PEDSIGMA`, `PEDT0`, `PEDPOLY`, `PEDEXPn`, `PEDSINn` | `pedestal_drift`: σ, t₀, then the drift terms (`'0, 0.05'`, `'A=20 tau=30'`, `'B=2 P=300 phi=0'`) |
| `PERDT`, `PERTRn` | Persistence: delay between two exposures (s) and traps, release rate (s⁻¹) and decay time (s) (`'f=1e-06 tau=80'`) |
| `FLUX`, `FLUXSEC` | Written by `nispSimulator`: input flux (ph/s/pixel, or file name) and the section of the detector it covers |

A file is recorded by its name only (a FITS string holds at most 68 characters). A simulator built on phot2pix adds its own inputs with `setInputCard()` and writes them with `writeInputCards()`; `PixelConfig::card(key)` gives the value or file name of any key.

---

## 7. Noise and detector effects

| Effect | Where | Model | Configured by |
|---|---|---|---|
| Shot noise (signal, dark) | `pixsim::integrate()` | Gaussian approximation of Poisson; accumulates | `QE`, `dark`, flux |
| Reset (kTC) noise, baseline | `pixsim::resetLevel()` | Drawn once per ramp, at the first frame | `KTCnoise`, `baseline` |
| Full well | `nisp::pixel::read()` | Charge limited at each read, per pixel | `fullWell` |
| Non-linearity | `nisp::pixel::read()` | Response model applied to the limited charge, before the read noise | `nlModel`, `nl0` – `nl2`, or `SCA::setNonLinearity()` |
| Inter-pixel capacitance | `nisp::SCA`, `phot2pix::PixelConvolution` | Signal shared with the neighbours at each read, by one kernel or one kernel per pixel; the charge is unchanged | `ipc` |
| Crosstalk | `nisp::SCA`, `phot2pix::PixelConvolution` | Ghosts added at each read to the same position in the other channels; limited by the full well of the source | `crosstalk`, `colPerChannel`, `nChannels` |
| Read noise | `pixsim::read()` | Independent at each read | `readNoise` |
| 1/f noise | `nisp::OneOverFNoise` | Per output channel, independent between channels, continuous in time over a ramp | `oneOverF` |
| Bias voltage noise | `nisp::OneOverFNoise` | 1/f series common to all channels, at the same pixel clock | `bias_voltage_noise` |
| Pedestal drift | `nisp::OneOverFNoise` | Offset per channel and per frame: drift equation + white noise | `pedestal_drift` |
| Reference pixels | `nisp::SCA`, `nisp::Geometry` | Physical border with QE = 0 | `refBorder`, `refDark` |
| Reference pixel correction | `nisp::SCA` | Eq. 3 of arXiv:2507.11326: medians of the top and bottom reference lines per channel, plus sliding medians of the left and right reference columns per line. Applied by the DPU, or on request to any cube | `mode` argument; the DPU uses `ALL`, switched by `refPixCorr` |
| Persistence | `nisp::persistence` | Charge released by traps filled during previous exposures | `persistence` |
| Digitisation | `nisp::SCA` | Rounding and clamping to `uint16`; `offset`/`scale` are applied in CDS, FLUX and DPU (not to its χ²), not to MACC/UTR layers | `offset`, `scale` |

### Reference pixels and the readout noise they measure

The reference pixels are read through the same channels as the active pixels, so a reference-pixel correction can only remove the readout noise the two share:

- **`pedestal_drift`**: one offset per channel and per frame. The top and bottom reference rows of a channel measure it.
- **`bias_voltage_noise`**: the same 1/f series in every channel at the same pixel clock, so horizontal bands across the full width. The left and right reference columns measure it row by row.
- **`oneOverF`**: independent in each channel. Only its slow part is seen by the reference rows; the rest stays in the data, as on the real detector.

The reference pixel correction ([§3.2](#reference-pixel-correction)) uses exactly these two measurements. T and B of each channel remove the pedestal, and $`L_y`$ and $`R_y`$ remove the bias voltage noise line by line.

### Inter-pixel capacitance, crosstalk and the noise

Both act on the signal before the read noise, which is measured with the IPC on a real detector and is therefore not smoothed again. They do act on the shot noise, whose statistics change:

- **IPC** correlates neighbouring pixels and lowers the variance of a uniform image by ΣK². With 1 % to each of the 4 neighbours, the variance is 0.92 of the Poisson one, so a photon transfer curve overestimates the gain by about 8 %.
- **Crosstalk** adds ghosts, and a ghost of the shot noise of its source, at the level of the couplings.

### Persistence

The `persistence::scheduler` stores the charge map of every exposure. During each later frame it computes the released signal from all of them, on half of the available threads, while the pixel charge is simulated on the other half, and adds it to the charge.

---

## 8. The `nispSimulator` tool

`src/utilities/nispSimulator.cxx` simulates one SCA from a JSON configuration (built when the DSTCore library is available; installed in `bin/`).

```bash
nispSimulator -o exposure.fits -d 1 my_config.json
```

| Option | Effect |
|---|---|
| `-o <file>` | Output FITS file. Without it, the simulation runs but nothing is written (a warning is printed). With a sequence (`"flux"` array): the directory the outputs are written to ([details](ScenerySequencer#5-nispsimulator-with-a-sequence)) |
| `-d <level>` | 0: warnings only; 1: run summary (matrix and its origin, mode, flux); 2: parameters, flux statistics, timings; 3+: full configuration dump |

The tool reads `width`, `height`, `flux`, `macc` or `exposureTime`, and `benchmark` ([§4.2](#42-matrix-flux-and-readout-pattern-nispsimulator)), passes the whole configuration to `SCA::initialize()`, reads a flux image with `SCA::loadImage()` (the file and HDU, as for a per-pixel map) and maps it onto the matrix with `SCA::toMatrix()` (matrix-sized, or detector-sized cropped at `detectorOrigin`), and adds `FLUX` and `FLUXSEC` to the output header.

In `DPU` mode, the three planes are written as three HDUs, `SCI` (slope), `CHI2` and `ERR2` (squared error of the slope, × 100). The other modes write one image or cube.

Ready-made configurations:

- `src/utilities/nispSimulator_default.json`: averaged NISP pixel values, with full well and non-linearity; readout noise present but set to 0.
- `src/utilities/testcases/`: one configuration per readout mode, and a sequence.
  - `CDS`: every detector effect off (gain 1, QE 1, non-linearity 0), 10 ph/s.
  - `MACC`: the same with the non-linearity on.
  - `UTR` and `DPU`: NISP pixel values, non-linearity and readout noise on, 10 ph/s.
  - `FLUX` and `SEQUENCE`: the same, with the flux read from FITS images (paths to adapt).

---

## 9. Validation

`src/test/gTest_nispReadoutModes.cxx` checks every readout mode on the 2040 × 2040 active area (origin [4, 4], no reference pixel), with MACC(4,16,4), 10 ph/s, gain 1, QE 1, and every other effect off. Each frame then adds S = 10 × 1.45408 = 14.5408 ADU, with Poisson variance S. Means and variances are checked over the 4.2 million pixels within 5σ.

| Mode | Expected | Flux recovered from |
|---|---|---|
| MACC | Layer 16g+f = (20g+f+1)·S | (G₃ − G₀) / (3·20·t_frame) |
| UTR | Read k = (k+1)·S, 76 reads | Least-squares slope |
| CDS | 75·S = 1090.56 ADU | CDS / (75·t_frame) |
| FLUX | 76·S = 1105.10 ADU | Value / (76·t_frame) |
| DPU | Slope ≈ 290.8 ADU per group, scatter ≈ 9.4 ADU; squared error ≈ 71 ADU²; loose check of the χ² | Slope / (20·t_frame) |

`src/test/gTest_nispSCA.cxx` checks the noise budget and the detector:

- Every MACC read carries RN² and read-to-read differences carry 2·RN² (read noise does not accumulate).
- A CDS carries 2·RN² while the reset level and kTC noise cancel.
- 1/f noise produces column-correlated banding.
- The readout noise is converted to ADU with the pixel gain; the bias voltage noise is identical in every channel; a pedestal drift survives in CDS.
- The reference pixels follow the physical border (default origin in the corner, 2044-wide matrix, `refBorder` 0); a detector-sized map is cropped at the origin; a matrix larger than the detector is rejected.
- The output header records the inputs.
- Persistence raises the signal of the next exposure, by the charge released during its frame after the delay `dT`.
- Full well and non-linearity, with a 20 ke⁻ well reached below the ADC limit:
  - the mean and variance of every MACC read follow P(q), then stay at the full well with the read noise only; FLUX, CDS and a window follow the same response;
  - saturated reads give exactly the full well of each pixel, read from a map;
  - the model can be a capturing functor; invalid full wells and models are rejected; the header records them;
  - the DPU slope is that of the chord of the response, and a plateau below 65535 is not flagged.

`src/test/gTest_persistence.cxx` checks the persistence model on its own: the release over a step, consecutive steps adding up to one step, a total release of f·τ, the scheduler clock and the delay between exposures.

`src/test/gTest_nispOneOverFNoise.cxx` checks the readout-noise components on their own: configuration parsing and errors, the placement of the matrix on the physical detector (partial channels line up with the physical ones), the spectral index set by `slope`, bias voltage noise shared by all channels (including the pixel-clock offset of a 2040-wide frame), one pedestal offset per channel and per frame with RMS σ, and the drift equation evaluated at the read time of each physical row, in windows and along a stream.

Three suites check the coupling of the pixels. The expected images come from a separate calculation, which spreads the signal of each pixel instead of collecting it:

- `src/test/gTest_pixelConvolution.cxx`: `phot2pix::PixelConvolution` alone, with cubes made in memory. Kernel forms, normalization, warning and errors; one kernel or one per pixel; borders, windows smaller than the kernel, number of threads; kernels that are not normalized.
- `src/test/gTest_nispIPC.cxx`: the `ipc` key (arrays, FITS cubes with sums below and above 1, matrix and detector sizes) and the IPC of a 2048 × 2048 image holding one lit pixel out of 64 along each axis, full frame and windows, with one kernel and with one per pixel.
- `src/test/gTest_nispCrosstalk.cxx`: the `crosstalk` key and the channels, the ghosts of sources placed so that they never overlap, and of a large Gaussian source, for several channel widths; the baseline and the full well of the source; the ghosts in MACC, FLUX, CDS and windowed exposures.

Set the environment variable `NISP_IPC_OUTPUT` to a directory to save the input and output images of the last two suites:

```bash
NISP_IPC_OUTPUT=images ./gTest_nispCrosstalk --gtest_filter='SCACrosstalk.*'
```

`gTest_nispSCA.cxx` also checks that `colPerChannel` and `nChannels` are read at the top level of the configuration only.

`src/test/gTest_detector.cxx` checks `phot2pix::Geometry` (placement, extraction of a matrix-sized or detector-sized map, derived geometries); `gTest_pixelconfig.cxx` checks `PixelConfig::card()` for every input form; `gTest_nispPixel.cxx` checks the pixel stages, the full well and the non-linearity of one read, and the `Ramp` estimator.

---

## 10. Known limitations

- With MACC parameters, the UTR mode produces the same output as MACC, and the MACC output always writes `MODE = MACC`. The MACC `EXPTIME` counts the drops after the last group.
- DPU:
  - the planes hold integers, so the χ² (about n_good − 2 on average for a linear pixel) and the squared error × 100 × `scale` are rounded to 1;
  - no test runs the DPU at a gain other than 1;
  - `nispSimulator` writes the three HDUs only when `exposureMode` is written exactly `DPU`, and for a single exposure; otherwise the output is the three-plane cube.
- The readout-noise model and the crosstalk read every channel in the same direction, and the readout noise has no row or frame overhead.
- Inter-pixel capacitance and crosstalk:
  - the IPC is applied before the full well. A saturated pixel gives its neighbours a part of all the charge it collected, and in a ramp they keep rising after it saturates. The crosstalk is limited by the full well;
  - a window only contains the signal of its own pixels: its border misses the IPC of the pixels around it, and it has no ghost of the sources outside of it (none at all if it is not wider than a channel);
  - the IPC is not recorded in the output header;
  - the couplings of the crosstalk are ≥ 0 and depend only on the distance between the channels: no negative ghost, and no difference between the two sides or between pairs of channels;
  - the crosstalk only gives ghosts of the first order: a ghost gives none;
  - each read is convolved, which costs time: the more weights in the kernel, the longer. A short crosstalk kernel is faster than 32 values.
- Persistence:
  - the charge stored after a MACC ramp uses the last integration interval instead of the total exposure;
  - the released charge, in electrons, is added to the charge in ADU without division by the gain;
  - the stored charge is not limited by the full well: a saturated star leaves more persistence than its pixels can hold.
- Full well and non-linearity:
  - the ADC limit (65535) is applied to each read in MACC, UTR and DPU, but only to the result in CDS, FLUX and windowed reads: a CDS can subtract a read above the ADC limit, and with `scale` ≠ 1 the limit moves;
  - `pixelresponse()` of a lone `nisp::pixel` is linear and has no full well: it calls `read(charge)` without the model. The SCA never uses it;
  - there is no count-rate non-linearity, no brighter-fatter effect and no blooming.
- Reference pixel correction:
  - it needs the full detector. A smaller matrix, such as the 2040 × 2040 active area used by `gTest_nispReadoutModes`, is not corrected, and a warning is printed;
  - in the DPU, `refPixCorr` only switches the correction on or off. The mode is always `ALL`;
  - no test covers the correction yet, nor the DPU with another number of channels than 32.

---

*The Phot2Pix logo has been generated by OpenAI's ChatGPT Image Generation tool (ImageGen).*

*This wiki page has been generated by AI (Claude Opus 5.5) by reading the code itself.*

<!-- This page is generated from wiki/NISP-example.md in the phot2pix repository. Edit it there: changes made in the web interface are overwritten on the next push. -->
