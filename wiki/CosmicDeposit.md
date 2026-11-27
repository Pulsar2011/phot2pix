# CosmicDeposit: cosmic rays and radioactivity

`phot2pix::CosmicRay::CosmicDeposit` deposits the charge of particle hits on the detector during each step of an exposure: cosmic rays in space (galactic and solar particles), muons on the ground, and particles from radioactive materials near the sensor.

It evaluates the impact of these particles **on the images**; it is not built to be an accurate representation of the physics of the interaction of cosmic rays with the detector material. The model is the simplest one that gets right what an image sees: the charge of a hit, the length of its track, where the track ends (Bragg peak), and the pixels it touches. Physics that does not change these quantities is left out and listed in the [known limitations](#9-known-limitations).

The user provides the spectrum of each species **at the detector**, that is, behind the shielding of the spacecraft or of the instrument, as computed with SPENVIS, CREME96 or Geant4. The model does not transport the particles through the shielding.

---

## 1. The model

The hits are drawn for one step of an exposure, on the whole physical detector.

### 1.1 Rate

The number of hits in a time step $\Delta t$ is drawn from a Poisson distribution of mean

```math
\lambda = \Delta t\; A \sum_s (1 - \tau_s)\, G_s \int J_s(E)\,\mathrm{d}E
```

| Symbol | Unit | Meaning |
|---|---|---|
| $\Delta t$ | s | Duration of the step: the `dt` of `generate()` |
| $A$ | m² | Area of the **physical** detector, whatever the size of the simulated matrix: width × height in pixels (`detectorSize`, by default the matrix) × `pitch_um`² |
| $s$ | — | A species: one entry of `"spectrum"` ([§2.2](#22-a-spectrum-entry)), that is a particle, its spectrum and its angular model |
| $E$ | MeV/nucleon or MeV | Kinetic energy: **per nucleon for a nucleus**, **per particle** for any other particle (e⁻, e⁺, muons) |
| $J_s(E)$ | m⁻² s⁻¹ sr⁻¹ per unit of $E$ | Differential intensity of the species at the detector, in the reference direction of its angular model: the table of its file ([§3.1](#31-columns-and-units)) |
| $G_s$ | sr | Geometric factor of its angular model (below) |
| $\tau_s$ | — | Share of its hits whose track is too long to be drawn ([grazing tracks](#12-directionality)) |

#### Geometric factor

In the direction $\Omega$, the intensity of the species is $I_s(E, \Omega) = J_s(E)\, f_s(\Omega)$, where $f_s$ is the shape of its angular model. A flat sensor receives $I_s\,\mu$ particles per unit of area, time, solid angle and energy, $\mu$ being the cosine of the angle between $\Omega$ and the sensor normal. The geometric factor is the integral of $f_s\,\mu$ over the directions that reach the faces of the species (`front`, `back` or `both`):

```math
G_s = \int f_s(\Omega)\,\mu\,\mathrm{d}\Omega
```

For an isotropic intensity ($f_s = 1$), $G_s = \pi$ per face. The other models are in the table of [§1.2](#12-directionality).

#### Hits not drawn

Particles that cross the sensor almost edge-on leave very long tracks, and those longer than `maxTrackPixels` are skipped: this is why the formula has the factor $(1 - \tau_s)$. The simulated images therefore contain slightly fewer hits than really cross the sensor. Two accessors give the true numbers, with the long tracks included ([§7](#7-accessors-without-simulating)):
- `eventRate(area)` is the number of hits per second;
- `expectedEvents(area, dt)` is the mean number of hits in a step of `dt` seconds.

The share of hits that are skipped, all species together, is given by `truncatedFraction()` and written in the output header as `CRTRUNC`.

#### Spectrum between table points

The spectrum of a species is obtained from a table of points $(E_k, J_k)$, read from a file ([§3](#3-spectrum-and-angular-files)) named in the JSON configuration of the simulation, so the spectrum is binned in energy. Inside each energy bin $[E_k, E_{k+1}]$, $J_s$ is the power law through both points:

```math
J_s(E) = J_k \left(\frac{E}{E_k}\right)^{-g_k},\qquad g_k = -\frac{\ln(J_{k+1}/J_k)}{\ln(E_{k+1}/E_k)}
```

Because $J_s$ is a power law inside each bin, the code needs no numerical integration and no tabulated random draw and uses two exact formulas, with $r_k = E_{k+1}/E_k$.

- **The flux** $F_k$ in the bin $k$ is the number of particles per square metre, per second and per steradian whose energy is contained in the bin. It sets how often the bin is chosen, as explained at the end of this section:

  ```math
  F_k = \int_{E_k}^{E_{k+1}} J_s\,\mathrm{d}E = J_k\,E_k\,\frac{r_k^{\,1-g_k} - 1}{1 - g_k}\qquad\left(J_k\,E_k \ln r_k \text{ if } g_k = 1\right)
  ```

- **The energy** $E$ of a hit in the bin comes from a random number $u$, drawn uniformly in $[0, 1)$. The link between them is that $u$ is the fraction of the particles of the bin whose energy is below $E$:

  ```math
  u = \frac{1}{F_k}\int_{E_k}^{E} J_s\,\mathrm{d}E'
  ```

  So $u = 0$ gives $E_k$, $u = 0.5$ gives the median energy of the bin, and $u$ close to 1 gives an energy close to $E_{k+1}$. Solving this equation for $E$ gives:

  ```math
  E = E_k\left[1 + u\left(r_k^{\,1-g_k} - 1\right)\right]^{1/(1-g_k)}\qquad\left(E_k\, r_k^{\,u} \text{ if } g_k = 1\right)
  ```

  Although $u$ is drawn from a uniform distribution, $E$ **follows the power law of the bin**. As a consequence, $E$ is drawn more often near $E_k$ when the spectrum falls ($g_k > 0$). $E$ is uniformly distributed in $\ln E$ only when the power-law index is $g_k = 1$, and uniformly distributed in $E$ only for a flat spectrum, $g_k = 0$.

The energy of a hit is therefore not one of the tabulated values and can be anywhere in its bin. $J_s$ is taken as $0$ outside the limits of the table.

Each hit is drawn in four steps. The code first picks a species and one of its bins, with a probability proportional to $(1 - \tau_s)\,G_s\,F_k$. It then draws an energy in that bin with the formula above, a direction as described in [§1.2](#12-directionality), and an entry point uniformly distributed over the physical detector.

### 1.2 Directionality

Particles do not all arrive from the same direction, and their direction matters twice. It sets how many of them cross the sensor, through the geometric factor $G_s$ of [§1.1](#11-rate), and it sets the length and the orientation of the track each one leaves. Each species therefore has an angular model, chosen with the key `"angular"` of its spectrum entry ([§2.2](#22-a-spectrum-entry)). The model says how the intensity $I_s(E, \Omega) = J_s(E)\,f_s(\Omega)$ changes with the direction $\Omega$. In the table below, $\theta$ is the angle between the direction and the sensor normal, $\mu = \cos\theta$, and $\theta_z$ is the angle between the direction and the zenith.

| `"model"` | Typical use | Intensity in a direction | How $\mu$ is drawn | $G$ per face |
|---|---|---|---|---|
| `isotropic` (default) | Cosmic rays in space | $J$, the same from every direction | $\mu^2$ uniform | $\pi$ |
| `surface` | Radioactive emitter on a surface parallel to the sensor ([§3.3](#33-from-a-surface-activity-to-j-surface-model)) | $J/\mu$, so that the same number of particles crosses the sensor in every solid angle | $\mu$ uniform | $2\pi$ |
| `cosn` | Muons on the ground ($n \approx 2$) | $J\cos^n\theta_z$ from above, 0 from below | Horizontal sensor: $\mu^{n+2}$ uniform. Otherwise by rejection | Horizontal sensor: $2\pi/(n+2)$, front face only. Otherwise numerical |
| `table` | A measured zenith distribution | $J\,f(\theta_z)$ from above, 0 from below, with $f$ read from a table ([§3.1](#31-columns-and-units)) | By rejection | Numerical |

The columns read as follows.

- The intensity column gives $I_s$ in one direction, $J$ being the intensity of the spectrum file ([§1.1](#11-rate)).
- The $\mu$ column says how the code draws the angle of each hit to the sensor normal. The number of particles that cross the sensor from a direction is $I_s\,\mu\,\mathrm{d}\Omega$ and not $I_s\,\mathrm{d}\Omega$, because a flat sensor seen at a grazing angle presents a smaller area. For an isotropic intensity, $\mu$ is therefore distributed as $2\mu\,\mathrm{d}\mu$, and the code draws it as $\mu = \sqrt{u}$ with $u$ uniform in $[0, 1)$: this is what "$\mu^2$ uniform" means. The other models follow from the same reasoning. "By rejection" means that the code draws directions uniformly over the sky and keeps each one with a probability proportional to $f_s\,\mu$. In the closed-form cases, the azimuth $\varphi$ of the track on the matrix is uniform in $[0, 2\pi)$.
- The last column is the geometric factor $G_s$ of one face. A species that reaches both faces (`"faces": "both"`) gets twice this value. "Numerical" means that $G_s$ is integrated numerically, to about 10⁻³ relative precision.

- The `isotropic` and `surface` models are defined in the frame of the sensor, so they do not depend on how the sensor is oriented. The `cosn` and `table` models are defined in the sky, relative to the zenith, so their directions and their $G$ depend on `tilt_deg`, the angle between the sensor normal and the zenith. The normal tilts towards $+x$ of the matrix. On a tilted sensor, and always with `table`, $G$ is integrated numerically and the directions are drawn by rejection.
- Each species reaches the `front` face of the sensor, its `back` face, or `both` (key `faces`). The face that a hit crosses is recorded with the hit, in the column `FACE`.
- A particle that crosses the sensitive depth $d$ at the angle $\theta$ travels $L = d/\mu$ inside it. Seen from above, its track on the matrix is $(d/\text{pitch})\tan\theta$ pixels long and points in the direction $\varphi$.
- Tracks longer than `maxTrackPixels` are not drawn. A track is that long when $\mu < \mu_{\min} = 1/\sqrt{1+q^2}$, with $q = \text{maxTrackPixels}\cdot\text{pitch}/d$. For an isotropic intensity, the share of hits that are skipped is $1/(1+q^2)$. `truncatedFraction()` gives this share for the configured species (see *Hits not drawn* in [§1.1](#11-rate)).

### 1.3 Energy deposited

A hit deposits an energy $\Delta E$ in the sensitive layer; its charge is $\Delta E$ / `eVPerElectron` ([§1.4](#14-charge-along-the-track)). $\Delta E$ depends on how the path through the layer, $L = d/\mu$ ([§1.2](#12-directionality)), compares with how far the particle can go and how much it loses on the way. Three quantities of the particle are used:

- $T$ is its kinetic energy at entry (MeV), for the whole particle; for a nucleus, it is $A$ times the energy per nucleon.
- $R(T)$ is its **CSDA range**, the path over which it would lose all of $T$ if it lost energy continuously at the mean rate (continuous slowing-down approximation, CSDA).
- $\langle\Delta E\rangle = T - T_{\text{out}}$ is its **CSDA mean loss** over $L$, $T_{\text{out}}$ being the energy left after the path $L$.

$R$ and $\langle\Delta E\rangle$ come from the range table of the species, built from its [mean stopping power](#mean-stopping-power). Each hit is then in one of three cases, tested in this order:

| Case | Condition | Path in the layer | $\Delta E$ |
|---|---|---|---|
| **Stops** | $R \le L$ | $R$ (the projected track is shortened by $R/L$) | $T$: all its kinetic energy |
| **Thick crossing** | $R > L$ and $\langle\Delta E\rangle > 0.1\,T$ | $L$ | [Gaussian](#thick-crossing-gaussian) of mean $\langle\Delta E\rangle$, kept in $[0, T]$ |
| **Thin crossing** | $R > L$ and $\langle\Delta E\rangle \le 0.1\,T$ | $L$ | [Landau](#thin-crossing-landau), kept in $[0, T]$ |

The three cases exist because the energy loss behaves differently in each of them:

- A particle that stops in the layer leaves all its energy there, and its charge rises to a Bragg peak at the end of its track ([§1.4](#14-charge-along-the-track)).
- In a thin crossing, the particle loses only a small share of its energy, and that loss fluctuates strongly from one hit to the next. Most hits lose a little less than the mean, and a few lose much more, when one collision gives an electron a large energy. This is the Landau distribution, whose most probable value $\Delta_p$ is below the mean loss.
- In a thick crossing, the particle loses more than 10 % of its energy (`kThickShare` in the code) and slows down along its path. Its loss is then the sum of many collisions, which is close to a Gaussian around the mean loss.

#### Thin crossing: Landau

```math
\Delta E = \Delta_p + \xi\,(\lambda - \lambda_{\text{mode}})
```

$\lambda$ is drawn from the standard Landau distribution (`DST::Math::standardLandau`) and $\lambda_{\text{mode}}$ is its mode, so the most frequent deposit is $\Delta_p$. The width $\xi$ and the most probable loss $\Delta_p$ are those of the PDG, at the velocity of the particle at entry:

```math
\xi = \frac{K}{2}\,\langle Z/A\rangle\, z^2\,\frac{x}{\beta^2},\qquad
\Delta_p = \xi\left[\ln\frac{2m_ec^2\beta^2\gamma^2}{I} + \ln\frac{\xi}{I} + j - \beta^2 - \delta(\beta\gamma)\right]
```

| Symbol | Meaning |
|---|---|
| $K$ | $4\pi N_A r_e^2 m_e c^2$ = 0.307075 MeV cm² mol⁻¹ |
| $\langle Z/A\rangle$, $I$ | Ratio of atomic number to mass number, and mean excitation energy, of the material |
| $z$ | Charge number of the particle |
| $x$ | Path $L$ in g/cm² ($L$ times the density) |
| $\beta$, $\gamma$ | Velocity (in units of $c$) and Lorentz factor of the particle at entry |
| $m_ec^2$ | Electron rest energy |
| $\delta(\beta\gamma)$ | Density effect (Sternheimer) |
| $j$ | 0.200 (Bichsel) |

A negative $\Delta_p$ is set to 0. A detector family can replace this distribution by overriding `energyDeposit_eV()` ([§6](#6-using-it-from-a-detector-family)). In the current implementation of Cosmic Ray hit, the Landau distribution is used even where the Vavilov distribution would apply ([§9](#9-known-limitations)).

#### Thick crossing: Gaussian

$\Delta E$ is drawn from a Gaussian of mean $\langle\Delta E\rangle$ with Bohr variance

```math
\sigma^2 = \xi\,W_{\max}\left(1 - \frac{\beta^2}{2}\right)
```

with $\xi$ defined as above over the path $L$, and $W_{\max}$ the largest energy the particle can give one electron in a collision. Because the particle slows down noticeably along $L$, $\beta$, $\xi$ and $W_{\max}$ are taken at the mid-path energy $(T + T_{\text{out}})/2$, not at entry.

#### Mean stopping power

The range $R$ and the mean loss $\langle\Delta E\rangle$ are integrals of the mean stopping power of the species in the material. For **heavy particles** (nuclei and muons), the stopping power is given by the Bethe-Bloch formula, with $W_{\max}$ computed from the mass of the particle and with the Sternheimer density effect $\delta$. For **electrons and positrons** (rest energy below 1 MeV), it is given by the ICRU 37 collision formula; bremsstrahlung is left out, because the emitted photons leave the layer. At **low energies**, below the energy at which the formula reaches its maximum, the stopping power is taken proportional to the velocity.

### 1.4 Charge along the track

Along the path, each pixel crossed receives the CSDA energy lost on its part of the path, $E(s_i) - E(s_{i+1})$, where $s_i$ and $s_{i+1}$ are the positions along the path at which the track enters and leaves the pixel. The charge is therefore uniform along a fast track, and it rises to a Bragg peak at the end of the track of a particle that stops in the layer.

The lateral spread of the species is then applied. A pixel at the offset $(a, b)$ from a pixel of the track receives the fraction

```math
(1 - f_h)\,g_{\sigma_c}(a)\,g_{\sigma_c}(b) + f_h\,g_{\sigma_h}(a)\,g_{\sigma_h}(b)
```

of the charge of that track pixel, where $g_\sigma$ is a Gaussian of standard deviation $\sigma$ integrated over the pixel and cut at $3\sigma$. Each Gaussian is renormalised after the cut, so the fractions around a pixel of the track add up to 1 and this operation conserves the charge. The only charge lost is the part that falls outside the physical detector. The width $\sigma_c$ of the core is the `sigmaDiffusion_px` of the spectrum entry or, when the entry gives none, that of `"cosmic"`. The width $\sigma_h$ and the share $f_h$ of the optional halo are given by the `"halo"` key of the entry; without this key, there is no halo ($f_h = 0$).

The halo reproduces the broad, round charge of a snowball, or of the δ-rays around a heavy ion. It does not model their physics, and it draws no δ-ray track.

The charge is the deposited energy divided by `eVPerElectron`. A hit can land on any pixel of the physical detector, reference pixels included.

### 1.5 Delta rays (`"deltaRays": true`)

When a heavy particle crosses the sensor, some of the electrons it knocks out of the atoms, the δ-rays, receive enough energy to reach the neighbouring pixels. With `"deltaRays": true` in a spectrum entry, the model puts a share $h$ of the charge of each hit into a round halo around each pixel of its track. The option is for heavy particles only (nuclei, muons), not for e⁻ or e⁺. Everything is computed from the velocity of the particle at entry ($\beta$, $\gamma$), with closed formulas only, in the four steps below.

#### Which δ-rays leave the pixel

One collision can give an electron at most the energy $W_{\max}$, the same value as in the Bethe formula. An electron of energy $T_0$ has a CSDA range of half a pixel, so a δ-ray below $T_0$ stays in the pixel of the track, and its charge is part of the track. Only the δ-rays between $T_0$ and $W_{\max}$ feed the halo. When $W_{\max} \le T_0$, there are none: slow particles, such as a 5 MeV α or a 10 MeV proton, have no halo.

#### Share of the energy loss carried by these δ-rays

The number of δ-rays per unit of their energy $T$ is $\mathrm{d}N/\mathrm{d}T \propto (1 - \beta^2 T/W_{\max})/T^2$. The energy they carry between $T_0$ and $W_{\max}$, divided by the Bethe mean loss of the particle, gives the share $f_\delta$ of the energy loss that goes into the halo:

```math
f_\delta = \frac{\ln(W_{\max}/T_0) - \beta^2\,(1 - T_0/W_{\max})}{\ln\!\left(2 m_e c^2 \beta^2\gamma^2 W_{\max}/I^2\right) - 2\beta^2 - \delta}
```

The numerator is the energy of these δ-rays and the denominator is the Bethe mean loss, both in units of $(K/2)\,z^2\,\langle Z/A\rangle/\beta^2$ per g/cm². $m_e c^2$, $I$ and $\delta$ are defined as in [§1.3](#thin-crossing-landau). The charge of the particle, $z^2$, cancels in the ratio, so $f_\delta$ depends only on the velocity. A proton and an iron nucleus at the same energy per nucleon therefore get the same share. The halo of the iron nucleus is brighter only because its total charge is $z^2$ times larger.

#### Where the δ-rays put their energy

Around a track, the dose of the δ-rays falls as $1/r^2$ with the distance $r$ (Butts & Katz). Each interval of $\ln r$ therefore receives the same energy, from $r_0$ = ½ pixel up to $r_{\max}$, the CSDA range of an electron of energy $W_{\max}$. The sensitive layer is only $d$ thick, though, so a δ-ray that travels farther than $d$ leaves the layer and deposits only about $d/r$ of its energy in it. The weight of each interval of $\ln r$ is therefore $w(r) = \min(1, d/r)$, with $r$ and $d$ in pixels.

#### Share of the charge in the halo

Because part of their energy escapes the layer, the δ-rays leave in the layer only the share

```math
f_{\text{kept}} = f_\delta\,\frac{\int_{r_0}^{r_{\max}} w(r)\,\mathrm{d}\ln r}{\ln(r_{\max}/r_0)}
```

of the energy loss, while the track keeps $1 - f_\delta$. Of the energy left in the layer, the halo therefore holds the share

```math
h = \frac{f_{\text{kept}}}{1 - f_\delta + f_{\text{kept}}}
```

which `deltaRayShare(i, E)` returns. The charge of the hit computed in [§1.3](#13-energy-deposited) does not change: it is split into $h$ for the halo and $1 - h$ for the track, which then gets its own lateral spread ([§1.4](#14-charge-along-the-track)).

#### Size of the halo

The halo around one pixel of the track holds $h\,N_e/L$ electrons, where $N_e$ is the charge of the hit and $L$ the length of the projected track in pixels (at least 1). The halo stops at the radius where it would put less than `deltaRayCut_e` electrons on a pixel, and never goes beyond $r_{\max}$. Its shares are renormalised within that radius, so the halo conserves its charge. If the halo is below `deltaRayCut_e` everywhere, the hit gets no halo and its track keeps all the charge. The halo is sampled on rings at most a quarter of a pixel wide, each in a multiple of 4 directions, so it keeps the symmetries of the pixel grid.

For 10 µm of HgCdTe, 18 µm pixels and `deltaRayCut_e` = 1, $h$ is a few per cent for any fast particle. The halo of a minimum-ionising proton then reaches about 3 pixels from the track, and that of an iron nucleus about 25.

This is a mean profile: the model draws no δ-ray track, and $h$ does not fluctuate from one hit to another.

---

## 2. Configuration reference

The source reads the `"cosmic"` object of the detector configuration (`configure()`). Absent keys keep their value. Everything is checked on a copy before anything is kept: a bad configuration leaves the source unchanged.

### 2.1 The `"cosmic"` keys

| Key | Unit | Default | Meaning |
|---|---|---|---|
| `depth_um` | µm | — | Sensitive depth $d$ |
| `pitch_um` | µm | — | Pixel pitch |
| `eVPerElectron` | eV | — | Mean energy to create one electron-hole pair in the sensor |
| `tilt_deg` | deg | 0 | Angle of the sensor normal (front face) to the zenith, in [0, 180]; used by the `cosn` and `table` models only |
| `sigmaDiffusion_px` | pixel | 0 | Lateral spread of the charge (σ); 0: none. Default of the entries that give none |
| `maxTrackPixels` | pixel | 64 | Longest projected track drawn |
| `deltaRayCut_e` | e⁻ | 1 | Smallest charge a δ-ray halo puts on a pixel, > 0 ([§1.5](#15-delta-rays-deltarays-true)) |
| `seed` | — | 0 | Seed of the random streams; 0: not reproducible |
| `spectrum` | — | — | One entry per species ([§2.2](#22-a-spectrum-entry)) |

`acceptance_sr` is no longer read: giving it is an error. The geometric factor comes from the angular model of each entry.

The **sensor material** is not a key: it is `CosmicDeposit::Params::material`, silicon by default, set in code by a detector family ([§4](#4-sensor-material), [§6](#6-using-it-from-a-detector-family)). NISP (`nisp::SCA`) uses `nisp::CosmicDeposit`, whose material is `HgCdTe(0.3)`.

### 2.2 A spectrum entry

| Key | Default | Meaning |
|---|---|---|
| `species` | — | Name of the species |
| `Z` | — | Charge number (−1 for an electron) |
| `A` | — | Mass number of a **nucleus**: energies per nucleon |
| `mc2` | — | Rest energy (MeV) of **any other particle**, **instead of** `A`: energies per particle |
| `file` | — | CSV file (extension `.csv`, any case), otherwise FITS file |
| `hdu` | — | FITS only, **required**: index of the table HDU (1 = primary, so 2 = first extension) or its EXTNAME |
| `filterColumn`, `filterValue` | every row | When the table holds several species, the rows of this one are those where `filterColumn == filterValue` (a string or a number). Give both or neither |
| `colEnergy` | — | **Required**: name of the energy column |
| `colFlux` | — | **Required**: name of the flux column |
| `angular` | isotropic | Angular model: `{"model": "isotropic"}`, `{"model": "surface"}`, `{"model": "cosn", "n": 2}`, or `{"model": "table", "file", "hdu", "colTheta", "colIntensity"[, "filterColumn", "filterValue"]}` |
| `faces` | `front` | Faces of the sensor reached: `front`, `back` or `both` |
| `sigmaDiffusion_px` | `sigmaDiffusion_px` of `"cosmic"` | σ of the core of the lateral spread of this species (pixel), ≥ 0 ([§1.4](#14-charge-along-the-track)) |
| `halo` | none | `{"sigma_px": σ > 0, "fraction": f in [0, 1]}`: a wider Gaussian that holds the share $f_h$ of the charge of this species ([§1.4](#14-charge-along-the-track)). Both keys are required |
| `deltaRays` | `false` | `true`: δ-ray halo of each hit, from its velocity ([§1.5](#15-delta-rays-deltarays-true)). Heavy particles only: an error for e⁻ or e⁺ |

Give `A` **or** `mc2`, not both. A proton given by `"mc2": 938.272` is valid, but then its energies are per particle.

To give two populations of the same particle different footprints, write two entries. For example, ordinary hits and snowballs each get their own table (or filter) and their own `sigmaDiffusion_px` and `halo`. Each entry is a species of its own: its own rate, its own `SPECIES` in the table of hits.

The halo cost grows with its width: each pixel of the track gives $(6\sigma_h + 1)^2$ deposits, about 360 for $\sigma_h = 3$.

### 2.3 Examples

#### Space: solar protons and electrons behind the shielding

```json
"cosmic": {
    "depth_um": 10.0, "pitch_um": 18.0, "eVPerElectron": 3.6, "seed": 7,
    "spectrum": [
        {"species": "p",  "Z": 1,  "A": 1,           "file": "spe_shielded.fits", "hdu": "SPECTRA",
         "filterColumn": "NAME", "filterValue": "p",  "colEnergy": "ENERGY", "colFlux": "FLUX"},
        {"species": "e-", "Z": -1, "mc2": 0.51099895, "file": "spe_shielded.fits", "hdu": "SPECTRA",
         "filterColumn": "NAME", "filterValue": "e-", "colEnergy": "ENERGY", "colFlux": "FLUX"}
    ]
}
```

#### Ground: sea-level muons on a vertical sensor

```json
"cosmic": {
    "depth_um": 10.0, "pitch_um": 18.0, "eVPerElectron": 3.6, "tilt_deg": 90.0,
    "spectrum": [
        {"species": "mu-", "Z": -1, "mc2": 105.658, "file": "muons_sea_level.csv",
         "colEnergy": "E", "colFlux": "J", "angular": {"model": "cosn", "n": 2}, "faces": "both"}
    ]
}
```

For `cosn`, the column `J` holds the intensity in the **vertical** direction.

#### Radioactivity: alphas from the package lid, betas from behind

```json
"cosmic": {
    "depth_um": 10.0, "pitch_um": 18.0, "eVPerElectron": 3.6,
    "spectrum": [
        {"species": "alpha", "Z": 2,  "A": 4,           "file": "lid_alphas.csv",
         "colEnergy": "E", "colFlux": "J", "angular": {"model": "surface"}, "faces": "front"},
        {"species": "beta",  "Z": -1, "mc2": 0.51099895, "file": "substrate_betas.csv",
         "colEnergy": "E", "colFlux": "J", "angular": {"model": "surface"}, "faces": "back"}
    ]
}
```

The alpha is a nucleus (`A`): its energies are **per nucleon**, so a 5 MeV alpha is at E = 1.25.

---

## 3. Spectrum and angular files

### 3.1 Columns and units

The values of the table are used as given: nothing is guessed or converted. The column named by `colEnergy` holds the kinetic energy, in MeV/nucleon for a nucleus and in MeV otherwise, strictly ascending over the rows of a species. The column named by `colFlux` holds the differential intensity $J$, in m⁻² s⁻¹ sr⁻¹ per unit of that energy, and its values must be positive. Any other column, such as `NAME`, is used only to select rows with `filterColumn`.

The energy and flux columns are real numbers (FITS `TFORM` E or D). Their `TUNIT` is not read: the values must be in the units above.

An **angular table** (`"model": "table"`) holds the zenith angle in degrees (`colTheta`, strictly ascending in [0, 90]) and the relative intensity $f$ (`colIntensity`, ≥ 0, not all 0). It is interpolated linearly, and $f = 0$ outside the table.

### 3.2 CSV layout

A CSV file is read by `DSL::FITStable::fromCSV()` in the layout `toCSV()` writes, then filtered like a FITS table (`DSL::FITStable::where()`):

```text
# solar protons and helium behind 2 g/cm2 of aluminium
NAME,ENERGY [MeV/nucleon],FLUX [m-2 s-1 sr-1 (MeV/nucleon)-1]
p,10,1
p,100,0.01
He,10,0.5
He,100,0.005
```

The first line holds the column names, each followed by ` [unit]` when the column has a unit. Fields are separated by commas; `"..."` quotes a field, and `""` is a quote inside a quoted field. An empty field is NULL, which becomes NaN in a numeric column. Blank lines and lines starting with `#` are skipped. A column in which every value is a number is read as real numbers, and any other column is read as strings.

### 3.3 From a surface activity to $J$ (`surface` model)

Consider an isotropic emitter on a plane parallel to the sensor, much larger than the sensor and close to it. Let $a$ be its activity per unit area (Bq m⁻²), $y$ its yield (particles of the species per decay), and $S(E)$ the normalised spectrum of the particles leaving the surface ($\int S\,\mathrm{d}E = 1$). Half of the particles go towards the sensor, uniformly in solid angle, so that

```math
J(E) = \frac{a\,y\,S(E)}{4\pi}\qquad(\text{m}^{-2}\,\text{s}^{-1}\,\text{sr}^{-1}\text{ per unit of }E)
```

With `"angular": {"model": "surface"}` the rate is then $2\pi\int J\,\mathrm{d}E\cdot A = a\,y\,A/2$ per face. $S(E)$ is the spectrum **leaving** the surface: the energy lost inside the source material must be accounted for by the user.

---

## 4. Sensor material

`CosmicDeposit::Params::material` (`target::material`) holds the mean $\langle Z/A\rangle$, the mean excitation energy $I$, the density, and the Sternheimer coefficients of the density effect.

By default, the material is **silicon**, with $\langle Z/A\rangle = 0.49848$, $I = 173$ eV and $\rho = 2.329$ g/cm³.

`material::compound(elements, density)` builds a solid from its composition. It computes $\langle Z/A\rangle$ from the elements and $I$ by Bragg additivity, the constant of the density effect $\bar C = 2\ln(I/\hbar\omega_p) + 1$, with $\hbar\omega_p$ the plasma energy of the material, and the Sternheimer coefficients $x_0$, $x_1$, $a$ and $m$ from the general formulas for solids of Sternheimer & Peierls (1971).

`material::HgCdTe(x)` builds Hg₁₋ₓCdₓTe, $x$ being the Cd fraction. The elements are taken from the PDG, with their mean excitation energy $I$ multiplied by 1.13 for atoms bound in a condensed compound (ICRU 37). The density is computed from the zinc-blende cell, with the lattice constant of Higgins et al. (1989). For $x = 1$, the material reproduces the PDG values of CdTe to 0.1 %.

```cpp
CosmicDeposit::Params p;
p.depth_um      = 10.0;
p.pitch_um      = 18.0;
p.eVPerElectron = pairEnergy_eV;   // mean energy per electron-hole pair of your sensor
p.material      = phot2pix::CosmicRay::target::material::HgCdTe(0.45);
```

---

## 5. Outputs

### 5.1 Header cards (`inputCards()`)

| Card | Unit | Value |
|---|---|---|
| `CRFLUX` | m⁻² s⁻¹ sr⁻¹ | $\sum_s \int J_s\,\mathrm{d}E$ |
| `CRACC` | sr | Geometric factor, weighted by the flux of each species |
| `CRTILT` | deg | `tilt_deg` |
| `CRDEPTH` | µm | `depth_um` |
| `CRPITCH` | µm | `pitch_um` |
| `CREVPE` | eV | `eVPerElectron` |
| `CRSIGD` | pixel | `sigmaDiffusion_px` of `"cosmic"` (the footprints of the entries are not written) |
| `CRMAXTRK` | pixel | `maxTrackPixels` |
| `CRDCUT` | e⁻ | `deltaRayCut_e` |
| `CRTRUNC` | — | Share of the hits not drawn (`truncatedFraction()`) |
| `CRSEED` | — | `seed` |

### 5.2 Table of hits (`toColumns()`, `columnInfo()`)

`generate(dt, geometry, out, &events)` also returns the hits. The hits are a complete record of the step: since `depositEvent()` is deterministic, the deposits can be recomputed from the hits and the configuration.

| Column | Unit | Meaning |
|---|---|---|
| `SPECIES` | — | Name of the species |
| `Z` | — | Charge number |
| `A` | — | Mass number |
| `ENERGY` | MeV | Kinetic energy, per nucleon for a nucleus, per particle otherwise |
| `MU` | — | $\cos\theta$ to the sensor normal |
| `PHI` | rad | Azimuth of the track on the matrix |
| `FACE` | — | Face crossed at entry: +1 front, −1 back |
| `X`, `Y` | pixel | Entry point (pixel $i$ covers $[i, i+1)$; FITS position: $x + 0.5$) |
| `DX`, `DY` | pixel | Projected displacement of the path in the layer |
| `PATH` | µm | Path in the layer: $d/\mu$, or the range if the particle stops |
| `NE` | electron | Charge deposited in the layer |

---

## 6. Using it from a detector family

`CosmicDeposit` is a `phot2pix::ChargeDeposit` (`ChargeDeposit.h`), the base of the sources of charge over several pixels:

| Member | Role |
|---|---|
| `configure(cfg)` | Reads the `"cosmic"` key; nothing if it is absent |
| `isConfigured()` | True when the spectrum and every required parameter are set |
| `inputCards()` | Header cards of the inputs ([§5.1](#51-header-cards-inputcards)) |
| `generate(dt, geometry, out)` | Appends the deposits of one step, as `{idx, electrons}`, `idx = Y * width + X`; throws `std::logic_error` if the source is not configured |
| `clone()` | Copy with its dynamic type |
| `setSeed(s)` | With a seed, the n-th call of `generate()` always draws the same stream, whatever the threads |

A detector family adapts the model by deriving from `CosmicDeposit`. Its constructor sets the default parameters `m_params` (the material, for instance) and the spectrum `m_spectrum`, and then calls **`buildTables()`**, which builds the range tables, the geometric factors, the draw weights and the lateral spreads; without this call, `drawEvents()` throws. To change the fluctuations of thin crossings, for instance to use a measured distribution, the family overrides `energyDeposit_eV()`. Particles that stop and thick crossings are always handled by `drawEvents()`.

A `phot2pix::DetectorSimulator` uses this source through `expose()`. When the configuration has a `"cosmic"` key, `initialize()` creates the source returned by the virtual function `_makeChargeDeposit()`, which is a silicon `CosmicDeposit` by default, configures it, and records its header cards. `expose()` then draws the deposits of the exposure time in a thread of its own while it simulates the pixels as `simulate()` does, adds the deposits that fall on the matrix (or on the window), divided by the gain of their pixel, and converts the image to the pixel type. `simulate()` itself never adds deposits. Since an exposure of a `DetectorSimulator` is a single read, every hit of the exposure time is in its image. A detector family chooses its source by overriding `_makeChargeDeposit()`, for instance to set another material. `nisp::SCA` returns `nullptr`: it owns its HgCdTe source (`nisp::CosmicDeposit`) and adds the hits at each step of its readout, so that they accumulate along the ramp and leave persistence ([Philosophy §5](Philosophy#example-a-source-of-cosmic-ray-hits)). `nisp::persistence::pixelArray` also returns `nullptr`, because it models the charge stored by past exposures and not a detector.

---

## 7. Accessors without simulating

| Function | Value |
|---|---|
| `eventRate(area)` | $\sum_s G_s\int J_s\,\mathrm{d}E\cdot A$ (s⁻¹), every track length |
| `expectedEvents(area, dt)`, `probability(n, area, dt)`, `probabilityAtLeastOne(area, dt)` | Poisson |
| `geometricFactor(i)`, `truncatedFraction(i)`, `truncatedFraction()` | Per species, and for all species |
| `stoppingPower(species, E)` | Mean collision stopping power (MeV cm²/g) |
| `stoppingTable(i).range(T)`, `.energy(R)`, `.energyAfter(T, s)` | CSDA range (g/cm²) and its inverse, $T$ in MeV |
| `energyLoss(species, E, path_um)` | Landau most probable loss and width (eV) |
| `deltaRayShare(i, E)` | Share $h$ of the charge of a hit in its δ-ray halo ([§1.5](#15-delta-rays-deltarays-true)); 0 without `deltaRays` |

### 7.1 Geometry of a sensor, without a source

The free functions of `phot2pix::CosmicRay` describe a flat sensor that nothing shields. They need no spectrum and no `CosmicDeposit`.

They use the following conventions. The sensor normal (front face) is tilted by `tilt_deg` from the zenith, towards $+x$ of the matrix. A sky direction is given by its zenith angle `thetaZ_deg`, which is 0 for a particle coming from above and 180 for a particle coming from below, and by its azimuth `phiZ_deg`, which is 0 towards $+x$ of the matrix. The incidence angle $\theta$ is measured from the sensor normal.

| Function | Value |
|---|---|
| `incidence(thetaZ_deg, phiZ_deg, tilt_deg)` | `Incidence{mu, phi, face}` of a particle coming from that direction: $\mu = \lvert\cos\theta\rvert$, azimuth of its track on the matrix, face reached |
| `effectiveArea(area_m2, thetaZ_deg, phiZ_deg, tilt_deg)` | Area presented to a parallel beam from that direction: $A\,\lvert\cos\theta\rvert$ (m²) |
| `geometricFactor(model, tilt_deg, thetaMin_deg, thetaMax_deg)` | $G = \int f\,\mu\,\mathrm{d}\Omega$ (sr) of an angular model ([§1.1](#geometric-factor)), over its faces and the incidence angles $\theta \in [\theta_{\min}, \theta_{\max}]$ (full sky by default) |
| `etendue(model, area_m2, tilt_deg, thetaMin_deg, thetaMax_deg)` | $A\,G$ (m² sr): hits per second $= A\,G\int J\,\mathrm{d}E$ |

For one face, $G$ has the closed forms below; in the other cases, it is integrated numerically to about 10⁻³ relative precision.

| Model | $G$ over $\theta \in [\theta_{\min}, \theta_{\max}]$ |
|---|---|
| `isotropic` | $\pi\,(\cos^2\theta_{\min} - \cos^2\theta_{\max})$ |
| `surface` | $2\pi\,(\cos\theta_{\min} - \cos\theta_{\max})$ |
| `cosn`, horizontal sensor | $\frac{2\pi}{n+2}\,(\cos^{n+2}\theta_{\min} - \cos^{n+2}\theta_{\max})$, front face only |

For example, cos² muons on a vertical sensor (`tilt_deg` = 90) give $G = \pi/8$ per face, a quarter of the $\pi/2$ of a horizontal sensor.

`CosmicDeposit` uses the same functions: its `geometricFactor(i)` is `geometricFactor(model, tilt_deg)` over the full sky.

---

## 8. Validation

`src/test/gTest_chargeDeposit.cxx` tests the `ChargeDeposit` base class. It checks the seeds and the reproducibility of the random streams, the copies and `clone()`, which continue the sequence, 64 threads drawing distinct streams, and `push()`, which drops the pixels outside the detector and the empty charges.

`src/test/gTest_cosmicDeposit.cxx` tests the source itself, as summarised in the table below.

| Part | Checks |
|---|---|
| Species | energy per nucleon or per particle; $\gamma$ of a 1 GeV/n proton (2.0735) and of a 1 MeV electron (2.957) |
| Spectrum | closed-form flux of a bin (g = 2, 1, 0); probabilities; KS test of the sampled energies; loading from CSV and FITS (filter, required columns, `hdu`, errors) |
| Materials | `compound()` of silicon; `HgCdTe(1)` against the PDG CdTe values |
| Configuration | every key; a bad load leaves the source unchanged; `acceptance_sr` rejected |
| Probabilities, generation | rate, Poisson, truncation; KS tests of the directions; same seed, same hits |
| Footprint | the length in each pixel; charge conserved, with and without lateral spread; detector edges; per-species core and halo (shares, charge conserved, configuration and its errors) |
| Delta rays | no halo below $T_0$; $h$ in (0, ½) and the same for p and He at the same energy per nucleon; charge kept; halo decreasing with the distance, with the symmetries of the pixel grid; radius growing with the charge, none under `deltaRayCut_e`; inclined tracks; configuration and its errors (not a boolean, e±, cut ≤ 0) |
| Stopping and range | `StoppingTable` against an analytical case; **Bethe against the PDG muon tables of Si, Ge and CdTe**: ionisation within 1 % at 10 MeV, 100 MeV and 1 GeV, CSDA range within 2–3 %; range scaling of helium; e± consistency |
| Bragg | a stopping particle deposits all its energy; charge rising to the Bragg peak; uniform charge of a fast particle; thick crossings average the CSDA loss |
| Angular models | geometric factors (π, 2π, π/2, π/8 for cos² on a vertical sensor); truncated shares; KS tests of $\mu$; faces; the rejection sampler ($\langle\mu\rangle = 8/15$ for cos² on a vertical sensor); configuration |
| Geometry helpers | `incidence()` for known beams (normal, 60°, from below, tilted sensors); `effectiveArea()`; `geometricFactor()` closed forms with windows of incidence angles (π/4, π, 15π/32), numerical cases, windows adding up to the full sky, the same value as the source; `etendue()` |

No independent table could be fetched for electrons (NIST ESTAR needs a form), so the e± stopping power is only checked for consistency.

---

## 9. Known limitations

- **Electrons follow straight tracks**, because the model has no multiple scattering. The effect is rarely seen on NISP. Adding it would require condensed-history steps with a Highland deflection, a radiation length in `target::material`, and a per-event seed in `Event`, so that `depositEvent()` stays deterministic.
- **Thin crossings use the Landau distribution even where the Vavilov distribution applies**, that is, where $\kappa = \xi/W_{\max} > 0.01$, for example for protons below about 100 MeV in 10 µm of HgCdTe. In that regime, the width and the upper tail of the deposit are too large (by about 1.7 at $\kappa = 1$), while its typical value is correct.
- **The low-energy stopping power is simplified.** It has no shell, Barkas or Bloch corrections and no effective charge, and it is taken proportional to the velocity below its maximum. Low-energy ranges are therefore too long in heavy materials: a 5 MeV alpha in HgCdTe travels about 40 µm instead of about 20 µm. Reading a stopping table given by the user (SRIM, NIST ASTAR/PSTAR) is left for a later version.
- **Some physics is missing.** Electrons have no radiative loss, δ-rays are not transported (`deltaRays` and `halo` spread a share of the charge, but they draw no δ-ray track), and there is no nuclear interaction and no secondary particle.
- **The δ-ray halo** ([§1.5](#15-delta-rays-deltarays-true)) has four simplifications. It is computed at the entry velocity, also for a particle that stops. It uses the CSDA ranges of the electrons, which are longer than their real reach in a heavy material such as HgCdTe, where the electrons scatter. Its escape from the layer is a geometric estimate, $d/r$. Finally, e⁻ and e⁺ primaries get no halo.
- **Photons are not followed.** For gamma rays, the user must give the spectrum of their secondary electrons.
- **Grazing tracks** longer than `maxTrackPixels` are not drawn (`CRTRUNC`).
- **The silicon defaults differ slightly from the PDG 2024 tables.** The defaults of `target::material` ($\bar C$ = 4.4351, $x_0$ = 0.2014, $x_1$ = 2.8715) differ from the PDG values (4.4355, 0.2015, 2.8716) by about 10⁻⁵ on dE/dx.

---

## 10. References

- PDG, "Passage of particles through matter" (D. E. Groom, S. R. Klein), in S. Navas *et al.* (Particle Data Group), *Phys. Rev. D* **110** (2024) 030001: Bethe equation, density effect, fluctuations in energy loss, energetic knock-on electrons (δ-rays).
- J. J. Butts, R. Katz, "Theory of RBE for heavy ion bombardment of dry enzymes and viruses", *Radiat. Res.* **30** (1967) 855 (radial dose of the δ-rays).
- PDG Atomic and Nuclear Properties (2024), muon tables: <https://pdg.lbl.gov/2024/AtomicNuclearProperties/>.
- R. M. Sternheimer, M. J. Berger, S. M. Seltzer, *At. Data Nucl. Data Tables* **30** (1984) 261; R. M. Sternheimer, *Phys. Rev.* **88** (1952) 851; R. M. Sternheimer, R. F. Peierls, *Phys. Rev. B* **3** (1971) 3681.
- H. Bichsel, "Straggling in thin silicon detectors", *Rev. Mod. Phys.* **60** (1988) 663.
- L. Landau, *J. Phys. USSR* **8** (1944) 201; K. S. Kölbig, B. Schorr, *Comput. Phys. Commun.* **31** (1984) 97 (CERNLIB G110).
- ICRU Report 37, "Stopping powers for electrons and positrons" (1984); F. Rohrlich, B. C. Carlson, *Phys. Rev.* **93** (1954) 38.
- S. M. Seltzer, M. J. Berger, *Int. J. Appl. Radiat. Isot.* **33** (1982) 1189.
- Higgins *et al.*, *J. Vac. Sci. Technol. A* **7** (1989) 271 (lattice constant of Hg₁₋ₓCdₓTe).

---

*The Phot2Pix logo has been generated by OpenAI's ChatGPT Image Generation tool (ImageGen).*

*This wiki page has been generated by AI (Claude Opus 5.5) by reading the code itself.*

<!-- This page is generated from wiki/CosmicDeposit.md in the phot2pix repository. Edit it there: changes made in the web interface are overwritten on the next push. -->
