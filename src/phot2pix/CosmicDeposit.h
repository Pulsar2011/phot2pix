// Created by gillard on 2026-11-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Written with the assistance of Claude Opus 5.5 (Anthropic's AI model) under guidance of William Gillard.
//       Claude Opus 5.5 was used to generate code snippets, comments, and documentation (Doxygen and wiki). It was used to generate unitary test cases, and to assist in code refactoring for optimisation in speed. It was particullarly used in the generation of the Delta ray and for the reading of the JSON file for the cosmic ray configuration. The AI model was guided by the developer's instructions and feedback to ensure the generated code met the project's requirements.
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
// You may share and adapt this code with attribution, Non Comercial use only.
//
// CosmicDeposit.h
//
/// @file CosmicDeposit.h
/// @brief phot2pix::CosmicDeposit, cosmic ray hits drawn from a cosmic ray spectrum.
///
/// Purpose: the impact of cosmic rays and radioactivity on the images, not cosmic-ray physics. The model is
/// the simplest that gets right what an image sees: the charge of a hit, the length of its track, where it
/// ends (Bragg peak) and the pixels it touches. The spectrum given is the one at the detector, behind the
/// shielding (from SPENVIS, CREME96, Geant4...): nothing is transported through the shielding here.
///
/// Model, for a charge deposition (step) of duration dt over the physical detector of area A:
///  1. Rate. The spectrum is the differential intensity J(E) of each species, E being the kinetic
///     energy per nucleon for a nucleus, per particle otherwise: [m-2 s-1 sr-1 (MeV/nucleon)-1] or
///     [m-2 s-1 sr-1 MeV-1]. The number of hits is Poisson, of mean
///        lambda = sum_species G * (∫ J dE) * A * dt,
///     G being the geometric factor of the species (sr), from its angular model (AngularModel): pi per
///     face for an isotropic field, 2 pi per face for a surface emitter, 2 pi / (n + 2) for cos^n(theta_z)
///     on a horizontal sensor, computed numerically otherwise. Which species and which energy a hit has
///     follows the cells of the spectrum, J being interpolated as a power law between two table points.
///  2. Direction, in the frame of the sensor: theta from its normal, mu = cos(theta), phi the azimuth on
///     the matrix, and the face crossed (Event::face). Isotropic: mu^2 uniform; surface: mu uniform;
///     cos^n on a horizontal sensor: mu^(n+2) uniform; otherwise drawn by rejection from the sky, the
///     sensor normal being tilted by tilt_deg from the zenith towards +x of the matrix.
///     The path through the sensitive depth d is d / mu, and its projection on the matrix
///     (d / pitch) * tan(theta) pixels. Tracks longer than maxTrackPixels are not drawn; their share
///     of the hits is truncatedFraction(), the share of the flux with mu < 1 / sqrt(1 + q^2),
///     q = maxTrackPixels * pitch / d (1 / (1 + q^2) for an isotropic field).
///  3. Energy deposited, from the CSDA range R of the particle (StoppingTable) and its path L = d / mu:
///     - stops (R <= L): the path is R, the projected track is shortened to R / L, and the particle
///       deposits all its kinetic energy;
///     - thick (CSDA mean loss over L above 10 % of its energy): Gaussian around the mean loss (Bohr);
///     - thin: Landau distribution of the thin layer (most probable value and width from the PDG,
///       "Passage of particles through matter"), limited to the kinetic energy of the particle.
///     Along the path, the charge follows the CSDA energy loss: uniform for a fast particle, rising to a
///     Bragg peak for one that stops.
///  4. Footprint. Exact length of the track in each pixel, then an optional lateral spread of the
///     species: a Gaussian core of sigma sigmaDiffusion_px (the entry's, else the source's) and an
///     optional wider Gaussian halo holding a share f of the charge,
///        (1 - f) G(sigma_core) + f G(sigma_halo),
///     each integrated over the pixel and cut at 3 sigma (the charge is kept). A hit is one pixel
///     (normal incidence, no lateral spread) or several (inclined track, lateral spread). The halo gives
///     the broad, round charge of a snowball or of the delta rays around a heavy ion, not their physics.
///     The deposit lies in the sensitive layer: no quantum efficiency, and a hit can land on any pixel
///     of the physical detector, reference pixels included.
///  5. Delta rays (entries with "deltaRays": true, heavy particles only). The knock-on electrons of a hit
///     carry a share h of its charge into a halo around each pixel of the track. h and the halo follow
///     from the velocity of the particle at entry (beta, gamma), with closed formulas only:
///      - Wmax, largest energy given to an electron; T0, energy of an electron whose CSDA range is half a
///        pixel: a delta ray below T0 stays in the pixel of the track (no halo when Wmax <= T0).
///      - Share of the mean energy loss carried by the delta rays from T0 to Wmax, their number being
///        dN/dT ∝ (1 - beta^2 T / Wmax) / T^2, against the Bethe mean loss (z^2 cancels):
///          f = [ln(Wmax / T0) - beta^2 (1 - T0 / Wmax)] / [ln(2 m c^2 beta^2 gamma^2 Wmax / I^2) - 2 beta^2 - delta].
///      - Radial profile, from the 1/r^2 radial dose of the delta rays (Butts & Katz): the same energy in each
///        interval of ln r, from r0 = 1/2 pixel to rMax, the CSDA range of an electron of energy Wmax. A delta
///        ray going farther than the depth d leaves the layer and keeps about d / r of its energy in it:
///        weight w(r) = min(1, d / r) per interval of ln r.
///      - Charge of the hit, kept: fk = f ∫ w dln r / ln(rMax / r0) is deposited by the delta rays, so the
///        halo holds h = fk / (1 - f + fk) of the charge and the track 1 - h (then its lateral spread).
///      - The halo of a pixel of the track (charge h * electrons / track length) stops at the radius where it
///        puts less than deltaRayCut_e on a pixel; its shares are normalised within that radius.
///     The halo is a mean profile: no delta-ray track and no fluctuation of h.
///
/// Configuration parameters of the "cosmic" object of a PixelConfig, read by configure(). There is
/// no default value but for tilt_deg, sigmaDiffusion_px, maxTrackPixels and deltaRayCut_e:
/// @verbatim
///   "depth_um"          : sensitive depth d (um)
///   "pitch_um"          : pixel pitch (um)
///   "eVPerElectron"     : energy per electron-hole pair (eV)
///   "tilt_deg"          : angle of the sensor normal to the zenith (deg), [0, 180]  [0]
///                         (zenith models only); "acceptance_sr" is no longer read: an error
///   "sigmaDiffusion_px" : lateral spread, sigma (pixel), 0: none; default of the species [0]
///   "maxTrackPixels"    : longest projected track drawn (pixel)  [64]
///   "deltaRayCut_e"     : smallest charge a delta-ray halo puts on a pixel (e-), > 0  [1]
///   "seed"              : seed of the random streams, 0: none    [0]
///   "spectrum"          : one entry per species, its table read from a CSV file or a FITS table (below)
/// @endverbatim
///
/// Spectrum. Each species names the file of its table:
/// @verbatim
///   "spectrum": [
///     {"species": "p",  "Z": 1, "A": 1, "file": "gcr_p.csv", "colEnergy": "ENERGY", "colFlux": "FLUX"},
///     {"species": "He", "Z": 2, "A": 4, "file": "gcr.fits", "hdu": "GCR",
///      "filterColumn": "NAME", "filterValue": "He", "colEnergy": "ENERGY", "colFlux": "FLUX"},
///     {"species": "e-", "Z": -1, "mc2": 0.51099895, "file": "electrons.csv", "colEnergy": "E", "colFlux": "J"},
///     {"species": "Fe", "Z": 26, "A": 56, "file": "gcr_fe.csv", "colEnergy": "ENERGY", "colFlux": "FLUX",
///      "sigmaDiffusion_px": 0.5, "deltaRays": true}
///   ]
///   "species"      : name of the species
///   "Z"            : charge number (-1 for an electron)
///   "A"            : mass number of a nucleus: energies per nucleon; or, instead of it,
///   "mc2"          : rest energy (MeV) of any other particle (e-, e+, muon...): energies per particle
///   "file"         : CSV file (extension .csv, any case), otherwise FITS file
///   "hdu"          : FITS only, required: HDU of the table, index (1 = primary, so 2 = first extension) or EXTNAME
///   "filterColumn" : when the table holds several species, the rows of this one are those where
///   "filterValue"  :   filterColumn == filterValue (a string or a number); both or neither         [every row]
///   "colEnergy"    : name of the energy column (required)
///   "colFlux"      : name of the flux column (required)
///   "angular"      : angular distribution (AngularModel)                                     [isotropic]
///                      {"model": "isotropic"} | {"model": "surface"} | {"model": "cosn", "n": 2} |
///                      {"model": "table", "file", "hdu", "colTheta", "colIntensity"[, "filterColumn",
///                       "filterValue"]}: zenith angle (deg) and relative intensity f, read as the spectra
///   "faces"        : faces of the sensor the particles reach, "front" | "back" | "both"       [front]
///   "sigmaDiffusion_px" : sigma of the core of the lateral spread of this species (pixel), >= 0
///                                                                    [sigmaDiffusion_px of "cosmic"]
///   "halo"         : wider Gaussian holding a share of the charge of this species           [none]
///                      {"sigma_px": sigma of the halo (pixel) > 0, "fraction": share of the charge in [0, 1]}
///   "deltaRays"    : true: delta-ray halo of each hit, from its velocity (5. above); not for e- or e+  [false]
/// @endverbatim
/// The table holds the spectrum in these columns, and it is up to the user to write them so:
///  - colEnergy: kinetic energy, per nucleon for a nucleus (MeV/nucleon), per particle otherwise (MeV),
///               strictly ascending over the rows of a species;
///  - colFlux:   differential intensity J (m-2 s-1 sr-1 per unit of that energy), > 0;
///  - any other column, such as NAME (the species), only used to select rows with filterColumn.
/// The energy and flux columns are real numbers (FITS TFORM E or D). Their units are not read (TUNIT is
/// ignored): the values must be in the units above.
/// A CSV file is read by DSL::FITStable::fromCSV(), in the layout toCSV() writes: a line of column names,
/// each followed by " [unit]" when it has one ("ENERGY [MeV/nucleon]"), then one line per row. Fields are
/// separated by ','; "..." quotes a field and "" is a quote inside it; an empty field is NULL (NaN in a
/// number column); blank lines and lines starting with '#' are skipped. A column whose every value is a
/// number is read as real numbers, any other as strings. The CSV file is turned into a DSL::FITStable,
/// so both kinds of file are read and filtered the same way (DSL::FITStable::where()).
///
/// Known limitations:
///  - Electrons follow straight tracks: no multiple scattering (rarely seen on NISP). A later version would
///    need condensed-history steps with a Highland deflection, a radiation length in target::material and a
///    per-event seed in Event so that depositEvent() stays deterministic.
///  - Thin crossings use Landau even where kappa = xi / Wmax > 0.01 (Vavilov regime, e.g. protons below about
///    100 MeV in 10 um of HgCdTe): the width and the upper tail of the deposit are too large (about 1.7 times
///    at kappa = 1), not its typical value. No Vavilov distribution.
///  - Stopping power without shell, Barkas or Bloch corrections nor effective charge, and proportional to the
///    velocity below the energy where the formula is largest: low-energy ranges are too long in heavy
///    materials (a 5 MeV alpha in HgCdTe: about 40 um instead of about 20 um). A stopping table given by the
///    user (SRIM, NIST ASTAR/PSTAR) is for later.
///  - No radiative loss of electrons (the bremsstrahlung photons leave the layer), no delta-ray transport
///    ("deltaRays" and "halo" spread a share of the charge, they draw no delta-ray track), no nuclear
///    interaction or secondary particle.
///  - Delta rays: mean halo from the entry velocity (also for a particle that stops), CSDA ranges of the
///    electrons (longer than their real reach in a heavy material, where they scatter), escape from the
///    layer estimated as d / r; none for e- and e+ primaries.
///  - Photons (gamma rays) are not followed: give the spectrum of their secondary electrons instead.
///  - Tracks longer than maxTrackPixels are not drawn (truncatedFraction()).
///  - The silicon defaults of target::material (Cbar 4.4351, x0 0.2014, x1 2.8715) differ slightly from the
///    PDG 2024 tables (Cbar 4.4355, x0 0.2015, x1 2.8716, a 0.1492): a change of about 1e-5 on dE/dx.

#ifndef P2P_COSMIC_DEPOSIT_H
#define P2P_COSMIC_DEPOSIT_H

#include "DetectorEffect.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace phot2pix
{
#pragma region - spectrum

    namespace CosmicRay
    {
        constexpr double MeV_per_u = 931.49410242;   ///< 1 u = 931.49410242 MeV/c2
        constexpr double N_A       = 6.02214076e23;  ///< Avogadro Cte  mol-1

        namespace target
        {
            /// @brief One element of a compound (see material::compound()).
            struct element
            {
                double Z;      ///< Atomic number
                double A;      ///< Molar mass (g/mol)
                double I_eV;   ///< Mean excitation energy of the atom in the compound (eV)
                double n;      ///< Atoms per formula unit, >= 0
            };

            /// @brief Sensor material, for the energy loss. Defaults: silicon.
            // Silicon, PDG Atomic and Nuclear Properties (2024), muon table header
            // (Groom, Mokhov & Striganov, At. Data Nucl. Data Tables 78 (2001) 183).
            // Density effect: Sternheimer, Berger & Seltzer, At. Data Nucl. Data Tables 30 (1984) 261.

            struct material
            {
                double ZoverA            = 0.49848;   ///< <Z/A>
                double meanExcitation_eV = 173.0;     ///< Mean excitation energy I (eV)
                double density_gcm3      = 2.329;     ///< Density (g/cm3)
                // Density effect, Sternheimer parametrisation. Silicon values from memory of the
                // Sternheimer-Berger-Seltzer table (PDG, atomic and nuclear properties): check them.
                double sternheimerC      = 4.4351;    ///< -C
                double sternheimerX0     = 0.2014;    ///< x0
                double sternheimerX1     = 2.8715;    ///< x1
                double sternheimerA      = 0.14921;   ///< a
                double sternheimerM      = 3.2546;    ///< m
                double sternheimerDelta0 = 0.14;      ///< delta0 (conductors; 0 for an insulator)
                double sternheimerDeltaMax = 0.059;   ///< deltaMax

                /**
                 * @brief Solid made of elements, from its composition and its density.
                 * @details <Z/A> = sum(n Z) / sum(n A).
                 *          I: Bragg additivity, ln I = sum(n Z ln I_i) / sum(n Z), I_i being the value of the atom
                 *          in the compound.
                 *          Density effect: Cbar = 2 ln(I / hbar omega_p) + 1, hbar omega_p = 28.816 sqrt(rho <Z/A>) eV.
                 *          No fitted x0, x1, a, m: general formulas for a solid of Sternheimer & Peierls,
                 *          Phys. Rev. B 3 (1971) 3681, a making delta continuous at x0; delta0 = 0 (insulator,
                 *          semiconductor). Tried on CdTe, delta differs from the PDG fit by ~0.2 at minimum
                 *          ionisation and ~0.4 at beta gamma = 10. sternheimerDeltaMax: NaN (no fit).
                 * @param elements     Elements and their number of atoms per formula unit.
                 * @param density_gcm3 Density (g/cm3).
                 * @throws std::invalid_argument if a value is not valid.
                 */
                static material compound(const std::vector<element>& elements, double density_gcm3);

                /**
                 * @brief Hg(1-x) Cd(x) Te, x being the Cd fraction.
                 * @details Elements from PDG (2024): Hg 200.592 g/mol, I = 800 eV; Cd 112.414 g/mol, 469 eV;
                 *          Te 127.60 g/mol, 485 eV. I multiplied by 1.13 for atoms bound in a condensed compound
                 *          (ICRU Report 37; Seltzer & Berger, Int. J. Appl. Radiat. Isot. 33 (1982) 1189).
                 *          Density: zinc blende, 4 formula units per cubic cell, lattice constant by Vegard's law
                 *          between HgTe (6.4614 A) and CdTe (6.481 A).
                 *          x = 1 gives the PDG CdTe values (<Z/A> 0.41665, I 539.3 eV, rho 5.850 g/cm3, Cbar 5.9677)
                 *          to 0.1 %. For another density, call compound() with the same elements.
                 * @param x Cd fraction, in [0, 1].
                 * @throws std::invalid_argument if x is not in [0, 1].
                 */
                static material HgCdTe(const double& x=0.3);
            };
        }

#pragma region - species
        /// @brief define the cosmic ray specie.
        class specie
        {
            private:
                std::string m_name;                ///< Name ("p", "He", ...)
                double      m_Z   = 0.0;           ///< Charge number
                double      m_A   = 0.0;           ///< Mass number, > 0
                double      m_mc2 = 0.0;           ///< Rest energy in unit of MeV
                bool        m_perNucleon = true;   ///< Energy of the spectrum per nucleon (nucleus), else per particle

                specie(const std::string& name, const double& Z, const double& A, const double& mc2, bool perNucleon)
                    : m_name(name), m_Z(Z), m_A(A), m_mc2(mc2), m_perNucleon(perNucleon) { validate(); }

#pragma region -- ctor/dtor
            public:
                /// @brief Nucleus from its mass number; rest energy A * MeV_per_u; energies per nucleon.
                specie(std::string name, double Z, double A)   : specie(std::move(name), Z, A, A * MeV_per_u, true) {}

                /// @brief Any particle from its rest energy (MeV); A = mc2 / MeV_per_u; energies per particle.
                static specie fromRestEnergy(const std::string& name, const double& Z, const double& mc2) { return specie(name, Z, mc2 / MeV_per_u, mc2, false); }
                
                specie() = default;
                virtual ~specie() = default;
#pragma endregion

#pragma region -- accessors
                const std::string& name() const { return m_name; }
                double Z()        const { return m_Z; }
                double A()        const { return m_A; }
                double restMass() const { return m_mc2; }
                double mc2()      const { return restMass(); }
#pragma endregion

#pragma region -- kinematics
                /// @brief True when the energies of the spectrum are per nucleon (nucleus), false when per particle.
                bool perNucleon() const { return m_perNucleon; }

                /// @brief Kinetic energy of the particle (MeV) for an energy E of the spectrum:
                ///        T = E * A for a nucleus (E per nucleon), T = E for any other particle (E per particle).
                double kineticEnergy(const double& E) const { return m_perNucleon ? E * m_A : E; }

                /// @brief Lorentz factor for an energy E of the spectrum: 1 + kineticEnergy(E) / mc2.
                ///        For a nucleus, 1 + E / MeV_per_u.
                double gamma(const double& E) const { return 1.0 + kineticEnergy(E) / m_mc2; }
#pragma endregion

#pragma region -- validation
            void validate() const
            {
                if (!(m_A > 0.0)) throw std::invalid_argument("specie \033[34m" + m_name + "\033[0m: A <= 0");
                if (std::abs(m_mc2 - m_A * MeV_per_u) > 1e-6)
                    std::cerr << "specie \033[34m" << m_name << "\033[0m: mc2 = \033[34m" << m_mc2
                              << "\033[0m != A * MeV_per_u = \033[34m" << m_A * MeV_per_u << "\033[0m" << std::endl;
                if (m_name.empty())
                    std::cerr << "specie (A = \033[34m" << m_A << "\033[0m, Z = \033[34m" << m_Z << "\033[0m): \033[31mname is empty\033[0m" << std::endl;
            }
#pragma endregion
        };

        inline const specie proton("p", 1.0, 1.0);
        inline const specie helium("He", 2.0, 4.0);
        inline const specie electron = specie::fromRestEnergy("e-", -1.0, 0.51099895069);
        inline const specie positron = specie::fromRestEnergy("e+",  1.0, 0.51099895069);
#pragma endregion

#pragma region - angular distribution
        /**
         * @class AngularModel
         * @brief Angular distribution of the particles of one species, and the faces of the sensor they reach.
         * @details J is the intensity of the spectrum (m-2 s-1 sr-1 per unit of energy). Models:
         *  - isotropic: intensity J in every direction. The flux through the sensor goes as mu = cos(theta):
         *    mu^2 uniform, geometric factor G = pi sr per face;
         *  - surface: isotropic emitter on a surface parallel to the sensor (package lid, contamination):
         *    flux uniform in solid angle, mu uniform, G = 2 pi sr per face;
         *  - cosn: zenith intensity I(theta_z) = J cos^n(theta_z), J the vertical intensity, downward-going
         *    particles only (sea-level muons: n = 2). Horizontal sensor: mu^(n+2) uniform, G = 2 pi / (n + 2) sr;
         *  - table: zenith intensity I(theta_z) = J f(theta_z), f tabulated in degrees (linear interpolation,
         *    0 outside the table), downward-going particles only.
         *  The sensor-frame models (isotropic, surface) do not depend on the orientation of the sensor; the zenith
         *  models (cosn, table) depend on its tilt (CosmicDeposit::Params::tilt_deg).
         */
        class AngularModel
        {
        public:
            enum class Kind  { isotropic, surface, cosn, table };   ///< Model
            enum class Faces { front, back, both };                 ///< Faces of the sensor the particles reach

            /// @brief Isotropic, through the front face.
            AngularModel() = default;

            static AngularModel isotropic(Faces faces = Faces::front);
            static AngularModel surface(Faces faces = Faces::front);
            /// @throws std::invalid_argument if n < 0.
            static AngularModel cosn(double n, Faces faces = Faces::front);
            /// @param theta_deg Zenith angles (deg), strictly ascending in [0, 90], at least 2.
            /// @param f         Relative intensity at each angle, finite, >= 0, not all 0.
            /// @throws std::invalid_argument on invalid values.
            static AngularModel table(std::vector<double> theta_deg, std::vector<double> f, Faces faces = Faces::front);

            /// @brief From the "angular" object and the "faces" key of a spectrum entry (see CosmicDeposit.h);
            ///        isotropic through the front face without them.
            /// @throws std::invalid_argument on a malformed value or a table that cannot be read.
            static AngularModel load(const PixelConfig& entry);

            Kind   kind()  const { return m_kind; }
            Faces  faces() const { return m_faces; }
            /// @brief Exponent of the cosn model.
            double n()     const { return m_n; }
            /// @brief True if the particles reach the front face, the back face.
            bool   front() const { return m_faces != Faces::back; }
            bool   back()  const { return m_faces != Faces::front; }
            /// @brief True for the models defined in the sky frame (cosn, table).
            bool   zenith() const { return m_kind == Kind::cosn || m_kind == Kind::table; }

            /// @brief I(theta_z) / J of a zenith model, cosZenith = cos(theta_z) >= 0 (downward-going); 1 otherwise.
            double relativeIntensity(double cosZenith) const;
            /// @brief Largest relativeIntensity().
            double maxRelativeIntensity() const;

        private:
            Kind   m_kind  = Kind::isotropic;
            Faces  m_faces = Faces::front;
            double m_n     = 0.0;
            std::vector<double> m_theta_deg;   ///< table: zenith angles (deg)
            std::vector<double> m_f;           ///< table: relative intensities
        };
#pragma endregion

#pragma region - geometric factor (no shielding)
        /**
         * Geometry of a flat sensor in the sky frame, nothing shielding it. The sensor normal (front face) is tilted
         * by tilt_deg from the zenith towards +x of the matrix; a sky direction is given by its zenith angle
         * thetaZ_deg (0: from above, 180: from below) and its azimuth phiZ_deg (0: towards +x of the matrix).
         * The incidence angle theta is the angle to the sensor normal, mu = |cos(theta)|.
         */

        /// @brief Incidence of a particle on the sensor.
        struct Incidence
        {
            double mu;     ///< |cos(theta)|, theta the angle to the sensor normal
            double phi;    ///< Azimuth of the track on the matrix (rad), in [0, 2 pi)
            int    face;   ///< Face reached: +1 front, -1 back
        };

        /// @brief Incidence of a particle coming from the sky direction (thetaZ_deg, phiZ_deg), its track running
        ///        away from that direction.
        /// @throws std::invalid_argument if tilt_deg or thetaZ_deg is not in [0, 180].
        Incidence incidence(double thetaZ_deg, double phiZ_deg, double tilt_deg = 0.0);

        /// @brief Area (m2) the sensor presents to a parallel beam from (thetaZ_deg, phiZ_deg): A |cos(theta)|.
        /// @throws std::invalid_argument if area_m2 <= 0, or see incidence().
        double effectiveArea(double area_m2, double thetaZ_deg, double phiZ_deg, double tilt_deg = 0.0);

        /// @brief Geometric factor G (sr) of an angular model on the sensor: ∫ I / J mu dOmega over the faces of the
        ///        model and the incidence angles theta in [thetaMin_deg, thetaMax_deg] (full sky by default).
        /// @details Closed forms: isotropic pi (cos^2 thetaMin - cos^2 thetaMax) per face; surface
        ///          2 pi (cos thetaMin - cos thetaMax) per face; cos^n on a horizontal sensor (front face only)
        ///          2 pi / (n + 2) (cos^(n+2) thetaMin - cos^(n+2) thetaMax). Otherwise numerical (midpoint rule,
        ///          400 x 720 sky directions, about 1e-3 relative).
        /// @throws std::invalid_argument if tilt_deg is not in [0, 180] or not 0 <= thetaMin < thetaMax <= 90.
        double geometricFactor(const AngularModel& model, double tilt_deg = 0.0,
                               double thetaMin_deg = 0.0, double thetaMax_deg = 90.0);

        /// @brief Etendue A G (m2 sr) of a sensor of area A for an angular model (see geometricFactor()): hits per
        ///        second = etendue * ∫ J dE.
        /// @throws std::invalid_argument if area_m2 <= 0, or see geometricFactor().
        double etendue(const AngularModel& model, double area_m2, double tilt_deg = 0.0,
                       double thetaMin_deg = 0.0, double thetaMax_deg = 90.0);
#pragma endregion

#pragma region - spectrum
        /**
         * @class CosmicSpectrum
         * @brief Differential intensity of the cosmic rays, per species: probabilities of a hit.
         * @details Each species is a table (E, J), E the kinetic energy (ascending), per nucleon for a nucleus
         *          (MeV/nucleon), per particle otherwise (MeV), and J > 0 the differential intensity
         *          (m-2 s-1 sr-1 per unit of E).
         *          J is taken to be a power law between two bin, both the bin's integral and the 
         *          sampling of energies have closed formulas. Nothing is integrated numerically
         *          or tabulated.
         *          The cells of every species, together, are the outcomes of a hit.
         */
        class CosmicSpectrum
        {
        public:
            /// @brief One interval between two table points of a species.
            struct EnergyBin
            {
                size_t species;   ///< Index of the species (see species())
                double E0, E1;    ///< Energy range (MeV/nucleon for a nucleus, MeV otherwise)
                double J0, J1;    ///< Differential intensity at E0 and E1 (m-2 s-1 sr-1 per unit of E)
                double flux;      ///< ∫ J dE over the cell (m-2 s-1 sr-1)
            };

            /// @brief A hit: species and kinetic energy.
            struct Sample
            {
                size_t species;   ///< Index of the species
                double E;         ///< Kinetic energy (MeV/nucleon for a nucleus, MeV otherwise)
            };

            /// @brief Lateral spread of the charge of a species: a Gaussian core and an optional Gaussian halo.
            /// @details Share of the charge of a pixel of the track in its neighbour (a, b):
            ///          (1 - haloFraction) g_core(a) g_core(b) + haloFraction g_halo(a) g_halo(b), g being a
            ///          Gaussian integrated over the pixel and cut at 3 sigma.
            struct Footprint
            {
                double sigma_px     = std::numeric_limits<double>::quiet_NaN();   ///< Core sigma (pixel), >= 0; NaN: sigmaDiffusion_px of the source
                double haloSigma_px = 0.0;   ///< Halo sigma (pixel), > 0 when haloFraction > 0
                double haloFraction = 0.0;   ///< Share of the charge in the halo, [0, 1]; 0: no halo
                bool   deltaRays    = false; ///< Delta-ray halo of each hit, from its velocity (heavy particles only)
            };

            /// @brief Add a species, with the lateral spread of the source and no halo (Footprint()).
            /// @details An overload rather than a default argument: Footprint() is not usable in a default argument
            ///          inside CosmicSpectrum (its member initializers are needed before the class is complete).
            void add(const specie& species, const std::vector<double>& E, const std::vector<double>& J,
                     const AngularModel& angular = AngularModel());

            /// @brief Add a species.
            /// @param species   Species.
            /// @param E         Energies (MeV/nucleon for a nucleus, MeV otherwise), strictly ascending, at least 2.
            /// @param J         Differential intensity at each energy (m-2 s-1 sr-1 per unit of E), > 0.
            /// @param angular   Angular distribution of the species and faces it reaches.
            /// @param footprint Lateral spread of its charge.
            /// @throws std::invalid_argument on invalid values; the spectrum is then unchanged.
            void add(const specie& species, const std::vector<double>& E, const std::vector<double>& J,
                     const AngularModel& angular, const Footprint& footprint);

            /// @brief Spectrum from the "spectrum" array of cfg: one entry per species, its ENERGY and FLUX read
            ///        from a CSV file or a FITS table, filtered on filterColumn == filterValue (see CosmicDeposit.h).
            /// @param cfg Configuration holding the "spectrum" key (the "cosmic" object).
            /// @throws std::out_of_range if cfg has no "spectrum"; std::invalid_argument on a malformed entry,
            ///         a file that cannot be read, a table without ENERGY and FLUX, no row for the species,
            ///         or values add() rejects.
            static CosmicSpectrum load(const PixelConfig& cfg);

            /// @brief True without any species.
            bool empty() const { return m_bins.empty(); }

            /// @brief Species.
            const std::vector<specie>& species() const { return m_species; }

            /// @brief Angular distribution of each species.
            const std::vector<AngularModel>& angular() const { return m_angular; }

            /// @brief Lateral spread of each species.
            const std::vector<Footprint>& footprints() const { return m_footprint; }

            /// @brief Bin of every species.
            const std::vector<EnergyBin>& bins() const { return m_bins; }

            /// @brief Energy (MeV/nucleon for a nucleus, MeV otherwise) with cumulated probability u in [0, 1] in a bin.
            /// @throws std::out_of_range if there is no such bin.
            double energyInBin(size_t bin, double u) const;

            /// @brief ∫ J dE over every species (m-2 s-1 sr-1).
            double integralIntensity() const { return m_total; }
            /// @brief ∫ J dE of one species (m-2 s-1 sr-1).
            double integralIntensity(size_t species) const;

            /// @brief Probability that a hit falls in each cell: flux / integralIntensity(). Sums to 1.
            std::vector<double> binProbabilities() const;
            /// @brief Probability that a hit is of each species. Sums to 1.
            std::vector<double> speciesProbabilities() const;

            /// @brief Species and energy of a hit, from two uniform numbers in [0, 1).
            /// @param uBin    Chooses the bin, with probability binProbabilities().
            /// @param uEnergy Chooses the energy in the cell, following its power law.
            /// @pre not empty().
            Sample sample(double uBin, double uEnergy) const;

        private:
            std::vector<specie>        m_species;        ///< Species
            std::vector<AngularModel>  m_angular;        ///< Angular distribution of each species
            std::vector<Footprint>     m_footprint;      ///< Lateral spread of each species
            std::vector<EnergyBin>     m_bins;           ///< Cells of every species
            std::vector<double>        m_cdf;            ///< Cumulated flux up to the end of each cell
            double                     m_total = 0.0;    ///< Total flux
        };
#pragma endregion

#pragma region - stopping power and range
        /**
         * @class StoppingTable
         * @brief Continuous slowing down approximation (CSDA) range of one species in one material, tabulated from its mean stopping power.
         * @details T is the kinetic energy of the particle (MeV), S the mean stopping power (MeV cm2/g), R the
         *          range (g/cm2) in the continuous slowing down approximation continuous slowing down approximation (CSDA).
         *          - Below Tc, \$ S(T) = S(Tc) \sqrt{\frac{T}{Tc}}\$: \$ R = 2 \frac{\sqrt{T Tc}{S(Tc)} \$, analytical.
         *          - Above Tc, R(T) = R(Tc) + ∫ dT / S on a grid in ln T (perDecade nodes per decade, Simpson on
         *            each interval), up to Tmax; above Tmax, extrapolated.
         *          - ln R is interpolated linearly in ln T, so energy() is the exact inverse of range().
         *          - Residual energy after a path s: E(s) = R^-1(R(T) - s) (energyAfter()).
         *          Known limitation: the low-energy part (below Tc) is a velocity-proportional stopping power
         *          joined to the formula given; it ignores shell corrections and the effective charge of ions.
         */
        class StoppingTable
        {
        public:
            /// @brief Empty table (empty() is true).
            StoppingTable() = default;

            /// @brief Tabulate the range from a stopping power.
            /// @param S         Mean stopping power (MeV cm2/g) of the kinetic energy (MeV); finite and > 0 from Tc to Tmax.
            /// @param Tc        Kinetic energy (MeV) below which S = S(Tc) sqrt(T / Tc), > 0.
            /// @param Tmax      Highest kinetic energy of the table (MeV).
            /// @param perDecade Nodes per decade of T, >= 1.
            /// @throws std::invalid_argument if Tc <= 0, perDecade = 0, or S is not finite and > 0 on the grid.
            StoppingTable(const std::function<double(double)>& S, double Tc, double Tmax, size_t perDecade = 100);

            /// @brief True without a table.
            bool empty() const { return !(m_Sc > 0.0); }

            /// @brief Kinetic energy (MeV) below which S is proportional to sqrt(T).
            double lowEnergyLimit() const { return m_Tc; }

            /// @brief CSDA range (g/cm2) of a particle of kinetic energy T (MeV); 0 if T <= 0.
            double range(double T) const;

            /// @brief Kinetic energy (MeV) of a particle whose CSDA range is R (g/cm2): inverse of range(); 0 if R <= 0.
            double energy(double R) const;

            /// @brief Kinetic energy (MeV) left after a path s (g/cm2) from T; 0 once the particle has stopped.
            double energyAfter(double T, double s) const
            {
                const double R = range(T);
                return (s >= R) ? 0.0 : energy(R - s);
            }

        private:
            double m_Tc = 0.0;            ///< End of the sqrt(T) part (MeV)
            double m_Sc = 0.0;            ///< S(Tc) (MeV cm2/g)
            double m_Rc = 0.0;            ///< R(Tc) = 2 Tc / S(Tc) (g/cm2)
            std::vector<double> m_lnT;    ///< ln T of the nodes, from ln Tc
            std::vector<double> m_lnR;    ///< ln R of the nodes
        };
#pragma endregion

#pragma region - cosmic deposit

        /**
         * @class CosmicDeposit
         * @brief Cosmic ray hits drawn from a CosmicSpectrum (see CosmicDeposit.h for the model).
         * @details The probabilities (expectedEvents(), probabilityAtLeastOne(), geometricFactor(), truncatedFraction())
         *          are available without simulating. A derived class adapts the model to one instrument:
         *          it sets defaults in its constructor (m_params), and may override energyDeposit_eV().
         */
        class CosmicDeposit : public ChargeDeposit
        {
        public:
            /// @brief Parameters of the model; NaN: not set yet (see CosmicDeposit.h).
            struct Params
            {
                double depth_um          = std::numeric_limits<double>::quiet_NaN();   ///< Sensitive depth d (um)
                double pitch_um          = std::numeric_limits<double>::quiet_NaN();   ///< Pixel pitch (um)
                double eVPerElectron     = std::numeric_limits<double>::quiet_NaN();   ///< Energy per pair (eV)
                double tilt_deg          = 0.0;                ///< Angle of the sensor normal (front face) to the zenith (deg), [0, 180]
                double sigmaDiffusion_px = 0.0;                ///< Lateral spread (pixel); 0: none
                double maxTrackPixels    = 64.0;               ///< Longest projected track drawn (pixel)
                double deltaRayCut_e     = 1.0;                ///< Smallest charge a delta-ray halo puts on a pixel (e-)
                target::material material;                     ///< Material of the sensitive depth
            };

            /// @brief One hit, in physical detector coordinates (pixel i covers [i, i+1)).
            struct Event
            {
                size_t species;     ///< Index in spectrum().species()
                double energy;      ///< Kinetic energy (MeV/nucleon for a nucleus, MeV otherwise)
                double mu;          ///< cos(theta) of the incidence, theta from the normal of the sensor
                double phi;         ///< Azimuth of the track on the matrix (rad)
                int    face;        ///< Face crossed at entry: +1 front, -1 back
                double x, y;        ///< Entry point of the track
                double dx, dy;      ///< Projected displacement of the path in the layer (pixel)
                double path_um;     ///< Path in the layer (um): depth / mu, or the range if the particle stops; 0: uniform charge
                double electrons;   ///< Charge deposited in the layer (e-)
            };

            /// @brief Most probable energy loss and width of the Landau distribution (eV).
            struct EnergyLoss { double mpv_eV; double xi_eV; };

            /// @brief Name, unit and meaning of one column of the table of hits (see toColumns()).
            struct ColumnInfo
            {
                std::string name;      ///< Column name, at most 8 characters (FITS TTYPE)
                std::string unit;      ///< Unit (FITS TUNIT); empty: none
                std::string comment;   ///< Meaning
            };

            /// @brief Hits as columns, one value per hit and per column, ready to be written as a table.
            /// @details Same order as columnInfo(). x and y are in physical detector coordinates with pixel i
            ///          covering [i, i+1): in the FITS convention (first pixel centred on 1.0) the position of
            ///          a hit is x + 0.5, y + 0.5.
            struct EventColumns
            {
                std::vector<std::string> species;   ///< SPECIES: name of the species
                std::vector<double> Z;              ///< Z: charge number
                std::vector<double> A;              ///< A: mass number
                std::vector<double> energy;         ///< ENERGY: kinetic energy (MeV/nucleon for a nucleus, MeV otherwise)
                std::vector<double> mu;             ///< MU: cos(theta) of the incidence
                std::vector<double> phi;            ///< PHI: azimuth of the track (rad)
                std::vector<double> face;           ///< FACE: face crossed at entry, +1 front, -1 back
                std::vector<double> x, y;           ///< X, Y: entry point (pixel)
                std::vector<double> dx, dy;         ///< DX, DY: projected displacement of the path in the layer (pixel)
                std::vector<double> path;           ///< PATH: path in the layer (um)
                std::vector<double> electrons;      ///< NE: charge deposited in the layer (e-)

                /// @brief Number of hits.
                size_t size() const { return electrons.size(); }
            };

            /// @brief Not configured: generate() throws until configure() succeeds.
            CosmicDeposit() = default;
            /// @brief From a spectrum and its parameters. @throws std::invalid_argument see load().
            CosmicDeposit(CosmicSpectrum spectrum, const Params& params);

            /// @brief Copy of this source with its dynamic type.
            std::shared_ptr<ChargeDeposit> clone() const override { return std::make_shared<CosmicDeposit>(*this); }

            /// @brief Read the "cosmic" key of cfg with load(); nothing if it is absent.
            void configure(const PixelConfig& cfg) override;

            /// @brief True when the spectrum and every required parameter are set.
            bool isConfigured() const override;

            /// @brief Parameters.
            const Params& params() const { return m_params; }
            /// @brief Spectrum.
            const CosmicSpectrum& spectrum() const { return m_spectrum; }

            // ── Probabilities ────────────────────────────────────────────────────

            /// @brief Hits per second through an area: sum over the species of G * ∫ J dE * area (s-1).
            double eventRate(double area_m2) const;
            /// @brief Mean number of hits in dt seconds on an area (m2), every track length included.
            double expectedEvents(double area_m2, double dt) const;
            /// @brief P(at least one hit) = 1 - exp(-expectedEvents()).
            double probabilityAtLeastOne(double area_m2, double dt) const;
            /// @brief P(n hits), Poisson with mean expectedEvents().
            double probability(size_t n, double area_m2, double dt) const;
            /// @brief Share of the hits whose track is longer than maxTrackPixels, so not drawn (all species).
            double truncatedFraction() const;
            /// @brief Same, for species i.
            /// @throws std::out_of_range if i is not a species of the spectrum.
            double truncatedFraction(size_t i) const { return m_acceptance.at(i).truncated; }
            /// @brief Geometric factor G of species i (sr): ∫ I / J * mu dOmega over the faces it reaches.
            /// @throws std::out_of_range if i is not a species of the spectrum.
            double geometricFactor(size_t i) const { return m_acceptance.at(i).G; }

            // ── Generation ───────────────────────────────────────────────────────

            /// @brief Hits of a step: entry points uniform on the physical detector.
            /// @param dt       Duration of the step (s).
            /// @param geometry Physical detector: no hit enters from outside it.
            /// @param rng      Random stream.
            std::vector<Event> drawEvents(double dt, const Geometry& geometry, std::mt19937_64& rng) const;

            /// @brief Append the deposits of one hit: in each pixel crossed, the CSDA energy lost on the part of the
            ///        path in it (uniform when event.path_um is 0), then the lateral spread of its species
            ///        (CosmicSpectrum::Footprint; sigmaDiffusion_px without halo for a species outside the spectrum).
            void depositEvent(const Event& event, const Geometry& geometry, std::vector<Deposit>& out) const;

            void generate(double dt, const Geometry& geometry, std::vector<Deposit>& out) const override;

            /// @brief generate() that also returns the hits drawn.
            /// @details The hits are the whole record of the step: depositEvent() is deterministic, so the
            ///          deposits are given again by the hits and the configuration. Hits are appended to
            ///          events, so one list can collect the steps of an exposure.
            /// @param dt       Duration of the step (s).
            /// @param geometry Physical detector.
            /// @param out      Deposits are appended to it.
            /// @param events   Hits are appended to it; null: not returned (as generate()).
            void generate(double dt, const Geometry& geometry, std::vector<Deposit>& out,
                          std::vector<Event>* events) const;

            /// @brief Columns of the table of hits: name, unit and meaning, in the order of EventColumns.
            static const std::vector<ColumnInfo>& columnInfo();

            /// @brief Hits as columns, with the species of this source (see EventColumns).
            EventColumns toColumns(const std::vector<Event>& events) const;

            /// @brief Header cards of the inputs (CRFLUX, CRACC, ..., see the wiki, CosmicDeposit §5.1).
            std::vector<HeaderCard> inputCards() const override;

            // ── Energy loss ──────────────────────────────────────────────────────

            /// @brief Most probable loss and width of the Landau distribution for a path in the material.
            /// @param species    Species.
            /// @param T          Kinetic energy (MeV/nucleon for a nucleus, MeV otherwise).
            /// @param path_um    Path in the material (um).
            EnergyLoss energyLoss(const specie& species, double T, double path_um) const;

            /// @brief Mean collision stopping power (MeV cm2/g) in the material (see collisionStopping()).
            /// @details Below the energy where the formula is largest (StoppingTable::lowEnergyLimit()), the
            ///          stopping power is proportional to the velocity.
            /// @param species Species.
            /// @param E       Kinetic energy (MeV/nucleon for a nucleus, MeV otherwise).
            double stoppingPower(const specie& species, double E) const;

            /// @brief Range table of species i of the spectrum (energies: kinetic energy of the particle, MeV).
            /// @throws std::out_of_range if i is not a species of the spectrum.
            const StoppingTable& stoppingTable(size_t i) const { return m_tables.at(i); }

            /// @brief Share h of the charge of a hit of species i in its delta-ray halo (see CosmicDeposit.h, 5.).
            /// @param i Species of the spectrum.
            /// @param E Kinetic energy at entry (MeV/nucleon for a nucleus, MeV otherwise).
            /// @return h in [0, 1); 0 for a species without "deltaRays" or when Wmax <= T0.
            /// @throws std::out_of_range if i is not a species of the spectrum.
            double deltaRayShare(size_t i, double E) const;

        protected:
            /// @brief Load the parameters and the spectrum from the "cosmic" configuration (see CosmicDeposit.h).
            /// @details Absent keys leave the value unchanged. Everything is checked on a copy first.
            ///          Protected: cfg is the content of the "cosmic" object, not a detector configuration;
            ///          from outside, configure() reads that object.
            /// @param cfg Configuration of the "cosmic" object.
            /// @throws std::invalid_argument if a value is not valid, or a required one is still missing;
            ///         the source is then unchanged.
            void load(const PixelConfig& cfg);

            /// @brief Rebuild what is derived from m_params and m_spectrum: range tables, geometric factors,
            ///        truncated shares, the weights of the draws and the lateral spread (one entry per species).
            /// @details The constructor and load() call it. A derived class that changes m_params or
            ///          m_spectrum directly calls it afterwards.
            /// @throws std::invalid_argument if a stopping power is not > 0 over a spectrum.
            void buildTables();

            /// @brief Energy deposited by a particle crossing a thin layer (eV): Landau, at most its kinetic energy.
            /// @details Used for the thin crossings only: drawEvents() handles a particle that stops or loses
            ///          more than 10 % of its energy. Override it to use another distribution (measured, Vavilov, ...).
            /// @param species Species.
            /// @param T       Kinetic energy (MeV/nucleon for a nucleus, MeV otherwise).
            /// @param path_um Path in the material (um).
            /// @param u       Uniform number in [0, 1).
            virtual double energyDeposit_eV(const specie& species, double T, double path_um, double u) const;

            Params         m_params;     ///< Parameters; a derived class sets its defaults in its constructor
            CosmicSpectrum m_spectrum;   ///< Spectrum

        private:
            /// @brief Density effect delta(beta*gamma), Sternheimer.
            double densityEffect(double betaGamma) const;
            /// @brief Mean collision stopping power (MeV cm2/g) of a particle of kinetic energy T (MeV), from the
            ///        formula only: Bethe-Bloch for heavy particles, collision formula for e- and e+.
            double collisionStopping(const specie& species, double T) const;
            /// @brief Direction of a hit of species i, in the frame of the sensor: mu, phi and face, grazing
            ///        tracks (mu below the truncation) excluded.
            void drawDirection(size_t i, std::mt19937_64& rng, double& mu, double& phi, int& face) const;

            /// @brief Geometric factor and truncation of one species.
            struct Acceptance
            {
                double G         = 0.0;   ///< Geometric factor (sr)
                double truncated = 0.0;   ///< Share of the flux through the sensor below muMin
                double muMin     = 0.0;   ///< Smallest mu drawn (maxTrackPixels)
            };

            /// @brief Lateral spread around a pixel: (2R+1)^2 shares of its charge, sum 1.
            struct Kernel
            {
                Kernel() = default;
                /// @brief (1 - haloFraction) core + haloFraction halo, each a Gaussian integrated over the pixel,
                ///        cut at 3 sigma; a sigma of 0 is the pixel alone.
                Kernel(double sigma_px, double haloSigma_px, double haloFraction);

                long R = 0;                  ///< Half width (pixel)
                std::vector<double> w{1.0};  ///< Share in pixel (a, b), a and b in [-R, R], at (b + R) (2R + 1) + a + R
            };

            /// @brief Delta-ray halo of a hit (see CosmicDeposit.h, 5.); radii in pixels.
            struct DeltaRays
            {
                double share = 0.0;   ///< h: share of the charge of the hit in the halo; 0: no halo
                double r0    = 0.0;   ///< Inner radius, 1/2 pixel: range of an electron of energy T0
                double rMax  = 0.0;   ///< CSDA range of an electron of energy Wmax
                double d     = 0.0;   ///< Depth of the layer: beyond it, a delta ray keeps d / r of its energy
                double a     = 0.0;   ///< d within [r0, rMax]
                double norm  = 0.0;   ///< ∫ min(1, d / r) dln r from r0 to rMax
            };

            /// @brief Delta-ray halo of a hit of species i at energy E; none for a species without "deltaRays".
            DeltaRays deltaRays(size_t i, double E) const;
            /// @brief Kernel of a delta-ray halo holding `electrons` around one pixel of the track, cut where it puts
            ///        less than deltaRayCut_e on a pixel; R = 0 (no halo) when that is everywhere.
            Kernel deltaRayKernel(const DeltaRays& delta, double electrons) const;
            /// @brief Deposit electrons in pixel (X, Y) and its neighbours, shared by kernel.
            void spread(long X, long Y, double electrons, const Kernel& kernel, const Geometry& geometry, std::vector<Deposit>& out) const;

            Kernel                     m_kernel;       ///< Lateral spread of sigmaDiffusion_px, for a species outside the spectrum
            std::vector<Kernel>        m_kernels;      ///< Lateral spread of each species of the spectrum
            std::vector<StoppingTable> m_tables;       ///< Range table of each species of the spectrum
            std::vector<Acceptance>    m_acceptance;   ///< Geometric factor of each species of the spectrum
            std::vector<double>        m_drawCdf;      ///< Cumulated flux * G * (1 - truncated) of the bins (m-2 s-1)
            StoppingTable              m_electronTable;   ///< CSDA range of the delta rays; empty without "deltaRays"
            double                     m_deltaT0 = 0.0;   ///< T0 (MeV): energy of an electron whose CSDA range is 1/2 pixel
        };
    }

#pragma endregion

} // namespace phot2pix

#endif // P2P_COSMIC_DEPOSIT_H