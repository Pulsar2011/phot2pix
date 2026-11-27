// Created by gillard on 2026-11-25.
// Written with the assistance of Claude Opus 5.5 (Anthropic's AI model) under guidance of William Gillard.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
// You may share and adapt this code with attribution, Non Comercial use only.
//
// CosmicDeposit.cxx

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>

#include <DSTfits/FITSmanager.h>      // DSL::FITSmanager, DSL::FITStable
#include <DSTmath/DSTstatistic.h>     // DST::Math::standardLandau

#include "CosmicDeposit.h"

namespace phot2pix
{
#pragma region Cosmic ray 
    namespace CosmicRay
    {

        constexpr double kPi              = std::numbers::pi;
        constexpr double kElectronMass_MeV = 0.510999;  // m_e c^2
        constexpr double kK_MeV_cm2_mol   = 0.307075;   // 4 pi N_A r_e^2 m_e c^2
        constexpr double kLandauJ         = 0.200;      // Bichsel, Landau most probable value
        constexpr double kThickShare      = 0.1;        // thick crossing: CSDA mean loss above this share of T

#pragma region target
        namespace target
        {
#pragma region -- building compound materials
            material material::compound(const std::vector<element>& elements, double density_gcm3)
            {
                if(elements.empty() || !(density_gcm3 > 0.0))
                    throw std::invalid_argument("[target::material::compound] No element, or density <= 0.");

                double nZ = 0.0, nA = 0.0, nZlnI = 0.0;
                for(const element& e : elements)
                {
                    if(!(e.Z > 0.0 && e.A > 0.0 && e.I_eV > 0.0 && e.n >= 0.0))
                        throw std::invalid_argument("[target::material::compound] Z, A, I must be > 0 and n >= 0.");
                    nZ    += e.n * e.Z;
                    nA    += e.n * e.A;
                    nZlnI += e.n * e.Z * std::log(e.I_eV);
                }
                if(!(nZ > 0.0))
                    throw std::invalid_argument("[target::material::compound] No atom.");

                material m;
                m.ZoverA            = nZ / nA;
                m.meanExcitation_eV = std::exp(nZlnI / nZ);
                m.density_gcm3      = density_gcm3;

                const double plasma_eV = 28.816 * std::sqrt(m.density_gcm3 * m.ZoverA);   // hbar omega_p
                m.sternheimerC = 2.0 * std::log(m.meanExcitation_eV / plasma_eV) + 1.0;

                // Sternheimer & Peierls, solids and liquids
                const bool lowI = m.meanExcitation_eV < 100.0;
                m.sternheimerX1 = lowI ? 2.0 : 3.0;
                m.sternheimerX0 = (m.sternheimerC < (lowI ? 3.681 : 5.215)) ? 0.2 : 0.326 * m.sternheimerC - (lowI ? 1.0 : 1.5);
                m.sternheimerM  = 3.0;
                m.sternheimerA  = (m.sternheimerC - 2.0 * std::numbers::ln10 * m.sternheimerX0)
                                / std::pow(m.sternheimerX1 - m.sternheimerX0, m.sternheimerM);
                m.sternheimerDelta0   = 0.0;
                m.sternheimerDeltaMax = std::numeric_limits<double>::quiet_NaN();

                return m;
            }
#pragma endregion

#pragma region -- HgCdTe
            /**
             * @brief Hg1−xCdxTe, from the Cd fraction x.
             * @details The default values of the material are for x = 0.3,
             *       which gives a density of 7.63 g/cm3 and a band gap of 0.1 eV at 77 K.
            */ 
            material material::HgCdTe(const double& x)
            {
                if(!(x >= 0.0 && x <= 1.0))
                    throw std::invalid_argument("[target::material::HgCdTe] The Cd fraction x must be in [0, 1].");

                constexpr double bound = 1.13;   // I of an atom bound in a condensed compound / I of the element
                const std::vector<element> elements = {
                    {80.0, 200.592, bound * 800.0, 1.0 - x},   // Hg
                    {48.0, 112.414, bound * 469.0, x      },   // Cd
                    {52.0, 127.60 , bound * 485.0, 1.0    }    // Te
                };

                double M = 0.0;   // g/mol of formula unit
                for(const element& e : elements)
                    M += e.n * e.A;

                // Lattice constant, room temperature (Higgins et al., J. Vac. Sci. Technol. A 7 (1989) 271):
                //   a(x) = 6.4614 + 0.0084 x + 0.0168 x^2 - 0.0057 x^3 (A). Thermal contraction at ~100 K (~0.1 % in a) neglected.
                const double a = (6.4614 + x * (0.0084 + x * (0.0168 - 0.0057 * x))) * 1e-8;   // cm

                return compound(elements, 4.0 * M / (N_A * a * a * a));
            }
        }   // namespace CosmicRay::target
#pragma endregion
#pragma endregion

#pragma region Helpers
        /// J = J0 (E/E0)^-g between (E0, J0) and (E1, J1).
        double powerLawIndex(double E0, double E1, double J0, double J1)
        {
            return -std::log(J1 / J0) / std::log(E1 / E0);
        }

        /// ∫ J dE over [E0, E1] of the power law.
        double powerLawIntegral(double E0, double E1, double J0, double J1)
        {
            const double r = E1 / E0;
            const double g = powerLawIndex(E0, E1, J0, J1);
            if (std::abs(1.0 - g) < 1.0e-9)
                return J0 * E0 * std::log(r);
            return J0 * E0 * (std::pow(r, 1.0 - g) - 1.0) / (1.0 - g);
        }

        /// Energy with cumulated probability u in [E0, E1] of the power law.
        double powerLawQuantile(double E0, double E1, double J0, double J1, double u)
        {
            const double r = E1 / E0;
            const double g = powerLawIndex(E0, E1, J0, J1);
            if (std::abs(1.0 - g) < 1.0e-9)
                return E0 * std::pow(r, u);
            return E0 * std::pow(1.0 + u * (std::pow(r, 1.0 - g) - 1.0), 1.0 / (1.0 - g));
        }
#pragma endregion

#pragma region -- spectrum files
        /// Table of the file of a species: a CSV file (".csv") through DSL::FITStable::fromCSV(),
        /// else the table in HDU "hdu" (required) of a FITS file.
        std::shared_ptr<DSL::FITStable> spectrumTable(const PixelConfig& entry, const std::string& path)
        {
            std::string ext = std::filesystem::path(path).extension().string();
            std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

            if (ext == ".csv")
            {
                if (entry.has("hdu"))
                    throw std::invalid_argument("\"hdu\" goes with a FITS file, not with a CSV file.");
                
                std::ifstream file(path);
                if (!file)
                    throw std::invalid_argument("cannot open the file.");
                
                return std::make_shared<DSL::FITStable>(DSL::FITStable::fromCSV(file));
            }

            DSL::FITSmanager mgr(path, true);   // read only; the table keeps its data once the file is closed
            
            if (!entry.has("hdu"))
                throw std::invalid_argument("\"hdu\" index or name should be provided.");

            double index;
            try
            { 
                index = entry.getAsDouble("hdu");
            }
            catch (const std::invalid_argument&)
            {
                return mgr.GetTable(entry.get<std::string>("hdu"));
            }   // EXTNAME
            
            if (!(index >= 2.0) || index != std::floor(index))
                throw std::invalid_argument("\"hdu\" is the index of a table HDU (2 = first extension) or its EXTNAME.");
            return mgr.GetTableAtIndex(static_cast<int>(index));
        }

        /// A real column (FITS E or D) of a table as doubles, on some rows.
        std::vector<double> columnAsDouble(const DSL::FITStable& table, const std::string& name, const DSL::RowSet& rows)
        {
            const DSL::dtype type = table.getColumn(name)->getType();   // std::out_of_range if absent
            if (type == DSL::tdouble)
                return table.column<double>(name).on(rows).data();
            if (type == DSL::tfloat)
            {
                const std::vector<float> v = table.column<float>(name).on(rows).data();
                return std::vector<double>(v.begin(), v.end());
            }
            throw std::invalid_argument("the column " + name + " must hold real numbers (FITS E or D).");
        }

        /// DSL::FITStable::where() expression selecting filterColumn == filterValue; empty without them (every row).
        std::string selection(const PixelConfig& entry)
        {
            if (!entry.has("filterColumn") && !entry.has("filterValue"))
                return std::string();

            if (!entry.has("filterColumn") || !entry.has("filterValue"))
                throw std::invalid_argument("give both \"filterColumn\" and \"filterValue\", or neither.");

            const json v = entry.toJson().at("filterValue");
            std::string value;
            
            if (v.is_string())
            {
                const std::string s = v.get<std::string>();
                const char q = (s.find('"') == std::string::npos) ? '"' : '\'';
                if (s.find(q) != std::string::npos)
                    throw std::invalid_argument("\"filterValue\" cannot hold both kinds of quote.");
                value = q + s + q;
            }
            else if (v.is_number())
                value = v.dump();
            else
                throw std::invalid_argument("\"filterValue\" must be a string or a number.");

            return entry.get<std::string>("filterColumn") + " == " + value;
        }

        /// Two real columns of the table of cfg ("file", "hdu"), named by the values of the keys keyX and keyY,
        /// on the rows where filterColumn == filterValue.
        void readColumns(const PixelConfig& cfg, const std::string& keyX, const std::string& keyY,
                         std::vector<double>& X, std::vector<double>& Y)
        {
            const std::string path = cfg.get<std::string>("file");
            try
            {
                const std::shared_ptr<DSL::FITStable> table = spectrumTable(cfg, path);
                const std::string where = selection(cfg);

                const DSL::RowSet rows = table->where(where);
                if (rows.empty())
                    throw std::invalid_argument(where.empty() ? std::string("the table has no row.") : "no row where " + where + ".");

                const std::string colX = cfg.get<std::string>(keyX);
                const std::string colY = cfg.get<std::string>(keyY);
                try
                {
                    X = columnAsDouble(*table, colX, rows);
                    Y = columnAsDouble(*table, colY, rows);
                }
                catch (const std::out_of_range&)
                {
                    throw std::invalid_argument("the table needs the columns " + colX + " and " + colY + ".");
                }
            }
            catch (const std::exception& ex)
            {
                throw std::invalid_argument("file '" + path + "': " + ex.what());
            }
        }

        /// Energy (colEnergy) and flux (colFlux) of a species, on the rows of its table where filterColumn == filterValue.
        void readSpectrum(const PixelConfig& entry, std::vector<double>& E, std::vector<double>& J)
        {
            readColumns(entry, "colEnergy", "colFlux", E, J);
        }

        /// Lateral spread of a species: "sigmaDiffusion_px" (core), "halo": {"sigma_px", "fraction"} and "deltaRays"
        /// (a boolean), all optional; CosmicSpectrum::add() checks the values.
        CosmicSpectrum::Footprint readFootprint(const PixelConfig& entry)
        {
            CosmicSpectrum::Footprint f;
            if (entry.has("sigmaDiffusion_px"))
                f.sigma_px = entry.getAsDouble("sigmaDiffusion_px");

            if (entry.has("deltaRays"))
            {
                try                           { f.deltaRays = entry.get<bool>("deltaRays"); }
                catch (const std::exception&) { throw std::invalid_argument("\"deltaRays\" must be true or false."); }
            }

            if (entry.has("halo"))
            {
                const json& h = entry.getJson("halo");
                if (!h.is_object())
                    throw std::invalid_argument("\"halo\" is an object {\"sigma_px\", \"fraction\"}.");
                const PixelConfig halo(h);
                if (!halo.has("sigma_px") || !halo.has("fraction"))
                    throw std::invalid_argument("\"halo\" needs \"sigma_px\" and \"fraction\".");
                f.haloSigma_px = halo.getAsDouble("sigma_px");
                f.haloFraction = halo.getAsDouble("fraction");
            }
            return f;
        }

#pragma endregion

#pragma region - AngularModel

    AngularModel AngularModel::isotropic(Faces faces)
    {
        AngularModel m;
        m.m_kind  = Kind::isotropic;
        m.m_faces = faces;
        return m;
    }

    AngularModel AngularModel::surface(Faces faces)
    {
        AngularModel m;
        m.m_kind  = Kind::surface;
        m.m_faces = faces;
        return m;
    }

    AngularModel AngularModel::cosn(double n, Faces faces)
    {
        if (!(std::isfinite(n) && n >= 0.0))
            throw std::invalid_argument("AngularModel: the exponent n of cosn must be finite and >= 0.");
        AngularModel m;
        m.m_kind  = Kind::cosn;
        m.m_faces = faces;
        m.m_n     = n;
        return m;
    }

    AngularModel AngularModel::table(std::vector<double> theta_deg, std::vector<double> f, Faces faces)
    {
        if (theta_deg.size() < 2 || theta_deg.size() != f.size())
            throw std::invalid_argument("AngularModel: a table needs at least 2 zenith angles, and as many intensities.");

        bool positive = false;
        for (size_t i = 0; i < theta_deg.size(); ++i)
        {
            if (!std::isfinite(theta_deg[i]) || theta_deg[i] < 0.0 || theta_deg[i] > 90.0 || (i > 0 && !(theta_deg[i] > theta_deg[i - 1])))
                throw std::invalid_argument("AngularModel: the zenith angles must be strictly ascending in [0, 90] deg.");
            if (!std::isfinite(f[i]) || f[i] < 0.0)
                throw std::invalid_argument("AngularModel: the relative intensities must be finite and >= 0.");
            positive = positive || f[i] > 0.0;
        }
        if (!positive)
            throw std::invalid_argument("AngularModel: the relative intensities are all 0.");

        AngularModel m;
        m.m_kind      = Kind::table;
        m.m_faces     = faces;
        m.m_theta_deg = std::move(theta_deg);
        m.m_f         = std::move(f);
        return m;
    }

    AngularModel AngularModel::load(const PixelConfig& entry)
    {
        Faces faces = Faces::front;
        if (entry.has("faces"))
        {
            std::string v;
            try                                    { v = entry.get<std::string>("faces"); }
            catch (const std::bad_variant_access&) { }
            if      (v == "front") faces = Faces::front;
            else if (v == "back")  faces = Faces::back;
            else if (v == "both")  faces = Faces::both;
            else throw std::invalid_argument("\"faces\" is \"front\", \"back\" or \"both\".");
        }

        if (!entry.has("angular"))
            return isotropic(faces);

        const json& a = entry.getJson("angular");   // std::invalid_argument if neither an object nor an array
        if (!a.is_object())
            throw std::invalid_argument("\"angular\" must be an object {\"model\": ...}.");
        const PixelConfig cfg(a);

        const std::string model = cfg.get<std::string>("model");
        if (model == "isotropic") return isotropic(faces);
        if (model == "surface")   return surface(faces);
        if (model == "cosn")      return cosn(cfg.getAsDouble("n"), faces);
        if (model == "table")
        {
            std::vector<double> theta, f;
            try
            {
                readColumns(cfg, "colTheta", "colIntensity", theta, f);
            }
            catch (const std::exception& ex)
            {
                throw std::invalid_argument(std::string("\"angular\" table: ") + ex.what());
            }
            return table(std::move(theta), std::move(f), faces);
        }
        throw std::invalid_argument("\"angular\": unknown model '" + model + "' (isotropic, surface, cosn or table).");
    }

    double AngularModel::relativeIntensity(double cosZenith) const
    {
        if (!zenith())
            return 1.0;
        if (!(cosZenith >= 0.0))
            return 0.0;   // upward-going
        const double cz = std::min(cosZenith, 1.0);
        if (m_kind == Kind::cosn)
            return std::pow(cz, m_n);

        const double theta = std::acos(cz) * 180.0 / kPi;
        if (theta < m_theta_deg.front() || theta > m_theta_deg.back())
            return 0.0;
        const size_t i = std::clamp<size_t>(static_cast<size_t>(std::upper_bound(m_theta_deg.begin(), m_theta_deg.end(), theta) - m_theta_deg.begin()),
                                            1, m_theta_deg.size() - 1);
        return m_f[i - 1] + (theta - m_theta_deg[i - 1]) / (m_theta_deg[i] - m_theta_deg[i - 1]) * (m_f[i] - m_f[i - 1]);
    }

    double AngularModel::maxRelativeIntensity() const
    {
        return (m_kind == Kind::table) ? *std::max_element(m_f.begin(), m_f.end()) : 1.0;   // cos^n <= 1
    }

#pragma endregion

#pragma region - geometric factor (no shielding)

    namespace
    {
        constexpr double kDeg = kPi / 180.0;   // rad per degree

        /// Incidence of a particle coming from u (unit vector towards the source, sky frame) on a sensor whose front
        /// normal is n = (sin a, 0, cos a): mu = |u.n|, the face (+1 when u.n >= 0), and the azimuth of its track
        /// -u on the matrix axes (cos a, 0, -sin a) and (0, 1, 0).
        Incidence incidenceOf(double ux, double uy, double uz, double sa, double ca)
        {
            const double dot = ux * sa + uz * ca;
            const double vx  = -(ux * ca - uz * sa), vy = -uy;
            double phi = std::atan2(vy, vx);
            if (phi < 0.0)
                phi += 2.0 * kPi;
            return {std::abs(dot), phi, (dot >= 0.0) ? 1 : -1};
        }

        /// G (sr) of a zenith model over mu in [muLo, muHi], the sensor tilted by alpha (rad): ∫ f(theta_z) mu dOmega
        /// over the downward-going directions that reach the faces of the model. Midpoint rule on cos(theta_z) in
        /// [0, 1] (400 steps) and on the azimuth (720 steps).
        double zenithFactor(const AngularModel& a, double alpha, double muLo, double muHi)
        {
            constexpr int nCos = 400, nPhi = 720;
            const double sa = std::sin(alpha), ca = std::cos(alpha);
            std::vector<double> cosPhi(nPhi);
            for (int j = 0; j < nPhi; ++j)
                cosPhi[j] = std::cos((j + 0.5) * 2.0 * kPi / nPhi);

            double G = 0.0;
            for (int i = 0; i < nCos; ++i)
            {
                const double cz = (i + 0.5) / nCos;
                const double f  = a.relativeIntensity(cz);
                if (!(f > 0.0))
                    continue;
                const double sz = std::sqrt(1.0 - cz * cz);
                for (int j = 0; j < nPhi; ++j)
                {
                    const double dot = sz * cosPhi[j] * sa + cz * ca;   // cos to the front normal
                    const double mu  = std::abs(dot);
                    if (((dot >= 0.0) ? a.front() : a.back()) && mu >= muLo && mu <= muHi)
                        G += f * mu;
                }
            }
            return G * (1.0 / nCos) * (2.0 * kPi / nPhi);
        }

        /// G (sr) of any model over mu in [muLo, muHi], the sensor tilted by alpha (rad); closed form where it exists.
        double factorInMu(const AngularModel& a, double alpha, double muLo, double muHi)
        {
            const double faces = (a.front() ? 1.0 : 0.0) + (a.back() ? 1.0 : 0.0);
            if (a.kind() == AngularModel::Kind::isotropic)                    // flux ∝ mu: d(mu^2)
                return kPi * (muHi * muHi - muLo * muLo) * faces;
            if (a.kind() == AngularModel::Kind::surface)                      // flux uniform in dOmega: d(mu)
                return 2.0 * kPi * (muHi - muLo) * faces;
            if (a.kind() == AngularModel::Kind::cosn && alpha == 0.0)         // horizontal sensor, from above
            {
                const double p = a.n() + 2.0;
                return a.front() ? 2.0 * kPi / p * (std::pow(muHi, p) - std::pow(muLo, p)) : 0.0;
            }
            return zenithFactor(a, alpha, muLo, muHi);
        }

        void checkTilt(double tilt_deg, const char* where)
        {
            if (!(tilt_deg >= 0.0 && tilt_deg <= 180.0))
                throw std::invalid_argument(std::string(where) + ": the tilt must be in [0, 180] deg.");
        }
    }

    Incidence incidence(double thetaZ_deg, double phiZ_deg, double tilt_deg)
    {
        checkTilt(tilt_deg, "incidence");
        if (!(thetaZ_deg >= 0.0 && thetaZ_deg <= 180.0) || !std::isfinite(phiZ_deg))
            throw std::invalid_argument("incidence: thetaZ_deg must be in [0, 180] and phiZ_deg finite.");
        const double t = thetaZ_deg * kDeg, p = phiZ_deg * kDeg, a = tilt_deg * kDeg;
        return incidenceOf(std::sin(t) * std::cos(p), std::sin(t) * std::sin(p), std::cos(t), std::sin(a), std::cos(a));
    }

    double effectiveArea(double area_m2, double thetaZ_deg, double phiZ_deg, double tilt_deg)
    {
        if (!(area_m2 > 0.0))
            throw std::invalid_argument("effectiveArea: the area must be > 0.");
        return area_m2 * incidence(thetaZ_deg, phiZ_deg, tilt_deg).mu;
    }

    double geometricFactor(const AngularModel& model, double tilt_deg, double thetaMin_deg, double thetaMax_deg)
    {
        checkTilt(tilt_deg, "geometricFactor");
        if (!(thetaMin_deg >= 0.0 && thetaMin_deg < thetaMax_deg && thetaMax_deg <= 90.0))
            throw std::invalid_argument("geometricFactor: the incidence angles must satisfy 0 <= thetaMin < thetaMax <= 90 deg.");
        const double muHi = (thetaMin_deg == 0.0)  ? 1.0 : std::cos(thetaMin_deg * kDeg);
        const double muLo = (thetaMax_deg == 90.0) ? 0.0 : std::cos(thetaMax_deg * kDeg);
        return factorInMu(model, tilt_deg * kDeg, muLo, muHi);
    }

    double etendue(const AngularModel& model, double area_m2, double tilt_deg, double thetaMin_deg, double thetaMax_deg)
    {
        if (!(area_m2 > 0.0))
            throw std::invalid_argument("etendue: the area must be > 0.");
        return area_m2 * geometricFactor(model, tilt_deg, thetaMin_deg, thetaMax_deg);
    }

#pragma endregion

#pragma region - CosmicSpectrum class implementation
#pragma region -- add species and EnergyBins
    
    void CosmicSpectrum::add(const specie& species, const std::vector<double>& E, const std::vector<double>& J,
                             const AngularModel& angular)
    {
        add(species, E, J, angular, Footprint());
    }

    void CosmicSpectrum::add(const specie& species, const std::vector<double>& E, const std::vector<double>& J,
                             const AngularModel& angular, const Footprint& footprint)
    {
        if (!(species.A() > 0.0))
            throw std::invalid_argument("CosmicSpectrum: species '" + species.name() + "' needs A > 0.");

        const double s = footprint.sigma_px;
        if (!(std::isnan(s) || (std::isfinite(s) && s >= 0.0)))
            throw std::invalid_argument("CosmicSpectrum: species '" + species.name() + "': sigmaDiffusion_px must be finite and >= 0.");
        if (!(footprint.haloFraction >= 0.0 && footprint.haloFraction <= 1.0))
            throw std::invalid_argument("CosmicSpectrum: species '" + species.name() + "': the halo fraction must be in [0, 1].");
        if (!(std::isfinite(footprint.haloSigma_px) && footprint.haloSigma_px >= 0.0) ||
            (footprint.haloFraction > 0.0 && !(footprint.haloSigma_px > 0.0)))
            throw std::invalid_argument("CosmicSpectrum: species '" + species.name() + "': the halo sigma_px must be finite and > 0.");
        if (footprint.deltaRays && species.mc2() < 1.0)   // e- and e+: rest energy below 1 MeV
            throw std::invalid_argument("CosmicSpectrum: species '" + species.name() + "': \"deltaRays\" is for heavy particles "
                                        "(nuclei, muons...), not for e- or e+.");

        if (E.size() < 2 || E.size() != J.size())
            throw std::invalid_argument("CosmicSpectrum: species '" + species.name() + "' needs at least 2 points, and as many E as J.");

        for (size_t i = 0; i < E.size(); ++i)
        {
            if (!(E[i] > 0.0) || !(J[i] > 0.0) || !std::isfinite(E[i]) || !std::isfinite(J[i]))
                throw std::invalid_argument("CosmicSpectrum: species '" + species.name() + "' needs finite E > 0 and J > 0 (power-law interpolation).");

            if (i > 0 && !(E[i] > E[i - 1]))
                throw std::invalid_argument("CosmicSpectrum: species '" + species.name() + "': E must be strictly ascending.");
        }

        std::vector<EnergyBin> bins = m_bins;
        for (size_t i = 0; i + 1 < E.size(); ++i)
            bins.push_back({m_species.size(), E[i], E[i + 1], J[i], J[i + 1], powerLawIntegral(E[i], E[i + 1], J[i], J[i + 1])});

        m_species.push_back(species);
        m_angular.push_back(angular);
        m_footprint.push_back(footprint);
        m_bins = std::move(bins);

        m_cdf.resize(m_bins.size());
        double sum = 0.0;
        for (size_t i = 0; i < m_bins.size(); ++i)
            m_cdf[i] = (sum += m_bins[i].flux);
        m_total = sum;
    }

    CosmicSpectrum CosmicSpectrum::load(const PixelConfig& cfg)
    {
        const json& spectra = cfg.getJson("spectrum");   // std::out_of_range if absent
        if (!spectra.is_array() || spectra.empty())
            throw std::invalid_argument("CosmicSpectrum: \"spectrum\" must be a non-empty array of species.");

        CosmicSpectrum s;
        for (const json& e : spectra)
        {
            if (!e.is_object())
                throw std::invalid_argument("CosmicSpectrum: a species is an object {\"species\", \"Z\", \"A\" or \"mc2\", \"file\", ...}.");
            const PixelConfig entry(e);

            std::string name;
            specie species;
            AngularModel angular;
            CosmicSpectrum::Footprint footprint;
            std::vector<double> E, J;
            try
            {
                name = entry.get<std::string>("species");
                const double Z = entry.getAsDouble("Z");

                // A nucleus (energies per nucleon) or any other particle from its rest energy (energies per particle)
                if (entry.has("A") == entry.has("mc2"))
                    throw std::invalid_argument("give \"A\" (a nucleus, energies per nucleon) or \"mc2\" (rest energy in MeV of "
                                                "any other particle, energies per particle), not both.");
                species = entry.has("mc2") ? specie::fromRestEnergy(name, Z, entry.getAsDouble("mc2"))
                                           : specie(name, Z, entry.getAsDouble("A"));

                readSpectrum(entry, E, J);
                angular = AngularModel::load(entry);
                footprint = readFootprint(entry);
            }
            catch (const std::exception& ex)
            {
                throw std::invalid_argument("CosmicSpectrum: species '" + (name.empty() ? std::string("?") : name) + "': " + ex.what());
            }

            s.add(species, E, J, angular, footprint);
        }
        return s;
    }

    double CosmicSpectrum::integralIntensity(size_t species) const
    {
        double sum = 0.0;
        for (const EnergyBin& c : m_bins)
            if (c.species == species) sum += c.flux;
        return sum;
    }

    std::vector<double> CosmicSpectrum::binProbabilities() const
    {
        std::vector<double> p(m_bins.size());
        for (size_t i = 0; i < p.size(); ++i) p[i] = m_bins[i].flux / m_total;
        return p;
    }

    std::vector<double> CosmicSpectrum::speciesProbabilities() const
    {
        std::vector<double> p(m_species.size(), 0.0);
        for (const EnergyBin& c : m_bins) p[c.species] += c.flux / m_total;
        return p;
    }

    double CosmicSpectrum::energyInBin(size_t bin, double u) const
    {
        const EnergyBin& c = m_bins.at(bin);
        return powerLawQuantile(c.E0, c.E1, c.J0, c.J1, u);
    }

    CosmicSpectrum::Sample CosmicSpectrum::sample(double uCell, double uEnergy) const
    {
        const auto it = std::upper_bound(m_cdf.begin(), m_cdf.end(), uCell * m_total);

        const size_t i = std::min<size_t>(static_cast<size_t>(it - m_cdf.begin()), m_bins.size() - 1);
        const EnergyBin& c = m_bins[i];

        return {c.species, powerLawQuantile(c.E0, c.E1, c.J0, c.J1, uEnergy)};
    }

#pragma endregion
#pragma region - StoppingTable

    namespace
    {
        /// Linear interpolation of ys(xs) at x, xs ascending with 2 values at least; the end intervals extrapolate.
        double interpolate(const std::vector<double>& xs, const std::vector<double>& ys, double x)
        {
            const size_t i = std::clamp<size_t>(static_cast<size_t>(std::upper_bound(xs.begin(), xs.end(), x) - xs.begin()),
                                                1, xs.size() - 1);
            return ys[i - 1] + (x - xs[i - 1]) / (xs[i] - xs[i - 1]) * (ys[i] - ys[i - 1]);
        }
    }

    StoppingTable::StoppingTable(const std::function<double(double)>& S, double Tc, double Tmax, size_t perDecade)
    {
        if (!(Tc > 0.0) || perDecade == 0)
            throw std::invalid_argument("StoppingTable: Tc must be > 0 and perDecade >= 1.");

        const auto stopping = [&S](double T)
        {
            const double s = S(T);
            if (!(std::isfinite(s) && s > 0.0))
                throw std::invalid_argument("StoppingTable: the stopping power is not > 0 at T = " + std::to_string(T) + " MeV.");
            return s;
        };

        const double Sc = stopping(Tc);
        const double Rc = 2.0 * Tc / Sc;   // ∫0^Tc dT / (Sc sqrt(T / Tc))

        std::vector<double> lnT{std::log(Tc)}, lnR{std::log(Rc)};
        const double step   = std::numbers::ln10 / static_cast<double>(perDecade);
        const double lnTmax = std::log(std::max(Tmax, Tc));
        double T0 = Tc, invS0 = 1.0 / Sc, R = Rc;
        for (size_t k = 1; lnT.size() < 2 || lnT.back() < lnTmax; ++k)
        {
            const double lnT1  = lnT.front() + step * static_cast<double>(k);
            const double T1    = std::exp(lnT1);
            const double invS1 = 1.0 / stopping(T1);
            R += (T1 - T0) / 6.0 * (invS0 + 4.0 / stopping(0.5 * (T0 + T1)) + invS1);   // Simpson
            lnT.push_back(lnT1);
            lnR.push_back(std::log(R));
            T0    = T1;
            invS0 = invS1;
        }

        m_Tc  = Tc;
        m_Sc  = Sc;
        m_Rc  = Rc;
        m_lnT = std::move(lnT);
        m_lnR = std::move(lnR);
    }

    double StoppingTable::range(double T) const
    {
        if (empty() || !(T > 0.0))
            return 0.0;
        if (T <= m_Tc)
            return 2.0 * std::sqrt(T * m_Tc) / m_Sc;
        return std::exp(interpolate(m_lnT, m_lnR, std::log(T)));
    }

    double StoppingTable::energy(double R) const
    {
        if (empty() || !(R > 0.0))
            return 0.0;
        if (R <= m_Rc)
        {
            const double x = 0.5 * R * m_Sc;   // sqrt(T Tc)
            return x * x / m_Tc;
        }
        return std::exp(interpolate(m_lnR, m_lnT, std::log(R)));
    }

#pragma endregion
#pragma region - CosmicDeposit: configuration

    CosmicDeposit::CosmicDeposit(CosmicSpectrum spectrum, const Params& params) : m_params(params), m_spectrum(std::move(spectrum))
    {
        if (!isConfigured())
            throw std::invalid_argument("CosmicDeposit: a spectrum, depth_um, pitch_um and eVPerElectron are required.");
        buildTables();
    }

    bool CosmicDeposit::isConfigured() const
    {
        const auto ok = [](double v) { return std::isfinite(v) && v > 0.0; };
        return !m_spectrum.empty() && 
                ok(m_params.depth_um) &&
                ok(m_params.pitch_um) &&
                ok(m_params.eVPerElectron) &&
                m_params.tilt_deg >= 0.0 && m_params.tilt_deg <= 180.0 &&
                m_params.sigmaDiffusion_px >= 0.0 &&
                ok(m_params.maxTrackPixels) &&
                ok(m_params.deltaRayCut_e);
    }

    void CosmicDeposit::configure(const PixelConfig& cfg)
    {
        if (!cfg.has("cosmic"))
            return;

        const json& j = cfg.getJson("cosmic");   // std::invalid_argument if neither an array nor an object
        
        if (!j.is_object())
            throw std::invalid_argument("CosmicDeposit: 'cosmic' must be an object.");
        
        load(PixelConfig(j));
    }

    void CosmicDeposit::load(const PixelConfig& cfg)
    {
        Params p = m_params;
        const auto number = [&cfg](const char* key, double& dst)
        {
            if (!cfg.has(key)) return;
            try                           { dst = cfg.getAsDouble(key); }
            catch (const std::exception&) { throw std::invalid_argument(std::string("CosmicDeposit: '") + key + "' must be a number."); }
        };

        if (cfg.has("acceptance_sr"))
            throw std::invalid_argument("CosmicDeposit: 'acceptance_sr' is no longer read: the geometric factor comes from the "
                                        "angular model of each spectrum entry ('angular', 'faces').");

        number("depth_um",          p.depth_um);
        number("pitch_um",          p.pitch_um);
        number("eVPerElectron",     p.eVPerElectron);
        number("tilt_deg",          p.tilt_deg);
        number("sigmaDiffusion_px", p.sigmaDiffusion_px);
        number("maxTrackPixels",    p.maxTrackPixels);
        number("deltaRayCut_e",     p.deltaRayCut_e);

        std::optional<uint64_t> seed;

        if (cfg.has("seed"))
        {
            // Signed or unsigned integer, as JSON or the caller stored it
            const json v = cfg.toJson().at("seed");

            if (!v.is_number_unsigned() && !(v.is_number_integer() && v.get<int64_t>() >= 0))
                throw std::invalid_argument("CosmicDeposit: 'seed' must be an integer >= 0.");
            
            seed = v.get<uint64_t>();
        }

        CosmicSpectrum spectrum = m_spectrum;
        if (cfg.has("spectrum"))
            spectrum = CosmicSpectrum::load(cfg);

        // Check on a copy, commit afterwards.
        CosmicDeposit check(*this);
        check.m_params   = p;
        check.m_spectrum = spectrum;
        
        if (!check.isConfigured())
            throw std::invalid_argument("CosmicDeposit: a spectrum, depth_um, pitch_um and eVPerElectron are required, "
                                        "and maxTrackPixels > 0, deltaRayCut_e > 0, sigmaDiffusion_px >= 0, tilt_deg in [0, 180].");
        check.buildTables();   // may throw: still nothing committed

        m_params     = p;
        m_spectrum   = std::move(spectrum);
        m_tables     = std::move(check.m_tables);
        m_acceptance = std::move(check.m_acceptance);
        m_drawCdf    = std::move(check.m_drawCdf);
        m_kernel     = std::move(check.m_kernel);
        m_kernels    = std::move(check.m_kernels);
        m_electronTable = std::move(check.m_electronTable);
        m_deltaT0       = check.m_deltaT0;
        if (seed) setSeed(*seed);
    }

    namespace
    {
        /// Gaussian of sigma s (pixel) integrated over the pixels k in [-R, R], R = ceil(3 s), normalised to 1
        /// (the charge is conserved within +-R); {1} for s = 0.
        std::vector<double> gaussianShares(double s)
        {
            if (!(s > 0.0))
                return {1.0};

            const long R = static_cast<long>(std::ceil(3.0 * s));
            std::vector<double> g(2 * R + 1, 0.0);
            double sum = 0.0;
            for (long k = -R; k <= R; ++k)
                sum += (g[k + R] = DST::Math::function::integrated_gauss<double>(static_cast<double>(k), 1.0, 0.0, s, 1.0));
            for (double& v : g) v /= sum;
            return g;
        }
    }

    CosmicDeposit::Kernel::Kernel(double sigma_px, double haloSigma_px, double haloFraction)
    {
        const bool halo = haloFraction > 0.0;
        const std::vector<double> core = gaussianShares(sigma_px);
        const std::vector<double> wide = halo ? gaussianShares(haloSigma_px) : std::vector<double>{};

        R = static_cast<long>(std::max(core.size(), wide.size()) / 2);
        const long n = 2 * R + 1;
        w.assign(static_cast<size_t>(n * n), 0.0);

        // share * g(a) g(b), centred in the (2R+1)^2 square
        const auto add = [this, n](const std::vector<double>& g, double share)
        {
            const long r = static_cast<long>(g.size() / 2);
            for (long b = -r; b <= r; ++b)
                for (long a = -r; a <= r; ++a)
                    w[static_cast<size_t>((b + R) * n + a + R)] += share * g[a + r] * g[b + r];
        };
        add(core, 1.0 - haloFraction);
        if (halo)
            add(wide, haloFraction);
    }

#pragma endregion
#pragma region - CosmicDeposit: probabilities

    double CosmicDeposit::eventRate(double area_m2) const
    {
        double sum = 0.0;   // sum over the species of G ∫ J dE (m-2 s-1)
        for (size_t i = 0; i < m_acceptance.size(); ++i)
            sum += m_acceptance[i].G * m_spectrum.integralIntensity(i);
        return sum * area_m2;
    }

    double CosmicDeposit::expectedEvents(double area_m2, double dt) const
    {
        return eventRate(area_m2) * dt;
    }

    double CosmicDeposit::probabilityAtLeastOne(double area_m2, double dt) const
    {
        return -std::expm1(-expectedEvents(area_m2, dt));
    }

    double CosmicDeposit::probability(size_t n, double area_m2, double dt) const
    {
        const double lambda = expectedEvents(area_m2, dt);
        if (lambda <= 0.0) return n == 0 ? 1.0 : 0.0;
        return std::exp(-lambda + static_cast<double>(n) * std::log(lambda) - std::lgamma(static_cast<double>(n) + 1.0));
    }

    double CosmicDeposit::truncatedFraction() const
    {
        double all = 0.0, cut = 0.0;   // hits of every track length, hits not drawn
        for (size_t i = 0; i < m_acceptance.size(); ++i)
        {
            const double w = m_acceptance[i].G * m_spectrum.integralIntensity(i);
            all += w;
            cut += w * m_acceptance[i].truncated;
        }
        return (all > 0.0) ? cut / all : 0.0;
    }

#pragma endregion
#pragma region - CosmicDeposit: energy loss

    /// Density effect, Sternheimer parametrisation, x = log10(beta gamma):
    ///   x >= x1      : delta = 2 ln(10) x - Cbar
    ///   x0 <= x < x1 : delta = 2 ln(10) x - Cbar + a (x1 - x)^m
    ///   x < x0       : delta = delta0 10^(2 (x - x0))   (conductors; 0 for an insulator)
    /// References:
    ///  - R. M. Sternheimer, M. J. Berger, S. M. Seltzer, "Density effect for the ionization loss of charged
    ///    particles in various substances", At. Data Nucl. Data Tables 30 (1984) 261.
    ///  - R. M. Sternheimer, "The density effect for the ionization loss in various materials",
    ///    Phys. Rev. 88 (1952) 851.
    ///  - R. M. Sternheimer, R. F. Peierls, "General expression for the density effect for the ionization loss
    ///    of charged particles", Phys. Rev. B 3 (1971) 3681 (parameters of target::material::compound()).
    ///  - PDG, "Passage of particles through matter" (D. E. Groom, S. R. Klein), section "Density effect",
    ///    in S. Navas et al. (Particle Data Group), Phys. Rev. D 110 (2024) 030001.
    double CosmicDeposit::densityEffect(double betaGamma) const
    {
        const target::material& m = m_params.material;
        const double x = std::log10(betaGamma);
        const double twoLn10x = 2.0 * std::numbers::ln10 * x;

        if (x >= m.sternheimerX1) return twoLn10x - m.sternheimerC;
        if (x >= m.sternheimerX0) return twoLn10x - m.sternheimerC + m.sternheimerA * std::pow(m.sternheimerX1 - x, m.sternheimerM);
        
        return m.sternheimerDelta0 * std::pow(10.0, 2.0 * (x - m.sternheimerX0));
    }

    /// Landau width and most probable energy loss in a thin layer of thickness x (g/cm2):
    ///   xi = (K/2) <Z/A> z^2 x / beta^2
    ///   Dp = xi [ln(2 m c^2 beta^2 gamma^2 / I) + ln(xi / I) + j - beta^2 - delta(beta gamma)],  j = 0.200
    ///   K  = 4 pi N_A r_e^2 m_e c^2 = 0.307075 MeV cm2 mol-1
    /// References:
    ///  - PDG, "Passage of particles through matter" (D. E. Groom, S. R. Klein), section "Fluctuations in
    ///    energy loss", in S. Navas et al. (Particle Data Group), Phys. Rev. D 110 (2024) 030001.
    ///  - H. Bichsel, "Straggling in thin silicon detectors", Rev. Mod. Phys. 60 (1988) 663 (Dp and j).
    ///  - L. Landau, "On the energy loss of fast particles by ionization", J. Phys. USSR 8 (1944) 201.
    CosmicDeposit::EnergyLoss CosmicDeposit::energyLoss(const specie& species, double T, double path_um) const
    {
        const target::material& m = m_params.material;
        const double gamma = species.gamma(T);   // T per nucleon for a nucleus, per particle for e±
        const double bg2   = gamma * gamma - 1.0;          // (beta gamma)^2
        const double beta2 = bg2 / (gamma * gamma);
        const double x     = path_um * 1.0e-4 * m.density_gcm3;   // g/cm2
        const double I     = m.meanExcitation_eV * 1.0e-6;        // MeV
        const double xi    = 0.5 * kK_MeV_cm2_mol * m.ZoverA * species.Z() * species.Z() / beta2 * x;   // MeV
        const double mpv   = xi * (std::log(2.0 * kElectronMass_MeV * bg2 / I) + std::log(xi / I) + kLandauJ
                                   - beta2 - densityEffect(std::sqrt(bg2)));
        return {std::max(mpv, 0.0) * 1.0e6, xi * 1.0e6};
    }

    /// Energy deposited: Landau of most probable value Dp and width xi, dE = Dp + xi (lambda - lambda_mode),
    /// lambda drawn from the standard Landau (DST::Math::standardLandau), capped at the kinetic energy.
    /// Known limitation: Landau is used even where kappa = xi / Wmax > 0.01 (Vavilov regime, e.g. protons below
    /// about 100 MeV in 10 um of HgCdTe). There the width and the upper tail of the deposit are too large
    /// (about 1.7 times at kappa = 1), not its typical value. No Vavilov distribution is implemented.
    /// Reference for the standard Landau density: K. S. Kolbig, B. Schorr, "A program package for the Landau
    /// distribution", Comput. Phys. Commun. 31 (1984) 97 (CERNLIB G110, DENLAN).
    double CosmicDeposit::energyDeposit_eV(const specie& species, double T, double path_um, double u) const
    {
        const EnergyLoss el = energyLoss(species, T, path_um);
        const double lambda = DST::Math::standardLandau::instance().quantile(u);
        const double dE     = el.mpv_eV + el.xi_eV * (lambda - DST::Math::standardLandau::lambda_mode);
        const double KE     = species.kineticEnergy(T) * 1.0e6;   // eV: a particle cannot deposit more than it has
        return std::clamp(dE, 0.0, KE);
    }

#pragma endregion
#pragma region - CosmicDeposit: stopping power and range

    namespace
    {
        /// e- and e+: rest energy below 1 MeV (their own collision formula, energies per particle).
        bool isElectronLike(const specie& s) { return s.mc2() < 1.0; }

        /// Largest energy (MeV) given to an atomic electron by a particle of kinetic energy T (MeV).
        double maxEnergyTransfer(const specie& s, double T)
        {
            if (isElectronLike(s))
                return (s.Z() < 0.0) ? 0.5 * T : T;   // e-: identical particles; e+: all of it
            const double gamma = 1.0 + T / s.mc2();
            const double r     = kElectronMass_MeV / s.mc2();
            return 2.0 * kElectronMass_MeV * (gamma * gamma - 1.0) / (1.0 + 2.0 * gamma * r + r * r);
        }

        /// Kinetic energy (MeV) where S is largest, searched from 1e-6 mc2 to 1e-2 mc2 (200 steps per decade):
        /// below it, the stopping power is taken proportional to the velocity. 0 if S is never > 0.
        double peakEnergy(const std::function<double(double)>& S, double mc2)
        {
            double bestT = 0.0, bestS = 0.0;
            for (int k = 0; k <= 800; ++k)
            {
                const double T = mc2 * std::pow(10.0, -6.0 + static_cast<double>(k) / 200.0);
                const double s = S(T);
                if (std::isfinite(s) && s > bestS)
                {
                    bestS = s;
                    bestT = T;
                }
            }
            return bestT;
        }
    }

    /// Mean collision stopping power (MeV cm2/g), T the kinetic energy of the particle (MeV):
    ///  - heavy particle (nucleus, muon, ...), Bethe-Bloch:
    ///      S = K z^2 <Z/A> / beta^2 [1/2 ln(2 m c^2 beta^2 gamma^2 Wmax / I^2) - beta^2 - delta / 2],
    ///      Wmax = 2 m c^2 beta^2 gamma^2 / (1 + 2 gamma m / M + (m / M)^2);
    ///  - e- and e+ (mc2 < 1 MeV), tau = T / m c^2:
    ///      S = K/2 <Z/A> / beta^2 [ln(tau^2 (tau + 2) / (2 (I / m c^2)^2)) + F(tau) - delta],
    ///      F-(tau) = 1 - beta^2 + [tau^2 / 8 - (2 tau + 1) ln 2] / (tau + 1)^2,
    ///      F+(tau) = 2 ln 2 - beta^2 / 12 [23 + 14 / (tau + 2) + 10 / (tau + 2)^2 + 4 / (tau + 2)^3].
    /// Radiative losses of electrons are left out: the bremsstrahlung photons leave a thin layer.
    /// Known limitation: no shell, Barkas or Bloch corrections, no effective charge of slow ions. Below the
    /// energy where the formula is largest, stoppingPower() and StoppingTable take a stopping power
    /// proportional to the velocity instead. Low-energy ranges are therefore too long in heavy materials:
    /// a 5 MeV alpha in HgCdTe gets about 40 um, against about 20 um (Bragg-Kleeman estimate). A stopping
    /// table given by the user (SRIM, NIST ASTAR/PSTAR) is for later.
    /// References:
    ///  - PDG, "Passage of particles through matter" (D. E. Groom, S. R. Klein), the Bethe equation, in
    ///    S. Navas et al. (Particle Data Group), Phys. Rev. D 110 (2024) 030001.
    ///  - ICRU Report 37, "Stopping powers for electrons and positrons" (1984).
    ///  - F. Rohrlich, B. C. Carlson, "Positron-electron differences in energy loss and multiple scattering",
    ///    Phys. Rev. 93 (1954) 38.
    double CosmicDeposit::collisionStopping(const specie& species, double T) const
    {
        const target::material& m = m_params.material;
        const double gamma = 1.0 + T / species.mc2();
        const double bg2   = gamma * gamma - 1.0;          // (beta gamma)^2
        const double beta2 = bg2 / (gamma * gamma);
        const double I     = m.meanExcitation_eV * 1.0e-6;   // MeV
        const double delta = densityEffect(std::sqrt(bg2));

        if (isElectronLike(species))
        {
            const double tau = T / kElectronMass_MeV;
            const double i   = I / kElectronMass_MeV;
            const double t2  = tau + 2.0;
            const double F   = (species.Z() < 0.0)
                ? 1.0 - beta2 + (tau * tau / 8.0 - (2.0 * tau + 1.0) * std::numbers::ln2) / ((tau + 1.0) * (tau + 1.0))
                : 2.0 * std::numbers::ln2 - beta2 / 12.0 * (23.0 + 14.0 / t2 + 10.0 / (t2 * t2) + 4.0 / (t2 * t2 * t2));
            return 0.5 * kK_MeV_cm2_mol * m.ZoverA / beta2 * (std::log(tau * tau * t2 / (2.0 * i * i)) + F - delta);
        }

        const double Wmax = maxEnergyTransfer(species, T);
        return kK_MeV_cm2_mol * species.Z() * species.Z() * m.ZoverA / beta2
             * (0.5 * std::log(2.0 * kElectronMass_MeV * bg2 * Wmax / (I * I)) - beta2 - 0.5 * delta);
    }

    double CosmicDeposit::stoppingPower(const specie& species, double E) const
    {
        const auto S  = [this, &species](double T) { return collisionStopping(species, T); };
        const double T  = species.kineticEnergy(E);
        const double Tc = peakEnergy(S, species.mc2());
        return (T < Tc) ? S(Tc) * std::sqrt(T / Tc) : S(T);
    }

    void CosmicDeposit::buildTables()
    {
        const std::vector<specie>& species = m_spectrum.species();
        std::vector<StoppingTable> tables;
        tables.reserve(species.size());
        double WmaxTop = 0.0;   // largest Wmax of the species with delta rays

        for (size_t i = 0; i < species.size(); ++i)
        {
            double Emax = 0.0;   // top of the spectrum of the species
            for (const CosmicSpectrum::EnergyBin& b : m_spectrum.bins())
                if (b.species == i)
                    Emax = std::max(Emax, b.E1);

            const specie& s = species[i];
            if (m_spectrum.footprints()[i].deltaRays)
                WmaxTop = std::max(WmaxTop, maxEnergyTransfer(s, s.kineticEnergy(Emax)));
            const auto S = [this, &s](double T) { return collisionStopping(s, T); };
            try
            {
                tables.emplace_back(S, peakEnergy(S, s.mc2()), 1.01 * s.kineticEnergy(Emax));
            }
            catch (const std::invalid_argument& ex)
            {
                throw std::invalid_argument("CosmicDeposit: species '" + s.name() + "': " + ex.what());
            }
        }

        // Geometric factor and truncated share of each species. Tracks longer than maxTrackPixels, i.e.
        // mu < muMin = 1 / sqrt(1 + q^2), are not drawn.
        const double q      = m_params.maxTrackPixels * m_params.pitch_um / m_params.depth_um;
        const double mu2min = 1.0 / (1.0 + q * q);
        const double muMin  = std::sqrt(mu2min);
        const double alpha  = m_params.tilt_deg * kPi / 180.0;

        std::vector<Acceptance> acceptance;
        acceptance.reserve(species.size());
        for (const AngularModel& a : m_spectrum.angular())
        {
            Acceptance acc;   // the same geometric factor as geometricFactor(), over the full sky and below muMin
            acc.muMin     = muMin;
            acc.G         = factorInMu(a, alpha, 0.0, 1.0);
            acc.truncated = (acc.G > 0.0) ? factorInMu(a, alpha, 0.0, muMin) / acc.G : 0.0;
            acceptance.push_back(acc);
        }

        // Weights of the draws: flux * G * (1 - truncated) of each bin
        std::vector<double> cdf;
        cdf.reserve(m_spectrum.bins().size());
        double sum = 0.0;
        for (const CosmicSpectrum::EnergyBin& b : m_spectrum.bins())
        {
            const Acceptance& acc = acceptance[b.species];
            sum += b.flux * acc.G * (1.0 - acc.truncated);
            cdf.push_back(sum);
        }

        // Lateral spread of each species: its own core sigma, else sigmaDiffusion_px, and its halo
        std::vector<Kernel> kernels;
        kernels.reserve(species.size());
        for (const CosmicSpectrum::Footprint& f : m_spectrum.footprints())
            kernels.emplace_back(std::isnan(f.sigma_px) ? m_params.sigmaDiffusion_px : f.sigma_px, f.haloSigma_px, f.haloFraction);

        // Delta rays: CSDA range of electrons up to the largest Wmax (1 MeV at least), and T0, the energy of an
        // electron whose range is half a pixel
        StoppingTable electronTable;
        double deltaT0 = 0.0;
        if (WmaxTop > 0.0)
        {
            const auto Se = [this](double T) { return collisionStopping(electron, T); };
            try
            {
                electronTable = StoppingTable(Se, peakEnergy(Se, electron.mc2()), std::max(1.01 * WmaxTop, 1.0));
            }
            catch (const std::invalid_argument& ex)
            {
                throw std::invalid_argument(std::string("CosmicDeposit: delta rays: ") + ex.what());
            }
            deltaT0 = electronTable.energy(0.5 * m_params.pitch_um * 1.0e-4 * m_params.material.density_gcm3);
        }

        m_tables     = std::move(tables);
        m_acceptance = std::move(acceptance);
        m_drawCdf    = std::move(cdf);
        m_kernel     = Kernel(m_params.sigmaDiffusion_px, 0.0, 0.0);
        m_kernels    = std::move(kernels);
        m_electronTable = std::move(electronTable);
        m_deltaT0       = deltaT0;
    }

#pragma endregion
#pragma region - CosmicDeposit: generation

    /// Hits of a step. For each hit, path L = depth / mu in the layer and CSDA range R of the particle:
    ///  - stops (R <= L): path R, projected track shortened to R / L, deposit = all its kinetic energy;
    ///  - thick (CSDA mean loss over L above kThickShare of T): Gaussian around the CSDA mean loss, Bohr variance
    ///    sigma^2 = xi Wmax (1 - beta^2 / 2) at the mid-path energy, between 0 and T;
    ///  - thin: energyDeposit_eV() (Landau).
    /// Reference for the Bohr variance: PDG, "Passage of particles through matter", section "Fluctuations in
    /// energy loss".
    std::vector<CosmicDeposit::Event> CosmicDeposit::drawEvents(double dt, const Geometry& g, std::mt19937_64& rng) const
    {
        if (!isConfigured())
            throw std::logic_error("CosmicDeposit: not configured (spectrum, depth_um, pitch_um, eVPerElectron).");
        if (m_tables.size() != m_spectrum.species().size() || m_acceptance.size() != m_spectrum.species().size() ||
            m_kernels.size() != m_spectrum.species().size() || m_drawCdf.size() != m_spectrum.bins().size())
            throw std::logic_error("CosmicDeposit: the tables do not match the spectrum (call buildTables()).");

        std::vector<Event> events;

        if (!(dt > 0.0))
            return events;

        // Entry points uniform on the whole physical detector; hits drawn: flux * G * (1 - truncated) of each bin.
        const double nx = static_cast<double>(g.width), ny = static_cast<double>(g.height);
        const double area_m2 = nx * ny * m_params.pitch_um * m_params.pitch_um * 1.0e-12;
        const double drawn   = m_drawCdf.empty() ? 0.0 : m_drawCdf.back();   // m-2 s-1
        const double lambda  = drawn * area_m2 * dt;
        if (!(lambda > 0.0))
            return events;

        const double um_to_gcm2 = 1.0e-4 * m_params.material.density_gcm3;

        std::uniform_real_distribution<double> U(0.0, 1.0);
        std::normal_distribution<double>       N(0.0, 1.0);
        const size_t n = std::poisson_distribution<size_t>(lambda)(rng);
        events.reserve(n);

        for (size_t i = 0; i < n; ++i)
        {
            Event e{};
            const auto   it  = std::upper_bound(m_drawCdf.begin(), m_drawCdf.end(), U(rng) * drawn);
            const size_t bin = std::min<size_t>(static_cast<size_t>(it - m_drawCdf.begin()), m_drawCdf.size() - 1);
            e.species = m_spectrum.bins()[bin].species;
            e.energy  = m_spectrum.energyInBin(bin, U(rng));

            const specie&        sp    = m_spectrum.species()[e.species];
            const StoppingTable& table = m_tables[e.species];
            const double         T     = sp.kineticEnergy(e.energy);   // MeV

            drawDirection(e.species, rng, e.mu, e.phi, e.face);
            const double L = m_params.depth_um / e.mu;            // geometric path through the depth (um)
            const double R = table.range(T) / um_to_gcm2;         // CSDA range (um)
            e.path_um = std::min(L, R);

            const double ell = m_params.depth_um / m_params.pitch_um * std::sqrt(1.0 - e.mu * e.mu) / e.mu * (e.path_um / L);   // pixels
            e.dx = ell * std::cos(e.phi);
            e.dy = ell * std::sin(e.phi);

            e.x = U(rng) * nx;
            e.y = U(rng) * ny;

            double dE_eV = 0.0;
            if (R <= L)                                                    // stops in the layer
                dE_eV = T * 1.0e6;
            else
            {
                const double Tout = table.energyAfter(T, L * um_to_gcm2);
                const double mean = T - Tout;                              // CSDA mean loss (MeV)
                if (mean > kThickShare * T)                                // thick: Gaussian, Bohr straggling
                {
                    const double Tm    = 0.5 * (T + Tout);
                    const double gm    = 1.0 + Tm / sp.mc2();
                    const double beta2 = 1.0 - 1.0 / (gm * gm);
                    const double xi    = 0.5 * kK_MeV_cm2_mol * m_params.material.ZoverA * sp.Z() * sp.Z() / beta2 * L * um_to_gcm2;
                    const double sigma = std::sqrt(xi * maxEnergyTransfer(sp, Tm) * (1.0 - 0.5 * beta2));
                    dE_eV = std::clamp(mean + sigma * N(rng), 0.0, T) * 1.0e6;
                }
                else                                                       // thin: Landau
                    dE_eV = energyDeposit_eV(sp, e.energy, L, U(rng));
            }

            e.electrons = dE_eV / m_params.eVPerElectron;
            events.push_back(e);
        }
        return events;
    }

    /// Direction of a hit in the frame of the sensor (mu = cos to the normal, phi on the matrix, face), mu >= muMin:
    ///  - isotropic: mu^2 uniform; surface: mu uniform; cos^n on a horizontal sensor: mu^(n+2) uniform; phi uniform;
    ///  - other zenith models: directions u (towards the source) uniform over the upper hemisphere, kept with
    ///    probability f(theta_z) |u.n| / max f, n = (sin a, 0, cos a) the front normal tilted by a = tilt_deg towards
    ///    +x of the matrix; the track runs along -u, projected on the matrix axes (cos a, 0, -sin a) and (0, 1, 0).
    void CosmicDeposit::drawDirection(size_t i, std::mt19937_64& rng, double& mu, double& phi, int& face) const
    {
        const AngularModel& a   = m_spectrum.angular()[i];
        const Acceptance&   acc = m_acceptance[i];
        std::uniform_real_distribution<double> U(0.0, 1.0);

        // Face of a sensor-frame model: both faces alike
        const auto sensorFace = [&a, &U, &rng]() { return (a.front() && a.back()) ? ((U(rng) < 0.5) ? 1 : -1) : (a.front() ? 1 : -1); };

        if (a.kind() == AngularModel::Kind::isotropic)
        {
            const double m2 = acc.muMin * acc.muMin;
            mu   = std::sqrt(m2 + U(rng) * (1.0 - m2));
            phi  = 2.0 * kPi * U(rng);
            face = sensorFace();
            return;
        }
        if (a.kind() == AngularModel::Kind::surface)
        {
            mu   = acc.muMin + U(rng) * (1.0 - acc.muMin);
            phi  = 2.0 * kPi * U(rng);
            face = sensorFace();
            return;
        }

        const double alpha = m_params.tilt_deg * kPi / 180.0;
        if (a.kind() == AngularModel::Kind::cosn && alpha == 0.0)   // horizontal sensor, from above
        {
            const double p   = a.n() + 2.0;
            const double low = std::pow(acc.muMin, p);
            mu   = std::pow(low + U(rng) * (1.0 - low), 1.0 / p);
            phi  = 2.0 * kPi * U(rng);
            face = 1;
            return;
        }

        const double sa = std::sin(alpha), ca = std::cos(alpha), fmax = a.maxRelativeIntensity();
        for (size_t k = 0; k < 100000000; ++k)
        {
            const double cz = U(rng), ph = 2.0 * kPi * U(rng);
            const double sz = std::sqrt(1.0 - cz * cz);
            const Incidence in = incidenceOf(sz * std::cos(ph), sz * std::sin(ph), cz, sa, ca);   // u towards the source
            if ((in.face == 1) ? !a.front() : !a.back())
                continue;
            if (in.mu < acc.muMin || U(rng) * fmax > a.relativeIntensity(cz) * in.mu)
                continue;

            mu   = in.mu;
            phi  = in.phi;
            face = in.face;
            return;
        }
        throw std::logic_error("CosmicDeposit: no direction drawn for species '" + m_spectrum.species()[i].name() +
                               "' (geometric factor about 0).");
    }

    /// Delta-ray halo (see CosmicDeposit.h, 5.), at the velocity of entry. In units of (K/2) z^2 <Z/A> / beta^2 per
    /// g/cm2, the delta rays from T0 to Wmax carry ln(Wmax / T0) - beta^2 (1 - T0 / Wmax) (dN/dT of a spin-0
    /// particle) and the Bethe mean loss is ln(2 m c^2 beta^2 gamma^2 Wmax / I^2) - 2 beta^2 - delta.
    /// References:
    ///  - PDG, "Passage of particles through matter" (D. E. Groom, S. R. Klein), sections "Maximum energy transfer
    ///    in a single collision" and "Energetic knock-on electrons (delta rays)", in S. Navas et al. (Particle
    ///    Data Group), Phys. Rev. D 110 (2024) 030001.
    ///  - J. J. Butts, R. Katz, "Theory of RBE for heavy ion bombardment of dry enzymes and viruses",
    ///    Radiat. Res. 30 (1967) 855 (radial dose of the delta rays, 1/r^2).
    CosmicDeposit::DeltaRays CosmicDeposit::deltaRays(size_t i, double E) const
    {
        DeltaRays d;
        if (i >= m_spectrum.species().size() || !m_spectrum.footprints()[i].deltaRays || m_electronTable.empty())
            return d;

        const specie& s = m_spectrum.species()[i];
        const target::material& m = m_params.material;
        const double gamma = s.gamma(E);
        const double bg2   = gamma * gamma - 1.0;          // (beta gamma)^2
        const double beta2 = bg2 / (gamma * gamma);
        const double Wmax  = maxEnergyTransfer(s, s.kineticEnergy(E));
        const double T0    = m_deltaT0;
        if (!(Wmax > T0) || !(T0 > 0.0))
            return d;   // every delta ray stays in the pixel of the track

        const double I     = m.meanExcitation_eV * 1.0e-6;   // MeV
        const double bethe = std::log(2.0 * kElectronMass_MeV * bg2 * Wmax / (I * I)) - 2.0 * beta2 - densityEffect(std::sqrt(bg2));
        if (!(bethe > 0.0))
            return d;
        const double f = std::clamp((std::log(Wmax / T0) - beta2 * (1.0 - T0 / Wmax)) / bethe, 0.0, 1.0);

        // Radii in pixels: 1/2 pixel (range of T0), range of Wmax, depth
        const double pxPerGcm2 = 1.0e4 / (m.density_gcm3 * m_params.pitch_um);
        d.r0   = 0.5;
        d.rMax = m_electronTable.range(Wmax) * pxPerGcm2;
        d.d    = m_params.depth_um / m_params.pitch_um;
        if (!(f > 0.0) || !(d.rMax > d.r0))
            return DeltaRays{};
        d.a    = std::clamp(d.d, d.r0, d.rMax);
        d.norm = std::log(d.a / d.r0) + d.d * (1.0 / d.a - 1.0 / d.rMax);

        const double kept = f * d.norm / std::log(d.rMax / d.r0);   // share of the loss the delta rays leave in the layer
        d.share = kept / (1.0 - f + kept);
        return d;
    }

    double CosmicDeposit::deltaRayShare(size_t i, double E) const
    {
        if (i >= m_spectrum.species().size())
            throw std::out_of_range("CosmicDeposit: no species " + std::to_string(i) + " in the spectrum.");
        return deltaRays(i, E).share;
    }

    CosmicDeposit::Kernel CosmicDeposit::deltaRayKernel(const DeltaRays& d, double electrons) const
    {
        Kernel k;
        if (!(d.share > 0.0) || !(electrons > 0.0))
            return k;

        // Charge per pixel at r: electrons w(r) / (2 pi norm r^2). Radius where it falls to deltaRayCut_e:
        // c / r^2 with w = 1 (r <= a), then c d / r^3 (r > a).
        const double c = electrons / (2.0 * kPi * d.norm * m_params.deltaRayCut_e);
        double R = std::sqrt(c);
        if (R > d.a)
            R = std::cbrt(c * d.d);
        R = std::min(R, d.rMax);
        if (!(R > d.r0))
            return k;   // below deltaRayCut_e everywhere: no halo

        // ∫ w dln r from r0 to r
        const auto W = [&d](double r) { return (r <= d.a) ? std::log(r / d.r0) : std::log(d.a / d.r0) + d.d * (1.0 / d.a - 1.0 / r); };

        k.R = static_cast<long>(std::floor(R + 0.5));
        const long n = 2 * k.R + 1;
        k.w.assign(static_cast<size_t>(n * n), 0.0);

        // Rings of width 0.25 min(1, r) pixel, each sampled every 1/4 pixel along its circle, in a multiple of 4
        // directions shifted by half a step: the kernel keeps the symmetries of the pixel grid.
        double sum = 0.0;
        for (double r1 = d.r0; r1 < R; )
        {
            const double r2     = std::min(R, r1 + 0.25 * std::min(1.0, r1));
            const double weight = W(r2) - W(r1);
            const double rm     = std::sqrt(r1 * r2);
            const size_t nPhi   = 4 * static_cast<size_t>(std::max(2.0, std::ceil(2.0 * kPi * rm)));
            for (size_t j = 0; j < nPhi; ++j)
            {
                const double phi = (static_cast<double>(j) + 0.5) * 2.0 * kPi / static_cast<double>(nPhi);
                const long a = static_cast<long>(std::floor(rm * std::cos(phi) + 0.5));   // pixel a covers [a - 1/2, a + 1/2)
                const long b = static_cast<long>(std::floor(rm * std::sin(phi) + 0.5));
                k.w[static_cast<size_t>((b + k.R) * n + a + k.R)] += weight / static_cast<double>(nPhi);
            }
            sum += weight;
            r1 = r2;
        }
        for (double& v : k.w) v /= sum;   // cut at R: the halo keeps its charge
        return k;
    }

    void CosmicDeposit::spread(long X, long Y, double electrons, const Kernel& k, const Geometry& g, std::vector<Deposit>& out) const
    {
        if (k.R == 0)
        {
            push(out, g, X, Y, electrons);
            return;
        }
        const long n = 2 * k.R + 1;
        for (long b = -k.R; b <= k.R; ++b)
            for (long a = -k.R; a <= k.R; ++a)
                push(out, g, X + a, Y + b, electrons * k.w[static_cast<size_t>((b + k.R) * n + a + k.R)]);
    }

    void CosmicDeposit::depositEvent(const Event& e, const Geometry& g, std::vector<Deposit>& out) const
    {
        // Parameters t in [0, 1] along the projected track where it crosses a pixel boundary.
        std::vector<double> t{0.0, 1.0};
        const auto crossings = [&t](double a0, double da)
        {
            if (da == 0.0) return;
            const double lo = std::min(a0, a0 + da), hi = std::max(a0, a0 + da);
            for (double k = std::floor(lo) + 1.0; k < hi; k += 1.0) t.push_back((k - a0) / da);
        };
        crossings(e.x, e.dx);
        crossings(e.y, e.dy);
        std::sort(t.begin(), t.end());

        // Charge along the path: the CSDA energy lost there, E(s_i) - E(s_i+1) with s = t * path (rises to a Bragg
        // peak for a particle that stops); uniform without a path (hit made by hand) or without a loss.
        const bool tabulated = e.path_um > 0.0 && e.species < m_tables.size() && !m_tables[e.species].empty();
        const StoppingTable* table = tabulated ? &m_tables[e.species] : nullptr;
        const double T      = tabulated ? m_spectrum.species()[e.species].kineticEnergy(e.energy) : 0.0;
        const double s_gcm2 = e.path_um * 1.0e-4 * m_params.material.density_gcm3;
        const double E0     = tabulated ? table->energyAfter(T, 0.0) : 0.0;
        const double lost   = tabulated ? E0 - table->energyAfter(T, s_gcm2) : 0.0;
        const bool   csda   = lost > 0.0;

        // Lateral spread of the species; delta rays: share h of the charge in their halo, around each pixel of
        // the track (the halo of one pixel holds h * electrons / track length)
        const Kernel&   kernel = (e.species < m_kernels.size()) ? m_kernels[e.species] : m_kernel;
        const DeltaRays delta  = deltaRays(e.species, e.energy);
        const Kernel    halo   = deltaRayKernel(delta, delta.share * e.electrons / std::max(1.0, std::hypot(e.dx, e.dy)));
        const double    h      = (halo.R > 0) ? delta.share : 0.0;   // none where the halo is below deltaRayCut_e

        for (size_t i = 0; i + 1 < t.size(); ++i)
        {
            if (!(t[i + 1] - t[i] > 0.0)) continue;
            const double weight = csda ? (table->energyAfter(T, t[i] * s_gcm2) - table->energyAfter(T, t[i + 1] * s_gcm2)) / lost
                                       : t[i + 1] - t[i];   // share of the charge in this pixel
            if (!(weight > 0.0)) continue;
            const double tm = 0.5 * (t[i] + t[i + 1]);
            const long   X  = static_cast<long>(std::floor(e.x + tm * e.dx));
            const long   Y  = static_cast<long>(std::floor(e.y + tm * e.dy));
            spread(X, Y, (1.0 - h) * e.electrons * weight, kernel, g, out);
            if (h > 0.0)
                spread(X, Y, h * e.electrons * weight, halo, g, out);
        }
    }

    void CosmicDeposit::generate(double dt, const Geometry& g, std::vector<Deposit>& out) const
    {
        generate(dt, g, out, nullptr);
    }

    void CosmicDeposit::generate(double dt, const Geometry& g, std::vector<Deposit>& out,
                                 std::vector<Event>* events) const
    {
        std::mt19937_64 rng = makeRng();
        const std::vector<Event> drawn = drawEvents(dt, g, rng);
        for (const Event& e : drawn)
            depositEvent(e, g, out);
        if (events)
            events->insert(events->end(), drawn.begin(), drawn.end());
    }

    const std::vector<CosmicDeposit::ColumnInfo>& CosmicDeposit::columnInfo()
    {
        static const std::vector<ColumnInfo> columns{
            {"SPECIES", "",         "Species of the particle"},
            {"Z",       "",         "Charge number"},
            {"A",       "",         "Mass number"},
            {"ENERGY",  "MeV",      "Kinetic energy, per nucleon for a nucleus, per particle otherwise"},
            {"MU",      "",         "cos(theta) of the incidence"},
            {"PHI",     "rad",      "Azimuth of the track on the matrix"},
            {"FACE",    "",         "Face crossed at entry: +1 front, -1 back"},
            {"X",       "pixel",    "Entry point, physical column (pixel i covers [i, i+1))"},
            {"Y",       "pixel",    "Entry point, physical row (pixel i covers [i, i+1))"},
            {"DX",      "pixel",    "Projected displacement of the path in the layer, columns"},
            {"DY",      "pixel",    "Projected displacement of the path in the layer, rows"},
            {"PATH",    "um",       "Path in the layer: depth / mu, or the range if the particle stops"},
            {"NE",      "electron", "Charge deposited in the layer"},
        };
        return columns;
    }

    CosmicDeposit::EventColumns CosmicDeposit::toColumns(const std::vector<Event>& events) const
    {
        EventColumns c;
        const size_t n = events.size();
        c.species.reserve(n);
        for (auto* v : {&c.Z, &c.A, &c.energy, &c.mu, &c.phi, &c.face, &c.x, &c.y, &c.dx, &c.dy, &c.path, &c.electrons})
            v->reserve(n);

        for (const Event& e : events)
        {
            const specie& s = m_spectrum.species().at(e.species);
            c.species.push_back(s.name());
            c.Z.push_back(s.Z());
            c.A.push_back(s.A());
            c.energy.push_back(e.energy);
            c.mu.push_back(e.mu);
            c.phi.push_back(e.phi);
            c.face.push_back(static_cast<double>(e.face));
            c.x.push_back(e.x);
            c.y.push_back(e.y);
            c.dx.push_back(e.dx);
            c.dy.push_back(e.dy);
            c.path.push_back(e.path_um);
            c.electrons.push_back(e.electrons);
        }
        return c;
    }

    std::vector<HeaderCard> CosmicDeposit::inputCards() const
    {
        // Geometric factor, weighted by the flux of the species
        const double flux = m_spectrum.integralIntensity();
        double G = 0.0;
        for (size_t i = 0; i < m_acceptance.size(); ++i)
            G += m_acceptance[i].G * m_spectrum.integralIntensity(i);
        G = (flux > 0.0) ? G / flux : 0.0;

        return {
            {"CRFLUX",   flux,                           "[m-2 s-1 sr-1] Cosmic ray integral intensity"},
            {"CRACC",    G,                              "[sr] Cosmic ray geometric factor, flux-weighted"},
            {"CRTILT",   m_params.tilt_deg,              "[deg] Tilt of the sensor normal from the zenith"},
            {"CRDEPTH",  m_params.depth_um,              "[um] Sensitive depth of the cosmic ray hits"},
            {"CRPITCH",  m_params.pitch_um,              "[um] Pixel pitch used for the cosmic ray hits"},
            {"CREVPE",   m_params.eVPerElectron,         "[eV] Energy per electron for the cosmic ray hits"},
            {"CRSIGD",   m_params.sigmaDiffusion_px,     "[pixel] Lateral spread of the cosmic ray charge"},
            {"CRMAXTRK", m_params.maxTrackPixels,        "[pixel] Longest cosmic ray track drawn"},
            {"CRDCUT",   m_params.deltaRayCut_e,         "[e-] Smallest charge per pixel of a delta-ray halo"},
            {"CRTRUNC",  truncatedFraction(),            "Share of the cosmic ray hits not drawn (longer tracks)"},
            {"CRSEED",   static_cast<size_t>(seed()),    "Seed of the cosmic ray draws, 0: none"},
        };
    }

#pragma endregion
#pragma endregion
    }
#pragma endregion

    std::shared_ptr<ChargeDeposit> defaultChargeDeposit()
    {
        return std::make_shared<CosmicRay::CosmicDeposit>();
    }
}