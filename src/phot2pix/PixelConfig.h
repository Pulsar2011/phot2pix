//
// Created by gillard on 2026-06-15.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// PixelConfig.h
//

#ifndef P2P_PIXEL_CONFIG_H
#define P2P_PIXEL_CONFIG_H

#include <cstdint>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <stdexcept>
#include <initializer_list>
#include <fstream>
#include <nlohmann/json.hpp>

namespace phot2pix
{
    using json = nlohmann::json;

    /// Value of a FITS header card: a number, a count, or a string.
    using CardValue = std::variant<double, size_t, std::string>;

    /// One FITS header card describing a simulation input: keyword (at most 8
    /// characters), value, and a comment that starts with the unit in square
    /// brackets, following the FITS convention ("[e-/ADU] Pixel gain").
    /// PixelConfig::card() builds the value of an input as the configuration gave
    /// it; DetectorSimulator records and writes the cards of its inputs.
    struct HeaderCard
    {
        std::string key;
        CardValue   value;
        std::string comment;
    };

    /**
     * @class PixelConfig
     * @brief Typed key-value configuration bag, JSON-compatible.
     *
     * Stores parameters as std::variant<bool, int{8..64}_t, uint{8..64}_t, size_t, float, double, std::string, json>
     * (size_t only where it differs from uint64_t, see Value).
     * Child-class configs simply call set() for their own keys before
     * passing the config up to the parent constructor.
     *
     * Usage:
     *   PixelConfig cfg = {{"dark", 0.001}, {"readNoise", 5.0}, {"gain", 1.5}};
     *   double rn = cfg.get<double>("readNoise");          // throws if missing
     *   double rn = cfg.get<double>("readNoise", 3.0);     // with default
     */
    class PixelConfig
    {
    public:
        // size_t is listed only when it is a type of its own. On macOS it is (unsigned long, while
        // uint64_t is unsigned long long); on 64-bit Linux both are unsigned long, and a variant
        // holding the same type twice rejects assignment, get<T>() and holds_alternative<T>() for it.
        using Value = std::conditional_t<std::is_same_v<size_t, uint64_t>,
            std::variant<bool, int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t,         float, double, std::string, json>,
            std::variant<bool, int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t, size_t, float, double, std::string, json>>;

        PixelConfig() = default;

        /// Brace-list construction: {{"dark", 0.001}, {"readNoise", 5.0}, ...}
        PixelConfig(std::initializer_list<std::pair<const std::string, Value>> params)
            : m_params(params)
        {}

        /// Build from an already-parsed nlohmann::json object.
        explicit PixelConfig(const json& j)
        {
            _fromJson(j);
        }

        /// Build from a JSON file path.
        explicit PixelConfig(const std::string& filepath)
        {
            std::ifstream file(filepath);
            if (!file.is_open())
                throw std::runtime_error("PixelConfig: cannot open '" + filepath + "'");
            _fromJson(json::parse(file));
        }

        // ------------------------------------------------------------------ //
        //  Accessors
        // ------------------------------------------------------------------ //

        /// Retrieve a value with the expected type T. Throws if key is absent
        /// or the stored type does not match T.
        template<typename T>
        T get(const std::string& key) const
        {
            auto it = m_params.find(key);
            if (it == m_params.end())
                throw std::out_of_range("PixelConfig: missing key '" + key + "'");
            if (!std::holds_alternative<T>(it->second))
                throw std::bad_variant_access{};
            return std::get<T>(it->second);
        }

        /// Retrieve a value, returning defaultVal if the key is absent.
        template<typename T>
        T get(const std::string& key, T defaultVal) const
        {
            auto it = m_params.find(key);
            if (it == m_params.end()) return defaultVal;
            if (!std::holds_alternative<T>(it->second))
                throw std::bad_variant_access{};
            return std::get<T>(it->second);
        }

        /// Store or overwrite a value.
        template<typename T>
        void set(const std::string& key, T value)
        {
            m_params[key] = Value(std::move(value));
        }

        bool has(const std::string& key) const
        {
            return m_params.count(key) != 0;
        }

        /// Merge entries from another PixelConfig (existing keys are NOT overwritten).
        void mergeDefaults(const PixelConfig& defaults)
        {
            for (const auto& [k, v] : defaults.m_params)
                m_params.emplace(k, v);   // emplace does nothing if key exists
        }

        /// Return the stored value cast to double, regardless of the stored numeric type.
        /// Works for all arithmetic types (double, float, int64_t, uint64_t, …).
        /// std::string is not convertible to double and causes std::invalid_argument.
        /// Throws std::out_of_range   if key is absent.
        /// Throws std::invalid_argument if the stored type is not numeric.
        double getAsDouble(const std::string& key) const
        {
            auto it = m_params.find(key);
            if (it == m_params.end())
                throw std::out_of_range("PixelConfig: missing key '" + key + "'");
            return std::visit([&key](const auto& v) -> double {
                using T = std::decay_t<decltype(v)>;
                if constexpr (std::is_convertible_v<T, double>)
                    return static_cast<double>(v);
                else
                    throw std::invalid_argument(
                        "PixelConfig: key '" + key + "' is not a numeric type.");
            }, it->second);
        }

        /// Return the raw json node for array/object keys. Throws if absent or not a json type.
        const json& getJson(const std::string& key) const
        {
            auto it = m_params.find(key);
            if (it == m_params.end())
                throw std::out_of_range("PixelConfig: missing key '" + key + "'");
            if (!std::holds_alternative<json>(it->second))
                throw std::invalid_argument("PixelConfig: key '" + key + "' is not a JSON array/object.");
            return std::get<json>(it->second);
        }

        /// Non-throwing overload: returns defaultVal if the key is absent or non-numeric.
        double getAsDouble(const std::string& key, double defaultVal) const noexcept
        {
            try { return getAsDouble(key); } catch (...) { return defaultVal; }
        }

        /// Serialize back to a nlohmann::json object.
        json toJson() const
        {
            json j;
            for (const auto& [k, v] : m_params)
            {
                std::visit([&](const auto& val){ j[k] = val; }, v);
            }
            return j;
        }

        /// Header-card value describing the input `key` as the configuration gave it,
        /// in the forms DetectorSimulator::initialize() accepts for a per-pixel parameter:
        ///   "key": number                        → the number
        ///   "key": "path" [+ "key.hdu": h]       → "file.fits" or "file.fits[h]"
        ///   "key": {"file": "path", "hdu": h}    → "file.fits[h]"
        ///   "key.file": "path" [+ "key.hdu": h]  → "file.fits[h]"
        ///   any other JSON value                 → its compact JSON text
        /// Only the file name is kept, not its directory: a FITS string value holds at
        /// most 68 characters. Returns nullopt when the key is absent.
        std::optional<CardValue> card(const std::string& key) const
        {
            const auto sibling = [this, &key]() -> json
            {
                const auto it = m_params.find(key + ".hdu");
                return (it == m_params.end()) ? json() : _toJson(it->second);
            };

            if (const auto it = m_params.find(key); it != m_params.end())
            {
                const Value& v = it->second;

                if (std::holds_alternative<std::string>(v))
                {
                    const std::string& s = std::get<std::string>(v);
                    if (!s.empty() && s.front() == '{')   // JSON-encoded object
                    {
                        try { return _fileCard(json::parse(s), s); }
                        catch (const json::exception&) { return s; }
                    }
                    return _fileName(s) + _hduSuffix(sibling());
                }

                if (std::holds_alternative<json>(v))
                {
                    const json& j = std::get<json>(v);
                    if (j.is_number())
                        return j.get<double>();
                    return _fileCard(j, j.dump());
                }

                return getAsDouble(key);   // any numeric type
            }

            if (const auto it = m_params.find(key + ".file"); it != m_params.end() && std::holds_alternative<std::string>(it->second))
                return _fileName(std::get<std::string>(it->second)) + _hduSuffix(sibling());

            return std::nullopt;
        }

    private:
        std::unordered_map<std::string, Value> m_params;

        static json _toJson(const Value& v)
        {
            return std::visit([](const auto& x) { return json(x); }, v);
        }

        static std::string _fileName(const std::string& path)
        {
            return std::filesystem::path(path).filename().string();
        }

        // "[h]" for an HDU given as an index or an EXTNAME, "" when absent.
        static std::string _hduSuffix(const json& hdu)
        {
            if (hdu.is_number())
                return "[" + std::to_string(hdu.get<long long>()) + "]";
            if (hdu.is_string())
                return "[" + hdu.get<std::string>() + "]";
            return "";
        }

        // {"file": path, "hdu": h} → "file.fits[h]"; anything else → fallback.
        static std::string _fileCard(const json& j, const std::string& fallback)
        {
            if (j.is_object() && j.contains("file") && j["file"].is_string())
                return _fileName(j["file"].get<std::string>()) + _hduSuffix(j.contains("hdu") ? j["hdu"] : json());
            return fallback;
        }

        void _fromJson(const json& j)
        {
            for (const auto& [k, v] : j.items())
            {
                if      (v.is_boolean())         m_params[k] = v.get<bool>();
                else if (v.is_number_unsigned())  m_params[k] = v.get<uint64_t>();
                else if (v.is_number_integer())   m_params[k] = v.get<int64_t>();
                else if (v.is_number_float())     m_params[k] = v.get<double>();
                else if (v.is_string())           m_params[k] = v.get<std::string>();
                else if (v.is_array() || v.is_object()) m_params[k] = v;
            }
        }
    };

} // namespace phot2pix

#endif // P2P_PIXEL_CONFIG_H
