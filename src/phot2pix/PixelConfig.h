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
#include <string>
#include <unordered_map>
#include <variant>
#include <stdexcept>
#include <initializer_list>
#include <fstream>
#include <nlohmann/json.hpp>

namespace phot2pix
{
    using json = nlohmann::json;

    /**
     * @class PixelConfig
     * @brief Typed key-value configuration bag, JSON-compatible.
     *
     * Stores parameters as std::variant<bool, int{8..64}_t, uint{8..64}_t, size_t, float, double, std::string>.
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
        using Value = std::variant<bool, int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t, size_t, float, double, std::string, json>;

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

    private:
        std::unordered_map<std::string, Value> m_params;

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
