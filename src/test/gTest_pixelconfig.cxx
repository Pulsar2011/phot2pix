//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_pixelconfig.cxx
//

#include <gtest/gtest.h>
#include <phot2pix/PixelConfig.h>
#include <filesystem>
#include <fstream>

using namespace phot2pix;
namespace fs = std::filesystem;

// ══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — construction
// ══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, DefaultConstructor)
{
    PixelConfig cfg;
    EXPECT_FALSE(cfg.has("anything"));
}

TEST(PixelConfig, InitializerList)
{
    PixelConfig cfg{
        {"gain",      2.5},
        {"readNoise", 5.0},
        {"label",     std::string("sensor")},
    };
    EXPECT_TRUE (cfg.has("gain"));
    EXPECT_TRUE (cfg.has("readNoise"));
    EXPECT_TRUE (cfg.has("label"));
    EXPECT_FALSE(cfg.has("dark"));
}

TEST(PixelConfig, FromJson)
{
    // bool → bool, unsigned → uint64_t, signed → int64_t, float → double, string → string.
    // Nested objects and arrays are stored as json and accessible via getJson().
    json j;
    j["flag"]   = true;
    j["uint_v"] = uint64_t{7};
    j["int_v"]  = int64_t{-3};
    j["dbl_v"]  = 1.5;
    j["str_v"]  = "hello";
    j["nested"] = { {"a", 1} };

    PixelConfig cfg(j);
    EXPECT_TRUE (cfg.get<bool>       ("flag"));
    EXPECT_EQ   (cfg.get<uint64_t>   ("uint_v"), uint64_t{7});
    EXPECT_EQ   (cfg.get<int64_t>    ("int_v"),  int64_t{-3});
    EXPECT_DOUBLE_EQ(cfg.get<double> ("dbl_v"),  1.5);
    EXPECT_EQ   (cfg.get<std::string>("str_v"),  "hello");
    EXPECT_TRUE (cfg.has("nested"));
}

TEST(PixelConfig, FromJsonFile)
{
    fs::path tmp = fs::temp_directory_path() / "p2p_test_pixelconfig.json";
    {
        std::ofstream f(tmp);
        f << R"({"gain":2.0,"dark":0.001,"label":"sensor"})";
    }

    PixelConfig cfg(tmp.string());
    EXPECT_DOUBLE_EQ(cfg.get<double>     ("gain"), 2.0);
    EXPECT_DOUBLE_EQ(cfg.get<double>     ("dark"), 0.001);
    EXPECT_EQ       (cfg.get<std::string>("label"), "sensor");

    fs::remove(tmp);
}

TEST(PixelConfig, FromJsonFileMissing)
{
    EXPECT_THROW(PixelConfig(std::string("/nonexistent/path/file.json")), std::runtime_error);
}

// ══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — has()
// ══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, HasPresentKey)
{
    PixelConfig cfg{{"x", 1.0}};
    EXPECT_TRUE(cfg.has("x"));
}

TEST(PixelConfig, HasAbsentKey)
{
    PixelConfig cfg{{"x", 1.0}};
    EXPECT_FALSE(cfg.has("y"));
}

// ══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — get<T>(key)
// ══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, GetTypedAllVariants)
{
    PixelConfig cfg;
    cfg.set<bool>       ("b",   true);
    cfg.set<int8_t>     ("i8",  int8_t{-1});
    cfg.set<int16_t>    ("i16", int16_t{-2});
    cfg.set<int32_t>    ("i32", int32_t{-3});
    cfg.set<int64_t>    ("i64", int64_t{-4});
    cfg.set<uint8_t>    ("u8",  uint8_t{5});
    cfg.set<uint16_t>   ("u16", uint16_t{6});
    cfg.set<uint32_t>   ("u32", uint32_t{7});
    cfg.set<uint64_t>   ("u64", uint64_t{8});
    cfg.set<size_t>     ("sz",  size_t{9});
    cfg.set<float>      ("f",   1.5f);
    cfg.set<double>     ("d",   2.5);
    cfg.set<std::string>("s",   std::string("hello"));

    EXPECT_EQ       (cfg.get<bool>       ("b"),   true);
    EXPECT_EQ       (cfg.get<int8_t>     ("i8"),  int8_t{-1});
    EXPECT_EQ       (cfg.get<int16_t>    ("i16"), int16_t{-2});
    EXPECT_EQ       (cfg.get<int32_t>    ("i32"), int32_t{-3});
    EXPECT_EQ       (cfg.get<int64_t>    ("i64"), int64_t{-4});
    EXPECT_EQ       (cfg.get<uint8_t>    ("u8"),  uint8_t{5});
    EXPECT_EQ       (cfg.get<uint16_t>   ("u16"), uint16_t{6});
    EXPECT_EQ       (cfg.get<uint32_t>   ("u32"), uint32_t{7});
    EXPECT_EQ       (cfg.get<uint64_t>   ("u64"), uint64_t{8});
    EXPECT_EQ       (cfg.get<size_t>     ("sz"),  size_t{9});
    EXPECT_FLOAT_EQ (cfg.get<float>      ("f"),   1.5f);
    EXPECT_DOUBLE_EQ(cfg.get<double>     ("d"),   2.5);
    EXPECT_EQ       (cfg.get<std::string>("s"),   "hello");
}

TEST(PixelConfig, GetTypedMissingKeyThrows)
{
    PixelConfig cfg;
    EXPECT_THROW(cfg.get<double>("absent"), std::out_of_range);
}

TEST(PixelConfig, GetTypedWrongTypeThrows)
{
    PixelConfig cfg{{"gain", 2.5}};           // stored as double
    EXPECT_THROW(cfg.get<int32_t>("gain"), std::bad_variant_access);
}

// ══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — get<T>(key, default)
// ══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, GetWithDefaultPresentKey)
{
    PixelConfig cfg{{"gain", 2.5}};
    EXPECT_DOUBLE_EQ(cfg.get<double>("gain", 1.0), 2.5);
}

TEST(PixelConfig, GetWithDefaultAbsentKey)
{
    PixelConfig cfg;
    EXPECT_DOUBLE_EQ(cfg.get<double>("gain", 1.0), 1.0);
}

TEST(PixelConfig, GetWithDefaultWrongTypeThrows)
{
    // Key exists but type mismatches — even with a default, it throws.
    PixelConfig cfg{{"gain", 2.5}};           // stored as double
    EXPECT_THROW(cfg.get<int32_t>("gain", int32_t{0}), std::bad_variant_access);
}

// ══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — set<T>()
// ══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, SetInsert)
{
    PixelConfig cfg;
    cfg.set<double>("x", 3.14);
    EXPECT_TRUE(cfg.has("x"));
    EXPECT_DOUBLE_EQ(cfg.get<double>("x"), 3.14);
}

TEST(PixelConfig, SetOverwriteSameType)
{
    PixelConfig cfg{{"x", 1.0}};
    cfg.set<double>("x", 9.9);
    EXPECT_DOUBLE_EQ(cfg.get<double>("x"), 9.9);
}

TEST(PixelConfig, SetOverwriteDifferentType)
{
    // Changing the stored type on an existing key is legal.
    PixelConfig cfg;
    cfg.set<double>     ("key", 1.0);
    cfg.set<std::string>("key", std::string("now a string"));
    EXPECT_EQ   (cfg.get<std::string>("key"), "now a string");
    EXPECT_THROW(cfg.get<double>("key"), std::bad_variant_access);
}

// ══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — mergeDefaults()
// ══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, MergeDefaultsAddsNewKeys)
{
    PixelConfig cfg {{"gain", 2.0}};
    PixelConfig defs{{"gain", 1.0}, {"dark", 0.001}};
    cfg.mergeDefaults(defs);

    EXPECT_DOUBLE_EQ(cfg.get<double>("gain"), 2.0);   // not overwritten
    EXPECT_DOUBLE_EQ(cfg.get<double>("dark"), 0.001); // inserted
}

TEST(PixelConfig, MergeDefaultsDoesNotOverwrite)
{
    PixelConfig cfg {{"a", 10.0}, {"b", 20.0}};
    PixelConfig defs{{"a", 99.0}, {"b", 99.0}, {"c", 99.0}};
    cfg.mergeDefaults(defs);

    EXPECT_DOUBLE_EQ(cfg.get<double>("a"), 10.0);
    EXPECT_DOUBLE_EQ(cfg.get<double>("b"), 20.0);
    EXPECT_DOUBLE_EQ(cfg.get<double>("c"), 99.0);
}

TEST(PixelConfig, MergeDefaultsIntoEmpty)
{
    PixelConfig cfg;
    PixelConfig defs{{"x", 5.0}, {"y", std::string("hello")}};
    cfg.mergeDefaults(defs);

    EXPECT_DOUBLE_EQ(cfg.get<double>     ("x"), 5.0);
    EXPECT_EQ       (cfg.get<std::string>("y"), "hello");
}

// ══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — getJson()
// ══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, GetJsonArray)
{
    json j;
    j["persistence"] = json::array({
        {{"f", 0.01},  {"t", 200.0}},
        {{"f", 0.005}, {"t", 1000.0}},
    });

    PixelConfig cfg(j);
    EXPECT_TRUE(cfg.has("persistence"));

    const json& traps = cfg.getJson("persistence");
    ASSERT_EQ(traps.size(), 2u);
    EXPECT_DOUBLE_EQ(traps[0].at("f").get<double>(), 0.01);
    EXPECT_DOUBLE_EQ(traps[0].at("t").get<double>(), 200.0);
    EXPECT_DOUBLE_EQ(traps[1].at("f").get<double>(), 0.005);
    EXPECT_DOUBLE_EQ(traps[1].at("t").get<double>(), 1000.0);
}

TEST(PixelConfig, GetJsonNestedObject)
{
    json j;
    j["nested"] = {{"a", 1}, {"b", 2}};

    PixelConfig cfg(j);
    EXPECT_TRUE(cfg.has("nested"));

    const json& obj = cfg.getJson("nested");
    EXPECT_EQ(obj.at("a").get<int>(), 1);
    EXPECT_EQ(obj.at("b").get<int>(), 2);
}

TEST(PixelConfig, GetJsonMissingKeyThrows)
{
    PixelConfig cfg;
    EXPECT_THROW(cfg.getJson("absent"), std::out_of_range);
}

TEST(PixelConfig, GetJsonOnScalarThrows)
{
    PixelConfig cfg{{"gain", 2.5}};
    EXPECT_THROW(cfg.getJson("gain"), std::invalid_argument);
}

TEST(PixelConfig, GetJsonFromFile)
{
    fs::path tmp = fs::temp_directory_path() / "p2p_test_pixelconfig_json.json";
    {
        std::ofstream f(tmp);
        f << R"({"gain":2.0,"persistence":[{"f":0.01,"t":200.0},{"f":0.005,"t":1000.0}]})";
    }

    PixelConfig cfg(tmp.string());
    EXPECT_DOUBLE_EQ(cfg.get<double>("gain"), 2.0);

    const json& traps = cfg.getJson("persistence");
    ASSERT_EQ(traps.size(), 2u);
    EXPECT_DOUBLE_EQ(traps[0].at("f").get<double>(), 0.01);
    EXPECT_DOUBLE_EQ(traps[1].at("t").get<double>(), 1000.0);

    fs::remove(tmp);
}

TEST(PixelConfig, GetJsonEmptyArray)
{
    json j;
    j["persistence"] = json::array();

    PixelConfig cfg(j);
    EXPECT_TRUE(cfg.has("persistence"));
    EXPECT_TRUE(cfg.getJson("persistence").empty());
}

// ══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — toJson() round-trip
// ══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, ToJsonRoundTrip)
{
    PixelConfig original;
    original.set<double>     ("gain",  2.5);
    original.set<int64_t>    ("count", int64_t{10});
    original.set<std::string>("label", std::string("test"));
    original.set<bool>       ("flag",  true);

    PixelConfig rebuilt(original.toJson());

    // After JSON serialisation int64_t re-parses as int64_t (signed integer path).
    EXPECT_DOUBLE_EQ(rebuilt.get<double>     ("gain"),  2.5);
    EXPECT_EQ       (rebuilt.get<int64_t>    ("count"), int64_t{10});
    EXPECT_EQ       (rebuilt.get<std::string>("label"), "test");
    EXPECT_TRUE     (rebuilt.get<bool>       ("flag"));
}

TEST(PixelConfig, ToJsonContainsAllKeys)
{
    PixelConfig cfg{{"a", 1.0}, {"b", std::string("x")}};
    json j = cfg.toJson();
    EXPECT_TRUE(j.contains("a"));
    EXPECT_TRUE(j.contains("b"));
    EXPECT_DOUBLE_EQ(j["a"].get<double>(), 1.0);
    EXPECT_EQ       (j["b"].get<std::string>(), "x");
}

// ══════════════════════════════════════════════════════════════════════════════
//  PixelConfig — card(): header-card value of an input
// ══════════════════════════════════════════════════════════════════════════════

TEST(PixelConfig, CardOfAbsentKeyIsEmpty)
{
    PixelConfig cfg;
    EXPECT_FALSE(cfg.card("dark").has_value());
}

TEST(PixelConfig, CardOfANumberIsTheNumber)
{
    PixelConfig cfg{{"dark", 0.01}, {"nPoles", 16}};
    ASSERT_TRUE(cfg.card("dark").has_value());
    EXPECT_DOUBLE_EQ(std::get<double>(*cfg.card("dark")),   0.01);
    EXPECT_DOUBLE_EQ(std::get<double>(*cfg.card("nPoles")), 16.0);
}

TEST(PixelConfig, CardOfAMapIsItsFileNameAndHdu)
{
    // Plain path: only the file name is kept.
    PixelConfig path{{"dark", std::string("/data/calib/dark_map.fits")}};
    EXPECT_EQ(std::get<std::string>(*path.card("dark")), "dark_map.fits");

    // Path + "key.hdu" (index or EXTNAME).
    PixelConfig withIndex{{"dark", std::string("/data/dark.fits")}, {"dark.hdu", 2}};
    EXPECT_EQ(std::get<std::string>(*withIndex.card("dark")), "dark.fits[2]");

    PixelConfig withName{{"gain", std::string("gain.fits")}, {"gain.hdu", std::string("GAIN")}};
    EXPECT_EQ(std::get<std::string>(*withName.card("gain")), "gain.fits[GAIN]");

    // Object form {"file", "hdu"}.
    PixelConfig object;
    object.set("QE", json{{"file", "/maps/qe.fits"}, {"hdu", "QE"}});
    EXPECT_EQ(std::get<std::string>(*object.card("QE")), "qe.fits[QE]");

    // Flat form "key.file" + "key.hdu".
    PixelConfig flat{{"readNoise.file", std::string("rn.fits")}, {"readNoise.hdu", 1}};
    EXPECT_EQ(std::get<std::string>(*flat.card("readNoise")), "rn.fits[1]");
}

TEST(PixelConfig, CardOfOtherJsonIsItsText)
{
    PixelConfig cfg;
    cfg.set("macc", json::array({4, 16, 4}));
    EXPECT_EQ(std::get<std::string>(*cfg.card("macc")), "[4,16,4]");
}
