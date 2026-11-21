//
// Created by gillard on 2026-06-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_detector.cxx
//

#include <filesystem>

#include <gtest/gtest.h>
#include <phot2pix/DetectorSimulator.h>
//#include <cfitsio/fitsio.h>

using namespace phot2pix;

// ═══════════════════════════════════════════════════════════════════════════════
//  DetectorSimulator — construction and geometry
// ═══════════════════════════════════════════════════════════════════════════════

TEST(DetectorSimulator, Construction)
{
    DetectorSimulator<> det(64, 32);
    EXPECT_EQ(det.getWidth(),  64u);
    EXPECT_EQ(det.getHeight(), 32u);
    EXPECT_EQ(det.size(),      64u * 32u);
}

TEST(DetectorSimulator, InvalidConstruction)
{
    EXPECT_THROW(DetectorSimulator<>(0,  32), std::invalid_argument); // zero width
    EXPECT_THROW(DetectorSimulator<>(64,  0), std::invalid_argument); // zero height
    EXPECT_THROW(DetectorSimulator<>(0,   0), std::invalid_argument); // both zero
}

// ── TestDet: promotes protected set*Map() to public for unit testing ─────────
// In production, callers use initialize(PixelConfig). This shim lets tests
// exercise the underlying per-pixel setter logic directly.
class TestDet : public DetectorSimulator<>
{
public:
    TestDet(size_t w, size_t h) : DetectorSimulator<>(w, h) {}
    using DetectorSimulator<>::setGainMap;
    using DetectorSimulator<>::setDarkCurrentMap;
    using DetectorSimulator<>::setReadNoiseMap;
    using DetectorSimulator<>::setQEMap;
    using DetectorSimulator<>::setBaselineMap;
    using DetectorSimulator<>::setKTCnoiseMap;
};

// ═══════════════════════════════════════════════════════════════════════════════
//  DetectorSimulator — construction of invalide datatype
// ═══════════════════════════════════════════════════════════════════════════════
TEST(DetectorSimulator, InvalidDataType)
{
    DetectorSimulator<> det(64, 32, 999); // invalid datatype
    EXPECT_THROW(det.simulate(0.0, 1.0), std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  DetectorSimulator — construction of different datatype
// ═══════════════════════════════════════════════════════════════════════════════
TEST(DetectorSimulator, Int8)
{
    const size_t W = 16;
    const size_t H = 8;
    DetectorSimulator<> det(W, H, SBYTE_IMG);
    std::valarray<double> flux(0.0, W * H);

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),SBYTE_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<int8_t> data = *result->GetData<int8_t>();
    EXPECT_EQ(data.size(), W * H);
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 0);
}

TEST(DetectorSimulator, UInt8)
{
    const size_t W = 16;
    const size_t H = 8;
    DetectorSimulator<> det(W, H, BYTE_IMG);
    std::valarray<double> flux(0.0, W * H);

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),BYTE_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<uint8_t> data = *result->GetData<uint8_t>();
    EXPECT_EQ(data.size(), W * H);
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 0);
}

TEST(DetectorSimulator, Int16)
{
    const size_t W = 16;
    const size_t H = 8;
    DetectorSimulator<> det(W, H, SHORT_IMG);
    std::valarray<double> flux(0.0, W * H);

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),SHORT_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<int16_t> data = *result->GetData<int16_t>();
    EXPECT_EQ(data.size(), W * H);
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 0);
}

TEST(DetectorSimulator, UInt16)
{
    const size_t W = 16;
    const size_t H = 8;
    DetectorSimulator<> det(W, H, USHORT_IMG);
    std::valarray<double> flux(0.0, W * H);

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),USHORT_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<uint16_t> data = *result->GetData<uint16_t>();
    EXPECT_EQ(data.size(), W * H);
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 0);
}

TEST(DetectorSimulator, Int32)
{
    const size_t W = 16;
    const size_t H = 8;
    DetectorSimulator<> det(W, H, LONG_IMG);
    std::valarray<double> flux(0.0, W * H);

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),LONG_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<int32_t> data = *result->GetData<int32_t>();
    EXPECT_EQ(data.size(), W * H);
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 0);
}

TEST(DetectorSimulator, UInt32)
{
    const size_t W = 16;
    const size_t H = 8;
    DetectorSimulator<> det(W, H, ULONG_IMG);
    std::valarray<double> flux(0.0, W * H);

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),ULONG_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<uint32_t> data = *result->GetData<uint32_t>();
    EXPECT_EQ(data.size(), W * H);
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 0);
}

TEST(DetectorSimulator, Int64)
{
    const size_t W = 16;
    const size_t H = 8;
    DetectorSimulator<> det(W, H, LONGLONG_IMG);
    std::valarray<double> flux(0.0, W * H);

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),LONGLONG_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<int64_t> data = *result->GetData<int64_t>();
    EXPECT_EQ(data.size(), W * H);
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 0);
}

TEST(DetectorSimulator, UInt64)
{
    const size_t W = 16;
    const size_t H = 8;
    DetectorSimulator<> det(W, H, ULONGLONG_IMG);
    std::valarray<double> flux(0.0, W * H);

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),ULONGLONG_IMG);


    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<uint64_t> data = *result->GetData<uint64_t>();
    EXPECT_EQ(data.size(), W * H);
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 0);
}

TEST(DetectorSimulator, Float)
{
    const size_t W = 16;
    const size_t H = 8;
    DetectorSimulator<> det(W, H, FLOAT_IMG);
    std::valarray<double> flux(0.0, W * H);

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),FLOAT_IMG);


    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<float> data = *result->GetData<float>();
    EXPECT_EQ(data.size(), W * H);
    for (size_t i = 0; i < data.size(); ++i)
        EXPECT_EQ(data[i], 0);
}



// ═══════════════════════════════════════════════════════════════════════════════
//  DetectorSimulator — construction with pixsim derived classe
// ═══════════════════════════════════════════════════════════════════════════════

class SimplePixel : public pixsim
{
public:
    SimplePixel() : pixsim(0.0, 0.0, 1.0, 1.0, 0.0, 0.0) {}
    
    double pixelresponse(const double& flux, const double& exposureTime) const override
    {
        return flux; // simple linear response
    }
};

TEST(SimplePixelSimulator, Int8)
{
    const size_t W = 23;
    const size_t H = 23;
    DetectorSimulator<SimplePixel> det(W, H, SBYTE_IMG);
    std::valarray<double> flux(0.0, W * H);
    for(size_t i = 0; i < flux.size(); ++i)
        flux[i] = static_cast<double>(i) - static_cast<double>(flux.size())/2.0; // some negative values

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),SBYTE_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<int8_t> data = *result->GetData<int8_t>();
    EXPECT_EQ(data.size(), W * H);

    for (size_t i = 0; i < flux.size(); ++i)
    {
        if(flux[i]+0.5 < static_cast<double>(std::numeric_limits<int8_t>::min()))
            EXPECT_EQ(data[i], std::numeric_limits<int8_t>::min());
        else if(flux[i]+0.5 > static_cast<double>(std::numeric_limits<int8_t>::max()))
            EXPECT_EQ(data[i], std::numeric_limits<int8_t>::max());
        else
            EXPECT_EQ(data[i], static_cast<int8_t>(flux[i]+0.5));
    }
}

TEST(SimplePixelSimulator, UInt8)
{
    const size_t W = 23;
    const size_t H = 23;
    DetectorSimulator<SimplePixel> det(W, H, BYTE_IMG);
    std::valarray<double> flux(0.0, W * H);
    for(size_t i = 0; i < flux.size(); ++i)
        flux[i] = static_cast<double>(i) - 50.0; // some negative values

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),BYTE_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<uint8_t> data = *result->GetData<uint8_t>();
    EXPECT_EQ(data.size(), W * H);

    for (size_t i = 0; i < flux.size(); ++i)
    {
        if(flux[i]+0.5 < 0.)
            EXPECT_EQ(data[i], 0);
        else if(flux[i]+0.5 > static_cast<double>(std::numeric_limits<uint8_t>::max()))
            EXPECT_EQ(data[i], std::numeric_limits<uint8_t>::max());
        else
            EXPECT_EQ(data[i], static_cast<uint8_t>(flux[i]+0.5));
    }
}

TEST(SimplePixelSimulator, Int16)
{
    const size_t W = 300;
    const size_t H = 300;
    DetectorSimulator<SimplePixel> det(W, H, SHORT_IMG);
    std::valarray<double> flux(0.0, W * H);
    for(size_t i = 0; i < flux.size(); ++i)
        flux[i] = static_cast<double>(i) - static_cast<double>(flux.size())/2.0; // some negative values

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),SHORT_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<int16_t> data = *result->GetData<int16_t>();
    EXPECT_EQ(data.size(), W * H);

    for (size_t i = 0; i < flux.size(); ++i)
    {
        if(flux[i]+0.5 < static_cast<double>(std::numeric_limits<int16_t>::min()))
            EXPECT_EQ(data[i], std::numeric_limits<int16_t>::min());
        else if(flux[i]+0.5 > static_cast<double>(std::numeric_limits<int16_t>::max()))
            EXPECT_EQ(data[i], std::numeric_limits<int16_t>::max());
        else
            EXPECT_EQ(data[i], static_cast<int16_t>(flux[i]+0.5));
    }
}

TEST(SimplePixelSimulator, UInt16)
{
    const size_t W = 23;
    const size_t H = 23;
    DetectorSimulator<SimplePixel> det(W, H, USHORT_IMG);
    std::valarray<double> flux(0.0, W * H);
    for(size_t i = 0; i < flux.size(); ++i)
        flux[i] = static_cast<double>(i) - 50.0; // some negative values

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),USHORT_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<uint16_t> data = *result->GetData<uint16_t>();
    EXPECT_EQ(data.size(), W * H);

    for (size_t i = 0; i < flux.size(); ++i)
    {
        if(flux[i]+0.5 < 0.)
            EXPECT_EQ(data[i], 0);
        else if(flux[i]+0.5 > static_cast<double>(std::numeric_limits<uint16_t>::max()))
            EXPECT_EQ(data[i], std::numeric_limits<uint16_t>::max());
        else
            EXPECT_EQ(data[i], static_cast<uint16_t>(flux[i]+0.5));
    }
}

TEST(SimplePixelSimulator, Int32)
{
    const size_t W = 300;
    const size_t H = 300;
    DetectorSimulator<SimplePixel> det(W, H, LONG_IMG);
    std::valarray<double> flux(0.0, W * H);
    for(size_t i = 0; i < flux.size(); ++i)
        flux[i] = static_cast<double>(i) - static_cast<double>(flux.size())/2.0; // some negative values

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),LONG_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<int32_t> data = *result->GetData<int32_t>();
    EXPECT_EQ(data.size(), W * H);

    for (size_t i = 0; i < flux.size(); ++i)
    {
        if(flux[i]+0.5 < static_cast<double>(std::numeric_limits<int32_t>::min()))
            EXPECT_EQ(data[i], std::numeric_limits<int32_t>::min());
        else if(flux[i]+0.5 > static_cast<double>(std::numeric_limits<int32_t>::max()))
            EXPECT_EQ(data[i], std::numeric_limits<int32_t>::max());
        else
            EXPECT_EQ(data[i], static_cast<int32_t>(flux[i]+0.5));
    }
}

TEST(SimplePixelSimulator, UInt32)
{
    const size_t W = 23;
    const size_t H = 23;
    DetectorSimulator<SimplePixel> det(W, H, ULONG_IMG);
    std::valarray<double> flux(0.0, W * H);
    for(size_t i = 0; i < flux.size(); ++i)
        flux[i] = static_cast<double>(i) - 50.0; // some negative values

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),ULONG_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<uint32_t> data = *result->GetData<uint32_t>();
    EXPECT_EQ(data.size(), W * H);

    for (size_t i = 0; i < flux.size(); ++i)
    {
        if(flux[i]+0.5 < 0.)
            EXPECT_EQ(data[i], 0);
        else if(flux[i]+0.5 > static_cast<double>(std::numeric_limits<uint32_t>::max()))
            EXPECT_EQ(data[i], std::numeric_limits<uint32_t>::max());
        else
            EXPECT_EQ(data[i], static_cast<uint32_t>(flux[i]+0.5));
    }
}

TEST(SimplePixelSimulator, Int64)
{
    const size_t W = 300;
    const size_t H = 300;
    DetectorSimulator<SimplePixel> det(W, H, LONGLONG_IMG);
    std::valarray<double> flux(0.0, W * H);
    for(size_t i = 0; i < flux.size(); ++i)
        flux[i] = static_cast<double>(i) - static_cast<double>(flux.size())/2.0; // some negative values

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),LONGLONG_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<int64_t> data = *result->GetData<int64_t>();
    EXPECT_EQ(data.size(), W * H);

    for (size_t i = 0; i < flux.size(); ++i)
    {
        if(flux[i]+0.5 < static_cast<double>(std::numeric_limits<int64_t>::min()))
            EXPECT_EQ(data[i], std::numeric_limits<int64_t>::min());
        else if(flux[i]+0.5 > static_cast<double>(std::numeric_limits<int64_t>::max()))
            EXPECT_EQ(data[i], std::numeric_limits<int64_t>::max());
        else
            EXPECT_EQ(data[i], static_cast<int64_t>(flux[i]+0.5));
    }
}

TEST(SimplePixelSimulator, UInt64)
{
    const size_t W = 23;
    const size_t H = 23;
    DetectorSimulator<SimplePixel> det(W, H, ULONGLONG_IMG);
    std::valarray<double> flux(0.0, W * H);
    for(size_t i = 0; i < flux.size(); ++i)
        flux[i] = static_cast<double>(i) - 50.0; // some negative values

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),ULONGLONG_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<uint64_t> data = *result->GetData<uint64_t>();
    EXPECT_EQ(data.size(), W * H);

    for (size_t i = 0; i < flux.size(); ++i)
    {
        if(flux[i]+0.5 < 0.)
            EXPECT_EQ(data[i], 0);
        else if(flux[i]+0.5 > static_cast<double>(std::numeric_limits<uint64_t>::max()))
            EXPECT_EQ(data[i], std::numeric_limits<uint64_t>::max());
        else
            EXPECT_EQ(data[i], static_cast<uint64_t>(flux[i]+0.5));
    }
}

TEST(SimplePixelSimulator, Float)
{
    const size_t W = 23;
    const size_t H = 23;
    DetectorSimulator<SimplePixel> det(W, H, FLOAT_IMG);
    std::valarray<double> flux(0.0, W * H);
    for(size_t i = 0; i < flux.size(); ++i)
        flux[i] = static_cast<double>(i) - W*H/2.0; // some negative values

    auto result = det.simulate(flux, 1.0);
    EXPECT_EQ(result->GetEqBitPerPixel(),FLOAT_IMG);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);

    std::valarray<float> data = *result->GetData<float>();
    EXPECT_EQ(data.size(), W * H);

    for (size_t i = 0; i < flux.size(); ++i)
    {
        EXPECT_NEAR(data[i], static_cast<float>(flux[i]), std::numeric_limits<float>::epsilon());
    }
}


// ═══════════════════════════════════════════════════════════════════════════════
//  Batch configuration
// ═══════════════════════════════════════════════════════════════════════════════

TEST(DetectorSimulator, Fill)
{
    DetectorSimulator<> det(8, 8);
    pixsim proto(0.5, 3.0, 2.0, 0.8, 1.0, 10.0);
    det.fill(proto);

    const auto& cdet = det;
    for (size_t i = 0; i < det.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(cdet.at(i).getDarkCurrent(),       0.5);
        EXPECT_DOUBLE_EQ(cdet.at(i).getReadNoise(),         3.0);
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(),              2.0);
        EXPECT_DOUBLE_EQ(cdet.at(i).getQuantumEfficiency(), 0.8);
        EXPECT_DOUBLE_EQ(cdet.at(i).getKTCnoise(),          1.0);
        EXPECT_DOUBLE_EQ(cdet.at(i).getBaseline(),         10.0);
    }
}

TEST(DetectorSimulator, MapSizeMismatch)
{
    TestDet det(8, 8); // 64 pixels
    std::vector<double> wrong(63, 1.0); // one too few

    EXPECT_THROW(det.setGainMap(wrong),        std::invalid_argument);
    EXPECT_THROW(det.setDarkCurrentMap(wrong), std::invalid_argument);
    EXPECT_THROW(det.setReadNoiseMap(wrong),   std::invalid_argument);
    EXPECT_THROW(det.setQEMap(wrong),          std::invalid_argument);
    EXPECT_THROW(det.setBaselineMap(wrong),    std::invalid_argument);
    EXPECT_THROW(det.setKTCnoiseMap(wrong),    std::invalid_argument);
}

TEST(DetectorSimulator, SetGainMap)
{
    TestDet det(4, 4); // 16 pixels
    std::vector<double> gains(16);
    for (size_t i = 0; i < 16; ++i) gains[i] = 1.0 + 0.1 * static_cast<double>(i);

    det.setGainMap(gains);
    const auto& cdet = det;

    for (size_t i = 0; i < 16; ++i)
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(), gains[i]);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Read-only pixel access
// ═══════════════════════════════════════════════════════════════════════════════

TEST(DetectorSimulator, ConstAccess)
{
    DetectorSimulator<> det(4, 4);
    pixsim proto(0.0, 0.0, 2.0, 1.0, 0.0, 0.0);
    det.fill(proto);
    const auto& cdet = det;

    // 1D index
    EXPECT_DOUBLE_EQ(cdet.at(5).getGain(), 2.0);

    // 2D index: col=1, row=1  →  flat index = 1 + 1*4 = 5
    EXPECT_DOUBLE_EQ(cdet.at(1, 1).getGain(), 2.0);

    // Row span — auto avoids PixelT which is a template param not in scope here
    auto row1 = cdet.Row(1);
    EXPECT_EQ(row1.size(), 4u);
    EXPECT_DOUBLE_EQ(row1[0].getGain(), 2.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  simulate() — shape and boundary checks
//  Note: simulate() requires DST::FITScube; these tests verify the interface
//  contracts (exceptions, dimensions) that do not depend on FITScube internals.
// ═══════════════════════════════════════════════════════════════════════════════

TEST(DetectorSimulator, WindowOutOfBounds)
{
    // simulate signature: (flux, t, c0, r0, ncols, nrows)
    DetectorSimulator<> det(32, 16);
    std::valarray<double> flux(0.0, 8 * 8);

    // row window exceeds detector height: r0=10, nrows=8 → 10+8=18 > 16
    EXPECT_THROW(det.simulate(flux, 1.0, 0, 10, 8, 8), std::out_of_range);

    // col window exceeds detector width: c0=26, ncols=8 → 26+8=34 > 32
    EXPECT_THROW(det.simulate(flux, 1.0, 26, 0, 8, 8), std::out_of_range);
}

TEST(DetectorSimulator, WindowFluxSizeMismatch)
{
    DetectorSimulator<> det(32, 16);
    std::valarray<double> flux(0.0, 63); // should be 8*8 = 64

    EXPECT_THROW(det.simulate(flux, 1.0, 0, 0, 8, 8), std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  _engine override — verify virtual dispatch reaches the derived hook
// ═══════════════════════════════════════════════════════════════════════════════

TEST(DetectorSimulator, EngineOverride)
{
    // Derived class that fills every output pixel with a sentinel value
    // instead of calling pixelresponse().
    class SentinelDetector : public DetectorSimulator<>
    {
    public:
        SentinelDetector(size_t w, size_t h) : DetectorSimulator<>(w, h) {}

        bool engineCalled = false;

    protected:
        void _engine(const std::valarray<double>&,
                     double,
                     size_t begin, size_t end,
                     std::valarray<double>& frame) const override
        {
            const_cast<SentinelDetector*>(this)->engineCalled = true;
            for (size_t i = begin; i < end; ++i)
                frame[i] = 42.0; // sentinel — never produced by pixelresponse()
        }
    };

    // Use a small detector so the sequential path is taken (_computeThreadCount → 1)
    // and we exercise the direct _engine() call in simulate().
    SentinelDetector det(4, 4);
    std::valarray<double> flux(0.0, 16);

    auto result = det.simulate(flux, 1.0);

    EXPECT_TRUE(det.engineCalled);
    // Verify every pixel was set to the sentinel
    const auto* data = result->GetData<double>();
    ASSERT_NE(data, nullptr);
    for (size_t i = 0; i < 16; ++i)
        EXPECT_DOUBLE_EQ((*data)[i], 42.0) << "pixel " << i;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  SimulateThreaded — smoke test on a large detector (exercises parallel path)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(DetectorSimulator, SimulateThreaded)
{
    // 4096×4096 = 16 777 216 pixels — well above the 100*nCores threshold
    // on any reasonable machine. This verifies the parallel path does not crash
    // and returns a result of the correct size.
    const size_t W = 4096, H = 4096;
    DetectorSimulator<> det(W, H);

    std::valarray<double> flux(10.0, W * H);
    auto result = det.simulate(flux, 1.0);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->Size(1), W);
    EXPECT_EQ(result->Size(2), H);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  DetectorConfig — initialize(PixelConfig)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(DetectorConfig, UniformDouble)
{
    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig{ {"gain", 3.0} });
    const auto& cdet = det;  // bind to const to force call to const accessors

    for (size_t i = 0; i < det.size(); ++i)
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(), 3.0);
}

TEST(DetectorConfig, UniformInt)
{
    // JSON integers land as int64_t — getAsDouble() must cast correctly
    DetectorSimulator<> det(4, 4);
    PixelConfig cfg;
    cfg.set<int64_t>("gain", 2);
    det.initialize(cfg);

    const auto& cdet = det;  // bind to const to force call to const accessors
    for (size_t i = 0; i < det.size(); ++i)
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(), 2.0);
}

TEST(DetectorConfig, UnknownKeyIgnored)
{
    DetectorSimulator<> det(4, 4);
    PixelConfig cfg{ {"unknown_field", 99.0} };
    EXPECT_NO_THROW(det.initialize(cfg));
}

TEST(DetectorConfig, MultipleKeys)
{
    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig{
        {"gain",      2.0},
        {"dark",      0.5},
        {"QE",        0.9},
        {"readNoise", 3.0}
    });

    const auto& cdet = det;  // bind to const to force call to const accessors
    
    for (size_t i = 0; i < det.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(),              2.0);
        EXPECT_DOUBLE_EQ(cdet.at(i).getDarkCurrent(),       0.5);
        EXPECT_DOUBLE_EQ(cdet.at(i).getQuantumEfficiency(), 0.9);
        EXPECT_DOUBLE_EQ(cdet.at(i).getReadNoise(),         3.0);
    }
}

TEST(DetectorConfig, FromFitsFile_form1)
{
    DSL::FITSimg<double> gain({4,4});
    for(size_t r=0; r<gain.Size(); ++r)
        gain.SetPixelValue(static_cast<double>(r+1),r);

    std::string url = "test_data";

    if(std::filesystem::exists(url+"/test_gain.fits"))
        std::filesystem::remove_all(url);
    
    std::filesystem::create_directory(url);
    

    gain.Write(url+"/test_gain.fits");

    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig{
        {"gain.file",url+"/test_gain.fits"},
        {"gain.hdu", 1},
        {"dark",      0.5},
        {"QE",        0.9},
        {"readNoise", 3.0}
    });

    const auto& cdet = det;  // bind to const to force call to const accessors
    
    for (size_t i = 0; i < det.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(),              static_cast<double>(i+1));
        EXPECT_DOUBLE_EQ(cdet.at(i).getDarkCurrent(),       0.5);
        EXPECT_DOUBLE_EQ(cdet.at(i).getQuantumEfficiency(), 0.9);
        EXPECT_DOUBLE_EQ(cdet.at(i).getReadNoise(),         3.0);
    }

    std::filesystem::remove_all(url+"/test_gain.fits");
}

TEST(DetectorConfig, FromFitsFile_form2)
{
    DSL::FITSimg<double> gain({4,4});
    for(size_t r=0; r<gain.Size(); ++r)
        gain.SetPixelValue(static_cast<double>(r+1),r);
    gain.SetName("GAIN");

    std::string url = "test_data";

    if(std::filesystem::exists(url+"/test_gain.fits"))
        std::filesystem::remove_all(url);
    
    std::filesystem::create_directory(url);
    

    gain.Write(url+"/test_gain.fits");

    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig{
        {"gain.file",url+"/test_gain.fits"},
        {"gain.hdu", std::string("GAIN")},
        {"dark",      0.5},
        {"QE",        0.9},
        {"readNoise", 3.0}
    });

    const auto& cdet = det;  // bind to const to force call to const accessors
    
    for (size_t i = 0; i < det.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(),              static_cast<double>(i+1));
        EXPECT_DOUBLE_EQ(cdet.at(i).getDarkCurrent(),       0.5);
        EXPECT_DOUBLE_EQ(cdet.at(i).getQuantumEfficiency(), 0.9);
        EXPECT_DOUBLE_EQ(cdet.at(i).getReadNoise(),         3.0);
    }

    std::filesystem::remove_all(url+"/test_gain.fits");
}

TEST(DetectorConfig, FromFitsFile_form3)
{
    DSL::FITSimg<double> gain({4,4});
    for(size_t r=0; r<gain.Size(); ++r)
        gain.SetPixelValue(static_cast<double>(r+1),r);

    std::string url = "test_data";

    if(std::filesystem::exists(url+"/test_gain.fits"))
        std::filesystem::remove_all(url);
    
    std::filesystem::create_directory(url);
    

    gain.Write(url+"/test_gain.fits");

    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig{
        {"gain","{\"file\":\""+url+"/test_gain.fits\",\"hdu\":" + std::to_string(1) + "}"},
        {"dark",      0.5},
        {"QE",        0.9},
        {"readNoise", 3.0}
    });

    const auto& cdet = det;  // bind to const to force call to const accessors
    
    for (size_t i = 0; i < det.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(),              static_cast<double>(i+1));
        EXPECT_DOUBLE_EQ(cdet.at(i).getDarkCurrent(),       0.5);
        EXPECT_DOUBLE_EQ(cdet.at(i).getQuantumEfficiency(), 0.9);
        EXPECT_DOUBLE_EQ(cdet.at(i).getReadNoise(),         3.0);
    }

    std::filesystem::remove_all(url+"/test_gain.fits");
}

TEST(DetectorConfig, FromFitsFile_form4)
{
    DSL::FITSimg<double> gain({4,4});
    for(size_t r=0; r<gain.Size(); ++r)
        gain.SetPixelValue(static_cast<double>(r+1),r);
    gain.SetName("GAIN");

    std::string url = "test_data";

    if(std::filesystem::exists(url+"/test_gain.fits"))
        std::filesystem::remove_all(url);
    
    std::filesystem::create_directory(url);
    
    gain.Write(url+"/test_gain.fits");

    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig{
        {"gain","{\"file\":\""+url+"/test_gain.fits\",\"hdu\":\"GAIN\"}"},
        {"dark",      0.5},
        {"QE",        0.9},
        {"readNoise", 3.0}
    });

    const auto& cdet = det;  // bind to const to force call to const accessors

    for (size_t i = 0; i < det.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(),              static_cast<double>(i+1));
        EXPECT_DOUBLE_EQ(cdet.at(i).getDarkCurrent(),       0.5);
        EXPECT_DOUBLE_EQ(cdet.at(i).getQuantumEfficiency(), 0.9);
        EXPECT_DOUBLE_EQ(cdet.at(i).getReadNoise(),         3.0);
    }

    std::filesystem::remove_all(url+"/test_gain.fits");
}

namespace
{
    // Two-HDU FITS file of 4x4 images: HDU 1 (primary) holds 10, 11, …; HDU 2, EXTNAME "SECOND", holds 1, 2, …
    std::string writeTwoHduFile(const std::string& path)
    {
        DSL::FITSimg<double> primary({4,4});
        DSL::FITSimg<double> second({4,4});
        for(size_t r=0; r<primary.Size(); ++r)
        {
            primary.SetPixelValue(static_cast<double>(r+10), r);
            second .SetPixelValue(static_cast<double>(r+1),  r);
        }
        second.SetName("SECOND");

        std::filesystem::create_directories(std::filesystem::path(path).parent_path());
        primary.Write(path, true);

        DSL::FITSmanager file(path, false);
        file.AppendImage(second);
        file.Close();
        return path;
    }
}

// The form a JSON configuration file gives: a JSON object {"file", "hdu"} (stored as a json
// value, not a JSON-encoded string), the HDU an index, an EXTNAME, or left out (primary HDU).
TEST(DetectorConfig, FromFitsFile_jsonObject)
{
    const std::string path = writeTwoHduFile("test_data/test_two_hdu.fits");

    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig(json{
        {"gain",      {{"file", path}, {"hdu", 2}}},
        {"baseline",  {{"file", path}, {"hdu", "SECOND"}}},
        {"dark",      {{"file", path}}},
        {"readNoise", 3.0}
    }));

    const auto& cdet = det;  // bind to const to force call to const accessors

    for (size_t i = 0; i < det.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(),        static_cast<double>(i+1));
        EXPECT_DOUBLE_EQ(cdet.at(i).getBaseline(),    static_cast<double>(i+1));
        EXPECT_DOUBLE_EQ(cdet.at(i).getDarkCurrent(), static_cast<double>(i+10));
        EXPECT_DOUBLE_EQ(cdet.at(i).getReadNoise(),   3.0);
    }

    std::filesystem::remove(path);
}

// loadImage() / loadMap() on a key initialize() does not know, as a tool reads its input flux.
TEST(DetectorConfig, LoadImageAndMap)
{
    const std::string path = writeTwoHduFile("test_data/test_two_hdu_load.fits");
    const DetectorSimulator<> det(4, 4);

    // Absent: no image, empty map.
    EXPECT_FALSE(det.loadImage(PixelConfig{}, "flux"));
    EXPECT_EQ(det.loadMap(PixelConfig{}, "flux").size(), 0u);

    // Number: no image, uniform map.
    const PixelConfig uniform{{"flux", 2.5}};
    EXPECT_FALSE(det.loadImage(uniform, "flux"));
    const std::valarray<double> flat = det.loadMap(uniform, "flux");
    ASSERT_EQ(flat.size(), det.size());
    EXPECT_DOUBLE_EQ(flat.min(), 2.5);
    EXPECT_DOUBLE_EQ(flat.max(), 2.5);

    // Every FITS form, with the HDU each one selects.
    const std::vector<std::pair<PixelConfig, double>> forms{
        {PixelConfig(json{{"flux", {{"file", path}, {"hdu", 2}}}}),        1.0},
        {PixelConfig(json{{"flux", {{"file", path}, {"hdu", "SECOND"}}}}), 1.0},
        {PixelConfig(json{{"flux", {{"file", path}}}}),                    10.0},
        {PixelConfig(json{{"flux", path}}),                                10.0},
        {PixelConfig(json{{"flux", path}, {"flux.hdu", 2}}),               1.0},
        {PixelConfig(json{{"flux.file", path}, {"flux.hdu", "SECOND"}}),   1.0},
    };
    for (const auto& [cfg, first] : forms)
    {
        const std::shared_ptr<DSL::FITScube> image = det.loadImage(cfg, "flux");
        ASSERT_TRUE(image) << cfg.toJson().dump();
        EXPECT_EQ(image->Size(1), 4u);
        EXPECT_EQ(image->Size(2), 4u);

        const std::valarray<double> map = det.loadMap(cfg, "flux");
        ASSERT_EQ(map.size(), det.size()) << cfg.toJson().dump();
        for (size_t i = 0; i < map.size(); ++i)
            EXPECT_DOUBLE_EQ(map[i], first + static_cast<double>(i)) << cfg.toJson().dump();
    }

    // Malformed entries.
    EXPECT_THROW(det.loadImage(PixelConfig(json{{"flux", {{"hdu", 2}}}}), "flux"),                    std::invalid_argument);
    EXPECT_THROW(det.loadImage(PixelConfig(json{{"flux", {{"file", path}, {"hdu", 0}}}}), "flux"),    std::invalid_argument);
    EXPECT_THROW(det.loadImage(PixelConfig(json{{"flux", {{"file", path}, {"hdu", 2.5}}}}), "flux"),  std::invalid_argument);
    EXPECT_THROW(det.loadImage(PixelConfig(json{{"flux", json::array({1, 2})}}), "flux"),             std::invalid_argument);

    std::filesystem::remove(path);
}

TEST(DetectorConfig, CustomConfigOverride)
{
    // Derived detector that overrides initialize() to handle an extra key
    // while delegating standard keys to the base class.
    class TrackingDetector : public DetectorSimulator<>
    {
    public:
        TrackingDetector(size_t w, size_t h) : DetectorSimulator<>(w, h) {}
        int    initCallCount = 0;
        double receivedCustomValue = -1.0;

        void initialize(const PixelConfig& cfg) override
        {
            DetectorSimulator<>::initialize(cfg);   // applies gain, dark, etc.
            ++initCallCount;
            receivedCustomValue = cfg.getAsDouble("custom", -1.0);
        }
    };

    TrackingDetector det(4, 4);
    det.initialize(PixelConfig{
        {"gain",   1.5},
        {"custom", 7.7}
    });

    EXPECT_EQ(det.initCallCount, 1);
    EXPECT_DOUBLE_EQ(det.receivedCustomValue, 7.7);

    const auto& cdet = det;
    for (size_t i = 0; i < det.size(); ++i)
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(), 1.5);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Geometry — physical detector, origin of the matrix, detector-sized maps
// ═══════════════════════════════════════════════════════════════════════════════

// Write a W × H double map whose pixel (x, y) holds value(x, y) to test_data/<name> and
// return its path. CFITSIO does not overwrite: a file left by a previous run is removed first.
template<typename F>
static std::string writeGeometryMap(const std::string& name, size_t W, size_t H, F value)
{
    const std::filesystem::path dir  = "test_data";
    const std::filesystem::path file = dir / name;

    std::filesystem::create_directories(dir);
    std::filesystem::remove(file);

    DSL::FITSimg<double> map({W, H});
    for (size_t y = 0; y < H; ++y)
        for (size_t x = 0; x < W; ++x)
            map.SetPixelValue(value(x, y), y * W + x);
    map.Write(file.string());
    return file.string();
}

TEST(DetectorGeometry, DefaultGeometryIsTheMatrix)
{
    DetectorSimulator<> det(8, 4);
    EXPECT_EQ(det.getGeometry().width,  8u);
    EXPECT_EQ(det.getGeometry().height, 4u);
    EXPECT_EQ(det.getOrigin(), std::make_pair(size_t(0), size_t(0)));
    EXPECT_FALSE(det.getGeometry().isReference(0, 0));

    det.initialize(PixelConfig{{"gain", 1.0}});   // no geometry key: unchanged
    EXPECT_EQ(det.getGeometry().width,  8u);
    EXPECT_EQ(det.getGeometry().height, 4u);
    EXPECT_EQ(det.getOrigin(), std::make_pair(size_t(0), size_t(0)));
}

TEST(DetectorGeometry, ExtractTakesTheMatrixSection)
{
    // 4 × 3 detector whose pixel (x, y) holds x + 4 y; 2 × 2 matrix at (1, 1).
    Geometry g(4, 3);
    g.origin = std::make_pair(size_t(1), size_t(1));

    std::valarray<double> detectorMap(12);
    for (size_t i = 0; i < detectorMap.size(); ++i)
        detectorMap[i] = static_cast<double>(i);

    const auto section = g.extract(detectorMap, 4, 3, 2, 2, "map");
    ASSERT_EQ(section.size(), 4u);
    EXPECT_DOUBLE_EQ(section[0], 5.0);    // (1, 1)
    EXPECT_DOUBLE_EQ(section[1], 6.0);    // (2, 1)
    EXPECT_DOUBLE_EQ(section[2], 9.0);    // (1, 2)
    EXPECT_DOUBLE_EQ(section[3], 10.0);   // (2, 2)

    const std::valarray<double> matrixMap(7.0, 4);
    const auto asIs = g.extract(matrixMap, 2, 2, 2, 2, "map");   // matrix-sized: used as is
    ASSERT_EQ(asIs.size(), 4u);
    EXPECT_DOUBLE_EQ(asIs[3], 7.0);

    EXPECT_THROW(g.extract(std::valarray<double>(0.0, 9), 3, 3, 2, 2, "map"), std::invalid_argument);   // neither size
    EXPECT_THROW(g.extract(detectorMap, 3, 4, 2, 2, "map"), std::invalid_argument);   // 12 values, wrong shape
    EXPECT_THROW(g.extract(detectorMap, 4, 2, 2, 2, "map"), std::invalid_argument);   // size and shape disagree
}

TEST(DetectorGeometry, DetectorSizedMapIsCroppedAtTheOrigin)
{
    // 8 × 8 gain map of an 8 × 8 detector: matrix pixel (x, y) takes map pixel (x0 + x, y0 + y).
    auto value = [](size_t x, size_t y) { return 1.0 + static_cast<double>(x + 8 * y); };
    const std::string path = writeGeometryMap("geometry_gain_8x8.fits", 8, 8, value);

    DetectorSimulator<> det(4, 4);
    const auto& cdet = det;   // const accessors

    det.initialize(PixelConfig{{"gain", path}, {"detectorSize", nlohmann::json::array({8, 8})}});
    for (size_t row = 0; row < 4; ++row)
        for (size_t col = 0; col < 4; ++col)
            ASSERT_DOUBLE_EQ(cdet.at(col, row).getGain(), value(col, row)) << "origin (0, 0), column " << col << ", row " << row;

    det.initialize(PixelConfig{{"gain", path}, {"detectorOrigin", nlohmann::json::array({2, 3})}});
    EXPECT_EQ(det.getOrigin(), std::make_pair(size_t(2), size_t(3)));
    for (size_t row = 0; row < 4; ++row)
        for (size_t col = 0; col < 4; ++col)
            ASSERT_DOUBLE_EQ(cdet.at(col, row).getGain(), value(2 + col, 3 + row)) << "origin (2, 3), column " << col << ", row " << row;

    std::filesystem::remove(path);
}

TEST(DetectorGeometry, MatrixSizedMapIsUsedAsIs)
{
    // A 4 × 4 map on a 4 × 4 matrix is in matrix coordinates, wherever the matrix is.
    auto value = [](size_t x, size_t y) { return 1.0 + static_cast<double>(x + 4 * y); };
    const std::string path = writeGeometryMap("geometry_gain_4x4.fits", 4, 4, value);

    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig{{"gain", path}, {"detectorSize", nlohmann::json::array({8, 8})},
                               {"detectorOrigin", nlohmann::json::array({4, 4})}});

    const auto& cdet = det;
    for (size_t row = 0; row < 4; ++row)
        for (size_t col = 0; col < 4; ++col)
            ASSERT_DOUBLE_EQ(cdet.at(col, row).getGain(), value(col, row)) << "column " << col << ", row " << row;

    std::filesystem::remove(path);
}

TEST(DetectorGeometry, MapOfAnotherSizeThrows)
{
    const std::string path = writeGeometryMap("geometry_gain_5x5.fits", 5, 5, [](size_t, size_t) { return 1.0; });

    DetectorSimulator<> det(4, 4);
    EXPECT_THROW(det.initialize(PixelConfig{{"gain", path}, {"detectorSize", nlohmann::json::array({8, 8})}}),
                 std::invalid_argument);

    std::filesystem::remove(path);
}

TEST(DetectorGeometry, MatrixThatDoesNotFitThrowsAndKeepsTheGeometry)
{
    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig{{"detectorSize",   nlohmann::json::array({8, 8})},
                               {"detectorOrigin", nlohmann::json::array({2, 2})}});

    EXPECT_THROW(det.initialize(PixelConfig{{"detectorOrigin", nlohmann::json::array({6, 0})}}), std::invalid_argument);
    EXPECT_EQ(det.getOrigin(), std::make_pair(size_t(2), size_t(2)));   // unchanged

    EXPECT_THROW(det.initialize(PixelConfig{{"detectorSize", nlohmann::json::array({3, 8})}}), std::invalid_argument);
    EXPECT_EQ(det.getGeometry().width, 8u);                              // unchanged

    EXPECT_THROW(det.initialize(PixelConfig{{"detectorOrigin", nlohmann::json::array({-1, 0})}}), std::invalid_argument);
    EXPECT_THROW(det.initialize(PixelConfig{{"detectorSize",   nlohmann::json::array({8})}}),     std::invalid_argument);
}

TEST(DetectorGeometry, DerivedGeometryIsConfiguredAndKept)
{
    // A detector family adds its own layout: here one reference column on the left of
    // an 8 × 8 detector, switched off by a "refColumn" key.
    class TestGeometry : public Geometry
    {
    public:
        bool refColumn = true;

        TestGeometry() : Geometry(8, 8) {}

        std::shared_ptr<Geometry> clone() const override { return std::make_shared<TestGeometry>(*this); }

        void configure(const PixelConfig& cfg) override
        {
            Geometry::configure(cfg);
            if (cfg.has("refColumn"))
                refColumn = cfg.getAsDouble("refColumn") != 0.0;
        }

        bool isReference(size_t px, size_t) const override { return refColumn && px == 0; }
    };

    DetectorSimulator<> det(4, 4, DOUBLE_IMG, std::make_shared<TestGeometry>());
    EXPECT_EQ(det.getGeometry().width, 8u);   // the geometry's own detector, not the matrix size
    EXPECT_TRUE(det.getGeometry().isReference(0, 3));

    det.initialize(PixelConfig{{"refColumn", 0.0}, {"detectorOrigin", nlohmann::json::array({4, 4})}});

    const auto* g = dynamic_cast<const TestGeometry*>(&det.getGeometry());
    ASSERT_NE(g, nullptr);                    // clone() kept the derived type
    EXPECT_FALSE(g->refColumn);               // the derived configure() was called...
    EXPECT_FALSE(det.getGeometry().isReference(0, 3));
    EXPECT_EQ(det.getOrigin(), std::make_pair(size_t(4), size_t(4)));   // ... and the base keys read
}

TEST(DetectorGeometry, CopiesKeepTheirOwnGeometry)
{
    DetectorSimulator<> det(4, 4);
    det.initialize(PixelConfig{{"detectorSize",   nlohmann::json::array({8, 8})},
                               {"detectorOrigin", nlohmann::json::array({2, 2})}});

    const DetectorSimulator<> copy = det;
    det.initialize(PixelConfig{{"detectorOrigin", nlohmann::json::array({4, 4})}});

    EXPECT_EQ(det.getOrigin(),  std::make_pair(size_t(4), size_t(4)));
    EXPECT_EQ(copy.getOrigin(), std::make_pair(size_t(2), size_t(2)));
}
