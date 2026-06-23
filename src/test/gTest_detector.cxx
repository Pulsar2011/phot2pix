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
    EXPECT_NO_THROW(det.initialize(cfg)); // _applyCustomConfig default: no-op
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

TEST(DetectorConfig, CustomConfigHook)
{
    // Derived detector that counts how many times _applyCustomConfig is called
    // and records the config it received.
    class TrackingDetector : public DetectorSimulator<>
    {
    public:
        TrackingDetector(size_t w, size_t h) : DetectorSimulator<>(w, h) {}
        int    hookCallCount = 0;
        double receivedCustomValue = -1.0;

    protected:
        void _applyCustomConfig(const PixelConfig& cfg) override
        {
            ++hookCallCount;
            receivedCustomValue = cfg.getAsDouble("custom", -1.0);
        }
    };

    TrackingDetector det(4, 4);
    det.initialize(PixelConfig{
        {"gain",   1.5},
        {"custom", 7.7}
    });

    EXPECT_EQ(det.hookCallCount, 1);
    EXPECT_DOUBLE_EQ(det.receivedCustomValue, 7.7);

    const auto& cdet = det;  // bind to const to force call to const accessors

    for (size_t i = 0; i < det.size(); ++i)
        EXPECT_DOUBLE_EQ(cdet.at(i).getGain(), 1.5);
}
