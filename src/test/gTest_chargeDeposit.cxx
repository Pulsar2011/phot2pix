//
// Created by gillard on 2026-11-26.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// gTest_chargeDeposit.cxx
//
// Unit tests of phot2pix::ChargeDeposit, the base class of the charge sources (CosmicDeposit, ...),
// through a minimal source defined here:
// (a) seed and reproducible random streams (makeRng());
// (b) copies: copy constructor, assignment, clone();
// (c) streams drawn from several threads;
// (d) push(): deposits outside the detector or without charge are dropped.
//
// gTest_cosmicDeposit.cxx tests the cosmic-ray source itself.
//

#include <gtest/gtest.h>
#include <phot2pix/DetectorEffect.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <thread>
#include <vector>

namespace
{
    using phot2pix::ChargeDeposit;
    using phot2pix::Deposit;
    using phot2pix::Geometry;
    using phot2pix::PixelConfig;

    // Minimal source: generate() puts a charge drawn from its stream in pixel (0, 0).
    // makeRng() and push() are exposed for the tests.
    class TestDeposit : public ChargeDeposit
    {
    public:
        std::shared_ptr<ChargeDeposit> clone() const override { return std::make_shared<TestDeposit>(*this); }

        void configure(const PixelConfig&) override {}

        void generate(double, const Geometry& g, std::vector<Deposit>& out) const override
        {
            std::mt19937_64 rng = makeRng();
            push(out, g, 0, 0, 1.0 + static_cast<double>(rng() % 1000));
        }

        // First values of the next stream (one call of makeRng()).
        std::vector<uint64_t> nextStream(size_t n = 4) const
        {
            std::mt19937_64 rng = makeRng();
            std::vector<uint64_t> v(n);
            for (uint64_t& x : v)
                x = rng();
            return v;
        }

        using ChargeDeposit::push;
    };
}

// ── (a) Seed and streams ─────────────────────────────────────────────────────────

TEST(ChargeDepositTest, SeedIsZeroByDefault)
{
    TestDeposit d;
    EXPECT_EQ(d.seed(), 0u);

    d.setSeed(42);
    EXPECT_EQ(d.seed(), 42u);
}

TEST(ChargeDepositTest, SameSeedGivesTheSameStreams)
{
    TestDeposit a, b;
    a.setSeed(1234);
    b.setSeed(1234);
    for (int n = 0; n < 5; ++n)
        EXPECT_EQ(a.nextStream(), b.nextStream()) << "call " << n;
}

TEST(ChargeDepositTest, SuccessiveStreamsDiffer)
{
    TestDeposit d;
    d.setSeed(1234);
    const std::vector<uint64_t> s0 = d.nextStream();
    const std::vector<uint64_t> s1 = d.nextStream();
    EXPECT_NE(s0, s1);
}

TEST(ChargeDepositTest, DifferentSeedsGiveDifferentStreams)
{
    TestDeposit a, b;
    a.setSeed(1);
    b.setSeed(2);
    EXPECT_NE(a.nextStream(), b.nextStream());
}

TEST(ChargeDepositTest, SetSeedRestartsTheSequence)
{
    TestDeposit d;
    d.setSeed(7);
    const std::vector<uint64_t> first = d.nextStream();
    d.nextStream();
    d.nextStream();

    d.setSeed(7);
    EXPECT_EQ(d.nextStream(), first);
}

TEST(ChargeDepositTest, SeedZeroIsNotReproducible)
{
    // std::random_device: equal streams have a negligible probability
    TestDeposit a, b;
    EXPECT_NE(a.nextStream(), b.nextStream());
    EXPECT_NE(a.nextStream(), a.nextStream());
}

TEST(ChargeDepositTest, GenerateDrawsTheNextStream)
{
    TestDeposit a, b;
    a.setSeed(55);
    b.setSeed(55);
    const Geometry g(4, 4);
    std::vector<Deposit> outA, outB;
    a.generate(1.0, g, outA);
    b.generate(1.0, g, outB);
    ASSERT_EQ(outA.size(), 1u);
    ASSERT_EQ(outB.size(), 1u);
    EXPECT_EQ(outA[0].idx, 0u);
    EXPECT_DOUBLE_EQ(outA[0].electrons, outB[0].electrons);
}

// ── (b) Copies ───────────────────────────────────────────────────────────────────

TEST(ChargeDepositTest, CopyContinuesTheSequence)
{
    TestDeposit d;
    d.setSeed(99);
    d.nextStream();   // call 0
    d.nextStream();   // call 1

    TestDeposit copy(d);
    EXPECT_EQ(copy.seed(), 99u);
    EXPECT_EQ(copy.nextStream(), d.nextStream());   // both draw call 2
}

TEST(ChargeDepositTest, AssignmentContinuesTheSequence)
{
    TestDeposit d;
    d.setSeed(99);
    d.nextStream();

    TestDeposit other;
    other.setSeed(5);
    other = d;
    EXPECT_EQ(other.seed(), 99u);
    EXPECT_EQ(other.nextStream(), d.nextStream());
}

TEST(ChargeDepositTest, CloneKeepsTheTypeTheSeedAndTheSequence)
{
    TestDeposit d;
    d.setSeed(3);
    d.nextStream();

    const std::shared_ptr<ChargeDeposit> c = d.clone();
    const auto* t = dynamic_cast<const TestDeposit*>(c.get());
    ASSERT_NE(t, nullptr);
    EXPECT_EQ(t->seed(), 3u);
    EXPECT_EQ(t->nextStream(), d.nextStream());
}

// ── (c) Threads ──────────────────────────────────────────────────────────────────

TEST(ChargeDepositTest, ThreadsGetDistinctStreamsOfTheSequence)
{
    constexpr size_t n = 64;

    TestDeposit d;
    d.setSeed(2026);
    std::vector<std::vector<uint64_t>> streams(n);
    std::vector<std::thread> threads;
    threads.reserve(n);
    for (size_t i = 0; i < n; ++i)
        threads.emplace_back([&d, &streams, i] { streams[i] = d.nextStream(); });
    for (std::thread& t : threads)
        t.join();

    std::sort(streams.begin(), streams.end());
    EXPECT_EQ(std::adjacent_find(streams.begin(), streams.end()), streams.end()) << "two threads drew the same stream";

    // Whatever the order of the threads, they drew the calls 0 .. n-1
    TestDeposit s;
    s.setSeed(2026);
    std::vector<std::vector<uint64_t>> sequential(n);
    for (std::vector<uint64_t>& v : sequential)
        v = s.nextStream();
    std::sort(sequential.begin(), sequential.end());
    EXPECT_EQ(streams, sequential);
}

// ── (d) push() ───────────────────────────────────────────────────────────────────

TEST(ChargeDepositTest, PushDropsDepositsOutsideTheDetectorOrWithoutCharge)
{
    const Geometry g(4, 3);
    std::vector<Deposit> out;

    TestDeposit::push(out, g, 1, 2, 5.0);                          // kept: idx = 2 * 4 + 1
    TestDeposit::push(out, g, 4, 0, 1.0);                          // X = width
    TestDeposit::push(out, g, 0, 3, 1.0);                          // Y = height
    TestDeposit::push(out, g, static_cast<size_t>(-1), 0, 1.0);    // a negative coordinate wraps around
    TestDeposit::push(out, g, 0, 0, 0.0);                          // no charge
    TestDeposit::push(out, g, 0, 0, -1.0);                         // negative charge

    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].idx, 9u);
    EXPECT_DOUBLE_EQ(out[0].electrons, 5.0);
}
