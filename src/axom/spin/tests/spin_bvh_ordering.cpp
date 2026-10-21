// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/spin/BVH.hpp"
#include "axom/spin/SpaceFillingCurve.hpp"
#include "axom/slic.hpp"
#include "gtest/gtest.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace
{
using axom::IndexType;
using axom::spin::SpaceFillingCurve;

// A test-only policy tag lets us exercise the public BVHPolicy extension point.
constexpr auto legacyPolicy = static_cast<axom::spin::BVHType>(-1);

// Policy with only the four-argument buildImpl interface.
class LegacyLinearBVH : public axom::spin::policy::LinearBVH<double, 2, axom::SEQ_EXEC>
{
public:
  template <typename BoxIndexable>
  void buildImpl(const BoxIndexable boxes, IndexType count, double scale, int allocator)
  {
    axom::spin::policy::LinearBVH<double, 2, axom::SEQ_EXEC>::buildImpl(boxes, count, scale, allocator);
  }
};
}  // namespace

namespace axom::spin
{
template <>
struct BVHPolicy<double, 2, SEQ_EXEC, legacyPolicy>
{
  using ImplType = LegacyLinearBVH;
};
}  // namespace axom::spin

namespace
{
// Check each query's expected input IDs, allowing any traversal order.
void expectCandidates(const char* queryKind,
                      const axom::Array<IndexType>& offsets,
                      const axom::Array<IndexType>& counts,
                      const axom::Array<IndexType>& candidates,
                      const std::vector<std::vector<IndexType>>& expected)
{
  SCOPED_TRACE(queryKind);
  const int host = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const axom::Array<IndexType> hostOffsets(offsets, host);
  const axom::Array<IndexType> hostCounts(counts, host);
  const axom::Array<IndexType> hostCandidates(candidates, host);
  ASSERT_EQ(hostOffsets.size(), expected.size());
  ASSERT_EQ(hostCounts.size(), expected.size());

  IndexType offset = 0;
  for(std::size_t q = 0; q < expected.size(); ++q)
  {
    SCOPED_TRACE(::testing::Message() << "query " << q);
    ASSERT_EQ(hostOffsets[q], offset);
    ASSERT_EQ(hostCounts[q], expected[q].size());
    ASSERT_LE(offset + hostCounts[q], hostCandidates.size());
    std::vector<IndexType> actual;
    for(IndexType i = 0; i < hostCounts[q]; ++i)
    {
      actual.push_back(hostCandidates[offset + i]);
    }
    std::sort(actual.begin(), actual.end());
    EXPECT_EQ(actual, expected[q]);
    offset += hostCounts[q];
  }
  EXPECT_EQ(offset, hostCandidates.size());
}

// Read original input IDs by leaf slot, independent of traversal visit order.
template <typename BVH>
std::vector<IndexType> leafOrder(const BVH& bvh, IndexType count)
{
  using ExecSpace = typename BVH::ExecSpaceType;
  using Box = typename BVH::BoxType;
  axom::Array<IndexType> ids(count, count, bvh.getAllocatorID());
  ids.fill(-1);
  const auto idsView = ids.view();
  const auto traverser = bvh.getTraverser();
  const auto bounds = bvh.getBounds();
  axom::for_all<ExecSpace>(1, [=] AXOM_HOST_DEVICE(IndexType) {
    auto leaf = [=](std::int32_t index, const std::int32_t* leafIDs) {
      idsView[index] = leafIDs[index];
    };
    traverser.traverse_tree(bounds, leaf, [](const Box&, const Box&) { return true; });
  });
  const axom::Array<IndexType> hostIds(ids, axom::execution_space<axom::SEQ_EXEC>::allocatorID());
  return std::vector<IndexType>(hostIds.begin(), hostIds.end());
}

// Run each test in 2D and 3D, with float/double and each enabled execution space.
template <typename BVH>
class BVHOrderingTest : public ::testing::Test
{ };

using BVHTypes = ::testing::Types<
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  axom::spin::BVH<2, axom::OMP_EXEC, float>,
  axom::spin::BVH<2, axom::OMP_EXEC, double>,
  axom::spin::BVH<3, axom::OMP_EXEC, float>,
  axom::spin::BVH<3, axom::OMP_EXEC, double>,
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA)
  axom::spin::BVH<2, axom::CUDA_EXEC<256>, float>,
  axom::spin::BVH<2, axom::CUDA_EXEC<256>, double>,
  axom::spin::BVH<3, axom::CUDA_EXEC<256>, float>,
  axom::spin::BVH<3, axom::CUDA_EXEC<256>, double>,
  #endif
  #if defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP)
  axom::spin::BVH<2, axom::HIP_EXEC<256>, float>,
  axom::spin::BVH<2, axom::HIP_EXEC<256>, double>,
  axom::spin::BVH<3, axom::HIP_EXEC<256>, float>,
  axom::spin::BVH<3, axom::HIP_EXEC<256>, double>,
  #endif
#endif
  axom::spin::BVH<2, axom::SEQ_EXEC, float>,
  axom::spin::BVH<2, axom::SEQ_EXEC, double>,
  axom::spin::BVH<3, axom::SEQ_EXEC, float>,
  axom::spin::BVH<3, axom::SEQ_EXEC, double>>;

TYPED_TEST_SUITE(BVHOrderingTest, BVHTypes);

// Both curves must return the original box IDs for hits, boundary hits, and misses.
AXOM_TYPED_TEST(BVHOrderingTest, query_candidates)
{
  using Box = typename TypeParam::BoxType;
  using Point = typename TypeParam::PointType;
  using Ray = typename TypeParam::RayType;
  using Vector = typename Ray::VectorType;

  TypeParam bvh;
  bvh.setScaleFactor(1);
  const int allocator = bvh.getAllocatorID();

  // Boxes 0 and 1 overlap; box 2 is nested in 0 with the same centroid; box 3 is isolated.
  // Point(value) sets every coordinate, so the same fixture works in 2D and 3D.
  const axom::Array<Box> boxes({Box(Point(0.), Point(2)),
                                Box(Point(1), Point(3)),
                                Box(Point(.5), Point(1.5)),
                                Box(Point(4), Point(5))},
                               allocator);
  const axom::Array<Point> points({Point(1), Point(3), Point(4.5), Point(6)}, allocator);
  const axom::Array<Box> queryBoxes({Box(Point(.75), Point(1.25)),
                                     Box(Point(2.5), Point(3.5)),
                                     Box(Point(4.25), Point(4.75)),
                                     Box(Point(6), Point(7))},
                                    allocator);
  const axom::Array<Ray> rays({Ray(Point(-1), Vector(1)),
                               Ray(Point(3.5), Vector(1)),
                               Ray(Point(6), Vector(1)),
                               Ray(Point(6), Vector(-1))},
                              allocator);
  axom::Array<IndexType> offsets(4, 4, allocator), counts(4, 4, allocator), candidates;

  for(auto ordering : {SpaceFillingCurve::Morton, SpaceFillingCurve::Hilbert})
  {
    SCOPED_TRACE(ordering == SpaceFillingCurve::Morton ? "Morton" : "Hilbert");
    bvh.setOrdering(ordering);
    ASSERT_EQ(bvh.initialize(boxes.view(), boxes.size()), axom::spin::BVH_BUILD_OK);

    bvh.findPoints(offsets, counts, candidates, points.size(), points.view());
    expectCandidates("points", offsets, counts, candidates, {{0, 1, 2}, {1}, {3}, {}});

    bvh.findBoundingBoxes(offsets, counts, candidates, queryBoxes.size(), queryBoxes.view());
    expectCandidates("boxes", offsets, counts, candidates, {{0, 1, 2}, {1}, {3}, {}});

    // Diagonal rays hit all boxes, only the isolated box, none, and all in reverse.
    bvh.findRays(offsets, counts, candidates, rays.size(), rays.view());
    expectCandidates("rays", offsets, counts, candidates, {{0, 1, 2, 3}, {3}, {}, {0, 1, 2, 3}});
  }
}

// Zero-width bounds and coincident boxes must still report every matching input ID.
AXOM_TYPED_TEST(BVHOrderingTest, degenerate_bounds)
{
  using Box = typename TypeParam::BoxType;
  using Point = typename TypeParam::PointType;
  TypeParam bvh;
  const int allocator = bvh.getAllocatorID();
  // Fill all but the last coordinate to make a segment in 2D or a square in 3D.
  const Box flat(Point(-1, Point::DIMENSION - 1), Point(1, Point::DIMENSION - 1));
  const Box point(Point(0.), Point(0.));
  const axom::Array<Point> queries({Point(0.), Point(2)}, allocator);
  axom::Array<IndexType> offsets(2, 2, allocator), counts(2, 2, allocator), candidates;
  for(const auto& box : {flat, point})
  {
    SCOPED_TRACE(::testing::Message() << "box " << box);
    const axom::Array<Box> boxes({box, box, box}, allocator);
    for(auto ordering : {SpaceFillingCurve::Morton, SpaceFillingCurve::Hilbert})
    {
      SCOPED_TRACE(ordering == SpaceFillingCurve::Morton ? "Morton" : "Hilbert");
      bvh.setOrdering(ordering);
      ASSERT_EQ(bvh.initialize(boxes.view(), boxes.size()), axom::spin::BVH_BUILD_OK);
      bvh.findPoints(offsets, counts, candidates, queries.size(), queries.view());
      expectCandidates("points", offsets, counts, candidates, {{0, 1, 2}, {}});
    }
  }
}

// Shrinking to one or zero boxes must remove old candidates; rebuilding must restore them.
AXOM_TYPED_TEST(BVHOrderingTest, reinitialize)
{
  using Box = typename TypeParam::BoxType;
  using Point = typename TypeParam::PointType;

  TypeParam bvh;
  const int allocator = bvh.getAllocatorID();
  const axom::Array<Box> boxes({Box(Point(0.), Point(1)), Box(Point(2), Point(3))}, allocator);
  const axom::Array<Point> queries({Point(.5), Point(2.5)}, allocator);
  axom::Array<IndexType> offsets(2, 2, allocator), counts(2, 2, allocator), candidates;
  const std::vector<std::vector<IndexType>> expectedByBoxCount[] = {{{}, {}}, {{0}, {}}, {{0}, {1}}};
  for(auto ordering : {SpaceFillingCurve::Morton, SpaceFillingCurve::Hilbert})
  {
    SCOPED_TRACE(ordering == SpaceFillingCurve::Morton ? "Morton" : "Hilbert");
    bvh.setOrdering(ordering);
    for(IndexType count : {2, 1, 0, 2})
    {
      SCOPED_TRACE(::testing::Message() << "box count " << count);
      ASSERT_EQ(bvh.initialize(boxes.view(), count), axom::spin::BVH_BUILD_OK);
      EXPECT_TRUE(bvh.isInitialized());
      bvh.findPoints(offsets, counts, candidates, queries.size(), queries.view());
      expectCandidates("points", offsets, counts, candidates, expectedByBoxCount[count]);
    }
  }
}

// Known cases distinguish the curves and check that equal keys retain input order.
AXOM_TYPED_TEST(BVHOrderingTest, leaf_order_and_stable_ties)
{
  using Box = typename TypeParam::BoxType;
  using Point = typename TypeParam::PointType;
  using Float = typename Point::CoordType;
  constexpr int DIM = Point::DIMENSION;
  constexpr IndexType corners = 1 << DIM;

  const int host = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  axom::Array<Box> boxes(2 * corners, 2 * corners, host);
  for(IndexType i = 0; i < corners; ++i)
  {
    // Reverse the binary corner numbering so Morton must sort too.
    Point p;
    for(int d = 0; d < DIM; ++d)
    {
      p[d] = Float(((corners - 1 - i) >> d) & 1);
    }
    boxes[i] = boxes[i + corners] = Box(p, p);
  }

  // Each corner occurs in both halves of the input; ties keep the earlier ID first.
  const std::vector<IndexType> morton2D {3, 7, 2, 6, 1, 5, 0, 4};
  const std::vector<IndexType> morton3D {7, 15, 6, 14, 5, 13, 4, 12, 3, 11, 2, 10, 1, 9, 0, 8};
  const std::vector<IndexType> hilbert2D {3, 7, 1, 5, 0, 4, 2, 6};
  const std::vector<IndexType> hilbert3D {7, 15, 3, 11, 1, 9, 5, 13, 4, 12, 0, 8, 2, 10, 6, 14};
  const auto& morton = DIM == 2 ? morton2D : morton3D;
  const auto& hilbert = DIM == 2 ? hilbert2D : hilbert3D;

  TypeParam bvh;
  const axom::Array<Box> deviceBoxes(boxes, bvh.getAllocatorID());
  ASSERT_EQ(bvh.initialize(deviceBoxes.data(), boxes.size()), axom::spin::BVH_BUILD_OK);
  EXPECT_EQ(bvh.getOrdering(), SpaceFillingCurve::Morton);
  EXPECT_EQ(leafOrder(bvh, boxes.size()), morton);

  bvh.setOrdering(SpaceFillingCurve::Hilbert);
  EXPECT_EQ(bvh.getOrdering(), SpaceFillingCurve::Hilbert);
  EXPECT_EQ(leafOrder(bvh, boxes.size()), morton);  // The setting takes effect on the next build.
  ASSERT_EQ(bvh.initialize(deviceBoxes.view(), boxes.size()), axom::spin::BVH_BUILD_OK);
  EXPECT_EQ(leafOrder(bvh, boxes.size()), hilbert);

  bvh.setOrdering(SpaceFillingCurve::Morton);
  ASSERT_EQ(bvh.initialize(deviceBoxes.view(), boxes.size()), axom::spin::BVH_BUILD_OK);
  EXPECT_EQ(leafOrder(bvh, boxes.size()), morton);
}

// Check clamping, truncation, and the fixed quantization depth before encoding each curve.
AXOM_TYPED_TEST(BVHOrderingTest, quantization)
{
  namespace lbvh = axom::spin::internal::linear_bvh;
  using ExecSpace = typename TypeParam::ExecSpaceType;
  using Point = typename TypeParam::PointType;
  using Float = typename Point::CoordType;
  using Vector = typename TypeParam::RayType::VectorType;
  using Cell = axom::primal::Point<std::int32_t, Point::DIMENSION>;
  constexpr int bits = Point::DIMENSION == 2 ? 16 : 10;
  constexpr int side = 1 << bits;
  struct Sample
  {
    Vector point;
    Cell cell;
  };

  // 2D uses the first two coordinates; 3D also checks the third.
  const Sample samples[] = {
    {Vector {-.25, 0, 1.25}, Cell {0, 0, side - 1}},         // Clamp below zero.
    {Vector {1, 1.25, -.25}, Cell {side - 1, side - 1, 0}},  // Clamp at and above one.
    {Vector {Float(2.5) / side, Float(3) / side, Float(3.5) / side}, Cell {2, 3, 3}},
    {Vector {Float(3.5) / side, .5, 0}, Cell {3, side / 2, 0}}};  // Truncate, don't round.
  const int allocator = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<std::uint32_t> keys(2, 2, allocator);
  const auto keyView = keys.view();
  for(const auto& sample : samples)
  {
    SCOPED_TRACE(::testing::Message() << "point " << sample.point);
    const auto point = sample.point;
    axom::for_all<ExecSpace>(1, [=] AXOM_HOST_DEVICE(IndexType) {
      keyView[0] = lbvh::encode_spatial_key(point, SpaceFillingCurve::Morton);
      keyView[1] = lbvh::encode_spatial_key(point, SpaceFillingCurve::Hilbert);
    });
    const axom::Array<std::uint32_t> hostKeys(keys,
                                              axom::execution_space<axom::SEQ_EXEC>::allocatorID());
    EXPECT_EQ(hostKeys[0], axom::spin::convertPointToMorton<std::uint32_t>(sample.cell));
    EXPECT_EQ(hostKeys[1], axom::spin::convertPointToHilbert<std::uint32_t>(sample.cell, bits));
  }
}
}  // namespace

// A four-argument policy accepts Morton and rejects Hilbert without replacing its tree.
TEST(spin_bvh_ordering, legacy_policy)
{
  using BVH = axom::spin::BVH<2, axom::SEQ_EXEC, double, legacyPolicy>;
  using Box = BVH::BoxType;
  using Point = BVH::PointType;
  const axom::Array<Box> boxes {Box(Point(0.), Point(1)), Box(Point(2), Point(3))};
  const axom::Array<Point> queries {Point(.5), Point(2.5)};
  axom::Array<IndexType> offsets(2), counts(2), candidates;
  BVH bvh;
  EXPECT_EQ(bvh.getOrdering(), SpaceFillingCurve::Morton);
  ASSERT_EQ(bvh.initialize(boxes.view(), boxes.size()), axom::spin::BVH_BUILD_OK);
  const auto bounds = bvh.getBounds();

  bvh.setOrdering(SpaceFillingCurve::Hilbert);
  for(IndexType count : {0, 1, 2})
  {
    EXPECT_EQ(bvh.initialize(boxes.view(), count), axom::spin::BVH_BUILD_FAILED);
    ASSERT_TRUE(bvh.isInitialized());
    EXPECT_EQ(bvh.getBounds(), bounds);
    bvh.findPoints(offsets, counts, candidates, queries.size(), queries.view());
    expectCandidates("preserved tree", offsets, counts, candidates, {{0}, {1}});
  }

  bvh.setOrdering(SpaceFillingCurve::Morton);
  const std::vector<std::vector<IndexType>> expectedByBoxCount[] = {{{}, {}}, {{0}, {}}, {{0}, {1}}};
  for(IndexType count : {0, 1, 2})
  {
    ASSERT_EQ(bvh.initialize(boxes.data(), count), axom::spin::BVH_BUILD_OK);
    bvh.findPoints(offsets, counts, candidates, queries.size(), queries.view());
    expectCandidates("legacy Morton", offsets, counts, candidates, expectedByBoxCount[count]);
  }

  BVH uninitialized;
  uninitialized.setOrdering(SpaceFillingCurve::Hilbert);
  EXPECT_EQ(uninitialized.initialize(boxes.data(), boxes.size()), axom::spin::BVH_BUILD_FAILED);
  EXPECT_FALSE(uninitialized.isInitialized());
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
