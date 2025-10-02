// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file spin_implicit_grid.cpp
 * \brief Unit tests for spin's ImplicitGrid class
 *
 * Uses gtest's TYPED_TESTS to test ImplicitGrid in 1,2 and 3 dimensions.
 */

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

#include <vector>
#include <unordered_set>
#include <algorithm>

/*!
 * Templated test fixture for ImplicitGrid tests
 *
 * \tparam T Is expected to wrap an integer, in its Value field
 */
template <typename T>
class ImplicitGridTest : public ::testing::Test
{
public:
  using GridT = T;
  using GridCell = typename GridT::GridCell;
  using SpacePt = typename GridT::SpacePoint;
  using BBox = typename GridT::SpatialBoundingBox;

  static const int DIM = SpacePt::DIMENSION;
};

/*! Type list for TypedTests on ImplicitGrid */
using MyTypes =
  ::testing::Types<axom::spin::ImplicitGrid<1>, axom::spin::ImplicitGrid<2>, axom::spin::ImplicitGrid<3>>;

TYPED_TEST_SUITE(ImplicitGridTest, MyTypes);

template <typename T>
class ImplicitGridExecTest;

template <int Dim, typename Exec>
class ImplicitGridExecTest<axom::spin::ImplicitGrid<Dim, Exec>> : public ::testing::Test
{
public:
  using GridT = axom::spin::ImplicitGrid<Dim, Exec>;
  using GridCell = typename GridT::GridCell;
  using SpacePt = typename GridT::SpacePoint;
  using BBox = typename GridT::SpatialBoundingBox;
  using ExecSpace = Exec;

  static const int DIM = SpacePt::DIMENSION;
};

/*! Type list for TypedTests on ImplicitGrid with execution space set*/
using MyTypesExec = ::testing::Types<
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  axom::spin::ImplicitGrid<1, axom::OMP_EXEC>,
  axom::spin::ImplicitGrid<2, axom::OMP_EXEC>,
  axom::spin::ImplicitGrid<3, axom::OMP_EXEC>,
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  axom::spin::ImplicitGrid<1, axom::CUDA_EXEC<256>>,
  axom::spin::ImplicitGrid<2, axom::CUDA_EXEC<256>>,
  axom::spin::ImplicitGrid<3, axom::CUDA_EXEC<256>>,
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  axom::spin::ImplicitGrid<1, axom::HIP_EXEC<256>>,
  axom::spin::ImplicitGrid<2, axom::HIP_EXEC<256>>,
  axom::spin::ImplicitGrid<3, axom::HIP_EXEC<256>>,
#endif
  axom::spin::ImplicitGrid<1, axom::SEQ_EXEC>,
  axom::spin::ImplicitGrid<2, axom::SEQ_EXEC>,
  axom::spin::ImplicitGrid<3, axom::SEQ_EXEC>>;

TYPED_TEST_SUITE(ImplicitGridExecTest, MyTypesExec);

TYPED_TEST(ImplicitGridExecTest, initialization_vectorized)
{
  const int DIM = TestFixture::DIM;
  using GridCell = typename TestFixture::GridCell;
  using BBox = typename TestFixture::BBox;
  using GridT = typename TestFixture::GridT;
  using SpacePt = typename TestFixture::SpacePt;
  using ExecSpace = typename TestFixture::ExecSpace;

  SLIC_INFO("Test ImplicitGrid constructor with " << axom::execution_space<ExecSpace>::name()
                                                  << " execution space and in " << DIM << "D");

  GridCell res(10);
  BBox bbox(SpacePt::zero(), SpacePt::ones());
  int numElts = 100;

  // Tests default constructor followed by initialize()
  GridT grid1;
  EXPECT_FALSE(grid1.isInitialized());

  grid1.initialize(bbox, &res, numElts);
  EXPECT_TRUE(grid1.isInitialized());
  EXPECT_EQ(grid1.gridResolution(), res);
  EXPECT_EQ(grid1.numIndexElements(), numElts);

  // Tests initializing constructor
  GridT grid2(bbox, &res, numElts);
  EXPECT_TRUE(grid2.isInitialized());
  EXPECT_EQ(grid2.gridResolution(), res);
  EXPECT_EQ(grid2.numIndexElements(), numElts);

  // Tests initializing from primitive types
  GridT grid3(bbox.getMin().data(), bbox.getMax().data(), res.data(), numElts);
  EXPECT_TRUE(grid3.isInitialized());
  EXPECT_EQ(grid3.gridResolution(), res);
  EXPECT_EQ(grid3.numIndexElements(), numElts);
}

TYPED_TEST(ImplicitGridTest, resolution)
{
  const int DIM = TestFixture::DIM;
  using GridCell = typename TestFixture::GridCell;
  using BBox = typename TestFixture::BBox;
  using GridT = typename TestFixture::GridT;
  using SpacePt = typename TestFixture::SpacePt;

  SLIC_INFO("Test ImplicitGrid resolution in " << DIM << "D");

  BBox bbox(SpacePt::zero(), SpacePt::ones());

  // Set the number of elements so that the DIM^th root is an integer
  const int dimRes = 8;
  int numElts = 1;
  for(int i = 0; i < DIM; ++i)
  {
    numElts *= dimRes;
  }

  // Tests explicitly set grid resolution
  {
    GridCell res = GridCell::make_point(12, 8, 16);
    GridT grid(bbox, &res, numElts);
    EXPECT_EQ(res, grid.gridResolution());
    EXPECT_EQ(numElts, grid.numIndexElements());
  }

  // Tests implicitly set grid resolution
  {
    GridT grid(bbox, nullptr, numElts);
    GridCell expRes(dimRes);
    EXPECT_EQ(expRes, grid.gridResolution());
    EXPECT_EQ(numElts, grid.numIndexElements());
  }

  // Test that grid resolution is at least one in each dim
  {
    // Default case -- using nullptr for resolution
    int zeroMeshElts = 0;
    GridT grid(bbox, nullptr, zeroMeshElts);
    GridCell expRes = GridCell::ones();
    EXPECT_EQ(expRes, grid.gridResolution());
    EXPECT_EQ(0, grid.numIndexElements());
  }
  {
    // Even when explicitly setting resolution to zero
    GridCell zeroRes = GridCell::zero();
    GridT grid(bbox, &zeroRes, numElts);
    GridCell expRes = GridCell::ones();
    EXPECT_EQ(expRes, grid.gridResolution());
    EXPECT_EQ(numElts, grid.numIndexElements());
  }
}

TYPED_TEST(ImplicitGridTest, insert_contains)
{
  const int DIM = TestFixture::DIM;
  using GridCell = typename TestFixture::GridCell;
  using BBox = typename TestFixture::BBox;
  using GridT = typename TestFixture::GridT;
  using SpacePt = typename TestFixture::SpacePt;
  using RangeBox = axom::primal::BoundingBox<int, DIM>;

  SLIC_INFO("Testing ImplicitGrid insert() and contains() in " << DIM << "D");

  // Note: A 10 x 10 x 10 implicit grid in the unit cube.
  //       Grid cells have a spacing of .1 along each dimension

  GridCell res(10);
  BBox bbox(SpacePt(0.), SpacePt(1.));
  const int maxElts = 10;
  const int numDefinedObjs = 7;

  GridT grid(bbox, &res, maxElts);

  /// Test insertion of several objects under assumption that we're
  /// in a unit cube with 10 grid cells per dimension resolution
  BBox objBox[maxElts];
  RangeBox ranges[maxElts];

  {
    // setup test data
    // Note: initializer list constructor ignores coordinates higher than its dimension
    SpacePt pts[] = {// Object 1 -- completely inside grid
                     SpacePt {.15, .25, .05},
                     SpacePt {.45, .25, .35},

                     // Object 2 -- extends past upper grid boundaries
                     SpacePt {.85, .85, .85},
                     SpacePt {1.25, 1.25, 1.25},

                     // Object 3 -- single point at center of a grid cell
                     SpacePt {.25, .35, .15},

                     // Object 4 -- single point on grid lines
                     SpacePt {.7, .7, .7},

                     // Object 5 -- completely outside grid
                     SpacePt {1.85, 1.85, 1.85},
                     SpacePt {2.25, 2.25, 2.25},

                     // Object 6 -- overlaps entire grid
                     SpacePt {-1.85, -1.85, -1.85},
                     SpacePt {2.25, 2.25, 2.25}};

    // Explicitly define expected ranges of cells
    // where we expect objects to be indexed
    GridCell rangePts[] = {GridCell {1, 2, 0},  // obj1
                           GridCell {4, 2, 3},

                           GridCell {8, 8, 8},  // obj2
                           GridCell {9, 9, 9},

                           GridCell {2, 3, 1},  // obj3

                           GridCell {6, 6, 6},  // obj4
                           GridCell {7, 7, 7}};

    objBox[0] = BBox();  // unused
    objBox[1] = BBox(pts[0], pts[1]);
    objBox[2] = BBox(pts[2], pts[3]);
    objBox[3] = BBox(pts[4]);
    objBox[4] = BBox(pts[5]);
    objBox[5] = BBox(pts[6], pts[7]);
    objBox[6] = BBox(pts[8], pts[9]);

    ranges[0] = RangeBox();  // unused
    ranges[1] = RangeBox(rangePts[0], rangePts[1]);
    ranges[2] = RangeBox(rangePts[2], rangePts[3]);
    ranges[3] = RangeBox(rangePts[4]);
    ranges[4] = RangeBox(rangePts[5], rangePts[6]);
    ranges[5] = RangeBox();  // unused
    ranges[6] = RangeBox();  // unused
  }

  // Insert objects into implicit grid.
  // Note: We do not insert an object with id 0
  grid.insert(numDefinedObjs - 1, objBox + 1, 1);

  // Test that points are contained in the expected grid cells
  const int i_max = DIM >= 1 ? grid.gridResolution()[0] : 1;
  const int j_max = DIM >= 2 ? grid.gridResolution()[1] : 1;
  const int k_max = DIM >= 3 ? grid.gridResolution()[2] : 1;
  for(int i = 0; i < i_max; ++i)
  {
    for(int j = 0; j < j_max; ++j)
    {
      for(int k = 0; k < k_max; ++k)
      {
        GridCell gridCell {i, j, k};

        // Item 0 -- was not inserted, so is never in grid
        {
          int idx = 0;
          EXPECT_FALSE(grid.contains(gridCell, idx));
          EXPECT_FALSE(grid.getCandidates(gridCell).test(idx));
        }

        // Objects 1,2,3 and 4 are expected to be in grid
        // when corresponding ranges contain the query point
        {
          for(int idx = 1; idx < 5; ++idx)
          {
            bool exp = ranges[idx].contains(gridCell);
            EXPECT_EQ(exp, grid.contains(gridCell, idx));
            EXPECT_EQ(exp, grid.getCandidates(gridCell).test(idx));
          }
        }

        // Object 5 is always outside grid
        {
          int idx = 5;
          EXPECT_FALSE(grid.contains(gridCell, idx));
          EXPECT_FALSE(grid.getCandidates(gridCell).test(idx));
        }

        // Object 6 is always inside grid
        {
          int idx = 6;
          EXPECT_TRUE(grid.contains(gridCell, idx));
          EXPECT_TRUE(grid.getCandidates(gridCell).test(idx));
        }

        // Test an index that is out of range
        {
          int idx = maxElts + 1;
          EXPECT_FALSE(grid.contains(gridCell, idx));
          EXPECT_FALSE(grid.getCandidates(gridCell).test(idx));
        }
      }
    }
  }
}

TYPED_TEST(ImplicitGridTest, get_candidates_pt)
{
  const int DIM = TestFixture::DIM;
  using GridCell = typename TestFixture::GridCell;
  using BBox = typename TestFixture::BBox;
  using GridT = typename TestFixture::GridT;
  using SpacePt = typename TestFixture::SpacePt;

  SLIC_INFO("Test ImplicitGrid getCandidates() for points in " << DIM << "D");

  using IndexType = typename GridT::IndexType;
  using CandidateBitset = typename GridT::BitsetType;
  using CandidateVector = std::vector<IndexType>;

  // Note: A 10 x 10 x 10 implicit grid in the unit cube.
  //       Grid cells have a spacing of .1 along each dimension
  GridCell res(10);
  BBox bbox(SpacePt(0.), SpacePt(1.));
  const int maxElts = 10;

  GridT grid(bbox, &res, maxElts);

  BBox objBox1(SpacePt {.15, .25, .05}, SpacePt {.45, .25, .35});
  grid.insert(objBox1, 1);

  {
    SpacePt queryPts[] = {
      // First three query points inside only obj1
      objBox1.getMin(),
      objBox1.getMax(),
      objBox1.getCentroid(),

      // Next two points are not in obj1, but in same grid cells
      SpacePt {0.11, 0.21, 0.01},
      SpacePt {0.49, 0.29, 0.39},

      // Next four points are not in obj1 or any of same cells
      SpacePt(0.55),
      SpacePt {.99, 0.25, 0.25},  // outside coord 0
      SpacePt {.35, 0.99, 0.25},  // outside coord 1
      SpacePt {.35, 0.25, 0.99},  // outside coord 2
    };

    // Test some points that are expected to match
    for(int i = 0; i < 5; ++i)
    {
      const std::size_t expSize = 1;
      const IndexType expIdx = 1;

      const SpacePt& queryPt = queryPts[i];

      // Test getCandidates() which returns a bitset
      CandidateBitset candidateBits = grid.getCandidates(queryPt);
      EXPECT_EQ(expSize, candidateBits.count());
      EXPECT_TRUE(candidateBits.test(expIdx));

      // Test getCandidatesAsArray() which returns a vector
      CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
      EXPECT_EQ(expSize, candidateVec.size());
      EXPECT_EQ(expIdx, candidateVec[0]);
    }

    // Test some points that are expected to not match
    for(int i = 5; i < 6 + DIM; ++i)
    {
      const std::size_t expSize = 0;

      const SpacePt& queryPt = queryPts[i];

      // Test getCandidates() which returns a bitset
      CandidateBitset candidateBits = grid.getCandidates(queryPt);
      EXPECT_EQ(expSize, candidateBits.count());

      // Test getCandidatesAsArray() which returns a vector
      CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
      EXPECT_EQ(expSize, candidateVec.size());
    }
  }

  BBox objBox2(SpacePt {.75, .85, .85}, SpacePt(.85));
  grid.insert(objBox2, 2);

  BBox objBox3(SpacePt {.85, .85, .75}, SpacePt(.95));
  grid.insert(objBox3, 3);

  // test a point that should contain only obj2
  {
    const std::size_t expSize = 1;
    SpacePt queryPt = SpacePt {.75, .85, .85};

    CandidateBitset candidateBits = grid.getCandidates(queryPt);
    EXPECT_EQ(expSize, candidateBits.count());
    EXPECT_FALSE(candidateBits.test(1));
    EXPECT_TRUE(candidateBits.test(2));
    EXPECT_FALSE(candidateBits.test(3));

    CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
    EXPECT_EQ(expSize, candidateVec.size());
    EXPECT_EQ(IndexType(2), candidateVec[0]);
  }

  // test a point that should contain only obj3
  {
    const std::size_t expSize = 1;
    SpacePt queryPt = SpacePt(.91);

    CandidateBitset candidateBits = grid.getCandidates(queryPt);
    EXPECT_EQ(expSize, candidateBits.count());
    EXPECT_FALSE(candidateBits.test(1));
    EXPECT_FALSE(candidateBits.test(2));
    EXPECT_TRUE(candidateBits.test(3));

    CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
    EXPECT_EQ(expSize, candidateVec.size());
    EXPECT_EQ(IndexType(3), candidateVec[0]);
  }

  // test a point that should contain obj2 and obj3, but not obj1
  {
    const std::size_t expSize = 2;
    SpacePt queryPt = SpacePt::make_point(.85, .85, .85);

    CandidateBitset candidateBits = grid.getCandidates(queryPt);
    EXPECT_EQ(expSize, candidateBits.count());
    EXPECT_FALSE(candidateBits.test(1));
    EXPECT_TRUE(candidateBits.test(2));
    EXPECT_TRUE(candidateBits.test(3));

    CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
    EXPECT_EQ(expSize, candidateVec.size());

    auto beg = candidateVec.cbegin();
    auto end = candidateVec.cend();

    bool has1 = end != std::find(beg, end, IndexType(1));
    EXPECT_FALSE(has1);

    bool has2 = end != std::find(beg, end, IndexType(2));
    EXPECT_TRUE(has2);

    bool has3 = end != std::find(beg, end, IndexType(3));
    EXPECT_TRUE(has3);
  }
}

TYPED_TEST(ImplicitGridExecTest, get_candidates_pt_vectorized)
{
  const int DIM = TestFixture::DIM;
  using GridCell = typename TestFixture::GridCell;
  using BBox = typename TestFixture::BBox;
  using GridT = typename TestFixture::GridT;
  using SpacePt = typename TestFixture::SpacePt;
  using ExecSpace = typename TestFixture::ExecSpace;

  SLIC_INFO("Test ImplicitGrid getCandidatesAsArray() with "
            << axom::execution_space<ExecSpace>::name() << " execution space for points in " << DIM
            << "D");

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Note: A 10 x 10 x 10 implicit grid in the unit cube.
  //       Grid cells have a spacing of .1 along each dimension
  GridCell res(10);
  BBox bbox(SpacePt(0.), SpacePt(1.));
  const int maxElts = 10;

  GridT grid(bbox, &res, maxElts);

  axom::Array<BBox> objBox1(1, 1, hostAllocID);
  objBox1[0] = BBox {SpacePt {.15, .25, .05}, SpacePt {.45, .25, .35}};

  axom::Array<BBox> objBox1_device = axom::Array<BBox>(objBox1, kernelAllocID);
  grid.insert(objBox1_device[0], 1);

  {
    axom::Array<SpacePt> queryPts(9, 9, hostAllocID);

    // First three query points inside only obj1
    queryPts[0] = objBox1[0].getMin();
    queryPts[1] = objBox1[0].getMax();
    queryPts[2] = objBox1[0].getCentroid();

    // Next two points are not in obj1, but in same grid cells
    queryPts[3] = SpacePt {0.11, 0.21, 0.01};
    queryPts[4] = SpacePt {0.49, 0.29, 0.39};

    // Next four points are not in obj1 or any of same cells
    queryPts[5] = SpacePt(0.55);
    queryPts[6] = SpacePt {.99, 0.25, 0.25};  // outside coord 0
    queryPts[7] = SpacePt {.35, 0.99, 0.25};  // outside coord 1
    queryPts[8] = SpacePt {.35, 0.25, 0.99};  // outside coord 2

    axom::Array<int> count, offset, candidates;
    {
      axom::Array<SpacePt> queryPtsDevice = axom::Array<SpacePt>(queryPts, kernelAllocID);
      axom::Array<int> countDevice(9, 9, kernelAllocID);
      axom::Array<int> offsetDevice(9, 9, kernelAllocID);
      axom::Array<int> candidatesDevice;

      // Run query against implicit grid
      grid.getCandidatesAsArray(9, queryPtsDevice.data(), offsetDevice, countDevice, candidatesDevice);

      // Copy results back to the host
      count = axom::Array<int>(countDevice, hostAllocID);
      offset = axom::Array<int>(offsetDevice, hostAllocID);
      candidates = axom::Array<int>(candidatesDevice, hostAllocID);
    }

    // Test some points that are expected to match
    for(int i = 0; i < 5; ++i)
    {
      std::unordered_set<int> expected {1};
      EXPECT_EQ(1, count[i]);
      EXPECT_EQ(1, candidates[i]);
    }

    // Test some points that are expected to not match
    for(int i = 5; i < 6 + DIM; ++i)
    {
      EXPECT_EQ(0, count[i]);
    }
  }

  axom::Array<BBox> objBox2(1, 1, hostAllocID);
  objBox2[0] = BBox {SpacePt {.75, .85, .85}, SpacePt(.85)};
  axom::Array<BBox> objBox2_device = axom::Array<BBox>(objBox2, kernelAllocID);
  grid.insert(objBox2_device[0], 2);

  axom::Array<BBox> objBox3(1, 1, hostAllocID);
  objBox3[0] = BBox {SpacePt {.85, .85, .75}, SpacePt(.95)};
  axom::Array<BBox> objBox3_device = axom::Array<BBox>(objBox3, kernelAllocID);
  grid.insert(objBox3_device[0], 3);

  {
    axom::Array<SpacePt> queryPts(3, 3, hostAllocID);

    // Should only be inside obj2
    queryPts[0] = SpacePt {.75, .85, .85};
    // Should only be inside obj3
    queryPts[1] = SpacePt(.91);
    // Should be inside obj2 and obj3, but not obj1
    queryPts[2] = SpacePt {.85, .85, .85};

    axom::Array<int> count, offset, candidates;
    {
      axom::Array<SpacePt> queryPtsDevice = axom::Array<SpacePt>(queryPts, kernelAllocID);
      axom::Array<int> countDevice(3, 3, kernelAllocID);
      axom::Array<int> offsetDevice(3, 3, kernelAllocID);
      axom::Array<int> candidatesDevice;
      // Run query against implicit grid
      grid.getCandidatesAsArray(3, queryPtsDevice.data(), offsetDevice, countDevice, candidatesDevice);

      // Copy results back to the host
      count = axom::Array<int>(countDevice, hostAllocID);
      offset = axom::Array<int>(offsetDevice, hostAllocID);
      candidates = axom::Array<int>(candidatesDevice, hostAllocID);
    }
    std::unordered_set<int> expected[] = {{2}, {3}, {2, 3}};

    for(int i = 0; i < 3; i++)
    {
      std::unordered_set<int> actual(candidates.data() + offset[i],
                                     candidates.data() + offset[i] + count[i]);
      EXPECT_EQ(expected[i], actual);
    }
  }
}

TYPED_TEST(ImplicitGridTest, get_candidates_box)
{
  const int DIM = TestFixture::DIM;
  using GridCell = typename TestFixture::GridCell;
  using BBox = typename TestFixture::BBox;
  using GridT = typename TestFixture::GridT;
  using SpacePt = typename TestFixture::SpacePt;

  SLIC_INFO("Test ImplicitGrid getCandidates() for boxes in " << DIM << "D");

  using IndexType = typename GridT::IndexType;
  using CandidateBitset = typename GridT::BitsetType;
  using CandidateVector = std::vector<IndexType>;

  // Note: A 10 x 10 x 10 implicit grid in the unit cube.
  //       Grid cells have a spacing of .1 along each dimension
  GridCell res(10);
  BBox bbox(SpacePt(0.), SpacePt(1.));
  const int maxElts = 30;

  GridT grid(bbox, &res, maxElts);
  const int i_max = DIM >= 1 ? grid.gridResolution()[0] : 1;
  const int j_max = DIM >= 2 ? grid.gridResolution()[1] : 1;
  const int k_max = DIM >= 3 ? grid.gridResolution()[2] : 1;

  // Add some boxes to the spatial index
  {
    // Assumes unit cube
    // Create 10 objects per dimension (up to D == 3)
    // Each object will be in one box in that dimension and all
    // boxes in the other dimensions
    //
    // Each grid cell will be covered by DIM objects

    for(int i = 0; i < res[0]; ++i)
    {
      grid.insert(BBox(SpacePt {0.025 + 0.1 * i, .05, .05}, SpacePt {0.075 + 0.1 * i, .95, .95}), i);
    }

    if(DIM >= 2)
    {
      for(int i = 0; i < res[1]; ++i)
      {
        grid.insert(BBox(SpacePt {0.05, 0.025 + 0.1 * i, .05}, SpacePt {0.95, 0.075 + 0.1 * i, .95}),
                    10 + i);
      }
    }

    if(DIM >= 3)
    {
      for(int i = 0; i < res[2]; ++i)
      {
        grid.insert(
          BBox(SpacePt {0.05, 0.05, 0.025 + 0.1 * i}, SpacePt {0.95, 0.95, 0.075 + 0.1 * i}),
          20 + i);
      }
    }

    // Check that each grid cell contains DIM objects
    for(int i = 0; i < i_max; ++i)
    {
      for(int j = 0; j < j_max; ++j)
      {
        for(int k = 0; k < k_max; ++k)
        {
          SpacePt queryPt {i * .1 + .05, j * .1 + .05, k * .1 + .05};
          EXPECT_EQ(DIM, grid.getCandidates(queryPt).count());
        }
      }
    }
  }

  //// Run some queries

  // Empty box -- covers no objects
  {
    BBox query;

    CandidateBitset bits = grid.getCandidates(query);
    EXPECT_EQ(0, bits.count());

    CandidateVector vec = grid.getCandidatesAsArray(query);
    EXPECT_EQ(0, vec.size());
  }

  // Box covers entire domain -- covers all objects
  {
    BBox query = bbox;

    CandidateBitset bits = grid.getCandidates(query);
    EXPECT_EQ(DIM * 10, bits.count());

    CandidateVector vec = grid.getCandidatesAsArray(query);
    EXPECT_EQ(DIM * 10, vec.size());
  }

  // Box is larger than domain -- covers all objects
  {
    BBox query(SpacePt(-1.), SpacePt(2.));

    CandidateBitset bits = grid.getCandidates(query);
    EXPECT_EQ(DIM * 10, bits.count());

    CandidateVector vec = grid.getCandidatesAsArray(query);
    EXPECT_EQ(DIM * 10, vec.size());
  }

  // Box only covers first quadrant/octant of domain
  {
    BBox query(SpacePt(-1.), SpacePt(0.45));

    int expCount = DIM * 5;  // covers five cells in each dimension

    CandidateBitset bits = grid.getCandidates(query);
    EXPECT_EQ(expCount, bits.count());

    CandidateVector vec = grid.getCandidatesAsArray(query);
    EXPECT_EQ(expCount, vec.size());
  }

  // Box covers a single cell
  {
    BBox query(SpacePt(0.525), SpacePt(0.575));

    int expCount = DIM;  // covers one cell per dimension

    CandidateBitset bits = grid.getCandidates(query);
    EXPECT_EQ(expCount, bits.count());

    CandidateVector vec = grid.getCandidatesAsArray(query);
    EXPECT_EQ(expCount, vec.size());

    EXPECT_EQ(DIM >= 1, bits.test(5));
    EXPECT_EQ(DIM >= 2, bits.test(15));
    EXPECT_EQ(DIM >= 3, bits.test(25));
  }

  // Box only covers last quadrant/octant of domain
  {
    BBox query(SpacePt(0.55), SpacePt(2.));

    int expCount = DIM * 5;  // covers five cells in each dimension

    CandidateBitset bits = grid.getCandidates(query);
    EXPECT_EQ(expCount, bits.count());

    CandidateVector vec = grid.getCandidatesAsArray(query);
    EXPECT_EQ(expCount, vec.size());
  }

  // Box covers middle of domain
  {
    BBox query(SpacePt(0.25), SpacePt(0.75));

    int expCount = DIM * 6;  // covers six cells in each dimension

    CandidateBitset bits = grid.getCandidates(query);
    EXPECT_EQ(expCount, bits.count());

    CandidateVector vec = grid.getCandidatesAsArray(query);
    EXPECT_EQ(expCount, vec.size());
  }

  // Box is inverted -- BoundingBox constructor fixes this
  {
    BBox query(SpacePt(2.), SpacePt(-1));

    int expCount = DIM * 10;  // covers all cells in each dimension

    CandidateBitset bits = grid.getCandidates(query);
    EXPECT_EQ(expCount, bits.count());

    CandidateVector vec = grid.getCandidatesAsArray(query);
    EXPECT_EQ(expCount, vec.size());
  }
}

TYPED_TEST(ImplicitGridExecTest, get_candidates_box_vectorized)
{
  const int DIM = TestFixture::DIM;
  using GridCell = typename TestFixture::GridCell;
  using BBox = typename TestFixture::BBox;
  using GridT = typename TestFixture::GridT;
  using SpacePt = typename TestFixture::SpacePt;
  using ExecSpace = typename TestFixture::ExecSpace;

  SLIC_INFO("Test ImplicitGrid getCandidatesAsArray() with "
            << axom::execution_space<ExecSpace>::name() << " execution space for boxes in " << DIM
            << "D");

  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
  int unifiedAllocID = axom::execution_space<ExecSpace>::allocatorID();

#if defined(AXOM_USE_GPU) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
    unifiedAllocID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  }
#endif

  // Note: A 10 x 10 x 10 implicit grid in the unit cube.
  //       Grid cells have a spacing of .1 along each dimension
  GridCell res(10);
  BBox bbox(SpacePt(0.), SpacePt(1.));
  const int maxElts = 30;

  // Initialize Implicit Grid in unified memory.
  // getCandidates() function not yet __host __device__ decorated
  GridT grid(bbox, &res, maxElts, unifiedAllocID);
  const int i_max = DIM >= 1 ? grid.gridResolution()[0] : 1;
  const int j_max = DIM >= 2 ? grid.gridResolution()[1] : 1;
  const int k_max = DIM >= 3 ? grid.gridResolution()[2] : 1;

  axom::Array<BBox> BBoxes(maxElts, maxElts, hostAllocID);
  axom::Array<SpacePt> queryPt(1, 1, hostAllocID);

  // Add some boxes to the spatial index
  {
    // Assumes unit cube
    // Create 10 objects per dimension (up to D == 3)
    // Each object will be in one box in that dimension and all
    // boxes in the other dimensions
    //
    // Each grid cell will be covered by DIM objects

    for(int i = 0; i < res[0]; ++i)
    {
      BBoxes[i] = BBox {SpacePt {0.025 + 0.1 * i, .05, .05}, SpacePt {0.075 + 0.1 * i, .95, .95}};
    }

    if(DIM >= 2)
    {
      for(int i = 0; i < res[1]; ++i)
      {
        BBoxes[10 + i] =
          BBox {SpacePt {0.05, 0.025 + 0.1 * i, .05}, SpacePt {0.95, 0.075 + 0.1 * i, .95}};
      }
    }

    if(DIM >= 3)
    {
      for(int i = 0; i < res[2]; ++i)
      {
        BBoxes[20 + i] =
          BBox {SpacePt {0.05, 0.05, 0.025 + 0.1 * i}, SpacePt {0.95, 0.95, 0.075 + 0.1 * i}};
      }
    }

    // Copy boxes to unified memory (for getCandidates() call)
    axom::Array<BBox> BBoxesDevice = axom::Array<BBox>(BBoxes, unifiedAllocID);

    grid.insert(10 * DIM, BBoxesDevice.data(), 0);
    // Check that each grid cell contains DIM objects
    for(int i = 0; i < i_max; ++i)
    {
      for(int j = 0; j < j_max; ++j)
      {
        for(int k = 0; k < k_max; ++k)
        {
          queryPt[0] = SpacePt {i * .1 + .05, j * .1 + .05, k * .1 + .05};
          axom::Array<SpacePt> queryPtDevice = axom::Array<SpacePt>(queryPt, unifiedAllocID);
          EXPECT_EQ(DIM, grid.getCandidates(queryPtDevice[0]).count());
        }
      }
    }
  }

  //// Run some queries
  constexpr int N_QUERIES = 8;
  axom::Array<BBox> queryBoxes(N_QUERIES, N_QUERIES, hostAllocID);

  // Empty box -- covers no objects
  queryBoxes[0] = BBox {};
  // Box covers entire domain -- covers all objects
  queryBoxes[1] = bbox;
  // Box is larger than domain -- covers all objects
  queryBoxes[2] = BBox {SpacePt(-1.), SpacePt(2.)};
  // Box only covers first quadrant/octant of domain
  queryBoxes[3] = BBox {SpacePt(-1.), SpacePt(0.45)};
  // Box covers a single cell
  queryBoxes[4] = BBox {SpacePt(0.525), SpacePt(0.575)};
  // Box only covers last quadrant/octant of domain
  queryBoxes[5] = BBox {SpacePt(0.55), SpacePt(2.)};
  // Box covers middle of domain
  queryBoxes[6] = BBox {SpacePt(0.25), SpacePt(0.75)};
  // Box is inverted -- BoundingBox constructor fixes this
  queryBoxes[7] = BBox {SpacePt(2.), SpacePt(-1)};

  std::vector<int> expectedVisits {0, DIM * 10, DIM * 10, DIM * 5, DIM, DIM * 5, DIM * 6, DIM * 10};

  axom::Array<int> offset, count, candidates;
  {
    axom::Array<BBox> queryBoxesDevice = axom::Array<BBox>(queryBoxes, kernelAllocID);
    axom::Array<int> countDevice(N_QUERIES, N_QUERIES, kernelAllocID);
    axom::Array<int> offsetDevice(N_QUERIES, N_QUERIES, kernelAllocID);
    axom::Array<int> candidatesDevice;

    // Run query against implicit grid
    grid.getCandidatesAsArray(N_QUERIES,
                              queryBoxesDevice.data(),
                              offsetDevice,
                              countDevice,
                              candidatesDevice);

    // Copy results back to the host
    count = axom::Array<int>(countDevice, hostAllocID);
    offset = axom::Array<int>(offsetDevice, hostAllocID);
    candidates = axom::Array<int>(candidatesDevice, hostAllocID);
  }

  for(int i = 0; i < N_QUERIES; i++)
  {
    EXPECT_EQ(expectedVisits[i], count[i]);
  }

  // check single cell results

  std::unordered_set<int> query_4(candidates.data() + offset[4],
                                  candidates.data() + offset[4] + count[4]);
  EXPECT_EQ(DIM >= 1, query_4.count(5) == 1);
  EXPECT_EQ(DIM >= 2, query_4.count(15) == 1);
  EXPECT_EQ(DIM >= 3, query_4.count(25) == 1);
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  result = RUN_ALL_TESTS();

  return result;
}
