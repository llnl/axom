// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <limits>

#include "gtest/gtest.h"

#include "axom/core/execution/runtime_policy.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/spin/UniformGrid.hpp"

//-----------------------------------------------------------------------------
TEST(spin_uniform_grid, array_constructor)
{
  const int DIM = 3;

  double p_max[DIM] = {10, 10, 10};
  double p_min[DIM] = {0, 0, 0};
  const int resolution = 4;
  int res[DIM] = {resolution, resolution, resolution};

  axom::spin::UniformGrid<int, DIM> valid(p_min, p_max, res);
  EXPECT_EQ(valid.getNumBins(), resolution * resolution * resolution);
  EXPECT_TRUE(valid.isBinEmpty(0));
}

TEST(spin_uniform_grid, bbox_constructor)
{
  const int DIM = 3;
  using CoordType = double;
  using QPoint = axom::primal::Point<CoordType, DIM>;

  QPoint pmax = QPoint::make_point(4, 2, 3);
  QPoint pmin = QPoint::make_point(6, 8, 6);
  int res[DIM] = {2, 3, 4};

  axom::primal::BoundingBox<double, DIM> theBbox(pmin, pmax);
  axom::spin::UniformGrid<int, DIM> valid(theBbox, res);
  EXPECT_EQ(valid.getNumBins(), res[0] * res[1] * res[2]);
  EXPECT_TRUE(valid.isBinEmpty(0));
}

TEST(spin_uniform_grid, indexing)
{
  const int DIM = 3;
  using CoordType = double;
  using QPoint = axom::primal::Point<CoordType, DIM>;

  double p_max[DIM] = {100, 100, 100};
  double p_min[DIM] = {0, 0, 0};
  const int resolution = 100;
  int res[DIM] = {resolution, resolution, resolution};

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  using execSpace = axom::HIP_EXEC<256>;
#elif defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  using execSpace = axom::CUDA_EXEC<256>;
#else
  using execSpace = axom::SEQ_EXEC;
#endif

  axom::spin::UniformGrid<int, DIM, execSpace> valid(p_min, p_max, res);

  // valid has 100 bins in each dimension, and each bin has a
  // width of 1.0.  The bins are laid out in row-major order.
  // So, each increment in the x-dimension increases the bin index by 1,
  // an increment in the y-dimension increases the bin by 100, and an
  // increment in the z-dimension increases the bin by 100*100 = 10000.

  QPoint pt1 = QPoint::make_point(1.5, 0, 0);
  int expectedBin = 1;  // The 0th z-slab, the 0th y-row, the 1st x-bin
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt1));

  QPoint pt2 = QPoint::make_point(0, 1.5, 0);
  expectedBin = 100;  // The 0th z-slab, the 1st y-row, the 0th x-bin
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt2));

  QPoint pt3 = QPoint::make_point(99.5, 0, 99.5);
  // The (0-based) 99th z-slab, the 0th y-row, the 99th x-bin:
  expectedBin = 99 * 100 * 100 + 0 + 99;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt3));

  QPoint pt4 = QPoint::make_point(16.1, 99.6, 89.2);
  // The 89th z-slab, the 99th y-row, the 16th x-bin:
  expectedBin = 89 * 100 * 100 + 99 * 100 + 16;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt4));

  // This is a set of upper-boundary cases.
  // In general, a bin is "half-open", including its lower-boundary
  // planes but excluding its upper boundary planes.  The bins located
  // on the max-x, max-y, and max-z planes are exceptions: they include
  // the upper boundary, so that points on all faces of the UniformGrid
  // are indexed in the UniformGrid.

  // The max point, in the top-most bucket.
  QPoint pt5 = QPoint::make_point(100., 100., 100.);
  expectedBin = 99 * 100 * 100 + 99 * 100 + 99;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt5));

  // A point on the max-x plane
  QPoint pt6 = QPoint::make_point(100., 55.5, 62.2);
  expectedBin = 62 * 100 * 100 + 55 * 100 + 99;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt6));
  // A point on the max-y plane
  QPoint pt7 = QPoint::make_point(81.4, 100., 2.1);
  expectedBin = 2 * 100 * 100 + 99 * 100 + 81;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt7));
  // A point on the max-z plane
  QPoint pt8 = QPoint::make_point(45.6, 1.8, 100.);
  expectedBin = 99 * 100 * 100 + 1 * 100 + 45;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt8));
  // A point on the max-xy edge
  QPoint pt9 = QPoint::make_point(100., 100., 62.2);
  expectedBin = 62 * 100 * 100 + 99 * 100 + 99;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt9));
  // A point on the max-yz edge
  QPoint pt10 = QPoint::make_point(81.4, 100., 100.);
  expectedBin = 99 * 100 * 100 + 99 * 100 + 81;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt10));
  // A point on the max-xz edge
  QPoint pt11 = QPoint::make_point(100., 1.8, 100.);
  expectedBin = 99 * 100 * 100 + 1 * 100 + 99;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt11));

  // Now go outside the grid, and get some invalid ones.

  QPoint pt12 = QPoint::make_point(12.5, 100.1, 0);
  // Above 100 is over the fence.
  expectedBin = axom::spin::UniformGrid<int, DIM>::INVALID_BIN_INDEX;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt12));

  QPoint pt13 = QPoint::make_point(-0.5, 12, 54.3);
  // Below 0 is over (under?) the fence.
  expectedBin = axom::spin::UniformGrid<int, DIM>::INVALID_BIN_INDEX;
  EXPECT_EQ(expectedBin, valid.getBinIndex(pt13));
}

// For insertions and deletions, we use these next few helper functions
// to track the count of objects in each bin.

// Verify the count in each bin against the map maintained "by hand".
template <typename T, int NDIMS>
void checkBinCounts(axom::spin::UniformGrid<T, NDIMS>& v, std::map<int, int>& bincounts)
{
  int bcount = v.getNumBins();
  for(int i = 0; i < bcount; ++i)
  {
    bool binAgrees = (bincounts.count(i) < 1 && v.isBinEmpty(i)) ||
      (bincounts[i] == ((int)v.getBinContents(i).size()));
    EXPECT_TRUE(binAgrees) << "Difference at bin " << i << ": v has " << v.getBinContents(i).size()
                           << " and bincounts has "
                           << (((int)bincounts.count(i)) < 1 ? 0 : bincounts[i]);
  }
}

// Increment the count for a bin
void incr(std::map<int, int>& m, int idx)
{
  int dat = 0;
  if(m.count(idx) > 0)
  {
    dat = m[idx];
  }
  m[idx] = dat + 1;
}

// Zero out the count for a bin
void zero(std::map<int, int>& m, int idx) { m.erase(idx); }

TEST(spin_uniform_grid, add_stuff_3D)
{
  const int DIM = 3;
  using CoordType = double;
  using QPoint = axom::primal::Point<CoordType, DIM>;
  using QBBox = axom::primal::BoundingBox<CoordType, DIM>;

  double origin[DIM] = {0, 0, 0};
  const int mpt = 6;
  double maxpoint[DIM] = {mpt, mpt, mpt};
  const int resolution = 6;
  int res[DIM] = {resolution, resolution, resolution};
  axom::spin::UniformGrid<int, DIM> valid(origin, maxpoint, res);

  std::map<int, int> check;

  QPoint pt1 = QPoint::make_point(2.5, 2.5, 2.5);
  QBBox bbox1(pt1);
  valid.insert(bbox1, 1);
  int index = valid.getBinIndex(pt1);
  incr(check, index);
  {
    SCOPED_TRACE("One bin");
    checkBinCounts(valid, check);
  }

  // Add a point with a bounding box overlapping the first one
  QPoint pt2 = QPoint::make_point(2.1, 2.1, 2.1);
  QPoint pt3 = QPoint::make_point(4.2, 2.9, 2.1);
  QBBox bbox2(pt2, pt3);
  valid.insert(bbox2, 2);
  incr(check, valid.getBinIndex(pt2));
  incr(check, valid.getBinIndex(bbox2.getCentroid()));
  incr(check, valid.getBinIndex(pt3));
  {
    SCOPED_TRACE("Added an overlapping bin");
    checkBinCounts(valid, check);
  }

  // Try inserting something partially outside the index
  QPoint pt4 = QPoint::make_point(2.7, 5.2, -3.1);
  QPoint pt5 = QPoint::make_point(3.1, 9.0, 2.8);
  QBBox bbox3(pt4, pt5);
  valid.insert(bbox3, 4);
  for(int k = 0; k < 3; ++k)
  {
    for(int i = 2; i < 4; ++i)
    {
      incr(check, valid.getBinIndex(QPoint::make_point(i + 0.5, 5.5, k + 0.5)));
    }
  }
  {
    SCOPED_TRACE("Bbox partially overlaps grid");
    checkBinCounts(valid, check);
  }

  // Try inserting something completely outside the index
  QPoint pt6 = QPoint::make_point(-2, 7.2, -3.1);
  QPoint pt7 = QPoint::make_point(-0.1, 9.0, 2.8);
  QBBox bbox4(pt6, pt7);
  valid.insert(bbox4, 6);
  {
    SCOPED_TRACE("Bbox completely misses grid");
    checkBinCounts(valid, check);
  }

  // Try inserting something completely overlapping the box
  QPoint pt8 = QPoint::make_point(-2, -0.2, -3.1);
  QPoint pt9 = QPoint::make_point(7, 6.2, 6.1);
  QBBox bbox5(pt8, pt9);
  valid.insert(bbox5, 8);
  for(int i = 0; i < valid.getNumBins(); ++i)
  {
    incr(check, i);
  }
  {
    SCOPED_TRACE("Bbox completely overlaps grid");
    checkBinCounts(valid, check);
  }
}

TEST(spin_uniform_grid, delete_stuff_3D)
{
  const int DIM = 3;
  using CoordType = double;
  using QPoint = axom::primal::Point<CoordType, DIM>;
  using QBBox = axom::primal::BoundingBox<CoordType, DIM>;

  double origin[DIM] = {0, 0, 0};
  const int mpt = 6;
  double maxpoint[DIM] = {mpt, mpt, mpt};
  const int resolution = 6;
  int res[DIM] = {resolution, resolution, resolution};
  axom::spin::UniformGrid<int, DIM> valid(origin, maxpoint, res);

  std::map<int, int> check;

  // Insert something small (one bin), then clear its bin
  QPoint pt1 = QPoint::make_point(1.1, 1.1, 1.1);
  QBBox bbox1(pt1);
  valid.insert(bbox1, 1);
  int index = valid.getBinIndex(pt1);
  incr(check, index);
  valid.clear(index);
  zero(check, index);
  {
    SCOPED_TRACE("Inserted one, then deleted it");
    checkBinCounts(valid, check);
  }

  // Insert some objects, clear some bins
  valid.insert(bbox1, 1);
  incr(check, valid.getBinIndex(pt1));
  QPoint pt2 = QPoint::make_point(0.1, 0.1, -1.1);
  QPoint pt3 = QPoint::make_point(1.9, 2.3, 7.1);
  QBBox bbox2(pt2, pt3);
  valid.insert(bbox2, 2);
  for(int k = 0; k < 6; ++k)
  {
    for(int j = 0; j < 3; ++j)
    {
      for(int i = 0; i < 2; ++i)
      {
        incr(check, valid.getBinIndex(QPoint::make_point(i + 0.5, j + 0.5, k + 0.5)));
      }
    }
  }
  {
    SCOPED_TRACE("Insert two, clear several bins");
    QPoint test = QPoint::make_point(0.5, 2.5, 5.5);
    int tidx = valid.getBinIndex(test);
    valid.clear(tidx);
    zero(check, tidx);
    test = QPoint::make_point(0.5, 1.5, 0.5);
    tidx = valid.getBinIndex(test);
    valid.clear(tidx);
    zero(check, tidx);
    test = QPoint::make_point(0.5, 2.5, 5.5);
    tidx = valid.getBinIndex(test);
    valid.clear(tidx);
    zero(check, tidx);

    checkBinCounts(valid, check);
  }
  {
    SCOPED_TRACE("Try clearing the invalid bin");
    valid.clear(axom::spin::UniformGrid<int, DIM>::INVALID_BIN_INDEX);
    checkBinCounts(valid, check);
  }
}

TEST(spin_uniform_grid, add_stuff_2D)
{
  const int DIM = 2;
  using CoordType = double;
  using QPoint = axom::primal::Point<CoordType, DIM>;
  using QBBox = axom::primal::BoundingBox<CoordType, DIM>;

  double origin[DIM] = {0, 0};
  const int mpt = 6;
  double maxpoint[DIM] = {mpt, mpt};
  const int resolution = 6;
  int res[DIM] = {resolution, resolution};
  axom::spin::UniformGrid<int, DIM> valid(origin, maxpoint, res);

  std::map<int, int> check;

  QPoint pt1 = QPoint::make_point(2.5, 2.5);
  QBBox bbox1(pt1);
  valid.insert(bbox1, 1);
  int index = valid.getBinIndex(pt1);
  incr(check, index);
  {
    SCOPED_TRACE("One bin");
    checkBinCounts(valid, check);
  }

  // Add a point with a bounding box overlapping the first one
  QPoint pt2 = QPoint::make_point(2.1, 2.1);
  QPoint pt3 = QPoint::make_point(4.2, 2.9);
  QBBox bbox2(pt2, pt3);
  valid.insert(bbox2, 2);
  incr(check, valid.getBinIndex(pt2));
  incr(check, valid.getBinIndex(bbox2.getCentroid()));
  incr(check, valid.getBinIndex(pt3));
  {
    SCOPED_TRACE("Added an overlapping bin");
    checkBinCounts(valid, check);
  }

  // Try inserting something partially outside the index
  QPoint pt4 = QPoint::make_point(2.7, 5.2);
  QPoint pt5 = QPoint::make_point(3.1, 9.0);
  QBBox bbox3(pt4, pt5);
  valid.insert(bbox3, 4);
  for(int i = 2; i < 4; ++i)
  {
    incr(check, valid.getBinIndex(QPoint::make_point(i + 0.5, 5.5)));
  }
  {
    SCOPED_TRACE("Bbox partially overlaps grid");
    checkBinCounts(valid, check);
  }

  // Try inserting something completely outside the index
  QPoint pt6 = QPoint::make_point(-2, 7.2);
  QPoint pt7 = QPoint::make_point(-0.1, 9.0);
  QBBox bbox4(pt6, pt7);
  valid.insert(bbox4, 6);
  {
    SCOPED_TRACE("Bbox completely misses grid");
    checkBinCounts(valid, check);
  }

  // Try inserting something completely overlapping the box
  QPoint pt8 = QPoint::make_point(-2, -0.2);
  QPoint pt9 = QPoint::make_point(7, 6.2);
  QBBox bbox5(pt8, pt9);
  valid.insert(bbox5, 8);
  for(int i = 0; i < valid.getNumBins(); ++i)
  {
    incr(check, i);
  }
  {
    SCOPED_TRACE("Bbox completely overlaps grid");
    checkBinCounts(valid, check);
  }
}

TEST(spin_uniform_grid, delete_stuff_2D)
{
  const int DIM = 2;
  using CoordType = double;
  using QPoint = axom::primal::Point<CoordType, DIM>;
  using QBBox = axom::primal::BoundingBox<CoordType, DIM>;

  double origin[DIM] = {0, 0};
  const int mpt = 6;
  double maxpoint[DIM] = {mpt, mpt};
  const int resolution = 6;
  int res[DIM] = {resolution, resolution};
  axom::spin::UniformGrid<int, DIM> valid(origin, maxpoint, res);

  std::map<int, int> check;

  // Insert something small (one bin), then clear its bin
  QPoint pt1 = QPoint::make_point(1.1, 1.1);
  QBBox bbox1(pt1);
  valid.insert(bbox1, 1);
  int index = valid.getBinIndex(pt1);
  incr(check, index);
  valid.clear(index);
  zero(check, index);
  {
    SCOPED_TRACE("Inserted one, then deleted it");
    checkBinCounts(valid, check);
  }

  // Insert some objects, clear some bins
  valid.insert(bbox1, 1);
  incr(check, valid.getBinIndex(pt1));
  QPoint pt2 = QPoint::make_point(0.1, 0.1);
  QPoint pt3 = QPoint::make_point(1.9, 2.3);
  QBBox bbox2(pt2, pt3);
  valid.insert(bbox2, 2);
  for(int j = 0; j < 3; ++j)
  {
    for(int i = 0; i < 2; ++i)
    {
      incr(check, valid.getBinIndex(QPoint::make_point(i + 0.5, j + 0.5)));
    }
  }
  {
    SCOPED_TRACE("Insert two, clear several bins");
    QPoint test = QPoint::make_point(0.5, 2.5);
    int tidx = valid.getBinIndex(test);
    valid.clear(tidx);
    zero(check, tidx);
    test = QPoint::make_point(0.5, 1.5);
    tidx = valid.getBinIndex(test);
    valid.clear(tidx);
    zero(check, tidx);
    test = QPoint::make_point(0.5, 2.5);
    tidx = valid.getBinIndex(test);
    valid.clear(tidx);
    zero(check, tidx);

    checkBinCounts(valid, check);
  }
  {
    SCOPED_TRACE("Try clearing the invalid bin");
    valid.clear(axom::spin::UniformGrid<int, DIM>::INVALID_BIN_INDEX);
    checkBinCounts(valid, check);
  }
}

namespace testing
{
template <typename ExecSpace>
struct ExecTraits
{
  static int getAllocatorId()
  {
#ifdef AXOM_USE_UMPIRE
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
#else
    return axom::getDefaultAllocatorID();
#endif
  }
  static int getUnifiedAllocatorId()
  {
#ifdef AXOM_USE_UMPIRE
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
#else
    return axom::getDefaultAllocatorID();
#endif
  }
};

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
template <int BLK_SZ>
struct ExecTraits<axom::CUDA_EXEC<BLK_SZ>>
{
  static int getAllocatorId()
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
  }
  static int getUnifiedAllocatorId()
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Unified);
  }
};
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
template <int BLK_SZ>
struct ExecTraits<axom::HIP_EXEC<BLK_SZ>>
{
  static int getAllocatorId()
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
  }
  static int getUnifiedAllocatorId()
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Unified);
  }
};
#endif
//------------------------------------------------------------------------------
//  This test harness defines some types that are useful for the tests below
//------------------------------------------------------------------------------
template <typename ExecSpace, int DIM>
class spin_uniform_grid_templated : public ::testing::Test
{
public:
  using PointType = axom::primal::Point<double, DIM>;
  using BoxType = axom::primal::BoundingBox<double, DIM>;

  using DynamicUniformGridType = axom::spin::UniformGrid<int, DIM, ExecSpace>;

  using FlatStorageType = axom::spin::policy::FlatGridStorage<int>;
  using FlatUniformGridType = axom::spin::UniformGrid<int, DIM, ExecSpace, FlatStorageType>;

  spin_uniform_grid_templated()
    : m_allocatorID(ExecTraits<ExecSpace>::getAllocatorId())
    , m_unifiedAllocatorID(ExecTraits<ExecSpace>::getUnifiedAllocatorId())
  { }

  void setNumObjects(int num_objects) { this->m_numObjects = num_objects; }

  void initialize()
  {
    this->m_boundingBoxes = axom::Array<BoxType>(0, this->m_numObjects, m_allocatorID);
    this->m_iota = axom::Array<int>(this->m_numObjects, this->m_numObjects, m_allocatorID);
    // generate n random objects to insert into the uniform grid.
    for(int i = 0; i < this->m_numObjects; i++)
    {
      PointType min, max;
      for(int idim = 0; idim < DIM; idim++)
      {
        double val1 = axom::utilities::random_real(0., 1.);
        double val2 = axom::utilities::random_real(0., 1.);
        if(val1 > val2)
        {
          axom::utilities::swap(val1, val2);
        }
        min[idim] = val1;
        max[idim] = val2;
      }
      this->m_boundingBoxes.push_back(BoxType(min, max));
    }

    const auto iota_v = m_iota.view();
    axom::for_all<ExecSpace>(this->m_numObjects, AXOM_LAMBDA(int idx) { iota_v[idx] = idx; });
  }

  std::unique_ptr<DynamicUniformGridType> constructUniformGridDynamic()
  {
    const int res_raw[3] = {2, 3, 4};
    const axom::NumericArray<int, DIM> resolution(res_raw);
    return std::make_unique<DynamicUniformGridType>(resolution,
                                                    m_boundingBoxes.view(),
                                                    m_iota.view(),
                                                    m_unifiedAllocatorID);
  }

  std::unique_ptr<FlatUniformGridType> constructUniformGridFlat()
  {
    const int res_raw[3] = {2, 3, 4};
    const axom::NumericArray<int, DIM> resolution(res_raw);
    return std::make_unique<FlatUniformGridType>(resolution,
                                                 m_boundingBoxes.view(),
                                                 m_iota.view(),
                                                 m_allocatorID);
  }

private:
  int m_allocatorID;
  int m_unifiedAllocatorID;

  int m_numObjects;
  axom::Array<BoxType> m_boundingBoxes;
  axom::Array<int> m_iota;
};

template <typename ExecSpace>
using UniformGrid2DTest = spin_uniform_grid_templated<ExecSpace, 2>;

template <typename ExecSpace>
using UniformGrid3DTest = spin_uniform_grid_templated<ExecSpace, 2>;

using MyTypes = ::testing::Types<
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  axom::OMP_EXEC,
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  axom::CUDA_EXEC<256>,
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  axom::HIP_EXEC<256>,
#endif
  axom::SEQ_EXEC>;

TYPED_TEST_SUITE(UniformGrid2DTest, MyTypes);
TYPED_TEST_SUITE(UniformGrid3DTest, MyTypes);

AXOM_TYPED_TEST(UniformGrid2DTest, initialize_tmpl)
{
  this->setNumObjects(1000);
  this->initialize();

  auto dynGrid = this->constructUniformGridDynamic();
  auto flatGrid = this->constructUniformGridFlat();
}

AXOM_TYPED_TEST(UniformGrid3DTest, initialize_tmpl)
{
  this->setNumObjects(1000);
  this->initialize();

  auto dynGrid = this->constructUniformGridDynamic();
  auto flatGrid = this->constructUniformGridFlat();
}

}  // namespace testing
