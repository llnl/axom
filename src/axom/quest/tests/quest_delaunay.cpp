// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/quest/Delaunay.hpp"
#include "axom/slic.hpp"

#include <cstddef>
#include <vector>

namespace
{

template <int DIM>
using DelaunayType = axom::quest::Delaunay<DIM>;

template <int DIM>
void insertPoints(DelaunayType<DIM>& dt,
                  const std::vector<typename DelaunayType<DIM>::PointType>& points)
{
  for(const auto& point : points)
  {
    dt.insertPoint(point);
  }
}

template <int DIM>
void expectValidDelaunay(DelaunayType<DIM>& dt,
                         const std::vector<typename DelaunayType<DIM>::PointType>& inserted_points,
                         int expected_num_elements = -1)
{
  dt.removeBoundary();

  EXPECT_TRUE(dt.getMeshData()->isValid(true));
  EXPECT_TRUE(dt.isValid(true));
  EXPECT_EQ(inserted_points.size(), static_cast<std::size_t>(dt.getMeshData()->vertices().size()));

  if(expected_num_elements >= 0)
  {
    EXPECT_EQ(expected_num_elements, static_cast<int>(dt.getMeshData()->elements().size()));
  }
  else
  {
    EXPECT_GT(static_cast<int>(dt.getMeshData()->elements().size()), 0);
  }
}

}  // namespace

TEST(quest_delaunay, cocircular_square_2d)
{
  using PointType = typename DelaunayType<2>::PointType;
  using BoundingBox = typename DelaunayType<2>::BoundingBox;

  DelaunayType<2> dt;
  dt.initializeBoundary(BoundingBox(PointType {-0.5, -0.5}, PointType {1.5, 1.5}));

  const std::vector<PointType> points {PointType {0.0, 0.0},
                                       PointType {1.0, 0.0},
                                       PointType {1.0, 1.0},
                                       PointType {0.0, 1.0}};

  insertPoints(dt, points);
  expectValidDelaunay(dt, points, 2);
}

TEST(quest_delaunay, regular_grid_2d)
{
  using PointType = typename DelaunayType<2>::PointType;
  using BoundingBox = typename DelaunayType<2>::BoundingBox;

  DelaunayType<2> dt;
  dt.initializeBoundary(BoundingBox(PointType {-1.0, -1.0}, PointType {3.0, 3.0}));

  std::vector<PointType> points;
  points.reserve(9);

  for(int y = 0; y < 3; ++y)
  {
    for(int x = 0; x < 3; ++x)
    {
      points.push_back(PointType {static_cast<double>(x), static_cast<double>(y)});
    }
  }

  insertPoints(dt, points);
  expectValidDelaunay(dt, points, 8);
}

TEST(quest_delaunay, boundary_location_regular_grid_2d)
{
  using PointType = typename DelaunayType<2>::PointType;
  using BoundingBox = typename DelaunayType<2>::BoundingBox;

  constexpr int NX = 20;
  constexpr int NY = 20;

  DelaunayType<2> dt;
  dt.initializeBoundary(BoundingBox(PointType {-1.0, -1.0}, PointType {2.0, 2.0}));

  std::vector<PointType> points;
  points.reserve(NX * NY);

  for(int y = 0; y < NY; ++y)
  {
    for(int x = 0; x < NX; ++x)
    {
      points.push_back(
        PointType {static_cast<double>(x) / (NX - 1), static_cast<double>(y) / (NY - 1)});
    }
  }

  insertPoints(dt, points);
  expectValidDelaunay(dt, points, 2 * (NX - 1) * (NY - 1));
}

TEST(quest_delaunay, cospherical_cube_3d)
{
  using PointType = typename DelaunayType<3>::PointType;
  using BoundingBox = typename DelaunayType<3>::BoundingBox;

  DelaunayType<3> dt;
  dt.initializeBoundary(BoundingBox(PointType {-0.5, -0.5, -0.5}, PointType {1.5, 1.5, 1.5}));

  std::vector<PointType> points;
  points.reserve(8);

  for(int z = 0; z < 2; ++z)
  {
    for(int y = 0; y < 2; ++y)
    {
      for(int x = 0; x < 2; ++x)
      {
        points.push_back(
          PointType {static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)});
      }
    }
  }

  insertPoints(dt, points);
  expectValidDelaunay(dt, points);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
