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

TEST(quest_delaunay, boundary_location_regular_grid_3d)
{
  using PointType = typename DelaunayType<3>::PointType;
  using BoundingBox = typename DelaunayType<3>::BoundingBox;

  const double one_third = 1. / 3.;
  const double two_thirds = 2. / 3.;

  DelaunayType<3> dt;
  dt.initializeBoundary(BoundingBox(PointType {-0.5, -0.5, -0.5}, PointType {1.5, 1.5, 1.5}));

  std::vector<PointType> inserted_points {PointType {0., two_thirds, 1.},
                                          PointType {0., two_thirds, one_third},
                                          PointType {two_thirds, two_thirds, one_third},
                                          PointType {two_thirds, two_thirds, two_thirds},
                                          PointType {one_third, one_third, one_third},
                                          PointType {0., 1., 1.},
                                          PointType {one_third, 1., one_third},
                                          PointType {two_thirds, 1., one_third},
                                          PointType {0., 0., two_thirds},
                                          PointType {one_third, two_thirds, two_thirds},
                                          PointType {1., 1., 0.},
                                          PointType {0., 0., 1.},
                                          PointType {one_third, 0., 1.},
                                          PointType {0., one_third, one_third},
                                          PointType {two_thirds, 0., one_third},
                                          PointType {one_third, two_thirds, one_third}};

  insertPoints(dt, inserted_points);

  const PointType query_pt {one_third, 0., two_thirds};
  EXPECT_NE(DelaunayType<3>::INVALID_INDEX, dt.findContainingElement(query_pt));

  inserted_points.push_back(query_pt);
  dt.insertPoint(query_pt);

  expectValidDelaunay(dt, inserted_points);
}

TEST(quest_delaunay, query_location_regular_grid_3d)
{
  using PointType = typename DelaunayType<3>::PointType;
  using BoundingBox = typename DelaunayType<3>::BoundingBox;

  DelaunayType<3> dt;
  dt.initializeBoundary(BoundingBox(PointType {-0.5, -0.5, -0.5}, PointType {1.5, 1.5, 1.5}));

  std::vector<PointType> points;
  points.reserve(4 * 4 * 4);
  for(int z = 0; z < 4; ++z)
  {
    for(int y = 0; y < 4; ++y)
    {
      for(int x = 0; x < 4; ++x)
      {
        points.push_back(PointType {x / 3., y / 3., z / 3.});
      }
    }
  }

  insertPoints(dt, points);
  dt.removeBoundary();

  const std::vector<PointType> queries {PointType {0.445693, 0.321066, 0.0526028},
                                        PointType {0.381996, 0.321531, 0.0811576},
                                        PointType {0.513459, 0.311519, 0.501756},
                                        PointType {0.499934, 0.32269, 0.0414115},
                                        PointType {0.541159, 0.307432, 0.546373},
                                        PointType {0.553746, 0.368503, 0.151107},
                                        PointType {0.514516, 0.372636, 0.575886},
                                        PointType {0.649638, 0.366708, 0.596617}};

  for(const auto& query : queries)
  {
    EXPECT_NE(DelaunayType<3>::INVALID_INDEX, dt.findContainingElement(query, false));
  }
}

TEST(quest_delaunay, query_outside_convex_hull_returns_invalid_3d)
{
  using PointType = typename DelaunayType<3>::PointType;
  using BoundingBox = typename DelaunayType<3>::BoundingBox;

  DelaunayType<3> dt;
  dt.initializeBoundary(BoundingBox(PointType {0., 0., 0.}, PointType {1., 1., 1.}));

  const std::vector<PointType> points {PointType {0.2, 0.2, 0.2},
                                       PointType {0.8, 0.2, 0.2},
                                       PointType {0.2, 0.8, 0.2},
                                       PointType {0.2, 0.2, 0.8}};

  insertPoints(dt, points);
  dt.removeBoundary();

  const PointType inside_bbox_outside_hull {0.8, 0.8, 0.8};
  EXPECT_EQ(DelaunayType<3>::INVALID_INDEX,
            dt.findContainingElement(inside_bbox_outside_hull, false));
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
