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

template <int DIM>
void expectConformingMesh(const DelaunayType<DIM>& dt)
{
  EXPECT_TRUE(dt.getMeshData()->isValid(true));
  EXPECT_TRUE(dt.getMeshData()->isConforming(true));
  EXPECT_TRUE(dt.isConforming(true));
}

template <int DIM>
void expectValidDelaunayWithBoundary(const DelaunayType<DIM>& dt)
{
  EXPECT_TRUE(dt.getMeshData()->isValid(true));
  EXPECT_TRUE(dt.getMeshData()->isConforming(true));
  EXPECT_TRUE(dt.isConforming(true));
  EXPECT_TRUE(dt.isValid(true));
}

}  // namespace

TEST(quest_delaunay, local_bad_geometry_pool_2d)
{
  using PointType = typename DelaunayType<2>::PointType;
  using BoundingBox = typename DelaunayType<2>::BoundingBox;
  using ValidationMode = typename DelaunayType<2>::InsertionValidationMode;

  const std::vector<PointType> points {
    PointType {0.72835559683934026, 0.82937230204996126},
    PointType {0.62600920572645069, 0.84649465468911222},
    PointType {0.057066612786023117, 0.98518528899105051},
    PointType {0.50420652587117898, 0.88131216460192052},
    PointType {0.64623860530040367, 0.84784930386894053},
    PointType {0.63047324496924584, 0.8517482773644065},
    PointType {0.61321189052614333, 0.8508361023097063},
    PointType {0.6153686267133367, 0.85060225303567505},
    PointType {0.57633479792557119, 0.86309485101808103},
    PointType {0.5734747290517247, 0.85999079345999363},
    PointType {0.59008692995299927, 0.85528988587964816},
    PointType {0.59200361573061988, 0.86005025488638098},
    PointType {0.60412589378177151, 0.85608666941414724},
    PointType {0.59008239905455229, 0.86082490771226516},
    PointType {0.59108093308168341, 0.85843584821320762},
    PointType {0.63956447880378053, 0.849279316704322},
    PointType {0.55233739946809535, 0.86640226672486664},
  };

  DelaunayType<2> dt;
  dt.initializeBoundary(BoundingBox(PointType {0., 0.}, PointType {1., 1.}));
  dt.setInsertionValidationMode(ValidationMode::Local);
  insertPoints(dt, std::vector<PointType>(points.begin(), points.end() - 1));
  expectValidDelaunayWithBoundary(dt);

  dt.insertPoint(points.back());
  expectValidDelaunayWithBoundary(dt);
}

TEST(quest_delaunay, local_bad_geometry_pool_3d)
{
  using PointType = typename DelaunayType<3>::PointType;
  using BoundingBox = typename DelaunayType<3>::BoundingBox;

  const std::vector<PointType> points {
    PointType {0.51599885628789743, 0.75470155608398004, 0.40855633126704233},
    PointType {0.50343513905278869, 0.76256644865428691, 0.41550002353843529},
    PointType {0.52181994461717562, 0.75507132048650449, 0.40434771160709904},
    PointType {0.51222305794382439, 0.76931259688762987, 0.41173021886913774},
    PointType {0.51893390889491309, 0.7608695730433962, 0.39023418626911949},
    PointType {0.49870963994503914, 0.73975835990074734, 0.40518385970598131},
    PointType {0.52105735259915775, 0.77213020101602725, 0.40722478047503557},
    PointType {0.49002686184837962, 0.75462235742107142, 0.41591466571185198},
    PointType {0.50501491031893275, 0.76598310556724625, 0.39817919534292223},
    PointType {0.48943926088386502, 0.75900280017153021, 0.40203090542081277},
    PointType {0.5113764391413268, 0.76710422391570976, 0.40773506111442864},
    PointType {0.51178826597709415, 0.77818416231550269, 0.3925790516462227},
    PointType {0.52301040315781211, 0.77096219982230552, 0.40733482938935728},
    PointType {0.53346303899079406, 0.77993659887621847, 0.39369144158113362},
    PointType {0.50029750641338455, 0.76954990769049825, 0.39799982472078405},
    PointType {0.52770121131451853, 0.76483515658193735, 0.40353420628419612},
    PointType {0.50952730554172287, 0.74599842639434555, 0.41708835434555841},
    PointType {0.52956607163853775, 0.75418313145478433, 0.40980177397908779},
    PointType {0.51225516061546172, 0.75892162524570184, 0.38822059123401836},
    PointType {0.51473015232418118, 0.74851193285704143, 0.39787104734305007},
    PointType {0.51352196009694051, 0.76305382863756566, 0.39295743334708788},
    PointType {0.51822520982029341, 0.76176127807719429, 0.42723370773918606},
    PointType {0.51440345815174604, 0.7705063813945725, 0.41385100134174907},
    PointType {0.50863912619096874, 0.76214178681356526, 0.39419055286657634},
    PointType {0.51550017071903309, 0.75494779413607704, 0.39942038195430957},
    PointType {0.51299271231867682, 0.76122513269993519, 0.4067726415350541},
  };

  DelaunayType<3> dt;
  dt.initializeBoundary(BoundingBox(PointType {0., 0., 0.}, PointType {1., 1., 1.}));
  insertPoints(dt, std::vector<PointType>(points.begin(), points.end() - 1));
  expectValidDelaunayWithBoundary(dt);

  dt.insertPoint(points.back());
  expectValidDelaunayWithBoundary(dt);
}

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

TEST(quest_delaunay, insertion_validation_regular_grid_3d)
{
  using PointType = typename DelaunayType<3>::PointType;
  using BoundingBox = typename DelaunayType<3>::BoundingBox;
  using ValidationMode = typename DelaunayType<3>::InsertionValidationMode;

  DelaunayType<3> dt;
  dt.initializeBoundary(BoundingBox(PointType {-0.5, -0.5, -0.5}, PointType {1.5, 1.5, 1.5}));
  dt.setInsertionValidationMode(ValidationMode::ConformingMesh);

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
  expectConformingMesh(dt);
  expectValidDelaunay(dt, points);
}

TEST(quest_delaunay, insertion_validation_full_small_grid_3d)
{
  using PointType = typename DelaunayType<3>::PointType;
  using BoundingBox = typename DelaunayType<3>::BoundingBox;
  using ValidationMode = typename DelaunayType<3>::InsertionValidationMode;

  DelaunayType<3> dt;
  dt.initializeBoundary(BoundingBox(PointType {-0.5, -0.5, -0.5}, PointType {1.5, 1.5, 1.5}));
  dt.setInsertionValidationMode(ValidationMode::Full);

  std::vector<PointType> points;
  points.reserve(3 * 3 * 3);
  for(int z = 0; z < 3; ++z)
  {
    for(int y = 0; y < 3; ++y)
    {
      for(int x = 0; x < 3; ++x)
      {
        points.push_back(PointType {x / 2., y / 2., z / 2.});
      }
    }
  }

  insertPoints(dt, points);
  expectConformingMesh(dt);
  expectValidDelaunay(dt, points);
}

TEST(quest_delaunay, stress_nearly_coplanar_insertions_3d)
{
  using PointType = typename DelaunayType<3>::PointType;
  using BoundingBox = typename DelaunayType<3>::BoundingBox;

  // A near-coplanar point set (most points lie close to a plane) tends to
  // generate sliver tetrahedra and near-zero orientation / in-sphere determinants
  // This is a good stress case for point-location and cavity predicates in 3D.
  DelaunayType<3> dt;
  dt.initializeBoundary(BoundingBox(PointType {-1., -1., -1.}, PointType {2., 2., 2.}));

  std::vector<PointType> points;
  points.reserve(5 * 5 + 2);

  // Keep the grid close to planar, but not so close that the geometric
  // circumsphere validator is dominated by sliver-conditioning noise.
  constexpr double eps = 1e-5;
  for(int y = 0; y < 5; ++y)
  {
    for(int x = 0; x < 5; ++x)
    {
      const double fx = static_cast<double>(x) / 4.;
      const double fy = static_cast<double>(y) / 4.;
      points.push_back(PointType {fx, fy, eps * (fx + 2. * fy)});
    }
  }

  // Add a couple of off-plane points to ensure a full 3D convex hull.
  points.push_back(PointType {0.2, 0.2, 0.4});
  points.push_back(PointType {0.8, 0.6, 0.7});

  insertPoints(dt, points);
  expectValidDelaunay(dt, points);
}

TEST(quest_delaunay, stress_large_coordinate_scale_grid_3d)
{
  using PointType = typename DelaunayType<3>::PointType;
  using BoundingBox = typename DelaunayType<3>::BoundingBox;
  using ValidationMode = typename DelaunayType<3>::InsertionValidationMode;

  // Large absolute coordinates with small relative spacing increase the risk
  // of cancellation in determinant-based predicates. This is a deterministic
  // stress case for the 3D orientation / in-sphere computations.
  constexpr double base = 1e9;
  constexpr double h = 1.0;

  DelaunayType<3> dt;
  dt.initializeBoundary(BoundingBox(PointType {base - 1., base - 1., base - 1.},
                                    PointType {base + 4., base + 4., base + 4.}));
  dt.setInsertionValidationMode(ValidationMode::Full);

  std::vector<PointType> points;
  points.reserve(4 * 4 * 4);
  for(int z = 0; z < 4; ++z)
  {
    for(int y = 0; y < 4; ++y)
    {
      for(int x = 0; x < 4; ++x)
      {
        points.push_back(PointType {base + h * x, base + h * y, base + h * z});
      }
    }
  }

  insertPoints(dt, points);
  expectConformingMesh(dt);
  expectValidDelaunay(dt, points);
}

TEST(quest_delaunay, stress_large_coordinate_scale_grid_2d)
{
  using PointType = typename DelaunayType<2>::PointType;
  using BoundingBox = typename DelaunayType<2>::BoundingBox;
  using ValidationMode = typename DelaunayType<2>::InsertionValidationMode;

  // Large absolute coordinates with small relative spacing increase the risk
  // of cancellation in determinant-based predicates. This stress case mirrors
  // the 3D test but exercises the 2D in-circle predicate.
  constexpr double base = 1e9;
  constexpr double h = 1.0;

  DelaunayType<2> dt;
  dt.initializeBoundary(
    BoundingBox(PointType {base - 1., base - 1.}, PointType {base + 4., base + 4.}));
  dt.setInsertionValidationMode(ValidationMode::Full);

  std::vector<PointType> points;
  points.reserve(4 * 4);
  for(int y = 0; y < 4; ++y)
  {
    for(int x = 0; x < 4; ++x)
    {
      points.push_back(PointType {base + h * x, base + h * y});
    }
  }

  insertPoints(dt, points);
  expectConformingMesh(dt);
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
