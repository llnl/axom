// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/OrientationResult.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/operators/orientation.hpp"

#include "axom/fmt.hpp"

#include <math.h>
#include <algorithm>  // std::next_permutation

namespace primal = axom::primal;

/// Test fixture for testing primal::Tetrahedron
class TetrahedronTest : public ::testing::Test
{
public:
  static const int DIM = 3;

  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QTet = primal::Tetrahedron<CoordType, DIM>;
  using QTri = primal::Triangle<CoordType, DIM>;
  using QSeg = primal::Segment<CoordType, DIM>;

protected:
  virtual void SetUp()
  {
    EPS = 1e-12;

    // Define coordinates for first tetrahedron
    qData0[0] = QPoint {0, 0, 0};
    qData0[1] = QPoint {1, 0, 0};
    qData0[2] = QPoint {1, 1, 0};
    qData0[3] = QPoint {1, 1, 1};

    // Define coordinates for second tetrahedron
    qData1[0] = QPoint {1, 0, 0};
    qData1[1] = QPoint {0, 1, 0};
    qData1[2] = QPoint {0, 0, 0};
    qData1[3] = QPoint {0, 0, 1};

    double angles[3];
    for(int i = 0; i < 3; ++i)
    {
      angles[i] = i * M_PI / 3.;
    }

    // Define coordinates for third tetrahedron
    {
      const double sc2 = .1;
      qData2[0] = QPoint {sc2 * std::cos(angles[0]), sc2 * std::sin(angles[0]), 0};
      qData2[1] = QPoint {sc2 * std::cos(angles[1]), sc2 * std::sin(angles[1]), 0};
      qData2[2] = QPoint {sc2 * std::cos(angles[2]), sc2 * std::sin(angles[2]), 0};
      qData2[3] = QPoint {0, 0, 100.};
    }

    // Define coordinates for fourth tetrahedron
    {
      const double sc3 = 100.;
      qData3[0] = QPoint {sc3 * std::cos(angles[0]), sc3 * std::sin(angles[0]), 0};
      qData3[1] = QPoint {sc3 * std::cos(angles[1]), sc3 * std::sin(angles[1]), 0};
      qData3[2] = QPoint {sc3 * std::cos(angles[2]), sc3 * std::sin(angles[2]), 0};
      qData3[3] = QPoint {0, 0, .1};
    }

    // A regular tetrahedron centered at the origin
    {
      // define some constants for regular tetrahedron
      const double c[3] = {std::sqrt(2) / 3, std::sqrt(6) / 3, 1. / 3};
      qData4[0] = QPoint {0, 0, 1};
      qData4[1] = QPoint {2 * c[0], 0, -c[2]};
      qData4[2] = QPoint {-c[0], -c[1], -c[2]};
      qData4[3] = QPoint {-c[0], c[1], -c[2]};
    }
  }

  int numTetrahedra() const { return 5; }

  QTet getTet(int idx)
  {
    EXPECT_TRUE(idx >= 0 && idx < numTetrahedra());

    QTet tet;

    switch(idx)
    {
    case 0:
      tet = QTet(qData0[0], qData0[1], qData0[2], qData0[3]);
      break;
    case 1:
      tet = QTet(qData1[0], qData1[1], qData1[2], qData1[3]);
      break;
    case 2:
      tet = QTet(qData2[0], qData2[1], qData2[2], qData2[3]);
      break;
    case 3:
      tet = QTet(qData3[0], qData3[1], qData3[2], qData3[3]);
      break;
    case 4:
      tet = QTet(qData4[0], qData4[1], qData4[2], qData4[3]);
      break;
    }

    return tet;
  }

  QPoint qData0[4];
  QPoint qData1[4];
  QPoint qData2[4];
  QPoint qData3[4];
  QPoint qData4[4];
  double EPS;
};

//------------------------------------------------------------------------------
TEST_F(TetrahedronTest, defaultConstructor)
{
  using QPoint = TetrahedronTest::QPoint;
  using QTet = TetrahedronTest::QTet;

  const QTet tet;

  // Test ostream operator
  SLIC_INFO("Empty tetrahedron coordinates: " << tet);

  // Check indirection operator
  EXPECT_EQ(QPoint::zero(), tet[0]);
  EXPECT_EQ(QPoint::zero(), tet[1]);
  EXPECT_EQ(QPoint::zero(), tet[2]);
  EXPECT_EQ(QPoint::zero(), tet[3]);
}

TEST_F(TetrahedronTest, constructFromPoints)
{
  using QPoint = TetrahedronTest::QPoint;
  using QTet = TetrahedronTest::QTet;

  // Access the test data
  const QPoint* pt = this->qData0;

  axom::Array<QPoint> ptArray({pt[0], pt[1], pt[2], pt[3]});

  QTet tet1(pt[0], pt[1], pt[2], pt[3]);

  QTet tet2(pt);

  QTet tet3(ptArray);

  QTet tet4({pt[0], pt[1], pt[2], pt[3]});

  // Test ostream operator
  SLIC_INFO("Tetrahedron 1 coordinates: " << tet1);
  SLIC_INFO("Tetrahedron 2 coordinates: " << tet2);
  SLIC_INFO("Tetrahedron 3 coordinates: " << tet3);
  SLIC_INFO("Tetrahedron 4 coordinates: " << tet4);

  // Check indirection operator
  EXPECT_EQ(pt[0], tet1[0]);
  EXPECT_EQ(pt[1], tet1[1]);
  EXPECT_EQ(pt[2], tet1[2]);
  EXPECT_EQ(pt[3], tet1[3]);

  EXPECT_EQ(pt[0], tet2[0]);
  EXPECT_EQ(pt[1], tet2[1]);
  EXPECT_EQ(pt[2], tet2[2]);
  EXPECT_EQ(pt[3], tet2[3]);

  EXPECT_EQ(pt[0], tet3[0]);
  EXPECT_EQ(pt[1], tet3[1]);
  EXPECT_EQ(pt[2], tet3[2]);
  EXPECT_EQ(pt[3], tet3[3]);

  EXPECT_EQ(pt[0], tet4[0]);
  EXPECT_EQ(pt[1], tet4[1]);
  EXPECT_EQ(pt[2], tet4[2]);
  EXPECT_EQ(pt[3], tet4[3]);
}

TEST_F(TetrahedronTest, volume)
{
  using QPoint = TetrahedronTest::QPoint;
  using QTet = TetrahedronTest::QTet;

  // Access the test data
  const QPoint* pt = this->qData0;
  QTet tet(pt[0], pt[1], pt[2], pt[3]);

  double expVolume = 1. / 6.;
  EXPECT_EQ(expVolume, tet.signedVolume());
  EXPECT_EQ(tet.signedVolume(), tet.volume());
}

TEST_F(TetrahedronTest, degenerate)
{
  using QPoint = TetrahedronTest::QPoint;
  using QTet = TetrahedronTest::QTet;

  // Access the test data
  const QPoint* pt = this->qData0;
  QTet tet(pt[0], pt[1], pt[2], pt[3]);

  EXPECT_FALSE(tet.degenerate());

  // Make the tet degenerate by identifying two vertices
  tet[0] = tet[1];
  EXPECT_TRUE(tet.degenerate());
}

TEST_F(TetrahedronTest, barycentric)
{
  using CoordType = TetrahedronTest::CoordType;
  using QPoint = TetrahedronTest::QPoint;
  using QTet = TetrahedronTest::QTet;
  using RPoint = primal::Point<CoordType, 4>;
  using TestVec = std::vector<std::pair<QPoint, RPoint>>;

  const QPoint* pt = this->qData1;
  QTet tet(pt[0], pt[1], pt[2], pt[3]);

  TestVec testData;

  // Test the four vertices
  testData.push_back(std::make_pair(pt[0], RPoint {1., 0., 0., 0.}));
  testData.push_back(std::make_pair(pt[1], RPoint {0., 1., 0., 0.}));
  testData.push_back(std::make_pair(pt[2], RPoint {0., 0., 1., 0.}));
  testData.push_back(std::make_pair(pt[3], RPoint {0., 0., 0., 1.}));

  // Test some of the edge midpoints
  testData.push_back(std::make_pair(QPoint::midpoint(pt[0], pt[1]), RPoint {0.5, 0.5, 0., 0.}));
  testData.push_back(std::make_pair(QPoint::midpoint(pt[1], pt[2]), RPoint {0., 0.5, 0.5, 0.}));
  testData.push_back(std::make_pair(QPoint::midpoint(pt[2], pt[3]), RPoint {0., 0., 0.5, 0.5}));
  testData.push_back(std::make_pair(QPoint::midpoint(pt[0], pt[2]), RPoint {0.5, 0., 0.5, 0.}));
  testData.push_back(std::make_pair(QPoint::midpoint(pt[0], pt[3]), RPoint {0.5, 0., 0., 0.5}));
  testData.push_back(std::make_pair(QPoint::midpoint(pt[1], pt[3]), RPoint {0., 0.5, 0., 0.5}));

  // Test the centroid
  testData.push_back(
    std::make_pair(QPoint(.25 * (pt[0].array() + pt[1].array() + pt[2].array() + pt[3].array())),
                   RPoint {.25, .25, .25, .25}));

  // Test a point outside the tetrahedron
  testData.push_back(
    std::make_pair(QPoint(-0.4 * pt[0].array() + 1.2 * pt[1].array() + 0.2 * pt[2].array()),
                   RPoint {-0.4, 1.2, 0.2, 0.}));

  // Now run the actual tests
  for(const auto& data : testData)
  {
    const QPoint& query = data.first;
    const RPoint& expBary = data.second;
    RPoint bary = tet.physToBarycentric(query);

    SLIC_DEBUG(
      axom::fmt::format("Computed barycentric coordinates for tetrahedron {} and point {} are {}",
                        tet,
                        query,
                        bary));

    EXPECT_NEAR(bary[0], expBary[0], this->EPS);
    EXPECT_NEAR(bary[1], expBary[1], this->EPS);
    EXPECT_NEAR(bary[2], expBary[2], this->EPS);
    EXPECT_NEAR(bary[3], expBary[3], this->EPS);
  }
}

//------------------------------------------------------------------------------
TEST_F(TetrahedronTest, barycentric_skipNormalization)
{
  using CoordType = TetrahedronTest::CoordType;
  using QPoint = TetrahedronTest::QPoint;
  using RPoint = primal::Point<CoordType, 4>;
  using TestVec = std::vector<std::pair<QPoint, RPoint>>;

  for(int i = 0; i < this->numTetrahedra(); ++i)
  {
    const auto tet = this->getTet(i);
    QPoint pt[4] = {tet[0], tet[1], tet[2], tet[3]};

    TestVec testData;

    // Test the four vertices
    testData.push_back(std::make_pair(pt[0], RPoint {1., 0., 0., 0.}));
    testData.push_back(std::make_pair(pt[1], RPoint {0., 1., 0., 0.}));
    testData.push_back(std::make_pair(pt[2], RPoint {0., 0., 1., 0.}));
    testData.push_back(std::make_pair(pt[3], RPoint {0., 0., 0., 1.}));

    // Test the edge midpoints
    testData.push_back(std::make_pair(QPoint::midpoint(pt[0], pt[1]), RPoint {0.5, 0.5, 0., 0.}));
    testData.push_back(std::make_pair(QPoint::midpoint(pt[1], pt[2]), RPoint {0., 0.5, 0.5, 0.}));
    testData.push_back(std::make_pair(QPoint::midpoint(pt[2], pt[3]), RPoint {0., 0., 0.5, 0.5}));
    testData.push_back(std::make_pair(QPoint::midpoint(pt[0], pt[2]), RPoint {0.5, 0., 0.5, 0.}));
    testData.push_back(std::make_pair(QPoint::midpoint(pt[0], pt[3]), RPoint {0.5, 0., 0., 0.5}));
    testData.push_back(std::make_pair(QPoint::midpoint(pt[1], pt[3]), RPoint {0., 0.5, 0., 0.5}));

    // Test the face midpoints
    constexpr double one_third = 1. / 3.;
    testData.push_back(
      std::make_pair(QPoint(one_third * (pt[0].array() + pt[1].array() + pt[2].array())),
                     RPoint {one_third, one_third, one_third, 0.}));
    testData.push_back(
      std::make_pair(QPoint(one_third * (pt[0].array() + pt[1].array() + pt[3].array())),
                     RPoint {one_third, one_third, 0., one_third}));
    testData.push_back(
      std::make_pair(QPoint(one_third * (pt[0].array() + pt[2].array() + pt[3].array())),
                     RPoint {one_third, 0., one_third, one_third}));
    testData.push_back(
      std::make_pair(QPoint(one_third * (pt[1].array() + pt[2].array() + pt[3].array())),
                     RPoint {0., one_third, one_third, one_third}));

    // Test the centroid
    testData.push_back(
      std::make_pair(QPoint(.25 * (pt[0].array() + pt[1].array() + pt[2].array() + pt[3].array())),
                     RPoint {.25, .25, .25, .25}));

    // Test a point outside the tetrahedron
    testData.push_back(
      std::make_pair(QPoint(-0.4 * pt[0].array() + 1.2 * pt[1].array() + 0.2 * pt[2].array()),
                     RPoint {-0.4, 1.2, 0.2, 0.}));

    // Now run the actual tests
    for(const auto& data : testData)
    {
      const QPoint& query = data.first;
      const RPoint& expBary = data.second;
      RPoint bary = tet.physToBarycentric(query, false);
      RPoint baryUnnormalized = tet.physToBarycentric(query, true);

      // The unnormalized weights are proportional to the actual barycentric weights
      // The factor of 6 is due to the use of parallelepiped volumes instead of tet volumes
      const double volumeScale = 6 * tet.signedVolume();

      SLIC_INFO(
        axom::fmt::format("For tet {} and point {} "
                          "-- barycentric coods {} (unnormalized barycentric {})"
                          "-- signed volume {}",
                          tet,
                          query,
                          bary,
                          baryUnnormalized,
                          tet.signedVolume()));

      EXPECT_NEAR(volumeScale, baryUnnormalized.array().sum(), this->EPS);
      for(int d = 0; d <= 3; ++d)
      {
        EXPECT_NEAR(bary[d] * volumeScale, baryUnnormalized[d], this->EPS);
        EXPECT_NEAR(expBary[d] * volumeScale, baryUnnormalized[d], this->EPS);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(TetrahedronTest, tetrahedron_roundtrip_bary_to_physical)
{
  const double EPS = 1e-12;

  using CoordType = TetrahedronTest::CoordType;
  using QPoint = TetrahedronTest::QPoint;
  using QTet = TetrahedronTest::QTet;
  using RPoint = primal::Point<CoordType, QTet::NUM_VERTS>;

  // Test tets
  std::vector<QTet> tets = {this->getTet(0),
                            this->getTet(1),
                            this->getTet(2),
                            this->getTet(3),
                            this->getTet(4)};

  // Compute circumsphere of test triangles and test some points
  for(const auto& tet : tets)
  {
    // test vertices
    {
      RPoint b_in[4] = {RPoint {1., 0., 0., 0.},
                        RPoint {0., 1., 0., 0.},
                        RPoint {0., 0., 1., 0.},
                        RPoint {0., 0., 0., 1.}};

      QPoint p_exp[4] = {tet[0], tet[1], tet[2], tet[3]};

      for(int i = 0; i < 4; ++i)
      {
        QPoint b2p = tet.baryToPhysical(b_in[i]);
        EXPECT_NEAR(0., primal::squared_distance(p_exp[i], b2p), EPS);

        RPoint p2b = tet.physToBarycentric(b2p);
        EXPECT_NEAR(0., primal::squared_distance(b_in[i], p2b), EPS);
      }

      // test edges
      {
        RPoint b_in[6] = {RPoint {.5, .5, 0., 0.},
                          RPoint {.5, 0., .5, 0.},
                          RPoint {.5, 0., 0., .5},
                          RPoint {0., .5, .5, 0.},
                          RPoint {0., .5, 0., .5},
                          RPoint {0., 0., .5, .5}};

        QPoint p_exp[6] = {QPoint::midpoint(tet[0], tet[1]),
                           QPoint::midpoint(tet[0], tet[2]),
                           QPoint::midpoint(tet[0], tet[3]),
                           QPoint::midpoint(tet[1], tet[2]),
                           QPoint::midpoint(tet[1], tet[3]),
                           QPoint::midpoint(tet[2], tet[3])};

        for(int i = 0; i < 6; ++i)
        {
          QPoint b2p = tet.baryToPhysical(b_in[i]);
          EXPECT_NEAR(0., primal::squared_distance(p_exp[i], b2p), EPS);

          RPoint p2b = tet.physToBarycentric(b2p);
          EXPECT_NEAR(0., primal::squared_distance(b_in[i], p2b), EPS);
        }
      }

      // test face barycenters
      {
        const double third = 1. / 3.;
        RPoint b_in[4] = {RPoint {third, third, third, 0.},
                          RPoint {third, third, 0., third},
                          RPoint {third, 0., third, third},
                          RPoint {0., third, third, third}};

        QPoint p_exp[4] = {QPoint(third * (tet[0].array() + tet[1].array() + tet[2].array())),
                           QPoint(third * (tet[0].array() + tet[1].array() + tet[3].array())),
                           QPoint(third * (tet[0].array() + tet[2].array() + tet[3].array())),
                           QPoint(third * (tet[1].array() + tet[2].array() + tet[3].array()))};

        for(int i = 0; i < 4; ++i)
        {
          QPoint b2p = tet.baryToPhysical(b_in[i]);
          EXPECT_NEAR(0., primal::squared_distance(p_exp[i], b2p), EPS);

          RPoint p2b = tet.physToBarycentric(b2p);
          EXPECT_NEAR(0., primal::squared_distance(b_in[i], p2b), EPS);
        }
      }

      // test tet barycenters
      {
        RPoint b_in[1] = {RPoint {.25, .25, .25, .25}};

        QPoint p_exp[1] = {
          QPoint(.25 * (tet[0].array() + tet[1].array() + tet[2].array() + tet[3].array()))};

        for(int i = 0; i < 1; ++i)
        {
          QPoint b2p = tet.baryToPhysical(b_in[i]);
          EXPECT_NEAR(0., primal::squared_distance(p_exp[i], b2p), EPS);

          RPoint p2b = tet.physToBarycentric(b2p);
          EXPECT_NEAR(0., primal::squared_distance(b_in[i], p2b), EPS);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(TetrahedronTest, tet_3D_circumsphere)
{
  using CoordType = TetrahedronTest::CoordType;
  using QPoint = TetrahedronTest::QPoint;
  using QTet = TetrahedronTest::QTet;
  using QSphere = primal::Sphere<CoordType, 3>;
  using RPoint = primal::Point<CoordType, 4>;
  const double EPS = 1e-9;

  using primal::ON_BOUNDARY;
  using primal::ON_NEGATIVE_SIDE;
  using primal::ON_POSITIVE_SIDE;

  // Compute circumsphere of test tetrahedra and test some points
  for(int ti = 0; ti < this->numTetrahedra(); ++ti)
  {
    QTet tet = this->getTet(ti);
    QSphere circumsphere = tet.circumsphere();

    SLIC_INFO("Circumsphere for tetrahedron: " << tet << " is " << circumsphere);

    // check that each vertex is on the sphere
    for(int i = 0; i < 4; ++i)
    {
      auto qpt = tet[i];
      EXPECT_NEAR(circumsphere.getRadius(),
                  sqrt(primal::squared_distance(qpt, circumsphere.getCenter())),
                  EPS);
    }

    // test vertices
    for(int i = 0; i < 4; ++i)
    {
      QPoint qpt = tet[i];
      EXPECT_EQ(ON_BOUNDARY, circumsphere.getOrientation(qpt, EPS));
    }

    // test edge centers
    {
      QPoint qpt[6] = {QPoint::midpoint(tet[0], tet[1]),
                       QPoint::midpoint(tet[0], tet[2]),
                       QPoint::midpoint(tet[0], tet[3]),
                       QPoint::midpoint(tet[1], tet[2]),
                       QPoint::midpoint(tet[1], tet[3]),
                       QPoint::midpoint(tet[2], tet[3])};
      for(int j = 0; j < 6; ++j)
      {
        EXPECT_EQ(ON_NEGATIVE_SIDE, circumsphere.getOrientation(qpt[j], EPS));
      }
    }

    // test face centers
    {
      constexpr CoordType third = 1. / 3.;
      constexpr CoordType zero {0};
      QPoint qpt[4] = {tet.baryToPhysical(RPoint {third, third, third, zero}),
                       tet.baryToPhysical(RPoint {third, third, zero, third}),
                       tet.baryToPhysical(RPoint {third, zero, third, third}),
                       tet.baryToPhysical(RPoint {zero, third, third, third})};
      for(int j = 0; j < 4; ++j)
      {
        EXPECT_EQ(ON_NEGATIVE_SIDE, circumsphere.getOrientation(qpt[j], EPS));
      }
    }

    // test tet center
    {
      QPoint qpt = tet.baryToPhysical(RPoint {.25, .25, .25, .25});
      EXPECT_EQ(ON_NEGATIVE_SIDE, circumsphere.getOrientation(qpt, EPS));
    }

    // test points that should be far outside tet
    {
      QPoint qpt[4] = {tet.baryToPhysical(RPoint {-1, 3, -1, 0}),
                       tet.baryToPhysical(RPoint {0, -1, 3, -1}),
                       tet.baryToPhysical(RPoint {-1, -1, 0, 3}),
                       tet.baryToPhysical(RPoint {3, -1, -1, 0})};
      for(int j = 0; j < 4; ++j)
      {
        EXPECT_EQ(ON_POSITIVE_SIDE, circumsphere.getOrientation(qpt[j], EPS));
      }
    }
  }
}

TEST_F(TetrahedronTest, regularTetrahedron)
{
  using QSeg = TetrahedronTest::QSeg;
  using QTri = TetrahedronTest::QTri;

  // get the regular tetrahedron
  auto tet = this->getTet(4);
  SLIC_INFO("Regular tetrahedron: " << tet);

  const double exp_edge_len = 2. * std::sqrt(6) / 3;
  const double exp_vol = 8 * sqrt(3) / 27;

  // check that all edge lengths are as expected
  for(const auto& seg : {QSeg {tet[0], tet[1]},
                         QSeg {tet[0], tet[2]},
                         QSeg {tet[0], tet[3]},
                         QSeg {tet[1], tet[2]},
                         QSeg {tet[1], tet[3]},
                         QSeg {tet[2], tet[3]}})
  {
    const double edgeLength = seg.length();
    EXPECT_NEAR(exp_edge_len, edgeLength, this->EPS);
  }

  // check that signed volume is as expected
  EXPECT_NEAR(exp_vol, tet.signedVolume(), this->EPS);

  // check that all face orientations are as expected
  for(int i = 0; i < 4; ++i)
  {
    const auto pt = tet[i];
    QTri tri;
    // clang-format off
    switch(i)
    {
    case 0: tri = QTri {tet[1], tet[2], tet[3]}; break;
    case 1: tri = QTri {tet[0], tet[3], tet[2]}; break; // note: swap since odd
    case 2: tri = QTri {tet[0], tet[1], tet[3]}; break;
    case 3: tri = QTri {tet[0], tet[2], tet[1]}; break; // note: swap since odd
    }
    // clang-format on
    EXPECT_EQ(primal::ON_NEGATIVE_SIDE, primal::orientation(pt, tri));

    // check that orientation changes if we permute a pair of triangle vertices
    axom::utilities::swap(tri[0], tri[1]);
    EXPECT_EQ(primal::ON_POSITIVE_SIDE, primal::orientation(pt, tri));
  }
}

TEST_F(TetrahedronTest, checkAndFixOrientation)
{
  using QTet = TetrahedronTest::QTet;

  int indices[] = {0, 1, 2, 3};

  for(int i = 0; i < this->numTetrahedra(); ++i)
  {
    QTet tet = this->getTet(i);
    double expVolume = tet.signedVolume();

    // Run sign check through all vertex permutations for the tetrahedron
    do
    {
      QTet tetPermuted = QTet(tet[indices[0]], tet[indices[1]], tet[indices[2]], tet[indices[3]]);

      double preCheckAbsoluteVolume = tetPermuted.volume();

      tetPermuted.checkAndFixOrientation();

      double postCheckAbsoluteVolume = tetPermuted.volume();

      EXPECT_NEAR(expVolume, postCheckAbsoluteVolume, this->EPS);

      // Verify absolute value of volume is still the same
      EXPECT_NEAR(preCheckAbsoluteVolume, postCheckAbsoluteVolume, this->EPS);

    } while(std::next_permutation(indices, indices + QTet::NUM_VERTS));
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}
