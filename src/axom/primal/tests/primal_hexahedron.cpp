// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Hexahedron.hpp"

#include <cmath>

namespace primal = axom::primal;

namespace
{
double signed_volume_tet_decomp(primal::Hexahedron<double, 3> hex)
{
  double retVol = 0.0;
  axom::StackArray<primal::Tetrahedron<double, 3>, 24> tets;

  hex.triangulate(tets);

  for(int i = 0; i < 24; i++)
  {
    // Get the signed volumes
    retVol += tets[i].signedVolume();
  }

  return retVol;
}

}  // namespace
/// Test fixture for testing primal::Hexahedron
class HexahedronTest : public ::testing::Test
{
public:
  static const int DIM = 3;

  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using Hex = primal::Hexahedron<CoordType, DIM>;

protected:
  virtual void SetUp()
  {
    EPS = 1e-8;

    /*
     * Define coordinates for first hexahedron:
     *
     * 7 +---------+ 6              +z
     *   |\        |\           +y
     *   |  \      |  \           <  ^
     *   | 4 + --------+ 5         \ |
     * 3 +---|-----+ 2 |            \|
     *   \   |     \   |             -----> +x
     *    \  |      \  |
     *     \ |       \ |
     *   0  +----------+ 1
     *
     */
    qData0[0] = QPoint {0, 0, 0};
    qData0[1] = QPoint {1, 0, 0};
    qData0[2] = QPoint {1, 1, 0};
    qData0[3] = QPoint {0, 1, 0};
    qData0[4] = QPoint {0, 0, 1};
    qData0[5] = QPoint {1, 0, 1};
    qData0[6] = QPoint {1, 1, 1};
    qData0[7] = QPoint {0, 1, 1};

    /*
     * Define coordinates for second hexahedron:
     *
     * 3 +---------+ 2              +y
     *   |\        |\
     *   |  \      |  \              ^
     *   | 7 + --------+ 6           |
     * 0 +---|-----+ 1 |             |
     *   \   |     \   |             -----> +x
     *    \  |      \  |              \
     *     \ |       \ |               \
     *   4  +----------+ 5              >
     *                                   +z
     */
    qData1[0] = QPoint {-1, 0, 0};
    qData1[1] = QPoint {0, 0, 0};
    qData1[2] = QPoint {0, 1, 0};
    qData1[3] = QPoint {-1, 1, 0};
    qData1[4] = QPoint {-1, 0, 1};
    qData1[5] = QPoint {0, 0, 1};
    qData1[6] = QPoint {0, 1, 1};
    qData1[7] = QPoint {-1, 1, 1};

    // Square frustum
    // (base side length 5, top side length 2, height 1)
    qData2[0] = QPoint {0, 0, 0};
    qData2[1] = QPoint {0, 0, 5};
    qData2[2] = QPoint {5, 0, 5};
    qData2[3] = QPoint {5, 0, 0};
    qData2[4] = QPoint {1.5, 1, 1.5};
    qData2[5] = QPoint {1.5, 1, 3.5};
    qData2[6] = QPoint {3.5, 1, 3.5};
    qData2[7] = QPoint {3.5, 1, 1.5};

    // Reproducer Test Case
    qData3[0] = QPoint {-70 - (5.0 / 6.0), -165.0000000000000, -238.0000000000000};
    qData3[1] = QPoint {-70 - (5.0 / 6.0), -143.0000000000000, -238.0000000000000};
    qData3[2] = QPoint {-52.0000000000000, -143.0000000000000, -238.0000000000000};
    qData3[3] = QPoint {-52.0000000000000, -165.0000000000000, -238.0000000000000};
    qData3[4] = QPoint {-70 - (5.0 / 6.0), -165.0000000000000, -221.0000000000000};
    qData3[5] = QPoint {-70 - (5.0 / 6.0), -143.0000000000000, -221.0000000000000};
    qData3[6] = QPoint {-52.0000000000000, -143.0000000000000, -221.0000000000000};
    qData3[7] = QPoint {-52.0000000000000, -165.0000000000000, -221.0000000000000};
  }

  QPoint qData0[8];
  QPoint qData1[8];
  QPoint qData2[8];
  QPoint qData3[8];
  double EPS;
};

//------------------------------------------------------------------------------
TEST_F(HexahedronTest, defaultConstructor)
{
  using QPoint = HexahedronTest::QPoint;
  using QHex = HexahedronTest::Hex;

  const QHex hex;

  // Test ostream operator
  SLIC_INFO("Empty Hexahedron coordinates: " << hex);

  // Check indirection operator
  EXPECT_EQ(QPoint::zero(), hex[0]);
  EXPECT_EQ(QPoint::zero(), hex[1]);
  EXPECT_EQ(QPoint::zero(), hex[2]);
  EXPECT_EQ(QPoint::zero(), hex[3]);
  EXPECT_EQ(QPoint::zero(), hex[4]);
  EXPECT_EQ(QPoint::zero(), hex[5]);
  EXPECT_EQ(QPoint::zero(), hex[6]);
  EXPECT_EQ(QPoint::zero(), hex[7]);
}

TEST_F(HexahedronTest, constructFromPoints)
{
  using QPoint = HexahedronTest::QPoint;
  using QHex = HexahedronTest::Hex;

  // Access the test data
  const QPoint* pt = this->qData0;

  axom::Array<QPoint> ptArray({pt[0], pt[1], pt[2], pt[3], pt[4], pt[5], pt[6], pt[7]});

  QHex hex1(pt[0], pt[1], pt[2], pt[3], pt[4], pt[5], pt[6], pt[7]);

  QHex hex2(pt);

  QHex hex3(ptArray);

  QHex hex4({pt[0], pt[1], pt[2], pt[3], pt[4], pt[5], pt[6], pt[7]});

  // Test ostream operator
  SLIC_INFO("Hexahedron 1 coordinates: " << hex1);
  SLIC_INFO("Hexahedron 2 coordinates: " << hex2);
  SLIC_INFO("Hexahedron 3 coordinates: " << hex3);
  SLIC_INFO("Hexahedron 4 coordinates: " << hex4);

  // Check indirection operator
  EXPECT_EQ(pt[0], hex1[0]);
  EXPECT_EQ(pt[1], hex1[1]);
  EXPECT_EQ(pt[2], hex1[2]);
  EXPECT_EQ(pt[3], hex1[3]);
  EXPECT_EQ(pt[4], hex1[4]);
  EXPECT_EQ(pt[5], hex1[5]);
  EXPECT_EQ(pt[6], hex1[6]);
  EXPECT_EQ(pt[7], hex1[7]);

  EXPECT_EQ(pt[0], hex2[0]);
  EXPECT_EQ(pt[1], hex2[1]);
  EXPECT_EQ(pt[2], hex2[2]);
  EXPECT_EQ(pt[3], hex2[3]);
  EXPECT_EQ(pt[4], hex2[4]);
  EXPECT_EQ(pt[5], hex2[5]);
  EXPECT_EQ(pt[6], hex2[6]);
  EXPECT_EQ(pt[7], hex2[7]);

  EXPECT_EQ(pt[0], hex3[0]);
  EXPECT_EQ(pt[1], hex3[1]);
  EXPECT_EQ(pt[2], hex3[2]);
  EXPECT_EQ(pt[3], hex3[3]);
  EXPECT_EQ(pt[4], hex3[4]);
  EXPECT_EQ(pt[5], hex3[5]);
  EXPECT_EQ(pt[6], hex3[6]);
  EXPECT_EQ(pt[7], hex3[7]);

  EXPECT_EQ(pt[0], hex4[0]);
  EXPECT_EQ(pt[1], hex4[1]);
  EXPECT_EQ(pt[2], hex4[2]);
  EXPECT_EQ(pt[3], hex4[3]);
  EXPECT_EQ(pt[4], hex4[4]);
  EXPECT_EQ(pt[5], hex4[5]);
  EXPECT_EQ(pt[6], hex4[6]);
  EXPECT_EQ(pt[7], hex4[7]);
}

TEST_F(HexahedronTest, volume)
{
  using QPoint = HexahedronTest::QPoint;
  using QHex = HexahedronTest::Hex;

  // Access the test data
  const QPoint* pt0 = this->qData0;
  const QPoint* pt1 = this->qData1;
  const QPoint* pt2 = this->qData2;
  const QPoint* pt3 = this->qData3;

  const QPoint non_planar_pt1 {0, -1, 0};
  const QPoint non_planar_pt2 {-0.5, -0.5, -0.5};
  const QPoint non_planar_pt3 {1.25, 1.25, 1.25};

  // Initialize hexahedrons
  QHex hex0(pt0[0], pt0[1], pt0[2], pt0[3], pt0[4], pt0[5], pt0[6], pt0[7]);
  QHex hex1(pt1[0], pt1[1], pt1[2], pt1[3], pt1[4], pt1[5], pt1[6], pt1[7]);
  QHex hex2(pt2[0], pt2[1], pt2[2], pt2[3], pt2[4], pt2[5], pt2[6], pt2[7]);
  QHex hex3(pt3[0], pt3[1], pt3[2], pt3[3], pt3[4], pt3[5], pt3[6], pt3[7]);

  // Hexahedron with one nonplanar side
  QHex hex4(non_planar_pt1, pt0[1], pt0[2], pt0[3], pt0[4], pt0[5], pt0[6], pt0[7]);

  // Hexahedron with three nonplanar sides
  QHex hex5(non_planar_pt2, pt0[1], pt0[2], pt0[3], pt0[4], pt0[5], pt0[6], pt0[7]);

  // Hexahedron with all nonplanar sides
  QHex hex6(non_planar_pt2, pt0[1], pt0[2], pt0[3], pt0[4], pt0[5], non_planar_pt3, pt0[7]);

  // Check volume
  EXPECT_DOUBLE_EQ(hex0.signedVolume(), 1);
  EXPECT_DOUBLE_EQ(hex1.signedVolume(), 1);
  EXPECT_DOUBLE_EQ(hex2.signedVolume(), 13);
  EXPECT_NEAR(hex3.signedVolume(), -7043.66666666, EPS);
  EXPECT_DOUBLE_EQ(hex4.signedVolume(), 1.25);
  EXPECT_DOUBLE_EQ(hex5.signedVolume(), 1.375);
  EXPECT_DOUBLE_EQ(hex6.signedVolume(), 1.5625);
  EXPECT_DOUBLE_EQ(hex0.volume(), 1);
  EXPECT_DOUBLE_EQ(hex1.volume(), 1);
  EXPECT_DOUBLE_EQ(hex2.volume(), 13);
  EXPECT_NEAR(hex3.volume(), 7043.66666666, EPS);
  EXPECT_DOUBLE_EQ(hex4.volume(), 1.25);
  EXPECT_DOUBLE_EQ(hex5.volume(), 1.375);
  EXPECT_DOUBLE_EQ(hex6.volume(), 1.5625);

  // Check hexahedron volume against 24-tetrahedron subvolumes
  EXPECT_DOUBLE_EQ(hex0.volume(), signed_volume_tet_decomp(hex0));
  EXPECT_DOUBLE_EQ(hex1.volume(), signed_volume_tet_decomp(hex1));
  EXPECT_DOUBLE_EQ(hex2.volume(), signed_volume_tet_decomp(hex2));

  // Negative volume expected
  EXPECT_DOUBLE_EQ(hex3.signedVolume(), signed_volume_tet_decomp(hex3));

  EXPECT_DOUBLE_EQ(hex4.volume(), signed_volume_tet_decomp(hex4));
  EXPECT_DOUBLE_EQ(hex5.volume(), signed_volume_tet_decomp(hex5));
  EXPECT_DOUBLE_EQ(hex6.volume(), signed_volume_tet_decomp(hex6));
}

TEST_F(HexahedronTest, equals)
{
  using QPoint = HexahedronTest::QPoint;
  using QHex = HexahedronTest::Hex;

  // Access the test data
  const QPoint* pt0 = this->qData0;
  const QPoint* pt1 = this->qData1;

  QHex hex0(pt0[0], pt0[1], pt0[2], pt0[3], pt0[4], pt0[5], pt0[6], pt0[7]);
  // Permute the points (but keep the same shape)
  QHex hex0a(pt0[4], pt0[5], pt0[6], pt0[7], pt0[0], pt0[1], pt0[2], pt0[3]);
  // Different shape
  QHex hex1(pt1[0], pt1[1], pt1[2], pt1[3], pt1[4], pt1[5], pt1[6], pt1[7]);

  // Hex should be equal to itself.
  EXPECT_TRUE(hex0.equals(hex0));
  // Hex should be equal to its points permuted.
  EXPECT_TRUE(hex0.equals(hex0a));
  // Hex should differ from an hex with different points.
  EXPECT_FALSE(hex0.equals(hex1));
  // Another hex should be equal to itself.
  EXPECT_TRUE(hex1.equals(hex1));
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
