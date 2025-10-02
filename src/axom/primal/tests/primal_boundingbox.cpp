// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/slic.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/NumericLimits.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_bb_policy()
{
  const int DIM = 3;
  using BoundingBoxType = primal::BoundingBox<double, DIM>;
  using PointType = primal::Point<double, DIM>;
  using VectorType = primal::Vector<double, DIM>;

  BoundingBoxType* box =
    axom::allocate<BoundingBoxType>(1, axom::execution_space<ExecSpace>::allocatorID());

  axom::for_all<ExecSpace>(
    1,
    AXOM_LAMBDA(int i) {
      box[i] = BoundingBoxType(PointType::zero(), PointType::ones());
      box[i].expand(4.5);
      box[i].dimension();
      box[i].getLongestDimension();
      box[i].shift(VectorType(4.5));
      box[i].isValid();
    });

  BoundingBoxType box_host;
  axom::copy(&box_host, box, sizeof(BoundingBoxType));

  EXPECT_EQ(box_host, BoundingBoxType(PointType(0.0), PointType(10.0)));

  axom::deallocate(box);
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_default_constructor)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  QBBox bbox;
  EXPECT_FALSE(bbox.isValid()) << "Default constructed bounding box is invalid";
  EXPECT_FALSE(bbox.contains(QPoint()))
    << "Default constructed bounding box should not contain any points";
  EXPECT_FALSE(bbox.contains(QPoint(1000)))
    << "Default constructed bounding box should not contain any points";
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_ctor_from_no_points)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  const QBBox empty_bbox;

  // Test on nullptr and zero points
  {
    QBBox bbox(nullptr, 0);

    EXPECT_FALSE(bbox.isValid());
    EXPECT_EQ(bbox, empty_bbox);
  }

  // Test with valid pointer and a size of zero
  {
    QPoint pt(2.);
    QBBox bbox(&pt, 0);

    EXPECT_FALSE(bbox.isValid());
    EXPECT_EQ(bbox, empty_bbox);
  }

  // Test with valid pointer and a size less than zero
  {
    QPoint pt(2.);
    QBBox bbox(&pt, -3);

    EXPECT_FALSE(bbox.isValid());
    EXPECT_EQ(bbox, empty_bbox);
  }
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_ctor_from_singlePt)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  QPoint pt1;
  QPoint pt2(2);

  QBBox bbox1(pt1);

  EXPECT_TRUE(bbox1.isValid());
  EXPECT_TRUE(bbox1.contains(pt1));
  EXPECT_FALSE(bbox1.contains(pt2));
  EXPECT_EQ(bbox1.getMin(), bbox1.getMax())
    << "BBox only has a single point, so bb.getMin()==bb.getMax()";

  QBBox bbox2(pt2);

  EXPECT_TRUE(bbox2.isValid());
  EXPECT_TRUE(bbox2.contains(pt2));
  EXPECT_FALSE(bbox2.contains(pt1));
  EXPECT_EQ(bbox2.getMin(), bbox2.getMax())
    << "BBox only has a single point, so bb.getMin()==bb.getMax()";
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_ctor_from_twoPoints)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  QPoint pt1(1);
  QPoint pt2(3);
  QPoint midPt = QPoint::midpoint(pt1, pt2);
  QPoint outPt(5);

  QBBox bbox1(pt1, pt2);

  EXPECT_TRUE(bbox1.isValid());
  EXPECT_TRUE(bbox1.contains(pt1));
  EXPECT_TRUE(bbox1.contains(pt2));
  EXPECT_TRUE(bbox1.contains(midPt));
  EXPECT_FALSE(bbox1.contains(outPt));

  //
  SLIC_INFO("\n** Testing from pairs of points that are " << "max and min of the bounding box");

  QBBox bbox2(pt2, pt1);

  EXPECT_TRUE(bbox2.isValid());
  EXPECT_TRUE(bbox2.contains(pt1));
  EXPECT_TRUE(bbox2.contains(pt2));
  EXPECT_TRUE(bbox2.contains(midPt));
  EXPECT_FALSE(bbox2.contains(outPt));

  //
  SLIC_INFO("\n** Testing from pairs of points that are "
            << "not the smallest and largest of the bounding box");

  const int val = 10;
  QPoint pt101 = QPoint::make_point(val, 0, val);
  QPoint pt010 = QPoint::make_point(0, val, 0);
  QPoint midPt2 = QPoint::midpoint(pt101, pt010);

  QBBox bbox3(pt101, pt010);

  EXPECT_TRUE(bbox3.isValid());
  EXPECT_TRUE(bbox3.contains(pt101));
  EXPECT_TRUE(bbox3.contains(pt010));
  EXPECT_TRUE(bbox3.contains(midPt2));
  EXPECT_TRUE(bbox3.contains(QPoint()));
  EXPECT_TRUE(bbox3.contains(QPoint(val)));
  EXPECT_TRUE(bbox3.contains(QPoint::make_point(0, val, val)));
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_ctor_from_many_points)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  // test single point
  QPoint pt1(0.);
  QBBox bbox1(&pt1, 1);

  QBBox bbox2(pt1);

  EXPECT_TRUE(bbox1 == bbox2);

  // test many points

  QPoint ptArr[10];

  for(int i = 0; i < 10; i++)
  {
    ptArr[i] = QPoint(i);
  }

  QBBox bbox3(&ptArr[0], 10);

  for(int i = 0; i < 10; i++)
  {
    EXPECT_TRUE(bbox3.contains(ptArr[i]));
  }
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_addPoint)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  QPoint pt1(1);
  QPoint pt2(3);
  QPoint pt10(10);
  QPoint pt20(20);
  QPoint pt30(30);

  QBBox bbox1(pt1);
  bbox1.addPoint(pt2);

  EXPECT_TRUE(bbox1.isValid());
  EXPECT_TRUE(bbox1.contains(pt1));
  EXPECT_TRUE(bbox1.contains(pt2));
  EXPECT_TRUE(bbox1.contains(QPoint::midpoint(pt1, pt2)));

  // Testing that if we add a missing point, it is within bounds
  EXPECT_FALSE(bbox1.contains(pt10));
  bbox1.addPoint(pt10);
  EXPECT_TRUE(bbox1.contains(pt10));

  // Testing that if we add a point, then points outside the bounds remain outside
  EXPECT_FALSE(bbox1.contains(pt30));
  bbox1.addPoint(pt20);  // note: adding 20, but testing 30
  EXPECT_FALSE(bbox1.contains(pt30));
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_test_clear)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  QPoint pt1(1);
  QPoint pt2(3);

  QBBox bbox(pt1, pt2);

  EXPECT_TRUE(bbox.isValid());
  EXPECT_TRUE(bbox.contains(pt1));
  EXPECT_TRUE(bbox.contains(pt2));
  EXPECT_TRUE(bbox.contains(QPoint::midpoint(pt1, pt2)));

  bbox.clear();
  EXPECT_FALSE(bbox.isValid()) << "After clear() the bounding box should be invalid";
  EXPECT_FALSE(bbox.contains(pt1))
    << "After clear() the bounding box should not contain any points";
  EXPECT_FALSE(bbox.contains(pt2))
    << "After clear() the bounding box should not contain any points";
  EXPECT_FALSE(bbox.contains(QPoint::midpoint(pt1, pt2)))
    << "After clear() the bounding box should not contain any points";
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_copy_and_assignment)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  QPoint pt1(1);
  QPoint pt2(3);
  QPoint midPt = QPoint::midpoint(pt1, pt2);

  QBBox bbox1(pt1, pt2);
  EXPECT_TRUE(bbox1.isValid());

  QBBox bbox2(bbox1);  // copy .ctor
  EXPECT_TRUE(bbox2.isValid());

  EXPECT_EQ(bbox1.getMin(), bbox2.getMin());
  EXPECT_EQ(bbox1.getMax(), bbox2.getMax());
  EXPECT_TRUE(bbox2.contains(midPt));

  QBBox bbox3(pt2);  // some initialization that we don't care about
  EXPECT_TRUE(bbox3.isValid());
  EXPECT_FALSE(bbox3.contains(pt1));
  bbox3 = bbox1;  // assignment operation

  EXPECT_EQ(bbox1.getMin(), bbox3.getMin());
  EXPECT_EQ(bbox1.getMax(), bbox3.getMax());
  EXPECT_TRUE(bbox3.contains(midPt));

  QBBox bbox4;
  QBBox bbox5;
  EXPECT_EQ(bbox4.getMin(), bbox5.getMin());
  EXPECT_EQ(bbox4.getMax(), bbox5.getMax());
  EXPECT_NE(bbox4.getMin(), bbox1.getMin())
    << "Empty bb should not be equal to an initialized one.";
  EXPECT_NE(bbox4.getMax(), bbox1.getMax())
    << "Empty bb should not be equal to an initialized one.";
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_test_equality)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  QPoint pt1(1);
  QPoint pt2(3);
  QPoint midPt = QPoint::midpoint(pt1, pt2);

  QBBox bbox1(pt1, pt2);
  EXPECT_TRUE(bbox1.isValid());
  EXPECT_TRUE(bbox1.contains(pt1));
  EXPECT_TRUE(bbox1.contains(pt2));

  QBBox bbox2;
  EXPECT_FALSE(bbox2.isValid()) << "Default constructed bbox is invalid";
  EXPECT_FALSE(bbox1 == bbox2)
    << "Default constructed bbox should differ from valid bbox (operator==)";
  EXPECT_TRUE(bbox1 != bbox2)
    << "Default constructed bbox should differ from valid bbox (operator!=)";

  bbox2.addPoint(midPt);
  EXPECT_TRUE(bbox2.isValid()) << "bbox should be valid after adding a point";
  EXPECT_FALSE(bbox1 == bbox2);
  EXPECT_TRUE(bbox1 != bbox2);

  bbox2.addPoint(pt2);
  EXPECT_TRUE(bbox2.isValid()) << "bbox should be valid after adding a point";
  EXPECT_FALSE(bbox1 == bbox2);
  EXPECT_TRUE(bbox1 != bbox2);

  bbox2.addPoint(pt1);
  EXPECT_TRUE(bbox1 == bbox2)
    << "after adding both endpoints, bounds should now be equal (operator==)";
  EXPECT_FALSE(bbox1 != bbox2)
    << "after adding both endpoints, bounds should now be equal (operator!=)";
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_add_box)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  //
  SLIC_INFO("Testing addBox() for two simple bounding boxes");
  QBBox bbox1(QPoint(1), QPoint(3));  // first box
  EXPECT_TRUE(bbox1.isValid());
  EXPECT_TRUE(bbox1.contains(QPoint(2)));
  EXPECT_FALSE(bbox1.contains(QPoint(4)));

  QBBox bbox2(QPoint(5), QPoint(7));  // second box
  EXPECT_TRUE(bbox2.isValid());
  EXPECT_TRUE(bbox2.contains(QPoint(6)));
  EXPECT_FALSE(bbox2.contains(QPoint(4)));

  bbox1.addBox(bbox2);  // adding the boxes (valid+valid)
  EXPECT_TRUE(bbox1.isValid());
  EXPECT_TRUE(bbox1.contains(QPoint(6)))
    << "After addBox() we should now contain points in the other box";
  EXPECT_TRUE(bbox1.contains(QPoint(4)))
    << "After addBox() we should now contain points in between the two "
       "bounding boxes";
  EXPECT_FALSE(bbox1.contains(QPoint(10))) << "Points outside both ranges are still outside";
  EXPECT_TRUE(bbox1.contains(bbox2));

  //
  SLIC_INFO("Testing addBox() for differently arranged boxes");
  QBBox bbox3(QPoint::make_point(1, 3, 3), QPoint::make_point(3, 1, 1));  // first box in 2nd test
  QBBox bbox4(QPoint::make_point(7, 5, 7), QPoint::make_point(5, 7, 7));  // second box in 2nd test
  EXPECT_TRUE(bbox3.isValid());
  EXPECT_TRUE(bbox3.contains(QPoint::make_point(2, 3, 1)));
  EXPECT_FALSE(bbox3.contains(QPoint(4)));

  bbox3.addBox(bbox4);
  EXPECT_TRUE(bbox3.isValid());
  EXPECT_TRUE(bbox3.contains(QPoint::make_point(2, 3, 1)));
  EXPECT_TRUE(bbox3.contains(QPoint::make_point(6, 6, 7)));
  EXPECT_TRUE(bbox3.contains(QPoint(4)));
  EXPECT_TRUE(bbox3.contains(bbox4));

  //
  SLIC_INFO("Testing addBox() to initially empty box");
  QBBox bbox5;
  EXPECT_FALSE(bbox5.isValid());
  EXPECT_NE(bbox5, bbox1);
  bbox5.addBox(bbox1);  // invalid + valid
  EXPECT_EQ(bbox5, bbox1);
  EXPECT_TRUE(bbox5.contains(bbox1));

  QBBox bbox6, bbox7;
  EXPECT_FALSE(bbox6.isValid());
  EXPECT_FALSE(bbox7.isValid());
  bbox6.addBox(bbox7);
  EXPECT_FALSE(bbox6.isValid());  // invalid + invalid

  QBBox bbox8(bbox1);
  EXPECT_EQ(bbox8, bbox1);
  EXPECT_TRUE(bbox8.isValid());
  bbox8.addBox(bbox7);  // valid + invalid
  EXPECT_EQ(bbox8, bbox1);
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_different_coord_types)
{
  constexpr int DIM = 3;
  using PointD = primal::Point<double, DIM>;
  using BBoxD = primal::BoundingBox<double, DIM>;

  using PointI = primal::Point<int, DIM>;
  using BBoxI = primal::BoundingBox<int, DIM>;

  // checking that an integer point is in the double bounding box
  BBoxD dBox(PointD(1.), PointD(3.));
  EXPECT_TRUE(dBox.contains(PointI(2)));
  EXPECT_FALSE(dBox.contains(PointI(4)));

  // Adding an integer point and testing
  EXPECT_FALSE(dBox.contains(PointD(3.5)));
  dBox.addPoint(PointI(4));
  EXPECT_TRUE(dBox.contains(PointD(3.5)));

  BBoxI iBox(PointI(1), PointI(3));
  EXPECT_TRUE(iBox.contains(PointD(2.5)));
  EXPECT_TRUE(iBox.contains(PointD(1.)));
  EXPECT_TRUE(iBox.contains(PointD(3.)));
  EXPECT_FALSE(iBox.contains(PointD(4.)));

  // Comparisons to double should be fine.
  EXPECT_FALSE(iBox.contains(PointD(3.5)));

  // TRICKY --- 4.5 will get rounded to 4,
  //    so 3.5 will be in the box, but 4.25 will not !!
  iBox.addPoint(PointD(4.5));
  EXPECT_TRUE(iBox.contains(PointD(3.5)));
  EXPECT_FALSE(iBox.contains(PointD(4.25)));
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_contains_bb)
{
  constexpr int DIM = 3;
  constexpr double ONE_THIRD = 1. / 3.;
  constexpr double TWO_THIRDS = 2. / 3.;

  using PointD = primal::Point<double, DIM>;
  using BBoxD = primal::BoundingBox<double, DIM>;

  BBoxD unit_box;
  unit_box.addPoint(PointD(0.));
  unit_box.addPoint(PointD(1.));

  const BBoxD empty_box;
  const BBoxD empty_box2;

  // check identity
  {
    EXPECT_TRUE(unit_box.contains(unit_box));
  }

  // Check that invalid bbox does not contain a point.
  EXPECT_FALSE(empty_box.contains(PointD(0.)));

  // Check that 1 point bbox contains the point.
  BBoxD bbox3;
  bbox3.addPoint(PointD(0.));
  EXPECT_TRUE(bbox3.contains(PointD(0.)));
  EXPECT_TRUE(bbox3.contains(bbox3));

  // check contains w/ empty boxes
  {
    EXPECT_TRUE(unit_box.contains(empty_box));
    EXPECT_FALSE(empty_box.contains(unit_box));

    EXPECT_TRUE(empty_box.contains(empty_box2));
  }

  // check full overlap
  {
    BBoxD other_box;
    other_box.addPoint(PointD(ONE_THIRD));
    other_box.addPoint(PointD(TWO_THIRDS));

    EXPECT_TRUE(unit_box.contains(other_box));
    EXPECT_FALSE(other_box.contains(unit_box));
  }

  // check full overlap, one point at boundary
  {
    BBoxD other_box;
    other_box.addPoint(PointD(ONE_THIRD));
    other_box.addPoint(PointD(1.));

    EXPECT_TRUE(unit_box.contains(other_box));
    EXPECT_FALSE(other_box.contains(unit_box));
  }

  // check partial overlap
  {
    BBoxD other_box;
    other_box.addPoint(PointD(0.5));
    other_box.addPoint(PointD(1.5));

    EXPECT_FALSE(unit_box.contains(other_box));
    EXPECT_FALSE(other_box.contains(unit_box));
  }
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_expand)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  //
  SLIC_INFO("Testing bounding box inflate");

  QBBox bbox(QPoint(1), QPoint(3));
  EXPECT_TRUE(bbox.isValid());

  // Expansion by positive number 0.5
  QBBox bbox1(bbox);
  bbox1.expand(0.5);
  EXPECT_TRUE(bbox1.isValid());
  EXPECT_EQ(bbox1.getMin(), QPoint(.5));
  EXPECT_EQ(bbox1.getMax(), QPoint(3.5));

  // Expansion by negative number -0.5 is contraction
  QBBox bbox2(bbox);
  bbox2.expand(-0.5);
  EXPECT_TRUE(bbox2.isValid());
  EXPECT_EQ(bbox2.getMin(), QPoint(1.5));
  EXPECT_EQ(bbox2.getMax(), QPoint(2.5));

  // Expand by too much 5 -- we are checking that bounds fix themselves
  QBBox bbox3(bbox);
  bbox3.expand(-5);
  EXPECT_TRUE(bbox3.isValid());
  EXPECT_EQ(bbox3.getMin(), QPoint(-2));  // 3 + -5 == -2
  EXPECT_EQ(bbox3.getMax(), QPoint(6));   // 1 - -5 == 6

  // Check that expanding an invalid bbox does nothing.
  QBBox bbox4;
  bbox4.expand(5.);
  EXPECT_FALSE(bbox4.isValid());
  EXPECT_EQ(bbox4, QBBox());
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_scale)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  //
  SLIC_INFO("Testing bounding box scale");

  QBBox bbox(QPoint(1), QPoint(3));
  EXPECT_TRUE(bbox.isValid());

  // Expansion by positive number 0.5
  QBBox bbox1(bbox);
  bbox1.scale(1.5);
  EXPECT_TRUE(bbox1.isValid());
  EXPECT_EQ(bbox1.getMin(), QPoint(.5));
  EXPECT_EQ(bbox1.getMax(), QPoint(3.5));

  // scale by a number less than 1
  QBBox bbox2(bbox);
  bbox2.scale(0.5);
  EXPECT_TRUE(bbox2.isValid());
  EXPECT_EQ(bbox2.getMin(), QPoint(1.5));
  EXPECT_EQ(bbox2.getMax(), QPoint(2.5));

  // Show that scaling by zero set the bounds to its midpoint
  QBBox bbox3(bbox);
  bbox3.scale(0.);
  EXPECT_EQ(bbox3.getMin(), bbox3.getMax());
  QPoint midpoint = bbox.getCentroid();
  EXPECT_EQ(bbox3.getMin(), midpoint);

  // Show that scaling by a negative is the same as a positive value
  QBBox bbox4(bbox);
  bbox4.scale(-1.);
  EXPECT_EQ(bbox, bbox4);

  // Show that scaling an invalid bounding box has no effect.
  QBBox bbox5, bbox6;
  EXPECT_FALSE(bbox5.isValid());
  EXPECT_FALSE(bbox6.isValid());
  bbox5.scale(1.5);
  EXPECT_EQ(bbox5, bbox6);
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_shift)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;
  using QVec = primal::Vector<CoordType, DIM>;

  //
  SLIC_INFO("Testing bounding box shift");

  QBBox bbox(QPoint(1), QPoint(3));
  EXPECT_TRUE(bbox.isValid());

  // Expansion by positive number 0.5
  QBBox bbox1(bbox);
  bbox1.shift(QVec(0.5));
  EXPECT_TRUE(bbox1.isValid());
  EXPECT_EQ(bbox1.getMin(), QPoint(1.5));
  EXPECT_EQ(bbox1.getMax(), QPoint(3.5));

  // Expansion by negative number -0.5 is contraction
  QBBox bbox2(bbox);
  bbox2.shift(QVec(-0.5));
  EXPECT_TRUE(bbox2.isValid());
  EXPECT_EQ(bbox2.getMin(), QPoint(0.5));
  EXPECT_EQ(bbox2.getMax(), QPoint(2.5));

  // Show that shifting an invalid bounding box has no effect.
  QBBox bbox3, bbox4;
  EXPECT_FALSE(bbox3.isValid());
  EXPECT_FALSE(bbox4.isValid());
  bbox3.shift(QVec(10.));
  EXPECT_EQ(bbox3, bbox4);
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, highest_lowest_values)
{
  constexpr int DIM = 3;

  // Testing that our type trait for highest and lowest values
  // is doing the right thing in our CXX11 and pre-CXX11 compilers

  // Test double
  double maxD = axom::numeric_limits<double>::max();
  double minD = -maxD;
  EXPECT_EQ(maxD, axom::numeric_limits<double>::max());
  EXPECT_EQ(minD, axom::numeric_limits<double>::lowest());

  // Test float
  double maxF = axom::numeric_limits<float>::max();
  double minF = -maxF;
  EXPECT_EQ(maxF, axom::numeric_limits<float>::max());
  EXPECT_EQ(minF, axom::numeric_limits<float>::lowest());

  // Test int
  int maxI = axom::numeric_limits<int>::max();
  int minI = axom::numeric_limits<int>::min();
  EXPECT_EQ(maxI, axom::numeric_limits<int>::max());
  EXPECT_EQ(minI, axom::numeric_limits<int>::lowest());

  // Test uint
  unsigned int maxU = axom::numeric_limits<unsigned int>::max();
  unsigned int minU = axom::numeric_limits<unsigned int>::min();
  EXPECT_EQ(maxU, axom::numeric_limits<unsigned int>::max());
  EXPECT_EQ(minU, axom::numeric_limits<unsigned int>::lowest());

  // Testing that our default constructor for bounding boxes is properly setting the range.

  // Note: The bounds are intentionally in the reverse order -- this is how we ensure
  // that adding a point to an empty bounding box always updates the bounds properly

  using BBoxD = primal::BoundingBox<double, DIM>;
  EXPECT_TRUE(BBoxD().getMin()[0] > 0);
  EXPECT_TRUE(BBoxD().getMax()[0] < 0);

  using BBoxF = primal::BoundingBox<float, DIM>;
  EXPECT_TRUE(BBoxF().getMin()[0] > 0);
  EXPECT_TRUE(BBoxF().getMax()[0] < 0);

  using BBoxI = primal::BoundingBox<int, DIM>;
  EXPECT_TRUE(BBoxI().getMin()[0] > 0);
  EXPECT_TRUE(BBoxI().getMax()[0] < 0);

  using BBoxU = primal::BoundingBox<unsigned int, DIM>;
  EXPECT_TRUE(BBoxU().getMin()[0] > 0);
  EXPECT_TRUE(BBoxU().getMax()[0] == 0);
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_longest_dimension)
{
  using PointType = primal::Point<double, 2>;
  using BoxType = primal::BoundingBox<double, 2>;

  const BoxType bbox(PointType::zero(), PointType::make_point(5.0, 10.0));
  const int longest_dimension = bbox.getLongestDimension();
  EXPECT_EQ(1, longest_dimension);
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_bisect)
{
  using PointType = primal::Point<double, 2>;
  using BoxType = primal::BoundingBox<double, 2>;

  BoxType bbox(PointType::zero(), PointType::ones());

  // Bisect along x
  BoxType right_x;
  BoxType left_x;
  bbox.bisect(right_x, left_x, 0);
  EXPECT_EQ(bbox.getMin(), right_x.getMin());
  EXPECT_EQ(PointType::make_point(0.5, 1.0), right_x.getMax());
  EXPECT_EQ(PointType::make_point(0.5, 0.0), left_x.getMin());
  EXPECT_EQ(bbox.getMax(), left_x.getMax());

  // Bisect along y
  BoxType bottom;
  BoxType top;
  bbox.bisect(bottom, top, 1);
  EXPECT_EQ(bbox.getMin(), bottom.getMin());
  EXPECT_EQ(PointType::make_point(1.0, 0.5), bottom.getMax());
  EXPECT_EQ(PointType::make_point(0.0, 0.5), top.getMin());
  EXPECT_EQ(bbox.getMax(), top.getMax());

  // Bisect longest dimension.
  const int longest = -1;
  BoxType bbox2(PointType::zero(), PointType::make_point(1., 10.));
  bbox2.bisect(bottom, top, longest);
  EXPECT_EQ(bbox2.getMin(), bottom.getMin());
  EXPECT_EQ(PointType::make_point(1.0, 5.), bottom.getMax());
  EXPECT_EQ(PointType::make_point(0.0, 5.), top.getMin());
  EXPECT_EQ(bbox2.getMax(), top.getMax());
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_get_centroid)
{
  using PointType = primal::Point<double, 2>;
  using BoxType = primal::BoundingBox<double, 2>;

  BoxType bbox(PointType::zero(), PointType::ones());
  PointType centroid = bbox.getCentroid();
  EXPECT_DOUBLE_EQ(0.5, centroid[0]);
  EXPECT_DOUBLE_EQ(0.5, centroid[1]);

  // The centroid of an invalid bbox is 0,0.
  BoxType bbox2;
  EXPECT_DOUBLE_EQ(bbox2.getCentroid()[0], 0.);
  EXPECT_DOUBLE_EQ(bbox2.getCentroid()[1], 0.);
}

//------------------------------------------------------------------------------
TEST(primal_boundingBox, bb_intersect)
{
  using PointType = primal::Point<double, 2>;
  using VectorType = primal::Vector<double, 2>;
  using BoxType = primal::BoundingBox<double, 2>;

  BoxType bbox(PointType::zero(), PointType::ones());
  BoxType bbox2(PointType::zero(), PointType::ones());
  bbox2.shift(VectorType(0.5));

  BoxType bbox_intersect = bbox.intersect(bbox2);
  EXPECT_DOUBLE_EQ(0.5, bbox_intersect.getMin()[0]);
  EXPECT_DOUBLE_EQ(0.5, bbox_intersect.getMin()[1]);
  EXPECT_DOUBLE_EQ(1., bbox_intersect.getMax()[0]);
  EXPECT_DOUBLE_EQ(1., bbox_intersect.getMax()[1]);

  // Empty bbox intersected with itself returns an invalid bbox.
  BoxType bbox_invalid;
  EXPECT_EQ(bbox_invalid, bbox_invalid.intersect(bbox_invalid));

  // Bounding box intersected with invalid bounding box returns an invalid bounding box.
  EXPECT_EQ(bbox_invalid, bbox_invalid.intersect(bbox));
  EXPECT_EQ(bbox_invalid, bbox.intersect(bbox_invalid));
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(primal_boundingBox, bb_check_policies)
{
  using seq_exec = axom::SEQ_EXEC;
  check_bb_policy<seq_exec>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using omp_exec = axom::OMP_EXEC;
  check_bb_policy<omp_exec>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_bb_policy<cuda_exec>();
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_bb_policy<hip_exec>();
#endif
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
