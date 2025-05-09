// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/spin/RectangularLattice.hpp"

#include "axom/slic.hpp"

// Define some helpful typedefs for 1D rectangular lattices
namespace lattice_1D
{
const int DIM = 1;

using LatticeT = axom::spin::RectangularLattice<DIM>;

using GridCell = LatticeT::GridCell;
using SpacePt = LatticeT::SpacePoint;
using SpaceVector = LatticeT::SpaceVector;
using BBox = LatticeT::SpatialBoundingBox;
using IntArray = axom::NumericArray<int, DIM>;
}  // namespace lattice_1D

// Define some helpful typedefs for 2D rectangular lattices
namespace lattice_2D
{
const int DIM = 2;

using LatticeT = axom::spin::RectangularLattice<DIM>;

using GridCell = LatticeT::GridCell;
using SpacePt = LatticeT::SpacePoint;
using SpaceVector = LatticeT::SpaceVector;
using BBox = LatticeT::SpatialBoundingBox;
using IntArray = axom::NumericArray<int, DIM>;
}  // namespace lattice_2D

// Define some helpful typedefs for 3D rectangular lattices
namespace lattice_3D
{
const int DIM = 3;

using LatticeT = axom::spin::RectangularLattice<DIM>;

using GridCell = LatticeT::GridCell;
using SpacePt = LatticeT::SpacePoint;
using SpaceVector = LatticeT::SpaceVector;
using BBox = LatticeT::SpatialBoundingBox;
using IntArray = axom::NumericArray<int, DIM>;
}  // namespace lattice_3D

TEST(spin_rectangle_lattice, lattice_ctor)
{
  SLIC_INFO("Testing lattice constructors in 1D, 2D and 3D");
  //1D
  {
    using namespace lattice_1D;
    SpacePt origin(1.1);
    SpaceVector spacing(SpacePt(.1));

    LatticeT defaultLattice;
    EXPECT_EQ(SpacePt::zero(), defaultLattice.origin());
    EXPECT_EQ(SpaceVector(SpacePt(1)), defaultLattice.spacing());

    LatticeT lattice(origin, spacing);
    EXPECT_EQ(origin, lattice.origin());
    EXPECT_EQ(spacing, lattice.spacing());
  }

  //2D
  {
    using namespace lattice_2D;
    SpacePt origin(1.1);
    SpaceVector spacing(SpacePt(.1));

    LatticeT defaultLattice;
    EXPECT_EQ(SpacePt::zero(), defaultLattice.origin());
    EXPECT_EQ(SpaceVector(SpacePt(1)), defaultLattice.spacing());

    LatticeT lattice(origin, spacing);
    EXPECT_EQ(origin, lattice.origin());
    EXPECT_EQ(spacing, lattice.spacing());
  }

  //3D
  {
    using namespace lattice_3D;
    SpacePt origin(1.1);
    SpaceVector spacing(SpacePt(.1));

    LatticeT defaultLattice;
    EXPECT_EQ(SpacePt::zero(), defaultLattice.origin());
    EXPECT_EQ(SpaceVector(SpacePt(1)), defaultLattice.spacing());

    LatticeT lattice(origin, spacing);
    EXPECT_EQ(origin, lattice.origin());
    EXPECT_EQ(spacing, lattice.spacing());
  }
}

TEST(spin_rectangle_lattice, lattice_ctor_degenerate_spacing)
{
  // The test checks that really small values for the grid spacing snap to zero
  // Here, really small is defined as EPS = 1E-50

  //2D, with small, but non-degenerate spacing in coordinate 1
  //This spacing value is retained.
  {
    using namespace lattice_2D;

    constexpr double SMALL_EPS = 1E-40;  // small, but above the threshold

    SpacePt origin(1.1);
    SpaceVector spacing(SpacePt::make_point(.1, SMALL_EPS));

    SpaceVector expSpacing;
    expSpacing[0] = .1;
    expSpacing[1] = SMALL_EPS;

    // Test from Point constructor
    {
      LatticeT lattice(origin, spacing);
      EXPECT_EQ(origin, lattice.origin());
      EXPECT_EQ(expSpacing, lattice.spacing());
    }

    // Test from pointer constructor
    {
      LatticeT lattice(origin.data(), spacing.data());
      EXPECT_EQ(origin, lattice.origin());
      EXPECT_EQ(expSpacing, lattice.spacing());
    }
  }

  //2D, with degenerate spacing in coordinate 1
  //This spacing value is snapped to zero
  {
    using namespace lattice_2D;

    constexpr double SMALL_EPS = 1E-100;  // small, and below the threshold

    SpacePt origin(1.1);
    SpaceVector spacing(SpacePt::make_point(.1, SMALL_EPS));

    SpaceVector expSpacing;
    expSpacing[0] = .1;
    expSpacing[1] = 0.;

    // Test from Point constructor
    {
      LatticeT lattice(origin, spacing);
      EXPECT_EQ(origin, lattice.origin());
      EXPECT_EQ(expSpacing, lattice.spacing());
    }

    // Test from pointer constructor
    {
      LatticeT lattice(origin.data(), spacing.data());
      EXPECT_EQ(origin, lattice.origin());
      EXPECT_EQ(expSpacing, lattice.spacing());
    }
  }
}

TEST(spin_rectangle_lattice, lattice_array_ctor)
{
  SLIC_INFO("Testing lattice constructors from arrays in 1D, 2D and 3D");

  //1D
  {
    using namespace lattice_1D;
    SpacePt origin(1.1);
    SpaceVector spacing(SpacePt(.1));

    LatticeT emptyLattice(nullptr, nullptr);
    EXPECT_EQ(SpacePt::zero(), emptyLattice.origin());
    EXPECT_EQ(SpaceVector(1.), emptyLattice.spacing());

    LatticeT latticeFromOrigin(origin.data(), nullptr);
    EXPECT_EQ(origin, latticeFromOrigin.origin());
    EXPECT_EQ(SpaceVector(1.), latticeFromOrigin.spacing());

    LatticeT latticeFromSpacing(nullptr, spacing.data());
    EXPECT_EQ(SpacePt::zero(), latticeFromSpacing.origin());
    EXPECT_EQ(spacing, latticeFromSpacing.spacing());

    LatticeT latticeFromArrays(origin.data(), spacing.data());
    EXPECT_EQ(origin, latticeFromArrays.origin());
    EXPECT_EQ(spacing, latticeFromArrays.spacing());
  }

  //2D
  {
    using namespace lattice_2D;
    SpacePt origin(1.1);
    SpaceVector spacing(SpacePt(.1));

    LatticeT emptyLattice(nullptr, nullptr);
    EXPECT_EQ(SpacePt::zero(), emptyLattice.origin());
    EXPECT_EQ(SpaceVector(1.), emptyLattice.spacing());

    LatticeT latticeFromOrigin(origin.data(), nullptr);
    EXPECT_EQ(origin, latticeFromOrigin.origin());
    EXPECT_EQ(SpaceVector(1.), latticeFromOrigin.spacing());

    LatticeT latticeFromSpacing(nullptr, spacing.data());
    EXPECT_EQ(SpacePt::zero(), latticeFromSpacing.origin());
    EXPECT_EQ(spacing, latticeFromSpacing.spacing());

    LatticeT latticeFromArrays(origin.data(), spacing.data());
    EXPECT_EQ(origin, latticeFromArrays.origin());
    EXPECT_EQ(spacing, latticeFromArrays.spacing());
  }

  //3D
  {
    using namespace lattice_3D;
    SpacePt origin(1.1);
    SpaceVector spacing(SpacePt(.1));

    LatticeT emptyLattice(nullptr, nullptr);
    EXPECT_EQ(SpacePt::zero(), emptyLattice.origin());
    EXPECT_EQ(SpaceVector(1.), emptyLattice.spacing());

    LatticeT latticeFromOrigin(origin.data(), nullptr);
    EXPECT_EQ(origin, latticeFromOrigin.origin());
    EXPECT_EQ(SpaceVector(1.), latticeFromOrigin.spacing());

    LatticeT latticeFromSpacing(nullptr, spacing.data());
    EXPECT_EQ(SpacePt::zero(), latticeFromSpacing.origin());
    EXPECT_EQ(spacing, latticeFromSpacing.spacing());

    LatticeT latticeFromArrays(origin.data(), spacing.data());
    EXPECT_EQ(origin, latticeFromArrays.origin());
    EXPECT_EQ(spacing, latticeFromArrays.spacing());
  }
}

TEST(spin_rectangle_lattice, operators)
{
  SLIC_INFO("Testing free operators in 1D, 2D and 3D");
  //1D
  {
    using namespace lattice_1D;

    SpacePt origin1(1.1);
    SpacePt origin2(2.1);

    SpaceVector spacing1(SpacePt(.1));
    SpaceVector spacing2(SpacePt(.2));

    // Test identity
    LatticeT lattice1(origin1, spacing1);
    EXPECT_EQ(lattice1, lattice1);

    // Test identify on different instances
    EXPECT_EQ(lattice1, LatticeT(origin1, spacing1));
    EXPECT_EQ(LatticeT(origin1, spacing1), lattice1);

    // Test inequality on different origins
    LatticeT lattice2(origin2, spacing1);
    EXPECT_NE(lattice1, lattice2);
    EXPECT_NE(lattice2, lattice1);

    // Test inequality on different spacings
    LatticeT lattice3(origin1, spacing2);
    EXPECT_NE(lattice1, lattice3);
    EXPECT_NE(lattice3, lattice1);

    // Test inequality on different origin and spacing
    LatticeT lattice4(origin2, spacing2);
    EXPECT_NE(lattice1, lattice4);
    EXPECT_NE(lattice4, lattice1);

    SLIC_INFO("Lattices " << (lattice1 == lattice4 ? "are" : "are not") << " equal."
                          << "\n\t First lattice: " << lattice1
                          << "\n\t Second lattice: " << lattice4);
  }

  //2D
  {
    using namespace lattice_2D;

    SpacePt origin1(1.1);
    SpacePt origin2(2.2);

    SpaceVector spacing1(SpacePt(.1));
    SpaceVector spacing2(SpacePt(.2));

    // Test identity
    LatticeT lattice1(origin1, spacing1);
    EXPECT_EQ(lattice1, lattice1);

    // Test identify on different instances
    EXPECT_EQ(lattice1, LatticeT(origin1, spacing1));
    EXPECT_EQ(LatticeT(origin1, spacing1), lattice1);

    // Test inequality on different origins
    LatticeT lattice2(origin2, spacing1);
    EXPECT_NE(lattice1, lattice2);
    EXPECT_NE(lattice2, lattice1);

    // Test inequality on different spacings
    LatticeT lattice3(origin1, spacing2);
    EXPECT_NE(lattice1, lattice3);
    EXPECT_NE(lattice3, lattice1);

    // Test inequality on different origin and spacing
    LatticeT lattice4(origin2, spacing2);
    EXPECT_NE(lattice1, lattice4);
    EXPECT_NE(lattice4, lattice1);

    SLIC_INFO("Lattices " << (lattice1 != lattice4 ? "are not" : "are") << " equal."
                          << "\n\t First lattice: " << lattice1
                          << "\n\t Second lattice: " << lattice4);
  }

  //3D
  {
    using namespace lattice_3D;

    SpacePt origin1(1.1);
    SpacePt origin2(2.2);

    SpaceVector spacing1(SpacePt(.1));
    SpaceVector spacing2(SpacePt(.2));

    // Test identity
    LatticeT lattice1(origin1, spacing1);
    EXPECT_EQ(lattice1, lattice1);

    // Test identify on different instances
    EXPECT_EQ(lattice1, LatticeT(origin1, spacing1));
    EXPECT_EQ(LatticeT(origin1, spacing1), lattice1);

    // Test inequality on different origins
    LatticeT lattice2(origin2, spacing1);
    EXPECT_NE(lattice1, lattice2);
    EXPECT_NE(lattice2, lattice1);

    // Test inequality on different spacings
    LatticeT lattice3(origin1, spacing2);
    EXPECT_NE(lattice1, lattice3);
    EXPECT_NE(lattice3, lattice1);

    // Test inequality on different origin and spacing
    LatticeT lattice4(origin2, spacing2);
    EXPECT_NE(lattice1, lattice4);
    EXPECT_NE(lattice4, lattice1);

    SLIC_INFO("Lattices " << (lattice1 == lattice4 ? "are" : "are not") << " equal."
                          << "\n\t First lattice: " << lattice1
                          << "\n\t Second lattice: " << lattice4);
  }
}

TEST(spin_rectangle_lattice, from_bounding_box)
{
  SLIC_INFO("Testing lattice creation from 1D, 2D and 3D bounding boxes");

  //1D
  {
    using namespace lattice_1D;

    BBox bbox(SpacePt(1.25), SpacePt(2.5));
    IntArray res(5);
    LatticeT lattice = axom::spin::rectangular_lattice_from_bounding_box(bbox, res);

    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[0]);
    EXPECT_DOUBLE_EQ(.25, lattice.spacing()[0]);
  }

  //2D
  {
    using namespace lattice_2D;

    BBox bbox(SpacePt(1.25), SpacePt(2.5));
    IntArray res(5);
    LatticeT lattice = axom::spin::rectangular_lattice_from_bounding_box(bbox, res);

    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[0]);
    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[1]);

    EXPECT_DOUBLE_EQ(.25, lattice.spacing()[0]);
    EXPECT_DOUBLE_EQ(.25, lattice.spacing()[1]);
  }

  //3D
  {
    using namespace lattice_3D;

    SpacePt bbMin = SpacePt::make_point(1.25, 2.5, 5.);
    SpacePt bbMax = SpacePt::make_point(2.5, 5., 10.);
    BBox bbox(bbMin, bbMax);

    int resData[3] = {5, 50, 500};
    IntArray res(resData);
    LatticeT lattice = axom::spin::rectangular_lattice_from_bounding_box(bbox, res);

    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[0]);
    EXPECT_DOUBLE_EQ(2.5, lattice.origin()[1]);
    EXPECT_DOUBLE_EQ(5., lattice.origin()[2]);

    EXPECT_DOUBLE_EQ(.25, lattice.spacing()[0]);
    EXPECT_DOUBLE_EQ(.05, lattice.spacing()[1]);
    EXPECT_DOUBLE_EQ(.01, lattice.spacing()[2]);
  }

  //2D, w/ bounding box that is degenerate in coordinate 1
  {
    using namespace lattice_2D;

    constexpr double EPS = 1E-100;
    BBox bbox(SpacePt::make_point(1.25, 1.25), SpacePt::make_point(2.5, 1.25 + EPS));
    IntArray res(5);
    LatticeT lattice = axom::spin::rectangular_lattice_from_bounding_box(bbox, res);

    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[0]);
    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[1]);

    EXPECT_DOUBLE_EQ(.25, lattice.spacing()[0]);
    EXPECT_DOUBLE_EQ(0., lattice.spacing()[1]);
  }

  //2D, w/ degenerate spacing in coordinate 1
  {
    using namespace lattice_2D;

    BBox bbox(SpacePt(1.25), SpacePt(2.5));
    IntArray res;
    res[0] = 5;
    res[1] = 0;

    LatticeT lattice = axom::spin::rectangular_lattice_from_bounding_box(bbox, res);

    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[0]);
    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[1]);

    EXPECT_DOUBLE_EQ(.25, lattice.spacing()[0]);
    EXPECT_DOUBLE_EQ(0., lattice.spacing()[1]);
  }
}

TEST(spin_rectangle_lattice, convert_point_cell_1D)
{
  SLIC_INFO("Testing point conversion of 1D lattice");

  using namespace lattice_1D;

  const double SPACING = 0.1;
  const double HALF_SPACING = SPACING / 2.;
  const double EPS = 1e-10;

  SpacePt origin = SpacePt::zero();
  SpaceVector spacing = SpaceVector(SpacePt(SPACING));

  LatticeT lattice(origin, spacing);

  EXPECT_EQ(origin, lattice.origin());
  EXPECT_EQ(spacing, lattice.spacing());

  // Test that we can map points to cells and cells to points
  for(int i = -10; i <= 10; ++i)
  {
    // Point near the lower bounds of a cell
    // Note: We have to add an epsilon to guarantee we find the right cell
    SpacePt lowerPoint(SPACING * i + EPS);

    GridCell lowerCell = lattice.gridCell(lowerPoint);
    EXPECT_EQ(i, lowerCell[0]);

    SpacePt ptFromLoweCell = lattice.spacePoint(lowerCell);
    EXPECT_DOUBLE_EQ(SPACING * i, ptFromLoweCell[0]);

    // Point near the middle of a cell
    SpacePt midPoint(SPACING * i + HALF_SPACING);

    GridCell midCell = lattice.gridCell(midPoint);
    EXPECT_EQ(i, midCell[0]);

    SpacePt ptFromMidCell = lattice.spacePoint(midCell);
    EXPECT_DOUBLE_EQ(SPACING * i, ptFromMidCell[0]);

    // Point near the upper bounds of a cell
    SpacePt upperPoint(SPACING * (i + 1) - EPS);
    GridCell upperCell = lattice.gridCell(upperPoint);
    EXPECT_EQ(i, upperCell[0]);

    SpacePt ptFromUpperCell = lattice.spacePoint(upperCell);
    EXPECT_DOUBLE_EQ(SPACING * i, ptFromUpperCell[0]);
  }
}

TEST(spin_rectangle_lattice, convert_point_cell_2D)
{
  SLIC_INFO("Testing point conversion of 2D lattice");

  using namespace lattice_2D;

  SpacePt origin = SpacePt::make_point(1.1, 2.2);
  SpaceVector spacing = SpaceVector::make_vector(.1, .2);
  LatticeT lattice(origin, spacing);

  // Test a few hand-selected points in 2D
  {
    SpacePt pt = lattice.spacePoint(GridCell::zero());
    EXPECT_DOUBLE_EQ(1.1, pt[0]);
    EXPECT_DOUBLE_EQ(2.2, pt[1]);
  }

  {
    SpacePt pt = lattice.spacePoint(GridCell::make_point(1, 0));
    EXPECT_DOUBLE_EQ(1.2, pt[0]);
    EXPECT_DOUBLE_EQ(2.2, pt[1]);
  }
  {
    SpacePt pt = lattice.spacePoint(GridCell::make_point(0, 1));
    EXPECT_DOUBLE_EQ(1.1, pt[0]);
    EXPECT_DOUBLE_EQ(2.4, pt[1]);
  }
  {
    SpacePt pt = lattice.spacePoint(GridCell::make_point(1, 1));
    EXPECT_DOUBLE_EQ(1.2, pt[0]);
    EXPECT_DOUBLE_EQ(2.4, pt[1]);
  }

  {
    SpacePt pt = lattice.spacePoint(GridCell::make_point(-1, 0));
    EXPECT_DOUBLE_EQ(1.0, pt[0]);
    EXPECT_DOUBLE_EQ(2.2, pt[1]);
  }
  {
    SpacePt pt = lattice.spacePoint(GridCell::make_point(0, -1));
    EXPECT_DOUBLE_EQ(1.1, pt[0]);
    EXPECT_DOUBLE_EQ(2.0, pt[1]);
  }
  {
    SpacePt pt = lattice.spacePoint(GridCell::make_point(-1, -1));
    EXPECT_DOUBLE_EQ(1.0, pt[0]);
    EXPECT_DOUBLE_EQ(2.0, pt[1]);
  }
}

TEST(spin_rectangle_lattice, negative_spacing_2D)
{
  SLIC_INFO("Testing point conversion of 2D lattice");

  using namespace lattice_2D;

  SpacePt origin = SpacePt::make_point(1.1, 2.2);
  SpaceVector spacing = SpaceVector::make_vector(-.1, -.2);
  LatticeT lattice(origin, spacing);

  // Test a few hand-selected points in 2D
  // Note: Spacing is negative
  {
    SpacePt pt = lattice.spacePoint(GridCell::zero());
    EXPECT_DOUBLE_EQ(1.1, pt[0]);
    EXPECT_DOUBLE_EQ(2.2, pt[1]);
  }

  {
    SpacePt pt = lattice.spacePoint(GridCell::make_point(1, 1));
    EXPECT_DOUBLE_EQ(1.0, pt[0]);
    EXPECT_DOUBLE_EQ(2.0, pt[1]);
  }

  {
    SpacePt pt = lattice.spacePoint(GridCell::make_point(-1, -1));
    EXPECT_DOUBLE_EQ(1.2, pt[0]);
    EXPECT_DOUBLE_EQ(2.4, pt[1]);
  }
}

TEST(spin_rectangle_lattice, zero_spacing_2D)
{
  SLIC_INFO("Testing point conversion of 2D lattice");

  using namespace lattice_2D;

  SpacePt origin = SpacePt::make_point(1.1, 2.2);
  SpaceVector spacing = SpaceVector::make_vector(.1, 0.);
  LatticeT lattice(origin, spacing);

  // Test a few hand-selected points in 2D
  // Note: Spacing for dim 1 is zero
  //       so, gridCell[1] must always be zero
  {
    SpacePt pt = lattice.spacePoint(GridCell::zero());
    EXPECT_DOUBLE_EQ(origin[0], pt[0]);
    EXPECT_DOUBLE_EQ(origin[1], pt[1]);

    GridCell cell = lattice.gridCell(pt);
    EXPECT_EQ(GridCell::zero(), cell);
  }

  {
    SpacePt pt = lattice.spacePoint(GridCell::make_point(1, 10));
    EXPECT_DOUBLE_EQ(1.2, pt[0]);
    EXPECT_DOUBLE_EQ(origin[1], pt[1]);

    GridCell cell = lattice.gridCell(pt);
    EXPECT_EQ(GridCell::make_point(1, 0), cell);
  }

  {
    GridCell testCell = GridCell::make_point(-1, -10);
    SpacePt pt = lattice.spacePoint(testCell);
    EXPECT_DOUBLE_EQ(1.0, pt[0]);
    EXPECT_DOUBLE_EQ(origin[1], pt[1]);  // goes to origin

    GridCell cell = lattice.gridCell(pt);
    EXPECT_EQ(GridCell::make_point(-2, 0), cell);  // cell[1] == 0

    SLIC_INFO("For lattice " << lattice << "\n\t Bounding box of cell " << testCell << " is "
                             << lattice.cellBounds(testCell));

    SLIC_INFO("For lattice " << lattice << "\n\t Bounding box of cell " << cell << " is "
                             << lattice.cellBounds(cell));
  }
}

TEST(spin_rectangle_lattice, cell_bounding_box_3D)
{
  SLIC_INFO("Testing cell bounding box in 3D");

  using namespace lattice_3D;

  SpacePt origin(1);
  SpaceVector spacing(SpacePt(2));
  LatticeT lattice(origin, spacing);

  // Get cell bounding box of cell at origin
  {
    SpacePt bbMin = origin;
    SpacePt bbMax = SpacePt(3);
    BBox expCellBBox(bbMin, bbMax);

    BBox cellBBox = lattice.cellBounds(GridCell::zero());

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(expCellBBox.getMin()[i], cellBBox.getMin()[i]);
      EXPECT_DOUBLE_EQ(expCellBBox.getMax()[i], cellBBox.getMax()[i]);
    }
  }

  // Get cell bounding box of cell at lattice point (5,4,3)
  {
    GridCell cell = GridCell::make_point(5, 4, 3);

    SpacePt bbMin = SpacePt::make_point(11., 9., 7.);
    SpacePt bbMax = SpacePt::make_point(13., 11., 9.);
    BBox expCellBBox(bbMin, bbMax);

    BBox cellBBox = lattice.cellBounds(cell);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(expCellBBox.getMin()[i], cellBBox.getMin()[i]);
      EXPECT_DOUBLE_EQ(expCellBBox.getMax()[i], cellBBox.getMax()[i]);
    }
  }
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  std::srand(105);

  result = RUN_ALL_TESTS();

  return result;
}
