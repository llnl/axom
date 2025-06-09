// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_trimming_curves.cpp
 * \brief This file tests primal's trimmed NURBS patch functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/NURBSPatch.hpp"

namespace primal = axom::primal;

// Test fixture to define a common untrimmed patch
class TrimmingCurveTest : public ::testing::Test
{
public:
  static const int DIM = 3;

  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  using ParameterPointType = primal::Point<CoordType, 2>;
  using TrimmingCurveType = primal::NURBSCurve<CoordType, 2>;

protected:
  virtual void SetUp()
  {
    const int degree_u = 1;
    const int degree_v = 1;

    const int npts_u = 3;
    const int npts_v = 3;

    // Construct from 1D axom::Array
    axom::Array<PointType> controlPointsArray({PointType {0.0, 0.0, 1.0},
                                               PointType {0.0, 1.0, 0.0},
                                               PointType {0.0, 2.0, 0.0},
                                               PointType {1.0, 0.0, 0.0},
                                               PointType {1.0, 1.0, -1.0},
                                               PointType {1.0, 2.0, 0.0},
                                               PointType {2.0, 0.0, 0.0},
                                               PointType {2.0, 1.0, 0.0},
                                               PointType {2.0, 2.0, 1.0}});

    nPatch = NURBSPatchType(controlPointsArray, npts_u, npts_v, degree_u, degree_v);
  }

  NURBSPatchType nPatch;
};

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, basic_operations)
{
  SLIC_INFO("Testing access of trimming curves");

  NURBSPatchType nPatch(this->nPatch);

  EXPECT_FALSE(nPatch.isTrimmed());
  EXPECT_EQ(nPatch.getNumTrimmingCurves(), 0);

  // Add simple trimming curves
  nPatch.makeTriviallyTrimmed();
  EXPECT_TRUE(nPatch.isTrimmed());
  EXPECT_EQ(nPatch.getNumTrimmingCurves(), 4);

  // Remove all trimming curves, while staying trimmed
  nPatch.clearTrimmingCurves();
  EXPECT_TRUE(nPatch.isTrimmed());
  EXPECT_EQ(nPatch.getNumTrimmingCurves(), 0);

  // Remove all trimming curves
  nPatch.makeUntrimmed();
  EXPECT_FALSE(nPatch.isTrimmed());

  // Mark as trimmed, but empty vector of trimming curves
  nPatch.markAsTrimmed();
  EXPECT_TRUE(nPatch.isTrimmed());
  EXPECT_EQ(nPatch.getNumTrimmingCurves(), 0);

  // Make untrimmed again
  nPatch.makeUntrimmed();
  EXPECT_FALSE(nPatch.isTrimmed());

  // Add a simple trimming curve from a NURBS curve
  axom::Array<ParameterPointType> trimmingCurveControlPoints {ParameterPointType {0.25, 0.25},
                                                              ParameterPointType {0.75, 0.25},
                                                              ParameterPointType {0.75, 0.75},
                                                              ParameterPointType {0.25, 0.75},
                                                              ParameterPointType {0.25, 0.25}};

  TrimmingCurveType trimmingCurve(trimmingCurveControlPoints, 2);

  // Adding trimming curves should mark patch as trimmed
  nPatch.addTrimmingCurve(trimmingCurve);
  EXPECT_TRUE(nPatch.isTrimmed());
  EXPECT_EQ(nPatch.getNumTrimmingCurves(), 1);

  // Copy the patch
  NURBSPatchType nPatchCopy(nPatch);

  // Check that the copy is the same
  EXPECT_EQ(nPatchCopy, nPatch);
  EXPECT_TRUE(nPatchCopy.isTrimmed());
  EXPECT_EQ(nPatchCopy.getNumTrimmingCurves(), 1);

  // Mark the copy as untrimmed
  nPatchCopy.makeUntrimmed();
  EXPECT_FALSE(nPatchCopy.isTrimmed());

  // Check that the two patches are not equal
  EXPECT_NE(nPatchCopy, nPatch);

  // Remove trimming curves from the original
  nPatch.clearTrimmingCurves();

  // Check that the two patches are still not equal
  EXPECT_NE(nPatchCopy, nPatch);

  // Mark the original one as untrimmed
  nPatch.makeUntrimmed();

  // Check that the two patches are now equal
  EXPECT_EQ(nPatchCopy, nPatch);
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, visibility_queries)
{
  SLIC_INFO("Testing queries for trimming curves");

  NURBSPatchType nPatch(this->nPatch);

  // Check visibility on the untrimmed patch
  EXPECT_TRUE(nPatch.isVisible(0.5, 0.5));
  EXPECT_FALSE(nPatch.isVisible(-0.5, 0.5));
  EXPECT_FALSE(nPatch.isVisible(0.5, -0.5));
  EXPECT_FALSE(nPatch.isVisible(1.5, 0.5));
  EXPECT_FALSE(nPatch.isVisible(0.5, 1.5));

  // Add a simple trimming curve from a NURBS curve
  axom::Array<ParameterPointType> trimmingCurveControlPoints {ParameterPointType {0.25, 0.25},
                                                              ParameterPointType {0.75, 0.25},
                                                              ParameterPointType {0.75, 0.75},
                                                              ParameterPointType {0.25, 0.75},
                                                              ParameterPointType {0.25, 0.25}};

  TrimmingCurveType trimmingCurve(trimmingCurveControlPoints, 2);
  nPatch.addTrimmingCurve(trimmingCurve);

  // Check that a point on the patch is visible
  EXPECT_TRUE(nPatch.isVisible(0.5, 0.5));
  EXPECT_TRUE(nPatch.isVisible(0.4, 0.5));
  EXPECT_TRUE(nPatch.isVisible(0.4, 0.6));

  // Check that a point outside the patch is invisible
  EXPECT_FALSE(nPatch.isVisible(0.0, 0.0));
  EXPECT_FALSE(nPatch.isVisible(0.0, 0.1));
  EXPECT_FALSE(nPatch.isVisible(0.9, 0.0));

  // Check that points outside the parameter space of the patch are invisible too
  EXPECT_FALSE(nPatch.isVisible(1.5, 0.5));
  EXPECT_FALSE(nPatch.isVisible(0.5, 1.5));
  EXPECT_FALSE(nPatch.isVisible(-0.5, 0.5));
  EXPECT_FALSE(nPatch.isVisible(0.5, -0.5));

  // Remove all trimming curves, but keep it trimmed,
  //  i.e. a fully invisible patch
  nPatch.clearTrimmingCurves();
  EXPECT_TRUE(nPatch.isTrimmed());

  // Check that all points are invisible
  EXPECT_FALSE(nPatch.isVisible(0.5, 0.5));
  EXPECT_FALSE(nPatch.isVisible(0.4, 0.5));
  EXPECT_FALSE(nPatch.isVisible(0.4, 0.6));

  EXPECT_FALSE(nPatch.isVisible(0.0, 0.0));
  EXPECT_FALSE(nPatch.isVisible(0.0, 0.1));
  EXPECT_FALSE(nPatch.isVisible(0.9, 0.0));
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, visibility_queries_open_curves)
{
  SLIC_INFO("Testing queries with open trimming curves");

  NURBSPatchType nPatch(this->nPatch);

  // Check visibility on the untrimmed patch
  EXPECT_TRUE(nPatch.isVisible(0.5, 0.5));
  EXPECT_FALSE(nPatch.isVisible(-0.5, 0.5));
  EXPECT_FALSE(nPatch.isVisible(0.5, -0.5));
  EXPECT_FALSE(nPatch.isVisible(1.5, 0.5));
  EXPECT_FALSE(nPatch.isVisible(0.5, 1.5));

  // Add an open trimming curve from a NURBS curve
  axom::Array<ParameterPointType> trimmingCurveControlPoints {ParameterPointType {0.30, 0.25},
                                                              ParameterPointType {0.75, 0.25},
                                                              ParameterPointType {0.75, 0.75},
                                                              ParameterPointType {0.25, 0.75},
                                                              ParameterPointType {0.25, 0.30}};

  TrimmingCurveType trimmingCurve(trimmingCurveControlPoints, 2);
  nPatch.addTrimmingCurve(trimmingCurve);

  // Check that a point on the patch is visible
  EXPECT_TRUE(nPatch.isVisible(0.5, 0.5));
  EXPECT_TRUE(nPatch.isVisible(0.4, 0.5));
  EXPECT_TRUE(nPatch.isVisible(0.4, 0.6));

  // Check that a point outside the patch is invisible
  EXPECT_FALSE(nPatch.isVisible(0.0, 0.0));
  EXPECT_FALSE(nPatch.isVisible(0.0, 0.1));
  EXPECT_FALSE(nPatch.isVisible(0.9, 0.0));
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, trimming_curve_orientation)
{
  SLIC_INFO("Testing queries for trimming curves");

  NURBSPatchType nPatch(this->nPatch);
  NURBSPatchType nPatchCopy(nPatch);

  // Add a simple trimming curve from a NURBS curve
  axom::Array<ParameterPointType> trimmingCurveControlPoints {ParameterPointType {0.25, 0.25},
                                                              ParameterPointType {0.75, 0.25},
                                                              ParameterPointType {0.75, 0.75},
                                                              ParameterPointType {0.25, 0.75},
                                                              ParameterPointType {0.25, 0.25}};

  TrimmingCurveType trimmingCurve(trimmingCurveControlPoints, 2);
  nPatch.addTrimmingCurve(trimmingCurve);

  trimmingCurve.reverseOrientation();
  nPatchCopy.addTrimmingCurve(trimmingCurve);

  // Current convention is that reversing the trimming curve doesn't change
  //  which points are visible, since doing so only negates the winding number
  //  field, and visibility is determined by an even-odd rule

  for(auto patch : {&nPatch, &nPatchCopy})
  {
    {
      EXPECT_TRUE(patch->isVisible(0.5, 0.5));
      EXPECT_TRUE(patch->isVisible(0.4, 0.5));
      EXPECT_TRUE(patch->isVisible(0.4, 0.6));

      EXPECT_FALSE(patch->isVisible(0.0, 0.0));
      EXPECT_FALSE(patch->isVisible(0.0, 0.1));
      EXPECT_FALSE(patch->isVisible(0.9, 0.0));
    }
  }

  // We also take the convention that the orientation of the trimmed surface
  //  is determined by the base surface, not the orientation of the curves

  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(0.0, 1.0, u_pts, npts);
  axom::numerics::linspace(0.0, 1.0, v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      // Note that `evaluate` returns the point on the untrimmed surface,
      //  and does not consider the patch's trimming curves
      PointType pt1 = nPatch.evaluate(u, v);
      PointType pt2 = nPatchCopy.evaluate(u, v);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, knot_vector_manipulation)
{
  SLIC_INFO("Testing position of trimming curves after changing surface orientation");

  NURBSPatchType nPatch(this->nPatch);

  // Add a simple, not symmetric trimming curve from a NURBS curve
  axom::Array<ParameterPointType> trimmingCurveControlPoints {ParameterPointType {0.20, 0.15},
                                                              ParameterPointType {0.70, 0.15},
                                                              ParameterPointType {0.70, 0.65},
                                                              ParameterPointType {0.20, 0.65},
                                                              ParameterPointType {0.20, 0.15}};
  TrimmingCurveType trimmingCurve(trimmingCurveControlPoints, 2);
  nPatch.addTrimmingCurve(trimmingCurve);

  // Copy the patch and its trimming curve
  NURBSPatchType nPatchCopy(nPatch);

  double min_u = 1.0, max_u = 2.0;
  double min_v = -1.0, max_v = 0.5;

  nPatchCopy.rescale_u(min_u, max_u);
  nPatchCopy.rescale_v(min_v, max_v);

  constexpr int npts = 11;
  double t_pts[npts];
  axom::numerics::linspace(0.0, 1.0, t_pts, npts);

  // Make sure that the trimming curves are in the right position after rescaling
  for(auto t : t_pts)
  {
    auto parameter_pt_1 = nPatch.getTrimmingCurve(0).evaluate(t);
    auto parameter_pt_2 = nPatchCopy.getTrimmingCurve(0).evaluate(t);

    auto space_pt_1 = nPatch.evaluate(parameter_pt_1[0], parameter_pt_1[1]);
    auto space_pt_2 = nPatchCopy.evaluate(parameter_pt_2[0], parameter_pt_2[1]);

    for(int N = 0; N < DIM; ++N)
    {
      EXPECT_NEAR(space_pt_1[N], space_pt_2[N], 1e-10);
    }
  }

  // Flip in the u-direction
  nPatchCopy.reverseOrientation(0);

  // Check that the trimming curve is in the right position
  for(auto t : t_pts)
  {
    auto parameter_pt_1 = nPatch.getTrimmingCurve(0).evaluate(t);
    auto parameter_pt_2 = nPatchCopy.getTrimmingCurve(0).evaluate(1.0 - t);

    auto space_pt_1 = nPatch.evaluate(parameter_pt_1[0], parameter_pt_1[1]);
    auto space_pt_2 = nPatchCopy.evaluate(parameter_pt_2[0], parameter_pt_2[1]);

    for(int N = 0; N < DIM; ++N)
    {
      EXPECT_NEAR(space_pt_1[N], space_pt_2[N], 1e-10);
    }
  }

  // Flip in the v-direction
  nPatchCopy.reverseOrientation(1);

  // Check that the trimming curve is in the right position
  for(auto t : t_pts)
  {
    auto parameter_pt_1 = nPatch.getTrimmingCurve(0).evaluate(t);
    auto parameter_pt_2 = nPatchCopy.getTrimmingCurve(0).evaluate(t);

    auto space_pt_1 = nPatch.evaluate(parameter_pt_1[0], parameter_pt_1[1]);
    auto space_pt_2 = nPatchCopy.evaluate(parameter_pt_2[0], parameter_pt_2[1]);

    for(int N = 0; N < DIM; ++N)
    {
      EXPECT_NEAR(space_pt_1[N], space_pt_2[N], 1e-10);
    }
  }

  // Copy the patch, reverse the axes
  nPatchCopy = nPatch;
  nPatchCopy.swapAxes();

  // Check that the trimming curve is in the right position
  for(auto t : t_pts)
  {
    auto parameter_pt_1 = nPatch.getTrimmingCurve(0).evaluate(t);
    auto parameter_pt_2 = nPatchCopy.getTrimmingCurve(0).evaluate(1.0 - t);

    EXPECT_NEAR(parameter_pt_1[0], parameter_pt_2[1], 1e-10);
    EXPECT_NEAR(parameter_pt_1[1], parameter_pt_2[0], 1e-10);

    auto space_pt_1 = nPatch.evaluate(parameter_pt_1[0], parameter_pt_1[1]);
    auto space_pt_2 = nPatchCopy.evaluate(parameter_pt_2[0], parameter_pt_2[1]);

    for(int N = 0; N < DIM; ++N)
    {
      EXPECT_NEAR(space_pt_1[N], space_pt_2[N], 1e-10);
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, trimming_disk_subdivision)
{
  SLIC_INFO("Testing disk subdivision of the curve");

  NURBSPatchType nPatch(this->nPatch);
  NURBSPatchType the_disk, the_rest;

  // Apply the simplest subdivision strategy
  bool clipDisk = true;
  nPatch.diskSplit(0.5, 0.5, 0.25, the_disk, the_rest, !clipDisk);

  constexpr int npts = 10;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(the_disk.getMinKnot_u(), the_disk.getMaxKnot_u(), u_pts, npts);
  axom::numerics::linspace(the_disk.getMinKnot_v(), the_disk.getMaxKnot_v(), v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      if((u - 0.5) * (u - 0.5) + (v - 0.5) * (v - 0.5) < 0.25 * 0.25)
      {
        EXPECT_TRUE(the_disk.isVisible(u, v));
        EXPECT_FALSE(the_rest.isVisible(u, v));
      }
      else
      {
        EXPECT_FALSE(the_disk.isVisible(u, v));
        EXPECT_TRUE(the_rest.isVisible(u, v));
      }
    }
  }

  // Punch out several nested disks from the surface, in order,
  //  and check for visibility in each of them
  NURBSPatchType disk1, disk2, disk3;
  the_rest = NURBSPatchType();
  nPatch.diskSplit(0.3, 0.6, 0.2, disk1, the_rest, !clipDisk);
  the_rest.diskSplit(0.6, 0.7, 0.2, disk2, the_rest, !clipDisk);
  the_rest.diskSplit(0.5, 0.4, 0.2, disk3, the_rest, !clipDisk);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      // Check points inside the first circle
      if((u - 0.3) * (u - 0.3) + (v - 0.6) * (v - 0.6) < 0.2 * 0.2)
      {
        EXPECT_FALSE(the_rest.isVisible(u, v));
        EXPECT_TRUE(disk1.isVisible(u, v));
        EXPECT_FALSE(disk2.isVisible(u, v));
        EXPECT_FALSE(disk3.isVisible(u, v));
      }
      else
      {
        // Check points inside the second circle, but not the first
        if((u - 0.6) * (u - 0.6) + (v - 0.7) * (v - 0.7) < 0.2 * 0.2)
        {
          EXPECT_FALSE(the_rest.isVisible(u, v));
          EXPECT_FALSE(disk1.isVisible(u, v));
          EXPECT_TRUE(disk2.isVisible(u, v));
          EXPECT_FALSE(disk3.isVisible(u, v));
        }
        else
        {
          // Check points inside the third, but not the first or second
          if((u - 0.5) * (u - 0.5) + (v - 0.4) * (v - 0.4) < 0.2 * 0.2)
          {
            EXPECT_FALSE(the_rest.isVisible(u, v));
            EXPECT_FALSE(disk1.isVisible(u, v));
            EXPECT_FALSE(disk2.isVisible(u, v));
            EXPECT_TRUE(disk3.isVisible(u, v));
          }
          // Check points inside none of the circles
          else
          {
            EXPECT_TRUE(the_rest.isVisible(u, v));
            EXPECT_FALSE(disk1.isVisible(u, v));
            EXPECT_FALSE(disk2.isVisible(u, v));
            EXPECT_FALSE(disk3.isVisible(u, v));
          }
        }
      }
    }
  }

  // Check that the disk is correctly clipped when clipDisk is true
  nPatch = this->nPatch;
  nPatch.diskSplit(0.5, 0.5, 0.25, the_disk, the_rest, clipDisk);

  EXPECT_NEAR(nPatch.getMinKnot_u(), the_rest.getMinKnot_u(), 1e-10);
  EXPECT_NEAR(nPatch.getMinKnot_v(), the_rest.getMinKnot_v(), 1e-10);
  EXPECT_NEAR(nPatch.getMaxKnot_u(), the_rest.getMaxKnot_u(), 1e-10);
  EXPECT_NEAR(nPatch.getMaxKnot_v(), the_rest.getMaxKnot_v(), 1e-10);

  EXPECT_NEAR(the_disk.getMinKnot_u(), 0.5 - 0.25, 1e-10);
  EXPECT_NEAR(the_disk.getMinKnot_v(), 0.5 - 0.25, 1e-10);
  EXPECT_NEAR(the_disk.getMaxKnot_u(), 0.5 + 0.25, 1e-10);
  EXPECT_NEAR(the_disk.getMaxKnot_v(), 0.5 + 0.25, 1e-10);
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, trimming_disk_subdivision_edge_cases)
{
  SLIC_INFO("Testing edge cases for disk subdivision");

  // Punch out some disks that are tangent to each other
  //  and the original patch boundaries
  NURBSPatchType nPatch(this->nPatch);
  NURBSPatchType the_rest, disks[7];

  double radius = 0.1;
  double centers[][2] = {{0.5, 0.5},
                         {0.7, 0.5},
                         {0.5, 0.3},
                         {0.4, 0.5 + 0.1 * std::sqrt(3.0)},
                         {0.1, 0.1},
                         {0.5, 0.9},
                         {0.9, 0.7}};

  bool clipDisk = true;

  // Get the first disk
  nPatch.diskSplit(centers[0][0], centers[0][1], radius, disks[0], the_rest, !clipDisk);

  // Get the rest of the disks
  for(int i = 1; i < 7; ++i)
  {
    the_rest.diskSplit(centers[i][0], centers[i][1], radius, disks[i], the_rest, !clipDisk);
  }

  constexpr int npts = 10;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(the_rest.getMinKnot_u(), the_rest.getMaxKnot_u(), u_pts, npts);

  axom::numerics::linspace(the_rest.getMinKnot_v(), the_rest.getMaxKnot_v(), v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      bool inADisk = false;
      for(int i = 0; i < 7; ++i)
      {
        // Ensure it's in the disk that was punctured out
        if((u - centers[i][0]) * (u - centers[i][0]) + (v - centers[i][1]) * (v - centers[i][1]) <
           radius * radius)
        {
          EXPECT_TRUE(disks[i].isVisible(u, v));
          EXPECT_FALSE(the_rest.isVisible(u, v));

          inADisk = true;
          break;
        }
      }

      // If it's not in any disk, then it should be in the rest of the patch
      if(!inADisk)
      {
        EXPECT_TRUE(the_rest.isVisible(u, v));
        for(int i = 0; i < 7; ++i)
        {
          EXPECT_FALSE(disks[i].isVisible(u, v));
        }
      }
    }
  }

  // Check very degenerate case, when two disks are punctured on top of each other
  nPatch = this->nPatch;
  NURBSPatchType disk1, disk2;

  nPatch.diskSplit(0.5, 0.5, radius, disk1, the_rest, !clipDisk);
  the_rest.diskSplit(0.5, 0.5, radius, disk2, the_rest, !clipDisk);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      // The second disk punctured out, disk2, should be empty
      if((u - 0.5) * (u - 0.5) + (v - 0.5) * (v - 0.5) < radius * radius)
      {
        EXPECT_TRUE(disk1.isVisible(u, v));
        EXPECT_FALSE(disk2.isVisible(u, v));
        EXPECT_FALSE(the_rest.isVisible(u, v));
      }
      else
      {
        EXPECT_FALSE(disk1.isVisible(u, v));
        EXPECT_FALSE(disk2.isVisible(u, v));
        EXPECT_TRUE(the_rest.isVisible(u, v));
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, trimming_edge_subdivision)
{
  SLIC_INFO("Testing disk subdivision of the curve");

  NURBSPatchType nPatch(this->nPatch);
  NURBSPatchType the_disk, the_rest;

  bool clipDisk = true, normalize = true;
  nPatch.diskSplit(0.65, 0.35, 0.18, the_disk, the_rest, !clipDisk);
  the_rest.diskSplit(0.35, 0.65, 0.18, the_disk, the_rest, !clipDisk);

  NURBSPatchType bottomleft, bottomright, topleft, topright;
  the_rest.split(0.5, 0.5, bottomleft, bottomright, topleft, topright, !normalize);

  constexpr int npts = 10;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(the_rest.getMinKnot_u(), the_rest.getMaxKnot_u(), u_pts, npts);

  axom::numerics::linspace(the_rest.getMinKnot_v(), the_rest.getMaxKnot_v(), v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      // Figure out if a point is inside either of the disks
      bool inDisks;
      if((u - 0.65) * (u - 0.65) + (v - 0.35) * (v - 0.35) < 0.18 * 0.18 ||
         (u - 0.35) * (u - 0.35) + (v - 0.65) * (v - 0.65) < 0.18 * 0.18)
      {
        inDisks = true;
      }
      else
      {
        inDisks = false;
      }

      // The point should only be in the patch if it's in the right parameter space,
      //  and isn't inside either disk
      EXPECT_EQ((u < 0.5) && (v < 0.5) && !inDisks, bottomleft.isVisible(u, v));
      EXPECT_EQ((u < 0.5) && (v > 0.5) && !inDisks, topleft.isVisible(u, v));
      EXPECT_EQ((u > 0.5) && (v < 0.5) && !inDisks, bottomright.isVisible(u, v));
      EXPECT_EQ((u > 0.5) && (v > 0.5) && !inDisks, topright.isVisible(u, v));
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, trimming_edge_subdivision_edge_cases)
{
  SLIC_INFO("Testing edge cases of disk subdivision of the curve");

  NURBSPatchType nPatch(this->nPatch);
  NURBSPatchType the_disk, the_rest;
  NURBSPatchType leftpatch, rightpatch;
  NURBSPatchType bottomleft, topleft;

  bool clipDisk = true, normalize = true;
  nPatch.diskSplit(0.5, 0.5, 0.2, the_disk, the_rest, !clipDisk);

  the_rest.split_u(0.7, leftpatch, rightpatch, !normalize);
  leftpatch.split_v(0.3, bottomleft, topleft, !normalize);

  constexpr int npts = 10;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(the_rest.getMinKnot_u(), the_rest.getMaxKnot_u(), u_pts, npts);
  axom::numerics::linspace(the_rest.getMinKnot_v(), the_rest.getMaxKnot_v(), v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      // Figure out if a point is inside either of the disks
      bool inDisk = (u - 0.5) * (u - 0.5) + (v - 0.5) * (v - 0.5) < 0.2 * 0.2;

      // The point should only be in the patch if it's in the right parameter space,
      //  and isn't inside the disk
      EXPECT_EQ((u < 0.7) && !inDisk, leftpatch.isVisible(u, v));
      EXPECT_EQ((u > 0.7) && !inDisk, rightpatch.isVisible(u, v));
      EXPECT_EQ((u < 0.7) && (v > 0.3) && !inDisk, topleft.isVisible(u, v));
      EXPECT_EQ((u < 0.7) && (v < 0.3) && !inDisk, bottomleft.isVisible(u, v));
    }
  }

  nPatch = this->nPatch;

  // Add trimming curves to the patch with a straight edge at u=0.5
  nPatch.addTrimmingCurve(TrimmingCurveType::make_linear_segment_nurbs({0.5, 0.75}, {0.5, 0.25}));

  nPatch.addTrimmingCurve(
    TrimmingCurveType::make_circular_arc_nurbs(-M_PI / 2.0, M_PI / 2.0, 0.5, 0.5, 0.25));

  // Split the patch along the same edge
  NURBSPatchType leftpatch2, rightpatch2;
  nPatch.split_u(0.5, leftpatch2, rightpatch2, !normalize);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      // Figure out if a point is inside either of the disks
      bool inHalfDisk = (u > 0.5) && (u - 0.5) * (u - 0.5) + (v - 0.5) * (v - 0.5) < 0.25 * 0.25;

      // The point should only be in the right patch if it's in the half disk,
      //  and should never be in the left patch
      EXPECT_EQ(inHalfDisk, rightpatch2.isVisible(u, v));
      EXPECT_FALSE(leftpatch2.isVisible(u, v));
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, trimming_edge_subdivision_tangent_edge_cases)
{
  SLIC_INFO("Testing edge cases of subdivision at tangent points to trimming curves");

  NURBSPatchType nPatch(this->nPatch);
  constexpr int npts = 10;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(nPatch.getMinKnot_u(), nPatch.getMaxKnot_u(), u_pts, npts);
  axom::numerics::linspace(nPatch.getMinKnot_v(), nPatch.getMaxKnot_v(), v_pts, npts);

  axom::Array<TrimmingCurveType> trimmingCurves;
  axom::Array<ParameterPointType> trimmingCurveControlPoints {ParameterPointType {1.0, 1.0},
                                                              ParameterPointType {0.75, 0.0},
                                                              ParameterPointType {0.25, 1.0},
                                                              ParameterPointType {0.0, 0.0}};

  TrimmingCurveType interestingCurve(trimmingCurveControlPoints, 3);

  trimmingCurves.push_back(interestingCurve);
  trimmingCurves.push_back(TrimmingCurveType::make_linear_segment_nurbs({0.0, 0.0}, {1.0, 0.0}));
  trimmingCurves.push_back(TrimmingCurveType::make_linear_segment_nurbs({1.0, 0.0}, {1.0, 1.0}));

  nPatch.addTrimmingCurves(trimmingCurves);

  bool normalize = true;
  NURBSPatchType toppatch, bottompatch;
  nPatch.split_v(0.5, bottompatch, toppatch, !normalize);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      if(v <= 0.5)
      {
        EXPECT_EQ(nPatch.isVisible(u, v), bottompatch.isVisible(u, v));
      }
      else
      {
        EXPECT_EQ(nPatch.isVisible(u, v), toppatch.isVisible(u, v));
      }
    }
  }

  // Reset the original patch
  nPatch = this->nPatch;
  nPatch.addTrimmingCurves(trimmingCurves);

  bool clipDisk = true;
  NURBSPatchType the_rest, the_disk;

  nPatch.diskSplit(0.5, -0.5, 1.0, the_disk, the_rest, !clipDisk);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      if(nPatch.isVisible(u, v))
      {
        // Figure out if a point is inside either of the disks
        bool inDisk = (u - 0.5) * (u - 0.5) + (v + 0.5) * (v + 0.5) < 1.0 * 1.0;

        EXPECT_EQ(inDisk, the_disk.isVisible(u, v));
        EXPECT_EQ(!inDisk, the_rest.isVisible(u, v));
      }
      else
      {
        EXPECT_EQ(the_disk.isVisible(u, v), false);
        EXPECT_EQ(the_rest.isVisible(u, v), false);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(TrimmingCurveTest, trimming_edge_subdivision_overlapping_edge)
{
  SLIC_INFO("Test a case where one of the trimming curves is entirely overlapping the edge");

  NURBSPatchType nPatch(this->nPatch);
  axom::Array<TrimmingCurveType> trimmingCurves;

  constexpr int npts = 10;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(nPatch.getMinKnot_u(), nPatch.getMaxKnot_u(), u_pts, npts);
  axom::numerics::linspace(nPatch.getMinKnot_v(), nPatch.getMaxKnot_v(), v_pts, npts);

  axom::Array<ParameterPointType> leftTrimmingControlPoints {ParameterPointType {1.0 / 2.0, 1.0},
                                                             ParameterPointType {0.75 / 2.0, 0.0},
                                                             ParameterPointType {0.25 / 2.0, 1.0},
                                                             ParameterPointType {0.0 / 2.0, 0.0}};

  TrimmingCurveType leftCurve(leftTrimmingControlPoints, 3), dummyCurve;
  leftCurve.split(0.5, dummyCurve, leftCurve);
  trimmingCurves.push_back(leftCurve);

  // Make a degree elevated curve between (0.75, 0.5) and (0.25, 0.5)
  TrimmingCurveType connecting_curve(5);
  for(int i = 0; i < 6; ++i)
  {
    connecting_curve[i] = ParameterPointType {0.75 - 0.5 * (i / 5.0) * (i / 5.0), 0.5};
  }
  trimmingCurves.push_back(connecting_curve);

  axom::Array<ParameterPointType> rightTrimmingControlPoints {
    ParameterPointType {0.5 + 1.0 / 2.0, 1.0},
    ParameterPointType {0.5 + 0.75 / 2.0, 0.0},
    ParameterPointType {0.5 + 0.25 / 2.0, 1.0},
    ParameterPointType {0.5 + 0.0 / 2.0, 0.0}};

  TrimmingCurveType rightCurve(rightTrimmingControlPoints, 3);
  rightCurve.split(0.5, rightCurve, dummyCurve);
  trimmingCurves.push_back(rightCurve);

  trimmingCurves.push_back(TrimmingCurveType::make_linear_segment_nurbs({1.0, 0.0}, {1.0, 1.0}));
  trimmingCurves.push_back(TrimmingCurveType::make_linear_segment_nurbs({0.0, 0.0}, {1.0, 0.0}));

  nPatch.addTrimmingCurves(trimmingCurves);

  bool normalize = true;
  NURBSPatchType toppatch, bottompatch;
  nPatch.split_v(0.5, bottompatch, toppatch, !normalize);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      if(v <= 0.5)
      {
        EXPECT_EQ(nPatch.isVisible(u, v), bottompatch.isVisible(u, v));
      }
      else
      {
        EXPECT_EQ(nPatch.isVisible(u, v), toppatch.isVisible(u, v));
      }
    }
  }
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
