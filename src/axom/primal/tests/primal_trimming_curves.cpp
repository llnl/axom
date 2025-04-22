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
    auto parameter_pt_2 = nPatchCopy.getTrimmingCurve(0).evaluate(t);

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
    auto parameter_pt_2 = nPatchCopy.getTrimmingCurve(0).evaluate(t);

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

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
