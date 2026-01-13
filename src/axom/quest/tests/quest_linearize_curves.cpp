// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mint.hpp"
#include "axom/slic.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/quest/LinearizeCurves.hpp"

#include "gtest/gtest.h"

#include <math.h>

using SegmentMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;

//------------------------------------------------------------------------------
/*!
 * \brief Compute the total length of segments in the segment mesh.
 *
 * \param mesh The input mesh whose segment length will be computed.
 * 
 * \return The total length of segments in the mesh.
 */
template <typename ExecPolicy = axom::SEQ_EXEC>
double totalSegmentLength(const SegmentMesh *mesh)
{
  axom::ReduceSum<ExecPolicy, double> totalSegmentLength(0.);
  axom::mint::for_all_cells<ExecPolicy, axom::mint::xargs::coords>(
    mesh,
    AXOM_LAMBDA(axom::IndexType AXOM_UNUSED_PARAM(cellID),
                const axom::numerics::Matrix<double> &coordsMatrix,
                const axom::IndexType *AXOM_UNUSED_PARAM(nodes)) {
      constexpr int xdim = 0;
      constexpr int ydim = 1;
      const double dx = coordsMatrix(xdim, 1) - coordsMatrix(xdim, 0);
      const double dy = coordsMatrix(ydim, 1) - coordsMatrix(ydim, 0);
      totalSegmentLength += sqrt(dx * dx + dy * dy);
    });
  return totalSegmentLength.get();
}

//------------------------------------------------------------------------------
/*!
 * \brief Make a circle using multiple NURBSCurves.
 *
 * \param[out] curves The array that will contain the curves.
 */
template <typename NURBSCurveType>
void makeCurves(axom::Array<NURBSCurveType> &curves, bool circle = true)
{
  const double center[] = {0., 0.};
  const double radius = 1.;
  if(circle)
  {
    const double theta[] = {0., M_PI, M_PI * 3. / 2., 2. * M_PI};
    curves.resize(3);
    curves[0] =
      NURBSCurveType::make_circular_arc_nurbs(theta[0], theta[1], center[0], center[1], radius);
    curves[1] =
      NURBSCurveType::make_circular_arc_nurbs(theta[1], theta[2], center[0], center[1], radius);
    curves[2] =
      NURBSCurveType::make_circular_arc_nurbs(theta[2], theta[3], center[0], center[1], radius);
  }
  else
  {
    // NOTE: keep the angles away from the more vertical parts of the curve to reduce the number of samples needed to approximate curve.
    const double theta[] = {M_PI / 8., M_PI * 7. / 8.};
    curves.resize(1);
    curves[0] =
      NURBSCurveType::make_circular_arc_nurbs(theta[0], theta[1], center[0], center[1], radius);
  }
}

//------------------------------------------------------------------------------
TEST(quest_linearize_curves, linearize_uniform)
{
  constexpr int DIM = 2;
  using NURBSCurveType = axom::primal::NURBSCurve<double, DIM>;

  // Define nurbs curves that represent a circle
  axom::Array<NURBSCurveType> curves;
  makeCurves(curves);

  // Linearize the curves uniformly.
  const double expectedLength = 2. * M_PI;
  const int segmentsPerKnotSpan = 30;
  axom::quest::LinearizeCurves lin;
  SegmentMesh *mesh = new SegmentMesh(DIM, axom::mint::SEGMENT);
  lin.getLinearMeshUniform(curves.view(), mesh, segmentsPerKnotSpan);
  const double actualLength = totalSegmentLength(mesh);
  EXPECT_NEAR(actualLength, expectedLength, 1.e-3);
  delete mesh;
}

//------------------------------------------------------------------------------
TEST(quest_linearize_curves, linearize_nonuniform)
{
  constexpr int DIM = 2;
  using NURBSCurveType = axom::primal::NURBSCurve<double, DIM>;

  // Define a nurbs curve that represents a circular arc from pi/8 to pi*7/8.
  axom::Array<NURBSCurveType> curves;
  makeCurves(curves, false);

  // Linearize the curves non-uniformly.
  const double expectedLength = M_PI * 3. / 4.;
  const double percentError = 0.01;
  axom::quest::LinearizeCurves lin;
  SegmentMesh *mesh = new SegmentMesh(DIM, axom::mint::SEGMENT);
  lin.getLinearMeshNonUniform(curves.view(), mesh, percentError);
  const double actualLength = totalSegmentLength(mesh);
  EXPECT_NEAR(((expectedLength - actualLength) / expectedLength), percentError, percentError);
  delete mesh;
}

//------------------------------------------------------------------------------
TEST(quest_linearize_curves, revolved_volume)
{
  constexpr int DIM = 2;
  using NURBSCurveType = axom::primal::NURBSCurve<double, DIM>;

  // Define a nurbs curve that represents a circular arc from pi/8 to pi*7/8.
  axom::Array<NURBSCurveType> curves;
  makeCurves(curves, false);

  // Compute the revolved volume
  const double expectedVolume = 4.15330715158103;
  axom::quest::LinearizeCurves lin;
  const auto transform = axom::numerics::Matrix<double>::identity(4);
  const double revolvedVolume = lin.getRevolvedVolume(curves.view(), transform);
  EXPECT_NEAR(revolvedVolume, expectedVolume, 3.e-3);
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
