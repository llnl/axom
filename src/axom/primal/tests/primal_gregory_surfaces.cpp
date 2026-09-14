// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file primal_gregory_surfaces.cpp
 * \brief This file tests primal's Gregory surface functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/BezierTriangle.hpp"
#include "axom/primal/geometry/GregoryPatch.hpp"
#include "axom/primal/geometry/GregoryTriangle.hpp"

#include <array>
#include <sstream>

namespace primal = axom::primal;

namespace
{
using CoordType = double;
using GregoryPatchType = primal::GregoryPatch<CoordType>;
using BezierPatchType = primal::BezierPatch<CoordType, 3>;
using GregoryTriangleType = primal::GregoryTriangle<CoordType>;
using BezierTriangleType = primal::BezierTriangle<CoordType, 3>;
using PointType = GregoryPatchType::PointType;
using VectorType = GregoryPatchType::VectorType;
using PatchCoordsVec = GregoryPatchType::CoordsVec;
using TriangleCoordsVec = GregoryTriangleType::CoordsVec;

BezierPatchType make_sample_bezier_patch()
{
  // clang-format off
  PointType bezierControlPoints[16] = {
    PointType {0, 0, 0}, PointType {0, 2,  1}, PointType {0, 4, -1}, PointType {0, 6, 0},
    PointType {2, 0, 1}, PointType {2, 2,  3}, PointType {2, 4,  2}, PointType {2, 6, 1},
    PointType {4, 0, 0}, PointType {4, 2,  2}, PointType {4, 4,  3}, PointType {4, 6, 0},
    PointType {6, 0, 0}, PointType {6, 2, -1}, PointType {6, 4,  1}, PointType {6, 6, 0}};
  // clang-format on

  return BezierPatchType(bezierControlPoints, 3, 3);
}

BezierTriangleType make_sample_bezier_triangle()
{
  BezierTriangleType bTri(4);

  for(int i = 0; i <= bTri.getOrder(); ++i)
  {
    for(int j = 0; j <= bTri.getOrder() - i; ++j)
    {
      const auto ii = static_cast<CoordType>(i);
      const auto jj = static_cast<CoordType>(j);
      bTri(i, j) = PointType {ii + 0.25 * jj, jj - 0.1 * ii, 0.5 * ii * ii - 0.25 * jj + ii * jj};
    }
  }

  return bTri;
}

void check_patch_control_points(const GregoryPatchType& patch, const PatchCoordsVec& expected)
{
  for(int i = 0; i < GregoryPatchType::NPTS; ++i)
  {
    EXPECT_EQ(patch.getControlPoints()[i], expected[i]);
  }
}

void check_triangle_control_points(const GregoryTriangleType& triangle,
                                   const TriangleCoordsVec& expected)
{
  for(int i = 0; i < GregoryTriangleType::NPTS; ++i)
  {
    EXPECT_EQ(triangle.getControlPoints()[i], expected[i]);
  }
}

}  // namespace

//------------------------------------------------------------------------------
TEST(primal_gregorypatch, array_constructors)
{
  SLIC_INFO("Testing Gregory patch point array constructors");

  auto controlPoints = GregoryPatchType(make_sample_bezier_patch()).getControlPoints();
  const auto constControlPoints = controlPoints;

  {
    SCOPED_TRACE("Testing C-array constructor");
    GregoryPatchType patch(controlPoints.data());
    check_patch_control_points(patch, controlPoints);
  }

  {
    SCOPED_TRACE("Testing const C-array constructor");
    GregoryPatchType patch(constControlPoints.data());
    check_patch_control_points(patch, controlPoints);
  }

  {
    SCOPED_TRACE("Testing StackArray constructor");
    GregoryPatchType patch(controlPoints);
    check_patch_control_points(patch, controlPoints);
  }

  {
    SCOPED_TRACE("Testing ArrayView constructor");
    axom::ArrayView<PointType> view(controlPoints.data(), GregoryPatchType::NPTS);
    GregoryPatchType patch(view);
    check_patch_control_points(patch, controlPoints);
  }

  {
    SCOPED_TRACE("Testing const ArrayView constructor");
    axom::ArrayView<const PointType> view(constControlPoints.data(), GregoryPatchType::NPTS);
    GregoryPatchType patch(view);
    check_patch_control_points(patch, controlPoints);
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorypatch, corner_vector_constructor)
{
  SLIC_INFO("Testing Gregory patch corner and vector constructor");

  PointType corners[4] = {PointType {0.0, 0.0, 0.0},
                          PointType {1.0, 0.0, 0.0},
                          PointType {1.0, 1.0, 0.0},
                          PointType {0.0, 1.0, 0.0}};

  VectorType normals[4] = {VectorType {0.0, 0.0, 1.0},
                           VectorType {0.0, 0.0, 1.0},
                           VectorType {0.0, 0.0, 1.0},
                           VectorType {0.0, 0.0, 1.0}};

  GregoryPatchType patch(axom::ArrayView<const PointType>(corners, 4),
                         axom::ArrayView<const VectorType>(normals, 4));

  for(int i = 0; i < 4; ++i)
  {
    EXPECT_EQ(patch.getCorner(i), corners[i]);
  }

  EXPECT_EQ(patch.evaluate(0.0, 0.0), corners[0]);
  EXPECT_EQ(patch.evaluate(1.0, 0.0), corners[1]);
  EXPECT_EQ(patch.evaluate(1.0, 1.0), corners[2]);
  EXPECT_EQ(patch.evaluate(0.0, 1.0), corners[3]);
}

//------------------------------------------------------------------------------
TEST(primal_gregorypatch, evaluate)
{
  SLIC_INFO("Testing Gregory patch evaluation");

  const BezierPatchType bPatch = make_sample_bezier_patch();
  const GregoryPatchType gPatch(bPatch);

  const double uValues[] = {0.0, 0.1, 0.25, 0.3, 0.5, 0.7, 0.75, 0.8, 1.0};
  const double vValues[] = {0.0, 0.1, 0.25, 0.3, 0.5, 0.7, 0.75, 0.8, 1.0};

  for(double u : uValues)
  {
    for(double v : vValues)
    {
      const PointType gregoryPoint = gPatch.evaluate(u, v);
      const PointType bezierPoint = bPatch.evaluate(u, v);
      for(int dim = 0; dim < 3; ++dim)
      {
        EXPECT_NEAR(gregoryPoint[dim], bezierPoint[dim], 1e-12);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorypatch, derivatives)
{
  SLIC_INFO("Testing Gregory patch derivative evaluation");

  const BezierPatchType bPatch = make_sample_bezier_patch();
  const GregoryPatchType gPatch(bPatch);

  const double u = 0.25;
  const double v = 0.75;

  PointType gEval, bEval;
  VectorType gDu, gDv, bDu, bDv;
  VectorType gDuDu, gDvDv, gDuDv, bDuDu, bDvDv, bDuDv;

  gPatch.evaluateSecondDerivatives(u, v, gEval, gDu, gDv, gDuDu, gDvDv, gDuDv);
  bPatch.evaluateSecondDerivatives(u, v, bEval, bDu, bDv, bDuDu, bDvDv, bDuDv);

  for(int dim = 0; dim < 3; ++dim)
  {
    EXPECT_NEAR(gEval[dim], bEval[dim], 1e-12);
    EXPECT_NEAR(gPatch.du(u, v)[dim], bDu[dim], 1e-12);
    EXPECT_NEAR(gPatch.dv(u, v)[dim], bDv[dim], 1e-12);
    EXPECT_NEAR(gPatch.dudu(u, v)[dim], bDuDu[dim], 1e-12);
    EXPECT_NEAR(gPatch.dvdv(u, v)[dim], bDvDv[dim], 1e-12);
    EXPECT_NEAR(gPatch.dudv(u, v)[dim], bDuDv[dim], 1e-12);
    EXPECT_NEAR(gDu[dim], bDu[dim], 1e-12);
    EXPECT_NEAR(gDv[dim], bDv[dim], 1e-12);
    EXPECT_NEAR(gDuDu[dim], bDuDu[dim], 1e-12);
    EXPECT_NEAR(gDvDv[dim], bDvDv[dim], 1e-12);
    EXPECT_NEAR(gDuDv[dim], bDuDv[dim], 1e-12);
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorypatch, bounding_box)
{
  SLIC_INFO("Testing Gregory patch bounding box");

  const auto controlPoints = GregoryPatchType(make_sample_bezier_patch()).getControlPoints();
  GregoryPatchType patch(controlPoints);
  const auto bbox = patch.boundingBox();
  const auto obb = patch.orientedBoundingBox();

  for(int i = 0; i < GregoryPatchType::NPTS; ++i)
  {
    EXPECT_TRUE(bbox.contains(controlPoints[i]));
    EXPECT_TRUE(obb.contains(controlPoints[i]));
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorypatch, print)
{
  SLIC_INFO("Testing Gregory patch output stream");

  GregoryPatchType patch(make_sample_bezier_patch());
  std::ostringstream oss;
  oss << patch;

  EXPECT_NE(std::string::npos, oss.str().find("GregoryPatch("));
}

//------------------------------------------------------------------------------
TEST(primal_gregorytriangle, array_constructors)
{
  SLIC_INFO("Testing Gregory triangle point array constructors");

  auto controlPoints = GregoryTriangleType(make_sample_bezier_triangle()).getControlPoints();
  const auto constControlPoints = controlPoints;

  {
    SCOPED_TRACE("Testing C-array constructor");
    GregoryTriangleType triangle(controlPoints.data());
    check_triangle_control_points(triangle, controlPoints);
  }

  {
    SCOPED_TRACE("Testing const C-array constructor");
    GregoryTriangleType triangle(constControlPoints.data());
    check_triangle_control_points(triangle, controlPoints);
  }

  {
    SCOPED_TRACE("Testing StackArray constructor");
    GregoryTriangleType triangle(controlPoints);
    check_triangle_control_points(triangle, controlPoints);
  }

  {
    SCOPED_TRACE("Testing ArrayView constructor");
    axom::ArrayView<PointType> view(controlPoints.data(), GregoryTriangleType::NPTS);
    GregoryTriangleType triangle(view);
    check_triangle_control_points(triangle, controlPoints);
  }

  {
    SCOPED_TRACE("Testing const ArrayView constructor");
    axom::ArrayView<const PointType> view(constControlPoints.data(), GregoryTriangleType::NPTS);
    GregoryTriangleType triangle(view);
    check_triangle_control_points(triangle, controlPoints);
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorytriangle, corner_vector_constructor)
{
  SLIC_INFO("Testing Gregory triangle corner and vector constructor");

  constexpr CoordType eps = 1e-12;

  const std::array<PointType, 3> corners = {
    PointType {0.0, 0.0, 0.0},
    PointType {1.0, 0.0, 0.0},
    PointType {0.0, 1.0, 0.0},
  };

  const std::array<VectorType, 3> normals = {
    VectorType {0.0, 0.0, 1.0},
    VectorType {0.0, 0.0, 1.0},
    VectorType {0.0, 0.0, 1.0},
  };

  GregoryTriangleType triangle(axom::ArrayView<const PointType>(corners.data(), 3),
                               axom::ArrayView<const VectorType>(normals.data(), 3));

  EXPECT_TRUE(triangle.evaluate(0.0, 0.0).isNearlyEqual(corners[0], eps));
  EXPECT_TRUE(triangle.evaluate(0.0, 1.0).isNearlyEqual(corners[1], eps));
  EXPECT_TRUE(triangle.evaluate(1.0, 0.0).isNearlyEqual(corners[2], eps));

  EXPECT_EQ(triangle.getBoundaryPoint(0, 0), corners[1]);
  EXPECT_EQ(triangle.getBoundaryPoint(0, 4), corners[2]);
  EXPECT_EQ(triangle.getBoundaryPoint(1, 0), corners[2]);
  EXPECT_EQ(triangle.getBoundaryPoint(1, 4), corners[0]);
  EXPECT_EQ(triangle.getBoundaryPoint(2, 0), corners[0]);
  EXPECT_EQ(triangle.getBoundaryPoint(2, 4), corners[1]);

  auto check_near = [&](const PointType& a, const PointType& b) {
    EXPECT_NEAR(a[0], b[0], eps);
    EXPECT_NEAR(a[1], b[1], eps);
    EXPECT_NEAR(a[2], b[2], eps);
  };

  const CoordType tvals[] = {0.0, 0.2, 0.6, 1.0};
  for(const CoordType t : tvals)
  {
    check_near(triangle.evaluate(0.0, t), PointType::lerp(corners[0], corners[1], t));
    check_near(triangle.evaluate(t, 0.0), PointType::lerp(corners[0], corners[2], t));
    check_near(triangle.evaluate(t, 1.0 - t), PointType::lerp(corners[1], corners[2], t));
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorytriangle, evaluate)
{
  SLIC_INFO("Testing Gregory triangle evaluation");

  const BezierTriangleType bTri = make_sample_bezier_triangle();
  const GregoryTriangleType gTri(bTri);

  const CoordType tValues[] = {0.0, 0.1, 0.25, 0.3, 0.5, 0.7, 0.75, 0.8, 1.0};
  for(CoordType u = 0.0; u <= 1.0; u += 0.25)
  {
    for(const CoordType v : tValues)
    {
      if(u + v <= 1.0)
      {
        const PointType gregoryPoint = gTri.evaluate(u, v);
        const PointType bezierPoint = bTri.evaluate(u, v);
        for(int dim = 0; dim < 3; ++dim)
        {
          EXPECT_NEAR(gregoryPoint[dim], bezierPoint[dim], 1e-12);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorytriangle, derivatives)
{
  SLIC_INFO("Testing Gregory triangle derivative evaluation");

  const BezierTriangleType bTri = make_sample_bezier_triangle();
  const GregoryTriangleType gTri(bTri);

  const CoordType u = 0.25;
  const CoordType v = 0.35;

  PointType gEval, bEval;
  VectorType gDu, gDv, bDu, bDv;
  VectorType gDuDu, gDvDv, gDuDv, bDuDu, bDvDv, bDuDv;

  gTri.evaluateSecondDerivatives(u, v, gEval, gDu, gDv, gDuDu, gDvDv, gDuDv);
  bTri.evaluateSecondDerivatives(u, v, bEval, bDu, bDv, bDuDu, bDvDv, bDuDv);

  for(int dim = 0; dim < 3; ++dim)
  {
    EXPECT_NEAR(gEval[dim], bEval[dim], 1e-12);
    EXPECT_NEAR(gTri.du(u, v)[dim], bDu[dim], 1e-12);
    EXPECT_NEAR(gTri.dv(u, v)[dim], bDv[dim], 1e-12);
    EXPECT_NEAR(gTri.dudu(u, v)[dim], bDuDu[dim], 1e-12);
    EXPECT_NEAR(gTri.dvdv(u, v)[dim], bDvDv[dim], 1e-12);
    EXPECT_NEAR(gTri.dudv(u, v)[dim], bDuDv[dim], 1e-12);
    EXPECT_NEAR(gDu[dim], bDu[dim], 1e-12);
    EXPECT_NEAR(gDv[dim], bDv[dim], 1e-12);
    EXPECT_NEAR(gDuDu[dim], bDuDu[dim], 1e-12);
    EXPECT_NEAR(gDvDv[dim], bDvDv[dim], 1e-12);
    EXPECT_NEAR(gDuDv[dim], bDuDv[dim], 1e-12);
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorytriangle, finite_difference_first_derivatives)
{
  SLIC_INFO("Testing Gregory triangle first derivatives with finite differences");

  constexpr CoordType h = 1e-7;
  constexpr CoordType tol = 5e-5;

  const std::array<PointType, 3> corners = {
    PointType {0.0, 0.0, 0.0},
    PointType {1.0, 0.0, 0.2},
    PointType {0.0, 1.0, -0.1},
  };

  const std::array<VectorType, 3> normals = {
    VectorType {0.0, 0.0, 1.0},
    VectorType {0.0, 0.0, 1.0},
    VectorType {0.0, 0.0, 1.0},
  };

  GregoryTriangleType triangle(axom::ArrayView<const PointType>(corners.data(), 3),
                               axom::ArrayView<const VectorType>(normals.data(), 3));

  const CoordType u0 = 0.21;
  const CoordType v0 = 0.27;

  PointType eval;
  VectorType Du, Dv;
  triangle.evaluateFirstDerivatives(u0, v0, eval, Du, Dv);

  const PointType fp_u = triangle.evaluate(u0 + h, v0);
  const PointType fm_u = triangle.evaluate(u0 - h, v0);
  const PointType fp_v = triangle.evaluate(u0, v0 + h);
  const PointType fm_v = triangle.evaluate(u0, v0 - h);

  const VectorType Du_fd(fm_u, fp_u);
  const VectorType Dv_fd(fm_v, fp_v);

  for(int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(Du[d], Du_fd[d] / (2 * h), tol);
    EXPECT_NEAR(Dv[d], Dv_fd[d] / (2 * h), tol);
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorytriangle, finite_difference_second_derivatives)
{
  SLIC_INFO("Testing Gregory triangle second derivatives with finite differences");

  constexpr CoordType h = 1e-5;
  constexpr CoordType tol = 1e-4;

  const std::array<PointType, 3> corners = {
    PointType {0.0, 0.0, 0.0},
    PointType {1.2, 0.1, 0.3},
    PointType {-0.1, 1.1, -0.2},
  };

  const std::array<VectorType, 3> normals = {
    VectorType {0.0, 0.0, 1.0},
    VectorType {0.2, 0.1, 1.0},
    VectorType {-0.1, 0.2, 1.0},
  };

  GregoryTriangleType triangle(axom::ArrayView<const PointType>(corners.data(), 3),
                               axom::ArrayView<const VectorType>(normals.data(), 3));

  const CoordType u0 = 0.23;
  const CoordType v0 = 0.31;

  PointType eval;
  VectorType Du, Dv, DuDu, DvDv, DuDv;
  triangle.evaluateSecondDerivatives(u0, v0, eval, Du, Dv, DuDu, DvDv, DuDv);

  PointType eval_pu, eval_mu, eval_pv, eval_mv;
  VectorType Du_pu, Dv_pu, Du_mu, Dv_mu;
  VectorType Du_pv, Dv_pv, Du_mv, Dv_mv;
  triangle.evaluateFirstDerivatives(u0 + h, v0, eval_pu, Du_pu, Dv_pu);
  triangle.evaluateFirstDerivatives(u0 - h, v0, eval_mu, Du_mu, Dv_mu);
  triangle.evaluateFirstDerivatives(u0, v0 + h, eval_pv, Du_pv, Dv_pv);
  triangle.evaluateFirstDerivatives(u0, v0 - h, eval_mv, Du_mv, Dv_mv);

  for(int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(triangle.dudu(u0, v0)[d], DuDu[d], tol);
    EXPECT_NEAR(triangle.dvdv(u0, v0)[d], DvDv[d], tol);
    EXPECT_NEAR(triangle.dudv(u0, v0)[d], DuDv[d], tol);

    EXPECT_NEAR(DuDu[d], (Du_pu[d] - Du_mu[d]) / (2 * h), tol);
    EXPECT_NEAR(DvDv[d], (Dv_pv[d] - Dv_mv[d]) / (2 * h), tol);
    EXPECT_NEAR(DuDv[d], (Du_pv[d] - Du_mv[d]) / (2 * h), tol);
    EXPECT_NEAR(DuDv[d], (Dv_pu[d] - Dv_mu[d]) / (2 * h), tol);
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorytriangle, bounding_box)
{
  SLIC_INFO("Testing Gregory triangle bounding box");

  const auto controlPoints = GregoryTriangleType(make_sample_bezier_triangle()).getControlPoints();
  GregoryTriangleType triangle(controlPoints);
  const auto bbox = triangle.boundingBox();
  const auto obb = triangle.orientedBoundingBox();

  for(int i = 0; i < GregoryTriangleType::NPTS; ++i)
  {
    EXPECT_TRUE(bbox.contains(controlPoints[i]));
    EXPECT_TRUE(obb.contains(controlPoints[i]));
  }
}

//------------------------------------------------------------------------------
TEST(primal_gregorytriangle, print)
{
  SLIC_INFO("Testing Gregory triangle output stream");

  GregoryTriangleType triangle(make_sample_bezier_triangle());
  std::ostringstream oss;
  oss << triangle;

  EXPECT_NE(std::string::npos, oss.str().find("GregoryTriangle("));
}
