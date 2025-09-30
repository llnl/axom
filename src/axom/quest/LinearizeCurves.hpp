// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_LINEARIZE_CURVES_HPP_
#define QUEST_LINEARIZE_CURVES_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"

namespace axom
{
namespace quest
{

/*!
 * \brief This class linearizes a vector of NURBSCurve objects into line segments stored in a mint mesh.
 */
class LinearizeCurves
{
public:
  using NURBSCurve = axom::primal::NURBSCurve<double, 2>;
  using CurveArrayView = axom::ArrayView<NURBSCurve>;
  using SegmentMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

public:
  LinearizeCurves() = default;

  /// Sets the threshold for welding vertices of adjacent Pieces of curves
  void setVertexWeldingThreshold(double thresh) { m_vertexWeldThreshold = thresh; }

  /*!
   * \brief Projects high-order NURBS contours onto a linear mesh using \a segmentsPerPiece 
   * linear segments per knot span of the contour
   * 
   * Knot spans are the sub-intervals within a spline
   *
   * \param[in] curves An array view of curves to linearize.
   * \param[in] mesh The mesh object that will contain the linearized line segments.
   * \param[in] segmentsPerKnotSpan The number of segments to make per knot span.
   */
  void getLinearMeshUniform(CurveArrayView curves, SegmentMesh *mesh, int segmentsPerKnotSpan) const;

  /*!
   * \brief Projects high-order NURBS contours onto a linear mesh using \a percentError 
   *        to decide when to stop refinement.
   * 
   * \param[in] curves An array view of curves to linearize.
   * \param[in] mesh The mesh object that will contain the linearized line segments.
   * \param[in] percentError A percent of error that is acceptable to stop refinement.
   */
  void getLinearMeshNonUniform(CurveArrayView curves, SegmentMesh *mesh, double percentError) const;

  /*!
   * \brief Compute the revolved volume of the curves using quadrature.
   *
   * \param[in] curves An array view of curves to linearize.
   * \param[in] transform A 4x4 matrix transform to apply to the shape before computing
   *                      the revolved volume.
   *
   * \note We compute revolved volume on the actual shapes so we can get a
   *       real revolved volume computed using the curve functions rather than
   *       relying on a linearized curve. The revolved volume is the volume
   *       enclosed by the surface of revolution when the shape is revolved
   *       about axis of revolution.
   *
   * \return The revolved volume.
   */
  double getRevolvedVolume(CurveArrayView curves, const numerics::Matrix<double> &transform) const;

protected:
  /*!
   * \brief Compute the revolved volume of a single curve using quadrature.
   *
   * \param[in] nurbs
   * \param[in] transform A 4x4 matrix transform to apply to the shape before computing
   *                      the revolved volume.
   *
   * \return The revolved volume.
   */
  double revolvedVolume(const NURBSCurve &nurbs, const numerics::Matrix<double> &transform) const;

protected:
  double m_vertexWeldThreshold {1E-9};
};

}  // end namespace quest
}  // end namespace axom

#endif
