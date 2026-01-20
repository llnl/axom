// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file winding_number_3d_memoization.hpp
 *
 * \brief Consists of data structures that accelerate GWN queries through "memoization," i.e.
 *  dynamically caching and reusing patch surface evaluations and tangents at quadrature points.
 */

#ifndef AXOM_PRIMAL_WINDING_NUMBER_3D_MEMOIZATION_HPP_
#define AXOM_PRIMAL_WINDING_NUMBER_3D_MEMOIZATION_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/KnotVector.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

#include "axom/primal/operators/is_convex.hpp"

#include "axom/core/numerics/quadrature.hpp"

#include <vector>
#include <ostream>
#include <math.h>

#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
namespace detail
{

template <typename T>
class NURBSPatchGWNCache;

/*!
 * \struct TrimmingCurveQuadratureData
 *
 * \brief Stores quadrature points and tangents for a trimming curve on a patch
 * for a Gaussian quadrature rule of a given number of nodes
 */
template <typename T>
struct TrimmingCurveQuadratureData
{
  TrimmingCurveQuadratureData() = default;

  /*!
   * \brief Constructor for quadrature data from a single trimming curve on a patch
   * 
   * \param [in] quad_npts The number of Gaussian nodes
   * \param [in] a_patch The 3D NURBS surface
   * \param [in] a_curve_index The index of the trimming curve
   * \param [in] a_refinementLevel How many subdivisions for the curve
   * \param [in] a_refinementSection Which subdivision for a given level
   */
  TrimmingCurveQuadratureData(const NURBSPatch<T, 3>& a_patch,
                              int a_curve_index,
                              int quad_npts,
                              int a_refinementLevel,
                              int a_refinementSection)
    : m_quad_npts(quad_npts)
  {
    // Generate the (cached) quadrature rules in parameter space
    const numerics::QuadratureRule& gl_rule = numerics::get_gauss_legendre(quad_npts);

    auto& the_curve = a_patch.getTrimmingCurve(a_curve_index);

    const T curve_min_knot = the_curve.getMinKnot();
    const T curve_max_knot = the_curve.getMaxKnot();

    // Find the right knot span based on the refinement level
    m_span_length = (curve_max_knot - curve_min_knot) / std::pow(2, a_refinementLevel);
    const T span_offset = m_span_length * a_refinementSection;

    m_quadrature_points.resize(m_quad_npts);
    m_quadrature_tangents.resize(m_quad_npts);
    for(int q = 0; q < m_quad_npts; ++q)
    {
      const T quad_x = gl_rule.node(q) * m_span_length + curve_min_knot + span_offset;

      Point<T, 2> c_eval;
      Vector<T, 2> c_Dt;
      the_curve.evaluateFirstDerivative(quad_x, c_eval, c_Dt);

      Point<T, 3> s_eval;
      Vector<T, 3> s_Du, s_Dv;
      a_patch.evaluateFirstDerivatives(c_eval[0], c_eval[1], s_eval, s_Du, s_Dv);

      m_quadrature_points[q] = s_eval;
      m_quadrature_tangents[q] = s_Du * c_Dt[0] + s_Dv * c_Dt[1];
    }
  }

  const Point<T, 3>& getQuadraturePoint(size_t idx) const { return m_quadrature_points[idx]; }
  const Vector<T, 3>& getQuadratureTangent(size_t idx) const { return m_quadrature_tangents[idx]; }
  double getQuadratureWeight(size_t idx) const
  {
    // Because the quadrature weights are identical for each trimming curve (up to a scaling factor),
    //  we query the static rule instead of storing redundant weights
    const numerics::QuadratureRule& gl_rule = numerics::get_gauss_legendre(m_quad_npts);
    return gl_rule.weight(idx) * m_span_length;
  }
  int getNumPoints() const { return m_quad_npts; }

private:
  axom::Array<Point<T, 3>> m_quadrature_points;
  axom::Array<Vector<T, 3>> m_quadrature_tangents;
  T m_span_length;
  int m_quad_npts;
};

/*!
 * \class NURBSPatchGWNCache
 *
 * \brief Represents a NURBS patch and associated data for GWN evaluation
 * \tparam T the coordinate type, e.g., double, float, etc.
 *
 * Stores an array of maps that associates subdivisions of each trimming
 * curve with quadrature data, i.e., nodes and surface tangents.
 * 
 * Once the cache is initialized, the patch and its trimming curves are const
 * 
 * \pre Assumes a 3D NURBS patch
 */
template <typename T>
class NURBSPatchGWNCache
{
public:
  NURBSPatchGWNCache() = default;

  /// \brief Initialize the cache with the data for a single NURBS patch
  NURBSPatchGWNCache(const NURBSPatch<T, 3>& a_patch) : m_alteredPatch(a_patch)
  {
    m_alteredPatch.normalizeBySpan();

    // Calculate the average normal for the untrimmed patch
    if(!m_alteredPatch.isTrimmed())
    {
      m_averageNormal = m_alteredPatch.calculateUntrimmedPatchNormal();
      m_alteredPatch.makeTriviallyTrimmed();
    }
    else
    {
      m_averageNormal = m_alteredPatch.calculateTrimmedPatchNormal();
    }

    m_pboxDiag = m_alteredPatch.getParameterSpaceDiagonal();

    // Make a bounding box by doing (trimmed) bezier extraction,
    //  splitting the resulting bezier patches in 4,
    //  and taking a union of those bounding boxes
    const auto split_patches = m_alteredPatch.extractTrimmedBezier();

    // Bounding boxes should be defined according to the *pre-expanded* surface,
    //  since the expanded portions are never visible
    m_oBox = m_alteredPatch.orientedBoundingBox();
    m_bBox.clear();
    for(int n = 0; n < split_patches.size(); ++n)
    {
      if(split_patches[n].getNumTrimmingCurves() == 0)
      {
        continue;  // Skip patches with no trimming curves
      }

      const auto the_patch = m_alteredPatch.isRational()
        ? BezierPatch<T, 3>(split_patches[n].getControlPoints(),
                            split_patches[n].getWeights(),
                            split_patches[n].getDegree_u(),
                            split_patches[n].getDegree_v())
        : BezierPatch<T, 3>(split_patches[n].getControlPoints(),
                            split_patches[n].getDegree_u(),
                            split_patches[n].getDegree_v());

      BezierPatch<T, 3> p1, p2, p3, p4;
      the_patch.split(0.5, 0.5, p1, p2, p3, p4);

      m_bBox.addBox(p1.boundingBox());
      m_bBox.addBox(p2.boundingBox());
      m_bBox.addBox(p3.boundingBox());
      m_bBox.addBox(p4.boundingBox());
    }

    m_alteredPatch.expandParameterSpace(0.05, 0.05);

    m_curveQuadratureMaps.resize(m_alteredPatch.getNumTrimmingCurves());
  }

  /// \brief Initialize the cache with the data for a single Bezier patch
  NURBSPatchGWNCache(const BezierPatch<T, 3>& a_patch)
    : NURBSPatchGWNCache(NURBSPatch<T, 3>(a_patch))
  { }

  ///@{
  //! \name Functions that mirror functionality of NURBSPatch so signatures match in GWN evaluation.
  //!
  //! By limiting access to these functions, we ensure memoized information is always accurate
  decltype(auto) getControlPoints() const { return m_alteredPatch.getControlPoints(); }
  int getNumControlPoints_u() const { return m_alteredPatch.getNumControlPoints_u(); }
  int getNumControlPoints_v() const { return m_alteredPatch.getNumControlPoints_v(); }
  decltype(auto) getWeights() const { return m_alteredPatch.getWeights(); }
  decltype(auto) getKnots_u() const { return m_alteredPatch.getKnots_u(); }
  decltype(auto) getKnots_v() const { return m_alteredPatch.getKnots_v(); }
  double getMinKnot_u() const { return m_alteredPatch.getMinKnot_u(); }
  double getMaxKnot_u() const { return m_alteredPatch.getMaxKnot_u(); }
  double getMinKnot_v() const { return m_alteredPatch.getMinKnot_v(); }
  double getMaxKnot_v() const { return m_alteredPatch.getMaxKnot_v(); }
  decltype(auto) getTrimmingCurves() const { return m_alteredPatch.getTrimmingCurves(); };
  int getNumTrimmingCurves() const { return m_alteredPatch.getNumTrimmingCurves(); }
  decltype(auto) getParameterSpaceDiagonal() const { return m_pboxDiag; }
  //@}

  ///@{
  //! \name Accessors for precomputed data
  const Vector<T, 3>& getAverageNormal() const { return m_averageNormal; }
  const BoundingBox<T, 3>& boundingBox() const { return m_bBox; }
  const OrientedBoundingBox<T, 3>& orientedBoundingBox() const { return m_oBox; }
  //@}

  /// \brief Creates or accesses the quadrature nodes for a given trimming curve
  TrimmingCurveQuadratureData<T>& getTrimmingCurveQuadratureData(int curveIndex,
                                                                 int quadNPts,
                                                                 int refinementLevel,
                                                                 int refinementIndex) const
  {
    // Check to see if we have already computed the quadrature data for this curve
    const auto hash_key = std::make_pair(refinementLevel, refinementIndex);

    if(m_curveQuadratureMaps[curveIndex].find(hash_key) == m_curveQuadratureMaps[curveIndex].end())
    {
      m_curveQuadratureMaps[curveIndex][hash_key] = TrimmingCurveQuadratureData<T>(m_alteredPatch,
                                                                                   curveIndex,
                                                                                   quadNPts,
                                                                                   refinementLevel,
                                                                                   refinementIndex);
    }

    return m_curveQuadratureMaps[curveIndex][hash_key];
  }

private:
  // The patch is private to prevent dirtying the cache by changing the patch,
  //  and because the stored internal patch is altered from the original input
  NURBSPatch<T, 3> m_alteredPatch;

  // Per patch data
  BoundingBox<T, 3> m_bBox;
  OrientedBoundingBox<T, 3> m_oBox;
  Vector<T, 3> m_averageNormal;
  double m_pboxDiag;

  // Per trimming curve data, keyed by (whichRefinementLevel, whichRefinementIndex)
  mutable axom::Array<std::map<std::pair<int, int>, TrimmingCurveQuadratureData<T>>> m_curveQuadratureMaps;
};

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_3D_MEMOIZATION_HPP_
