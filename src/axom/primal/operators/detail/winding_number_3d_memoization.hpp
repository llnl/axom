// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file NURBSCurve.hpp
 *
 * \brief A NURBS curve primitive
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

#include <vector>
#include <ostream>
#include <math.h>

#include "axom/fmt.hpp"

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#else
  #error "3D GWN evaluation requires mfem library."
#endif

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
   * \param [in] a_curve The 2D parameter space NURBS curve on \a a_patch
   * \param [in] a_patch The 3D NURBS surface
   * \param [in] a_refinementLevel How mnay subdivisions for the curve
   * \param [in] a_refinementSection Which subdivision for a given level
   */
  TrimmingCurveQuadratureData(int quad_npts,
                              const NURBSCurve<T, 2>& a_curve,
                              const NURBSPatch<T, 3>& a_patch,
                              int a_refinementLevel,
                              int a_refinementSection)
    : m_quad_npts(quad_npts)
  {
    // Generate the (cached) quadrature rules in parameter space
    static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
    const mfem::IntegrationRule& quad_rule =
      my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * quad_npts - 1);

    const T curve_min_knot = a_curve.getMinKnot();
    const T curve_max_knot = a_curve.getMaxKnot();

    // Find the right knot span based on the refinement level
    m_span_length = (curve_max_knot - curve_min_knot) / std::pow(2, a_refinementLevel);
    const T span_offset = m_span_length * a_refinementSection;

    m_quadrature_points.resize(m_quad_npts);
    m_quadrature_tangents.resize(m_quad_npts);
    for(int q = 0; q < m_quad_npts; ++q)
    {
      T quad_x = quad_rule.IntPoint(q).x * m_span_length + curve_min_knot + span_offset;

      Point<T, 2> c_eval;
      Vector<T, 2> c_Dt;
      a_curve.evaluateFirstDerivative(quad_x, c_eval, c_Dt);

      Point<T, 3> s_eval;
      Vector<T, 3> s_Du, s_Dv;
      a_patch.evaluateFirstDerivatives(c_eval[0], c_eval[1], s_eval, s_Du, s_Dv);

      m_quadrature_points[q] = s_eval;
      m_quadrature_tangents[q] = s_Du * c_Dt[0] + s_Dv * c_Dt[1];
    }
  }

  Point<T, 3> getQuadraturePoint(size_t idx) const { return m_quadrature_points[idx]; }
  Vector<T, 3> getQuadratureTangent(size_t idx) const { return m_quadrature_tangents[idx]; }
  double getQuadratureWeight(size_t idx) const
  {
    // Is this efficient because it's cached? More or less efficient than storing a bunch of duplicated weights for each curve?
    static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
    const mfem::IntegrationRule& quad_rule =
      my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * m_quad_npts - 1);
    return quad_rule.IntPoint(idx).weight * m_span_length;
  }
  double getNumPoints() const { return m_quad_npts; }

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

    // Make a bounding box by doing bezier extraction, then splitting the resulting bezier patches in 4,
    //  and taking a union of those bounding boxes
    axom::Array<T> knot_vals_u = m_alteredPatch.getKnots_u().getUniqueKnots();
    axom::Array<T> knot_vals_v = m_alteredPatch.getKnots_v().getUniqueKnots();

    const auto num_knot_span_u = knot_vals_u.size() - 1;
    const auto num_knot_span_v = knot_vals_v.size() - 1;

    axom::Array<NURBSPatch<T, 3>, 2> split_patches(num_knot_span_u, num_knot_span_v);
    split_patches(0, 0) = m_alteredPatch;
    for(int i = 0; i < num_knot_span_u - 1; ++i)
    {
      split_patches(i, 0).split_u(knot_vals_u[i + 1], split_patches(i, 0), split_patches(i + 1, 0));
    }

    for(int i = 0; i < num_knot_span_u; ++i)
    {
      for(int j = 0; j < num_knot_span_v - 1; ++j)
      {
        split_patches(i, j).split_v(knot_vals_v[j + 1], split_patches(i, j), split_patches(i, j + 1));
      }
    }

    // Bounding boxes should be defined according to the *pre-expanded* surface,
    //  since the expanded portions are never visibile
    m_oBox = m_alteredPatch.orientedBoundingBox();
    m_bBox.clear();
    for(int i = 0; i < num_knot_span_u; ++i)
    {
      for(int j = 0; j < num_knot_span_v; ++j)
      {
        if(split_patches(i, j).getNumTrimmingCurves() == 0)
        {
          continue;  // Skip patches with no trimming curves
        }

        BezierPatch<T, 3> the_patch;
        if(m_alteredPatch.isRational())
        {
          the_patch = BezierPatch<T, 3>(split_patches(i, j).getControlPoints(),
                                        split_patches(i, j).getWeights(),
                                        split_patches(i, j).getDegree_u(),
                                        split_patches(i, j).getDegree_v());
        }
        else
        {
          the_patch = BezierPatch<T, 3>(split_patches(i, j).getControlPoints(),
                                        split_patches(i, j).getDegree_u(),
                                        split_patches(i, j).getDegree_v());
        }

        BezierPatch<T, 3> p1, p2, p3, p4;
        the_patch.split(0.5, 0.5, p1, p2, p3, p4);

        m_bBox.addBox(p1.boundingBox());
        m_bBox.addBox(p2.boundingBox());
        m_bBox.addBox(p3.boundingBox());
        m_bBox.addBox(p4.boundingBox());
      }
    }

    m_alteredPatch.expandParameterSpace(0.05, 0.05);

    m_curveQuadratureMaps.resize(m_alteredPatch.getNumTrimmingCurves());
  }

  NURBSPatchGWNCache(const BezierPatch<T, 3> a_patch)
    : NURBSPatchGWNCache(NURBSPatch<T, 3>(a_patch))
  { }

  // Mirror the functionality of NURBSPatch so signatures match in GWN evaluation.
  // Allowing only access ensures the memoized information is always accurate
  auto getControlPoints() const { return m_alteredPatch.getControlPoints(); }
  auto getNumControlPoints_u() const { return m_alteredPatch.getNumControlPoints_u(); }
  auto getNumControlPoints_v() const { return m_alteredPatch.getNumControlPoints_v(); }
  auto getWeights() const { return m_alteredPatch.getWeights(); }
  auto getKnots_u() const { return m_alteredPatch.getKnots_u(); }
  auto getKnots_v() const { return m_alteredPatch.getKnots_v(); }
  auto getMinKnot_u() const { return m_alteredPatch.getMinKnot_u(); }
  auto getMaxKnot_u() const { return m_alteredPatch.getMaxKnot_u(); }
  auto getMinKnot_v() const { return m_alteredPatch.getMinKnot_v(); }
  auto getMaxKnot_v() const { return m_alteredPatch.getMaxKnot_v(); }
  auto getTrimmingCurves() const { return m_alteredPatch.getTrimmingCurves(); };
  auto getNumTrimmingCurves() const { return m_alteredPatch.getNumTrimmingCurves(); }
  auto getParameterSpaceDiagonal() const { return m_pboxDiag; }

  // Access precomputed data
  auto getAverageNormal() const { return m_averageNormal; }
  auto boundingBox() const { return m_bBox; }
  auto orientedBoundingBox() const { return m_oBox; }

  /// \brief Creates or accesses the quadrature nodes for a given trimming curve
  TrimmingCurveQuadratureData<T>& getTrimmingCurveQuadratureData(int curveIndex,
                                                                 int quadNPts,
                                                                 int refinementLevel,
                                                                 int refinementIndex) const
  {
    // Check to see if we have already computed the quadrature data for this curve
    auto hash_key = std::make_pair(refinementLevel, refinementIndex);

    if(m_curveQuadratureMaps[curveIndex].find(hash_key) == m_curveQuadratureMaps[curveIndex].end())
    {
      m_curveQuadratureMaps[curveIndex][hash_key] =
        TrimmingCurveQuadratureData<T>(quadNPts,
                                       m_alteredPatch.getTrimmingCurves()[curveIndex],
                                       m_alteredPatch,
                                       refinementLevel,
                                       refinementIndex);
    }

    return m_curveQuadratureMaps[curveIndex][hash_key];
  }

private:
  // The patch is private to present dirtying the cachce by changing the patch,
  //  and because the stored internal patch is altered from the original input
  NURBSPatch<T, 3> m_alteredPatch;

  // Per patch data
  BoundingBox<T, 3> m_bBox;
  OrientedBoundingBox<T, 3> m_oBox;
  axom::primal::Vector<T, 3> m_averageNormal;
  double m_pboxDiag;

  // Per trimming curve data, keyed by (whichRefinementLevel, whichRefinementIndex)
  mutable axom::Array<std::map<std::pair<int, int>, TrimmingCurveQuadratureData<T>>> m_curveQuadratureMaps;
};

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_3D_MEMOIZATION_HPP_
