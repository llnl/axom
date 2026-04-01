// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
#include <unordered_map>
#include <cstdint>
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
    const numerics::QuadratureRule gl_rule = numerics::get_gauss_legendre(quad_npts);
    const auto quad_nodes = gl_rule.nodes();
    const auto quad_weights = gl_rule.weights();

    auto& the_curve = a_patch.getTrimmingCurve(a_curve_index);

    const T curve_min_knot = the_curve.getMinKnot();
    const T curve_max_knot = the_curve.getMaxKnot();

    // Find the right knot span based on the refinement level
    m_span_length = (curve_max_knot - curve_min_knot) / std::pow(2, a_refinementLevel);
    const T span_offset = m_span_length * a_refinementSection;

    m_quadrature_scaled_weights.resize(m_quad_npts);
    m_quadrature_points.resize(m_quad_npts);
    m_quadrature_tangents.resize(m_quad_npts);
    for(int q = 0; q < m_quad_npts; ++q)
    {
      m_quadrature_scaled_weights[q] = quad_weights[q] * m_span_length;

      const T quad_x = quad_nodes[q] * m_span_length + curve_min_knot + span_offset;

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
  double getQuadratureWeight(size_t idx) const { return m_quadrature_scaled_weights[idx]; }
  int getNumPoints() const { return m_quad_npts; }

  axom::ArrayView<const Point<T, 3>> getQuadraturePoints() const
  {
    return m_quadrature_points.view();
  }
  axom::ArrayView<const Vector<T, 3>> getQuadratureTangents() const
  {
    return m_quadrature_tangents.view();
  }
  axom::ArrayView<const double> getQuadratureWeights() const
  {
    return m_quadrature_scaled_weights.view();
  }

private:
  axom::Array<Point<T, 3>> m_quadrature_points;
  axom::Array<Vector<T, 3>> m_quadrature_tangents;
  axom::Array<double> m_quadrature_scaled_weights;
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

    // Cast direction is set to average normal, unless it is near zero
    if(m_averageNormal.norm() < 1e-10)
    {
      // ...unless the average direction is zero
      double theta = axom::utilities::random_real(0.0, 2 * M_PI);
      double u = axom::utilities::random_real(-1.0, 1.0);
      m_castDirection = Vector<T, 3> {sin(theta) * sqrt(1 - u * u), cos(theta) * sqrt(1 - u * u), u};
    }
    else
    {
      m_castDirection = m_averageNormal.unitVector();
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
  const Vector<T, 3>& getCastDirection() const { return m_castDirection; }
  const BoundingBox<T, 3>& boundingBox() const { return m_bBox; }
  const OrientedBoundingBox<T, 3>& orientedBoundingBox() const { return m_oBox; }
  //@}

  /// \brief Creates or accesses the quadrature nodes for a given trimming curve
  TrimmingCurveQuadratureData<T>& getTrimmingCurveQuadratureData(int curveIndex,
                                                                 int quadNPts,
                                                                 int refinementLevel,
                                                                 int refinementIndex) const
  {
    // Cache quadrature data per trimming curve keyed by (refinementLevel, refinementIndex).
    // Note: `quadNPts` is fixed in the 3D winding-number implementation (currently 15).
    const auto make_key = [](int level, int index) -> std::uint64_t {
      return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(level)) << 32) |
        static_cast<std::uint64_t>(static_cast<std::uint32_t>(index));
    };

    const std::uint64_t key = make_key(refinementLevel, refinementIndex);
    auto& curve_map = m_curveQuadratureMaps[curveIndex];

    // Single lookup for both hit and miss; avoids multiple tree traversals and operator[].
    auto [it, inserted] =
      curve_map.try_emplace(key, m_alteredPatch, curveIndex, quadNPts, refinementLevel, refinementIndex);
    return it->second;
  }

private:
  // The patch is private to prevent dirtying the cache by changing the patch,
  //  and because the stored internal patch is altered from the original input
  NURBSPatch<T, 3> m_alteredPatch;

  // Per patch data
  BoundingBox<T, 3> m_bBox;
  OrientedBoundingBox<T, 3> m_oBox;
  Vector<T, 3> m_averageNormal, m_castDirection;
  double m_pboxDiag;

  // Per trimming curve data, keyed by (whichRefinementLevel, whichRefinementIndex)
  mutable axom::Array<std::unordered_map<std::uint64_t, TrimmingCurveQuadratureData<T>>> m_curveQuadratureMaps;
};

}  // namespace detail

/*!
 * \brief Manage an array of NURBSPatchGWNCache<double>
 */
class NURBSPatchCacheManager
{
  using NURBSCache = axom::primal::detail::NURBSPatchGWNCache<double>;
  using NURBSCacheArray = axom::Array<NURBSCache>;
  using NURBSCacheArrayView = axom::ArrayView<NURBSCache>;

  using PatchArrayView = axom::ArrayView<axom::primal::NURBSPatch<double, 3>>;

public:
  NURBSPatchCacheManager() = default;

  NURBSPatchCacheManager(PatchArrayView patchs)
  {
    for(auto& patch : patchs)
    {
      m_nurbs_caches.push_back(NURBSCache(patch));
    }
  }

  /// A view of the manager object.
  struct View
  {
    NURBSCacheArrayView m_view;

    /// Return the NURBSCacheArrayView.
    NURBSCacheArrayView caches() const { return m_view; }
  };

  /// Return a view of this manager to pass into a device function.
  View view() { return View {m_nurbs_caches.view()}; }

  /// Return if the underlying array is empty
  bool empty() const { return m_nurbs_caches.empty(); }

private:
  NURBSCacheArray m_nurbs_caches;
};

template <typename ExecSpace>
struct nurbs_cache_3d_traits
{
  using type = NURBSPatchCacheManager;
};

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
/*!
 * \brief Manage per-thread arrays of NURBSPatchGWNCache<double>
 */
class NURBSPatchCacheManagerOMP
{
  using NURBSCache = axom::primal::detail::NURBSPatchGWNCache<double>;
  using NURBSCachePerThreadArray = axom::Array<axom::Array<NURBSCache>>;
  using NURBSCachePerThreadArrayView = axom::ArrayView<axom::Array<NURBSCache>>;
  using NURBSCacheArrayView = axom::ArrayView<NURBSCache>;

  using PatchArrayView = axom::ArrayView<axom::primal::NURBSPatch<double, 3>>;

public:
  NURBSPatchCacheManagerOMP() = default;

  NURBSPatchCacheManagerOMP(PatchArrayView patches)
  {
    const int nt = omp_get_max_threads();
    m_nurbs_caches.resize(nt);
    auto nurbs_caches_view = m_nurbs_caches.view();

    // Make the first one
    nurbs_caches_view[0].resize(patches.size());
    axom::for_all<axom::OMP_EXEC>(
      patches.size(),
      AXOM_LAMBDA(axom::IndexType i) { nurbs_caches_view[0][i] = NURBSCache(patches[i]); });
    SLIC_INFO("Finished the first construction");
    // Copy the constructed cache to the other threads' copies (less work than construction)
    axom::for_all<axom::OMP_EXEC>(
      1,
      nt,
      AXOM_LAMBDA(axom::IndexType t) { nurbs_caches_view[t].resize(nurbs_caches_view[0].size()); });
    axom::for_all<axom::OMP_EXEC>(
      patches.size(),
      AXOM_LAMBDA(axom::IndexType i) {
        for(int t = 0; t < nt; t++)
        {
          nurbs_caches_view[t][i] = nurbs_caches_view[0][i];
        }
      });
  }

  /// A view of the manager object.
  struct View
  {
    NURBSCachePerThreadArrayView m_views;

    /// Return the NURBSCacheArrayView for the current OMP thread.
    NURBSCacheArrayView caches() const { return m_views[omp_get_thread_num()].view(); }
  };

  /// Return a view of this manager to pass into a device function.
  View view() { return View {m_nurbs_caches.view()}; }

  /// Return if the underlying array is empty
  bool empty() const { return m_nurbs_caches.empty(); }

private:
  NURBSCachePerThreadArray m_nurbs_caches;
};

template <>
struct nurbs_cache_3d_traits<axom::OMP_EXEC>
{
  using type = NURBSPatchCacheManagerOMP;
};
#endif

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_3D_MEMOIZATION_HPP_
