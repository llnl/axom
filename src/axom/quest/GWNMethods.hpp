// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file GWNMethods.hpp
 *
 * \brief Helper classes and type traits for GWN Evaluation methods
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/spin.hpp"
#include "axom/primal.hpp"

#include "axom/quest/FastApproximateGWN.hpp"
#include "axom/quest/util/mesh_helpers.hpp"

#include "axom/fmt.hpp"

#include "mfem.hpp"

namespace axom::quest
{

//------------------------------------------------------------------------------
// Query-mesh helpers
//------------------------------------------------------------------------------

/// Helper function to set up the mesh and associated winding and inout fields.
/// Uses an mfem::DataCollection to hold everything together.
void setup_gwn_mesh(mfem::DataCollection& dc, mfem::Mesh* query_mesh, int queryOrder)
{
  AXOM_ANNOTATE_SCOPE("setup_mesh");

  dc.SetOwnData(true);
  dc.SetMesh(query_mesh);

  const int dim = query_mesh->Dimension();

  // Create grid functions for the winding field; will take care of fes and fec memory via MakeOwner()
  auto* winding_fec = new mfem::H1Pos_FECollection(queryOrder, dim);
  auto* winding_fes = new mfem::FiniteElementSpace(query_mesh, winding_fec, 1);
  mfem::GridFunction* winding = new mfem::GridFunction(winding_fes);
  winding->MakeOwner(winding_fec);

  // Create grid functions for the inout field; will take care of fes and fec memory via MakeOwner()
  auto* inout_fec = new mfem::H1Pos_FECollection(queryOrder, dim);
  auto* inout_fes = new mfem::FiniteElementSpace(query_mesh, inout_fec, 1);
  mfem::GridFunction* inout = new mfem::GridFunction(inout_fes);
  inout->MakeOwner(inout_fec);

  dc.RegisterField("winding", winding);
  dc.RegisterField("inout", inout);
}

/// Query grid setup; has some dimension-specific types;
///  if user did not provide a bounding box, use shape bounding box scaled by 10%
template <int NDIMS>
void generate_gwn_query_mesh(mfem::DataCollection& dc,
                             const axom::primal::BoundingBox<double, NDIMS>& shape_bbox,
                             const std::vector<double>& boxMins,
                             const std::vector<double>& boxMaxs,
                             const std::vector<int>& boxResolution,
                             int queryOrder)
{
  AXOM_ANNOTATE_SCOPE("generate_query_mesh");

  const int query_dim = static_cast<int>(boxResolution.size());
  const bool has_query_box = !boxMins.empty();

  constexpr double scale_factor = 1.1;
  if(query_dim == 2)
  {
    using Point2D = primal::Point<double, 2>;
    using BoundingBox2D = primal::BoundingBox<double, 2>;

    const auto query_res = axom::NumericArray<int, 2>(boxResolution.data());
    const auto query_box = has_query_box
      ? BoundingBox2D(Point2D(boxMins.data()), Point2D(boxMaxs.data()))
      : BoundingBox2D(Point2D({shape_bbox.getMin()[0], shape_bbox.getMin()[1]}),
                      Point2D({shape_bbox.getMax()[0], shape_bbox.getMax()[1]}))
          .scale(scale_factor);

    SLIC_INFO(
      axom::fmt::format("Query grid resolution {} within bounding box {}", query_res, query_box));

    mfem::Mesh* query_mesh =
      axom::quest::util::make_cartesian_mfem_mesh_2D(query_box, query_res, queryOrder);

    setup_gwn_mesh(dc, query_mesh, queryOrder);
  }
  else
  {
    using Point3D = primal::Point<double, 3>;
    using BoundingBox3D = primal::BoundingBox<double, 3>;

    const auto query_res = axom::NumericArray<int, 3>(boxResolution.data());
    const auto query_box = has_query_box
      ? BoundingBox3D(Point3D(boxMins.data()), Point3D(boxMaxs.data()))
      : BoundingBox3D(
          Point3D({shape_bbox.getMin()[0], shape_bbox.getMin()[1], shape_bbox.getMin()[2]}),
          Point3D({shape_bbox.getMax()[0], shape_bbox.getMax()[1], shape_bbox.getMax()[2]}))
          .scale(scale_factor);

    SLIC_INFO(
      axom::fmt::format("Query grid resolution {} within bounding box {}", query_res, query_box));

    mfem::Mesh* query_mesh =
      axom::quest::util::make_cartesian_mfem_mesh_3D(query_box, query_res, queryOrder);

    setup_gwn_mesh(dc, query_mesh, queryOrder);
  }
}

//------------------------------------------------------------------------------
// Query classes
//------------------------------------------------------------------------------

namespace detail
{

/*!
 * \class FastGWNQuery
 *
 * \brief Shared acceleration state and evaluation interface for fast GWN queries.
 *
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam NDIMS The spatial dimension of the query and BVH.
 * \tparam ORDER The order of the Taylor expansion used for agglomeration.
 *
 * This class is the common core for the geometry-specific GWN query classes, 
 *   and manages the BVH and internal moment data shared by fast evaluation.
 * Each geometry-specific adapter handles preprocessing its input geometry and maps mesh points
 *   to evaluation points (e.g., accounting for scaling transforms).
 */
template <typename ExecSpace, int NDIMS, int ORDER = 2>
class FastGWNQuery
{
public:
  using PointType = primal::Point<double, NDIMS>;
  using BoundingBoxType = primal::BoundingBox<double, NDIMS>;
  using MomentType = GWNMomentData<double, NDIMS, ORDER>;
  using BVHType = axom::spin::BVH<NDIMS, ExecSpace>;

  FastGWNQuery() = default;

  /// Returns whether this query has acceleration data available.
  bool isInitialized() const { return m_bvh.isInitialized(); }

  /// \brief Clear all acceleration data.
  void clear()
  {
    m_bvh = BVHType {};
    m_internal_moments.clear();
    m_source_primitive_count = 0;
  }

  /*!
   * \brief Build acceleration data from the source geometry used for moments.
   *
   * \tparam SourceGeometry The geometry type used to construct bounding boxes
   *         and `GWNMomentData`.
   * \param [in] source_primitives The processed source primitives, in BVH order.
   *
   * Builds bounding boxes, initializes `m_bvh`, and reduces source-primtiive
   * moments into `m_internal_moments`.
   */
  template <typename SourceGeometry>
  void initialize(const axom::ArrayView<const SourceGeometry>& source_primitives)
  {
    clear();

    if(source_primitives.empty())
    {
      SLIC_WARNING("Quest: Input shape contains no primitives; skipping fast GWN initialization.");
      return;
    }

    m_source_primitive_count = source_primitives.size();
    const auto source_primitive_count = m_source_primitive_count;
    axom::Array<BoundingBoxType> geometry_boxes(axom::ArrayOptions::Uninitialized {},
                                                source_primitive_count);
    auto geometry_boxes_view = geometry_boxes.view();

    {
      AXOM_ANNOTATE_SCOPE("bvh_initialization");
      axom::for_all<ExecSpace>(source_primitive_count, [=] AXOM_HOST_DEVICE(axom::IndexType i) {
        geometry_boxes_view[i] = source_primitives[i].boundingBox();
      });
      m_bvh.initialize(geometry_boxes_view, source_primitive_count);
    }

    {
      AXOM_ANNOTATE_SCOPE("moment_precomputation");
      const auto compute_moments =
        [source_primitives](std::int32_t source_primitive_index) -> MomentType {
        return source_primitive_index >= 0 && source_primitive_index < source_primitives.size()
          ? MomentType(source_primitives[source_primitive_index])
          : MomentType {};
      };
      const auto traverser = m_bvh.getTraverser();
      m_internal_moments = traverser.template reduce_leaves<ExecSpace, MomentType>(compute_moments);
    }
  }

  /*!
   * \brief Evaluate the GWN at every nodal degree of freedom of an MFEM mesh.
   *
   * \tparam EvaluationGeometry An indexed collection of primitives used for direct evaluation.
   *           Note: The type EvaluationGeometry may differ from SourceGeometry, as is the case
   *           with, for example, NURBSCurve vs NURBSCurveGWNCache
   * \tparam QueryPointFunction A callable that maps a nodal DOF index to a
   *         point in the GWN coordinate system.
   * \param [in,out] dc The data collection containing `winding` and `inout`
   *        fields.
   * \param [in] evaluation_primitives Primitive objects corresponding to the source-geometry order.
   * \param [in] tolerances Evaluation tolerances.
   * \param [in] query_point Maps a nodal DOF index to a query point.
   */
  template <typename EvaluationGeometry, typename QueryPointFunction>
  void evaluate(mfem::DataCollection& dc,
                const EvaluationGeometry& evaluation_primitives,
                const primal::WindingTolerances& tolerances,
                QueryPointFunction query_point)
  {
    if(!dc.HasField("winding") || !dc.HasField("inout"))
    {
      SLIC_WARNING(
        axom::fmt::format("Quest: Input data collection has no field `{}`. Exiting early.",
                          dc.HasField("winding") ? "inout" : "winding"));
      return;
    }

    const auto* query_mesh = dc.GetMesh();
    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");
    const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();
    const auto tolerances_copy = tolerances;

    AXOM_ANNOTATE_SCOPE("query");
    if(isInitialized())
    {
      // Use fast approximation
      SLIC_ERROR_IF(
        evaluation_primitives.size() != m_source_primitive_count,
        axom::fmt::format("FastGWNQuery: {} evaluation primitives, but acceleration data was "
                          "initialized for {} source primitives.",
                          evaluation_primitives.size(),
                          m_source_primitive_count));
      const auto traverser = m_bvh.getTraverser();
      const auto internal_moments_view = m_internal_moments.view();
      axom::for_all<ExecSpace>(num_query_points, [=, &winding, &inout](axom::IndexType index) {
        const double gwn = fast_approximate_winding_number(query_point(index),
                                                           traverser,
                                                           evaluation_primitives,
                                                           internal_moments_view,
                                                           tolerances_copy);
        winding[static_cast<int>(index)] = gwn;
        inout[static_cast<int>(index)] = std::lround(gwn);
      });
    }
    else
    {
      // Use direct evaluation
      axom::for_all<ExecSpace>(num_query_points, [=, &winding, &inout](axom::IndexType index) {
        const auto query = query_point(index);
        double gwn = 0.0;
        for(axom::IndexType i = 0; i < evaluation_primitives.size(); ++i)
        {
          gwn += primal::winding_number(query, evaluation_primitives[i], tolerances_copy);
        }
        winding[static_cast<int>(index)] = gwn;
        inout[static_cast<int>(index)] = std::lround(gwn);
      });
    }
  }

private:
  BVHType m_bvh;
  axom::Array<MomentType> m_internal_moments;
  axom::IndexType m_source_primitive_count {0};
};

}  // namespace detail

/*!
 * \class DipoleGWNQuery
 *
 * \brief Evaluates a 2D or 3D dipole field on an MFEM query mesh.
 *
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam NDIMS The dimension of the dipole field, either 2 or 3.
 * \tparam ORDER The order of the Taylor expansion used for agglomeration.
 *
 * Dipoles are already BVH leaves, so preprocessing only retains the input
 * view and optionally constructs the shared fast-GWN acceleration data.
 */
template <typename ExecSpace, int NDIMS, int ORDER = 2>
class DipoleGWNQuery
{
public:
  static_assert(NDIMS == 2 || NDIMS == 3, "Dipole GWN queries are defined only for 2D or 3D.");

  using DipoleType = primal::Dipole<double, NDIMS>;
  using DipoleArrayType = axom::Array<DipoleType>;
  using FastQuery = detail::FastGWNQuery<ExecSpace, NDIMS, ORDER>;

  DipoleGWNQuery() = default;

  /*!
   * \brief Retain dipoles and optionally initialize BVH acceleration data.
   *
   * \param [in] dipoles The source dipoles. They must remain alive while this query is used.
   * \param [in] use_direct_eval If false, construct the BVH and moment data.
   */
  void preprocess(const DipoleArrayType& dipoles, bool use_direct_eval = false)
  {
    // Clear stale data
    m_fast_query.clear();

    m_dipoles_view = dipoles.view();
    if(m_dipoles_view.empty())
    {
      SLIC_WARNING("Quest: Input shape contains no dipoles; skipping preprocessing.");
      return;
    }

    if(!use_direct_eval)
    {
      // Dipoles require no conversion, subdivision, or cache construction.
      m_fast_query.initialize(m_dipoles_view);
    }
  }

  /*!
   * \brief Evaluate the dipole GWN field at the DOFs of a query mesh.
   *
   * \param [in] dc A query grid to be evaluated at the DOFs.
   * \param [in] tolerances A collection of possible tolerances for GWN evaluation.
   * \param [in] slice_z The z offset for a 2D query mesh evaluating a 3D field.
   */
  void query(mfem::DataCollection& dc,
             const primal::WindingTolerances& tolerances,
             const double slice_z = 0.0)
  {
    if(m_dipoles_view.empty())
    {
      SLIC_WARNING("Quest: Skipping query; Input shape not properly initialized.");
      return;
    }

    // Dispatch the query point location based on whether the query is of a 2D domain,
    //  or a 2D slice of a 3D domain.
    const auto* query_mesh = dc.GetMesh();
    const auto query_point = [query_mesh,
                              slice_z](axom::IndexType index) -> primal::Point<double, NDIMS> {
      if constexpr(NDIMS == 2)
      {
        AXOM_UNUSED_VAR(slice_z);  // unused in 2D
        primal::Point<double, 2> point;
        query_mesh->GetNode(static_cast<int>(index), point.data());
        return point;
      }
      else
      {
        primal::Point<double, 3> point {0.0, 0.0, slice_z};
        query_mesh->GetNode(static_cast<int>(index), point.data());
        return point;
      }
    };

    m_fast_query.evaluate(dc, m_dipoles_view, tolerances, query_point);
  }

private:
  axom::ArrayView<const DipoleType> m_dipoles_view;
  FastQuery m_fast_query;
};

///@{
/// \name Query methods for 2D GWN applications

/*!
 * \class NURBSCurveGWNQuery
 *
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam ORDER If agglomeration is used, this is the order of the Taylor expansion.
 *
 * \brief Preprocesses NURBSCurve geometry for GWN evaluation, 
 *         and performs the calculation on the DOFs of an input MFEM mesh.
 * 
 * Possible evaluation modes are
 *  `use_direct_eval` : If true, evaluation is done curve-by-curve.
 *                      If false, evaluation is sped up with agglomeration via Taylor-expansion
 *  `use_memoization` : Caches and re-uses subdivision data for curve evaluations
 */
template <typename ExecSpace, int ORDER = 2>
class NURBSCurveGWNQuery
{
public:
  using CurveType = axom::primal::NURBSCurve<double, 2>;
  using CurveArrayType = axom::Array<CurveType>;
  using NURBSCacheManager = typename axom::primal::nurbs_cache_2d_traits<ExecSpace>::type;
  using FastQuery = detail::FastGWNQuery<ExecSpace, 2, ORDER>;

  NURBSCurveGWNQuery() = default;

  /// \brief Returns whether this query uses memoization during evaluation.
  bool is_memoized() const { return !m_nurbs_cache_mgr.empty(); }

  ///@{
  /// \name Setters for misc algorithm parameters
  void setSubdivisionBboxThreshold(double subdivision_bbox_threshold)
  {
    m_subdivision_bbox_threshold = subdivision_bbox_threshold;
  }

  void setSubdivisionMaxPasses(int subdivision_max_passes)
  {
    m_subdivision_max_passes = subdivision_max_passes;
  }

  void setSubdivisionMaxNumCurves(int subdivision_max_curves)
  {
    m_subdivision_max_curves = subdivision_max_curves;
  }
  ///@}

  /*!
   * \brief Process input curves, optionally building a BVH
   *
   * \param [in] input_curves A view to the input curves
   * \param [in] use_direct_eval If false, use accelerated agglomeration algorithm via BVH
   * \param [in] use_memoization If true, allocate a per-thread cache for each curve
   */
  void preprocess(const CurveArrayType& input_curves,
                  bool use_direct_eval = false,
                  bool use_memoization = true)
  {
    // Clear stale data
    m_fast_query.clear();
    m_nurbs_cache_mgr = NURBSCacheManager {};
    m_subdivided_curves.clear();
    m_processed_curves_view = {};

    m_input_curves_view = input_curves.view();
    if(m_input_curves_view.size() <= 0)
    {
      SLIC_WARNING("Quest: Input shape contains no curves; skipping preprocessing.");
      return;
    }

    AXOM_ANNOTATE_SCOPE("preprocessing");

    if(!use_direct_eval)
    {
      {
        AXOM_ANNOTATE_SCOPE("subdivision");
        m_subdivided_curves = subdivide_curves(m_input_curves_view,
                                               m_subdivision_bbox_threshold,
                                               m_subdivision_max_curves,
                                               m_subdivision_max_passes);
        m_processed_curves_view = m_subdivided_curves.view();
      }

      // Build boxes and moments from the processed curves. A cache, if requested,
      // is only a direct-evaluation representation of these same leaves.
      m_fast_query.initialize(m_processed_curves_view);
    }
    else
    {
      // Without fast-approximation, processing is unnecessary
      m_processed_curves_view = m_input_curves_view;
    }

    if(use_memoization)
    {
      AXOM_ANNOTATE_SCOPE("cache_initialization");
      m_nurbs_cache_mgr = NURBSCacheManager(m_processed_curves_view);
    }
  }

  /*!
   * \brief Evaluate the GWN for a query grid at the DOFs of the \a dc query mesh
   *
   * \param [in] dc A query grid to be evaluated at the DOFs
   * \param [in] tol A collection of possible tolerances for GWN evaluation
   */
  void query(mfem::DataCollection& dc, const primal::WindingTolerances& tol)
  {
    if(m_input_curves_view.empty())
    {
      SLIC_WARNING("Quest: Skipping query; Input shape not properly initialized.");
      return;
    }

    const auto* query_mesh = dc.GetMesh();
    const auto query_point = [query_mesh](int idx) -> axom::primal::Point<double, 2> {
      axom::primal::Point<double, 2> pt;
      query_mesh->GetNode(idx, pt.data());
      return pt;
    };

    if(is_memoized())
    {
      const auto cache_manager_view = m_nurbs_cache_mgr.view();

      // `FastGWNQuery` indexes direct leaves with the BVH's processed-curve
      // indices. The cache manager creates exactly one cache per processed
      // curve, and selects each OMP thread's cache during calls to operator[].
      SLIC_ASSERT(cache_manager_view.size() == m_processed_curves_view.size());
      m_fast_query.evaluate(dc, cache_manager_view, tol, query_point);
    }
    else
    {
      m_fast_query.evaluate(dc, m_processed_curves_view, tol, query_point);
    }
  }

private:
  // For the input curves/BVH leaf nodes
  axom::ArrayView<const CurveType> m_input_curves_view;
  NURBSCacheManager m_nurbs_cache_mgr;

  // For preprocessed curves
  axom::Array<CurveType> m_subdivided_curves;
  axom::ArrayView<const CurveType> m_processed_curves_view;
  FastQuery m_fast_query;

  // Additional algorithm parameters
  double m_subdivision_bbox_threshold {1.0};
  int m_subdivision_max_passes {10};
  int m_subdivision_max_curves {1000000};
};

/*!
 * \class PolylineGWNQuery
 *
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam ORDER If agglomeration is used, this is the order of the Taylor expansion.
 *
 * \brief Preprocesses a linear mesh for GWN evaluation, 
 *         and performs the calculation on the DOFs of an input MFEM mesh.
 * 
 * Possible evaluation modes are
 *  `use_direct_eval` : If true, evaluation is done segment-by-segment.
 *                      If false, evaluation is sped up with agglomeration via Taylor-expansion
 */
template <typename ExecSpace, int ORDER = 2>
class PolylineGWNQuery
{
public:
  using Point2D = axom::primal::Point<double, 2>;
  using SegmentType = axom::primal::Segment<double, 2>;
  using FastQuery = detail::FastGWNQuery<ExecSpace, 2, ORDER>;

  PolylineGWNQuery() = default;

  /*!
   * \brief Process mint::mesh into axom::Segments, optionally building a BVH
   *
   * \param [in] poly_mesh The input mesh
   * \param [in] use_direct_eval If false, use accelerated agglomeration algorithm via BVH
   */
  void preprocess(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* poly_mesh,
                  bool use_direct_eval)
  {
    // Clear stale data
    m_fast_query.clear();
    m_segments.clear();

    if(poly_mesh == nullptr || poly_mesh->getNumberOfCells() <= 0)
    {
      SLIC_WARNING("Quest: Input mesh contains no segments; skipping preprocessing.");
      return;
    }

    AXOM_ANNOTATE_SCOPE("preprocessing");

    {
      AXOM_ANNOTATE_SCOPE("extract_segments");

      m_segments.resize(poly_mesh->getNumberOfCells());
      auto segments_view = m_segments.view();

      axom::mint::for_all_cells<ExecSpace, axom::mint::xargs::coords>(
        poly_mesh,
        [=] AXOM_HOST_DEVICE(axom::IndexType cellIdx,
                             const axom::numerics::Matrix<double>& coords,
                             const axom::IndexType* AXOM_UNUSED_PARAM(nodeIds)) {
          segments_view[cellIdx] =
            SegmentType {Point2D {coords(0, 0), coords(1, 0)}, Point2D {coords(0, 1), coords(1, 1)}};
        });
    }

    // If direct evaluation is preferred, skip BVH initialization
    if(!use_direct_eval)
    {
      const axom::Array<SegmentType>& const_segments = m_segments;
      m_fast_query.initialize(const_segments.view());
    }
  }

  /*!
   * \brief Evaluate the GWN for a query grid at the DOFs of the \a dc query mesh
   *
   * \param [in] dc A query grid to be evaluated at the DOFs
   * \param [in] tol A collection of possible tolerances for GWN evaluation
   */
  void query(mfem::DataCollection& dc, const primal::WindingTolerances& tol)
  {
    if(m_segments.empty())
    {
      SLIC_WARNING("Quest: Skipping query; Input shape not properly initialized.");
      return;
    }

    const auto* query_mesh = dc.GetMesh();
    const auto query_point = [query_mesh](axom::IndexType idx) -> Point2D {
      axom::primal::Point<double, 2> pt({0., 0.});
      query_mesh->GetNode(static_cast<int>(idx), pt.data());
      return pt;
    };

    const axom::Array<SegmentType>& const_segments = m_segments;
    m_fast_query.evaluate(dc, const_segments.view(), tol, query_point);
  }

private:
  // For the procsesed input curves/BVH leaf nodes
  axom::Array<SegmentType> m_segments;
  FastQuery m_fast_query;
};
///@}

///@{
/// \name Query methods for 3D GWN applications

/*!
 * \class NURBSPatchGWNQuery
 *
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam ORDER If agglomeration is used, this is the order of the Taylor expansion.
 *
 * \brief Preprocesses NURBSPatch geometry for GWN evaluation, 
 *         and performs the calculation on the DOFs of an input MFEM mesh.
 * 
 * Possible evaluation modes are
 *  `use_direct_eval` : If true, evaluation is done patch-by-patch.
 *                      If false, evaluation is sped up with agglomeration via Taylor-expansion
 *  `use_memoization` : Caches and re-uses trimming curve quadrature data for patch evaluation
 */
template <typename ExecSpace, int ORDER = 2>
class NURBSPatchGWNQuery
{
public:
  using PatchType = axom::primal::NURBSPatch<double, 3>;
  using PatchArrayType = axom::Array<PatchType>;
  using NURBSCacheManager = typename axom::primal::nurbs_cache_3d_traits<ExecSpace>::type;
  using FastQuery = detail::FastGWNQuery<ExecSpace, 3, ORDER>;

  NURBSPatchGWNQuery() = default;

  /// \brief Returns whether this query uses memoization during evaluation.
  bool is_memoized() const { return !m_nurbs_cache_mgr.empty(); }

  ///@{
  /// \name Setters for misc algorithm parameters
  void setSubdivisionBboxThreshold(double subdivision_bbox_threshold)
  {
    m_subdivision_bbox_threshold = subdivision_bbox_threshold;
  }

  void setSubdivisionMaxPasses(int subdivision_max_passes)
  {
    m_subdivision_max_passes = subdivision_max_passes;
  }

  void setSubdivisionMaxNumPatches(int subdivision_max_patches)
  {
    m_subdivision_max_patches = subdivision_max_patches;
  }
  ///@}

  /*!
   * \brief Process input patches, optionally building a BVH
   *
   * Processing involves "cleaning" input surfaces for more robust GWN evaluation by 
   *  - Normalizing the parameter space of each surface according to the number of knot spans
   *  - Ensuring each input is represented as a trimmed surface
   *  - Clipping the parameter space of each surface to the visible (i.e. trimmed) portion 
   * 
   * \param [in] input_patches A view to the input trimmed NURBS surfaces
   * \param [in] use_direct_eval If false, use accelerated agglomeration algorithm via BVH
   * \param [in] use_memoization If true, allocate a per-thread cache for each patch
   */
  void preprocess(const PatchArrayType& input_patches,
                  bool use_direct_eval = false,
                  bool use_memoization = true)
  {
    // Clear stale data
    m_fast_query.clear();
    m_nurbs_cache_mgr = NURBSCacheManager {};
    m_subdivided_patches.clear();
    m_processed_patches_view = {};

    m_input_patches_view = input_patches.view();
    if(m_input_patches_view.size() <= 0)
    {
      SLIC_WARNING("Quest: Input shape contains no patches; skipping preprocessing.");
      return;
    }

    AXOM_ANNOTATE_SCOPE("preprocessing");

    if(!use_direct_eval)
    {
      {
        AXOM_ANNOTATE_SCOPE("subdivision");
        m_subdivided_patches = subdivide_patches(m_input_patches_view,
                                                 m_subdivision_bbox_threshold,
                                                 m_subdivision_max_patches,
                                                 m_subdivision_max_passes);
        m_processed_patches_view = m_subdivided_patches.view();
      }

      m_fast_query.initialize(m_processed_patches_view);
    }
    else
    {
      // Without fast-approximation, processing is unnecessary, but we still clean the
      //  the trimmed representation for more precise moment calculation
      for(auto& surf : m_input_patches_view)
      {
        auto cleaned = surf.cleanedTrimmedRepresentation();
        if(cleaned.getNumTrimmingCurves() == 0)
        {
          continue;
        }

        m_subdivided_patches.push_back(std::move(cleaned));
      }
      m_processed_patches_view = m_subdivided_patches.view();
    }

    if(use_memoization)
    {
      AXOM_ANNOTATE_SCOPE("cache_initialization");

      // TODO: We now recompute normals and surface areas already considered as moment
      //  data. It would be cheaper to reuse them, but it's a small part of an already
      //  small preprocessing cost that greatly increases the complexity of FastGWNQuery
      // For now, pass views to empty arrays to NURBSCacheManager
      axom::Array<primal::Vector<double, 3>> no_precomputed_normals;
      axom::Array<double> no_precomputed_surface_areas;
      m_nurbs_cache_mgr = NURBSCacheManager(m_processed_patches_view,
                                            no_precomputed_normals.view(),
                                            no_precomputed_surface_areas.view());
    }
  }

  /*!
   * \brief Evaluate the GWN for a query grid at the DOFs of the \a dc query mesh
   *
   * \param [in] dc A query grid to be evaluated at the DOFs
   * \param [in] tol A collection of possible tolerances for GWN evaluation
   * \param [in] slice_z If the dc mesh is 2D, the GWN will be evaluated on a slice 
   *                      parallel to the x-y plane with this offset on the z-axis
   */
  void query(mfem::DataCollection& dc, const primal::WindingTolerances& tol, const double slice_z = 0.0)
  {
    if(m_input_patches_view.empty())
    {
      SLIC_WARNING("Quest: Skipping query; Input shape not properly initialized.");
      return;
    }

    const auto* query_mesh = dc.GetMesh();
    const auto query_point = [query_mesh, slice_z](int idx) -> axom::primal::Point<double, 3> {
      axom::primal::Point<double, 3> pt({0., 0., slice_z});
      query_mesh->GetNode(idx, pt.data());
      return pt;
    };

    if(is_memoized())
    {
      const auto cache_manager_view = m_nurbs_cache_mgr.view();

      // `FastGWNQuery` indexes direct leaves with the BVH's processed-patch
      // indices. The cache manager creates exactly one cache per processed
      // patch, and selects each OMP thread's cache during calls to operator[].
      SLIC_ASSERT(cache_manager_view.size() == m_processed_patches_view.size());
      m_fast_query.evaluate(dc, cache_manager_view, tol, query_point);
    }
    else
    {
      SLIC_WARNING(
        "Quest warning: Patch GWN evaluation is prohibitively slow without memoization.");
      m_fast_query.evaluate(dc, m_processed_patches_view, tol, query_point);
    }
  }

private:
  // For the input patches and direct evaluation leaves.
  axom::ArrayView<const PatchType> m_input_patches_view;
  NURBSCacheManager m_nurbs_cache_mgr;

  // For preprocessed patches
  axom::Array<PatchType> m_subdivided_patches;
  axom::ArrayView<const PatchType> m_processed_patches_view;
  FastQuery m_fast_query;

  // Additional algorithm parameters
  double m_subdivision_bbox_threshold {1.0};
  int m_subdivision_max_passes {10};
  int m_subdivision_max_patches {10000};
};

/*!
 * \class TriangleGWNQuery
 * 
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam ORDER If agglomeration is used, this is the order of the Taylor expansion.
 *
 * \brief Preprocesses a triangle mesh for GWN evaluation, 
 *         and performs the calculation on the DOFs of an input MFEM mesh.
 * 
 * Possible evaluation modes are
 *  `use_direct_eval` : If true, evaluation is done triangle-by-triangle.
 *                      If false, evaluation is sped up with agglomeration via Taylor-expansion
 */
template <typename ExecSpace, int ORDER = 2>
class TriangleGWNQuery
{
public:
  using Point3D = axom::primal::Point<double, 3>;
  using BoxType = axom::primal::BoundingBox<double, 3>;
  using TriangleType = axom::primal::Triangle<double, 3>;
  using FastQuery = detail::FastGWNQuery<ExecSpace, 3, ORDER>;

  TriangleGWNQuery() = default;

  /*!
   * \brief Load triangles from mint::mesh to primal::Triangles, optionally building a BVH
   *
   * \param [in] tri_mesh The input mesh
   * \param [in] use_direct_eval If false, use accelerated agglomeration algorithm via BVH
   */
  void preprocess(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* tri_mesh,
                  bool use_direct_eval)
  {
    // Clear stale data
    m_fast_query.clear();
    m_triangles.clear();
    m_shape_center = Point3D {};
    m_scale = 1.0;

    AXOM_ANNOTATE_SCOPE("preprocessing");

    if(tri_mesh == nullptr || tri_mesh->getNumberOfCells() <= 0)
    {
      SLIC_WARNING("Quest: Input mesh contains no triangles; skipping preprocessing.");
      return;
    }

    {
      AXOM_ANNOTATE_SCOPE("extract_triangles");

      // Iterate over mesh nodes and get a bounding box for the shape
      BoxType shape_bbox;
      BoxType* shape_bbox_ptr = &shape_bbox;
      axom::mint::for_all_nodes<ExecSpace, axom::mint::xargs::xyz>(
        tri_mesh,
        [=] AXOM_HOST_DEVICE(axom::IndexType, double x, double y, double z) {
          shape_bbox_ptr->addPoint(Point3D {x, y, z});
        });

      // Extract the triangles from the mesh into axom primitives,
      //  scaled and translated so that `shape_box` is centered at the origin
      //  and has roughly unit volume. Otherwise, small triangles introduce numerical issues
      m_shape_center = shape_bbox.getCentroid();
      const auto longest_dim = shape_bbox.getLongestDimension();
      m_scale = shape_bbox.getMax()[longest_dim] - shape_bbox.getMin()[longest_dim];
      const auto shape_center = m_shape_center;
      const double scale = m_scale;

      m_triangles.resize(tri_mesh->getNumberOfCells());
      auto triangles_view = m_triangles.view();
      axom::mint::for_all_cells<ExecSpace, axom::mint::xargs::coords>(
        tri_mesh,
        [=] AXOM_HOST_DEVICE(axom::IndexType cellIdx,
                             const axom::numerics::Matrix<double>& coords,
                             const axom::IndexType* AXOM_UNUSED_PARAM(nodeIds)) {
          triangles_view[cellIdx] =
            TriangleType {Point3D {(coords(0, 0) - shape_center[0]) / scale,
                                   (coords(1, 0) - shape_center[1]) / scale,
                                   (coords(2, 0) - shape_center[2]) / scale},
                          Point3D {(coords(0, 1) - shape_center[0]) / scale,
                                   (coords(1, 1) - shape_center[1]) / scale,
                                   (coords(2, 1) - shape_center[2]) / scale},
                          Point3D {(coords(0, 2) - shape_center[0]) / scale,
                                   (coords(1, 2) - shape_center[1]) / scale,
                                   (coords(2, 2) - shape_center[2]) / scale}};
        });
    }

    // If direct evaluation is preferred, skip BVH initialization
    if(!use_direct_eval)
    {
      const axom::Array<TriangleType>& const_triangles = m_triangles;
      m_fast_query.initialize(const_triangles.view());
    }
  }

  /*!
   * \brief Evaluate the GWN for a query grid at the DOFs of the \a dc query mesh
   *
   * \param [in] dc A query grid to be evaluated at the DOFs
   * \param [in] tol A collection of possible tolerances for GWN evaluation
   * \param [in] slice_z If the dc mesh is 2D, the GWN will be evaluated on a slice 
   *                      parallel to the x-y plane with this offset on the z-axis
   */
  void query(mfem::DataCollection& dc, const primal::WindingTolerances& tol, const double slice_z = 0.0)
  {
    if(m_triangles.empty())
    {
      SLIC_WARNING("Quest: Skipping query; Input shape not properly initialized.");
      return;
    }

    const auto* query_mesh = dc.GetMesh();

    // Get the query point from the mesh, scaled to the proper normalization
    const auto& ctr = m_shape_center;
    const auto& scl = m_scale;
    auto scaled_query_point =
      [query_mesh, slice_z, ctr, scl](axom::IndexType idx) -> axom::primal::Point<double, 3> {
      axom::primal::Point<double, 3> pt({0., 0., slice_z});
      query_mesh->GetNode(static_cast<int>(idx), pt.data());
      pt.array() = (pt.array() - ctr.array()) / scl;
      return pt;
    };

    const axom::Array<TriangleType>& const_triangles = m_triangles;
    m_fast_query.evaluate(dc, const_triangles.view(), tol, scaled_query_point);
  }

private:
  // For the processed input curves/BVH leaf nodes
  axom::Array<TriangleType> m_triangles;

  FastQuery m_fast_query;

  // Parameters for normalization
  axom::primal::Point<double, 3> m_shape_center;
  double m_scale {1.0};
};
///@}

//------------------------------------------------------------------------------
// Type Traits
//------------------------------------------------------------------------------

enum class GWNInputType
{
  Curve,
  Polyline,
  Surface,
  Triangulation
};

template <typename GWNQueryType>
struct gwn_input_traits;

template <typename ExecSpace, int ORDER>
struct gwn_input_traits<axom::quest::PolylineGWNQuery<ExecSpace, ORDER>>
  : std::integral_constant<GWNInputType, GWNInputType::Polyline>
{ };

template <typename ExecSpace, int ORDER>
struct gwn_input_traits<axom::quest::NURBSCurveGWNQuery<ExecSpace, ORDER>>
  : std::integral_constant<GWNInputType, GWNInputType::Curve>
{ };

template <typename ExecSpace, int ORDER>
struct gwn_input_traits<axom::quest::TriangleGWNQuery<ExecSpace, ORDER>>
  : std::integral_constant<GWNInputType, GWNInputType::Triangulation>
{ };

template <typename ExecSpace, int ORDER>
struct gwn_input_traits<axom::quest::NURBSPatchGWNQuery<ExecSpace, ORDER>>
  : std::integral_constant<GWNInputType, GWNInputType::Surface>
{ };

template <typename GWNQueryType>
inline constexpr GWNInputType gwn_input_type_v = gwn_input_traits<GWNQueryType>::value;

//------------------------------------------------------------------------------
// Compute postprocessing stats
//------------------------------------------------------------------------------

struct FieldStats
{
  double dof_l2 {};
  double dof_linf {};
  double l2 {};
  double min {};
  double max {};
};

FieldStats compute_field_stats(const mfem::GridFunction& gf)
{
  FieldStats s {};

  s.dof_l2 = gf.Norml2();
  s.dof_linf = gf.Normlinf();
  s.min = gf.Min();
  s.max = gf.Max();

  // Compute L2 norm over the physical domain: sqrt( Integral gf^2 dOmega )
  // We do this by assembling b_i = Integral phi_i * gf dOmega, then taking dot(gf, b).
  auto* fes = const_cast<mfem::FiniteElementSpace*>(gf.FESpace());
  auto* gf_ptr = const_cast<mfem::GridFunction*>(&gf);
  mfem::GridFunctionCoefficient gf_coeff(gf_ptr);
  mfem::LinearForm gf_sq_form(fes);
  gf_sq_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(gf_coeff));
  gf_sq_form.Assemble();

  const double integral_gf_sq = gf * gf_sq_form;
  s.l2 = std::sqrt(std::max(0.0, integral_gf_sq));

  return s;
}

struct IntegralStats
{
  double integral {};
  double domain_volume {};
};

IntegralStats compute_integrals(const mfem::GridFunction& gf)
{
  IntegralStats s {};

  auto* fes = const_cast<mfem::FiniteElementSpace*>(gf.FESpace());
  mfem::ConstantCoefficient one(1.0);
  mfem::LinearForm vol_form(fes);
  vol_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
  vol_form.Assemble();

  s.integral = gf * vol_form;

  mfem::GridFunction unity(fes);
  unity.ProjectCoefficient(one);
  s.domain_volume = unity * vol_form;

  return s;
}

}  // namespace axom::quest
