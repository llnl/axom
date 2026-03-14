// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file GWNMethods.hpp
 *
 * \brief Helper classes and type traits for GWN Evaluation methods
 */

#ifndef AXOM_QUEST_GWN_METHODS_HPP_
#define AXOM_QUEST_GWN_METHODS_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/spin.hpp"
#include "axom/primal.hpp"

#include "axom/quest/FastApproximateGWN.hpp"
#include "axom/quest/util/mesh_helpers.hpp"

#include "axom/fmt.hpp"

#include "mfem.hpp"

namespace axom
{
namespace quest
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

///@{
/// \name Query methods for 2D GWN applications

class DirectGWN2D
{
public:
  using CurveArrayType = axom::Array<axom::primal::NURBSCurve<double, 2>>;
  using NURBSCacheArray = axom::Array<axom::primal::detail::NURBSCurveGWNCache<double>>;

  DirectGWN2D() = default;

  /// \brief Define view for NURBS data.
  ///    If memoization is used, allocate a cache for each curve.
  void preprocess(const CurveArrayType& input_curves, bool use_memoization = true)
  {
    m_input_curves_view = input_curves.view();
    if(m_input_curves_view.size() <= 0)
    {
      SLIC_WARNING("Quest: Input shape contains no curves; skipping preprocessing.");
      return;
    }

    axom::utilities::Timer timer(true);
    {
      AXOM_ANNOTATE_SCOPE("preprocessing");
      if(use_memoization)
      {
        m_nurbs_caches.reserve(input_curves.size());
        for(const auto& curv : input_curves)
        {
          m_nurbs_caches.emplace_back(curv);
        }
      }
    }
    timer.stop();
    AXOM_ANNOTATE_METADATA("preprocessing_time", timer.elapsed(), "");
    SLIC_INFO(axom::fmt::format("Direct query preprocessing (loading curves{}): {} s",
                                use_memoization ? " and memoization caches" : "",
                                timer.elapsedTimeInSec()));
  }

  /*!
   * \brief Evaluate the GWN for a query grid at the DOFs of the \a dc query mesh
   *
   * \param [in] dc A query grid to be evaluated at the DOFs
   * \param [in] tol A collection of possible tolerances for GWN evaluation
   */
  void query(mfem::DataCollection& dc, const primal::WindingTolerances& tol)
  {
    if(!dc.HasField("winding") || !dc.HasField("inout"))
    {
      SLIC_WARNING(
        axom::fmt::format("Quest: Input data collection has no field `{}`. Exiting Early.",
                          dc.HasField("winding") ? "inout" : "winding"));
    }

    if(m_input_curves_view.empty())
    {
      SLIC_WARNING("Quest: Skipping query; Input shape not properly initialized.");
      return;
    }

    auto* query_mesh = dc.GetMesh();
    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

    auto query_point = [&query_mesh](int idx) -> axom::primal::Point<double, 2> {
      axom::primal::Point<double, 2> pt;
      query_mesh->GetNode(idx, pt.data());
      return pt;
    };

    axom::utilities::Timer query_timer(true);
    {
      AXOM_ANNOTATE_SCOPE("query");
      const primal::WindingTolerances tol_copy = tol;

      // Use non-memoized form
      if(m_nurbs_caches.empty())
      {
        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType nidx) {
          const auto q = query_point(static_cast<int>(nidx));
          double wn {};
          for(const auto& cache : m_input_curves_view)
          {
            wn += axom::primal::winding_number(q, cache, tol_copy.edge_tol, tol_copy.EPS);
          }
          winding[static_cast<int>(nidx)] = wn;
          inout[static_cast<int>(nidx)] = std::lround(wn);
        });
      }
      else  // Use memoized form
      {
        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType nidx) {
          const auto q = query_point(static_cast<int>(nidx));
          double wn {};
          for(const auto& cache : m_nurbs_caches)
          {
            wn += axom::primal::winding_number(q, cache, tol_copy.edge_tol, tol_copy.EPS);
          }
          winding[static_cast<int>(nidx)] = wn;
          inout[static_cast<int>(nidx)] = std::lround(wn);
        });
      }
    }
    query_timer.stop();

    const double query_time_s = query_timer.elapsed();
    const double ms_per_query = query_timer.elapsedTimeInMilliSec() / num_query_points;
    SLIC_INFO(axom::fmt::format(
      axom::utilities::locale(),
      "Querying {:L} samples in winding number field with{} memoization took {:.3Lf} seconds"
      " (@ {:.0Lf} queries per second; {:.6Lf} ms per query)",
      num_query_points,
      m_nurbs_caches.empty() ? "out" : "",
      query_time_s,
      num_query_points / query_time_s,
      ms_per_query));
    AXOM_ANNOTATE_METADATA("query_points", num_query_points, "");
    AXOM_ANNOTATE_METADATA("query_time", query_time_s, "");
  }

private:
  axom::ArrayView<const axom::primal::NURBSCurve<double, 2>> m_input_curves_view;
  NURBSCacheArray m_nurbs_caches;
};

template <int ORDER>
class PolylineGWN2D
{
public:
  using Point2D = axom::primal::Point<double, 2>;
  using BoxType = axom::primal::BoundingBox<double, 2>;
  using CurveArrayType = axom::Array<axom::primal::NURBSCurve<double, 2>>;
  using SegmentType = axom::primal::Segment<double, 2>;
  using GWNMoments = axom::quest::GWNMomentData<double, 2, ORDER>;

  PolylineGWN2D() = default;

  /// \brief Load polyline data into primal::Segments.
  ///    If fast-approximation is used, construct BVH
  void preprocess(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* poly_mesh,
                  bool useDirectEval)
  {
    if(poly_mesh == nullptr || poly_mesh->getNumberOfCells() <= 0)
    {
      SLIC_WARNING("Quest: Input mesh contains no segments; skipping preprocessing.");
      return;
    }

    axom::utilities::Timer timer(true);
    axom::utilities::Timer stage_timer(false);

    AXOM_ANNOTATE_SCOPE("preprocessing");

    stage_timer.start();
    {
      AXOM_ANNOTATE_SCOPE("extract_segments");

      m_segments.resize(poly_mesh->getNumberOfCells());
      auto segments_view = m_segments.view();

      axom::mint::for_all_cells<axom::SEQ_EXEC, axom::mint::xargs::coords>(
        poly_mesh,
        AXOM_LAMBDA(axom::IndexType cellIdx,
                    const axom::numerics::Matrix<double>& coords,
                    [[maybe_unused]] const axom::IndexType* nodeIds) {
          segments_view[cellIdx] =
            SegmentType {Point2D {coords(0, 0), coords(1, 0)}, Point2D {coords(0, 1), coords(1, 1)}};
        });
    }
    stage_timer.stop();
    SLIC_INFO(axom::fmt::format("  Preprocessing stage (extract_segments): {} s",
                                stage_timer.elapsedTimeInSec()));

    // If direct evaluation is preferred, skip BVH initialization
    if(!useDirectEval)
    {
      stage_timer.reset();
      stage_timer.start();
      {
        AXOM_ANNOTATE_SCOPE("bvh_init");
        const int nlines = m_segments.size();
        axom::Array<BoxType> aabbs(nlines, nlines);
        auto aabbs_view = aabbs.view();
        const auto segments_view = m_segments.view();

        axom::for_all<axom::SEQ_EXEC>(
          nlines,
          AXOM_LAMBDA(axom::IndexType i) {
            aabbs_view[i] = BoxType {segments_view[i].source(), segments_view[i].target()};
          });
        m_bvh.initialize(aabbs_view, nlines);
      }
      stage_timer.stop();
      SLIC_INFO(
        axom::fmt::format("  Preprocessing stage (bvh): {} s", stage_timer.elapsedTimeInSec()));

      stage_timer.reset();
      stage_timer.start();
      {
        AXOM_ANNOTATE_SCOPE("moments");
        const auto segments_view = m_segments.view();

        auto compute_moments = [segments_view](std::int32_t currentNode,
                                               const std::int32_t* leafNodes) -> GWNMoments {
          const auto idx = leafNodes[currentNode];
          return GWNMoments(segments_view[idx]);
        };

        const auto traverser = m_bvh.getTraverser();
        m_internal_moments = traverser.reduce_tree<axom::SEQ_EXEC, GWNMoments>(compute_moments);
      }
      stage_timer.stop();
      SLIC_INFO(
        axom::fmt::format("  Preprocessing stage (moments): {} s", stage_timer.elapsedTimeInSec()));
    }
    timer.stop();
    SLIC_INFO(axom::fmt::format("Total preprocessing: {} s", timer.elapsedTimeInSec()));
  }

  /*!
   * \brief Evaluate the GWN for a query grid at the DOFs of the \a dc query mesh
   *
   * \param [in] dc A query grid to be evaluated at the DOFs
   * \param [in] tol A collection of possible tolerances for GWN evaluation
   */
  void query(mfem::DataCollection& dc, const primal::WindingTolerances& tol)
  {
    if(!dc.HasField("winding") || !dc.HasField("inout"))
    {
      SLIC_WARNING(
        axom::fmt::format("Quest: Input data collection has no field `{}`. Exiting Early.",
                          dc.HasField("winding") ? "inout" : "winding"));
    }

    if(m_segments.empty())
    {
      SLIC_WARNING("Quest: Skipping query; Input shape not properly initialized.");
      return;
    }

    const auto* query_mesh = dc.GetMesh();
    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

    auto query_point = [query_mesh](axom::IndexType idx) -> Point2D {
      axom::primal::Point<double, 2> pt({0., 0.});
      query_mesh->GetNode(static_cast<int>(idx), pt.data());
      return pt;
    };

    axom::utilities::Timer query_timer(true);
    {
      AXOM_ANNOTATE_SCOPE("query");

      const auto segments_view = m_segments.view();
      const primal::WindingTolerances tol_copy = tol;

      // Use fast approximation
      if(m_bvh.isInitialized())
      {
        const auto traverser = m_bvh.getTraverser();
        const auto internal_moments_view = m_internal_moments.view();

        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType index) {
          const double wn = axom::quest::fast_approximate_winding_number(query_point(index),
                                                                         traverser,
                                                                         segments_view,
                                                                         internal_moments_view,
                                                                         tol_copy);

          winding[static_cast<int>(index)] = wn;
          inout[static_cast<int>(index)] = std::lround(wn);
        });
      }
      // Use direct formula
      else
      {
        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType index) {
          const axom::primal::Point<double, 2> q = query_point(static_cast<int>(index));
          double wn {};
          for(const auto& seg : m_segments)
          {
            wn += axom::primal::winding_number(q, seg, tol_copy.edge_tol);
          }

          winding[static_cast<int>(index)] = wn;
          inout[static_cast<int>(index)] = std::lround(wn);
        });
      }
    }
    query_timer.stop();

    const double query_time_s = query_timer.elapsed();
    const double ms_per_query = query_timer.elapsedTimeInMilliSec() / num_query_points;
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Querying {:L} samples in winding number field took {:.3Lf} seconds"
                                " (@ {:.0Lf} queries per second; {:.5Lf} ms per query)",
                                num_query_points,
                                query_time_s,
                                num_query_points / query_time_s,
                                ms_per_query));
    AXOM_ANNOTATE_METADATA("query_points", num_query_points, "");
    AXOM_ANNOTATE_METADATA("query_time", query_time_s, "");
  }

private:
  // For the procsesed input curves/BVH leaf nodes
  axom::Array<SegmentType> m_segments;

  // Only needed for fast approximation method
  axom::Array<GWNMoments> m_internal_moments;
  axom::spin::BVH<2, axom::SEQ_EXEC> m_bvh;
};
///@}

///@{
/// \name Query methods for 3D GWN applications
class DirectGWN3D
{
public:
  using PatchArrayType = axom::Array<axom::primal::NURBSPatch<double, 3>>;
  using NURBSCacheArray = axom::Array<axom::primal::detail::NURBSPatchGWNCache<double>>;

  DirectGWN3D() = default;

  /// \brief Define view for NURBS data.
  ///    If memoization is used, allocate a cache for each patch.
  void preprocess(const PatchArrayType& input_patches, bool use_memoization = true)
  {
    m_input_patches_view = input_patches.view();
    if(m_input_patches_view.size() <= 0)
    {
      SLIC_WARNING("Quest: Input shape contains no patches; skipping preprocessing.");
      return;
    }

    axom::utilities::Timer timer(true);
    {
      AXOM_ANNOTATE_SCOPE("preprocessing");
      if(use_memoization)
      {
        m_nurbs_caches.reserve(input_patches.size());
        for(const auto& patch : input_patches)
        {
          m_nurbs_caches.emplace_back(patch);
        }
      }
    }
    timer.stop();
    AXOM_ANNOTATE_METADATA("preprocessing_time", timer.elapsed(), "");
    SLIC_INFO(axom::fmt::format("Direct query preprocessing (loading surfaces{}): {} s",
                                use_memoization ? " and memoization caches" : "",
                                timer.elapsedTimeInSec()));
  }

  /*!
   * \brief Evaluate the GWN for a query grid at the DOFs of the \a dc query mesh
   *
   * \param [in] dc A query grid to be evaluated at the DOFs
   * \param [in] tol A collection of possible tolerances for GWN evaluation
   * \param [in] slice_z If the dc mesh is 2D, the GWN will be evaluated on a slice 
   *                      parallel to the x-y plane with this offset on the z-axis
   */
  void query(mfem::DataCollection& dc,
             const primal::WindingTolerances& tol,
             const double slice_z = 0.0) const
  {
    if(!dc.HasField("winding") || !dc.HasField("inout"))
    {
      SLIC_WARNING(
        axom::fmt::format("Quest: Input data collection has no field `{}`. Exiting Early.",
                          dc.HasField("winding") ? "inout" : "winding"));
    }

    if(m_input_patches_view.empty())
    {
      SLIC_WARNING("Quest: Skipping query; Input shape not properly initialized.");
      return;
    }

    auto* query_mesh = dc.GetMesh();
    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

    auto query_point = [query_mesh, slice_z](int idx) -> axom::primal::Point<double, 3> {
      axom::primal::Point<double, 3> pt({0., 0., slice_z});
      query_mesh->GetNode(idx, pt.data());
      return pt;
    };

    axom::utilities::Timer query_timer(true);
    {
      AXOM_ANNOTATE_SCOPE("query");
      const primal::WindingTolerances tol_copy = tol;

      // Use non-memoized form
      if(m_nurbs_caches.empty())
      {
        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType nidx) {
          const auto q = query_point(static_cast<int>(nidx));
          double wn {};
          for(const auto& patch : m_input_patches_view)
          {
            wn += axom::primal::winding_number(q,
                                               patch,
                                               tol_copy.edge_tol,
                                               tol_copy.ls_tol,
                                               tol_copy.quad_tol,
                                               tol_copy.disk_size,
                                               tol_copy.EPS);
          }
          winding[static_cast<int>(nidx)] = wn;
          inout[static_cast<int>(nidx)] = std::lround(wn);
        });
      }
      else  // Use memoized form
      {
        const auto nurbs_patches_view = m_nurbs_caches.view();
        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType nidx) {
          const auto q = query_point(static_cast<int>(nidx));
          double wn {};
          for(const auto& cache : nurbs_patches_view)
          {
            wn += axom::primal::winding_number(q,
                                               cache,
                                               tol_copy.edge_tol,
                                               tol_copy.ls_tol,
                                               tol_copy.quad_tol,
                                               tol_copy.disk_size,
                                               tol_copy.EPS);
          }
          winding[static_cast<int>(nidx)] = wn;
          inout[static_cast<int>(nidx)] = std::lround(wn);
        });
      }
    }
    query_timer.stop();

    const double query_time_s = query_timer.elapsed();
    const double ms_per_query = query_timer.elapsedTimeInMilliSec() / num_query_points;
    SLIC_INFO(axom::fmt::format(
      axom::utilities::locale(),
      "Querying {:L} samples in winding number field with{} memoization took {:.3Lf} seconds"
      " (@ {:.0Lf} queries per second; {:.6Lf} ms per query)",
      num_query_points,
      m_nurbs_caches.empty() ? "out" : "",
      query_time_s,
      num_query_points / query_time_s,
      ms_per_query));
    AXOM_ANNOTATE_METADATA("query_points", num_query_points, "");
    AXOM_ANNOTATE_METADATA("query_time", query_time_s, "");
  }

private:
  axom::ArrayView<const axom::primal::NURBSPatch<double, 3>> m_input_patches_view;
  NURBSCacheArray m_nurbs_caches;
};

template <int ORDER>
class TriangleGWN3D
{
public:
  using Point3D = axom::primal::Point<double, 3>;
  using BoxType = axom::primal::BoundingBox<double, 3>;
  using TriangleType = axom::primal::Triangle<double, 3>;
  using GWNMoments = axom::quest::GWNMomentData<double, 3, ORDER>;

  TriangleGWN3D() = default;

  /// \brief Load mesh data into primal::Triangles.
  ///    If fast-approximation is used, construct BVH
  void preprocess(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* tri_mesh, bool useDirectEval)
  {
    axom::utilities::Timer timer(true);
    axom::utilities::Timer stage_timer(false);

    AXOM_ANNOTATE_SCOPE("preprocessing");

    const auto ntris = tri_mesh->getNumberOfCells();
    if(ntris <= 0)
    {
      SLIC_WARNING("Quest: Input mesh contains no triangles; skipping preprocessing.");
      return;
    }

    stage_timer.start();
    {
      AXOM_ANNOTATE_SCOPE("extract_triangles");

      // Iterate over mesh nodes and get a bounding box for the shape
      BoxType shape_bbox;
      BoxType* shape_bbox_ptr = &shape_bbox;
      axom::mint::for_all_nodes<axom::SEQ_EXEC, axom::mint::xargs::xyz>(
        tri_mesh,
        AXOM_LAMBDA(axom::IndexType, double x, double y, double z) {
          shape_bbox_ptr->addPoint(Point3D {x, y, z});
        });

      // Extract the triangles from the mesh into axom primitives,
      //  scaled and translated so that `shape_box` is centered at the origin
      //  and has roughly unit volume. Otherwise, small triangles introduce numerical issues
      m_shape_center = shape_bbox.getCentroid();
      const auto longest_dim = shape_bbox.getLongestDimension();
      m_scale = shape_bbox.getMax()[longest_dim] - shape_bbox.getMin()[longest_dim];

      m_triangles.resize(ntris);
      auto triangles_view = m_triangles.view();
      axom::mint::for_all_cells<axom::SEQ_EXEC, axom::mint::xargs::coords>(
        tri_mesh,
        AXOM_LAMBDA(axom::IndexType cellIdx,
                    const axom::numerics::Matrix<double>& coords,
                    [[maybe_unused]] const axom::IndexType* nodeIds) {
          const auto& ctr = m_shape_center;
          const auto& scl = m_scale;
          triangles_view[cellIdx] = TriangleType {Point3D {(coords(0, 0) - ctr[0]) / scl,
                                                           (coords(1, 0) - ctr[1]) / scl,
                                                           (coords(2, 0) - ctr[2]) / scl},
                                                  Point3D {(coords(0, 1) - ctr[0]) / scl,
                                                           (coords(1, 1) - ctr[1]) / scl,
                                                           (coords(2, 1) - ctr[2]) / scl},
                                                  Point3D {(coords(0, 2) - ctr[0]) / scl,
                                                           (coords(1, 2) - ctr[1]) / scl,
                                                           (coords(2, 2) - ctr[2]) / scl}};
        });
    }
    stage_timer.stop();
    SLIC_INFO(axom::fmt::format("  Preprocessing stage (extract_triangles): {} s",
                                stage_timer.elapsedTimeInSec()));

    // If direct evaluation is preferred, skip BVH initialization
    if(!useDirectEval)
    {
      stage_timer.reset();
      stage_timer.start();
      {
        AXOM_ANNOTATE_SCOPE("bvh_init");
        axom::Array<BoxType> aabbs(ntris, ntris);
        auto aabbs_view = aabbs.view();
        const auto triangles_view = m_triangles.view();

        axom::for_all<axom::SEQ_EXEC>(
          ntris,
          AXOM_LAMBDA(axom::IndexType i) {
            aabbs_view[i] =
              BoxType {triangles_view[i][0], triangles_view[i][1], triangles_view[i][2]};
          });
        m_bvh.initialize(aabbs_view, ntris);
      }
      stage_timer.stop();
      SLIC_INFO(
        axom::fmt::format("  Preprocessing stage (bvh): {} s", stage_timer.elapsedTimeInSec()));

      stage_timer.reset();
      stage_timer.start();
      {
        AXOM_ANNOTATE_SCOPE("moments");
        const auto triangles_view = m_triangles.view();

        auto compute_moments = [triangles_view](std::int32_t currentNode,
                                                const std::int32_t* leafNodes) -> GWNMoments {
          const auto idx = leafNodes[currentNode];
          return GWNMoments(triangles_view[idx]);
        };

        const auto traverser = m_bvh.getTraverser();
        m_internal_moments = traverser.reduce_tree<axom::SEQ_EXEC, GWNMoments>(compute_moments);
      }
      stage_timer.stop();
      SLIC_INFO(
        axom::fmt::format("  Preprocessing stage (moments): {} s", stage_timer.elapsedTimeInSec()));
    }
    timer.stop();

    SLIC_INFO(axom::fmt::format("Total preprocessing: {} s", timer.elapsedTimeInSec()));
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
    if(!dc.HasField("winding") || !dc.HasField("inout"))
    {
      SLIC_WARNING(
        axom::fmt::format("Quest: Skipping query; Input data collection has no field `{}`.",
                          dc.HasField("winding") ? "inout" : "winding"));
      return;
    }

    if(m_triangles.empty())
    {
      SLIC_WARNING("Quest: Skipping query; Input shape not properly initialized.");
      return;
    }

    const auto* query_mesh = dc.GetMesh();
    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

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

    axom::utilities::Timer query_timer(true);
    {
      AXOM_ANNOTATE_SCOPE("query");

      const auto triangles_view = m_triangles.view();
      const primal::WindingTolerances tol_copy = tol;

      // Use fast approximation
      if(m_bvh.isInitialized())
      {
        const auto traverser = m_bvh.getTraverser();
        const auto internal_moments_view = m_internal_moments.view();

        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType index) {
          const double wn = axom::quest::fast_approximate_winding_number(scaled_query_point(index),
                                                                         traverser,
                                                                         triangles_view,
                                                                         internal_moments_view,
                                                                         tol_copy);

          winding[static_cast<int>(index)] = wn;
          inout[static_cast<int>(index)] = std::lround(wn);
        });
      }
      // Use direct formula
      else
      {
        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType index) {
          const auto q = scaled_query_point(static_cast<int>(index));
          double wn {};
          for(const auto& tri : m_triangles)
          {
            wn += axom::primal::winding_number(q, tri, tol_copy.edge_tol, tol_copy.EPS);
          }

          winding[static_cast<int>(index)] = wn;
          inout[static_cast<int>(index)] = std::lround(wn);
        });
      }
    }
    query_timer.stop();

    const double query_time_s = query_timer.elapsed();
    const double ms_per_query = query_timer.elapsedTimeInMilliSec() / num_query_points;
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Querying {:L} samples in winding number field took {:.3Lf} seconds"
                                " (@ {:.0Lf} queries per second; {:.5Lf} ms per query)",
                                num_query_points,
                                query_time_s,
                                num_query_points / query_time_s,
                                ms_per_query));
    AXOM_ANNOTATE_METADATA("query_points", num_query_points, "");
    AXOM_ANNOTATE_METADATA("query_time", query_time_s, "");
  }

private:
  // For the procsesed input curves/BVH leaf nodes
  axom::Array<TriangleType> m_triangles;

  // Only needed for fast approximation method
  axom::Array<GWNMoments> m_internal_moments;
  axom::spin::BVH<3, axom::SEQ_EXEC> m_bvh;

  // Parameters for normalization
  axom::primal::Point<double, 3> m_shape_center;
  double m_scale;
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

template <int ORDER>
struct gwn_input_traits<axom::quest::PolylineGWN2D<ORDER>>
  : std::integral_constant<GWNInputType, GWNInputType::Polyline>
{ };

template <>
struct gwn_input_traits<axom::quest::DirectGWN2D>
  : std::integral_constant<GWNInputType, GWNInputType::Curve>
{ };

template <int ORDER>
struct gwn_input_traits<axom::quest::TriangleGWN3D<ORDER>>
  : std::integral_constant<GWNInputType, GWNInputType::Triangulation>
{ };

template <>
struct gwn_input_traits<axom::quest::DirectGWN3D>
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

}  // namespace quest
}  // namespace axom
#endif
