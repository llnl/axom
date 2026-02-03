// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_winding_number.cpp
 * \brief Example that computes the winding number of a grid of points
 * against a collection of 2D parametric rational curves.
 * Supports MFEM meshes in the cubic positive Bernstein basis or the (rational)
 * NURBS basis.
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/quest.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "mfem.hpp"

namespace primal = axom::primal;
using Point2D = primal::Point<double, 2>;
using NURBSCurve2D = primal::NURBSCurve<double, 2>;
using CurveGWNCache = primal::detail::NURBSCurveGWNCache<double>;
using BoundingBox2D = primal::BoundingBox<double, 2>;

//------------------------------------------------------------------------------
// Query-mesh helpers
//------------------------------------------------------------------------------

/// Helper function to set up the mesh and associated winding and inout fields.
/// Uses an mfem::DataCollection to hold everything together.
void setup_mesh(mfem::DataCollection& dc, mfem::Mesh* query_mesh, int queryOrder)
{
  AXOM_ANNOTATE_SCOPE("setup_mesh");

  dc.SetOwnData(true);
  dc.SetMesh(query_mesh);

  constexpr int DIM = 2;

  // Create grid functions for the winding field; will take care of fes and fec memory via MakeOwner()
  auto* winding_fec = new mfem::H1_FECollection(queryOrder, DIM);
  auto* winding_fes = new mfem::FiniteElementSpace(query_mesh, winding_fec, 1);
  mfem::GridFunction* winding = new mfem::GridFunction(winding_fes);
  winding->MakeOwner(winding_fec);

  // Create grid functions for the inout field; will take care of fes and fec memory via MakeOwner()
  auto* inout_fec = new mfem::H1_FECollection(queryOrder, DIM);
  auto* inout_fes = new mfem::FiniteElementSpace(query_mesh, inout_fec, 1);
  mfem::GridFunction* inout = new mfem::GridFunction(inout_fes);
  inout->MakeOwner(inout_fec);

  dc.RegisterField("winding", winding);
  dc.RegisterField("inout", inout);
}

/// Query grid setup; if user did not provide a bounding box, user input bounding box scaled by 10%
void generate_query_mesh(mfem::DataCollection& dc,
                         const BoundingBox2D& bbox,
                         const std::vector<double>& boxMins,
                         const std::vector<double>& boxMaxs,
                         const std::vector<int>& boxResolution,
                         int queryOrder)
{
  const bool has_query_box = !boxMins.empty();
  constexpr double scale_factor = 1.1;

  using Point2D = primal::Point<double, 2>;
  const auto query_res = axom::NumericArray<int, 2>(boxResolution.data());
  const auto query_box = has_query_box
    ? BoundingBox2D(Point2D(boxMins.data()), Point2D(boxMaxs.data()))
    : BoundingBox2D(Point2D({bbox.getMin()[0], bbox.getMin()[1]}),
                    Point2D({bbox.getMax()[0], bbox.getMax()[1]}))
        .scale(scale_factor);

  SLIC_INFO(
    axom::fmt::format("Query grid resolution {} within bounding box {}", query_res, query_box));

  mfem::Mesh* query_mesh =
    axom::quest::util::make_cartesian_mfem_mesh_2D(bbox, query_res, queryOrder);
  dc.SetMesh(query_mesh);
}

//------------------------------------------------------------------------------
// Query classes
//------------------------------------------------------------------------------

class DirectGWN2D
{
public:
  using CurveArrayType = axom::Array<NURBSCurve2D>;
  using NURBSCacheArray = axom::Array<CurveGWNCache>;

  DirectGWN2D() = default;

  void preprocess(const CurveArrayType& input_curves, bool use_memoization)
  {
    m_input_curves_view = input_curves.view();
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
    SLIC_INFO(axom::fmt::format("Direct query preprocessing (memoization caches): {} s",
                                timer.elapsedTimeInSec()));
  }

  void query(mfem::DataCollection& dc, const primal::WindingTolerances& tol)
  {
    auto* query_mesh = dc.GetMesh();
    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

    auto query_point = [&query_mesh](int idx) -> Point2D {
      Point2D pt;
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
          const Point2D q = query_point(static_cast<int>(nidx));
          double wn {};
          for(const auto& cache : m_input_curves_view)
          {
            wn += axom::primal::winding_number(q, cache, tol_copy.edge_tol, tol_copy.EPS);
          }
          winding[static_cast<int>(nidx)] = wn;
          inout[static_cast<int>(nidx)] = std::round(wn);
        });
      }
      else  // Use memoized form
      {
        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType nidx) {
          const Point2D q = query_point(static_cast<int>(nidx));
          double wn {};
          for(const auto& cache : m_nurbs_caches)
          {
            wn += axom::primal::winding_number(q, cache, tol_copy.edge_tol, tol_copy.EPS);
          }
          winding[static_cast<int>(nidx)] = wn;
          inout[static_cast<int>(nidx)] = std::round(wn);
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
  axom::ArrayView<const NURBSCurve2D> m_input_curves_view;
  NURBSCacheArray m_nurbs_caches;
};

class PolylineGWN2D
{
public:
  using BoxType = axom::primal::BoundingBox<double, 2>;
  using CurveArrayType = axom::Array<NURBSCurve2D>;
  using SegmentType = axom::primal::Segment<double, 2>;
  using GWNMoments = axom::quest::GWNMomentData<double, 2, 2>;

  PolylineGWN2D() = default;

  void preprocess(const CurveArrayType& input_curves, int segmentsPerKnotSpan)
  {
    axom::utilities::Timer timer(true);
    axom::utilities::Timer stage_timer(false);

    AXOM_ANNOTATE_SCOPE("preprocessing");

    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> poly_mesh(2, axom::mint::SEGMENT);

    stage_timer.start();
    {
      AXOM_ANNOTATE_SCOPE("linearization");
      axom::quest::LinearizeCurves lc;
      lc.getLinearMeshUniform(input_curves.view(), &poly_mesh, segmentsPerKnotSpan);
    }
    stage_timer.stop();

    const auto nlines = poly_mesh.getNumberOfCells();
    {
      SLIC_INFO(axom::fmt::format(
        axom::utilities::locale(),
        "Discretized {} curves with {} segments in each knot span for a total of {:L} segments.",
        input_curves.size(),
        segmentsPerKnotSpan,
        nlines));
      SLIC_INFO(axom::fmt::format("Preprocessing stage (linearization): {} s",
                                  stage_timer.elapsedTimeInSec()));
    }

    stage_timer.reset();
    stage_timer.start();
    {
      AXOM_ANNOTATE_SCOPE("extract_segments");

      m_segments.resize(nlines);
      auto segments_view = m_segments.view();

      axom::mint::for_all_cells<axom::SEQ_EXEC, axom::mint::xargs::coords>(
        &poly_mesh,
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

    stage_timer.reset();
    stage_timer.start();
    {
      AXOM_ANNOTATE_SCOPE("bvh_init");
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
    SLIC_INFO(axom::fmt::format("  Preprocessing stage (bvh): {} s", stage_timer.elapsedTimeInSec()));

    stage_timer.reset();
    stage_timer.start();
    {
      AXOM_ANNOTATE_SCOPE("moments");
      const auto triangles_view = m_segments.view();

      auto compute_moments = [triangles_view](std::int32_t currentNode,
                                              const std::int32_t* leafNodes) -> GWNMoments {
        const auto idx = leafNodes[currentNode];
        return GWNMoments(triangles_view[idx]);
      };

      m_internal_moments = m_bvh.template reduceTree<GWNMoments>(compute_moments);
    }
    stage_timer.stop();
    SLIC_INFO(
      axom::fmt::format("  Preprocessing stage (moments): {} s", stage_timer.elapsedTimeInSec()));

    timer.stop();
    SLIC_INFO(axom::fmt::format("Total preprocessing: {} s", timer.elapsedTimeInSec()));
  }

  void query(mfem::DataCollection& dc, const primal::WindingTolerances& tol, const double slice_z = 0.0)
  {
    if(m_segments.empty())
    {
      SLIC_WARNING("Skipping query; segment data is empty.");
      return;
    }

    const auto* query_mesh = dc.GetMesh();
    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

    auto query_point = [query_mesh](axom::IndexType idx) -> Point2D {
      Point2D pt({0., 0.});
      query_mesh->GetNode(static_cast<int>(idx), pt.data());
      return pt;
    };

    axom::utilities::Timer query_timer(true);
    {
      AXOM_ANNOTATE_SCOPE("query");

      const auto traverser = m_bvh.getTraverser();
      const auto segments_view = m_segments.view();
      const auto internal_moments_view = m_internal_moments.view();
      const primal::WindingTolerances tol_copy = tol;

      axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType index) {
        const double wn = axom::quest::fast_approximate_winding_number(query_point(index),
                                                                       traverser,
                                                                       segments_view,
                                                                       internal_moments_view,
                                                                       tol_copy);

        winding[static_cast<int>(index)] = wn;
        inout[static_cast<int>(index)] = std::round(wn);
      });
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
  axom::Array<SegmentType> m_segments;
  axom::Array<GWNMoments> m_internal_moments;
  axom::spin::BVH<2, axom::SEQ_EXEC> m_bvh;
};

using WNQueryType = std::variant<DirectGWN2D, PolylineGWN2D>;

// This framework will make more sense when there are more than two types available
WNQueryType make_wn_query(bool linearize_curves)
{
  if(linearize_curves)
  {
    return DirectGWN2D {};
  }

  return PolylineGWN2D {};
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger raii_logger;

  axom::CLI::App app {
    "Load mesh containing collection of curves"
    " and optionally generate a query mesh of winding numbers."};

  std::string inputFile;
  std::string outputPrefix = {"winding"};

  bool verbose {false};
  std::string annotationMode {"none"};
  bool memoized {true};
  bool vis {true};

  double segmentsPerKnotSpan {10};

  // Query mesh parameters
  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int queryOrder {1};

  primal::WindingTolerances tol;

  app.add_option("-i,--input", inputFile)
    ->description("MFEM mesh containing contours (1D segments)")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_option("-o,--output-prefix", outputPrefix)
    ->description(
      "Prefix for output 2D query mesh (in MFEM format) mesh containing "
      "winding number calculations")
    ->capture_default_str();

  app.add_flag("-v,--verbose", verbose, "verbose output")->capture_default_str();
  app.add_flag("--memoized,!--no-memoized", memoized, "Cache geometric data during query?")
    ->capture_default_str();
  app.add_flag("--vis,!--no-vis", vis, "Should we write out the results for visualization?")
    ->capture_default_str();

  // Options for query tolerances
  app.add_option("--edge-tol", tol.edge_tol)
    ->description("Relative edge tolerance for queries")
    ->check(axom::CLI::PositiveNumber)
    ->capture_default_str();
  app.add_option("--eps-tol", tol.EPS)
    ->description("Additional generic tolerance parameter")
    ->check(axom::CLI::PositiveNumber)
    ->capture_default_str();

#ifdef AXOM_USE_CALIPER
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif

  // Options for triangulation of the input STEP file
  auto* linearize_curves_subcommand = app.add_subcommand("linearize_curves")
                                        ->description("Options for linearizing NURBS curves")
                                        ->fallthrough();

  linearize_curves_subcommand->add_option("--num-segments", segmentsPerKnotSpan)
    ->description("Number of segments for each knot span of each input curve.")
    ->check(axom::CLI::PositiveNumber)
    ->capture_default_str();

  auto* query_mesh_subcommand =
    app.add_subcommand("query_mesh")->description("Options for setting up a query mesh")->fallthrough();
  auto* minbb = query_mesh_subcommand->add_option("--min", boxMins)
                  ->description("Min bounds for box mesh (x,y)")
                  ->expected(2);
  auto* maxbb = query_mesh_subcommand->add_option("--max", boxMaxs)
                  ->description("Max bounds for box mesh (x,y)")
                  ->expected(2);
  query_mesh_subcommand->add_option("--res", boxResolution)
    ->description("Resolution of the box mesh (i,j)")
    ->expected(2)
    ->required();
  query_mesh_subcommand->add_option("--order", queryOrder)
    ->description("polynomial order of the query mesh")
    ->check(axom::CLI::PositiveNumber);

  // add some requirements -- if user provides minbb or maxbb, we need both
  minbb->needs(maxbb);
  maxbb->needs(minbb);

  CLI11_PARSE(app, argc, argv);

  axom::utilities::raii::AnnotationsWrapper annotation_raii_wrapper(annotationMode);
  AXOM_ANNOTATE_SCOPE("winding number example");

  axom::Array<NURBSCurve2D> curves;
  {
    AXOM_ANNOTATE_SCOPE("read_mesh");

    axom::quest::MFEMReader mfem_reader;
    mfem_reader.setFileName(inputFile);

    const int ret = mfem_reader.read(curves);
    if(ret != axom::quest::MFEMReader::READ_SUCCESS)
    {
      return 1;
    }
  }

  // Early return if user didn't set up a query mesh
  if(boxResolution.empty())
  {
    return 0;
  }

  // Extract the curves and compute their bounding boxes along the way
  BoundingBox2D bbox;
  axom::Array<CurveGWNCache> memoized_curves;
  for(const auto& cur : curves)
  {
    bbox.addBox(cur.boundingBox());
  }
  SLIC_INFO(axom::fmt::format("Curve mesh bounding box: {}", bbox));

  // Query grid setup;
  // if user did not provide a bounding box, user input bounding box scaled by 10%
  mfem::DataCollection dc("winding_query");
  {
    // Create the desired winding number query instance
    auto wn_query = make_wn_query(app.got_subcommand("linearize_curves"));

    // Generate teh query grid and fields
    generate_query_mesh(dc, bbox, boxMins, boxMaxs, boxResolution, queryOrder);

    // Run the preprocess
    std::visit(
      [&](auto& wn) {
        using T = std::decay_t<decltype(wn)>;
        if constexpr(std::is_same_v<T, DirectGWN2D>)
        {
          wn.preprocess(curves, memoized);
        }
        else if constexpr(std::is_same_v<T, PolylineGWN2D>)
        {
          wn.preprocess(curves, segmentsPerKnotSpan);
        }
      },
      wn_query);

    // Run the query
    std::visit([&](auto& wn) { wn.query(dc, tol); }, wn_query);
  }

  // Save the query mesh and fields to disk using a format that can be viewed in VisIt
  if(vis)
  {
    AXOM_ANNOTATE_SCOPE("dump_mesh");

    mfem::VisItDataCollection windingDC(outputPrefix, dc.GetMesh());
    windingDC.RegisterField("winding", dc.GetField("winding"));
    windingDC.RegisterField("inout", dc.GetField("inout"));
    windingDC.Save();

    SLIC_INFO(axom::fmt::format("Outputting generated mesh '{}' to '{}'",
                                windingDC.GetCollectionName(),
                                axom::utilities::filesystem::getCWD()));
  }

  return 0;
}
