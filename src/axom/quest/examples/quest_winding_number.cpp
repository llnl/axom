// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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

/// Helper function to set up the mesh and associated winding and inout fields; uses an mfem::DataCollection to hold everything together
void setup_mesh(mfem::DataCollection& dc,
                const BoundingBox2D& query_box,
                const axom::NumericArray<int, 2>& query_res,
                int queryOrder)
{
  AXOM_ANNOTATE_SCOPE("setup_mesh");

  constexpr int DIM = 2;

  dc.SetOwnData(true);

  mfem::Mesh* query_mesh =
    axom::quest::util::make_cartesian_mfem_mesh_2D(query_box, query_res, queryOrder);
  dc.SetMesh(query_mesh);

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

template <typename CurveArray>
void run_query(mfem::DataCollection& dc, const CurveArray& curves)
{
  AXOM_ANNOTATE_SCOPE("run_query");

  auto* query_mesh = dc.GetMesh();
  auto& winding = *dc.GetField("winding");
  auto& inout = *dc.GetField("inout");

  // Utility function to get query point from query index
  const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();
  auto query_point = [&query_mesh](int idx) -> Point2D {
    Point2D pt;
    query_mesh->GetNode(idx, pt.data());
    return pt;
  };

  // Query the winding numbers at each degree of freedom (DoF) of the query mesh.
  // The loop below independently checks every curve for each query point.
  for(int nidx = 0; nidx < num_query_points; ++nidx)
  {
    const Point2D q = query_point(nidx);
    double wn {};
    for(const auto& c : curves)
    {
      wn += axom::primal::winding_number(q, c);
    }

    winding[nidx] = wn;
    inout[nidx] = std::round(wn);
  }
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

  // Query mesh parameters
  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int queryOrder {1};

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

#ifdef AXOM_USE_CALIPER
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif

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
  // Extract the curves and compute their bounding boxes along the way
  BoundingBox2D bbox;
  axom::Array<CurveGWNCache> memoized_curves;
  {
    AXOM_ANNOTATE_SCOPE("preprocessing");

    axom::utilities::Timer preproc_timer(true);
    int count {0};
    for(const auto& cur : curves)
    {
      SLIC_INFO_IF(verbose, axom::fmt::format("Element {}: {}", count++, cur));

      bbox.addBox(cur.boundingBox());

      // Add curves to GWN Cache objects that dynamically store intermediate
      //  curve subdivisions to be reused across query points
      if(memoized)
      {
        memoized_curves.emplace_back(CurveGWNCache(cur));
      }
    }
    preproc_timer.stop();

    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Preprocessing curves took {:.4Lf} seconds",
                                preproc_timer.elapsed()));
    AXOM_ANNOTATE_METADATA("preprocessing_time", preproc_timer.elapsed(), "");
  }
  SLIC_INFO(axom::fmt::format("Curve mesh bounding box: {}", bbox));

  // Early return if user didn't set up a query mesh
  if(boxResolution.empty())
  {
    return 0;
  }

  // Generate a Cartesian (high order) mesh for the query points
  const bool has_query_box = boxMins.size() > 0;
  const auto query_res = axom::NumericArray<int, 2>(boxResolution.data());
  const auto query_box = has_query_box
    ? BoundingBox2D(Point2D(boxMins.data()), Point2D(boxMaxs.data()))
    : BoundingBox2D(bbox.getMin(), bbox.getMax()).scale(1.05);

  mfem::DataCollection dc("winding_query");
  setup_mesh(dc, query_box, query_res, queryOrder);

  // Run the query (optionally, with memoization)
  {
    axom::utilities::Timer query_timer(true);
    if(memoized)
    {
      run_query(dc, memoized_curves);
    }
    else
    {
      run_query(dc, curves);
    }
    query_timer.stop();

    const int ndofs = dc.GetField("winding")->FESpace()->GetNDofs();
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Querying {:L} samples in winding number field took {:.3Lf} seconds"
                                " (@ {:.0Lf} queries per second; {:.5Lf} ms per query)",
                                ndofs,
                                query_timer.elapsed(),
                                ndofs / query_timer.elapsed(),
                                query_timer.elapsedTimeInMilliSec() / ndofs));
    AXOM_ANNOTATE_METADATA("query_points", ndofs, "");
    AXOM_ANNOTATE_METADATA("query_time", query_timer.elapsed(), "");
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
