// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_winding_number_3d.cpp
 * \brief Example that computes the generalized winding number (GWN) of query points
 * against a collection of 3D trimmed NURBS patches loaded from a STEP file.
 *
 * The query points are taken from either:
 *  - a 3D Cartesian MFEM mesh, or
 *  - a 2D Cartesian MFEM mesh interpreted as a z-constant slice (default z=0).
 *
 * \note This example requires Axom to be configured with both MFEM and OpenCascade enabled.
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "mfem.hpp"

namespace primal = axom::primal;

using Point3D = primal::Point<double, 3>;
using BoundingBox2D = primal::BoundingBox<double, 2>;
using BoundingBox3D = primal::BoundingBox<double, 3>;

using NURBSPatch3D = axom::quest::STEPReader::NURBSPatch;
using PatchGWNCache = primal::detail::NURBSPatchGWNCache<double>;

struct WindingTolerances
{
  double edge_tol {1e-6};
  double ls_tol {1e-6};
  double quad_tol {1e-6};
  double disk_size {0.01};
  double EPS {1e-8};
};

void setup_mesh(mfem::DataCollection& dc, mfem::Mesh* query_mesh, int queryOrder)
{
  AXOM_ANNOTATE_SCOPE("setup_mesh");

  dc.SetOwnData(true);
  dc.SetMesh(query_mesh);

  const int dim = query_mesh->Dimension();

  // Create grid functions for the winding field; will take care of fes and fec memory via MakeOwner()
  auto* winding_fec = new mfem::H1_FECollection(queryOrder, dim);
  auto* winding_fes = new mfem::FiniteElementSpace(query_mesh, winding_fec, 1);
  mfem::GridFunction* winding = new mfem::GridFunction(winding_fes);
  winding->MakeOwner(winding_fec);

  // Create grid functions for the inout field; will take care of fes and fec memory via MakeOwner()
  auto* inout_fec = new mfem::H1_FECollection(queryOrder, dim);
  auto* inout_fes = new mfem::FiniteElementSpace(query_mesh, inout_fec, 1);
  mfem::GridFunction* inout = new mfem::GridFunction(inout_fes);
  inout->MakeOwner(inout_fec);

  dc.RegisterField("winding", winding);
  dc.RegisterField("inout", inout);
}

template <typename PatchArrayType>
void run_query(mfem::DataCollection& dc,
               const PatchArrayType& patches,
               const WindingTolerances& tol,
               double slice_z = 0.0)
{
  AXOM_ANNOTATE_SCOPE("run_query");

  auto* query_mesh = dc.GetMesh();
  auto& winding = *dc.GetField("winding");
  auto& inout = *dc.GetField("inout");

  const int space_dim = query_mesh->SpaceDimension();
  const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

  auto query_point = [&](int idx) -> Point3D {
    double coords[3] = {0., 0., 0.};
    query_mesh->GetNode(idx, coords);

    if(space_dim == 2)
    {
      return Point3D {coords[0], coords[1], slice_z};
    }

    return Point3D {coords[0], coords[1], coords[2]};
  };

  for(int nidx = 0; nidx < num_query_points; ++nidx)
  {
    const Point3D q = query_point(nidx);

    double wn {};
    for(const auto& patch : patches)
    {
      wn += axom::primal::winding_number(q,
                                         patch,
                                         tol.edge_tol,
                                         tol.ls_tol,
                                         tol.quad_tol,
                                         tol.disk_size,
                                         tol.EPS);
    }

    winding[nidx] = wn;
    inout[nidx] = std::round(wn);
  }
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger raii_logger;

  axom::CLI::App app {
    "Load a STEP file containing trimmed NURBS patches and optionally generate a "
    "query mesh of 3D generalized winding numbers."};

  std::string inputFile;
  std::string outputPrefix {"winding3d"};

  bool verbose {false};
  bool validate {false};
  std::string annotationMode {"none"};
  bool memoized {true};

  // Query mesh parameters (either 2D slice or full 3D grid)
  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int queryOrder {1};
  double sliceZ {0.0};

  // Winding-number tolerances
  WindingTolerances tol;

  app.add_option("-i,--input", inputFile)
    ->description("Input STEP file containing a trimmed NURBS BRep")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_option("-o,--output-prefix", outputPrefix)
    ->description(
      "Prefix for output query mesh (in MFEM/VisIt format) containing winding number results")
    ->capture_default_str();

  app.add_flag("-v,--verbose", verbose, "verbose output")->capture_default_str();
  app.add_flag("--validate", validate, "Run STEP model validation checks")->capture_default_str();
  app.add_flag("--memoized,!--no-memoized", memoized, "Cache geometric data during query?")
    ->capture_default_str();

  app.add_option("--edge-tol", tol.edge_tol)
    ->description("Tolerance for boundary proximity checks")
    ->check(axom::CLI::PositiveNumber)
    ->capture_default_str();
  app.add_option("--ls-tol", tol.ls_tol)
    ->description("Tolerance for line-surface intersection")
    ->check(axom::CLI::PositiveNumber)
    ->capture_default_str();
  app.add_option("--quad-tol", tol.quad_tol)
    ->description("Relative error tolerance for quadrature")
    ->check(axom::CLI::PositiveNumber)
    ->capture_default_str();
  app.add_option("--disk-size", tol.disk_size)
    ->description("Extracted disk size as fraction of parameter-space bbox diagonal")
    ->check(axom::CLI::PositiveNumber)
    ->capture_default_str();
  app.add_option("--eps", tol.EPS)
    ->description("Miscellaneous numerical tolerance")
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

  auto* query_mesh_subcommand =
    app.add_subcommand("query_mesh")->description("Options for setting up a query mesh")->fallthrough();

  query_mesh_subcommand->add_option("--min", boxMins)
    ->description("Min bounds for box mesh (x,y[,z])")
    ->expected(2, 3)
    ->required();
  query_mesh_subcommand->add_option("--max", boxMaxs)
    ->description("Max bounds for box mesh (x,y[,z])")
    ->expected(2, 3)
    ->required();
  query_mesh_subcommand->add_option("--res", boxResolution)
    ->description("Resolution of the box mesh (i,j[,k])")
    ->expected(2, 3)
    ->required();
  query_mesh_subcommand->add_option("--order", queryOrder)
    ->description("polynomial order of the query mesh")
    ->check(axom::CLI::PositiveNumber);
  query_mesh_subcommand->add_option("--slice-z", sliceZ)
    ->description("Z value for 2D slice query meshes (when --min/--max are 2D)")
    ->capture_default_str();

  CLI11_PARSE(app, argc, argv);

  axom::utilities::raii::AnnotationsWrapper annotation_raii_wrapper(annotationMode);
  AXOM_ANNOTATE_SCOPE("3D winding number example");

  axom::quest::STEPReader step_reader;
  step_reader.setFileName(inputFile);
  step_reader.setVerbosity(verbose);

  {
    AXOM_ANNOTATE_SCOPE("read_step");

    const int ret = step_reader.read(validate);
    if(ret != 0)
    {
      SLIC_ERROR(axom::fmt::format("Failed to read STEP file '{}'", inputFile));
      return 1;
    }
  }

  const auto& patches = step_reader.getPatchArray();
  SLIC_INFO(step_reader.getBRepStats());
  SLIC_INFO(axom::fmt::format("STEP file units: {}", step_reader.getFileUnits()));
  SLIC_INFO(axom::fmt::format("Loaded {} trimmed NURBS patches.", patches.size()));

  // Extract the patches and compute their bounding boxes along the way
  BoundingBox3D bbox;
  axom::Array<PatchGWNCache> memoized_patches;
  {
    AXOM_ANNOTATE_SCOPE("preprocessing");

    int count {0};
    for(const auto& patch : patches)
    {
      auto pbox = patch.boundingBox();
      bbox.addBox(pbox);

      SLIC_INFO_IF(verbose, axom::fmt::format("Patch {} bbox: {}", count++, pbox));

      if(memoized)
      {
        memoized_patches.emplace_back(PatchGWNCache(patch));
      }
    }
  }
  SLIC_INFO(axom::fmt::format("Patch collection bounding box: {}", bbox));

  // Early return if user didn't set up a query mesh
  if(boxMins.empty())
  {
    return 0;
  }

  if((boxMins.size() != boxMaxs.size()) || (boxMins.size() != boxResolution.size()))
  {
    SLIC_ERROR("--min/--max/--res must have the same number of entries (2 or 3).");
    return 1;
  }

  if(boxMins.size() != 2 && boxMins.size() != 3)
  {
    SLIC_ERROR("--min/--max/--res must have 2 entries (slice) or 3 entries (3D grid).");
    return 1;
  }

  mfem::DataCollection dc("winding_query");

  if(boxMins.size() == 2)
  {
    using Point2D = primal::Point<double, 2>;
    const auto query_res = axom::NumericArray<int, 2>(boxResolution.data());
    const auto query_box = BoundingBox2D(Point2D(boxMins.data()), Point2D(boxMaxs.data()));

    mfem::Mesh* query_mesh =
      axom::quest::util::make_cartesian_mfem_mesh_2D(query_box, query_res, queryOrder);

    setup_mesh(dc, query_mesh, queryOrder);

    if(memoized)
    {
      run_query(dc, memoized_patches, tol, sliceZ);
    }
    else
    {
      run_query(dc, patches, tol, sliceZ);
    }
  }
  else
  {
    using Point3 = primal::Point<double, 3>;
    const auto query_res = axom::NumericArray<int, 3>(boxResolution.data());
    const auto query_box = BoundingBox3D(Point3(boxMins.data()), Point3(boxMaxs.data()));

    mfem::Mesh* query_mesh =
      axom::quest::util::make_cartesian_mfem_mesh_3D(query_box, query_res, queryOrder);

    setup_mesh(dc, query_mesh, queryOrder);

    if(memoized)
    {
      run_query(dc, memoized_patches, tol);
    }
    else
    {
      run_query(dc, patches, tol);
    }
  }

  // Save the query mesh and fields to disk using a format that can be viewed in VisIt
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
