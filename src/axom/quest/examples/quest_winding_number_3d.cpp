// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace primal = axom::primal;

using Point3D = primal::Point<double, 3>;
using BoundingBox2D = primal::BoundingBox<double, 2>;
using BoundingBox3D = primal::BoundingBox<double, 3>;

using NURBSPatch3D = axom::quest::STEPReader::NURBSPatch;
using PatchGWNCache = primal::detail::NURBSPatchGWNCache<double>;

struct WindingTolerances
{
  double ls_tol {1e-6};
  double quad_tol {1e-6};
  double disk_size {0.01};
  double edge_tol {1e-8};
  double EPS {1e-8};
};

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

/// Helper function to set up the mesh and associated winding and inout fields.
/// Uses an mfem::DataCollection to hold everything together.
void setup_mesh(mfem::DataCollection& dc, mfem::Mesh* query_mesh, int queryOrder)
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

template <typename PatchArrayType>
void run_query(mfem::DataCollection& dc,
               const PatchArrayType& patches,
               const WindingTolerances& tol,
               const double slice_z = 0.0)
{
  AXOM_ANNOTATE_SCOPE("run_query");

  auto* query_mesh = dc.GetMesh();
  auto& winding = *dc.GetField("winding");
  auto& inout = *dc.GetField("inout");

  const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

  auto query_point = [&](int idx) -> Point3D {
    Point3D pt({0., 0., slice_z});
    query_mesh->GetNode(idx, pt.data());
    return pt;
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

  std::string inputFile;
  std::string outputPrefix {"winding3d"};

  bool verbose {false};
  std::string annotationMode {"none"};
  bool memoized {true};
  bool vis {true};
  bool validate {false};
  bool stats {false};

  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int queryOrder {1};
  double sliceZ {0.0};
  WindingTolerances tol;

  axom::CLI::App app {
    "Load a STEP file containing trimmed NURBS patches "
    "and optionally generate a query grid of generalized winding numbers."};

  // Command line options and validation
  {
    app.add_option("-i,--input", inputFile)
      ->description("Input STEP file containing a trimmed NURBS BRep")
      ->required()
      ->check(axom::CLI::ExistingFile);

    app.add_option("-o,--output-prefix", outputPrefix)
      ->description(
        "Prefix for output query grid (in MFEM format) containing winding number results")
      ->capture_default_str();

    app.add_flag("-v,--verbose", verbose, "verbose output")->capture_default_str();
    app.add_flag("--validate", validate, "Run STEP model validation checks")->capture_default_str();
    app.add_flag("--memoized,!--no-memoized", memoized, "Cache geometric data during query?")
      ->capture_default_str();
    app.add_flag("--vis,!--no-vis", vis, "Should we write out the results for visualization?")
      ->capture_default_str();
    app.add_flag("--stats,!--no-stats", stats, "Compute summary stats for query fields?")
      ->capture_default_str();

    // Options for query tolerances; for now, only expose the line search and quadrature tolerances
    app.add_option("--ls-tol", tol.ls_tol)
      ->description("Tolerance for line-surface intersection")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();
    app.add_option("--quad-tol", tol.quad_tol)
      ->description("Relative error tolerance for quadrature")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();
    app.add_option("--disk-size", tol.disk_size)
      ->description("Relative disk size for winding number edge cases")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();
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

    auto* query_mesh_subcommand =
      app.add_subcommand("query_mesh")->description("Options for setting up a query mesh")->fallthrough();

    auto* minbb = query_mesh_subcommand->add_option("--min", boxMins)
                    ->description("Min bounds for box mesh (x,y[,z])")
                    ->expected(2, 3);
    auto* maxbb = query_mesh_subcommand->add_option("--max", boxMaxs)
                    ->description("Max bounds for box mesh (x,y[,z])")
                    ->expected(2, 3);
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

    // add some requirements -- if user provides minbb or maxbb, we need both
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    // let's also check that they're consistently sized w/ each other and with the resolution
    query_mesh_subcommand->callback([&]() {
      if(const bool have_box = (minbb->count() > 0 || maxbb->count() > 0); have_box)
      {
        if(boxMins.size() != boxMaxs.size())
        {
          throw axom::CLI::ValidationError(
            "--min/--max",
            axom::fmt::format("must have the same number of values (2 for 2D or 3 for 3D). "
                              "Got --min={}, --max={}",
                              boxMins.size(),
                              boxMaxs.size()));
        }

        for(size_t d = 0; d < boxMins.size(); ++d)
        {
          if(boxMins[d] >= boxMaxs[d])
          {
            throw axom::CLI::ValidationError(
              "--min/--max",
              axom::fmt::format(
                "must satisfy min < max in every dimension; failed at index {}, mins: {}, maxs: {}",
                d,
                boxMins,
                boxMaxs));
          }
        }

        if(boxResolution.size() != boxMins.size())
        {
          throw axom::CLI::ValidationError(
            "--res",
            axom::fmt::format(
              "must have the same number of values as --min/--max. Got --res={}, --min/--max={}",
              boxResolution.size(),
              boxMins.size()));
        }
      }
    });
  }

  CLI11_PARSE(app, argc, argv);

  axom::utilities::raii::AnnotationsWrapper annotation_raii_wrapper(annotationMode);
  AXOM_ANNOTATE_SCOPE("3D winding number example");

  axom::utilities::Timer step_read_timer(true);
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
  step_read_timer.stop();
  SLIC_INFO(step_reader.getBRepStats());
  SLIC_INFO(axom::fmt::format("STEP file units: {}", step_reader.getFileUnits()));
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Loaded {} trimmed NURBS patches in {:.3Lf} seconds",
                              patches.size(),
                              step_read_timer.elapsed()));

  // Early return if user didn't set up a query mesh
  if(boxResolution.empty())
  {
    return 0;
  }

  // Preprocessing: Extract the patches and compute their bounding boxes along the way
  BoundingBox3D bbox;
  axom::Array<PatchGWNCache> memoized_patches(0, memoized ? patches.size() : 0);
  {
    AXOM_ANNOTATE_SCOPE("preprocessing");

    axom::utilities::Timer preproc_timer(true);
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
    preproc_timer.stop();

    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Preprocessing patches took {:.3Lf} seconds",
                                preproc_timer.elapsed()));
    AXOM_ANNOTATE_METADATA("preprocessing_time", preproc_timer.elapsed(), "");
  }
  SLIC_INFO(axom::fmt::format("Patch collection bounding box: {}", bbox));

  // Query grid setup; has some dimension-specific types;
  // if user did not provide a bounding box, user input bounding box scaled by 10%
  mfem::DataCollection dc("winding_query");

  const int query_dim = boxResolution.size();
  const bool has_query_box = boxMins.size() > 0;
  constexpr double scale_factor = 1.1;
  if(query_dim == 2)
  {
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
      axom::quest::util::make_cartesian_mfem_mesh_2D(query_box, query_res, queryOrder);

    setup_mesh(dc, query_mesh, queryOrder);
    AXOM_ANNOTATE_METADATA("query_dimension", 2, "");
  }
  else
  {
    using Point3D = primal::Point<double, 3>;
    const auto query_res = axom::NumericArray<int, 3>(boxResolution.data());
    const auto query_box = has_query_box
      ? BoundingBox3D(Point3D(boxMins.data()), Point3D(boxMaxs.data()))
      : BoundingBox3D(bbox.getMin(), bbox.getMax()).scale(scale_factor);

    SLIC_INFO(
      axom::fmt::format("Query grid resolution {} within bounding box {}", query_res, query_box));

    mfem::Mesh* query_mesh =
      axom::quest::util::make_cartesian_mfem_mesh_3D(query_box, query_res, queryOrder);

    setup_mesh(dc, query_mesh, queryOrder);
    AXOM_ANNOTATE_METADATA("query_dimension", 3, "");
  }

  // Run the query
  double query_time_s = 0.0;
  int ndofs = 0;
  double ms_per_query = 0.0;
  {
    axom::utilities::Timer query_timer(true);
    if(memoized)
    {
      run_query(dc, memoized_patches, tol, sliceZ);
    }
    else
    {
      run_query(dc, patches, tol, sliceZ);
    }
    query_timer.stop();

    ndofs = dc.GetField("winding")->FESpace()->GetNDofs();
    query_time_s = query_timer.elapsed();
    ms_per_query = query_timer.elapsedTimeInMilliSec() / ndofs;
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Querying {:L} samples in winding number field took {:.3Lf} seconds"
                                " (@ {:.0Lf} queries per second; {:.4Lf} ms per query)",
                                ndofs,
                                query_time_s,
                                ndofs / query_time_s,
                                ms_per_query));
    AXOM_ANNOTATE_METADATA("query_points", ndofs, "");
    AXOM_ANNOTATE_METADATA("query_time", query_time_s, "");
  }

  // Postprocess query results: norms, ranges, and integral statistics
  if(stats)
  {
    AXOM_ANNOTATE_SCOPE("postprocess");

    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto winding_stats = compute_field_stats(winding);
    const auto inout_stats = compute_field_stats(inout);
    const auto inout_integrals = compute_integrals(inout);

    std::int64_t pos_inout_dofs = 0;
    std::int64_t neg_inout_dofs = 0;
    for(int i = 0; i < inout.Size(); ++i)
    {
      const double v = inout[i];
      if(v > 0.0)
      {
        ++pos_inout_dofs;
      }
      else if(v < 0.0)
      {
        ++neg_inout_dofs;
      }
    }

    SLIC_INFO(
      axom::fmt::format("WN_STATS: dof_l2={:.6e} dof_linf={:.6e} l2={:.6e} min={:.6e} max={:.6e}",
                        winding_stats.dof_l2,
                        winding_stats.dof_linf,
                        winding_stats.l2,
                        winding_stats.min,
                        winding_stats.max));

    SLIC_INFO(axom::fmt::format(
      "INOUT_STATS: dof_l2={:.6e} dof_linf={:.6e} l2={:.6e} min={:.6e} max={:.6e} volume={:.6e} "
      "domain_volume={:.6e} vol_frac={:.6e} pos_dofs={} neg_dofs={}",
      inout_stats.dof_l2,
      inout_stats.dof_linf,
      inout_stats.l2,
      inout_stats.min,
      inout_stats.max,
      inout_integrals.integral,
      inout_integrals.domain_volume,
      (inout_integrals.domain_volume > 0.0 ? inout_integrals.integral / inout_integrals.domain_volume
                                           : 0.0),
      pos_inout_dofs,
      neg_inout_dofs));
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
