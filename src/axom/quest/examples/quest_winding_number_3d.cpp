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

#include "axom/core/utilities/StringUtilities.hpp"

#include "mfem.hpp"

#include <vector>

namespace primal = axom::primal;
namespace quest = axom::quest;

using Point3D = primal::Point<double, 3>;
using BoundingBox3D = primal::BoundingBox<double, 3>;

using NURBSPatch3D = quest::STEPReader::NURBSPatch;
using Triangle3D = primal::Triangle<double, 3>;

using RuntimePolicy = axom::runtime_policy::Policy;

//------------------------------------------------------------------------------
// CLI input
//------------------------------------------------------------------------------

class Input
{
public:
  std::string inputFile;
  std::string outputPrefix {"winding3d"};

  bool verbose {false};
  std::string annotationMode {"none"};
  bool memoized {true};
  bool vis {true};
  bool validate {false};
  bool stats {false};

  axom::runtime_policy::Policy policy = RuntimePolicy::seq;

  const std::array<std::string, 2> valid_algorithms {"direct", "fast-approximation"};
  std::string algorithm {valid_algorithms[1]};  // fast-approximation

  bool triangulate {false};
  double linear_deflection {0.1};
  double angular_deflection {0.5};
  bool deflection_is_relative {false};
  int approximation_order {2};

  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int queryOrder {1};
  double sliceZ {0.0};

  primal::WindingTolerances tol;

  void parse(int argc, char** argv, axom::CLI::App& app)
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

    // Options for triangulation of the input STEP file
    auto* triangulate_step_subcommand = app.add_subcommand("triangulate_step")
                                          ->description("Options for triangulating NURBS surfaces")
                                          ->fallthrough();

    triangulate_step_subcommand->add_option("--linear-deflection", linear_deflection)
      ->description(
        "Maximum allowed deviation between the original geometry and the triangulation.")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();
    triangulate_step_subcommand->add_option("--angular-deflection", angular_deflection)
      ->description(
        "Maximum allowed angular deviation (in radians) between normals of adjacent triangles.")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();
    triangulate_step_subcommand
      ->add_flag(
        "--is-relative",
        deflection_is_relative,
        "Is linear deflection in relative to local edge lengths (true) or mesh units (false)")
      ->capture_default_str();

    triangulate_step_subcommand->add_option("--algorithm", algorithm)
      ->description(
        "Use direct evaluation instead of fast, heirarchical approximation? (significantly "
        "slower, slightly more precise)")
      ->capture_default_str()
      ->check(axom::CLI::IsMember(valid_algorithms));
    triangulate_step_subcommand
      ->add_option("--expansion-order",
                   approximation_order,
                   "The order of the Taylor expansion (lower is faster, less precise)")
      ->expected(0, 2)
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
    std::stringstream pol_sstr;
    pol_sstr << "Set MIR runtime policy method.";
    pol_sstr << "\nSet to 'seq' or 0 to use the RAJA sequential policy.";
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
    pol_sstr << "\nSet to 'omp' or 1 to use the RAJA OpenMP policy.";
#endif

    app.add_option("-p, --policy", policy, pol_sstr.str())
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

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

    app.parse(argc, argv);

    triangulate = app.got_subcommand("triangulate_step");
  }
};

using GWNQueryType = std::variant<axom::quest::DirectGWN3D,
                                  axom::quest::TriangleGWN3D<axom::SEQ_EXEC, 0>,
                                  axom::quest::TriangleGWN3D<axom::SEQ_EXEC, 1>,
                                  axom::quest::TriangleGWN3D<axom::SEQ_EXEC, 2>,
                                  axom::quest::TriangleGWN3D<axom::OMP_EXEC, 0>,
                                  axom::quest::TriangleGWN3D<axom::OMP_EXEC, 1>,
                                  axom::quest::TriangleGWN3D<axom::OMP_EXEC, 2>>;

template <typename ExecSpace>
GWNQueryType pick_gwn_method(bool triangulate, int approximation_order)
{
  if(triangulate)
  {
    if(approximation_order == 0)
    {
      return axom::quest::TriangleGWN3D<ExecSpace, 0> {};
    }
    else if(approximation_order == 1)
    {
      return axom::quest::TriangleGWN3D<ExecSpace, 1> {};
    }
    else  // approximation_order == 2
    {
      return axom::quest::TriangleGWN3D<ExecSpace, 2> {};
    }
  }

  return axom::quest::DirectGWN3D {};
}

GWNQueryType make_gwn_query(axom::runtime_policy::Policy policy,
                            bool triangulate,
                            int approximation_order)
{
  if(policy == RuntimePolicy::omp)
  {
    SLIC_INFO(axom::fmt::format("Using policy omp with {} threads", omp_get_max_threads()));
    return pick_gwn_method<axom::OMP_EXEC>(triangulate, approximation_order);
  }

  SLIC_INFO("Using policy seq");
  return pick_gwn_method<axom::SEQ_EXEC>(triangulate, approximation_order);
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger raii_logger;

  // Parse command line arguments into input
  Input input;
  axom::CLI::App app {
    "Load a STEP file containing trimmed NURBS patches "
    "and optionally generate a query grid of generalized winding numbers."};

  try
  {
    input.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  axom::utilities::raii::AnnotationsWrapper annotation_raii_wrapper(input.annotationMode);
  AXOM_ANNOTATE_SCOPE("3D winding number example");

  // Bounding box for input shape
  BoundingBox3D shape_bbox;

  // Declare possible geometry input types
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tri_mesh(3, axom::mint::TRIANGLE);
  axom::Array<NURBSPatch3D> patches;

  if(axom::utilities::string::endsWith(input.inputFile, ".stl"))
  {
    AXOM_ANNOTATE_SCOPE("read_stl");

    axom::quest::STLReader stl_reader;
    stl_reader.setFileName(input.inputFile);

    axom::utilities::Timer read_timer(true);
    const int ret = stl_reader.read();

    if(ret != 0)
    {
      SLIC_ERROR(axom::fmt::format("Failed to read STL file '{}'", input.inputFile));
      return 1;
    }

    stl_reader.getMesh(&tri_mesh);
    read_timer.stop();

    BoundingBox3D* shape_bbox_ptr = &shape_bbox;
    axom::mint::for_all_nodes<axom::SEQ_EXEC, axom::mint::xargs::xyz>(
      &tri_mesh,
      AXOM_LAMBDA(axom::IndexType, double x, double y, double z) {
        shape_bbox_ptr->addPoint(Point3D {x, y, z});
      });

    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Loaded {} triangles in {:.3Lf} seconds",
                                stl_reader.getNumFaces(),
                                read_timer.elapsed()));
  }
  else if(axom::utilities::string::endsWith(input.inputFile, ".step"))
  {
    AXOM_ANNOTATE_SCOPE("read_step");

    axom::quest::STEPReader step_reader;
    step_reader.setFileName(input.inputFile);
    step_reader.setVerbosity(input.verbose);

    axom::utilities::Timer read_timer(true);
    const int ret = step_reader.read(input.validate);
    if(ret != 0)
    {
      SLIC_ERROR(axom::fmt::format("Failed to read STEP file '{}'", input.inputFile));
      return 1;
    }
    read_timer.stop();

    shape_bbox = step_reader.getBRepBoundingBox();

    int num_trimming_curves = 0;
    for(const auto& patch : step_reader.getPatchArray())
    {
      num_trimming_curves += patch.getNumTrimmingCurves();
    }

    SLIC_INFO(step_reader.getBRepStats());
    SLIC_INFO(axom::fmt::format("STEP file units: {}", step_reader.getFileUnits()));
    SLIC_INFO(axom::fmt::format(
      axom::utilities::locale(),
      "Loaded {} trimmed NURBS patches (with {} trimming curves) in {:.3Lf} seconds",
      patches.size(),
      num_trimming_curves,
      read_timer.elapsed()));

    if(input.triangulate)
    {
      read_timer.reset();
      read_timer.start();
      AXOM_ANNOTATE_SCOPE("triangulation");
      const int tc = step_reader.getTriangleMesh(&tri_mesh,
                                                 input.linear_deflection,
                                                 input.angular_deflection,
                                                 input.deflection_is_relative,
                                                 /* trimmed */ true);
      if(tc != 0)
      {
        SLIC_ERROR("Failed to triangulate STEP geometry.");
        return 1;
      }
      read_timer.stop();

      SLIC_INFO(
        axom::fmt::format(axom::utilities::locale(),
                          "Triangulated geometry with deflection {} and angular deflection {}"
                          " containing {:L} triangles in {:.3Lf} seconds",
                          input.linear_deflection,
                          input.angular_deflection,
                          tri_mesh.getNumberOfCells(),
                          read_timer.elapsed()));
    }
    else
    {
      patches = step_reader.getPatchArray();
    }
  }
  else
  {
    SLIC_WARNING(axom::fmt::format("Unsupported file type for input {}", input.inputFile));
  }

  // Early return if user didn't set up a query mesh
  if(input.boxResolution.empty())
  {
    return 0;
  }

  // Query grid setup; has some dimension-specific types;
  // if user did not provide a bounding box, user input bounding box scaled by 10%
  mfem::DataCollection dc("winding_query");
  {
    // Create the desired winding number query instance
    auto wn_query = make_gwn_query(input.policy, input.triangulate, input.approximation_order);

    // Generate the query grid and fields
    quest::generate_gwn_query_mesh(dc,
                                   shape_bbox,
                                   input.boxMins,
                                   input.boxMaxs,
                                   input.boxResolution,
                                   input.queryOrder);

    // Run the preprocess
    std::visit(
      [&](auto& wn) {
        using T = std::decay_t<decltype(wn)>;
        if constexpr(quest::gwn_input_type_v<T> == quest::GWNInputType::Surface)
        {
          wn.preprocess(patches, input.memoized);
        }
        else if constexpr(quest::gwn_input_type_v<T> == quest::GWNInputType::Triangulation)
        {
          wn.preprocess(&tri_mesh, input.algorithm == "direct");
        }
      },
      wn_query);

    // Run the query
    std::visit([&](auto& wn) { wn.query(dc, input.tol, input.sliceZ); }, wn_query);
  }

  // Postprocess query results: norms, ranges, and integral statistics
  if(input.stats)
  {
    AXOM_ANNOTATE_SCOPE("postprocess");

    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto winding_stats = axom::quest::compute_field_stats(winding);
    const auto inout_stats = axom::quest::compute_field_stats(inout);
    const auto inout_integrals = axom::quest::compute_integrals(inout);

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
  if(input.vis)
  {
    AXOM_ANNOTATE_SCOPE("dump_mesh");

    mfem::VisItDataCollection windingDC(input.outputPrefix, dc.GetMesh());
    windingDC.RegisterField("winding", dc.GetField("winding"));
    windingDC.RegisterField("inout", dc.GetField("inout"));
    windingDC.Save();

    SLIC_INFO(axom::fmt::format("Outputting generated mesh '{}' to '{}'",
                                windingDC.GetCollectionName(),
                                axom::utilities::filesystem::getCWD()));
  }

  return 0;
}
