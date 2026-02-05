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

using Triangle3D = primal::Triangle<double, 3>;

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
///  if user did not provide a bounding box, user input bounding box scaled by 10%
void generate_query_mesh(mfem::DataCollection& dc,
                         const BoundingBox3D& bbox,
                         const std::vector<double>& boxMins,
                         const std::vector<double>& boxMaxs,
                         const std::vector<int>& boxResolution,
                         int queryOrder)
{
  AXOM_ANNOTATE_SCOPE("generate_query_mesh");

  const int query_dim = static_cast<int>(boxResolution.size());
  const bool has_query_box = !boxMins.empty();

  SLIC_INFO(axom::fmt::format("Patch collection bounding box: {}", bbox));

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
  }
}

//------------------------------------------------------------------------------
// Query classes
//------------------------------------------------------------------------------

class DirectGWN3D
{
public:
  using PatchArrayType = axom::Array<NURBSPatch3D>;
  using NURBSCacheArray = axom::Array<PatchGWNCache>;

  DirectGWN3D() = default;

  void preprocess(const PatchArrayType& input_patches, bool use_memoization)
  {
    m_input_patches_view = input_patches.view();
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
    SLIC_INFO(axom::fmt::format("Direct query preprocessing (memoization caches): {} s",
                                timer.elapsedTimeInSec()));
  }

  void query(mfem::DataCollection& dc,
             const primal::WindingTolerances& tol,
             const double slice_z = 0.0) const
  {
    auto* query_mesh = dc.GetMesh();
    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

    auto query_point = [query_mesh, slice_z](int idx) -> Point3D {
      Point3D pt({0., 0., slice_z});
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
          const Point3D q = query_point(static_cast<int>(nidx));
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
          inout[static_cast<int>(nidx)] = std::round(wn);
        });
      }
      else  // Use memoized form
      {
        const auto nurbs_patches_view = m_nurbs_caches.view();
        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType nidx) {
          const Point3D q = query_point(static_cast<int>(nidx));
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
  axom::ArrayView<const NURBSPatch3D> m_input_patches_view;
  NURBSCacheArray m_nurbs_caches;
};

class TriangleGWN3D
{
public:
  using BoxType = axom::primal::BoundingBox<double, 3>;
  using TriangleType = axom::primal::Triangle<double, 3>;
  using GWNMoments = axom::quest::GWNMomentData<double, 3, 2>;

  TriangleGWN3D() = default;

  void preprocess(axom::quest::STEPReader& step_reader,
                  double linear_deflection,
                  double angular_deflection,
                  bool is_relative,
                  bool useDirectEval)
  {
    axom::utilities::Timer timer(true);
    axom::utilities::Timer stage_timer(false);

    AXOM_ANNOTATE_SCOPE("preprocessing");

    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tri_mesh(3, axom::mint::TRIANGLE);
    stage_timer.start();
    {
      AXOM_ANNOTATE_SCOPE("triangulation");
      const int rc =
        step_reader.getTriangleMesh(&tri_mesh, linear_deflection, angular_deflection, is_relative, true);
      if(rc != 0)
      {
        SLIC_ERROR("Failed to triangulate STEP geometry.");
        return;
      }
    }
    stage_timer.stop();

    const auto ntris = tri_mesh.getNumberOfCells();
    {
      SLIC_INFO(
        axom::fmt::format(axom::utilities::locale(),
                          "Triangulated geometry with deflection {} and angular deflection {}"
                          " containing {:L} triangles",
                          linear_deflection,
                          angular_deflection,
                          ntris));
      SLIC_INFO(axom::fmt::format("  Preprocessing stage (triangulation): {} s",
                                  stage_timer.elapsedTimeInSec()));
    }

    if(ntris <= 0)
    {
      SLIC_WARNING("Triangle mesh contains no cells; skipping preprocessing.");
      return;
    }

    stage_timer.reset();
    stage_timer.start();
    {
      AXOM_ANNOTATE_SCOPE("extract_triangles");

      m_triangles.resize(ntris);
      auto triangles_view = m_triangles.view();

      axom::mint::for_all_cells<axom::SEQ_EXEC, axom::mint::xargs::coords>(
        &tri_mesh,
        AXOM_LAMBDA(axom::IndexType cellIdx,
                    const axom::numerics::Matrix<double>& coords,
                    [[maybe_unused]] const axom::IndexType* nodeIds) {
          triangles_view[cellIdx] =
            TriangleType {Point3D {coords(0, 0), coords(1, 0), coords(2, 0)},
                          Point3D {coords(0, 1), coords(1, 1), coords(2, 1)},
                          Point3D {coords(0, 2), coords(1, 2), coords(2, 2)}};
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

        m_internal_moments = m_bvh.template reduceTree<GWNMoments>(compute_moments);
      }
      stage_timer.stop();
      SLIC_INFO(
        axom::fmt::format("  Preprocessing stage (moments): {} s", stage_timer.elapsedTimeInSec()));
    }

    timer.stop();
    SLIC_INFO(axom::fmt::format("Total preprocessing: {} s", timer.elapsedTimeInSec()));
  }

  void query(mfem::DataCollection& dc, const primal::WindingTolerances& tol, const double slice_z = 0.0)
  {
    if(m_triangles.empty())
    {
      SLIC_WARNING("Skipping query; triangle data is empty.");
      return;
    }

    const auto* query_mesh = dc.GetMesh();
    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

    auto query_point = [query_mesh, slice_z](axom::IndexType idx) -> Point3D {
      Point3D pt({0., 0., slice_z});
      query_mesh->GetNode(static_cast<int>(idx), pt.data());
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
          const double wn = axom::quest::fast_approximate_winding_number(query_point(index),
                                                                         traverser,
                                                                         triangles_view,
                                                                         internal_moments_view,
                                                                         tol_copy);

          winding[static_cast<int>(index)] = wn;
          inout[static_cast<int>(index)] = std::round(wn);
        });
      }
      // Use direct formula
      else
      {
        axom::for_all<axom::SEQ_EXEC>(num_query_points, [=, &winding, &inout](axom::IndexType index) {
          const Point3D q = query_point(static_cast<int>(index));
          double wn {};
          for(const auto& tri : m_triangles)
          {
            wn += axom::primal::winding_number(q, tri, tol_copy.edge_tol, tol_copy.EPS);
          }

          winding[static_cast<int>(index)] = wn;
          inout[static_cast<int>(index)] = std::round(wn);
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
};

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

  bool triangulate {false};
  double linear_deflection {0.1};
  double angular_deflection {0.5};
  bool deflection_is_relative {false};
  bool directTriangleEval {false};

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

    triangulate_step_subcommand
      ->add_flag("--direct,!--fast-approximation",
                 directTriangleEval,
                 "Use direct evaluation instead of fast, heirarchical approximation? "
                 "(significantly slower, slightly more precise)")
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

    app.parse(argc, argv);

    triangulate = app.got_subcommand("triangulate_step");
  }
};

using WNQueryType = std::variant<DirectGWN3D, TriangleGWN3D>;

// This framework will make more sense when there are more than two types available
WNQueryType make_wn_query(const Input& input)
{
  if(input.triangulate)
  {
    return TriangleGWN3D {};
  }

  return DirectGWN3D {};
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

  // Load the Step file
  axom::utilities::Timer step_read_timer(true);
  axom::quest::STEPReader step_reader;
  step_reader.setFileName(input.inputFile);
  step_reader.setVerbosity(input.verbose);

  {
    AXOM_ANNOTATE_SCOPE("read_step");

    const int ret = step_reader.read(input.validate);
    if(ret != 0)
    {
      SLIC_ERROR(axom::fmt::format("Failed to read STEP file '{}'", input.inputFile));
      return 1;
    }
  }
  step_read_timer.stop();

  // Read the patches/triangles and, if a query mesh is specified, get the bbox from the reader
  const auto& patches = step_reader.getPatchArray();
  step_read_timer.stop();
  {
    int num_trimming_curves = 0;
    for(const auto& patch : patches)
    {
      num_trimming_curves += patch.getNumTrimmingCurves();
    }
    AXOM_ANNOTATE_METADATA("num_original_trimming_curves", num_trimming_curves, "");
    AXOM_ANNOTATE_METADATA("num_original_patches", patches.size(), "");
    AXOM_ANNOTATE_METADATA("units", step_reader.getFileUnits(), "");

    SLIC_INFO(step_reader.getBRepStats());
    SLIC_INFO(axom::fmt::format("STEP file units: {}", step_reader.getFileUnits()));
    SLIC_INFO(axom::fmt::format(
      axom::utilities::locale(),
      "Loaded {} trimmed NURBS patches (with {} trimming curves) in {:.3Lf} seconds",
      patches.size(),
      num_trimming_curves,
      step_read_timer.elapsed()));
  }

  // Early return if user didn't set up a query mesh
  if(input.boxResolution.empty())
  {
    return 0;
  }

  BoundingBox3D bbox = step_reader.getBRepBoundingBox();

  // Query grid setup; has some dimension-specific types;
  // if user did not provide a bounding box, user input bounding box scaled by 10%
  mfem::DataCollection dc("winding_query");
  {
    // Create the desired winding number query instance
    auto wn_query = make_wn_query(input);

    // Generate the query grid and fields
    generate_query_mesh(dc, bbox, input.boxMins, input.boxMaxs, input.boxResolution, input.queryOrder);

    // Run the preprocess
    std::visit(
      [&](auto& wn) {
        using T = std::decay_t<decltype(wn)>;
        if constexpr(std::is_same_v<T, DirectGWN3D>)
        {
          wn.preprocess(patches, input.memoized);
        }
        else if constexpr(std::is_same_v<T, TriangleGWN3D>)
        {
          wn.preprocess(step_reader,
                        input.linear_deflection,
                        input.angular_deflection,
                        input.deflection_is_relative,
                        input.directTriangleEval);
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
