// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_gwn_benchmark.cpp
 * \brief Whole-grid GWN benchmarks for 2D contours and 3D STEP geometry.
 *
 * Build with ENABLE_BENCHMARKS, MFEM, OpenCascade, RAJA, and OpenMP enabled,
 * then run, for example:
 *
 *   OMP_NUM_THREADS=8 ./build/tests/quest_gwn_benchmark \
 *     --benchmark_min_time=0.1s --benchmark_repetitions=3 \
 *     --benchmark_display_aggregates_only=true --benchmark_time_unit=ms
 *
 * The case matrix is defined in CASES below. --benchmark_filter selects from
 * those cases; filter on 2d/curved, 2d/discrete, 3d/curved, or 3d/discrete,
 * and add /query/ or /preprocess/ to narrow the selection. --policy selects
 * seq, omp (the default), or both. OMP_NUM_THREADS controls the OMP runs.
 * Query cases report points_per_second using wall time and query_points.
 * 
 * Preprocess cases time subdivision or primitive extraction, BVH and moment
 * construction, and cache setup. Preprocess cases do not construct or query
 * the grid. Input reading, curve linearization or STEP triangulation, and
 * query-grid construction are outside both timed phases.
 *
 * --mesh-file and --step-file select the 2D and 3D inputs. Both files are
 * loaded before benchmark filtering, so both are required when AXOM_DATA_DIR
 * is unavailable, even if the filter selects only one dimension. 
 * 
 * The default grids contain 315^2 and 45^3 cells, producing 99,856 and 97,336
 * query points; --grid-resolution-2d and --grid-resolution-3d change these
 * values. 
 * 
 * Curved cases use memoization and order-2 fast approximation.
 */

#include "axom/config.hpp"

#include "axom/CLI11.hpp"
#include "axom/core.hpp"
#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/fmt.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/quest/GWNMethods.hpp"
#include "axom/quest/LinearizeCurves.hpp"
#include "axom/quest/io/MFEMReader.hpp"
#include "axom/quest/io/STEPReader.hpp"
#include "axom/slic.hpp"

#include "benchmark/benchmark.h"
#include "mfem.hpp"

#include <array>
#include <cmath>
#include <string>
#include <type_traits>
#include <vector>

namespace quest = axom::quest;
namespace spin = axom::spin;
namespace slic = axom::slic;

namespace
{

using MintMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
using Curve = axom::primal::NURBSCurve<double, 2>;

struct Problem3D
{
  quest::STEPReader reader;
  MintMesh triangles {3, axom::mint::TRIANGLE};
  axom::primal::BoundingBox<double, 3> bounds;
  int grid_resolution {45};

  bool load(const std::string& step_file)
  {
    reader.setFileName(step_file);
    if(reader.read(false) != 0 || reader.getPatchArray().empty())
    {
      return false;
    }
    bounds = reader.getBRepBoundingBox();
    return reader.getTriangleMesh(&triangles, 0.1, 0.5, false, true) == 0 &&
      triangles.getNumberOfCells() > 0;
  }
};

struct Problem2D
{
  axom::Array<Curve> curves;
  MintMesh segments {2, axom::mint::SEGMENT};
  axom::primal::BoundingBox<double, 2> bounds;
  int grid_resolution {315};

  bool load(const std::string& mesh_file)
  {
    quest::MFEMReader reader;
    reader.setFileName(mesh_file);
    if(reader.read(curves) != quest::MFEMReader::READ_SUCCESS || curves.empty())
    {
      return false;
    }
    for(const auto& curve : curves)
    {
      bounds.addBox(curve.boundingBox());
    }
    quest::LinearizeCurves linearizer;
    linearizer.getLinearMeshUniform(curves.view(), &segments, 10);
    return segments.getNumberOfCells() > 0;
  }
};

struct Case
{
  const char* name;
  axom::IndexType bucket_size;
  spin::SpaceFillingCurve ordering;
  spin::BucketedBVHGrouping grouping;
};

constexpr std::array<Case, 10> CASES {{
  {"morton/b1/none", 1, spin::SpaceFillingCurve::Morton, spin::BucketedBVHGrouping::Sequential},
  {"hilbert/b1/none", 1, spin::SpaceFillingCurve::Hilbert, spin::BucketedBVHGrouping::Sequential},
  {"morton/b2/sequential", 2, spin::SpaceFillingCurve::Morton, spin::BucketedBVHGrouping::Sequential},
  {"morton/b2/shared_prefix", 2, spin::SpaceFillingCurve::Morton, spin::BucketedBVHGrouping::SharedPrefix},
  {"hilbert/b2/sequential", 2, spin::SpaceFillingCurve::Hilbert, spin::BucketedBVHGrouping::Sequential},
  {"hilbert/b2/shared_prefix", 2, spin::SpaceFillingCurve::Hilbert, spin::BucketedBVHGrouping::SharedPrefix},
  {"morton/b8/sequential", 8, spin::SpaceFillingCurve::Morton, spin::BucketedBVHGrouping::Sequential},
  {"morton/b8/shared_prefix", 8, spin::SpaceFillingCurve::Morton, spin::BucketedBVHGrouping::SharedPrefix},
  {"hilbert/b8/sequential", 8, spin::SpaceFillingCurve::Hilbert, spin::BucketedBVHGrouping::Sequential},
  {"hilbert/b8/shared_prefix", 8, spin::SpaceFillingCurve::Hilbert, spin::BucketedBVHGrouping::SharedPrefix},
}};

enum class Phase
{
  Preprocess,
  Query
};

template <typename Query, typename Preprocess>
void measureCase(benchmark::State& state,
                 Query& query,
                 Preprocess preprocess,
                 mfem::DataCollection& dc,
                 Phase phase,
                 axom::IndexType input_primitives)
{
  preprocess();  // Warm up construction; query cases also need a built query.
  if(phase == Phase::Preprocess)
  {
    for(auto _ : state)
    {
      preprocess();
    }
    state.counters["input_primitives"] = input_primitives;
    return;
  }

  const axom::primal::WindingTolerances tolerances;
  query.query(dc, tolerances);  // Warm caches and populate the output field.
  for(auto _ : state)
  {
    query.query(dc, tolerances);
  }

  const auto& winding = *dc.GetField("winding");
  for(int i = 0; i < winding.Size(); ++i)
  {
    if(!std::isfinite(winding[i]))
    {
      state.SkipWithError("GWN query returned a non-finite value");
      return;
    }
  }
  const auto query_points = winding.Size();
  state.counters["query_points"] = query_points;
  state.counters["input_primitives"] = input_primitives;
  state.counters["points_per_second"] =
    benchmark::Counter(static_cast<double>(query_points) * state.iterations(),
                       benchmark::Counter::kIsRate);
}

template <typename ExecSpace, bool Curved>
void run3DCase(benchmark::State& state, Problem3D& problem, Case config, Phase phase)
{
  using Query =
    std::conditional_t<Curved, quest::NURBSPatchGWNQuery<ExecSpace, 2>, quest::TriangleGWNQuery<ExecSpace, 2>>;

  Query query;
  query.setBucketSize(config.bucket_size);
  query.setBVHOrdering(config.ordering);
  query.setBucketGrouping(config.grouping);

  mfem::DataCollection dc("gwn_benchmark");
  if(phase == Phase::Query)
  {
    quest::generate_gwn_query_mesh(dc,
                                   problem.bounds,
                                   std::vector<double> {},
                                   std::vector<double> {},
                                   std::vector<int>(3, problem.grid_resolution),
                                   1);
  }

  const auto preprocess = [&]() {
    if constexpr(Curved)
    {
      query.preprocess(problem.reader.getPatchArray(), false, true);
    }
    else
    {
      query.preprocess(&problem.triangles, false);
    }
  };

  const auto input_primitives =
    Curved ? problem.reader.getPatchArray().size() : problem.triangles.getNumberOfCells();
  measureCase(state, query, preprocess, dc, phase, input_primitives);
}

template <typename ExecSpace, bool Curved>
void run2DCase(benchmark::State& state, Problem2D& problem, Case config, Phase phase)
{
  using Query =
    std::conditional_t<Curved, quest::NURBSCurveGWNQuery<ExecSpace, 2>, quest::PolylineGWNQuery<ExecSpace, 2>>;

  Query query;
  query.setBucketSize(config.bucket_size);
  query.setBVHOrdering(config.ordering);
  query.setBucketGrouping(config.grouping);

  mfem::DataCollection dc("gwn_benchmark");
  if(phase == Phase::Query)
  {
    quest::generate_gwn_query_mesh(dc,
                                   problem.bounds,
                                   std::vector<double> {},
                                   std::vector<double> {},
                                   std::vector<int>(2, problem.grid_resolution),
                                   1);
  }

  const auto preprocess = [&]() {
    if constexpr(Curved)
    {
      query.preprocess(problem.curves, false, true);
    }
    else
    {
      query.preprocess(&problem.segments, false);
    }
  };

  const auto input_primitives = Curved ? problem.curves.size() : problem.segments.getNumberOfCells();
  measureCase(state, query, preprocess, dc, phase, input_primitives);
}

template <typename ExecSpace>
void registerPolicy(const char* policy, Problem2D& problem2d, Problem3D& problem3d)
{
  for(const auto phase : {Phase::Preprocess, Phase::Query})
  {
    const char* phase_name = phase == Phase::Query ? "query" : "preprocess";
    for(const auto& config : CASES)
    {
      const std::string curved2d_name = axom::fmt::format("gwn/2d/curved/{}/{}/{}/grid{}",
                                                          phase_name,
                                                          policy,
                                                          config.name,
                                                          problem2d.grid_resolution);
      benchmark::RegisterBenchmark(curved2d_name.c_str(),
                                   [&problem2d, config, phase](benchmark::State& state) {
                                     run2DCase<ExecSpace, true>(state, problem2d, config, phase);
                                   })
        ->UseRealTime();

      const std::string discrete2d_name = axom::fmt::format("gwn/2d/discrete/{}/{}/{}/grid{}",
                                                            phase_name,
                                                            policy,
                                                            config.name,
                                                            problem2d.grid_resolution);
      benchmark::RegisterBenchmark(discrete2d_name.c_str(),
                                   [&problem2d, config, phase](benchmark::State& state) {
                                     run2DCase<ExecSpace, false>(state, problem2d, config, phase);
                                   })
        ->UseRealTime();

      const std::string curved3d_name = axom::fmt::format("gwn/3d/curved/{}/{}/{}/grid{}",
                                                          phase_name,
                                                          policy,
                                                          config.name,
                                                          problem3d.grid_resolution);
      benchmark::RegisterBenchmark(curved3d_name.c_str(),
                                   [&problem3d, config, phase](benchmark::State& state) {
                                     run3DCase<ExecSpace, true>(state, problem3d, config, phase);
                                   })
        ->UseRealTime();

      const std::string discrete3d_name = axom::fmt::format("gwn/3d/discrete/{}/{}/{}/grid{}",
                                                            phase_name,
                                                            policy,
                                                            config.name,
                                                            problem3d.grid_resolution);
      benchmark::RegisterBenchmark(discrete3d_name.c_str(),
                                   [&problem3d, config, phase](benchmark::State& state) {
                                     run3DCase<ExecSpace, false>(state, problem3d, config, phase);
                                   })
        ->UseRealTime();
    }
  }
}

}  // namespace

int main(int argc, char** argv)
{
  slic::SimpleLogger logger;
  slic::setLoggingMsgLevel(slic::message::Warning);

  std::string step_file;
  std::string mesh_file;
#ifdef AXOM_DATA_DIR
  // Select defaults from axom_data if present. Otherwise they must be provided via CLI
  step_file = axom::utilities::filesystem::joinPath(AXOM_DATA_DIR, "quest/step/connector.step");
  mesh_file = axom::utilities::filesystem::joinPath(AXOM_DATA_DIR, "contours/svg/axom_logo.mesh");
#endif
  int grid_resolution_2d = 315;
  int grid_resolution_3d = 45;
  std::string policy = "omp";
  axom::CLI::App app {"2D and 3D generalized winding number benchmarks"};
  app.add_option("--mesh-file", mesh_file, "2D MFEM contour mesh")->check(axom::CLI::ExistingFile);
  app.add_option("--step-file", step_file, "STEP input file")->check(axom::CLI::ExistingFile);
  app.add_option("--grid-resolution-2d", grid_resolution_2d, "2D cells per grid axis")
    ->check(axom::CLI::PositiveNumber);
  app.add_option("--grid-resolution-3d", grid_resolution_3d, "3D cells per grid axis")
    ->check(axom::CLI::PositiveNumber);
  app.add_option("--policy", policy, "Execution policy: seq, omp, or both (default: omp)")
    ->check(axom::CLI::IsMember({"seq", "omp", "both"}));
  app.allow_extras();
  CLI11_PARSE(app, argc, argv);
  if(step_file.empty() || mesh_file.empty())
  {
    SLIC_ERROR("Pass --mesh-file and --step-file when Axom was configured without AXOM_DATA_DIR");
    return 1;
  }

  std::vector<std::string> benchmark_args = app.remaining_for_passthrough();
  benchmark_args.insert(benchmark_args.begin(), argv[0]);
  std::vector<char*> benchmark_argv;
  benchmark_argv.reserve(benchmark_args.size());
  for(auto& arg : benchmark_args)
  {
    benchmark_argv.push_back(arg.data());
  }
  int benchmark_argc = static_cast<int>(benchmark_argv.size());
  benchmark::Initialize(&benchmark_argc, benchmark_argv.data());
  if(benchmark::ReportUnrecognizedArguments(benchmark_argc, benchmark_argv.data()))
  {
    return 1;
  }

  Problem2D problem2d;
  problem2d.grid_resolution = grid_resolution_2d;
  if(!problem2d.load(mesh_file))
  {
    SLIC_ERROR("Failed to read or linearize 2D mesh input " << mesh_file);
    return 1;
  }

  Problem3D problem3d;
  problem3d.grid_resolution = grid_resolution_3d;
  if(!problem3d.load(step_file))
  {
    SLIC_ERROR("Failed to read or triangulate STEP input " << step_file);
    return 1;
  }

  if(policy == "seq" || policy == "both")
  {
    registerPolicy<axom::SEQ_EXEC>("seq", problem2d, problem3d);
  }
  if(policy == "omp" || policy == "both")
  {
    registerPolicy<axom::OMP_EXEC>("omp", problem2d, problem3d);
  }
  benchmark::RunSpecifiedBenchmarks();
  return 0;
}
