// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file spin_bvh_benchmark.cpp
 * \brief Measure regular, single-primitive BVH construction and point queries.
 *
 * Build with ENABLE_BENCHMARKS=ON, then run from the build directory:
 * \code{.sh}
 * OMP_NUM_THREADS=8 ./tests/spin_bvh_benchmark \
 *   --benchmark_filter='spin_bvh/(build|query)/(seq|omp)/' \
 *   --benchmark_min_time=1s --benchmark_repetitions=5 \
 *   --benchmark_display_aggregates_only=true --benchmark_time_unit=ms |
 *   awk '$1 ~ /_median$/ { print $1, $2, $3 }'
 * \endcode
 * Use --benchmark_list_tests to see the available cases; omit the filter to run
 * every configured policy. CUDA and HIP cases require a corresponding GPU.
 * The command prints each case name and its median Time in milliseconds. Remove
 * the awk filter to see the full Google Benchmark table. Add
 * --benchmark_out=spin-bvh.json --benchmark_out_format=json to save full results;
 * Google Benchmark still prints the console table.
 *
 * The fixed workload uses 50,000 boxes and 10,000 point queries, with uniform
 * and clustered distributions. Timings exclude input generation, device
 * transfers, warmup, and result validation; GPU work is synchronized before
 * each timer stops. Compare matching case names from separate builds.
 */

#include "axom/config.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/execution/synchronize.hpp"
#include "axom/slic.hpp"
#include "axom/spin/BVH.hpp"

#include "benchmark/benchmark.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using Box = axom::primal::BoundingBox<double, 3>;
using Point = axom::primal::Point<double, 3>;
using Clock = std::chrono::steady_clock;

bool s_benchmark_failed = false;

// Scramble an integer deterministically for both input generation and result hashing.
std::uint64_t mix(std::uint64_t x)
{
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// Use 53 mixed bits to make a value in [0, 1), with a distinct key for each axis.
double unit(std::uint64_t index, int axis)
{
  return static_cast<double>(mix(index * 3 + axis) >> 11) * (1.0 / 9007199254740992.0);
}

// Spread centers through the domain or place them in four dense xy groups.
Point center(axom::IndexType index, bool clustered)
{
  if(clustered)
  {
    const double base = (index % 4) < 2 ? 0.25 : 0.65;
    const double other = (index % 4) % 2 == 0 ? 0.25 : 0.65;
    return Point {base + 0.15 * unit(index, 0),
                  other + 0.15 * unit(index, 1),
                  0.3 + 0.15 * unit(index, 2)};
  }
  return Point {0.05 + 0.9 * unit(index, 0), 0.05 + 0.9 * unit(index, 1), 0.05 + 0.9 * unit(index, 2)};
}

// Validate the query output and summarize its candidate IDs for run comparisons.
std::uint64_t fingerprint(const axom::Array<axom::IndexType>& offsets,
                          const axom::Array<axom::IndexType>& counts,
                          const axom::Array<axom::IndexType>& candidates,
                          axom::IndexType box_count)
{
  std::uint64_t hash = 1469598103934665603ULL;
  axom::IndexType expected_offset = 0;
  for(axom::IndexType q = 0; q < counts.size(); ++q)
  {
    const auto offset = offsets[q];
    const auto count = counts[q];
    if(offset != expected_offset || count < 0 || offset + count > candidates.size())
    {
      throw std::runtime_error("Invalid candidate range");
    }
    std::vector<axom::IndexType> ids;
    ids.reserve(count);
    for(axom::IndexType j = 0; j < count; ++j)
    {
      const auto id = candidates[offset + j];
      if(id < 0 || id >= box_count)
      {
        throw std::runtime_error("Candidate ID outside input range");
      }
      ids.push_back(id);
    }
    // Candidate order can vary by execution policy, so hash each query in sorted order.
    std::sort(ids.begin(), ids.end());
    // Even-numbered queries use a box center and must find that box.
    if(q % 2 == 0 && !std::binary_search(ids.begin(), ids.end(), q % box_count))
    {
      throw std::runtime_error("Box-center query did not return its own box");
    }
    // Include the count so the digest preserves each query's candidate boundary.
    hash = mix(hash ^ static_cast<std::uint64_t>(count));
    for(const auto id : ids)
    {
      hash = mix(hash ^ static_cast<std::uint64_t>(id));
    }
    expected_offset += count;
  }
  if(expected_offset != candidates.size())
  {
    throw std::runtime_error("Candidate array length does not match counts");
  }
  return hash;
}

template <typename ExecSpace>
void runBenchmark(benchmark::State& state, bool clustered, bool build)
{
  using BVH = axom::spin::BVH<3, ExecSpace, double, axom::spin::BVHType::LinearBVH>;
  const auto box_count = static_cast<axom::IndexType>(state.range(0));
  const auto query_count = static_cast<axom::IndexType>(state.range(1));

  try
  {
    const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    const int exec_allocator = axom::execution_space<ExecSpace>::allocatorID();
    if(exec_allocator == axom::INVALID_ALLOCATOR_ID)
    {
      throw std::runtime_error("Execution-space allocator is unavailable");
    }

    axom::Array<Box> host_boxes(box_count, box_count, host_allocator);
    axom::Array<Point> host_queries(query_count, query_count, host_allocator);
    const double half_width = clustered ? 0.006 : 0.004;
    for(axom::IndexType i = 0; i < box_count; ++i)
    {
      const Point c = center(i, clustered);
      host_boxes[i] = Box(Point {c[0] - half_width, c[1] - half_width, c[2] - half_width},
                          Point {c[0] + half_width, c[1] + half_width, c[2] + half_width});
    }
    for(axom::IndexType i = 0; i < query_count; ++i)
    {
      host_queries[i] =
        i % 2 == 0 ? center(i % box_count, clustered) : center(box_count + i, clustered);
    }

    // Input creation and host/device copies are outside the measured operation.
    axom::Array<Box> boxes(host_boxes, exec_allocator);
    axom::Array<Point> queries(host_queries, exec_allocator);
    axom::Array<axom::IndexType> offsets(query_count, query_count, exec_allocator);
    axom::Array<axom::IndexType> counts(query_count, query_count, exec_allocator);
    axom::Array<axom::IndexType> candidates;
    BVH bvh;
    bvh.setScaleFactor(1.0);

    // Warm up both paths before measuring either one.
    if(bvh.initialize(boxes.view(), box_count) != axom::spin::BVH_BUILD_OK)
    {
      throw std::runtime_error("BVH build failed");
    }
    axom::synchronize<ExecSpace>();
    bvh.findPoints(offsets.view(), counts.view(), candidates, query_count, queries.view());
    axom::synchronize<ExecSpace>();

    for(auto _ : state)
    {
      const auto start = Clock::now();
      if(build)
      {
        if(bvh.initialize(boxes.view(), box_count) != axom::spin::BVH_BUILD_OK)
        {
          throw std::runtime_error("BVH build failed");
        }
      }
      else
      {
        bvh.findPoints(offsets.view(), counts.view(), candidates, query_count, queries.view());
      }
      axom::synchronize<ExecSpace>();
      const auto end = Clock::now();
      state.SetIterationTime(std::chrono::duration<double>(end - start).count());
    }

    if(build)
    {
      bvh.findPoints(offsets.view(), counts.view(), candidates, query_count, queries.view());
      axom::synchronize<ExecSpace>();
    }

    // Validate after timing, including a host copy for device results.
    const axom::Array<axom::IndexType> host_offsets(offsets, host_allocator);
    const axom::Array<axom::IndexType> host_counts(counts, host_allocator);
    const axom::Array<axom::IndexType> host_candidates(candidates, host_allocator);
    const auto hash = fingerprint(host_offsets, host_counts, host_candidates, box_count);
    state.SetLabel("fingerprint=" + std::to_string(hash));
    state.counters["candidates"] = static_cast<double>(host_candidates.size());
    state.SetItemsProcessed(state.iterations() * (build ? box_count : query_count));
  }
  catch(const std::exception& error)
  {
    s_benchmark_failed = true;
    state.SkipWithError(error.what());
  }
}

template <typename ExecSpace>
void registerPolicy(const std::string& policy)
{
  for(const bool clustered : {false, true})
  {
    const std::string distribution = clustered ? "clustered" : "uniform";
    for(const bool build : {true, false})
    {
      const std::string operation = build ? "build" : "query";
      const std::string name = "spin_bvh/" + operation + "/" + policy + "/" + distribution;
      benchmark::RegisterBenchmark(
        name,
        [=](benchmark::State& state) { runBenchmark<ExecSpace>(state, clustered, build); })
        ->Args({50000, 10000})
        ->UseManualTime();
    }
  }
}
}  // namespace

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;
  benchmark::Initialize(&argc, argv);
  if(benchmark::ReportUnrecognizedArguments(argc, argv))
  {
    return 1;
  }

  registerPolicy<axom::SEQ_EXEC>("seq");
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  registerPolicy<axom::OMP_EXEC>("omp");
#endif
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  registerPolicy<axom::CUDA_EXEC<256>>("cuda");
#endif
#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  registerPolicy<axom::HIP_EXEC<256>>("hip");
#endif

  benchmark::RunSpecifiedBenchmarks();
  return s_benchmark_failed ? 1 : 0;
}
