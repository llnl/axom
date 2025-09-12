// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file core_flatmap_perf.cpp
 *  \brief This example measures performance of the FlatMap batched construction
 *   interface and FlatMapView, demonstrating portability between CPU, GPU, and OpenMP.
 */

#include "axom/core/Array.hpp"
#include "axom/core/FlatMap.hpp"
#include "axom/core/FlatMapView.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/core/utilities/Timer.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

struct InputParams
{
  using RuntimePolicy = axom::runtime_policy::Policy;

  axom::IndexType num_elems = 10000;
  RuntimePolicy runtime_policy = RuntimePolicy::seq;
  axom::IndexType rep_count = 100;

public:
  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-n, --numElems", num_elems)->description("Number of elements to insert");

    app.add_option("-p, --policy", runtime_policy)
      ->description("Set runtime policy for test")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    app.add_option("-r, --repCount", rep_count)->description("Number of repetitions to run");

    app.get_formatter()->column_width(60);

    app.parse(argc, argv);
  }
};

/*!
 * \brief Sample an RNG with lookahead.
 *  Based on PCG implementation: https://www.pcg-random.org
 */
AXOM_HOST_DEVICE uint64_t SampleRNG(uint64_t seed, int distance)
{
  uint64_t a = 6364136223846793005ULL;
  uint64_t output = seed;

  // LCG recurrence relation:
  //  x_n+1 = a*x_n mod m (m is just size of integer seed, i.e. 2^64)
  // Closed-form solution is:
  //  x_n+1 = a^n*x_0 mod m

  // Compute modular exponent a^n mod 2^64
  if(distance > 0)
  {
    uint64_t a_n = 1;
    int n = distance;
    while(n > 0)
    {
      if(n % 2 == 1)
      {
        a_n *= a;
      }
      a = a * a;
      n /= 2;
    }
    output *= a_n;
  }

  uint32_t xorshifted = ((output >> 18u) ^ output) >> 27u;
  uint32_t rot = output >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

template <typename ExecPolicy, typename T>
void test_flatmap_init_and_query(axom::IndexType num_elems, axom::IndexType rep_count)
{
  int allocatorID = axom::execution_space<ExecPolicy>::allocatorID();

  axom::utilities::Timer initTimer(false);
  axom::utilities::Timer findTimer(false);
  axom::utilities::Timer rehashTimer(false);
  std::random_device rnd_dev {};
  for(int i = 0; i < rep_count; i++)
  {
    axom::Array<T> keys_vec(num_elems, num_elems, allocatorID);
    axom::Array<T> values_vec(num_elems, num_elems, allocatorID);
    const auto keys = keys_vec.view();
    const auto values = values_vec.view();

    // Get a random seed
    uint64_t rnd = rand();
    rnd <<= 32;
    rnd += rand();

    // Generate random keys and values
    axom::for_all<ExecPolicy>(
      num_elems,
      AXOM_LAMBDA(axom::IndexType index) {
        keys[index] = SampleRNG(rnd, 2 * index);
        values[index] = SampleRNG(rnd, 2 * index + 1);
      });

    // Construct a flat map in a single batch.
    initTimer.start();
    auto new_map =
      axom::FlatMap<T, T>::template create<ExecPolicy>(keys, values, axom::Allocator {allocatorID});
    initTimer.stop();

    // Get a view of the map.
    auto map_view = new_map.view();

    findTimer.start();
    // Test use of FlatMap::find() within a kernel.
    axom::for_all<ExecPolicy>(
      num_elems,
      AXOM_LAMBDA(axom::IndexType index) {
        auto it = map_view.find(keys[index]);
        if(it != map_view.end())
        {
          T value = it->second;
          values[index] = value * 2;
        }
      });
    findTimer.stop();

    rehashTimer.start();
    new_map.rehash(num_elems);
    rehashTimer.stop();
  }

  axom::fmt::print(" - Average construction time: {:.8f} seconds\n",
                   initTimer.elapsedTimeInSec() / rep_count);
  axom::fmt::print(" - Average construction throughput: {:.3f} keys/s\n",
                   (num_elems * rep_count) / initTimer.elapsedTimeInSec());
  axom::fmt::print(" - Average rehash time: {:.8f} seconds\n",
                   rehashTimer.elapsedTimeInSec() / rep_count);
  axom::fmt::print(" - Average rehash throughput: {:.3f} keys/s\n",
                   (num_elems * rep_count) / rehashTimer.elapsedTimeInSec());
  axom::fmt::print(" - Average query time: {:.8f} seconds\n",
                   findTimer.elapsedTimeInSec() / rep_count);
  axom::fmt::print(" - Average query throughput: {:.3f} keys/s\n",
                   (num_elems * rep_count) / findTimer.elapsedTimeInSec());
}

int main(int argc, char** argv)
{
  axom::CLI::App app {"Driver for flat-map performance tests"};

  InputParams params;
  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    auto retval = app.exit(e);
    exit(retval);
  }

  axom::fmt::print("Runtime policy: {}\n", axom::runtime_policy::policyToName(params.runtime_policy));
  axom::fmt::print("Number of key-value pairs: {}\n", params.num_elems);
  axom::fmt::print("Repetition count: {}\n", params.rep_count);

  using RuntimePolicy = axom::runtime_policy::Policy;

  if(params.runtime_policy == RuntimePolicy::seq)
  {
    test_flatmap_init_and_query<axom::SEQ_EXEC, int>(params.num_elems, params.rep_count);
  }
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  else if(params.runtime_policy == RuntimePolicy::omp)
  {
    test_flatmap_init_and_query<axom::OMP_EXEC, int>(params.num_elems, params.rep_count);
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  else if(params.runtime_policy == RuntimePolicy::cuda)
  {
    test_flatmap_init_and_query<axom::CUDA_EXEC<256>, int>(params.num_elems, params.rep_count);
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(params.runtime_policy == RuntimePolicy::hip)
  {
    test_flatmap_init_and_query<axom::HIP_EXEC<256>, int>(params.num_elems, params.rep_count);
  }
#endif
}
