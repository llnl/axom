// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"

#include <conduit/conduit.hpp>

namespace utils = axom::bump::utilities;

namespace
{
AXOM_HOST_DEVICE inline std::int32_t expected_value(axom::IndexType i)
{
  std::int32_t x = static_cast<std::int32_t>(i);
  for(int iter = 0; iter < 2048; iter++)
  {
    x = x * 1664525 + 1013904223;
  }
  return x;
}
}  // namespace

template <typename ExecSpace>
void run_unified_copy_sync_test()
{
  constexpr axom::IndexType N = 1 << 18;

  conduit::Node n_dev;
  utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
  n_dev["values"].set_allocator(c2a.getConduitAllocatorID());
  n_dev["values"].set(conduit::DataType::int32(N));

  auto values = utils::make_array_view<std::int32_t>(n_dev["values"]);
  values.fill(std::int32_t {-1});

  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType i) { values[i] = expected_value(i); });

  conduit::Node n_host;
  utils::copy<axom::SEQ_EXEC>(n_host, n_dev);

  const auto host_values = n_host["values"].as_int32_accessor();
  EXPECT_EQ(host_values[0], expected_value(0));
  EXPECT_EQ(host_values[N / 2], expected_value(N / 2));
  EXPECT_EQ(host_values[N - 1], expected_value(N - 1));
}

TEST(bump_conduit_memory, CopySynchronizesUnifiedMemory)
{
#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  constexpr int HIP_BLOCK_SIZE = 64;
  using hip_exec_async = axom::HIP_EXEC<HIP_BLOCK_SIZE, axom::ASYNC>;
  run_unified_copy_sync_test<hip_exec_async>();
#elif defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  constexpr int CUDA_BLOCK_SIZE = 256;
  using cuda_exec_async = axom::CUDA_EXEC<CUDA_BLOCK_SIZE, axom::ASYNC>;
  run_unified_copy_sync_test<cuda_exec_async>();
#else
  GTEST_SKIP() << "Requires a GPU execution space (HIP or CUDA) with RAJA+Umpire.";
#endif
}

