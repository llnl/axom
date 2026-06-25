// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"

#include <cstdlib>

namespace
{
template <typename ExecSpace>
struct CopySlices
{
  static void run(int* dst, const int* src, axom::IndexType nslices, axom::IndexType sliceSize)
  {
    axom::for_all<ExecSpace>(nslices, [=](axom::IndexType slice) {
      axom::copy(dst + slice * sliceSize, src + slice * sliceSize, sliceSize * sizeof(int));
    });
  }
};
}  // namespace

TEST(core_copy_openmp_race, first_use_inside_for_all)
{
  constexpr axom::IndexType nslices = 128;
  constexpr axom::IndexType sliceSize = 64;
  constexpr axom::IndexType size = nslices * sliceSize;

  auto* src = static_cast<int*>(std::malloc(size * sizeof(int)));
  auto* dst = static_cast<int*>(std::malloc(size * sizeof(int)));

  ASSERT_NE(src, nullptr);
  ASSERT_NE(dst, nullptr);

  for(axom::IndexType i = 0; i < size; ++i)
  {
    src[i] = static_cast<int>(i);
    dst[i] = -1;
  }

  CopySlices<axom::OMP_EXEC>::run(dst, src, nslices, sliceSize);

  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(dst[i], src[i]);
  }

  std::free(dst);
  std::free(src);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  return RUN_ALL_TESTS();
}
