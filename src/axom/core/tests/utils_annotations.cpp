// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/RAII.hpp"
#include "axom/core/utilities/Annotations.hpp"
#include "axom/core/AnnotationMacros.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#ifdef AXOM_USE_CALIPER
  #include "caliper/cali.h"
#endif

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

#include <map>

namespace
{
std::string s_annotation_mode {"none"};

}  // namespace

TEST(utils_annotations, initialize_finalize)
{
#ifdef AXOM_USE_MPI
  axom::utilities::annotations::initialize(MPI_COMM_WORLD, "none");
#else
  axom::utilities::annotations::initialize("none");
#endif

  axom::utilities::annotations::finalize();

  SUCCEED();
}

TEST(utils_annotations, print_adiak_metadata)
{
#ifdef AXOM_USE_MPI
  axom::utilities::annotations::initialize(MPI_COMM_WORLD, "none");
#else
  axom::utilities::annotations::initialize("none");
#endif

  const auto default_metadata = axom::utilities::annotations::retrieve_metadata();

#ifndef AXOM_USE_ADIAK
  EXPECT_EQ(std::size_t {0}, default_metadata.size());
#else
  // some default metadata is present, e.g. user
  EXPECT_TRUE(default_metadata.find("user") != default_metadata.end());

  // but other metadata is not present before we register it
  EXPECT_TRUE(default_metadata.find("custom_int_metadata") == default_metadata.end());
  EXPECT_TRUE(default_metadata.find("custom_str_metadata") == default_metadata.end());
#endif

  // register some additional metadata
  const std::string category {"test_params"};
  AXOM_ANNOTATE_METADATA("custom_str_metadata", "arbitrary string", category);
  AXOM_ANNOTATE_METADATA("custom_int_metadata", 42, category);
  AXOM_ANNOTATE_METADATA("test_annotation_mode", s_annotation_mode, category);

  // check that registered metadata is now present
  auto updated_metadata = axom::utilities::annotations::retrieve_metadata();
#ifndef AXOM_USE_ADIAK
  EXPECT_EQ(std::size_t {0}, updated_metadata.size());
#else
  // some default metadata is present
  EXPECT_TRUE(updated_metadata.find("user") != updated_metadata.end());

  // but other metadata is not present before we register it
  EXPECT_TRUE(updated_metadata.find("custom_int_metadata") != updated_metadata.end());
  EXPECT_TRUE(updated_metadata.find("custom_str_metadata") != updated_metadata.end());

  EXPECT_EQ(std::to_string(42), updated_metadata["custom_int_metadata"]);
  EXPECT_EQ("arbitrary string", updated_metadata["custom_str_metadata"]);

  // print the key-value pairs on rank 0
  {
    int my_rank {0};
  #ifdef AXOM_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  #endif
    if(my_rank == 0)
    {
      std::cout << "Adiak metadata: \n";
      for(const auto &kv : updated_metadata)
      {
        std::cout << axom::fmt::format("- {}: {}\n", kv.first, kv.second);
      }
    }
  }
#endif  // AXOM_USE_ADIAK

  axom::utilities::annotations::finalize();
}

TEST(utils_annotations, check_modes)
{
  EXPECT_TRUE(axom::utilities::annotations::detail::is_mode_valid("none"));

  for(const auto &m : {"counts", "file", "trace", "report", "gputx", "nvprof", "nvtx", "roctx"})
  {
#ifdef AXOM_USE_CALIPER
    EXPECT_TRUE(axom::utilities::annotations::detail::is_mode_valid(m));
#else
    EXPECT_FALSE(axom::utilities::annotations::detail::is_mode_valid(m));
#endif
  }

  EXPECT_FALSE(axom::utilities::annotations::detail::is_mode_valid("_other_"));
}

TEST(utils_annotations, modes)
{
  std::cout << "Testing caliper service '" << s_annotation_mode << "'\n";

#ifdef AXOM_USE_MPI
  axom::utilities::annotations::initialize(MPI_COMM_WORLD, s_annotation_mode);
#else
  axom::utilities::annotations::initialize(s_annotation_mode);
#endif

#ifdef AXOM_USE_CALIPER
  {
    AXOM_ANNOTATE_BEGIN("my region");

    for(int i = 0; i < 10; ++i)
    {
      AXOM_ANNOTATE_SCOPE("inner1");
    }

    for(int i = 0; i < 100; ++i)
    {
      AXOM_ANNOTATE_SCOPE("inner2");
    }
    AXOM_ANNOTATE_END("my region");
  }
#endif

  axom::utilities::annotations::finalize();

  SUCCEED();
}

TEST(utils_annotations, print_help)
{
  int my_rank {0}, num_ranks {1};

#ifdef AXOM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#endif

  // This prints a lot, so let's only print for the "none" mode test
  // when running w/ a single rank
  if(s_annotation_mode == "none" && num_ranks == 1 && my_rank == 0)
  {
    std::cout << "Caliper help string: \n"
              << axom::utilities::annotations::detail::mode_help_string() << std::endl;
  }

  SUCCEED();
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  // Initialize MPI if axom is configured w/ MPI
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);

  // Parse annotation mode for the current test invocation
  axom::CLI::App app {"Axom annotation tests"};
  app.add_option("-m,--mode", s_annotation_mode, "Annotation mode")->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  // run tests
  int result = RUN_ALL_TESTS();

  return result;
}
