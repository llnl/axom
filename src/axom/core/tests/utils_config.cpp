// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"

#ifdef AXOM_USE_OPENMP
  #include <omp.h>
#endif

#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

TEST(utils_config, axom_version)
{
  const int AXOM_MAJOR = AXOM_VERSION_MAJOR;
  const int AXOM_MINOR = AXOM_VERSION_MINOR;
  const int AXOM_PATCH = AXOM_VERSION_PATCH;
  EXPECT_TRUE(AXOM_MAJOR >= 0);
  EXPECT_TRUE(AXOM_MINOR >= 0);
  EXPECT_TRUE(AXOM_PATCH >= 0);

  const std::string AXOM_FULL = AXOM_VERSION_FULL;
  EXPECT_FALSE(AXOM_FULL.empty());

  std::cout << "AXOM_VERSION_FULL: " << AXOM_VERSION_FULL << std::endl;
  std::cout << "AXOM_VERSION_MAJOR: " << AXOM_VERSION_MAJOR << std::endl;
  std::cout << "AXOM_VERSION_MINOR: " << AXOM_VERSION_MINOR << std::endl;
  std::cout << "AXOM_VERSION_PATCH: " << AXOM_VERSION_PATCH << std::endl;
}

TEST(utils_config, config_openmp)
{
  // This test checks that the per-target OpenMP guards
  // in our configuration file 'axom/config.hpp' are working properly

#ifdef AXOM_USE_OPENMP
  std::cout << "OpenMP is available in this configuration." << std::endl;

  #pragma omp parallel
  {
    int thId = omp_get_thread_num();
    int thNum = omp_get_num_threads();
    int thMax = omp_get_max_threads();

  #pragma omp critical
    std::cout << "\tMy thread id is: " << thId << "\tNum threads is: " << thNum
              << "\tMax threads is: " << thMax << std::endl;
  }

#else
  std::cout << "OpenMP not available in this configuration." << std::endl;
#endif

  // Tests a simple reduction over integers
  // Note: Non-openmp configurations are compiled with a flag
  //       to ignore usage of unknown openmp pragmas
  const int N = 100;
  int sum = 0;
#pragma omp parallel for reduction(+ : sum)
  for(int i = 1; i <= N; ++i)
  {
    sum += i;
  }
  EXPECT_EQ(N * (N + 1) / 2, sum) << "Bad reduction on first " << N << " integers";
  std::cout << "Sum of first " << N << " numbers is: " << sum << std::endl;
}

#ifdef AXOM_USE_MFEM
TEST(utils_config, mfem_configuration)
{
  #ifdef MFEM_VERSION
  std::cout << "Using mfem version " << MFEM_VERSION_MAJOR << "."  // major version
            << MFEM_VERSION_MINOR << "."                           // minor version
            << MFEM_VERSION_PATCH                                  // patch level
            << std::endl;
  #endif  // MFEM_VERSION

  // Verify that this copy of mfem is configured appropriately with respect to MPI
  {
    bool axomHasMPI = false;
  #ifdef AXOM_USE_MPI
    axomHasMPI = true;
  #endif
    std::cout << "Axom is built " << (axomHasMPI ? "with" : "without") << " MPI" << std::endl;

    bool mfemHasMPI = false;
  #ifdef MFEM_USE_MPI
    mfemHasMPI = true;
  #endif
    std::cout << "mfem is built " << (mfemHasMPI ? "with" : "without") << " MPI" << std::endl;

    if(!axomHasMPI)
    {
      EXPECT_FALSE(mfemHasMPI) << "Axom expects mfem to be built without MPI "
                                  "when it is not built with MPI";
    }
  }

  // Verify that this copy of mfem is configured without Sidre
  {
    bool mfemHasSidre = false;
  #ifdef MFEM_USE_SIDRE
    mfemHasSidre = true;
  #endif

    EXPECT_FALSE(mfemHasSidre) << "Axom expects mfem to be built without Sidre";
  }
}
#endif  // AXOM_USE_MFEM
