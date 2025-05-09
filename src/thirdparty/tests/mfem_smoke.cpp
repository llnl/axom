// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: mfem_smoke.cpp
///
//-----------------------------------------------------------------------------
#include "axom/config.hpp"  // for compile-time definitions
#include "mfem.hpp"
#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
TEST(mfem_smoke, check_version)
{
  std::cout << "Using mfem version: "     //
            << MFEM_VERSION_MAJOR << "."  //
            << MFEM_VERSION_MINOR << "."  //
            << MFEM_VERSION_PATCH << std::endl;

  EXPECT_TRUE(MFEM_VERSION_MAJOR >= 0);
  EXPECT_TRUE(MFEM_VERSION_MINOR >= 0);
  EXPECT_TRUE(MFEM_VERSION_PATCH >= 0);
}

//-----------------------------------------------------------------------------
TEST(mfem_smoke, basic_use)
{
  // Simple usage of a basic mfem type
  mfem::Element* el = new mfem::Quadrilateral(0, 1, 2, 3);

  EXPECT_EQ(mfem::Element::QUADRILATERAL, el->GetType());
  EXPECT_EQ(4, el->GetNVertices());

  delete el;
}

TEST(mfem_smoke, basic_device_use)
{
  mfem::Device mfem_device;

  // Check that Backend::CPU is enabled (default constructor)
  EXPECT_TRUE(mfem_device.IsDisabled());

#if defined(AXOM_USE_OPENMP) && defined(MFEM_USE_OPENMP)
  mfem_device.Configure("omp");
  EXPECT_TRUE(mfem_device.Allows(mfem::Backend::OMP));
  mfem_device.Print();
#endif

#if defined(AXOM_USE_CUDA) && defined(MFEM_USE_CUDA)
  mfem_device.Configure("cuda");
  EXPECT_TRUE(mfem_device.Allows(mfem::Backend::CUDA_MASK));
  mfem_device.Print();
#endif

#if defined(AXOM_USE_HIP) && defined(MFEM_USE_HIP)
  mfem_device.Configure("hip");
  EXPECT_TRUE(mfem_device.Allows(mfem::Backend::HIP_MASK));
  mfem_device.Print();
#endif
}
