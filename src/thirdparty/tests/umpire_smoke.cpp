// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Umpire includes
#include "umpire/config.hpp"
#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"

#include "gtest/gtest.h"

#include <iostream>

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(umpire_smoke, check_version)
{
  // Umpire version must be 2.1.0 or greater
  EXPECT_TRUE(UMPIRE_VERSION_MAJOR >= 2);
  if(UMPIRE_VERSION_MAJOR == 2)
  {
    EXPECT_TRUE(UMPIRE_VERSION_MINOR >= 1);
  }
  else
  {
    EXPECT_TRUE(UMPIRE_VERSION_MINOR >= 0);
  }
  EXPECT_TRUE(UMPIRE_VERSION_PATCH >= 0);

  std::cout << "Umpire version: "           //
            << UMPIRE_VERSION_MAJOR << "."  //
            << UMPIRE_VERSION_MINOR << "."  //
            << UMPIRE_VERSION_PATCH << std::endl;
}

//------------------------------------------------------------------------------
TEST(umpire_smoke, basic_use)
{
  constexpr int N = 100;
  constexpr int BYTESIZE = N * sizeof(int);
  constexpr int MAGIC_VAL = 42;

  auto& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator("HOST");
  EXPECT_EQ(allocator.getName(), "HOST");
  EXPECT_TRUE(allocator.getId() >= 0);

  int* data = static_cast<int*>(allocator.allocate(BYTESIZE));
  EXPECT_TRUE(data != nullptr);
  EXPECT_EQ(BYTESIZE, allocator.getCurrentSize());

  for(int i = 0; i < N; ++i)
  {
    data[i] = MAGIC_VAL;
  }

  auto found_allocator = rm.getAllocator(data);
  EXPECT_EQ(found_allocator.getId(), allocator.getId());
  EXPECT_EQ(found_allocator.getName(), allocator.getName());

  allocator.deallocate(data);
  data = nullptr;
}

TEST(umpire_smoke, available_allocators)
{
  auto& rm = umpire::ResourceManager::getInstance();

  std::cout << "Available Umpire allocators: \n";
  for(const auto& name : rm.getAllocatorNames())
  {
    auto alloc = rm.getAllocator(name);
    std::cout << "  {name: " << alloc.getName()             //
              << ", id: " << alloc.getId()                  //
              << ", strategy: " << alloc.getStrategyName()  //
    // umpire@2022.10 introduced a function to convert a Platform to a string
#if (UMPIRE_VERSION_MAJOR > 2022) || (UMPIRE_VERSION_MAJOR == 2022 && UMPIRE_VERSION_MINOR > 3)
              << ", platform: " << umpire::platform_to_string(alloc.getPlatform())  //
#endif
              << "}\n";
  }
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  return (result);
}
