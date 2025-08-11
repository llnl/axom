// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/ShapeSet.hpp"

#include "gtest/gtest.h"

#include <stdexcept>

namespace axom
{
namespace klee
{
namespace
{

TEST(ShapeSetTest, dimensions_getAndSet)
{
  ShapeSet shapeSet;

  // Can't query the ShapeSet dimensions until calling setDimensions
  EXPECT_THROW(shapeSet.getDimensions(), std::logic_error);

  {
    shapeSet.setDimensions(Dimensions::Two);

    EXPECT_EQ(Dimensions::Two, shapeSet.getDimensions());
    EXPECT_NE(Dimensions::Three, shapeSet.getDimensions());
  }

  {
    shapeSet.setDimensions(Dimensions::Three);

    EXPECT_NE(Dimensions::Two, shapeSet.getDimensions());
    EXPECT_EQ(Dimensions::Three, shapeSet.getDimensions());
  }
}

}  // namespace
}  // namespace klee
}  // namespace axom
