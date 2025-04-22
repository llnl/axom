// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Shape.hpp"

#include <stdexcept>

#include "gtest/gtest.h"

namespace axom
{
namespace klee
{
namespace
{
Geometry createTestGeometry()
{
  return Geometry {TransformableGeometryProperties {
                     Dimensions::Three,
                     LengthUnit::miles,
                   },
                   "test format",
                   "test path",
                   nullptr};
}

TEST(ShapeTest, replaces_no_lists_given)
{
  Shape shape {"name", "material", {}, {}, createTestGeometry()};
  EXPECT_TRUE(shape.replaces("some_material"));
}

TEST(ShapeTest, replaces_replacement_list_given)
{
  Shape shape {"name", "material", {"mat1", "mat2"}, {}, createTestGeometry()};
  EXPECT_TRUE(shape.replaces("mat1"));
  EXPECT_TRUE(shape.replaces("mat2"));
  EXPECT_FALSE(shape.replaces("mat3"));
}

TEST(ShapeTest, replaces_non_replacement_list_given)
{
  Shape shape {"name", "material", {}, {"mat1", "mat2"}, createTestGeometry()};
  EXPECT_FALSE(shape.replaces("mat1"));
  EXPECT_FALSE(shape.replaces("mat2"));
  EXPECT_TRUE(shape.replaces("mat3"));
}

TEST(ShapeTest, both_replacement_lists_given)
{
  EXPECT_THROW(Shape("name", "material", {"replaced"}, {"not replaced"}, createTestGeometry()),
               std::logic_error);
}
}  // namespace
}  // namespace klee
}  // namespace axom
