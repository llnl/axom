// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/Shape.hpp"
#include "axom/quest/DiscreteShape.hpp"

#include "gtest/gtest.h"

#include <utility>

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;

TEST(quest_discrete_shape, appliesNonAffinePointTransform)
{
  klee::TransformableGeometryProperties properties {klee::Dimensions::Three, klee::LengthUnit::cm};
  auto transform = std::make_shared<klee::PointTransform>(
    [](const klee::Geometry::Point3D& p) {
      return klee::Geometry::Point3D {p[0] + 10., p[1] * 2., p[2] - 1.};
    },
    properties,
    "operators/1/transform",
    "Error evaluating callback for 'transform'");

  primal::Tetrahedron<double, 3> tet {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
  klee::Geometry geometry {properties, tet, transform};
  klee::Shape shape {"warped_tet", "steel", {}, {}, std::move(geometry)};

  quest::DiscreteShape discreteShape {shape, nullptr};
  auto mesh = discreteShape.createMeshRepresentation();
  ASSERT_TRUE(mesh);
  ASSERT_EQ(4, mesh->getNumberOfNodes());

  const double* x = mesh->getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh->getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh->getCoordinateArray(axom::mint::Z_COORDINATE);

  EXPECT_DOUBLE_EQ(10., x[0]);
  EXPECT_DOUBLE_EQ(0., y[0]);
  EXPECT_DOUBLE_EQ(-1., z[0]);

  EXPECT_DOUBLE_EQ(11., x[1]);
  EXPECT_DOUBLE_EQ(0., y[1]);
  EXPECT_DOUBLE_EQ(-1., z[1]);

  EXPECT_DOUBLE_EQ(10., x[2]);
  EXPECT_DOUBLE_EQ(2., y[2]);
  EXPECT_DOUBLE_EQ(-1., z[2]);

  EXPECT_DOUBLE_EQ(10., x[3]);
  EXPECT_DOUBLE_EQ(0., y[3]);
  EXPECT_DOUBLE_EQ(0., z[3]);
}
