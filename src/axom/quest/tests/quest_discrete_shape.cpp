// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/KleeError.hpp"
#include "axom/klee/Shape.hpp"
#include "axom/quest/DiscreteShape.hpp"

#include "gtest/gtest.h"

#include <utility>

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;

TEST(quest_discrete_shape, appliesNonAffinePointTransform)
{
  // Test that Geometry correctly reports and applies non-affine transforms
  klee::TransformableGeometryProperties properties {klee::Dimensions::Three, klee::LengthUnit::cm};
  auto transform = std::make_shared<klee::PointTransform>(
    [](const klee::Geometry::Point3D& p) {
      return klee::Geometry::Point3D {p[0] + 10., p[1] * 2., p[2] - 1.};
    },
    properties,
    "operators/1/transform",
    "Error evaluating callback for 'transform'");

  // Use a file-based format that doesn't require Conduit serialization
  klee::Geometry geometry {properties, "stl", "dummy.stl", transform};

  // Verify the geometry has non-affine operators
  EXPECT_TRUE(geometry.hasNonAffineOperators());

  // Verify the transform is applied correctly
  klee::Geometry::Point3D input {1., 2., 3.};
  klee::Geometry::Point3D expected {11., 4., 2.};  // x+10, y*2, z-1
  klee::Geometry::Point3D result = geometry.applyTransform(input);

  EXPECT_DOUBLE_EQ(expected[0], result[0]);
  EXPECT_DOUBLE_EQ(expected[1], result[1]);
  EXPECT_DOUBLE_EQ(expected[2], result[2]);
}

// The following guard is necessary since ShapeMesh.hpp is included in MeshClipperStrategy.hpp
// and has a hard dependency on Sidre
#ifdef AXOM_USE_SIDRE
  #include "axom/quest/MeshClipperStrategy.hpp"

TEST(quest_mesh_clipper, rejectsNonAffineOperators)
{
  // MeshClipperStrategy requires affine transformations and should reject
  // geometries with non-affine operators like PointTransform
  klee::TransformableGeometryProperties properties {klee::Dimensions::Three, klee::LengthUnit::cm};
  auto transform = std::make_shared<klee::PointTransform>(
    [](const klee::Geometry::Point3D& p) { return klee::Geometry::Point3D {p[0] * 2., p[1], p[2]}; },
    properties,
    "operators/1/transform",
    "Error evaluating callback for 'transform'");

  // Use a file-based format that doesn't require Conduit serialization
  klee::Geometry geometry {properties, "stl", "dummy.stl", transform};

  // MeshClipperStrategy constructor should throw for non-affine operators
  try
  {
    quest::experimental::MeshClipperStrategy strategy(geometry);
    FAIL() << "MeshClipperStrategy should reject non-affine operators";
  }
  catch(const klee::KleeError& err)
  {
    // The error should mention non-affine operators
    const std::string error_msg = err.what();
    EXPECT_NE(error_msg.find("non-affine"), std::string::npos) << "Error was: " << error_msg;
  }
}
#endif  // AXOM_USE_SIDRE
