// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/io/IO.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/KleeError.hpp"

#include "axom/slic.hpp"

#include "KleeMatchers.hpp"

#include "gtest/gtest.h"

#include <fstream>
#include <sstream>

namespace axom
{
namespace klee
{
namespace
{
using primal::Vector3D;
using test::AlmostEqVector;
using ::testing::Contains;
using ::testing::HasSubstr;
using ::testing::Truly;

ShapeSet readShapeSetFromString(const std::string &input)
{
  std::istringstream istream(input);
  return readShapeSet(istream);
}

TEST(IOTest, readShapeSet_noShapes)
{
  auto shapeSet = readShapeSetFromString(R"(
        dimensions: 2
        shapes: [])");
  EXPECT_TRUE(shapeSet.getShapes().empty());
  EXPECT_EQ(Dimensions::Two, shapeSet.getDimensions());
  EXPECT_NE(Dimensions::Three, shapeSet.getDimensions());
  EXPECT_NE(Dimensions::Unspecified, shapeSet.getDimensions());
}

TEST(IOTest, readShapeSet_invalidDimensions)
{
  try
  {
    readShapeSetFromString(R"(
      dimensions: 5
      shapes: [])");
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &err)
  {
    EXPECT_THAT(err.what(), HasSubstr("dimensions"));
  }
}

TEST(IOTest, readShapeSet_shapeWithNoReplacementLists)
{
  auto shapeSet = readShapeSetFromString(R"(
        dimensions: 2
        shapes:
          - name: wheel
            material: steel
            geometry:
              format: test_format
              path: path/to/file.format
    )");

  auto &shapes = shapeSet.getShapes();
  ASSERT_EQ(1u, shapes.size());

  auto &shape = shapes[0];
  EXPECT_TRUE(shape.replaces("mat1"));
  EXPECT_TRUE(shape.replaces("mat2"));
  EXPECT_EQ("wheel", shape.getName());
  EXPECT_EQ("steel", shape.getMaterial());

  auto &geometry = shape.getGeometry();
  EXPECT_EQ("test_format", geometry.getFormat());
  EXPECT_EQ("path/to/file.format", geometry.getPath());
  EXPECT_FALSE(geometry.getGeometryOperator());

  EXPECT_EQ(geometry.getInputDimensions(), Dimensions::Two);
  EXPECT_EQ(geometry.getOutputDimensions(), Dimensions::Two);
}

TEST(IOTest, readShapeSet_shapeWithReplacesList)
{
  auto shapeSet = readShapeSetFromString(R"(
        dimensions: 2
        shapes:
          - name: wheel
            material: steel
            replaces: [mat1, mat2]
            geometry:
              format: test_format
              path: path/to/file.format
)");

  auto &shapes = shapeSet.getShapes();
  ASSERT_EQ(1u, shapes.size());

  auto &shape = shapes[0];
  EXPECT_TRUE(shape.replaces("mat1"));
  EXPECT_TRUE(shape.replaces("mat2"));
  EXPECT_FALSE(shape.replaces("material_not_in_list"));
}

TEST(IOTest, readShapeSet_shapeWithDoesNotReplaceList)
{
  auto shapeSet = readShapeSetFromString(R"(
        dimensions: 2
        shapes:
          - name: wheel
            material: steel
            does_not_replace: [mat1, mat2]
            geometry:
              format: test_format
              path: path/to/file.format
    )");

  auto &shapes = shapeSet.getShapes();
  ASSERT_EQ(1u, shapes.size());

  auto &shape = shapes[0];
  EXPECT_FALSE(shape.replaces("mat1"));
  EXPECT_FALSE(shape.replaces("mat2"));
  EXPECT_TRUE(shape.replaces("material_not_in_list"));
}

TEST(IOTest, readShapeSet_missingName)
{
  auto input = R"(
    dimensions: 2

    shapes:
      - material: steel
        geometry:
          format: test_format
          path: path/to/my.file
  )";

  try
  {
    readShapeSetFromString(input);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &err)
  {
    EXPECT_THAT(err.what(), HasSubstr("name"));
  }
}

TEST(IOTest, readShapeSet_missingMaterial)
{
  auto input = R"(
    dimensions: 2

    shapes:
      - name: wheel
        geometry:
          format: test_format
          path: path/to/my.file
  )";

  try
  {
    readShapeSetFromString(input);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &err)
  {
    EXPECT_THAT(err.what(), HasSubstr("material"));
  }
}

TEST(IOTest, readShapeSet_missingGeometryPath)
{
  // Expect a validation error when 'geometry/path' is missing
  // and 'geometry/format' is not "none"
  {
    auto input = R"(
      dimensions: 2

      shapes:
        - name: wheel
          material: steel
          geometry:
            format: test_format
    )";

    try
    {
      readShapeSetFromString(input);
      FAIL() << "Should have thrown";
    }
    catch(const KleeError &err)
    {
      EXPECT_THAT(err.what(), HasSubstr("Provided format"));
    }
  }

  // Valid when 'geometry/path' is missing and 'geometry/format' is none
  {
    auto input = R"(
      dimensions: 2

      shapes:
        - name: wheel
          material: steel
          geometry:
            format: none
    )";

    try
    {
      readShapeSetFromString(input);
      SUCCEED();
    }
    catch(const KleeError &err)
    {
      FAIL() << "Should not have thrown. Error message: " << err.what();
    }
  }
}

TEST(IOTest, readShapeSet_formatGeometryFormat)
{
  auto input = R"(
    dimensions: 2

    shapes:
      - name: wheel
        material: steel
        geometry:
          path: my/file.format
  )";

  try
  {
    readShapeSetFromString(input);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &err)
  {
    EXPECT_THAT(err.what(), HasSubstr("format"));
  }
}

TEST(IOTest, readShapeSet_file)
{
  std::string fileName = "testFile.yaml";

  std::string fileContents = R"(
    dimensions: 2

    shapes:
      - name: wheel
        material: steel
        geometry:
          format: test_format
          path: path/to/file.format)";
  std::ofstream fout {fileName};
  fout << fileContents;
  fout.close();

  auto shapeSet = readShapeSet(fileName);
  EXPECT_EQ(1u, shapeSet.getShapes().size());
  EXPECT_EQ("testFile.yaml", shapeSet.getPath());
}

TEST(IOTest, readShapeSet_shapeWithReplacesAndDoesNotReplaceLists)
{
  auto input = R"(
      dimensions: 2
      shapes:
        - name: wheel
          material: steel
          replaces: [mat1, mat2]
          does_not_replace: [mat1, mat2]
          geometry:
            format: test_format
            path: path/to/file.format
  )";
  try
  {
    readShapeSetFromString(input);
  }
  catch(const KleeError &error)
  {
    EXPECT_THAT(error.getErrors(), Contains(Truly([](const inlet::VerificationError &err) {
                  return err.path == axom::Path {"shapes/_inlet_collection/0"} &&
                    err.messageContains("replaces") && err.messageContains("does_not_replace");
                })));
  }
}

TEST(IOTest, readShapeSet_geometryOperators)
{
  auto shapeSet = readShapeSetFromString(R"(
      dimensions: 2
      shapes:
        - name: wheel
          material: steel
          geometry:
            format: test_format
            path: path/to/file.format
            units: m
            operators:
              - rotate: 90
              - translate: [10, 20]
    )");
  auto &shapes = shapeSet.getShapes();
  ASSERT_EQ(1u, shapes.size());
  auto &shape = shapes[0];
  auto &geometryOperator = shape.getGeometry().getGeometryOperator();
  ASSERT_TRUE(geometryOperator);
  auto composite = std::dynamic_pointer_cast<const CompositeOperator>(geometryOperator);
  ASSERT_TRUE(composite);

  EXPECT_EQ(2u, composite->getOperators().size());

  auto rotation = dynamic_cast<const Rotation *>(composite->getOperators()[0].get());
  ASSERT_NE(rotation, nullptr);
  EXPECT_EQ(rotation->getAngle(), 90);

  auto translation = dynamic_cast<const Translation *>(composite->getOperators()[1].get());
  ASSERT_NE(translation, nullptr);
  EXPECT_THAT(translation->getOffset(), AlmostEqVector(Vector3D {10, 20, 0}));
  EXPECT_EQ(LengthUnit::m, translation->getEndProperties().units);
  EXPECT_EQ(shapeSet.getDimensions(), translation->getEndProperties().dimensions);
}

TEST(IOTest, readShapeSet_geometryOperatorsWithoutUnits)
{
  try
  {
    readShapeSetFromString(R"(
      dimensions: 2
      shapes:
        - name: wheel
          material: steel
          geometry:
            format: test_format
            path: path/to/file.format
            operators:
              - rotate: 90
              - translate: [10, 20]
    )");
    FAIL() << "Expected a failure";
  }
  catch(const KleeError &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("operator"));
    EXPECT_THAT(ex.what(), HasSubstr("units"));
  }
}

TEST(IOTest, readShapeSet_geometryOperatorsWithUnits)
{
  auto shapeSet = readShapeSetFromString(R"(
      dimensions: 2
      shapes:
        - name: wheel
          material: steel
          geometry:
            format: test_format
            path: path/to/file.format
            units: mm
            operators:
              - rotate: 90
              - translate: [10, 20]
    )");
  auto &shapes = shapeSet.getShapes();
  ASSERT_EQ(1u, shapes.size());

  auto &geometryOperator = shapes[0].getGeometry().getGeometryOperator();
  ASSERT_TRUE(geometryOperator);

  auto composite = std::dynamic_pointer_cast<const CompositeOperator>(geometryOperator);
  ASSERT_TRUE(composite);

  EXPECT_EQ(2u, composite->getOperators().size());
  auto translation = dynamic_cast<const Translation *>(composite->getOperators()[1].get());
  ASSERT_NE(translation, nullptr);
  EXPECT_THAT(translation->getOffset(), AlmostEqVector(Vector3D {10, 20, 0}));
}

TEST(IOTest, readShapeSet_differentDimensions)
{
  auto shapeSet = readShapeSetFromString(R"(
      dimensions: 2
      shapes:
        - name: wheel
          material: steel
          geometry:
            format: test_format
            path: path/to/file.format
            units: cm
            start_dimensions: 3
            operators:
              - slice:
                 x: 10
    )");
  auto &shapes = shapeSet.getShapes();
  ASSERT_EQ(1u, shapes.size());
  auto &shape = shapes[0];
  auto &geometry = shape.getGeometry();
  TransformableGeometryProperties expectedStartProperties {Dimensions::Three, LengthUnit::cm};
  TransformableGeometryProperties expectedEndProperties {Dimensions::Two, LengthUnit::cm};
  EXPECT_EQ(expectedStartProperties, geometry.getStartProperties());
  EXPECT_EQ(expectedEndProperties, geometry.getEndProperties());
  EXPECT_EQ(shapeSet.getDimensions(), geometry.getEndProperties().dimensions);
  auto &geometryOperator = geometry.getGeometryOperator();
  ASSERT_TRUE(geometryOperator);
  auto composite = std::dynamic_pointer_cast<const CompositeOperator>(geometryOperator);
  ASSERT_TRUE(composite);
  EXPECT_EQ(1u, composite->getOperators().size());
  auto slice = std::dynamic_pointer_cast<const SliceOperator>(composite->getOperators()[0]);
  EXPECT_TRUE(slice);
}

TEST(IOTest, readShapeSet_explicitDimensions)
{
  // let's start w/ a global dimension of 2
  {
    auto shapeSet = readShapeSetFromString(R"(
      dimensions: 2
      shapes:
        - name: no_explicit_dims
          material: mat0
          geometry:
            format: test_format
            path: path/to/file0.format
            units: cm
        - name: explicit_dims_same_as_global
          material: mat1
          geometry:
            format: test_format
            path: path/to/file1.format
            units: cm
            dimensions: 2
        - name: explicit_dims_different_from_global
          material: mat2
          geometry:
            format: test_format
            path: path/to/file2.format
            units: cm
            dimensions: 3
        - name: explicit_dims_with_start_dim
          material: mat3
          geometry:
            format: test_format
            path: path/to/file3.format
            units: cm
            start_dimensions: 3
            dimensions: 2
            operators:
              - slice:
                 x: 10
    )");

    auto &shapes = shapeSet.getShapes();
    ASSERT_EQ(4u, shapes.size());

    const Dimensions exp_global_dims {Dimensions::Two};
    EXPECT_EQ(shapeSet.getDimensions(), exp_global_dims);

    // no_explicit_dims -- should be same as global dims
    {
      auto &geometry = shapes[0].getGeometry();
      const Dimensions exp_start_dims {Dimensions::Two};
      const Dimensions exp_end_dims {Dimensions::Two};
      EXPECT_EQ(exp_start_dims, geometry.getInputDimensions());
      EXPECT_EQ(exp_end_dims, geometry.getOutputDimensions());
      EXPECT_EQ(geometry.getOutputDimensions(), exp_global_dims);
    }

    // explicit_dims_same_as_global -- should be same as global dims
    {
      auto &geometry = shapes[1].getGeometry();
      const Dimensions exp_start_dims {Dimensions::Two};
      const Dimensions exp_end_dims {Dimensions::Two};
      EXPECT_EQ(exp_start_dims, geometry.getInputDimensions());
      EXPECT_EQ(exp_end_dims, geometry.getOutputDimensions());
    }

    // explicit_dims_different_from_global -- differs from global dims
    {
      auto &geometry = shapes[2].getGeometry();
      const Dimensions exp_start_dims {Dimensions::Three};
      const Dimensions exp_end_dims {Dimensions::Three};
      EXPECT_EQ(exp_start_dims, geometry.getInputDimensions());
      EXPECT_EQ(exp_end_dims, geometry.getOutputDimensions());
    }

    // explicit_dims_with_start_dim -- changes dimension
    {
      auto &geometry = shapes[3].getGeometry();
      const Dimensions exp_start_dims {Dimensions::Three};
      const Dimensions exp_end_dims {Dimensions::Two};
      EXPECT_EQ(exp_start_dims, geometry.getInputDimensions());
      EXPECT_EQ(exp_end_dims, geometry.getOutputDimensions());
    }
  }

  // next, we'll test a global dimension of 3
  {
    auto shapeSet = readShapeSetFromString(R"(
      dimensions: 3
      shapes:
        - name: no_explicit_dims
          material: mat0
          geometry:
            format: test_format
            path: path/to/file0.format
            units: cm
        - name: explicit_dims_same_as_global
          material: mat1
          geometry:
            format: test_format
            path: path/to/file1.format
            units: cm
            dimensions: 3
        - name: explicit_dims_different_from_global
          material: mat2
          geometry:
            format: test_format
            path: path/to/file2.format
            units: cm
            dimensions: 2
        - name: explicit_dims_with_start_dim
          material: mat3
          geometry:
            format: test_format
            path: path/to/file3.format
            units: cm
            start_dimensions: 3
            dimensions: 2
            operators:
              - slice:
                 x: 10
    )");

    auto &shapes = shapeSet.getShapes();
    ASSERT_EQ(4u, shapes.size());

    const Dimensions exp_global_dims {Dimensions::Three};
    EXPECT_EQ(shapeSet.getDimensions(), exp_global_dims);

    // no_explicit_dims -- should be same as global dims
    {
      auto &geometry = shapes[0].getGeometry();
      const Dimensions exp_start_dims {Dimensions::Three};
      const Dimensions exp_end_dims {Dimensions::Three};
      EXPECT_EQ(exp_start_dims, geometry.getInputDimensions());
      EXPECT_EQ(exp_end_dims, geometry.getOutputDimensions());
      EXPECT_EQ(geometry.getOutputDimensions(), exp_global_dims);
    }

    // explicit_dims_same_as_global -- should be same as global dims
    {
      auto &geometry = shapes[1].getGeometry();
      const Dimensions exp_start_dims {Dimensions::Three};
      const Dimensions exp_end_dims {Dimensions::Three};
      EXPECT_EQ(exp_start_dims, geometry.getInputDimensions());
      EXPECT_EQ(exp_end_dims, geometry.getOutputDimensions());
    }

    // explicit_dims_different_from_global -- differs from global dims
    {
      auto &geometry = shapes[2].getGeometry();
      const Dimensions exp_start_dims {Dimensions::Two};
      const Dimensions exp_end_dims {Dimensions::Two};
      EXPECT_EQ(exp_start_dims, geometry.getInputDimensions());
      EXPECT_EQ(exp_end_dims, geometry.getOutputDimensions());
    }

    // explicit_dims_with_start_dim -- changes dimension
    {
      auto &geometry = shapes[3].getGeometry();
      const Dimensions exp_start_dims {Dimensions::Three};
      const Dimensions exp_end_dims {Dimensions::Two};
      EXPECT_EQ(exp_start_dims, geometry.getInputDimensions());
      EXPECT_EQ(exp_end_dims, geometry.getOutputDimensions());
    }
  }
}

TEST(IOTest, readShapeSet_wrongEndDimensions)
{
  try
  {
    readShapeSetFromString(R"(
        dimensions: 3
        shapes:
          - name: wheel
            material: steel
            geometry:
              format: test_format
              path: path/to/file.format
              start_dimensions: 2
      )");
    FAIL() << "Expected an error";
  }
  catch(const KleeError &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("dimensions"));
  }
}

TEST(IOTest, readShapeSet_namedGeometryOperators)
{
  auto shapeSet = readShapeSetFromString(R"(
      dimensions: 2

      shapes:
        - name: wheel
          material: steel
          geometry:
            format: test_format
            path: path/to/file.format
            units: m
            operators:
              - ref: my_operation

      named_operators:
        - name: my_operation
          units: m
          value:
            - rotate: 90
            - translate: [10, 20]
    )");
  auto &shapes = shapeSet.getShapes();
  ASSERT_EQ(1u, shapes.size());
  auto &shape = shapes[0];
  auto &geometryOperator = shape.getGeometry().getGeometryOperator();
  ASSERT_TRUE(geometryOperator);
  auto composite = std::dynamic_pointer_cast<const CompositeOperator>(geometryOperator);
  ASSERT_TRUE(composite);
  EXPECT_EQ(1u, composite->getOperators().size());
  auto referenced = dynamic_cast<const CompositeOperator *>(composite->getOperators()[0].get());
  EXPECT_EQ(2u, referenced->getOperators().size());
  auto translation = dynamic_cast<const Translation *>(referenced->getOperators()[1].get());
  ASSERT_NE(translation, nullptr);
  EXPECT_THAT(translation->getOffset(), AlmostEqVector(Vector3D {10, 20, 0}));
}

}  // namespace
}  // namespace klee
}  // namespace axom

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  int result = RUN_ALL_TESTS();
  return result;
}
