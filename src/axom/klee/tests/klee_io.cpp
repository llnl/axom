// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/io/IO.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/KleeError.hpp"

#include "axom/config.hpp"
#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/slic.hpp"

#include "KleeMatchers.hpp"

#include "gtest/gtest.h"

#include <fstream>
#include <sstream>

namespace klee = axom::klee;
namespace inlet = axom::inlet;
namespace test = axom::klee::test;
namespace primal = axom::primal;

using klee::CompositeOperator;
using klee::Dimensions;
using klee::InputFormat;
using klee::KleeError;
using klee::LengthUnit;
using klee::Rotation;
using klee::Scale;
using klee::ShapeSet;
using klee::SliceOperator;
using klee::TransformableGeometryProperties;
using klee::Translation;
using primal::Point3D;
using primal::Vector3D;
using test::AlmostEqPoint;
using test::AlmostEqVector;
using ::testing::Contains;
using ::testing::HasSubstr;
using ::testing::Truly;

namespace
{
ShapeSet readShapeSetFromString(const std::string &input)
{
  std::istringstream istream(input);
  return klee::readShapeSet(istream);
}

ShapeSet readShapeSetFromString(const std::string &input, InputFormat format)
{
  std::istringstream istream(input);
  return klee::readShapeSet(istream, format);
}
}  // end namespace

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

  auto shapeSet = klee::readShapeSet(fileName);
  EXPECT_EQ(1u, shapeSet.getShapes().size());
  EXPECT_EQ("testFile.yaml", shapeSet.getPath());
}

TEST(IOTest, readShapeSet_explicitYamlStreamFormat)
{
  auto shapeSet = readShapeSetFromString(R"(
    dimensions: 2
    shapes:
      - name: wheel
        material: steel
        geometry:
          format: test_format
          path: path/to/file.format
  )",
                                         InputFormat::YAML);

  ASSERT_EQ(1u, shapeSet.getShapes().size());
  EXPECT_EQ("wheel", shapeSet.getShapes()[0].getName());
}

TEST(IOTest, readShapeSet_fileWithoutExtensionDefaultsToYaml)
{
  const std::string fileName = "missingKleeInputWithoutExtension";
  ASSERT_FALSE(axom::utilities::filesystem::pathExists(fileName));
  try
  {
    klee::readShapeSet(fileName);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &error)
  {
    ASSERT_EQ(1u, error.getErrors().size());
    EXPECT_EQ(axom::Path {fileName}, error.getErrors()[0].path);
    EXPECT_THAT(error.what(), HasSubstr("Failed to parse YAML Klee input"));
  }
}

TEST(IOTest, readShapeSet_explicitYamlOverridesFileExtension)
{
  axom::utilities::filesystem::TempFile input {"explicitYaml", "lua"};
  input.write(R"(
    dimensions: 2
    shapes: [])");

  auto shapeSet = klee::readShapeSet(input.getPath(), InputFormat::YAML);
  EXPECT_EQ(Dimensions::Two, shapeSet.getDimensions());
  EXPECT_EQ(input.getPath(), shapeSet.getPath());
}

TEST(IOTest, readShapeSet_explicitFormatReportsParseFailure)
{
  axom::utilities::filesystem::TempFile input {"invalidExplicitYaml", "lua"};
  input.write("dimensions: [");
  try
  {
    klee::readShapeSet(input.getPath(), InputFormat::YAML);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &error)
  {
    ASSERT_EQ(1u, error.getErrors().size());
    EXPECT_EQ(axom::Path {input.getPath()}, error.getErrors()[0].path);
    EXPECT_THAT(error.what(), HasSubstr("Failed to parse YAML Klee input"));
  }
}

TEST(IOTest, readShapeSet_emptyStreamReportsParseFailure)
{
  try
  {
    readShapeSetFromString("", InputFormat::YAML);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &error)
  {
    ASSERT_EQ(1u, error.getErrors().size());
    EXPECT_EQ(axom::Path {"<stream>"}, error.getErrors()[0].path);
    EXPECT_STREQ("Failed to parse YAML Klee input from stream.", error.what());
  }
}

TEST(IOTest, readShapeSet_missingFileReportsParseFailure)
{
  const std::string fileName = "missingKleeInput.yaml";
  try
  {
    klee::readShapeSet(fileName);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &error)
  {
    ASSERT_EQ(1u, error.getErrors().size());
    EXPECT_EQ(axom::Path {fileName}, error.getErrors()[0].path);
    EXPECT_STREQ("Failed to parse YAML Klee input from file 'missingKleeInput.yaml'.",
                 error.what());
  }
}

TEST(IOTest, readShapeSet_unsupportedFileExtension)
{
  try
  {
    klee::readShapeSet("testFile.json");
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &err)
  {
    EXPECT_THAT(err.what(), HasSubstr("Unsupported Klee input file extension '.json'"));
    EXPECT_THAT(err.what(), HasSubstr(".yaml, .yml, and .lua"));
  }
}

TEST(IOTest, readShapeSet_streamDefaultsToYaml)
{
  try
  {
    readShapeSetFromString(R"(
      dimensions = 2
      shapes = {}
    )");
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &err)
  {
    EXPECT_THAT(err.what(), HasSubstr("dimensions"));
  }
}

#ifndef AXOM_USE_LUA
TEST(IOTest, readShapeSet_luaUnavailableDiagnostic)
{
  try
  {
    readShapeSetFromString(R"(
      dimensions = 2
      shapes = {}
    )",
                           InputFormat::Lua);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &err)
  {
    EXPECT_STREQ(
      "Lua input decks require Axom configured with AXOM_ENABLE_LUA=ON and Sol library "
      "support. Rebuild Axom with Lua enabled or convert deck to YAML.",
      err.what());
  }
}
#endif

#ifdef AXOM_USE_LUA
TEST(IOTest, readShapeSet_explicitLuaOverridesFileExtension)
{
  axom::utilities::filesystem::TempFile input {"explicitLua", "yaml"};
  input.write(R"(
    dimensions = 2
    shapes = {})");

  auto shapeSet = klee::readShapeSet(input.getPath(), InputFormat::Lua);
  EXPECT_EQ(Dimensions::Two, shapeSet.getDimensions());
  EXPECT_EQ(input.getPath(), shapeSet.getPath());
}

TEST(IOTest, readShapeSet_malformedLuaReportsParseFailure)
{
  try
  {
    readShapeSetFromString("dimensions =", InputFormat::Lua);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &error)
  {
    ASSERT_EQ(1u, error.getErrors().size());
    EXPECT_EQ(axom::Path {"<stream>"}, error.getErrors()[0].path);
    EXPECT_THAT(error.what(), HasSubstr("Failed to parse Lua Klee input from stream"));
  }
}

TEST(IOTest, readShapeSet_luaStreamMinimalShapeList)
{
  auto shapeSet = readShapeSetFromString(R"(
    dimensions = 2
    shapes = {
      {
        name = "wheel",
        material = "steel",
        geometry = {
          format = "test_format",
          path = "path/to/file.format"
        }
      }
    }
  )",
                                         InputFormat::Lua);

  ASSERT_EQ(1u, shapeSet.getShapes().size());
  const auto &shape = shapeSet.getShapes()[0];
  EXPECT_EQ("wheel", shape.getName());
  EXPECT_EQ("steel", shape.getMaterial());
  EXPECT_EQ("test_format", shape.getGeometry().getFormat());
  EXPECT_EQ("path/to/file.format", shape.getGeometry().getPath());
  EXPECT_EQ(Dimensions::Two, shapeSet.getDimensions());
}

TEST(IOTest, readShapeSet_luaFileExtension)
{
  std::string fileName = "testFile.lua";

  std::string fileContents = R"(
    dimensions = 2
    shapes = {
      {
        name = "wheel",
        material = "steel",
        geometry = {
          format = "test_format",
          path = "relative/path.format"
        }
      }
    }
  )";
  std::ofstream fout {fileName};
  fout << fileContents;
  fout.close();

  auto shapeSet = klee::readShapeSet(fileName);
  ASSERT_EQ(1u, shapeSet.getShapes().size());
  EXPECT_EQ("testFile.lua", shapeSet.getPath());
  EXPECT_EQ("relative/path.format", shapeSet.getShapes()[0].getGeometry().getPath());
}

TEST(IOTest, readShapeSet_luaReplacementRules)
{
  auto replaces = readShapeSetFromString(R"(
    dimensions = 2
    shapes = {
      {
        name = "wheel",
        material = "steel",
        replaces = {"mat1", "mat2"},
        geometry = {
          format = "test_format",
          path = "path/to/file.format"
        }
      }
    }
  )",
                                         InputFormat::Lua);
  ASSERT_EQ(1u, replaces.getShapes().size());
  EXPECT_TRUE(replaces.getShapes()[0].replaces("mat1"));
  EXPECT_FALSE(replaces.getShapes()[0].replaces("mat3"));

  auto doesNotReplace = readShapeSetFromString(R"(
    dimensions = 2
    shapes = {
      {
        name = "wheel",
        material = "steel",
        does_not_replace = {"mat1", "mat2"},
        geometry = {
          format = "test_format",
          path = "path/to/file.format"
        }
      }
    }
  )",
                                               InputFormat::Lua);
  ASSERT_EQ(1u, doesNotReplace.getShapes().size());
  EXPECT_FALSE(doesNotReplace.getShapes()[0].replaces("mat1"));
  EXPECT_TRUE(doesNotReplace.getShapes()[0].replaces("mat3"));
}

TEST(IOTest, readShapeSet_luaGeometryOperators)
{
  auto shapeSet = readShapeSetFromString(R"(
    dimensions = 3
    shapes = {
      {
        name = "windshield",
        material = "glass",
        geometry = {
          format = "stl",
          path = "windshield.stl",
          start_units = "m",
          end_units = "cm",
          operators = {
            { rotate = 90, axis = {0, 1, 0}, center = {0, 0, -10} },
            { translate = {10, 20, 30} },
            { scale = {1.5, 2.5, 3.5}, center = {1, 2, 3} },
            { convert_units_to = "cm" }
          }
        }
      },
      {
        name = "slice",
        material = "steel",
        geometry = {
          format = "stl",
          path = "slice.stl",
          start_dimensions = 3,
          dimensions = 2,
          units = "cm",
          operators = {
            { slice = { x = 10 } }
          }
        }
      }
    }
  )",
                                         InputFormat::Lua);

  ASSERT_EQ(2u, shapeSet.getShapes().size());
  const auto &geometryOperator = shapeSet.getShapes()[0].getGeometry().getGeometryOperator();
  ASSERT_TRUE(geometryOperator);
  auto composite = std::dynamic_pointer_cast<const CompositeOperator>(geometryOperator);
  ASSERT_TRUE(composite);
  ASSERT_EQ(4u, composite->getOperators().size());

  auto rotation = dynamic_cast<const Rotation *>(composite->getOperators()[0].get());
  ASSERT_NE(rotation, nullptr);
  EXPECT_EQ(rotation->getAngle(), 90);

  auto translation = dynamic_cast<const Translation *>(composite->getOperators()[1].get());
  ASSERT_NE(translation, nullptr);
  EXPECT_THAT(translation->getOffset(), AlmostEqVector(Vector3D {10, 20, 30}));

  auto scale = dynamic_cast<const Scale *>(composite->getOperators()[2].get());
  ASSERT_NE(scale, nullptr);
  EXPECT_DOUBLE_EQ(1.5, scale->getXFactor());
  EXPECT_DOUBLE_EQ(2.5, scale->getYFactor());
  EXPECT_DOUBLE_EQ(3.5, scale->getZFactor());
  EXPECT_THAT(scale->getCenter(), AlmostEqPoint(Point3D {1, 2, 3}));
  EXPECT_EQ(LengthUnit::cm, composite->getEndProperties().units);

  auto sliceComposite = std::dynamic_pointer_cast<const CompositeOperator>(
    shapeSet.getShapes()[1].getGeometry().getGeometryOperator());
  ASSERT_TRUE(sliceComposite);
  ASSERT_EQ(1u, sliceComposite->getOperators().size());
  EXPECT_TRUE(std::dynamic_pointer_cast<const SliceOperator>(sliceComposite->getOperators()[0]));
}

TEST(IOTest, readShapeSet_luaNamedGeometryOperatorsWithNestedRef)
{
  auto shapeSet = readShapeSetFromString(R"(
    dimensions = 2

    shapes = {
      {
        name = "wheel",
        material = "steel",
        geometry = {
          format = "test_format",
          path = "path/to/file.format",
          units = "m",
          operators = {
            { ref = "outer_operation" }
          }
        }
      }
    }

    named_operators = {
      {
        name = "inner_operation",
        units = "m",
        value = {
          { rotate = 90 }
        }
      },
      {
        name = "outer_operation",
        units = "m",
        value = {
          { ref = "inner_operation" },
          { translate = {10, 20} }
        }
      }
    }
  )",
                                         InputFormat::Lua);

  ASSERT_EQ(1u, shapeSet.getShapes().size());
  auto composite = std::dynamic_pointer_cast<const CompositeOperator>(
    shapeSet.getShapes()[0].getGeometry().getGeometryOperator());
  ASSERT_TRUE(composite);
  ASSERT_EQ(1u, composite->getOperators().size());
  auto referenced = dynamic_cast<const CompositeOperator *>(composite->getOperators()[0].get());
  ASSERT_NE(referenced, nullptr);
  ASSERT_EQ(2u, referenced->getOperators().size());
  auto nested = dynamic_cast<const CompositeOperator *>(referenced->getOperators()[0].get());
  ASSERT_NE(nested, nullptr);
  EXPECT_EQ(1u, nested->getOperators().size());
  auto translation = dynamic_cast<const Translation *>(referenced->getOperators()[1].get());
  ASSERT_NE(translation, nullptr);
  EXPECT_THAT(translation->getOffset(), AlmostEqVector(Vector3D {10, 20, 0}));
}

TEST(IOTest, readShapeSet_luaDifferentDimensions)
{
  auto shapeSet = readShapeSetFromString(R"(
    dimensions = 2
    shapes = {
      {
        name = "flat",
        material = "steel",
        geometry = {
          format = "stl",
          path = "flat.stl",
          dimensions = 3,
          units = "cm"
        }
      },
      {
        name = "sliced",
        material = "glass",
        geometry = {
          format = "stl",
          path = "sliced.stl",
          start_dimensions = 3,
          dimensions = 2,
          units = "cm",
          operators = {
            { slice = { z = 0 } }
          }
        }
      }
    }
  )",
                                         InputFormat::Lua);

  ASSERT_EQ(2u, shapeSet.getShapes().size());
  EXPECT_EQ(Dimensions::Three, shapeSet.getShapes()[0].getGeometry().getInputDimensions());
  EXPECT_EQ(Dimensions::Three, shapeSet.getShapes()[0].getGeometry().getOutputDimensions());
  EXPECT_EQ(Dimensions::Three, shapeSet.getShapes()[1].getGeometry().getInputDimensions());
  EXPECT_EQ(Dimensions::Two, shapeSet.getShapes()[1].getGeometry().getOutputDimensions());
}

TEST(IOTest, readShapeSet_luaGeneratedOrdinaryTableValues)
{
  auto shapeSet = readShapeSetFromString(R"(
    local dim = 2
    local r = 4.0
    local z = 8.0
    local x = 1.0
    local y = 2.0

    dimensions = dim

    shapes = {
      {
        name = "part",
        material = "steel",
        geometry = {
          format = "stl",
          path = "part.stl",
          units = "cm",
          operators = {
            { translate = (dim == 2) and {r, z} or {x, y, z} }
          }
        }
      }
    }
  )",
                                         InputFormat::Lua);

  ASSERT_EQ(1u, shapeSet.getShapes().size());
  auto composite = std::dynamic_pointer_cast<const CompositeOperator>(
    shapeSet.getShapes()[0].getGeometry().getGeometryOperator());
  ASSERT_TRUE(composite);
  ASSERT_EQ(1u, composite->getOperators().size());
  auto translation = dynamic_cast<const Translation *>(composite->getOperators()[0].get());
  ASSERT_NE(translation, nullptr);
  EXPECT_THAT(translation->getOffset(), AlmostEqVector(Vector3D {4, 8, 0}));
}

TEST(IOTest, readShapeSet_luaUnexpectedGlobalDiagnostic)
{
  try
  {
    readShapeSetFromString(R"(
      dimensions = 2
      unexpected_global = {
        nested_value = 42
      }
      shapes = {}
    )",
                           InputFormat::Lua);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &err)
  {
    ASSERT_EQ(1u, err.getErrors().size());
    EXPECT_EQ(axom::Path {"unexpected_global"}, err.getErrors()[0].path);
    EXPECT_THAT(err.what(), HasSubstr("unexpected_global"));
  }
}

TEST(IOTest, readShapeSet_luaNestedUnexpectedFieldsMatchYamlValidation)
{
  auto shapeSet = readShapeSetFromString(R"(
    dimensions = 2
    shapes = {
      {
        name = "wheel",
        material = "steel",
        extra_shape_value = true,
        geometry = {
          format = "stl",
          path = "wheel.stl",
          extra_geometry_table = {
            nested_value = 42
          }
        }
      }
    }
  )",
                                         InputFormat::Lua);

  ASSERT_EQ(1u, shapeSet.getShapes().size());
  EXPECT_EQ("wheel", shapeSet.getShapes()[0].getName());
}
#endif

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

TEST(IOTest, readShapeSet_geometryOperators_scaleWithCenter)
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
              - scale: [1.5, 2.5]
                center: [10, 20]
    )");
  auto &shapes = shapeSet.getShapes();
  ASSERT_EQ(1u, shapes.size());
  auto &geometryOperator = shapes[0].getGeometry().getGeometryOperator();
  ASSERT_TRUE(geometryOperator);

  auto composite = std::dynamic_pointer_cast<const CompositeOperator>(geometryOperator);
  ASSERT_TRUE(composite);
  ASSERT_EQ(1u, composite->getOperators().size());

  auto scale = dynamic_cast<const Scale *>(composite->getOperators()[0].get());
  ASSERT_NE(scale, nullptr);
  EXPECT_DOUBLE_EQ(1.5, scale->getXFactor());
  EXPECT_DOUBLE_EQ(2.5, scale->getYFactor());
  EXPECT_DOUBLE_EQ(1.0, scale->getZFactor());
  EXPECT_THAT(scale->getCenter(), AlmostEqPoint(Point3D {10, 20, 0}));
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

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  int result = RUN_ALL_TESTS();
  return result;
}
