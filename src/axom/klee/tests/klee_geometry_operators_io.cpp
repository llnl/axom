// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/inlet.hpp"

#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/KleeError.hpp"
#include "axom/klee/io/GeometryOperatorsIO.hpp"

#include "KleeMatchers.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <memory>
#include <unordered_map>

namespace axom
{
namespace klee
{
namespace internal
{
namespace
{
using primal::Point3D;
using primal::Vector3D;
using test::AlmostEqMatrix;
using test::AlmostEqPoint;
using test::AlmostEqSlice;
using test::AlmostEqVector;

using ::testing::HasSubstr;

using OperatorPointer = CompositeOperator::OpPtr;

/**
 * Read geometry operators from a string
 * \param startProperties the starting properties
 * \param input the operators expressed in yaml
 * \param an optional map of named operators
 * \return the operators that were read.
 */
OperatorPointer readOperators(const TransformableGeometryProperties &startProperties,
                              const std::string &input,
                              const NamedOperatorMap &namedOperators = NamedOperatorMap {})
{
  auto reader = std::unique_ptr<inlet::YAMLReader>(new inlet::YAMLReader());
  std::string wrappedInput = "test_list:\n  ";
  wrappedInput += input;
  reader->parseString(wrappedInput);

  sidre::DataStore dataStore;
  inlet::Inlet doc(std::move(reader), dataStore.getRoot());
  GeometryOperatorData::defineSchema(doc.getGlobalContainer(), "test_list", "The test list");

  std::vector<inlet::VerificationError> errors;
  if(!doc.verify(&errors))
  {
    throw KleeError(errors);
  }
  auto opData = doc["test_list"].get<GeometryOperatorData>();
  return opData.makeOperator(startProperties, namedOperators);
}

/**
 * Read named operators from the given input string.
 *
 * \param startingDimensions the starting dimensions
 * \param input the operators expressed in yaml
 * \return the named operators that were read
 */
NamedOperatorMap readNamedOperators(Dimensions startingDimensions, const std::string &input)
{
  auto reader = std::unique_ptr<inlet::YAMLReader>(new inlet::YAMLReader());
  std::string wrappedInput = "op_list:\n";
  wrappedInput += input;
  reader->parseString(wrappedInput);

  sidre::DataStore dataStore;
  inlet::Inlet doc(std::move(reader), dataStore.getRoot());
  NamedOperatorMapData::defineSchema(doc.getGlobalContainer(), "op_list");
  std::vector<inlet::VerificationError> errors;
  if(!doc.verify(&errors))
  {
    throw KleeError(errors);
  }
  auto operators = doc["op_list"].get<NamedOperatorMapData>();
  return operators.makeNamedOperatorMap(startingDimensions);
}

/**
 * Copy an operator from a pointer.
 *
 * \tparam T the expected type of the operator
 * \param ptr the pointer to the operator. Must not be null after a dynamic_cast to T.
 * \return a copy of the operator
 * \throws std::logic error if the pointer is of the wrong type
 */
template <typename T>
T copyOperator(const OperatorPointer &ptr)
{
  auto desired = dynamic_cast<const T *>(ptr.get());
  if(desired == nullptr)
  {
    throw KleeError {inlet::VerificationError {Path {"<no path>"}, "Did not get expected type"}};
  }
  return *desired;
}

/**
 * Get a single operator from a CompositeOperator.
 *
 * \tparam T the expected type of the operator
 * \param ptr the pointer to the operator. Must not be null after a dynamic_cast to CompositeOperator.
 * \return a copy of the operator
 * \throws std::logic error if the pointer is of the wrong type
 */
template <typename T>
T getSingleOperatorFromComposite(const OperatorPointer &ptr)
{
  auto composite = dynamic_cast<const CompositeOperator *>(ptr.get());
  if(composite == nullptr)
  {
    throw KleeError {inlet::VerificationError {Path {"<no path>"}, "Did not get CompositeOperator"}};
  }
  if(composite->getOperators().size() != 1u)
  {
    throw KleeError {
      inlet::VerificationError {Path {"<no path>"}, "Did not have exactly one operator"}};
  }
  return copyOperator<T>(composite->getOperators()[0]);
}

/**
 * Read a single operator.
 *
 * \tparam T the expected type of the operator
 * \param startProperties the starting properties
 * \param input the operator as yaml.
 * \param an optional map of named operators
 * \return the read operator
 * \throws KleeError if there isn't exactly one operator of
 * the specified type.
 */
template <typename T>
T readSingleOperator(const TransformableGeometryProperties &startProperties,
                     const std::string &input,
                     const std::unordered_map<std::string, OperatorPointer> &namedOperators =
                       std::unordered_map<std::string, OperatorPointer> {})
{
  std::string wrappedInput {"-\n"};
  wrappedInput += input;
  OperatorPointer genericOperator = readOperators(startProperties, wrappedInput, namedOperators);
  return getSingleOperatorFromComposite<T>(genericOperator);
}

/**
 * Get a named operator from the given map. The composite operator in the
 * map must have a single entry of the specified type.
 *
 * \tparam T the expected type of the operator
 * \param operators the map from which to get the operator
 * \param name the name of the operator
 * \return the requested operator
 * \throws KleeError if any of the operators are of the wrong type
 */
template <typename T>
T getSingleNamedOperator(const NamedOperatorMap &operators, std::string const &name)
{
  auto iter = operators.find(name);
  if(iter == operators.end())
  {
    std::string message = "No operator named ";
    throw KleeError {inlet::VerificationError {Path {"<no path>"}, message + name}};
  }
  return getSingleOperatorFromComposite<T>(iter->second);
}

/// Helper function to create a klee::SliceOperator from its constituent parts
SliceOperator make_slice(Point3D origin,
                         Vector3D normal,
                         Vector3D up,
                         TransformableGeometryProperties startProperties)
{
  return SliceOperator {origin, normal, up, startProperties};
}

TEST(GeometryOperatorsIO, readMultipleOperatorsIncluded)
{
  try
  {
    readOperators({Dimensions::Two, LengthUnit::cm}, R"(
      - translate: [10, 20]
        rotate: 90
    )");
    FAIL() << "Should have thrown";
  }
  catch(KleeError const &error)
  {
    EXPECT_THAT(error.what(), HasSubstr("translate"));
    EXPECT_THAT(error.what(), HasSubstr("rotate"));
  }
}

TEST(GeometryOperatorsIO, readTranslation_2D)
{
  auto translation = readSingleOperator<Translation>({Dimensions::Two, LengthUnit::cm}, R"(
      translate: [10, 20]
    )");
  EXPECT_THAT(translation.getOffset(), AlmostEqVector(Vector3D {10, 20, 00}));
  TransformableGeometryProperties expectedProperties {Dimensions::Two, LengthUnit::cm};
  EXPECT_EQ(expectedProperties, translation.getStartProperties());
  EXPECT_EQ(expectedProperties, translation.getEndProperties());
}

TEST(GeometryOperatorsIO, readTranslation_3D)
{
  auto translation = readSingleOperator<Translation>({Dimensions::Three, LengthUnit::cm}, R"(
      translate: [10, 20, 30]
    )");
  EXPECT_THAT(translation.getOffset(), AlmostEqVector(Vector3D {10, 20, 30}));
  TransformableGeometryProperties expectedProperties {Dimensions::Three, LengthUnit::cm};
  EXPECT_EQ(expectedProperties, translation.getStartProperties());
  EXPECT_EQ(expectedProperties, translation.getEndProperties());
}

TEST(GeometryOperatorsIO, readTranslation_unknownKeys)
{
  try
  {
    readSingleOperator<Translation>({Dimensions::Two, LengthUnit::cm},
                                    R"(
        translate: [10, 20]
        UNKNOWN_KEY: UNKNOWN_VALUE
      )");
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &err)
  {
    EXPECT_THAT(err.what(), HasSubstr("UNKNOWN_KEY"));
  }
}

TEST(GeometryOperatorsIO, readRotation_2D_requiredOnly)
{
  auto rotation = readSingleOperator<Rotation>({Dimensions::Two, LengthUnit::cm}, R"(
      rotate: 45
    )");
  EXPECT_DOUBLE_EQ(45, rotation.getAngle());
  EXPECT_THAT(rotation.getCenter(), AlmostEqPoint(Point3D {0, 0, 0}));
  EXPECT_THAT(rotation.getAxis(), AlmostEqVector(Vector3D {0, 0, 1}));
  TransformableGeometryProperties expectedProperties {Dimensions::Two, LengthUnit::cm};
  EXPECT_EQ(expectedProperties, rotation.getStartProperties());
  EXPECT_EQ(expectedProperties, rotation.getEndProperties());
}

TEST(GeometryOperatorsIO, readRotation_2D_optionalFields)
{
  auto rotation = readSingleOperator<Rotation>({Dimensions::Two, LengthUnit::cm}, R"(
      rotate: 45
      center: [10, 20]
    )");
  EXPECT_DOUBLE_EQ(45, rotation.getAngle());
  EXPECT_THAT(rotation.getCenter(), AlmostEqPoint(Point3D {10, 20, 0}));
  EXPECT_THAT(rotation.getAxis(), AlmostEqVector(Vector3D {0, 0, 1}));
  TransformableGeometryProperties expectedProperties {Dimensions::Two, LengthUnit::cm};
  EXPECT_EQ(expectedProperties, rotation.getStartProperties());
  EXPECT_EQ(expectedProperties, rotation.getEndProperties());
}

TEST(GeometryOperatorsIO, readRotation_2D_axisNotAllowed)
{
  try
  {
    readSingleOperator<Rotation>({Dimensions::Two, LengthUnit::cm}, R"(
            rotate: 45
            axis: [1, 2, 3]
        )");
    FAIL() << "Should have thrown an exception";
  }
  catch(const KleeError &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("rotate"));
    EXPECT_THAT(ex.what(), HasSubstr("axis"));
  }
}

TEST(GeometryOperatorsIO, readRotation_3D_requiredOnly)
{
  auto rotation = readSingleOperator<Rotation>({Dimensions::Three, LengthUnit::cm}, R"(
      rotate: 45
      axis: [1, 2, 3]
    )");
  EXPECT_DOUBLE_EQ(45, rotation.getAngle());
  EXPECT_THAT(rotation.getCenter(), AlmostEqPoint(Point3D {0, 0, 0}));
  EXPECT_THAT(rotation.getAxis(), AlmostEqVector(Vector3D {1, 2, 3}));
  TransformableGeometryProperties expectedProperties {Dimensions::Three, LengthUnit::cm};
  EXPECT_EQ(expectedProperties, rotation.getStartProperties());
  EXPECT_EQ(expectedProperties, rotation.getEndProperties());
}

TEST(GeometryOperatorsIO, readRotation_3D_optionalFields)
{
  auto rotation = readSingleOperator<Rotation>({Dimensions::Three, LengthUnit::cm}, R"(
      rotate: 45
      axis: [1, 2, 3]
      center: [4, 5, 6]
    )");
  EXPECT_DOUBLE_EQ(45, rotation.getAngle());
  EXPECT_THAT(rotation.getCenter(), AlmostEqPoint(Point3D {4, 5, 6}));
  EXPECT_THAT(rotation.getAxis(), AlmostEqVector(Vector3D {1, 2, 3}));
  TransformableGeometryProperties expectedProperties {Dimensions::Three, LengthUnit::cm};
  EXPECT_EQ(expectedProperties, rotation.getStartProperties());
  EXPECT_EQ(expectedProperties, rotation.getEndProperties());
}

TEST(GeometryOperatorsIO, readRotation_3D_axisMissing)
{
  try
  {
    readOperators({Dimensions::Three, LengthUnit::cm}, R"(
          - rotate: 45
        )");
    FAIL() << "Should not have parsed";
  }
  catch(const KleeError &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("axis"));
  }
}

TEST(GeometryOperatorsIO, readScale_singleValue)
{
  Dimensions all_dims[] = {Dimensions::Two, Dimensions::Three};
  for(Dimensions dims : all_dims)
  {
    auto scale = readSingleOperator<Scale>({dims, LengthUnit::cm}, R"(
          scale: 1.2
        )");
    EXPECT_DOUBLE_EQ(1.2, scale.getXFactor());
    EXPECT_DOUBLE_EQ(1.2, scale.getYFactor());
    EXPECT_DOUBLE_EQ(1.2, scale.getZFactor());
    TransformableGeometryProperties expectedProperties {dims, LengthUnit::cm};
    EXPECT_EQ(expectedProperties, scale.getStartProperties());
    EXPECT_EQ(expectedProperties, scale.getEndProperties());
  }
}

TEST(GeometryOperatorsIO, readScale_2d_array)
{
  auto scale = readSingleOperator<Scale>({Dimensions::Two, LengthUnit::cm}, R"(
      scale: [1.2, 3.4]
    )");
  EXPECT_DOUBLE_EQ(1.2, scale.getXFactor());
  EXPECT_DOUBLE_EQ(3.4, scale.getYFactor());
  EXPECT_DOUBLE_EQ(1.0, scale.getZFactor());
  TransformableGeometryProperties expectedProperties {Dimensions::Two, LengthUnit::cm};
  EXPECT_EQ(expectedProperties, scale.getStartProperties());
  EXPECT_EQ(expectedProperties, scale.getEndProperties());
}

TEST(GeometryOperatorsIO, readScale_3d_array)
{
  auto scale = readSingleOperator<Scale>({Dimensions::Three, LengthUnit::cm},
                                         R"(
      scale: [1.2, 3.4, 5.6]
  )");
  EXPECT_DOUBLE_EQ(1.2, scale.getXFactor());
  EXPECT_DOUBLE_EQ(3.4, scale.getYFactor());
  EXPECT_DOUBLE_EQ(5.6, scale.getZFactor());
  TransformableGeometryProperties expectedProperties {Dimensions::Three, LengthUnit::cm};
  EXPECT_EQ(expectedProperties, scale.getStartProperties());
  EXPECT_EQ(expectedProperties, scale.getEndProperties());
}

TEST(GeometryOperatorsIO, readConvertUnits)
{
  auto converter = readSingleOperator<UnitConverter>({Dimensions::Three, LengthUnit::inches}, R"(
      convert_units_to: cm
  )");
  TransformableGeometryProperties expectedStartProperties {Dimensions::Three, LengthUnit::inches};
  TransformableGeometryProperties expectedEndProperties {Dimensions::Three, LengthUnit::cm};
  EXPECT_EQ(expectedStartProperties, converter.getStartProperties());
  EXPECT_EQ(expectedEndProperties, converter.getEndProperties());
}

TEST(GeometryOperatorsIO, readSlice_specifyAll)
{
  auto slice = readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, R"(
      slice:
        origin: [1, 2, 3]
        normal: [4, 5, 6]
        up: [-5, 4, 0]
    )");
  EXPECT_THAT(
    slice,
    AlmostEqSlice(make_slice({1, 2, 3}, {4, 5, 6}, {-5, 4, 0}, {Dimensions::Three, LengthUnit::cm})));
}

TEST(GeometryOperatorsIO, readSlice_x_defaults)
{
  auto slice = readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, R"(
      slice:
        x: 10
    )");
  EXPECT_THAT(
    slice,
    AlmostEqSlice(make_slice({10, 0, 0}, {1, 0, 0}, {0, 0, 1}, {Dimensions::Three, LengthUnit::cm})));
}

TEST(GeometryOperatorsIO, readSlice_x_optionals)
{
  auto slice = readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, R"(
      slice:
        x: 10
        origin: [10, 20, 30]
        normal: [40, 0, 0]
        up: [0, -50, 60]
    )");
  EXPECT_THAT(
    slice,
    AlmostEqSlice(
      make_slice({10, 20, 30}, {40, 0, 0}, {0, -50, 60}, {Dimensions::Three, LengthUnit::cm})));
}

TEST(GeometryOperatorsIO, readSlice_y_defaults)
{
  auto slice = readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, R"(
      slice:
        y: 20
    )");
  EXPECT_THAT(
    slice,
    AlmostEqSlice(make_slice({0, 20, 0}, {0, 1, 0}, {1, 0, 0}, {Dimensions::Three, LengthUnit::cm})));
}

TEST(GeometryOperatorsIO, readSlice_y_optionals)
{
  auto slice = readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, R"(
      slice:
        y: 20
        origin: [10, 20, 30]
        normal: [0, 40, 0]
        up: [-50, 0, 60]
    )");
  EXPECT_THAT(
    slice,
    AlmostEqSlice(
      make_slice({10, 20, 30}, {0, 40, 0}, {-50, 0, 60}, {Dimensions::Three, LengthUnit::cm})));
}

TEST(GeometryOperatorsIO, readSlice_z_defaults)
{
  auto slice = readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, R"(
      slice:
        z: 30
    )");
  EXPECT_THAT(
    slice,
    AlmostEqSlice(make_slice({0, 0, 30}, {0, 0, 1}, {0, 1, 0}, {Dimensions::Three, LengthUnit::cm})));
}

TEST(GeometryOperatorsIO, readSlice_z_optionals)
{
  auto slice = readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, R"(
      slice:
        z: 30
        origin: [10, 20, 30]
        normal: [0, 0, 40]
        up: [-50, 60, 0]
    )");
  EXPECT_THAT(
    slice,
    AlmostEqSlice(
      make_slice({10, 20, 30}, {0, 0, 40}, {-50, 60, 0}, {Dimensions::Three, LengthUnit::cm})));
}

TEST(GeometryOperatorsIO, readSlice_zeroNormal)
{
  try
  {
    readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, R"(
          slice:
            origin: [10, 20, 30]
            normal: [0, 0, 0]
            up: [0, -50, 60]
        )");
    FAIL() << "Should have thrown a message about the normal being zero";
  }
  catch(KleeError &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("normal"));
    EXPECT_THAT(ex.what(), HasSubstr("zero"));
  }
}

TEST(GeometryOperatorsIO, readSlice_upAndNormalNotNormal)
{
  try
  {
    readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, R"(
          slice:
            origin: [10, 20, 30]
            normal: [10, 20, 30]
            up: [10, 0, 0]
        )");
    FAIL() << "Should have thrown a message about the normal and up "
              "vectors not being perpendicular";
  }
  catch(KleeError &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("normal"));
    EXPECT_THAT(ex.what(), HasSubstr("up"));
  }
}

TEST(GeometryOperatorsIO, readSlice_badPlaneValues)
{
  auto badOriginInput = R"(
      slice:
        x: 10
        origin: [20, 0, 0]
  )";
  EXPECT_THROW(readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, badOriginInput),
               KleeError)
    << "Bad origin";

  auto badNormalInput = R"(
      slice:
        x: 10
        normal: [1, 2, 3]
  )";
  EXPECT_THROW(readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, badNormalInput),
               KleeError)
    << "Bad normal";

  auto badUpInput = R"(
      slice:
        x: 10
        up: [1, 2, 3]
  )";
  EXPECT_THROW(readSingleOperator<SliceOperator>({Dimensions::Three, LengthUnit::cm}, badUpInput), KleeError)
    << "Bad up";
}

TEST(GeometryOperatorsIO, readSlice_start2D)
{
  auto input = R"(
      slice:
        x: 10
        origin: [20, 0, 0]
  )";
  EXPECT_THROW(readSingleOperator<SliceOperator>({Dimensions::Two, LengthUnit::cm}, input),
               KleeError);
}

TEST(GeometryOperatorsIO, readMultiple_matchingDimensions)
{
  auto op = readOperators({Dimensions::Three, LengthUnit::cm}, R"(
      - translate: [10, 20, 30]
      - translate: [40, 50, 60]
    )");
  auto composite = std::dynamic_pointer_cast<const CompositeOperator>(op);
  ASSERT_TRUE(composite);
  EXPECT_EQ(2u, composite->getOperators().size());
}

TEST(GeometryOperatorsIO, readMultiple_nonMatchingDimensions)
{
  auto input = R"(
      - translate: [10, 20, 30]
      - translate: [40, 50]
  )";
  EXPECT_THROW(readOperators({Dimensions::Three, LengthUnit::cm}, input), KleeError);
}

TEST(GeometryOperatorsIO, readMultiple_unknownOperator)
{
  try
  {
    readOperators({Dimensions::Three, LengthUnit::cm}, R"(
          - UNKNOWN_OPERATOR: [10, 20, 30]
         )");
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &error)
  {
    EXPECT_THAT(error.what(), HasSubstr("UNKNOWN_OPERATOR"));
  }
}

TEST(GeometryOperatorsIO, readRef_missing)
{
  try
  {
    readOperators({Dimensions::Two, LengthUnit::cm}, R"(
       - ref: MISSING
    )");
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("MISSING"));
  }
}

TEST(GeometryOperatorsIO, readRef_present)
{
  TransformableGeometryProperties props {Dimensions::Three, LengthUnit::cm};
  NamedOperatorMap namedOperators = {
    {"op1", std::make_shared<Translation>(Vector3D {10, 20, 30}, props)},
    {"op2", std::make_shared<Scale>(1.5, 2.0, 2.5, props)},
  };
  auto translation = readSingleOperator<Translation>({Dimensions::Three, LengthUnit::cm},
                                                     R"(
      ref: op1
    )",
                                                     namedOperators);
  EXPECT_THAT(translation.getOffset(), AlmostEqVector(Vector3D {10, 20, 30}));
}

class RefIoUnitsMismatchTest : public ::testing::Test
{
protected:
  void SetUp() override;
  std::shared_ptr<GeometryOperator> referencedOperator;
  CompositeOperator readOperator(const TransformableGeometryProperties &startProperties) const;
};

void RefIoUnitsMismatchTest::SetUp()
{
  TransformableGeometryProperties startProps {Dimensions::Three, LengthUnit::cm};
  auto composite = std::make_shared<CompositeOperator>(startProps);
  composite->addOperator(
    std::make_shared<Translation>(Vector3D {10, 20, 30}, composite->getEndProperties()));
  composite->addOperator(
    std::make_shared<UnitConverter>(LengthUnit::mm, composite->getEndProperties()));
  composite->addOperator(
    std::make_shared<Translation>(Vector3D {40, 50, 60}, composite->getEndProperties()));
  referencedOperator = composite;
}

CompositeOperator RefIoUnitsMismatchTest::readOperator(
  const TransformableGeometryProperties &startProperties) const
{
  NamedOperatorMap namedOperators = {
    {"op1", referencedOperator},
  };

  return readSingleOperator<CompositeOperator>(startProperties,
                                               R"(
      ref: op1
  )",
                                               namedOperators);
}

TEST_F(RefIoUnitsMismatchTest, startUnitMismatch)
{
  TransformableGeometryProperties startProperties {Dimensions::Three, LengthUnit::mm};
  auto composite = readOperator(startProperties);

  EXPECT_EQ(2u, composite.getOperators().size());

  auto automaticInitialConversion = copyOperator<UnitConverter>(composite.getOperators()[0]);
  EXPECT_EQ(LengthUnit::mm, automaticInitialConversion.getStartProperties().units);
  EXPECT_EQ(LengthUnit::cm, automaticInitialConversion.getEndProperties().units);

  EXPECT_EQ(referencedOperator, composite.getOperators()[1]);
}

TEST_F(RefIoUnitsMismatchTest, endUnitMismatch)
{
  TransformableGeometryProperties startProperties {Dimensions::Three, LengthUnit::cm};
  auto composite = readOperator(startProperties);

  EXPECT_EQ(2u, composite.getOperators().size());

  EXPECT_EQ(referencedOperator, composite.getOperators()[0]);

  auto automaticFinalConversion = copyOperator<UnitConverter>(composite.getOperators()[1]);
  EXPECT_EQ(LengthUnit::mm, automaticFinalConversion.getStartProperties().units);
  EXPECT_EQ(LengthUnit::cm, automaticFinalConversion.getEndProperties().units);
}

TEST_F(RefIoUnitsMismatchTest, allUnitMismatch)
{
  TransformableGeometryProperties startProperties {Dimensions::Three, LengthUnit::inches};
  auto composite = readOperator(startProperties);

  EXPECT_EQ(3u, composite.getOperators().size());

  auto automaticInitialConversion = copyOperator<UnitConverter>(composite.getOperators()[0]);
  EXPECT_EQ(LengthUnit::inches, automaticInitialConversion.getStartProperties().units);
  EXPECT_EQ(LengthUnit::cm, automaticInitialConversion.getEndProperties().units);

  EXPECT_EQ(referencedOperator, composite.getOperators()[1]);

  auto automaticFinalConversion = copyOperator<UnitConverter>(composite.getOperators()[2]);
  EXPECT_EQ(LengthUnit::mm, automaticFinalConversion.getStartProperties().units);
  EXPECT_EQ(LengthUnit::inches, automaticFinalConversion.getEndProperties().units);
}

TEST(GeometryOperatorsIO, readNamedOperators_basic)
{
  NamedOperatorMap readOperators = readNamedOperators(Dimensions::Two, R"(
      - name: op1
        units: cm
        value:
          - translate: [10, 20]
      - name: op2
        units: in
        value:
          - scale: 1.5
    )");
  EXPECT_EQ(2, readOperators.size());

  auto translation = getSingleNamedOperator<Translation>(readOperators, "op1");
  TransformableGeometryProperties expectedTranslationProperties {Dimensions::Two, LengthUnit::cm};
  EXPECT_EQ(expectedTranslationProperties, translation.getStartProperties());
  EXPECT_EQ(expectedTranslationProperties, translation.getEndProperties());
  EXPECT_THAT(translation.getOffset(), AlmostEqVector(Vector3D {10, 20, 0}));

  auto scale = getSingleNamedOperator<Scale>(readOperators, "op2");
  TransformableGeometryProperties expectedScaleProperties {Dimensions::Two, LengthUnit::inches};
  EXPECT_EQ(expectedScaleProperties, scale.getStartProperties());
  EXPECT_EQ(expectedScaleProperties, scale.getEndProperties());
  EXPECT_EQ(1.5, scale.getXFactor());
  EXPECT_EQ(1.5, scale.getYFactor());
}

TEST(GeometryOperatorsIO, readNamedOperators_invalidDimensions)
{
  auto input = R"(
      - name: op1
        units: cm
        value:
          - translate: [10, 20, 30]
  )";
  EXPECT_THROW(readNamedOperators(Dimensions::Two, input), KleeError);
}

TEST(GeometryOperatorsIO, readNamedOperators_differentInitialDimensions)
{
  NamedOperatorMap readOperators = readNamedOperators(Dimensions::Two, R"(
      - name: op1
        units: cm
        start_dimensions: 3
        value:
          - translate: [10, 20, 30]
    )");
  ASSERT_EQ(1, readOperators.size());

  auto translation = getSingleNamedOperator<Translation>(readOperators, "op1");
  EXPECT_THAT(translation.getOffset(), AlmostEqVector(Vector3D {10, 20, 30}));
  TransformableGeometryProperties expectedProperties {Dimensions::Three, LengthUnit::cm};
  EXPECT_EQ(expectedProperties, translation.getStartProperties());
  EXPECT_EQ(expectedProperties, translation.getEndProperties());
}

TEST(GeometryOperatorsIO, readNamedOperators_differentStartAndEnd)
{
  NamedOperatorMap readOperators = readNamedOperators(Dimensions::Three, R"(
      - name: op1
        start_units: mm
        end_units: cm
        value:
          - translate: [10, 20, 30]
          - convert_units_to: cm
          - translate: [40, 50, 60]
    )");
  ASSERT_EQ(1, readOperators.size());

  auto composite = copyOperator<CompositeOperator>(readOperators["op1"]);
  ASSERT_EQ(3u, composite.getOperators().size());
  TransformableGeometryProperties expectedStartProperties {Dimensions::Three, LengthUnit::mm};
  TransformableGeometryProperties expectedEndProperties {Dimensions::Three, LengthUnit::cm};
  EXPECT_EQ(expectedStartProperties, composite.getStartProperties());
  EXPECT_EQ(expectedEndProperties, composite.getEndProperties());

  auto firstTranslation = copyOperator<Translation>(composite.getOperators()[0]);
  EXPECT_THAT(firstTranslation.getOffset(), AlmostEqVector(Vector3D {10, 20, 30}));
  EXPECT_EQ(expectedStartProperties, firstTranslation.getStartProperties());
  EXPECT_EQ(expectedStartProperties, firstTranslation.getEndProperties());

  auto conversion = copyOperator<UnitConverter>(composite.getOperators()[1]);
  EXPECT_EQ(expectedStartProperties, conversion.getStartProperties());
  EXPECT_EQ(expectedEndProperties, conversion.getEndProperties());

  auto secondTranslation = copyOperator<Translation>(composite.getOperators()[2]);
  EXPECT_THAT(secondTranslation.getOffset(), AlmostEqVector(Vector3D {40, 50, 60}));
  EXPECT_EQ(expectedEndProperties, secondTranslation.getStartProperties());
  EXPECT_EQ(expectedEndProperties, secondTranslation.getEndProperties());
}

TEST(GeometryOperatorsIO, readNamedOperators_errorInEndUnits)
{
  try
  {
    readNamedOperators(Dimensions::Three, R"(
      - name: op1
        start_units: mm
        end_units: in
        value:
          - translate: [10, 20, 30]
          - convert_units_to: cm
          - translate: [40, 50, 60]
    )");
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("units"));
  }
}

TEST(GeometryOperatorsIO, readNamedOperators_ref)
{
  NamedOperatorMap readOperators = readNamedOperators(Dimensions::Two, R"(
      - name: op1
        units: in
        value:
          - translate: [10, 20]
      - name: op2
        units: in
        value:
          - scale: 1.5
          - ref: op1
    )");
  EXPECT_EQ(2, readOperators.size());

  auto translation = getSingleNamedOperator<Translation>(readOperators, "op1");
  EXPECT_THAT(translation.getOffset(), AlmostEqVector(Vector3D {10, 20, 0}));

  auto op2 = readOperators["op2"];
  auto composite = dynamic_cast<const CompositeOperator *>(op2.get());
  ASSERT_NE(nullptr, composite);
  TransformableGeometryProperties expectedOp2Units {Dimensions::Two, LengthUnit::inches};
  EXPECT_EQ(expectedOp2Units, composite->getStartProperties());
  EXPECT_EQ(expectedOp2Units, composite->getEndProperties());
  ASSERT_EQ(2u, composite->getOperators().size());

  auto scale = copyOperator<Scale>(composite->getOperators()[0]);
  EXPECT_EQ(1.5, scale.getXFactor());
  EXPECT_EQ(1.5, scale.getYFactor());

  auto referencedOperator = copyOperator<CompositeOperator>(composite->getOperators()[1]);
  ASSERT_EQ(1u, referencedOperator.getOperators().size());
  auto referencedTranslation = copyOperator<Translation>(referencedOperator.getOperators()[0]);
  EXPECT_THAT(referencedTranslation.getOffset(), AlmostEqVector(Vector3D {10, 20, 0}));
}

}  // namespace
}  // namespace internal
}  // namespace klee
}  // namespace axom

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  int result = RUN_ALL_TESTS();
  return result;
}
