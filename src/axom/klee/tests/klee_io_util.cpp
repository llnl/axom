// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/inlet.hpp"

#include "axom/klee/io/IOUtil.hpp"
#include "axom/klee/KleeError.hpp"

#include "KleeMatchers.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <memory>

namespace axom
{
namespace klee
{
namespace internal
{
using primal::Point3D;
using primal::Vector3D;

using test::AlmostEqPoint;
using test::AlmostEqVector;
using ::testing::ElementsAre;

static std::unique_ptr<inlet::Reader> readYaml(const std::string &input)
{
  auto reader = std::unique_ptr<inlet::YAMLReader>(new inlet::YAMLReader());
  reader->parseString(input);
  return reader;
}

class InletTestData
{
public:
  template <typename DefOp>
  InletTestData(const std::string &input, DefOp defOp);

private:
  sidre::DataStore m_store;

public:
  inlet::Inlet doc;
};

template <typename DefOp>
InletTestData::InletTestData(const std::string &input, DefOp defOp)
  : m_store {}
  , doc {readYaml(input), m_store.getRoot()}
{
  defOp(doc.getGlobalContainer());
  // Typically this would be done by the main parsing functions. Since
  // we're testing small pieces out of context, we need to check here
  // if everything is fine.
  std::vector<inlet::VerificationError> errors;
  if(!doc.verify(&errors))
  {
    throw KleeError(errors);
  }
}

std::vector<double> parseDoubleVector(const std::string &vectorInput, Dimensions dims)
{
  std::string fullInput = "values: ";
  fullInput += vectorInput;
  InletTestData data {fullInput, [](inlet::Container &c) { c.addDoubleArray("values"); }};
  return toDoubleVector(data.doc["values"], dims, "values");
}

TEST(io_util, toDoubleVector)
{
  EXPECT_THAT(parseDoubleVector("[1.2, 3.4]", Dimensions::Two), ElementsAre(1.2, 3.4));
  EXPECT_THAT(parseDoubleVector("[1, 2]", Dimensions::Two), ElementsAre(1.0, 2.0));
  EXPECT_THROW(parseDoubleVector("[1, 2]", Dimensions::Three), KleeError) << "Wrong length";
  EXPECT_THROW(parseDoubleVector("[a, b]", Dimensions::Three), KleeError) << "Wrong type";
}

Dimensions defineAndParseDimension(const char *input)
{
  std::string fullInput = "dims: ";
  fullInput += input;
  InletTestData data {fullInput, [](inlet::Container &c) {
                        defineDimensionsField(c, "dims", "some description");
                      }};
  return toDimensions(data.doc["dims"]);
}

TEST(io_util, defineAndConvertDimensions)
{
  EXPECT_EQ(Dimensions::Two, defineAndParseDimension("2"));
  EXPECT_EQ(Dimensions::Three, defineAndParseDimension("3"));
  EXPECT_THROW(defineAndParseDimension("4"), KleeError);
}

/**
 * Wrapper around defineUnitsSchema which calls it with default descriptions
 * to make it easy to test.
 *
 * @param container the Container on which to define the units fields
 */
void defineUnitsSchemaWithDefaults(inlet::Container &container) { defineUnitsSchema(container); }

TEST(io_util, getOptionalStartAndEndUnits_nothingSpecified)
{
  // Random input or Inlet issues a warning about blank input
  InletTestData data {"foo: 123", defineUnitsSchemaWithDefaults};
  auto units = getOptionalStartAndEndUnits(data.doc.getGlobalContainer());
  EXPECT_EQ(LengthUnit::unspecified, std::get<0>(units));
  EXPECT_EQ(LengthUnit::unspecified, std::get<1>(units));
}

TEST(io_util, getOptionalStartAndEndUnits_unitsSpecified)
{
  InletTestData data {"units: cm", defineUnitsSchemaWithDefaults};
  auto units = getOptionalStartAndEndUnits(data.doc.getGlobalContainer());
  EXPECT_EQ(LengthUnit::cm, std::get<0>(units));
  EXPECT_EQ(LengthUnit::cm, std::get<1>(units));
}

TEST(io_util, getOptionalStartAndEndUnits_startAndEndSpecified)
{
  InletTestData data {R"(
        start_units: cm
        end_units: in
    )",
                      defineUnitsSchemaWithDefaults};
  auto units = getOptionalStartAndEndUnits(data.doc.getGlobalContainer());
  EXPECT_EQ(LengthUnit::cm, std::get<0>(units));
  EXPECT_EQ(LengthUnit::inches, std::get<1>(units));
}

TEST(io_util, getOptionalStartAndEndUnits_partialSpecification)
{
  InletTestData startOnly {"start_units: cm", defineUnitsSchemaWithDefaults};
  EXPECT_THROW(getOptionalStartAndEndUnits(startOnly.doc.getGlobalContainer()), KleeError);
  InletTestData endOnly {"end_units: cm", defineUnitsSchemaWithDefaults};
  EXPECT_THROW(getOptionalStartAndEndUnits(endOnly.doc.getGlobalContainer()), KleeError);
}

TEST(io_util, getOptionalStartAndEndUnits_startEndAndUnits)
{
  InletTestData data {R"(
        start_units: cm
        end_units: cm
        units: cm
    )",
                      defineUnitsSchemaWithDefaults};
  EXPECT_THROW(getOptionalStartAndEndUnits(data.doc.getGlobalContainer()), KleeError);
}

TEST(io_util, getStartAndEndUnits_unitsPresent)
{
  InletTestData data {"units: cm", defineUnitsSchemaWithDefaults};
  auto units = getStartAndEndUnits(data.doc.getGlobalContainer());
  EXPECT_EQ(LengthUnit::cm, std::get<0>(units));
  EXPECT_EQ(LengthUnit::cm, std::get<1>(units));
}

TEST(io_util, getStartAndEndUnits_nothingSpecified)
{
  // Random input or Inlet issues a warning about blank input
  InletTestData data {"foo: 123", defineUnitsSchemaWithDefaults};
  EXPECT_THROW(getStartAndEndUnits(data.doc.getGlobalContainer()), KleeError);
}

template <typename T, typename Op>
T parseArray(const char *value, Dimensions dims, Op op)
{
  std::string input = "value: ";
  input += value;
  InletTestData data {input, [](inlet::Container &c) { c.addDoubleArray("value"); }};
  return op(data.doc.getGlobalContainer(), "value", dims);
}

template <typename T, typename Op>
T parseArray(const char *value, Dimensions dims, const T &defaultValue, Op op)
{
  std::string input;
  if(value != nullptr)
  {
    input = "value: ";
    input += value;
  }
  else
  {
    // avoid warning about empty input
    input = "foo: bar";
  }
  InletTestData data {input, [](inlet::Container &c) { c.addDoubleArray("value"); }};
  return op(data.doc.getGlobalContainer(), "value", dims, defaultValue);
}

Point3D parsePoint(const char *value, Dimensions dims)
{
  return parseArray<Point3D>(
    value,
    dims,
    static_cast<Point3D (*)(const inlet::Container &, char const *, Dimensions)>(toPoint));
}

Point3D parsePoint(const char *value, Dimensions dims, Point3D defaultValue)
{
  return parseArray<Point3D>(
    value,
    dims,
    defaultValue,
    static_cast<Point3D (*)(const inlet::Container &, char const *, Dimensions, const Point3D &)>(
      toPoint));
}

Vector3D parseVector(const char *value, Dimensions dims)
{
  return parseArray<Vector3D>(
    value,
    dims,
    static_cast<Vector3D (*)(const inlet::Container &, char const *, Dimensions)>(toVector));
}

Vector3D parseVector(const char *value, Dimensions dims, Vector3D defaultValue)
{
  return parseArray<Vector3D>(
    value,
    dims,
    defaultValue,
    static_cast<Vector3D (*)(const inlet::Container &, char const *, Dimensions, const Vector3D &)>(
      toVector));
}

TEST(io_util, toPoint)
{
  EXPECT_THAT(parsePoint("[1, 2]", Dimensions::Two), AlmostEqPoint(Point3D {1, 2, 0}));
  EXPECT_THAT(parsePoint("[1, 2, 3]", Dimensions::Three), AlmostEqPoint(Point3D {1, 2, 3}));
  EXPECT_THROW(parsePoint("[1, 2]", Dimensions::Three), KleeError);
  EXPECT_THROW(parsePoint("[1, 2, 3]", Dimensions::Two), KleeError);
}

TEST(io_util, toVector)
{
  EXPECT_THAT(parseVector("[1, 2]", Dimensions::Two), AlmostEqVector(Vector3D {1, 2, 0}));
  EXPECT_THAT(parseVector("[1, 2, 3]", Dimensions::Three), AlmostEqVector(Vector3D {1, 2, 3}));
  EXPECT_THROW(parseVector("[1, 2]", Dimensions::Three), KleeError);
  EXPECT_THROW(parseVector("[1, 2, 3]", Dimensions::Two), KleeError);
}

TEST(io_util, toPoint_default)
{
  EXPECT_THAT(parsePoint("[1, 2]", Dimensions::Two, Point3D {4, 5, 6}),
              AlmostEqPoint(Point3D {1, 2, 0}));
  EXPECT_THAT(parsePoint("[1, 2, 3]", Dimensions::Three, Point3D {4, 5, 6}),
              AlmostEqPoint(Point3D {1, 2, 3}));
  EXPECT_THAT(parsePoint(nullptr, Dimensions::Two, Point3D {4, 5, 0}),
              AlmostEqPoint(Point3D {4, 5, 0}));
  EXPECT_THAT(parsePoint(nullptr, Dimensions::Three, Point3D {4, 5, 6}),
              AlmostEqPoint(Point3D {4, 5, 6}));
}

TEST(io_util, toVector_default)
{
  EXPECT_THAT(parseVector("[1, 2]", Dimensions::Two, Vector3D {4, 5, 6}),
              AlmostEqVector(Vector3D {1, 2, 0}));
  EXPECT_THAT(parseVector("[1, 2, 3]", Dimensions::Three, Vector3D {4, 5, 6}),
              AlmostEqVector(Vector3D {1, 2, 3}));
  EXPECT_THAT(parseVector(nullptr, Dimensions::Two, Vector3D {4, 5, 0}),
              AlmostEqVector(Vector3D {4, 5, 0}));
  EXPECT_THAT(parseVector(nullptr, Dimensions::Three, Vector3D {4, 5, 6}),
              AlmostEqVector(Vector3D {4, 5, 6}));
}

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
