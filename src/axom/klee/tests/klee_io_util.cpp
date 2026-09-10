// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/inlet.hpp"

#include "axom/klee/io/IOUtil.hpp"
#include "axom/klee/KleeError.hpp"

#include "gtest/gtest.h"

#include <memory>

namespace axom
{
namespace klee
{
namespace internal
{
static std::unique_ptr<inlet::Reader> readYaml(const std::string& input)
{
  auto reader = std::unique_ptr<inlet::YAMLReader>(new inlet::YAMLReader());
  reader->parseString(input);
  return reader;
}

class InletTestData
{
public:
  template <typename DefOp>
  InletTestData(const std::string& input, DefOp defOp);

private:
  sidre::DataStore m_store;

public:
  inlet::Inlet doc;
};

template <typename DefOp>
InletTestData::InletTestData(const std::string& input, DefOp defOp)
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

Dimensions defineAndParseDimension(const char* input)
{
  std::string fullInput = "dims: ";
  fullInput += input;
  InletTestData data {fullInput, [](inlet::Container& c) {
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
void defineUnitsSchemaWithDefaults(inlet::Container& container) { defineUnitsSchema(container); }

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

}  // namespace internal
}  // namespace klee
}  // namespace axom

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  int result = RUN_ALL_TESTS();
  return result;
}
