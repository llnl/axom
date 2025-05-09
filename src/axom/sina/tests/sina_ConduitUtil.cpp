// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <stdexcept>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/core/ConduitUtil.hpp"
#include "axom/sina/tests/SinaMatchers.hpp"

namespace axom
{
namespace sina
{
namespace testing
{
namespace
{

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::ElementsAre;
using ::testing::HasSubstr;

TEST(ConduitUtil, getRequiredField_present)
{
  conduit::Node parent;
  parent["fieldName"] = "field value";
  auto &field = getRequiredField("fieldName", parent, "parent name");
  EXPECT_TRUE(field.dtype().is_string());
  EXPECT_EQ("field value", field.as_string());
}

TEST(ConduitUtil, getRequiredField_missing)
{
  conduit::Node parent;
  try
  {
    auto &field = getRequiredField("fieldName", parent, "parent name");
    FAIL() << "Should not have found field, but got " << field.name();
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr("fieldName"));
    EXPECT_THAT(expected.what(), HasSubstr("parent name"));
  }
}

TEST(ConduitUtil, getRequiredField_slashes)
{
  conduit::Node parent;
  // Conduit by default parses /, creating parent["some"]["name"]
  parent["some/name"] = 24;
  // This is how we provide a literal name with slashes
  parent.add_child("some/name") = 42;
  EXPECT_EQ(42, getRequiredField("some/name", parent, "parent name").to_int64());
}

TEST(ConduitUtil, getRequiredString_valid)
{
  conduit::Node parent;
  parent["fieldName"] = "field value";
  EXPECT_EQ("field value", getRequiredString("fieldName", parent, "parent name"));
}

TEST(ConduitUtil, getRequiredString_missing)
{
  conduit::Node parent;
  try
  {
    auto value = getRequiredString("fieldName", parent, "parent name");
    FAIL() << "Should not have found string, but got " << value;
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr("fieldName"));
    EXPECT_THAT(expected.what(), HasSubstr("parent name"));
  }
}

TEST(ConduitUtil, getRequiredString_wrongType)
{
  conduit::Node parent;
  parent["fieldName"] = 123;
  try
  {
    auto value = getRequiredString("fieldName", parent, "parent name");
    FAIL() << "Should not have found string, but got " << value;
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr("fieldName"));
    EXPECT_THAT(expected.what(), HasSubstr("parent name"));
    EXPECT_THAT(expected.what(), HasSubstr("string"));
  }
}

TEST(ConduitUtil, getRequiredString_slashes)
{
  conduit::Node parent;
  parent["some/name"] = "undesired value";
  parent.add_child("some/name") = "desired value";
  EXPECT_EQ("desired value", getRequiredString("some/name", parent, "parent name"));
}

TEST(ConduitUtil, getRequiredDouble_valid)
{
  conduit::Node parent;
  parent["fieldName"] = 3.14;
  EXPECT_THAT(3.14, DoubleEq(getRequiredDouble("fieldName", parent, "parent name")));
}

TEST(ConduitUtil, getRequiredDouble_missing)
{
  conduit::Node parent;
  try
  {
    auto value = getRequiredDouble("fieldName", parent, "parent name");
    FAIL() << "Should not have found double, but got " << value;
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr("fieldName"));
    EXPECT_THAT(expected.what(), HasSubstr("parent name"));
  }
}

TEST(ConduitUtil, getRequiredDouble_wrongType)
{
  conduit::Node parent;
  parent["fieldName"] = "field value";
  try
  {
    auto value = getRequiredDouble("fieldName", parent, "parent name");
    FAIL() << "Should not have found double, but got " << value;
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr("fieldName"));
    EXPECT_THAT(expected.what(), HasSubstr("parent name"));
    EXPECT_THAT(expected.what(), HasSubstr("double"));
  }
}

TEST(ConduitUtil, getOptionalString_valid)
{
  conduit::Node parent;
  parent["fieldName"] = "the value";
  EXPECT_EQ("the value", getOptionalString("fieldName", parent, "parent name"));
}

TEST(ConduitUtil, getOptionalString_missing)
{
  conduit::Node parent;
  EXPECT_EQ("", getOptionalString("fieldName", parent, "parent name"));
}

TEST(ConduitUtil, getOptionalString_explicitNullValue)
{
  conduit::Node parent;
  parent["fieldName"];
  EXPECT_EQ("", getOptionalString("fieldName", parent, "parent name"));
}

TEST(ConduitUtil, getOptionalString_wrongType)
{
  conduit::Node parent;
  parent["fieldName"] = 123;
  try
  {
    auto value = getOptionalString("fieldName", parent, "parent name");
    FAIL() << "Should not have found string, but got " << value;
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr("fieldName"));
    EXPECT_THAT(expected.what(), HasSubstr("parent name"));
    EXPECT_THAT(expected.what(), HasSubstr("string"));
  }
}

TEST(ConduitUtil, getOptionalField_slashes)
{
  conduit::Node parent;
  parent["some/name"] = "undesired value";
  EXPECT_EQ("", getOptionalString("some/name", parent, "parent name"));
}

TEST(ConduitUtil, toDoubleVector_empty)
{
  conduit::Node emptyList = parseJsonValue("[]");
  EXPECT_EQ(std::vector<double> {}, toDoubleVector(emptyList, "testNode"));
}

TEST(ConduitUtil, toDoubleVector_validValues)
{
  conduit::Node nonEmptyList = parseJsonValue("[1.0, 2.0, 3.0, 4.0]");
  EXPECT_THAT(toDoubleVector(nonEmptyList, "testNode"), ElementsAre(1.0, 2.0, 3.0, 4.0));
}

TEST(ConduitUtil, toDoubleVector_NotList)
{
  conduit::Node notList = parseJsonValue("\"this is not a list of doubles\"");
  try
  {
    toDoubleVector(notList, "someName");
    FAIL() << "Should have thrown an exception";
  }
  catch(std::invalid_argument const &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("someName"));
  }
}

TEST(ConduitUtil, toStringVector_empty)
{
  conduit::Node emptyList = parseJsonValue("[]");
  EXPECT_THAT(toStringVector(emptyList, "testNode"), ContainerEq(std::vector<std::string> {}));
}

TEST(ConduitUtil, toStringVector_validValues)
{
  conduit::Node nonEmptyList = parseJsonValue(R"(["s1", "s2", "s3"])");
  EXPECT_THAT(toStringVector(nonEmptyList, "testNode"), ElementsAre("s1", "s2", "s3"));
}

TEST(ConduitUtil, toStringVector_NotList)
{
  conduit::Node notList = parseJsonValue("\"this is not a list of doubles\"");
  try
  {
    toStringVector(notList, "someName");
    FAIL() << "Should have thrown an exception.";
  }
  catch(std::invalid_argument const &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("someName"));
  }
}

TEST(ConduitUtil, toStringVector_NotListOfStrings)
{
  conduit::Node notList = parseJsonValue(R"([1, 2, "a string"])");
  try
  {
    toStringVector(notList, "someName");
    FAIL() << "Should have thrown an exception.";
  }
  catch(std::invalid_argument const &ex)
  {
    EXPECT_THAT(ex.what(), HasSubstr("someName"));
  }
}

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
