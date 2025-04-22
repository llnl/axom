// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/tests/inlet_test_utils.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

using axom::inlet::Container;
using axom::inlet::Field;
using axom::inlet::Inlet;
using axom::inlet::InletType;
using axom::inlet::Proxy;
using axom::inlet::VerificationError;

using ::testing::Contains;
using ::testing::Truly;

template <typename InletReader>
Inlet createBasicInlet(const std::string& luaString, bool enableDocs = true)
{
  std::unique_ptr<InletReader> reader(new InletReader());
  reader->parseString(axom::inlet::detail::fromLuaTo<InletReader>(luaString));

  return Inlet(std::move(reader), enableDocs);
}

template <typename InletReader>
class inlet_Inlet_basic : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_Inlet_basic, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_Inlet_basic, getTopLevelBools)
{
  std::string testString = "foo = true; bar = false";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addBool("foo", "foo's description");
  inlet.addBool("bar", "bar's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addBool("non/existant", "nothing");

  //
  // Check stored values from get
  //

  bool value = false;

  // Check for existing fields
  value = inlet["foo"];
  EXPECT_TRUE(value);

  value = inlet.get<bool>("foo");
  EXPECT_TRUE(value);

  value = inlet["bar"];
  EXPECT_FALSE(value);

  value = inlet.get<bool>("bar");
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto container = inlet["non/existant"];
  EXPECT_EQ(container.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getNestedBools)
{
  std::string testString = "foo = { bar = true; baz = false }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");
  inlet.addBool("foo/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addBool("foo/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  bool value = false;

  // Check for existing fields
  value = inlet["foo/bar"];
  EXPECT_TRUE(value);

  value = inlet["foo/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto container = inlet["foo/nonexistant"];
  EXPECT_EQ(container.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getDoublyNestedBools)
{
  std::string testString = "foo = { quux = { bar = true; baz = false } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addBool("foo/quux/bar", "bar's description");
  inlet.addBool("foo/quux/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addBool("foo/quux/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  bool value = false;

  // Check for existing fields
  value = inlet["foo/quux/bar"];
  EXPECT_TRUE(value);

  value = inlet["foo/quux/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto container = inlet["foo/quux/nonexistant"];
  EXPECT_EQ(container.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getDeeplyNestedBools)
{
  std::string testString =
    "foo = { quux = { corge = { quuz = { grault = { bar = true; baz = false } "
    "} } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addBool("foo/quux/corge/quuz/grault/bar", "bar's description");
  inlet.addBool("foo/quux/corge/quuz/grault/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addBool("foo/quux/corge/quuz/grault/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  bool value = false;

  // Check for existing fields
  value = inlet["foo/quux/corge/quuz/grault/bar"];
  EXPECT_TRUE(value);

  value = inlet["foo/quux/corge/quuz/grault/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto container = inlet["foo/quux/corge/quuz/grault/nonexistant"];
  EXPECT_EQ(container.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getNestedBoolsThroughContainer)
{
  std::string testString = "foo = { bar = true; baz = false }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");
  inlet.addBool("foo/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addBool("foo/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  bool value = false;

  // Grab the subcontainer
  auto foo = inlet["foo"];

  // Check for existing fields
  value = foo["bar"];
  EXPECT_TRUE(value);

  value = foo["baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = foo["nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getDeeplyNestedBoolsThroughContainer)
{
  std::string testString =
    "foo = { quux = { corge = { quuz = { grault = { bar = true; baz = false } "
    "} } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addBool("foo/quux/corge/quuz/grault/bar", "bar's description");
  inlet.addBool("foo/quux/corge/quuz/grault/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addBool("foo/quux/corge/quuz/grault/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  bool value = false;

  auto container = inlet["foo/quux/corge"];

  // Check for existing fields
  value = container["quuz/grault/bar"];
  EXPECT_TRUE(value);

  value = container["quuz/grault/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = container["quuz/grault/nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getDeeplyNestedBoolsThroughField)
{
  std::string testString =
    "foo = { quux = { corge = { quuz = { grault = { bar = true; baz = false } "
    "} } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addBool("foo/quux/corge/quuz/grault/bar", "bar's description");
  inlet.addBool("foo/quux/corge/quuz/grault/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addBool("foo/quux/corge/quuz/grault/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  bool value = false;

  auto container = inlet["foo/quux/corge"];

  // Check for existing fields
  value = container["quuz/grault/bar"];
  EXPECT_TRUE(value);

  value = container["quuz/grault/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto nonexistant_field = container["quuz/grault/nonexistant"];
  EXPECT_EQ(nonexistant_field.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getTopLevelDoubles)
{
  std::string testString = "foo = 5.05; bar = 15.1";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addDouble("foo", "foo's description");
  inlet.addDouble("bar", "bar's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addDouble("nonexistant", "nothing");

  //
  // Check stored values from get
  //

  double value = -1;

  // Check for existing fields
  value = inlet["foo"];
  EXPECT_EQ(value, 5.05);

  value = inlet["bar"];
  EXPECT_EQ(value, 15.1);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = inlet["nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getNestedDoubles)
{
  std::string testString = "foo = { bar = 200.5; baz = 100.987654321 }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addDouble("foo/bar", "bar's description");
  inlet.addDouble("foo/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addDouble("foo/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  double value = -1;

  // Check for existing fields
  value = inlet["foo/bar"];
  EXPECT_EQ(value, 200.5);

  value = inlet["foo/baz"];
  EXPECT_EQ(value, 100.987654321);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = inlet["foo/nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getTopLevelInts)
{
  std::string testString = "foo = 5; bar = 15";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addInt("foo", "foo's description");
  inlet.addInt("bar", "bar's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addInt("nonexistant", "nothing");

  //
  // Check stored values from get
  //

  int value = -1;

  // Check for existing fields
  value = inlet["foo"];
  EXPECT_EQ(value, 5);

  value = inlet["bar"];
  EXPECT_EQ(value, 15);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = inlet["nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getNestedInts)
{
  std::string testString = "foo = { bar = 200; baz = 100 }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addInt("foo/bar", "bar's description");
  inlet.addInt("foo/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addInt("foo/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  int value = -1;

  // Check for existing fields
  value = inlet["foo/bar"];
  EXPECT_EQ(value, 200);

  value = inlet["foo/baz"];
  EXPECT_EQ(value, 100);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = inlet["foo/nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getTopLevelStrings)
{
  std::string testString = "foo = 'test string'; bar = 'other test string'";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addString("foo", "foo's description");
  inlet.addString("bar", "bar's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addString("nonexistant", "nothing");

  //
  // Check stored values from get
  //

  std::string value = "";

  // Check for existing fields
  value = inlet["foo"].get<std::string>();
  EXPECT_EQ(value, "test string");

  value = inlet["bar"].get<std::string>();
  EXPECT_EQ(value, "other test string");

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = inlet["nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getNestedStrings)
{
  std::string testString = "foo = { bar = 'yet another string'; baz = 'string 2' }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addString("foo/bar", "bar's description");
  inlet.addString("foo/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addString("foo/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  std::string value = "";

  // Check for existing fields
  value = inlet["foo/bar"].get<std::string>();
  EXPECT_EQ(value, "yet another string");

  value = inlet["foo/baz"].get<std::string>();
  EXPECT_EQ(value, "string 2");

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = inlet["foo/nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TYPED_TEST(inlet_Inlet_basic, getNestedValuesAddedUsingContainer)
{
  std::string strVal = "";
  int intVal = 0;
  double doubleVal = 0;
  bool boolVal = false;

  std::string testString = "foo = { bar = 'yet another string'; so = 3.5; re = 9; mi = true }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Check for existing fields
  Container& container = inlet.addStruct("foo", "A container called foo");
  container.required(true);

  container.addString("bar", "bar's description").required(true);

  container.addDouble("so", "so's description");
  container.addInt("re", "re's description");
  container.addBool("mi", "mi's description");

  //
  // Check stored values from get
  //

  strVal = inlet["foo/bar"].get<std::string>();
  EXPECT_EQ(strVal, "yet another string");

  boolVal = inlet["foo/mi"];
  EXPECT_EQ(boolVal, true);

  doubleVal = inlet["foo/so"];
  EXPECT_EQ(doubleVal, 3.5);

  intVal = inlet["foo/re"];
  EXPECT_EQ(intVal, 9);
}

template <typename InletReader>
class inlet_Inlet_views : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_Inlet_views, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_Inlet_views, NestedContainerViewCheck1)
{
  std::string testString =
    "field1 = true; field2 = 5632; NewContainer = { str = 'hello'; integer = "
    "32 }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);
  inlet.addBool("field1", "this is field #1, a boolean value").required(true);
  inlet.addInt("field2", "this is field #2, an integer").required(false);
  Container& t = inlet.addStruct("NewContainer", "It's blue").required(false);
  t.addString("str", "str's description").required(true);
  t.addInt("integer", "a whole number").required(false);

  axom::sidre::Group* sidreGroup = inlet.sidreGroup();

  EXPECT_EQ(inlet["field1"].type(), InletType::Bool);
  EXPECT_EQ(inlet["field2"].type(), InletType::Integer);
  EXPECT_EQ(inlet["NewContainer/str"].type(), InletType::String);
  EXPECT_EQ(inlet["NewContainer/integer"].type(), InletType::Integer);

  EXPECT_TRUE(sidreGroup->hasView("field1/required"));
  EXPECT_TRUE(sidreGroup->hasView("field2/required"));
  EXPECT_TRUE(sidreGroup->hasView("NewContainer/str/required"));
  EXPECT_TRUE(sidreGroup->hasView("NewContainer/integer/required"));

  EXPECT_TRUE(sidreGroup->hasView("field1/description"));
  EXPECT_TRUE(sidreGroup->hasView("field2/description"));
  EXPECT_TRUE(sidreGroup->hasView("NewContainer/str/description"));
  EXPECT_TRUE(sidreGroup->hasView("NewContainer/integer/description"));
}

TYPED_TEST(inlet_Inlet_views, NestedContainerViewCheck2)
{
  std::string testString =
    "foo = false; bar = true; Container1 = { float1 = 3.14; Container11 = { "
    "Container111 = "
    "{ x = 4 } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);
  inlet.addBool("foo", "foo's description").required(true);
  inlet.addBool("bar", "bar's description").required(false);

  auto& t = inlet.addStruct("Container1", "The first container").required(false);
  t.addDouble("float1", "floating point number within container 1").required(true);
  auto& t2 = t.addStruct("Container11", "Container within Container 1");
  auto& t3 = t2.addStruct("Container111", "Container within Container 11");
  t3.addInt("x", "A variable");

  axom::sidre::Group* sidreGroup = inlet.sidreGroup();

  EXPECT_EQ(inlet["foo"].type(), InletType::Bool);
  EXPECT_EQ(inlet["bar"].type(), InletType::Bool);
  EXPECT_EQ(inlet["Container1/float1"].type(), InletType::Double);

  EXPECT_TRUE(sidreGroup->hasView("foo/required"));
  EXPECT_TRUE(sidreGroup->hasView("bar/required"));
  EXPECT_TRUE(sidreGroup->hasView("Container1/float1/required"));

  EXPECT_TRUE(sidreGroup->hasView("foo/description"));
  EXPECT_TRUE(sidreGroup->hasView("bar/description"));
  EXPECT_TRUE(sidreGroup->hasView("Container1/float1/description"));
  EXPECT_TRUE(sidreGroup->hasView("Container1/Container11/Container111/x/description"));
}

TYPED_TEST(inlet_Inlet_views, NestedContainerViewCheck3)
{
  std::string testString =
    "Container1 = { float1 = 5.6 }; Container2 = { int1 = 95 }; Container3 = { "
    "bool1 = "
    "true }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& t = inlet.addStruct("Container1", "The first container");
  t.addDouble("float1", " A floating point number in Container 1");
  auto& t2 = inlet.addStruct("Container2", "The second container");
  t2.addInt("int1", "An integer in Container 2");
  auto& t3 = inlet.addStruct("Container3", "The third container");
  t3.addBool("bool1", "A boolean value in Container 3");

  axom::sidre::Group* sidreGroup = inlet.sidreGroup();

  EXPECT_EQ(inlet["Container1/float1"].type(), InletType::Double);
  EXPECT_EQ(inlet["Container2/int1"].type(), InletType::Integer);
  EXPECT_EQ(inlet["Container3/bool1"].type(), InletType::Bool);

  EXPECT_TRUE(sidreGroup->hasView("Container1/float1/description"));
  EXPECT_TRUE(sidreGroup->hasView("Container2/int1/description"));
  EXPECT_TRUE(sidreGroup->hasView("Container3/bool1/description"));
}

template <typename InletReader>
class inlet_Inlet_classes : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_Inlet_classes, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_Inlet_classes, mixLevelContainers)
{
  // regression test for the lua stack being at the wrong
  // level with nested containers with non-nested values

  std::string input =
    "thermal_solver={"
    "   u0 = { type = 'function', func = 'BoundaryTemperature'},"
    "   kappa = { type = 'constant', constant = 0.5},"
    "   solver = {"
    "     rel_tol = 1.0e-6,"
    "     abs_tol = 1.0e-12,"
    "     print_level = 0,"
    "     max_iter = 100,"
    "     dt = 1.0,"
    "     steps = 1 "
    "   }"
    "}";

  Inlet inlet = createBasicInlet<TypeParam>(input);

  //
  // Define input file schema
  //

  // fields that are expected to not be present in above string
  inlet.addString("thermal_solver/mesh/filename", "filename");
  inlet.addInt("thermal_solver/mesh/serial", "serial");
  inlet.addInt("thermal_solver/mesh/parallel", "parallel");
  inlet.addInt("thermal_solver/order", "order");
  inlet.addString("thermal_solver/timestepper", "timestepper");

  // Rest of the fields are present
  inlet.addString("thermal_solver/u0/type", "type");
  inlet.addString("thermal_solver/u0/func", "func");

  inlet.addString("thermal_solver/kappa/type", "type");
  inlet.addDouble("thermal_solver/kappa/constant", "constant");

  inlet.addDouble("thermal_solver/solver/rel_tol", "rel_tol");
  inlet.addDouble("thermal_solver/solver/abs_tol", "abs_tol");
  inlet.addInt("thermal_solver/solver/print_level", "print_level");
  inlet.addInt("thermal_solver/solver/max_iter", "max_iter");
  inlet.addDouble("thermal_solver/solver/dt", "dt");
  inlet.addInt("thermal_solver/solver/steps", "steps");

  //
  //  Verify values found in input file
  //
  std::string strVal;
  int intVal;
  double doubleVal;

  // u0
  strVal = inlet["thermal_solver/u0/type"].get<std::string>();
  EXPECT_EQ(strVal, "function");

  strVal = inlet["thermal_solver/u0/func"].get<std::string>();
  EXPECT_EQ(strVal, "BoundaryTemperature");

  // kappa
  strVal = inlet["thermal_solver/kappa/type"].get<std::string>();
  EXPECT_EQ(strVal, "constant");

  doubleVal = inlet["thermal_solver/kappa/constant"];
  EXPECT_EQ(doubleVal, 0.5);

  // solver
  doubleVal = inlet["thermal_solver/solver/rel_tol"];
  EXPECT_EQ(doubleVal, 1.e-6);

  doubleVal = inlet["thermal_solver/solver/abs_tol"];
  EXPECT_EQ(doubleVal, 1.e-12);

  intVal = inlet["thermal_solver/solver/print_level"];
  EXPECT_EQ(intVal, 0);

  intVal = inlet["thermal_solver/solver/max_iter"];
  EXPECT_EQ(intVal, 100);

  doubleVal = inlet["thermal_solver/solver/dt"];
  EXPECT_EQ(doubleVal, 1.0);

  intVal = inlet["thermal_solver/solver/steps"];
  EXPECT_EQ(intVal, 1);
}

TYPED_TEST(inlet_Inlet_classes, defaultValuesDocsEnabled)
{
  std::string testString =
    "Container1 = { float1 = 5.6; Container2 = { int1 = 95; Container4 = { "
    "str1= 'hi' } } "
    "}; Container3 = { bool1 = true }";
  // Creating basic inlet with documentation enabled (indicated by the last param)
  Inlet inlet = createBasicInlet<TypeParam>(testString, true);

  // new fields
  inlet.addDouble("field1").defaultValue(2);  // int argument will get casted to double
  inlet.addInt("Container1/Container2/field2").defaultValue(5);
  inlet.addBool("Container1/field3").defaultValue(true);
  inlet.addString("Container3/field4").defaultValue("default for new string");

  // existing fields
  inlet.addInt("Container1/Container2/int1").defaultValue(100);
  inlet.addBool("Container3/bool1").defaultValue(false);
  inlet.addDouble("Container1/float1").defaultValue(3.14);
  inlet.addString("Container1/Container2/Container4/str1").defaultValue("default for old string");

  axom::sidre::Group* sidreGroup = inlet.sidreGroup();

  EXPECT_TRUE(sidreGroup->hasView("field1/defaultValue"));
  double doubleVal = sidreGroup->getView("field1/defaultValue")->getScalar();
  EXPECT_EQ(doubleVal, 2.0);

  EXPECT_TRUE(sidreGroup->hasView("Container1/Container2/field2/defaultValue"));
  int intVal = sidreGroup->getView("Container1/Container2/field2/defaultValue")->getScalar();
  EXPECT_EQ(intVal, 5);

  EXPECT_TRUE(sidreGroup->hasView("Container1/field3/defaultValue"));
  int8_t boolVal = sidreGroup->getView("Container1/field3/defaultValue")->getScalar();
  EXPECT_EQ(boolVal, 1);

  EXPECT_TRUE(sidreGroup->hasView("Container3/field4/defaultValue"));
  std::string strVal = sidreGroup->getView("Container3/field4/defaultValue")->getString();
  EXPECT_EQ(strVal, "default for new string");

  EXPECT_TRUE(sidreGroup->hasView("Container1/Container2/int1/defaultValue"));
  intVal = sidreGroup->getView("Container1/Container2/int1/defaultValue")->getScalar();
  EXPECT_EQ(intVal, 100);

  EXPECT_TRUE(sidreGroup->hasView("Container3/bool1/defaultValue"));
  boolVal = sidreGroup->getView("Container3/bool1/defaultValue")->getScalar();
  EXPECT_EQ(boolVal, 0);

  EXPECT_TRUE(sidreGroup->hasView("Container1/float1/defaultValue"));
  doubleVal = sidreGroup->getView("Container1/float1/defaultValue")->getScalar();
  EXPECT_EQ(doubleVal, 3.14);

  EXPECT_TRUE(sidreGroup->hasView("Container1/Container2/Container4/str1/defaultValue"));
  strVal = sidreGroup->getView("Container1/Container2/Container4/str1/defaultValue")->getString();
  EXPECT_EQ(strVal, "default for old string");

  doubleVal = sidreGroup->getView("field1/value")->getScalar();
  EXPECT_EQ(doubleVal, 2.0);

  intVal = sidreGroup->getView("Container1/Container2/field2/value")->getScalar();
  EXPECT_EQ(intVal, 5);

  boolVal = sidreGroup->getView("Container1/field3/value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  strVal = sidreGroup->getView("Container3/field4/value")->getString();
  EXPECT_EQ(strVal, "default for new string");

  intVal = sidreGroup->getView("Container1/Container2/int1/value")->getScalar();
  EXPECT_EQ(intVal, 95);

  boolVal = sidreGroup->getView("Container3/bool1/value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  doubleVal = sidreGroup->getView("Container1/float1/value")->getScalar();
  EXPECT_EQ(doubleVal, 5.6);

  strVal = sidreGroup->getView("Container1/Container2/Container4/str1/value")->getString();
  EXPECT_EQ(strVal, "hi");
}

TYPED_TEST(inlet_Inlet_classes, defaultValuesDocsDisabled)
{
  std::string testString =
    "Container1 = { float1 = 5.6; Container2 = { int1 = 95; Container4 = { "
    "str1= 'hi' } } "
    "}; Container3 = { bool1 = true }";
  // Creating basic inlet with documentation disabled (indicated by the last param)
  Inlet inlet = createBasicInlet<TypeParam>(testString, false);

  // new fields
  inlet.addDouble("field1").defaultValue(2.0);
  inlet.addInt("Container1/Container2/field2").defaultValue(5);
  inlet.addBool("Container1/field3").defaultValue(true);
  inlet.addString("Container3/field4").defaultValue("default for new string");

  // existing fields
  inlet.addInt("Container1/Container2/int1").defaultValue(100);
  inlet.addBool("Container3/bool1").defaultValue(false);
  inlet.addDouble("Container1/float1").defaultValue(3.14);
  inlet.addString("Container1/Container2/Container4/str1").defaultValue("default for old string");

  axom::sidre::Group* sidreGroup = inlet.sidreGroup();

  EXPECT_FALSE(sidreGroup->hasView("field1/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Container1/Container2/field2/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Container1/field3/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Container3/field4/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Container1/Container2/int1/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Container3/bool1/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Container1/float1/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Container1/Container2/Container4/str1/defaultValue"));

  double doubleVal = sidreGroup->getView("field1/value")->getScalar();
  EXPECT_EQ(doubleVal, 2.0);

  int intVal = sidreGroup->getView("Container1/Container2/field2/value")->getScalar();
  EXPECT_EQ(intVal, 5);

  int8_t boolVal = sidreGroup->getView("Container1/field3/value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  std::string strVal = sidreGroup->getView("Container3/field4/value")->getString();
  EXPECT_EQ(strVal, "default for new string");

  intVal = sidreGroup->getView("Container1/Container2/int1/value")->getScalar();
  EXPECT_EQ(intVal, 95);

  boolVal = sidreGroup->getView("Container3/bool1/value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  doubleVal = sidreGroup->getView("Container1/float1/value")->getScalar();
  EXPECT_EQ(doubleVal, 5.6);

  strVal = sidreGroup->getView("Container1/Container2/Container4/str1/value")->getString();
  EXPECT_EQ(strVal, "hi");
}

TYPED_TEST(inlet_Inlet_classes, ranges)
{
  std::string testString =
    "Container1 = { float1 = 5.6; Container2 = { int1 = 95; Container4 = { "
    "str1= 'hi' } } "
    "}; Container3 = { bool1 = true }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  axom::sidre::Group* sidreGroup = inlet.sidreGroup();

  inlet.addInt("Container1/set").validValues({2, 4, 6});
  EXPECT_TRUE(sidreGroup->hasView("Container1/set/validValues"));
  int* bufferArr1 = sidreGroup->getView("Container1/set/validValues")->getArray();
  EXPECT_TRUE(bufferArr1);
  EXPECT_EQ(sidreGroup->getView("Container1/set/validValues")->getBuffer()->getNumElements(), 3);
  EXPECT_EQ(bufferArr1[0], 2);
  EXPECT_EQ(bufferArr1[1], 4);
  EXPECT_EQ(bufferArr1[2], 6);

  inlet.addDouble("Container1/Container2/Container4/double1").range(2.0, 5.0);
  EXPECT_TRUE(sidreGroup->hasView("Container1/Container2/Container4/double1/range"));
  double* bufferArr2 =
    sidreGroup->getView("Container1/Container2/Container4/double1/range")->getArray();
  EXPECT_TRUE(bufferArr2);
  EXPECT_EQ(
    sidreGroup->getView("Container1/Container2/Container4/double1/range")->getBuffer()->getNumElements(),
    2);
  EXPECT_EQ(bufferArr2[0], 2.0);
  EXPECT_EQ(bufferArr2[1], 5.0);

  inlet.addInt("Container1/Container2/int1").range(1, 50);
  EXPECT_TRUE(sidreGroup->hasView("Container1/Container2/int1/range"));
  int* bufferArr3 = sidreGroup->getView("Container1/Container2/int1/range")->getArray();
  EXPECT_TRUE(bufferArr3);
  EXPECT_EQ(sidreGroup->getView("Container1/Container2/int1/range")->getBuffer()->getNumElements(),
            2);
  EXPECT_EQ(bufferArr3[0], 1);
  EXPECT_EQ(bufferArr3[1], 50);
}

template <typename InletReader>
class inlet_Inlet_verify : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_Inlet_verify, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_Inlet_verify, verifyRequired)
{
  std::string testString =
    "field1 = true; field2 = 5632; NewContainer = { str = 'hello'; integer = "
    "32 }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addString("NewContainer/str").required(true);
  inlet.addInt("NewContainer/int").required(false);
  inlet.addBool("field1").required(true);
  inlet.addInt("field2").required(true);
  EXPECT_TRUE(inlet.verify());

  inlet.addBool("NewContainer/field3").required(true);
  inlet.addStruct("NewContainer/LastContainer").addDouble("field4").required(true);
  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_Inlet_verify, verifyDoubleRange)
{
  // For checking values
  std::string testString =
    "field1 = true; field2 = 56.32; NewContainer = { str = 'hello'; field4 = "
    "22.19 "
    "}";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addDouble("field2").range(1.0, 57.2);
  EXPECT_TRUE(inlet.verify());

  inlet.addDouble("field3").range(-10.5, 13.23);
  EXPECT_TRUE(inlet.verify());

  inlet.addDouble("NewContainer/field4").range(1.0, 22.0);
  EXPECT_FALSE(inlet.verify());

  // For checking default values
  std::string testString1 =
    "field1 = true; field2 = 56.32; NewContainer = { str = 'hello'; field4 = "
    "22.19 "
    "}";
  Inlet inlet1 = createBasicInlet<TypeParam>(testString1);

  inlet1.addDouble("field2").defaultValue(5).range(0,
                                                   60);  // ints will get casted to double
  EXPECT_TRUE(inlet1.verify());

  inlet1.addDouble("field3").defaultValue(-12).range(-10.5, 13.23);
  EXPECT_FALSE(inlet1.verify());
}

TYPED_TEST(inlet_Inlet_verify, verifyIntRange)
{
  // For checking values
  std::string testString =
    "field1 = true; field2 = 56; NewContainer = { field4 = 22; field5 = 48 }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("field2").range(0, 56);
  EXPECT_TRUE(inlet.verify());

  inlet.addInt("field3").range(-12, 13);
  EXPECT_TRUE(inlet.verify());

  inlet.addInt("NewContainer/field4").range(1, 23);
  EXPECT_TRUE(inlet.verify());

  inlet.addInt("NewContainer/field5").range(1, 7);
  EXPECT_FALSE(inlet.verify());

  // For checking default values
  std::string testString1 =
    "field1 = true; field2 = 56; NewContainer = { field4 = 22; field5 = 48 }";
  Inlet inlet1 = createBasicInlet<TypeParam>(testString1);
  inlet1.addInt("field2").range(0, 56).defaultValue(32);
  EXPECT_TRUE(inlet1.verify());

  inlet1.addInt("field3").range(-12, 13).defaultValue(-3);
  EXPECT_TRUE(inlet1.verify());

  inlet1.addInt("NewContainer/field4").range(1, 23).defaultValue(24);
  EXPECT_FALSE(inlet1.verify());
}

TYPED_TEST(inlet_Inlet_verify, verifyValidIntValues)
{
  // check values
  std::string testString =
    "field1 = true; field2 = 56; NewContainer = { field4 = 22; field5 = 48 }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("field2").validValues({1, 2, 3, 56, 57, 58});
  EXPECT_TRUE(inlet.verify());

  std::vector<int> nums = {-1, -2, -6, -18, 21};

  inlet.addInt("field3").validValues(nums);
  EXPECT_TRUE(inlet.verify());

  inlet.addInt("NewContainer/field4").validValues({44, 40, 48, 22});
  EXPECT_TRUE(inlet.verify());

  inlet.addInt("NewContainer/field5").validValues(nums);
  EXPECT_FALSE(inlet.verify());

  // check default values
  std::string testString1 =
    "field1 = true; field2 = 56; NewContainer = { field4 = 22; field5 = 48 }";
  Inlet inlet1 = createBasicInlet<TypeParam>(testString1);

  inlet1.addInt("field2").validValues({1, 2, 3, 56, 57, 58}).defaultValue(2);
  EXPECT_TRUE(inlet1.verify());

  inlet1.addInt("field3").validValues(nums).defaultValue(21);
  EXPECT_TRUE(inlet1.verify());

  inlet1.addInt("NewContainer/field4").validValues({44, 40, 48, 22}).defaultValue(48);
  EXPECT_TRUE(inlet1.verify());

  inlet1.addInt("NewContainer/field5").validValues(nums).defaultValue(90);
  EXPECT_FALSE(inlet1.verify());
}

TYPED_TEST(inlet_Inlet_verify, verifyValidDoubleValues)
{
  // check values
  std::string testString =
    "field1 = true; field2 = 56.0; NewContainer = { field4 = 22.0; field5 = "
    "48.23 "
    "}";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addDouble("field2").validValues({1, 2, 3, 56, 57, 58});
  EXPECT_TRUE(inlet.verify());

  std::vector<int> nums = {-1, -2, -6, -18, 21};

  inlet.addDouble("field3").validValues(nums);
  EXPECT_TRUE(inlet.verify());

  inlet.addDouble("NewContainer/field4").validValues({44.1, 40., 48., 22.});
  EXPECT_TRUE(inlet.verify());

  inlet.addDouble("NewContainer/field5").validValues(nums);
  EXPECT_FALSE(inlet.verify());

  // check default values
  std::string testString1 =
    "field1 = true; field2 = 56.0; NewContainer = { field4 = 22.0; field5 = "
    "48.23 "
    "}";
  Inlet inlet1 = createBasicInlet<TypeParam>(testString1);

  inlet1.addDouble("field2").validValues({1, 2, 3, 56, 57, 58}).defaultValue(2.);
  EXPECT_TRUE(inlet1.verify());

  inlet1.addDouble("field3").validValues(nums).defaultValue(21);
  EXPECT_TRUE(inlet1.verify());

  inlet1.addDouble("NewContainer/field4").validValues({44.05, 40.13, 48.28, 22.}).defaultValue(48.28);
  EXPECT_TRUE(inlet1.verify());

  inlet1.addDouble("NewContainer/field5").validValues(nums).defaultValue(90.1);
  EXPECT_FALSE(inlet1.verify());
}

TYPED_TEST(inlet_Inlet_verify, verifyValidStringValues)
{
  // check values
  std::string testString =
    "field1 = true; field2 = 'abc'; NewContainer = { field3 = 'xyz'; field4 = "
    "'yes' }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addString("field2").validValues({"abc", "defg", "hijk", "lm"});
  EXPECT_TRUE(inlet.verify());

  std::vector<std::string> strs = {"nop", "qrstuv", "xyz", "wx"};

  inlet.addString("NewContainer/field3").validValues(strs);
  EXPECT_TRUE(inlet.verify());

  inlet.addString("Container1/field5").validValues(strs);
  EXPECT_TRUE(inlet.verify());

  inlet.addString("NewContainer/field4").validValues(strs);
  EXPECT_FALSE(inlet.verify());

  // check default values
  std::string testString1 = "field1 = true; NewContainer = { field5 = 'nop' }";
  Inlet inlet1 = createBasicInlet<TypeParam>(testString1);

  inlet1.addString("field2").validValues({"abc", "defg", "hijk", "lm"}).defaultValue("defg");
  EXPECT_TRUE(inlet1.verify());

  inlet1.addString("field3").validValues(strs).defaultValue("wx");
  EXPECT_TRUE(inlet1.verify());

  inlet1.addString("NewContainer/field4").validValues(strs);
  EXPECT_TRUE(inlet1.verify());

  inlet1.addString("NewContainer/field5").validValues(strs).defaultValue("zyx");
  EXPECT_FALSE(inlet1.verify());
}

TYPED_TEST(inlet_Inlet_verify, verifyFieldLambda)
{
  std::string testString =
    "field1 = true; field2 = 'abc'; NewContainer = { field3 = 'xyz'; field4 = "
    "'yes' }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addBool("field1");
  auto& field2 = inlet.addString("field2");
  auto& field3 = inlet.addString("NewContainer/field3");

  field2.registerVerifier([](const Field& field) {
    auto str = field.get<std::string>();
    return (str.size() >= 1 && str[0] == 'a');
  });
  EXPECT_TRUE(inlet.verify());

  field3.registerVerifier([](const Field& field) {
    auto str = field.get<std::string>();
    return (str.size() >= 1 && str[0] == 'a');
  });
  EXPECT_FALSE(inlet.verify());

  // Check that verifiers that take a list of errors are given a proper list
  // they can modify and that the errors are returned to the caller.
  field2.registerVerifier([](const Field&, std::vector<VerificationError>* errors) {
    INLET_VERIFICATION_WARNING("<base>", "bad thing", errors);
    return false;
  });
  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  EXPECT_THAT(
    errors,
    Contains(Truly([](const VerificationError& error) { return error.message == "bad thing"; })));
}

TYPED_TEST(inlet_Inlet_verify, verifyContainerLambda1)
{
  std::string testString =
    "field1 = true; field2 = 'abc'; NewContainer = { field3 = 'xyz'; field4 = "
    "'yes' }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addBool("field1");
  auto& field2 = inlet.addString("field2");
  Container& container1 = inlet.addStruct("NewContainer");
  auto& field3 = container1.addString("field3");

  field2.registerVerifier([](const Field& field) {
    auto str = field.get<std::string>();
    return (str.size() >= 1 && str[0] == 'a');
  });
  EXPECT_TRUE(inlet.verify());

  field3.registerVerifier([](const Field& field) {
    auto str = field.get<std::string>();
    return (str.size() >= 1 && str[0] == 'x');
  });
  EXPECT_TRUE(inlet.verify());

  EXPECT_TRUE(container1.contains("field3"));

  container1.registerVerifier([](const Container& container) { return container.contains("field3"); });
  EXPECT_TRUE(inlet.verify());

  container1.registerVerifier(
    [](const Container& container) { return container.contains("field22"); });
  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_Inlet_verify, verifyContainerLambda3)
{
  std::string testString = "dimensions = 2; vector = { x = 1; y = 2; z = 3; }";
  auto myInlet = createBasicInlet<TypeParam>(testString);
  myInlet.addInt("dimensions").required(true);
  auto& v = myInlet.addStruct("vector").required(true);
  v.addInt("x");

  v.registerVerifier([&myInlet](const Container& container) {
    int dim = myInlet["dimensions"];
    bool x_present = container.contains("x") && (container["x"].type() == InletType::Integer);
    bool y_present = container.contains("y") && (container["y"].type() == InletType::Integer);
    bool z_present = container.contains("z") && (container["z"].type() == InletType::Integer);
    if(dim == 1 && x_present)
    {
      return true;
    }
    else if(dim == 2 && x_present && y_present)
    {
      return true;
    }
    else if(dim == 3 && x_present && y_present && z_present)
    {
      return true;
    }
    return false;
  });

  EXPECT_FALSE(myInlet.verify());

  v.addInt("y");
  v.addInt("z");

  EXPECT_TRUE(myInlet.verify());
}

template <typename InletReader>
class inlet_Inlet_array : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_Inlet_array, axom::inlet::detail::ReaderTypes);

// Checks all of the Container::getArray functions
TYPED_TEST(inlet_Inlet_array, getArray)
{
  std::string testString =
    "luaArrays = { arr1 = { [0] = 4}, "
    "              arr2 = {[0] = true, [1] = false}, "
    "              arr3 = {[0] = 'hello', [1] = 'bye'}, "
    "              arr4 = { [0] = 2.4 } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addIntArray("luaArrays/arr1");
  inlet.addBoolArray("luaArrays/arr2");
  inlet.addStringArray("luaArrays/arr3");
  inlet.addDoubleArray("luaArrays/arr4");

  std::unordered_map<int, int> expectedInts {{0, 4}};
  std::unordered_map<int, bool> expectedBools {{0, true}, {1, false}};
  std::unordered_map<int, double> expectedDoubles {{0, 2.4}};
  std::unordered_map<int, std::string> expectedStrs {{0, "hello"}, {1, "bye"}};

  std::unordered_map<int, int> intMap = inlet["luaArrays/arr1"];
  std::unordered_map<int, bool> boolMap = inlet["luaArrays/arr2"];
  std::unordered_map<int, std::string> strMap = inlet["luaArrays/arr3"];
  std::unordered_map<int, double> doubleMap = inlet["luaArrays/arr4"];

  EXPECT_EQ(intMap, expectedInts);

  EXPECT_EQ(boolMap, expectedBools);

  EXPECT_EQ(strMap, expectedStrs);

  EXPECT_EQ(doubleMap, expectedDoubles);
}

// Checks the underlying Sidre representation of the arrays added from Lua
TYPED_TEST(inlet_Inlet_array, inletArraysInSidre)
{
  std::string testString =
    "luaArrays = { arr1 = { [0] = 4.0, [1] = 5.0, [2] = 6.0, [3] = 2.4}, "
    "              arr2 = { [0] = true, [1] = false}, "
    "              arr3 = { [0] = 'hello', [1] = 'bye'}, "
    "              arr4 = { [0] = 2.4 } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addIntArray("luaArrays/arr1");

  const axom::sidre::Group* group =
    inlet.sidreGroup()->getGroup("luaArrays/arr1/_inlet_collection");
  auto idx = group->getGroup("0");
  EXPECT_TRUE(idx);
  int val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 4);

  idx = group->getGroup("3");
  EXPECT_TRUE(idx);
  val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 2);

  idx = group->getGroup("1");
  EXPECT_TRUE(idx);
  val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 5);

  idx = group->getGroup("2");
  EXPECT_TRUE(idx);
  val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 6);

  inlet.addBoolArray("luaArrays/arr2");
  group = inlet["luaArrays/arr2/_inlet_collection"].sidreGroup();

  idx = group->getGroup("0");
  EXPECT_TRUE(idx);
  int8_t boolVal = idx->getView("value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  idx = group->getGroup("1");
  EXPECT_TRUE(idx);
  boolVal = idx->getView("value")->getScalar();
  EXPECT_EQ(boolVal, 0);

  inlet.addStringArray("luaArrays/arr3");
  group = inlet["luaArrays/arr3/_inlet_collection"].sidreGroup();

  idx = group->getGroup("0");
  EXPECT_TRUE(idx);
  std::string str = idx->getView("value")->getString();
  EXPECT_EQ(str, "hello");

  idx = group->getGroup("1");
  EXPECT_TRUE(idx);
  str = idx->getView("value")->getString();
  EXPECT_EQ(str, "bye");

  inlet.addDoubleArray("luaArrays/arr4");
  group = inlet["luaArrays/arr4/_inlet_collection"].sidreGroup();

  idx = group->getGroup("0");
  EXPECT_TRUE(idx);
  double doubleVal = idx->getView("value")->getScalar();
  EXPECT_EQ(doubleVal, 2.4);
}

#ifdef AXOM_USE_SOL
// Using integer literals in strings is lua-specific
TEST(inlet_Inlet_basic_lua, getTopLevelStrings)
{
  std::string testString = "foo = 'test string'; bar = '15'";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addString("foo", "foo's description");
  inlet.addString("bar", "bar's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addString("nonexistant", "nothing");

  //
  // Check stored values from get
  //

  std::string value = "";

  // Check for existing fields
  value = inlet["foo"].get<std::string>();
  EXPECT_EQ(value, "test string");

  value = inlet["bar"].get<std::string>();
  EXPECT_EQ(value, "15");

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = inlet["nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

// Empty string literals are lua-specific
TEST(inlet_Inlet_basic_lua, getNestedStrings)
{
  std::string testString = "foo = { bar = 'yet another string'; baz = '' }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  //
  // Define schema
  //

  // Check for existing fields
  inlet.addString("foo/bar", "bar's description");
  inlet.addString("foo/baz", "baz's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addString("foo/nonexistant", "nothing");

  //
  // Check stored values from get
  //

  std::string value = "";

  // Check for existing fields
  value = inlet["foo/bar"].get<std::string>();
  EXPECT_EQ(value, "yet another string");

  value = inlet["foo/baz"].get<std::string>();
  EXPECT_EQ(value, "");

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = inlet["foo/nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TEST(inlet_Inlet_verify_lua, verifyContainerLambda2)
{
  std::string testString =
    "thermal_solver={}\n"
    "thermal_solver.order = 2\n"
    "thermal_solver.timestepper = 'quasistatic'\n"
    "solid_solver={}\n"
    "solid_solver.order = 3\n"
    "solid_solver.timestepper = 'BackwardEuler'\n"
    "material={}\n"
    "material.attribute = 1\n"
    "material.thermalview = 'isotropic'\n"
    "material.solidview = 'hyperelastic'";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);
  inlet.addInt("thermal_solver/order");
  inlet.addString("thermal_solver/timestepper");
  inlet.addInt("solid_solver/order");
  inlet.addString("solid_solver/timestepper");
  inlet.addInt("material/attribute");
  auto& globalContainer = inlet.getGlobalContainer();
  auto material = globalContainer["material"];

  globalContainer.registerVerifier([](const Container& container) {
    bool verifySuccess = true;
    if(container.contains("thermal_solver") && !container["material"].contains("thermalview"))
    {
      verifySuccess = false;
    }
    if(container.contains("solid_solver") && !container["material"].contains("solidview"))
    {
      verifySuccess = false;
    }
    return verifySuccess;
  });

  EXPECT_FALSE(inlet.verify());

  inlet.addString("material/thermalview");
  inlet.addString("material/solidview");

  EXPECT_TRUE(material.contains("solidview"));
  EXPECT_TRUE(material.contains("thermalview"));

  EXPECT_TRUE(inlet.verify());

  auto& thing = inlet.addStruct("test");
  auto& thing2 = thing.addStruct("test2");
  thing2.addStruct("test3/test4/test5");

  EXPECT_EQ(globalContainer["test"].name(), "test");
  EXPECT_EQ(thing["test2"].name(), "test/test2");
  EXPECT_EQ(thing2["test3/test4/test5"].name(), "test/test2/test3/test4/test5");
}

TEST(inlet_Inlet_verify_lua, requiredContainer)
{
  std::string testString =
    "thermal_solver={}\n"
    "thermal_solver.order = 2\n"
    "thermal_solver.timestepper = 'quasistatic'\n"
    "solid_solver={}\n"
    "solid_solver.order = 3\n"
    "solid_solver.timestepper = 'BackwardEuler'\n"
    "material={}\n"
    "material.attribute = 1\n"
    "material.thermalview = 'isotropic'\n"
    "material.solidview = 'hyperelastic'";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);
  auto& material = inlet.addStruct("material");
  material.required(true);

  EXPECT_FALSE(inlet.verify());

  inlet.addString("material/thermalview");
  inlet.addString("material/solidview");

  EXPECT_TRUE(inlet.contains("material/thermalview"));
  EXPECT_TRUE(inlet.contains("material/solidview"));

  EXPECT_TRUE(inlet.verify());
}

// Check discontiguous arrays specifically
TEST(inlet_Inlet_array_lua, getArray)
{
  std::string testString =
    "luaArrays = { arr1 = { [1] = 4}, "
    "              arr2 = {[4] = true, [8] = false}, "
    "              arr3 = {[33] = 'hello', [2] = 'bye'}, "
    "              arr4 = { [12] = 2.4 } }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  inlet.addIntArray("luaArrays/arr1");
  inlet.addBoolArray("luaArrays/arr2");
  inlet.addStringArray("luaArrays/arr3");
  inlet.addDoubleArray("luaArrays/arr4");

  std::unordered_map<int, int> expectedInts {{1, 4}};
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  std::unordered_map<int, double> expectedDoubles {{12, 2.4}};
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"}, {2, "bye"}};

  std::unordered_map<int, int> intMap = inlet["luaArrays/arr1"];
  std::unordered_map<int, bool> boolMap = inlet["luaArrays/arr2"];
  std::unordered_map<int, std::string> strMap = inlet["luaArrays/arr3"];
  std::unordered_map<int, double> doubleMap = inlet["luaArrays/arr4"];

  EXPECT_EQ(intMap, expectedInts);

  EXPECT_EQ(boolMap, expectedBools);

  EXPECT_EQ(strMap, expectedStrs);

  EXPECT_EQ(doubleMap, expectedDoubles);
}

// Checks the underlying Sidre representation of the arrays added from Lua
TEST(inlet_Inlet_array_lua, inletArraysInSidre)
{
  std::string testString =
    "luaArrays = { arr1 = { [1] = 4, [2] = 5, [3] = 6 , [12] = 2.4}, "
    "              arr2 = { [4] = true, [8] = false}, "
    "              arr3 = { [33] = 'hello', [2] = 'bye'}, "
    "              arr4 = { [12] = 2.4 } }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  inlet.addIntArray("luaArrays/arr1");

  const axom::sidre::Group* group =
    inlet.sidreGroup()->getGroup("luaArrays/arr1/_inlet_collection");
  auto idx = group->getGroup("1");
  EXPECT_TRUE(idx);
  int val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 4);

  idx = group->getGroup("12");
  EXPECT_TRUE(idx);
  val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 2);

  idx = group->getGroup("2");
  EXPECT_TRUE(idx);
  val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 5);

  idx = group->getGroup("3");
  EXPECT_TRUE(idx);
  val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 6);

  inlet.addBoolArray("luaArrays/arr2");
  group = inlet["luaArrays/arr2/_inlet_collection"].sidreGroup();

  idx = group->getGroup("4");
  EXPECT_TRUE(idx);
  int8_t boolVal = idx->getView("value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  idx = group->getGroup("8");
  EXPECT_TRUE(idx);
  boolVal = idx->getView("value")->getScalar();
  EXPECT_EQ(boolVal, 0);

  inlet.addStringArray("luaArrays/arr3");
  group = inlet["luaArrays/arr3/_inlet_collection"].sidreGroup();

  idx = group->getGroup("33");
  EXPECT_TRUE(idx);
  std::string str = idx->getView("value")->getString();
  EXPECT_EQ(str, "hello");

  idx = group->getGroup("2");
  EXPECT_TRUE(idx);
  str = idx->getView("value")->getString();
  EXPECT_EQ(str, "bye");

  inlet.addDoubleArray("luaArrays/arr4");
  group = inlet["luaArrays/arr4/_inlet_collection"].sidreGroup();

  idx = group->getGroup("12");
  EXPECT_TRUE(idx);
  double doubleVal = idx->getView("value")->getScalar();
  EXPECT_EQ(doubleVal, 2.4);
}

#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
