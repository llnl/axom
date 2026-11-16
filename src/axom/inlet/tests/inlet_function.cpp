// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"

#include "axom/core/utilities/FileUtilities.hpp"

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/YAMLReader.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/JSONSchemaWriter.hpp"
#include "axom/inlet/SphinxWriter.hpp"

#include "gtest/gtest.h"

#include <array>
#include <functional>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <iostream>

using axom::inlet::FunctionTag;
using axom::inlet::FunctionType;
using axom::inlet::Inlet;
using axom::inlet::InletType;
using axom::inlet::LuaReader;
using axom::inlet::VerificationError;

#include "axom/sol.hpp"

class SolStateReader : public LuaReader
{
public:
  using LuaReader::solState;
};

Inlet createBasicInlet(const std::string& luaString, bool enableDocs = true)
{
  auto lr = std::make_unique<LuaReader>();
  lr->parseString(luaString);
  return Inlet(std::move(lr), enableDocs);
}

/// Collects SLIC error messages; registered with the logger in main()
std::ostringstream& slicErrors()
{
  static std::ostringstream errors;
  return errors;
}

/// Checks whether a SLIC error containing \a expected was logged, then clears the log
::testing::AssertionResult slicErrorContains(const std::string& expected)
{
  const std::string errors = slicErrors().str();
  slicErrors().str("");
  if(errors.find(expected) != std::string::npos)
  {
    return ::testing::AssertionSuccess();
  }
  return ::testing::AssertionFailure()
    << "Expected a SLIC error containing '" << expected << "', got: '" << errors << "'";
}

TEST(inlet_function, simple_vec3_to_double_raw)
{
  std::string testString = "function foo (v) return v.x + v.y + v.z end";
  auto inlet = createBasicInlet(testString);

  auto func = inlet.reader().getFunction("foo", FunctionTag::Double, {FunctionTag::Vector});

  EXPECT_TRUE(func);
  auto result = func.call<double>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_raw)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {FunctionTag::Vector});

  EXPECT_TRUE(func);
  auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, vector_function_accepts_lua_table_returns)
{
  auto inlet = createBasicInlet(R"(
    function make_vector (dim)
      local result = {}
      for i = 1, dim do result[i] = i end
      return result
    end
  )");

  auto function =
    inlet.reader().getFunction("make_vector", FunctionTag::Vector, {FunctionTag::Double});
  ASSERT_TRUE(function);
  for(int dim = 1; dim <= 3; ++dim)
  {
    const auto result = function.call<FunctionType::Vector>(static_cast<double>(dim));
    ASSERT_EQ(result.dim, dim);
    for(int component = 0; component < result.dim; ++component)
    {
      EXPECT_DOUBLE_EQ(result[component], component + 1.0);
    }
  }
}

TEST(inlet_function, vector_function_accepts_equivalent_table_forms)
{
  // Explicit keys, in any order, and float keys with integral values are the same vector
  for(const std::string table : {"{4.5, -2}", "{[2] = -2, [1] = 4.5}", "{[1.0] = 4.5, [2.0] = -2}"})
  {
    auto inlet = createBasicInlet("function foo () return " + table + " end");
    inlet.addFunction("foo", FunctionTag::Vector, {});
    ASSERT_TRUE(inlet.verify()) << table;

    // Through Proxy::call and through a std::function copy
    const auto viaCall = inlet["foo"].call<FunctionType::Vector>();
    const auto viaCopy = inlet["foo"].get<std::function<FunctionType::Vector()>>()();
    for(const auto& result : {viaCall, viaCopy})
    {
      ASSERT_EQ(result.dim, 2) << table;
      EXPECT_DOUBLE_EQ(result[0], 4.5) << table;
      EXPECT_DOUBLE_EQ(result[1], -2.0) << table;
    }
  }
}

TEST(inlet_function, vector_function_rejects_malformed_table_returns)
{
  axom::slic::ScopedAbortToThrow abortGuard;

  // Lua vectors must be dense numeric sequences with a supported dimension.
  // Each case fails exactly one check; match its message to confirm which one.
  const std::array<std::pair<std::string, std::string>, 13> cases {{
    {"2.0", "return types possibly incorrect"},
    {"'text'", "return types possibly incorrect"},
    {"{}", "empty table"},
    {"{1, 2, 3, 4}", "integers between 1 and 3"},
    {"{[1] = 1, [3] = 3}", "contiguous starting at 1"},
    {"{1, 'two'}", "components must be numeric"},
    {"{1, true}", "components must be numeric"},
    {"{1, {2}}", "components must be numeric"},
    {"{1, 2, label = 3}", "only contain numeric indices"},
    {"{[0] = 0, 1, 2}", "integers between 1 and 3"},
    {"{[1.5] = 1}", "integers between 1 and 3"},
    {"{[-1] = 1}", "integers between 1 and 3"},
    {"{[math.huge] = 1}", "integers between 1 and 3"},
  }};

  for(const auto& [result, message] : cases)
  {
    auto inlet = createBasicInlet("function foo () return " + result + " end");
    auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {});
    ASSERT_TRUE(func);
    EXPECT_THROW(func.call<FunctionType::Vector>(), axom::slic::SlicAbortException) << result;
    EXPECT_TRUE(slicErrorContains(message)) << result;
  }
}

TEST(inlet_function, lua_callback_failures_report_slic_errors)
{
  axom::slic::ScopedAbortToThrow abortGuard;

  auto inlet = createBasicInlet(R"(
    function runtime_error () error('callback failed') end
    function wrong_type () return 'not a number' end
  )");

  auto wrongType = inlet.reader().getFunction("wrong_type", FunctionTag::Double, {});
  ASSERT_TRUE(wrongType);
  EXPECT_THROW(wrongType.call<FunctionType::Double>(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("return types possibly incorrect"));

  // The SLIC error includes Lua's error message
  auto runtimeError = inlet.reader().getFunction("runtime_error", FunctionTag::Double, {});
  ASSERT_TRUE(runtimeError);
  EXPECT_THROW(runtimeError.call<FunctionType::Double>(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("callback failed"));
}

TEST(inlet_function, lua_callback_failures_report_slic_errors_through_every_access_path)
{
  axom::slic::ScopedAbortToThrow abortGuard;

  auto inlet = createBasicInlet(R"(
    function fail_void (x) error('void failed') end
    function add_one (x) return x + 1 end
  )");
  inlet.addFunction("fail_void", FunctionTag::Void, {FunctionTag::Double});
  inlet.addFunction("add_one", FunctionTag::Double, {FunctionTag::String});
  ASSERT_TRUE(inlet.verify());

  // A void callback still reports an execution error
  EXPECT_THROW(inlet["fail_void"].call<void>(1.0), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("void failed"));

  // An invalid input reports an execution error through Proxy::call.
  EXPECT_THROW(inlet["add_one"].call<double>(std::string {"text"}), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("Lua function call failed:"));

  // Copies made with get<std::function> report errors the same way
  auto addOne = inlet["add_one"].get<std::function<double(std::string)>>();
  EXPECT_DOUBLE_EQ(addOne("2"), 3.0);
  EXPECT_THROW(addOne("text"), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("Lua function call failed:"));
}

TEST(inlet_function, lua_callback_failure_in_verifier_reports_slic_error)
{
  axom::slic::ScopedAbortToThrow abortGuard;

  auto inlet = createBasicInlet("function scale () error('scale failed') end");
  inlet.addFunction("scale", FunctionTag::Double, {})
    .registerVerifier([](const axom::inlet::Function& func) { return func.call<double>() > 0.0; });

  EXPECT_THROW(inlet.verify(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("scale failed"));
}

TEST(inlet_function, simple_vec3_to_vec3_raw_partial_init)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {FunctionTag::Vector});

  EXPECT_TRUE(func);

  auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 2});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 0);

  result = func.call<FunctionType::Vector>(FunctionType::Vector {1});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 0);
  EXPECT_FLOAT_EQ(result[2], 0);
}

TEST(inlet_function, simple_vec3_to_double_through_container)
{
  std::string testString = "function foo (v) return v.x + v.y + v.z end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::Vector}, "foo's description");

  auto callable = inlet["foo"].get<std::function<double(FunctionType::Vector)>>();
  auto result = callable({1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_through_container)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Vector(FunctionType::Vector)>>();
  auto result = callable({1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_double_to_double_through_container)
{
  std::string testString = "function foo (a) return (a * 3.4) + 9.64 end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::Double}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Double(FunctionType::Double)>>();
  double arg = -6.37;
  double result = callable(arg);
  EXPECT_FLOAT_EQ(result, (arg * 3.4) + 9.64);
}

TEST(inlet_function, returned_function_keeps_lua_state_alive)
{
  // An extracted callback must retain its Lua state after Inlet is destroyed.
  std::function<double(double)> callback;
  {
    auto inlet = createBasicInlet("offset = 3.0; function foo (value) return value + offset end");
    inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::Double});
    callback = inlet["foo"].get<std::function<double(double)>>();
  }

  EXPECT_DOUBLE_EQ(callback(4.0), 7.0);
}

TEST(inlet_function, returned_functions_share_lua_state_after_inlet_is_destroyed)
{
  // Callbacks from one input share its Lua state, including after the Inlet is gone.
  // Copying, calling, and destroying them in any order must leave the state usable.
  std::function<double()> increment;
  std::function<double()> current;
  {
    auto inlet = createBasicInlet(R"(
      count = 0
      function increment () count = count + 1; return count end
      function current () return count end
    )");
    inlet.addFunction("increment", FunctionTag::Double, {});
    inlet.addFunction("current", FunctionTag::Double, {});
    increment = inlet["increment"].get<std::function<double()>>();
    current = inlet["current"].get<std::function<double()>>();
  }

  EXPECT_DOUBLE_EQ(increment(), 1.0);
  auto incrementCopy = increment;
  increment = nullptr;
  EXPECT_DOUBLE_EQ(incrementCopy(), 2.0);
  incrementCopy = nullptr;
  EXPECT_DOUBLE_EQ(current(), 2.0);
}

TEST(inlet_function, reader_function_outlives_reader)
{
  axom::inlet::FunctionVariant function;
  {
    LuaReader reader;
    reader.parseString("function foo (value) return 2 * value end");
    function = reader.getFunction("foo", FunctionTag::Double, {FunctionTag::Double});
  }

  ASSERT_TRUE(function);
  EXPECT_DOUBLE_EQ(function.call<double>(4.0), 8.0);
}

TEST(inlet_function, function_declaration_is_idempotent_across_containers)
{
  // Redeclaring a function returns the existing one, whichever Container declared it first
  for(const bool firstOnRoot : {true, false})
  {
    for(const bool secondOnRoot : {true, false})
    {
      auto inlet = createBasicInlet("group = { f = function() return 2.0 end }");
      auto& group = inlet.addStruct("group");

      auto declare = [&](bool onRoot) -> axom::inlet::Verifiable<axom::inlet::Function>& {
        return onRoot ? inlet.addFunction("group/f", FunctionTag::Double, {})
                      : group.addFunction("f", FunctionTag::Double, {});
      };

      auto& first = declare(firstOnRoot);
      auto& second = declare(secondOnRoot);
      EXPECT_EQ(&first, &second) << firstOnRoot << secondOnRoot;
      EXPECT_EQ(1u, group.getChildFunctions().size());
      EXPECT_TRUE(inlet.getGlobalContainer().getChildFunctions().empty());
      EXPECT_TRUE(inlet.verify());
      EXPECT_DOUBLE_EQ(group["f"].call<double>(), 2.0);
    }
  }
}

TEST(inlet_function, or_function_selects_the_supplied_form)
{
  auto inlet = createBasicInlet(R"(
    label = function() return 'computed' end
    scale = 4.0
  )");
  inlet.addString("label").orFunction();
  inlet.addDouble("scale").orFunction();
  EXPECT_TRUE(inlet.verify());

  // type() tells the caller which form the input supplied
  EXPECT_EQ(inlet["label"].type(), InletType::Function);
  EXPECT_EQ(inlet["scale"].type(), InletType::Double);

  // evaluate() handles either form
  EXPECT_EQ(inlet["label"].evaluate<std::string>(), "computed");
  EXPECT_DOUBLE_EQ(inlet["scale"].evaluate<double>(), 4.0);

  // The concrete form is still an ordinary value
  EXPECT_DOUBLE_EQ(inlet["scale"].get<double>(), 4.0);
}

TEST(inlet_function, or_function_evaluate_passes_arguments_to_the_function_only)
{
  auto inlet = createBasicInlet(R"(
    doubled = function(x) return 2 * x end
    constant = 4.0
    greeting = function(name) return 'hello ' .. name end
  )");
  inlet.addDouble("doubled").orFunction({FunctionTag::Double});
  inlet.addDouble("constant").orFunction({FunctionTag::Double});
  inlet.addString("greeting").orFunction({FunctionTag::String});
  EXPECT_TRUE(inlet.verify());

  // Numbers are passed as double, so an int literal or an lvalue also matches
  const double x = 3.0;
  EXPECT_DOUBLE_EQ(inlet["doubled"].evaluate<double>(3.0), 6.0);
  EXPECT_DOUBLE_EQ(inlet["doubled"].evaluate<double>(3), 6.0);
  EXPECT_DOUBLE_EQ(inlet["doubled"].evaluate<double>(x), 6.0);

  // A concrete value acts as a constant function of the arguments
  EXPECT_DOUBLE_EQ(inlet["constant"].evaluate<double>(3.0), 4.0);

  // Character strings are passed as std::string
  EXPECT_EQ(inlet["greeting"].evaluate<std::string>("inlet"), "hello inlet");
}

TEST(inlet_function, or_function_entry_is_present_in_either_form)
{
  auto inlet = createBasicInlet(R"(
    from_function = function() return 2.0 end
    from_value = 3.0
  )");
  inlet.addDouble("from_function").orFunction();
  inlet.addDouble("from_value").orFunction();
  inlet.addDouble("absent").orFunction();
  EXPECT_TRUE(inlet.verify());

  for(const std::string name : {"from_function", "from_value"})
  {
    EXPECT_TRUE(inlet.contains(name)) << name;
    EXPECT_TRUE(inlet.isUserProvided(name)) << name;
  }
  EXPECT_FALSE(inlet.contains("absent"));
  EXPECT_FALSE(inlet.isUserProvided("absent"));

  const auto& fields = inlet.getGlobalContainer().getChildFields();
  EXPECT_TRUE(fields.at("from_function")->isFunction());
  EXPECT_FALSE(fields.at("from_value")->isFunction());
  EXPECT_FALSE(fields.at("absent")->isFunction());
}

TEST(inlet_function, or_function_satisfies_required_in_either_form)
{
  for(const std::string input : {"scale = function() return 2.0 end", "scale = 2.0"})
  {
    auto inlet = createBasicInlet(input);
    inlet.addDouble("scale").orFunction().required();
    EXPECT_TRUE(inlet.verify()) << input;
  }

  auto missing = createBasicInlet("");
  missing.addDouble("scale").orFunction().required();
  EXPECT_FALSE(missing.verify());
}

TEST(inlet_function, or_function_default_applies_only_when_neither_form_is_supplied)
{
  // The default may be declared before or after orFunction()
  for(const bool defaultFirst : {true, false})
  {
    auto declare = [defaultFirst](Inlet& inlet) {
      auto& scale = inlet.addDouble("scale");
      if(defaultFirst)
      {
        scale.defaultValue(1.0).orFunction();
      }
      else
      {
        scale.orFunction().defaultValue(1.0);
      }
    };

    auto callback = createBasicInlet("scale = function() return 3.0 end");
    declare(callback);
    EXPECT_TRUE(callback.verify());
    EXPECT_EQ(callback["scale"].type(), InletType::Function) << defaultFirst;
    EXPECT_DOUBLE_EQ(callback["scale"].evaluate<double>(), 3.0) << defaultFirst;
    {
      // The default does not linger behind the function
      axom::slic::ScopedAbortToThrow abortGuard;
      EXPECT_THROW(callback.getGlobalContainer().getChildFields().at("scale")->get<double>(),
                   axom::slic::SlicAbortException)
        << defaultFirst;
      EXPECT_TRUE(slicErrorContains("was supplied as a function"));
    }

    auto concrete = createBasicInlet("scale = 2.0");
    declare(concrete);
    EXPECT_TRUE(concrete.verify());
    EXPECT_DOUBLE_EQ(concrete["scale"].evaluate<double>(), 2.0) << defaultFirst;

    auto absent = createBasicInlet("");
    declare(absent);
    EXPECT_TRUE(absent.verify());
    EXPECT_FALSE(absent.isUserProvided("scale"));
    EXPECT_DOUBLE_EQ(absent["scale"].evaluate<double>(), 1.0) << defaultFirst;
  }
}

TEST(inlet_function, or_function_reading_a_supplied_function_as_a_value_is_an_error)
{
  auto inlet = createBasicInlet("scale = function() return 2.0 end");
  inlet.addDouble("scale").orFunction();
  EXPECT_TRUE(inlet.verify());

  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(inlet["scale"].get<double>(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("was supplied as a function"));
  EXPECT_THROW(inlet.getGlobalContainer().getChildFields().at("scale")->get<double>(),
               axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("was supplied as a function"));
}

TEST(inlet_function, or_function_rejects_an_unrelated_wrong_type)
{
  auto inlet = createBasicInlet("scale = 'neither a number nor a function'");
  inlet.addDouble("scale").orFunction();
  EXPECT_FALSE(inlet.verify());
  EXPECT_FALSE(inlet.contains("scale"));
  EXPECT_FALSE(inlet.getGlobalContainer().getChildFields().at("scale")->isFunction());
}

TEST(inlet_function, or_function_rejects_unsupported_field_types)
{
  // FunctionTag has no Bool or Int return type
  auto inlet = createBasicInlet("");
  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(inlet.addBool("flag").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("but Field 'flag' has type"));
  EXPECT_THROW(inlet.addInt("count").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("but Field 'count' has type"));
}

TEST(inlet_function, or_function_redeclaration_must_match_the_signature)
{
  auto inlet = createBasicInlet("group = { scale = function(x) return 2 * x end }");
  auto& group = inlet.addStruct("group");
  group.addDouble("scale").orFunction({FunctionTag::Double});

  // Same signature: no-op. Redeclaring the Field returns the existing Field
  group.addDouble("scale").orFunction({FunctionTag::Double});
  EXPECT_TRUE(inlet.verify());
  EXPECT_DOUBLE_EQ(inlet["group/scale"].evaluate<double>(2.0), 4.0);

  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(group.addDouble("scale").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("different signature"));
}

TEST(inlet_function, or_function_is_recognized_by_strict_containers)
{
  auto inlet = createBasicInlet(R"(
    group = { scale = function() return 2.0 end, offset = 1.0 }
  )");
  auto& group = inlet.addStruct("group").strict();
  group.addDouble("scale").orFunction();
  group.addDouble("offset").orFunction();
  EXPECT_TRUE(inlet.verify());
  EXPECT_TRUE(inlet.unexpectedNames().empty());
  EXPECT_TRUE(group.exists());
  EXPECT_TRUE(group.isUserProvided());
}

TEST(inlet_function, or_function_constraints_apply_to_the_value_only)
{
  // Inlet does not call the function during verification
  auto callback = createBasicInlet("scale = function() return 5.0 end");
  callback.addDouble("scale").range(0.0, 1.0).orFunction();
  EXPECT_TRUE(callback.verify());

  auto concrete = createBasicInlet("scale = 5.0");
  concrete.addDouble("scale").range(0.0, 1.0).orFunction();
  EXPECT_FALSE(concrete.verify());
}

TEST(inlet_function, or_function_with_a_reader_that_has_no_functions)
{
  // A reader without function support only ever supplies the value
  auto reader = std::make_unique<axom::inlet::YAMLReader>();
  reader->parseString("scale: 2.5");
  Inlet inlet(std::move(reader));
  inlet.addDouble("scale").orFunction();
  EXPECT_TRUE(inlet.verify());
  EXPECT_EQ(inlet["scale"].type(), InletType::Double);
  EXPECT_DOUBLE_EQ(inlet["scale"].evaluate<double>(), 2.5);
}

TEST(inlet_function, or_function_array_selects_the_supplied_form)
{
  auto inlet = createBasicInlet(R"(
    from_function = function() return {1.0, 2.0, 3.0} end
    from_value = {4.0, 5.0}
  )");
  inlet.addDoubleArray("from_function").orFunction();
  inlet.addDoubleArray("from_value").orFunction();
  EXPECT_TRUE(inlet.verify());

  EXPECT_EQ(inlet["from_function"].type(), InletType::Function);
  EXPECT_EQ(inlet["from_value"].type(), InletType::Collection);

  EXPECT_EQ(inlet["from_function"].evaluate<std::vector<double>>(),
            (std::vector<double> {1.0, 2.0, 3.0}));
  EXPECT_EQ(inlet["from_value"].evaluate<std::vector<double>>(), (std::vector<double> {4.0, 5.0}));
  EXPECT_EQ(inlet["from_value"].get<std::vector<double>>(), (std::vector<double> {4.0, 5.0}));
}

TEST(inlet_function, or_function_array_evaluate_passes_arguments_to_the_function_only)
{
  auto inlet = createBasicInlet(R"(
    from_function = function(t) return {t, 2 * t} end
    from_value = {4.0, 5.0}
  )");
  inlet.addDoubleArray("from_function").orFunction({FunctionTag::Double});
  inlet.addDoubleArray("from_value").orFunction({FunctionTag::Double});
  EXPECT_TRUE(inlet.verify());

  EXPECT_EQ(inlet["from_function"].evaluate<std::vector<double>>(3),
            (std::vector<double> {3.0, 6.0}));
  EXPECT_EQ(inlet["from_value"].evaluate<std::vector<double>>(3), (std::vector<double> {4.0, 5.0}));
}

TEST(inlet_function, or_function_array_is_present_in_either_form)
{
  for(const std::string input :
      {"values = function() return {1.0} end", "values = {1.0}", "group = { other = 1.0 }"})
  {
    const bool supplied = input.rfind("values", 0) == 0;
    auto inlet = createBasicInlet(input);
    inlet.getGlobalContainer().strict();
    inlet.addDoubleArray("values").orFunction().required(supplied);
    inlet.addStruct("group").addDouble("other");
    EXPECT_TRUE(inlet.verify()) << input;
    EXPECT_TRUE(inlet.unexpectedNames().empty()) << input;
    EXPECT_EQ(inlet.contains("values"), supplied) << input;
    EXPECT_EQ(inlet.isUserProvided("values"), supplied) << input;
  }

  auto missing = createBasicInlet("");
  missing.addDoubleArray("values").orFunction().required();
  EXPECT_FALSE(missing.verify());
}

TEST(inlet_function, or_function_array_rejects_an_unrelated_wrong_type)
{
  auto inlet = createBasicInlet("values = 3.0");
  inlet.addDoubleArray("values").orFunction();
  EXPECT_FALSE(inlet.verify());
  EXPECT_NE(inlet["values"].type(), InletType::Function);
}

TEST(inlet_function, or_function_array_reading_a_supplied_function_as_values_is_an_error)
{
  auto inlet = createBasicInlet("values = function() return {1.0, 2.0} end");
  inlet.addDoubleArray("values").orFunction();
  EXPECT_TRUE(inlet.verify());

  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(inlet["values"].get<std::vector<double>>(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("'values' was supplied as a function"));
  // Reading through the Container, bypassing operator[], is caught as well
  EXPECT_THROW(
    inlet.getGlobalContainer().getChildContainers().at("values")->get<std::vector<double>>(),
    axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("'values' was supplied as a function"));
}

TEST(inlet_function, or_function_array_returns_at_most_three_components)
{
  // A function for an array returns FunctionTag::Vector
  auto inlet = createBasicInlet("values = function() return {1.0, 2.0, 3.0, 4.0} end");
  inlet.addDoubleArray("values").orFunction();
  EXPECT_TRUE(inlet.verify());

  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(inlet["values"].evaluate<std::vector<double>>(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("between 1 and 3"));
}

TEST(inlet_function, or_function_rejects_unsupported_collections)
{
  auto inlet = createBasicInlet("");
  axom::slic::ScopedAbortToThrow abortGuard;
  const std::string expected = "supports double fields, string fields, and arrays of doubles";

  EXPECT_THROW(inlet.addIntArray("ints").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains(expected));
  EXPECT_THROW(inlet.addStringArray("strings").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains(expected));
  EXPECT_THROW(inlet.addDoubleDictionary("dict").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains(expected));
  EXPECT_THROW(inlet.addVariantArray("mixed").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains(expected));
  EXPECT_THROW(inlet.addStruct("group").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains(expected));
  EXPECT_THROW(inlet.addFunction("func", FunctionTag::Double, {}).orFunction(),
               axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains(expected));
}

TEST(inlet_function, or_function_array_with_a_reader_that_has_no_functions)
{
  auto reader = std::make_unique<axom::inlet::YAMLReader>();
  reader->parseString("values: [1.0, 2.0]");
  Inlet inlet(std::move(reader));
  inlet.addDoubleArray("values").orFunction();
  EXPECT_TRUE(inlet.verify());
  EXPECT_EQ(inlet["values"].evaluate<std::vector<double>>(), (std::vector<double> {1.0, 2.0}));
}

TEST(inlet_function, or_function_is_documented_on_the_field)
{
  axom::utilities::filesystem::TempFile doc_file("inlet_or_function_docs", ".rst");
  {
    auto inlet = createBasicInlet(R"(
      scale = function(t) return 2.0 * t end
      label = 'fixed'
    )");
    inlet.addDouble("scale", "A scale").orFunction({FunctionTag::Double});
    inlet.addString("label", "A label.").orFunction();
    inlet.addDouble("plain", "A plain value");
    EXPECT_TRUE(inlet.verify());
    inlet.write(axom::inlet::SphinxWriter(doc_file.getPath()));
  }

  // Document the function form on the Field's row, regardless of which form the input supplied
  const std::string contents = doc_file.getFileContents();
  EXPECT_NE(std::string::npos,
            contents.find("A scale. Accepts a function with signature Double(Double)."));
  EXPECT_NE(std::string::npos, contents.find("A label. Accepts a function with signature String()."));
  EXPECT_NE(std::string::npos, contents.find("A plain value\n"));
  EXPECT_EQ(std::string::npos, contents.find("Functions"));
}

TEST(inlet_function, simple_void_to_double_through_container)
{
  std::string testString = "function foo () return 9.64 end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Double()>>();
  double result = callable();
  EXPECT_FLOAT_EQ(result, 9.64);
}

TEST(inlet_function, simple_double_to_void_through_container)
{
  // Test a function that returns nothing by using it to modify a global
  std::string testString = "bar = 19.9; function foo (a) bar = a end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Void, {FunctionTag::Double}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Void(FunctionType::Double)>>();
  double arg = -6.37;
  callable(arg);

  inlet.addDouble("bar", "bar's description");
  double result = inlet["bar"];
  EXPECT_FLOAT_EQ(result, arg);
}

TEST(inlet_function, simple_string_to_double_through_container)
{
  std::string testString =
    "function foo(s) "
    "  if s == 'a' then return 9.1 "
    "  elseif s == 'b' then return -6.3 "
    "  else return 66.5 end "
    "end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::String}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Double(FunctionType::String)>>();

  EXPECT_FLOAT_EQ(callable("a"), 9.1);
  EXPECT_FLOAT_EQ(callable("b"), -6.3);
  EXPECT_FLOAT_EQ(callable("c"), 66.5);
}

TEST(inlet_function, simple_double_to_string_through_container)
{
  std::string testString =
    "function foo(d) "
    "  if d == 1 then return 'a' "
    "  elseif d == 2 then return 'b' "
    "  else return 'c' end "
    "end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::String, {FunctionTag::Double}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::String(FunctionType::Double)>>();
  EXPECT_EQ(callable(1), "a");
  EXPECT_EQ(callable(2), "b");
  EXPECT_EQ(callable(3), "c");
}

TEST(inlet_function, simple_vec3_to_double_through_container_call)
{
  std::string testString = "function foo (v) return v.x + v.y + v.z end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::Vector}, "foo's description");

  auto result = inlet["foo"].call<double>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_through_container_call)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description");

  auto result = inlet["foo"].call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_double_to_double_through_container_call)
{
  std::string testString = "function foo (v, t) return t * (v.x + v.y + v.z) end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo",
                    FunctionTag::Double,
                    {FunctionTag::Vector, FunctionTag::Double},
                    "foo's description");

  auto result = inlet["foo"].call<double>(FunctionType::Vector {1, 2, 3}, 2.0);
  EXPECT_FLOAT_EQ(result, 12);
}

TEST(inlet_function, simple_vec3_double_to_vec3_through_container_call)
{
  std::string testString = "function foo (v, t) return t*v end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo",
                    FunctionTag::Vector,
                    {FunctionTag::Vector, FunctionTag::Double},
                    "foo's description");

  auto result = inlet["foo"].call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3}, 2.0);
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_pass)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto& func =
    inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description").required();
  func.registerVerifier([](const axom::inlet::Function& func) {
    auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  EXPECT_TRUE(inlet.verify());
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_fail)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto& func =
    inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description").required();
  func.registerVerifier([](const axom::inlet::Function& func) {
    auto result = func.call<FunctionType::Vector>(FunctionType::Vector {2, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  EXPECT_FALSE(inlet.verify());
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_with_errors_fail)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto& func =
    inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description").required();
  func.registerVerifier([](const axom::inlet::Function& func, std::vector<VerificationError>* errors) {
    INLET_VERIFICATION_WARNING("foo", "Something bad happened", errors);
    auto result = func.call<FunctionType::Vector>(FunctionType::Vector {2, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  ASSERT_FALSE(errors.empty());
  ASSERT_EQ(axom::Path("foo"), errors[0].path);
  ASSERT_EQ("Something bad happened", errors[0].message);
}

struct Foo
{
  bool bar;
  std::function<FunctionType::Vector(FunctionType::Vector)> baz;
};

template <>
struct FromInlet<Foo>
{
  Foo operator()(const axom::inlet::Container& base)
  {
    Foo f {base["bar"], base["baz"]};
    return f;
  }
};

TEST(inlet_function, simple_vec3_to_vec3_struct)
{
  std::string testString = "foo = { bar = true; baz = function (v) return 2*v end }";
  auto inlet = createBasicInlet(testString);

  // Define schema
  inlet.addBool("foo/bar", "bar's description");
  inlet.addFunction("foo/baz", FunctionTag::Vector, {FunctionTag::Vector}, "baz's description")
    .required();
  Foo foo = inlet["foo"].get<Foo>();
  EXPECT_TRUE(foo.bar);
  auto result = foo.baz({4, 5, 6});
  EXPECT_FLOAT_EQ(result[0], 8);
  EXPECT_FLOAT_EQ(result[1], 10);
  EXPECT_FLOAT_EQ(result[2], 12);
}

TEST(inlet_function, simple_vec3_to_vec3_array_of_struct)
{
  std::string testString =
    "foo = { [7] = { bar = true, "
    "                baz = function (v) return 2*v end }, "
    "       [12] = { bar = false, "
    "                baz = function (v) return 3*v end } "
    "}";
  auto inlet = createBasicInlet(testString);

  auto& arr_container = inlet.addStructArray("foo");

  // Define schema
  arr_container.addBool("bar", "bar's description");
  arr_container.addFunction("baz", FunctionTag::Vector, {FunctionTag::Vector}, "baz's description")
    .required();

  auto foos = inlet["foo"].get<std::unordered_map<int, Foo>>();
  EXPECT_TRUE(foos[7].bar);
  auto first_result = foos[7].baz({4, 5, 6});
  EXPECT_FLOAT_EQ(first_result[0], 8);
  EXPECT_FLOAT_EQ(first_result[1], 10);
  EXPECT_FLOAT_EQ(first_result[2], 12);

  EXPECT_FALSE(foos[12].bar);
  auto second_result = foos[12].baz({4, 5, 6});
  EXPECT_FLOAT_EQ(second_result[0], 12);
  EXPECT_FLOAT_EQ(second_result[1], 15);
  EXPECT_FLOAT_EQ(second_result[2], 18);
}

struct OperatorWithOrFunction
{
  double rotate;
};

template <>
struct FromInlet<OperatorWithOrFunction>
{
  OperatorWithOrFunction operator()(const axom::inlet::Container& base)
  {
    return {base["rotate"].evaluate<double>()};
  }
};

struct OperatorsWithOrFunctionDictionary
{
  std::unordered_map<std::string, OperatorWithOrFunction> operators;
};

template <>
struct FromInlet<OperatorsWithOrFunctionDictionary>
{
  OperatorsWithOrFunctionDictionary operator()(const axom::inlet::Container& base)
  {
    return {base["ops"].get<std::unordered_map<std::string, OperatorWithOrFunction>>()};
  }
};

TEST(inlet_function, or_function_applies_to_every_element_of_a_struct_array)
{
  // Elements may use different forms
  auto inlet = createBasicInlet(R"(
    operators = {
      { rotate = 30.0 },
      { rotate = function() return 45.0 end },
      { }
    }
  )");
  auto& ops = inlet.addStructArray("operators").strict();
  ops.addDouble("rotate").orFunction().defaultValue(0.0);
  EXPECT_TRUE(inlet.verify());
  EXPECT_TRUE(inlet.unexpectedNames().empty());

  const auto values = inlet["operators"].get<std::vector<OperatorWithOrFunction>>();
  ASSERT_EQ(values.size(), 3u);
  EXPECT_DOUBLE_EQ(values[0].rotate, 30.0);
  EXPECT_DOUBLE_EQ(values[1].rotate, 45.0);
  EXPECT_DOUBLE_EQ(values[2].rotate, 0.0);
}

TEST(inlet_function, or_function_in_nested_dictionary_of_struct)
{
  auto inlet = createBasicInlet(R"(
    groups = {
      [0] = {
        ops = {
          first = {rotate = 2},
          second = {rotate = function () return 3 end}
        }
      }
    }
  )");

  auto& groups = inlet.addStructArray("groups");
  auto& ops = groups.addStructDictionary("ops");
  ops.addDouble("rotate").orFunction().required();

  EXPECT_TRUE(inlet.verify());
  const auto values =
    inlet["groups"].get<std::unordered_map<int, OperatorsWithOrFunctionDictionary>>();
  EXPECT_DOUBLE_EQ(values.at(0).operators.at("first").rotate, 2.0);
  EXPECT_DOUBLE_EQ(values.at(0).operators.at("second").rotate, 3.0);
}

struct TranslateWithOrFunction
{
  std::vector<double> translate;
};

template <>
struct FromInlet<TranslateWithOrFunction>
{
  TranslateWithOrFunction operator()(const axom::inlet::Container& base)
  {
    return {base["translate"].evaluate<std::vector<double>>()};
  }
};

TEST(inlet_function, or_function_array_applies_to_every_element_of_a_struct_array)
{
  auto inlet = createBasicInlet(R"(
    operators = {
      { translate = {1.0, 2.0} },
      { translate = function() return {3.0, 4.0, 5.0} end }
    }
  )");
  inlet.addStructArray("operators").strict().addDoubleArray("translate").orFunction().required();
  EXPECT_TRUE(inlet.verify());
  EXPECT_TRUE(inlet.unexpectedNames().empty());

  const auto values = inlet["operators"].get<std::vector<TranslateWithOrFunction>>();
  ASSERT_EQ(values.size(), 2u);
  EXPECT_EQ(values[0].translate, (std::vector<double> {1.0, 2.0}));
  EXPECT_EQ(values[1].translate, (std::vector<double> {3.0, 4.0, 5.0}));
}

TEST(inlet_function, dimension_dependent_result)
{
  std::string testString =
    "foo = function (v) "
    "first = 2 * v.x "
    "last = 2 * v.y "
    "if v.dim == 2 then "
    "return Vector.new(first, last) "
    "else "
    "return Vector.new(first, 0, last) "
    "end "
    "end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Vector(FunctionType::Vector)>>();

  FunctionType::Vector input_3d({3.5, 0.5, 7.5});
  auto result = callable(input_3d);

  EXPECT_EQ(result.dim, 3);
  EXPECT_FLOAT_EQ(result.vec[0], 7);
  EXPECT_FLOAT_EQ(result.vec[1], 0);
  EXPECT_FLOAT_EQ(result.vec[2], 1);

  FunctionType::Vector input_2d({3.5, 0.5});
  result = callable(input_2d);

  EXPECT_EQ(result.dim, 2);
  EXPECT_FLOAT_EQ(result.vec[0], 7);
  EXPECT_FLOAT_EQ(result.vec[1], 1);
}

struct FooWithScalarFunc
{
  std::function<double(double)> bar;
};

template <>
struct FromInlet<FooWithScalarFunc>
{
  FooWithScalarFunc operator()(const axom::inlet::Container& base) { return {base["foo/bar"]}; }
};

TEST(inlet_function, nested_function_in_struct)
{
  std::string testString =
    "quux = { [0] = { foo = { bar = function (x) return x + 1 end } }, "
    "         [1] = { foo = { bar = function (x) return x + 3 end } } }";
  auto inlet = createBasicInlet(testString);

  auto& quux_schema = inlet.addStructArray("quux");
  auto& foo_schema = quux_schema.addStruct("foo");

  foo_schema.addFunction("bar", FunctionTag::Double, {FunctionTag::Double}, "bar's description");

  auto foos = inlet["quux"].get<std::vector<FooWithScalarFunc>>();
  EXPECT_EQ(foos.size(), 2);

  auto& first_func = foos[0].bar;
  // Check that the function object contains a valid target
  EXPECT_TRUE(static_cast<bool>(first_func));
  EXPECT_DOUBLE_EQ(first_func(4.0), 5.0);

  auto& second_func = foos[1].bar;
  // Check that the function object contains a valid target
  EXPECT_TRUE(static_cast<bool>(second_func));
  EXPECT_DOUBLE_EQ(second_func(4.0), 7.0);
}

template <typename Ret, typename... Args>
Ret checkedCall(const axom::sol::protected_function& func, Args&&... args)
{
  auto tentative_result = func(std::forward<Args>(args)...);
  EXPECT_TRUE(tentative_result.valid());
  return tentative_result;
}

/**
 * The inlet_function_usertype suite is intended to verify the correctness of the
 * definition of the correspondence between the FunctionType::Vector type and its
 * lua usertype equivalent.  Instead of using the Inlet interface to define and
 * access functions, the axom::sol::state member is interrogated directly via the
 * derived class SolStateReader to avoid mixing concerns in these tests.
 * 
 * Each entry in the Lua table/metatable for this usertype has a corresponding
 * test, i.e., one for each operator overload/constructor/member variable.
 */
TEST(inlet_function_usertype, lua_usertype_basic)
{
  std::string testString = "function func(vec) return 7 end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  int result = checkedCall<int>(func, vec);
  EXPECT_EQ(result, 7);
}

TEST(inlet_function_usertype, lua_usertype_basic_ret)
{
  std::string testString = "function func(x, y, z) return Vector.new(x, y, z) end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1, 2, 3);
  EXPECT_EQ(vec, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_ret_2d)
{
  std::string testString = "function func(x, y, z) return Vector.new(x, y) end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1, 2, 3);
  EXPECT_EQ(vec, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_ret_default)
{
  std::string testString = "function func(x, y, z) return Vector.new() end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {0, 0, 0};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1, 2, 3);
  EXPECT_EQ(vec, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_add)
{
  std::string testString = "function func(vec1, vec2) return vec1 + vec2 end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  const axom::inlet::FunctionType::Vector sum {5, 7, 9};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec1, vec2);
  EXPECT_EQ(result, sum);
}

TEST(inlet_function_usertype, lua_usertype_basic_sub)
{
  std::string testString = "function func(vec1, vec2) return vec1 - vec2 end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  const axom::inlet::FunctionType::Vector difference {-3, -3, -3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec1, vec2);
  EXPECT_EQ(result, difference);
}

TEST(inlet_function_usertype, lua_usertype_basic_negate)
{
  std::string testString = "function func(vec) return -vec end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const axom::inlet::FunctionType::Vector negated {-1, -2, -3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec);
  EXPECT_EQ(result, negated);
}

TEST(inlet_function_usertype, lua_usertype_basic_scalar_mult)
{
  std::string testString =
    "function func1(vec, x) return vec * x end; function func2(vec, x) return "
    "x * vec end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func1 = (*reader.solState())["func1"];
  axom::sol::protected_function func2 = (*reader.solState())["func2"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const axom::inlet::FunctionType::Vector doubled {2, 4, 6};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func1, vec, 2.0);
  EXPECT_EQ(result, doubled);
  const axom::inlet::FunctionType::Vector tripled {3, 6, 9};
  result = checkedCall<axom::inlet::FunctionType::Vector>(func2, vec, 3.0);
  EXPECT_EQ(result, tripled);
}

TEST(inlet_function_usertype, lua_usertype_basic_index_get)
{
  std::string testString = "function func(vec, idx) return vec[idx] end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  // Use 1-based indexing in these tests as lua is 1-indexed
  auto result = checkedCall<double>(func, vec, 1);
  EXPECT_FLOAT_EQ(1, result);
  result = checkedCall<double>(func, vec, 2);
  EXPECT_FLOAT_EQ(2, result);
  result = checkedCall<double>(func, vec, 3);
  EXPECT_FLOAT_EQ(3, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_index_set)
{
  std::string testString =
    "function func(idx) vec = Vector.new(1,1,1); vec[idx] = -1; return vec end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1);
  EXPECT_FLOAT_EQ(-1, result[0]);
  result = checkedCall<axom::inlet::FunctionType::Vector>(func, 2);
  EXPECT_FLOAT_EQ(-1, result[1]);
  result = checkedCall<axom::inlet::FunctionType::Vector>(func, 3);
  EXPECT_FLOAT_EQ(-1, result[2]);
}

TEST(inlet_function_usertype, lua_usertype_basic_norm)
{
  std::string testString = "function func(vec) return vec:norm() end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const double l2_norm = std::sqrt((1 * 1) + (2 * 2) + (3 * 3));
  auto result = checkedCall<double>(func, vec);
  EXPECT_FLOAT_EQ(l2_norm, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_squared_norm)
{
  std::string testString = "function func(vec) return vec:squared_norm() end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const double squared_l2_norm = (1 * 1) + (2 * 2) + (3 * 3);
  auto result = checkedCall<double>(func, vec);
  EXPECT_FLOAT_EQ(squared_l2_norm, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_unit_vec)
{
  std::string testString = "function func(vec) return vec:unitVector() end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const double l2_norm = std::sqrt((1 * 1) + (2 * 2) + (3 * 3));
  const axom::inlet::FunctionType::Vector unit {1 / l2_norm, 2 / l2_norm, 3 / l2_norm};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec);
  EXPECT_EQ(unit, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_dot)
{
  std::string testString = "function func(vec1, vec2) return vec1:dot(vec2) end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  const double dot = (1 * 4) + (2 * 5) + (3 * 6);
  auto result = checkedCall<double>(func, vec1, vec2);
  EXPECT_EQ(dot, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_cross)
{
  std::string testString = "function func(vec1, vec2) return vec1:cross(vec2) end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  const double i = (2 * 6) - (3 * 5);
  const double j = (3 * 4) - (1 * 6);
  const double k = (1 * 5) - (2 * 4);
  const axom::inlet::FunctionType::Vector cross {i, j, k};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec1, vec2);
  EXPECT_EQ(cross, result);
}

TEST(inlet_function_usertype, lua_usertype_check_dim)
{
  std::string testString = "function func(vec) return vec.dim end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5};
  auto result = checkedCall<double>(func, vec1);
  EXPECT_EQ(result, 3);

  result = checkedCall<double>(func, vec2);
  EXPECT_EQ(result, 2);
}

TEST(inlet_function_usertype, lua_usertype_named_access)
{
  std::string testString =
    "function func(vec, comp) if comp == 1 then return vec.x elseif comp == 2 "
    "then return vec.y else return vec.z end end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {4, 5, 6};
  auto result = checkedCall<double>(func, vec1, 1);
  EXPECT_EQ(result, 4);

  result = checkedCall<double>(func, vec1, 2);
  EXPECT_EQ(result, 5);

  result = checkedCall<double>(func, vec1, 3);
  EXPECT_EQ(result, 6);
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  axom::slic::addStreamToMsgLevel(new axom::slic::GenericOutputStream(&slicErrors(), "<MESSAGE>\n"),
                                  axom::slic::message::Error);

  result = RUN_ALL_TESTS();

  return result;
}
