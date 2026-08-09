// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/JSONSchemaWriter.hpp"
#include "axom/inlet/SphinxWriter.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

using axom::inlet::FunctionTag;
using axom::inlet::FunctionType;
using axom::inlet::Inlet;
using axom::inlet::InletType;
using axom::inlet::InputPath;
using axom::inlet::JSONSchemaWriter;
using axom::inlet::LuaReader;
using axom::inlet::SphinxWriter;
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

TEST(inlet_function, simple_vec3_to_vec3_raw_table_return)
{
  std::string testString = "function foo (v) return {v.x + 1, v.y + 2, v.z + 3} end";
  auto inlet = createBasicInlet(testString);

  auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {FunctionTag::Vector});

  EXPECT_TRUE(func);
  auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3});
  EXPECT_EQ(result.dim, 3);
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, vector_function_accepts_one_and_two_entry_table_returns)
{
  auto inlet = createBasicInlet(
    "function one () return {4.0} end\n"
    "function two () return {4.0, 5.0} end");

  auto one = inlet.reader().getFunction("one", FunctionTag::Vector, {});
  ASSERT_TRUE(one);
  const auto oneResult = one.call<FunctionType::Vector>();
  EXPECT_EQ(oneResult.dim, 1);
  EXPECT_FLOAT_EQ(oneResult[0], 4.0);

  auto two = inlet.reader().getFunction("two", FunctionTag::Vector, {});
  ASSERT_TRUE(two);
  const auto twoResult = two.call<FunctionType::Vector>();
  EXPECT_EQ(twoResult.dim, 2);
  EXPECT_FLOAT_EQ(twoResult[0], 4.0);
  EXPECT_FLOAT_EQ(twoResult[1], 5.0);
}

TEST(inlet_function, vector_function_rejects_scalar_return)
{
  auto inlet = createBasicInlet("function foo () return 2.0 end");
  auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {});

  EXPECT_TRUE(func);
  EXPECT_THROW(func.call<FunctionType::Vector>(), std::runtime_error);
}

TEST(inlet_function, vector_function_rejects_malformed_table_returns)
{
  // Lua vectors must be dense numeric sequences with a supported dimension.
  const std::array<std::string, 5> inputs {{
    "function foo () return {} end",
    "function foo () return {1, 2, 3, 4} end",
    "function foo () return {[1] = 1, [3] = 3} end",
    "function foo () return {1, 'two'} end",
    "function foo () return {1, 2, label = 3} end",
  }};

  for(const auto& input : inputs)
  {
    auto inlet = createBasicInlet(input);
    auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {});
    ASSERT_TRUE(func);
    EXPECT_THROW(func.call<FunctionType::Vector>(), std::runtime_error);
  }
}

TEST(inlet_function, scalar_and_string_functions_reject_wrong_return_types)
{
  auto inlet = createBasicInlet(
    "function scalar () return 'not a number' end\n"
    "function string () return {} end");

  auto scalar = inlet.reader().getFunction("scalar", FunctionTag::Double, {});
  ASSERT_TRUE(scalar);
  EXPECT_THROW(scalar.call<FunctionType::Double>(), std::runtime_error);

  auto string = inlet.reader().getFunction("string", FunctionTag::String, {});
  ASSERT_TRUE(string);
  EXPECT_THROW(string.call<FunctionType::String>(), std::runtime_error);
}

TEST(inlet_function, lua_callback_runtime_error_is_catchable)
{
  auto inlet = createBasicInlet("function foo () error('callback failed') end");
  auto func = inlet.reader().getFunction("foo", FunctionTag::Double, {});
  ASSERT_TRUE(func);

  try
  {
    func.call<FunctionType::Double>();
    FAIL() << "Expected the Lua callback to throw";
  }
  catch(const std::runtime_error& error)
  {
    EXPECT_NE(std::string(error.what()).find("callback failed"),
              std::string::npos);
  }
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

TEST(inlet_function, function_path_override)
{
  auto inlet = createBasicInlet("function public_name (x) return x + 2 end");

  auto& schema = inlet.addStruct("internal_group");
  schema.addFunction("internal_name",
                     FunctionTag::Double,
                     {FunctionTag::Double},
                     "",
                     "public_name");

  auto callback = inlet["internal_group/internal_name"].get<std::function<double(double)>>();
  EXPECT_DOUBLE_EQ(callback(3.0), 5.0);
}

TEST(inlet_function, explicit_exact_function_input_path)
{
  auto inlet = createBasicInlet("function public_name (x) return x + 2 end");

  inlet.addFunction("internal_name",
                    FunctionTag::Double,
                    {FunctionTag::Double},
                    InputPath::exact("public_name"));

  auto callback = inlet["internal_name"].get<std::function<double(double)>>();
  EXPECT_DOUBLE_EQ(callback(3.0), 5.0);
}

TEST(inlet_function, function_value_alternative_is_schema_order_independent)
{
  // Both schema entries inspect "foo"; declaration order must not select one.
  const auto addSchema = [](Inlet& inlet, bool functionFirst) {
    if(!functionFirst)
    {
      inlet.addDouble("foo");
    }
    inlet.addFunctionAsValueAlternative(
      "foo_callback",
      FunctionTag::Double,
      {},
      "foo");
    if(functionFirst)
    {
      inlet.addDouble("foo");
    }
  };

  for(const bool functionFirst : {true, false})
  {
    auto inlet = createBasicInlet("function foo () return 2.0 end");
    addSchema(inlet, functionFirst);

    EXPECT_TRUE(inlet.verify());
    EXPECT_FALSE(inlet.contains("foo"));
    ASSERT_TRUE(inlet.contains("foo_callback"));
    EXPECT_DOUBLE_EQ(inlet["foo_callback"].call<double>(), 2.0);
  }
}

TEST(inlet_function, function_value_alternative_preserves_concrete_value)
{
  auto inlet = createBasicInlet("foo = 4.0");
  inlet.addFunctionAsValueAlternative(
    "foo_callback",
    FunctionTag::Double,
    {},
    "foo");
  inlet.addDouble("foo");

  EXPECT_TRUE(inlet.verify());
  EXPECT_FALSE(inlet.contains("foo_callback"));
  ASSERT_TRUE(inlet.contains("foo"));
  EXPECT_DOUBLE_EQ(inlet["foo"].get<double>(), 4.0);
}

TEST(inlet_function, auto_named_function_value_alternative_uses_input_path)
{
  // set and access function in alternative
  {
    auto inlet = createBasicInlet("function foo () return 2.0 end");
    inlet.addFunctionAsValueAlternative(FunctionTag::Double, {}, "foo");
    inlet.addDouble("foo");

    EXPECT_TRUE(inlet.verify());
    auto& container = inlet.getGlobalContainer();
    EXPECT_TRUE(container.containsFunctionValueAlternative("foo"));

    const auto names = container.getFunctionValueAlternativeNames();
    ASSERT_EQ(names.size(), 1u);
    EXPECT_EQ(names[0], "foo");
    EXPECT_DOUBLE_EQ(
      container.getFunctionValueAlternative("foo").call<double>(),
      2.0);
  }

  // set and access value in alternative
  {
    auto concreteInlet = createBasicInlet("foo = 4.0");
    concreteInlet.addFunctionAsValueAlternative(FunctionTag::Double, {}, "foo");
    concreteInlet.addDouble("foo");

    EXPECT_TRUE(concreteInlet.verify());
    auto& concreteContainer = concreteInlet.getGlobalContainer();
    EXPECT_FALSE(concreteContainer.containsFunctionValueAlternative("foo"));

    EXPECT_TRUE(concreteContainer.getFunctionValueAlternativeNames().empty());
    EXPECT_DOUBLE_EQ(concreteInlet["foo"].get<double>(), 4.0);
  }
}

TEST(inlet_function, required_function_value_alternative_missing)
{
  auto inlet = createBasicInlet("");
  inlet.addDouble("foo");
  inlet
    .addFunctionAsValueAlternative(
      "foo_callback",
      FunctionTag::Double,
      {},
      "foo")
    .required();

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  EXPECT_FALSE(errors.empty());
  EXPECT_FALSE(inlet.contains("foo"));
  EXPECT_FALSE(inlet.contains("foo_callback"));
  EXPECT_FALSE(inlet.getGlobalContainer().exists());
}

TEST(inlet_function, function_value_alternative_rejects_unrelated_wrong_type)
{
  const auto addSchema = [](Inlet& inlet, bool functionFirst) {
    if(!functionFirst)
    {
      inlet.addDouble("foo");
    }
    inlet.addFunctionAsValueAlternative(
      "foo_callback",
      FunctionTag::Double,
      {},
      "foo");
    if(functionFirst)
    {
      inlet.addDouble("foo");
    }
  };

  for(const bool functionFirst : {true, false})
  {
    auto inlet = createBasicInlet("foo = 'not a number or function'");
    addSchema(inlet, functionFirst);

    EXPECT_FALSE(inlet.verify());
    EXPECT_FALSE(inlet.contains("foo"));
    EXPECT_FALSE(inlet.contains("foo_callback"));
    // The input exists even though neither schema entry accepts its type.
    EXPECT_TRUE(inlet.isUserProvided("foo"));
    EXPECT_FALSE(inlet.isUserProvided("foo_callback"));
    EXPECT_FALSE(inlet.getGlobalContainer().exists());
  }
}

TEST(inlet_function, function_value_alternative_is_valid_in_strict_container)
{
  for(const bool useFunction : {true, false})
  {
    auto inlet =
      createBasicInlet(useFunction ? "function foo () return 2.0 end" : "foo = 4.0");
    inlet.getGlobalContainer().strict();
    inlet.addFunctionAsValueAlternative(
      "foo_callback",
      FunctionTag::Double,
      {},
      "foo");
    inlet.addDouble("foo");

    EXPECT_TRUE(inlet.verify());
    EXPECT_TRUE(inlet.unexpectedNames().empty());
    EXPECT_EQ(inlet.contains("foo_callback"), useFunction);
    EXPECT_EQ(inlet.contains("foo"), !useFunction);
    EXPECT_TRUE(inlet.getGlobalContainer().exists());
  }
}

TEST(inlet_function, function_value_alternative_rejects_duplicate_callback_alternative)
{
  auto inlet = createBasicInlet("function scale () return 2.0 end");
  auto& callback = inlet.addFunctionAsValueAlternative(
    "scale_callback",
    FunctionTag::Double,
    {},
    "scale");
  auto& repeatedCallback = inlet.addFunctionAsValueAlternative(
    "scale_callback",
    FunctionTag::Double,
    {},
    "scale");

  // Like other Inlet schema entries, repeating the same name is idempotent.
  EXPECT_EQ(&callback, &repeatedCallback);

  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(inlet.addFunctionAsValueAlternative(
                 "another_scale_callback",
                 FunctionTag::Double,
                 {},
                 "scale"),
               axom::slic::SlicAbortException);
  EXPECT_EQ(inlet.getGlobalContainer().getChildFunctions().size(), 1u);
}

TEST(inlet_function, function_value_alternative_sidre_and_generated_documentation)
{
  const std::string sphinxFile = "inlet_function_value_alternative.rst";
  const std::string jsonFile = "inlet_function_value_alternative.json";
  auto inlet = createBasicInlet("scale = 2.0");
  inlet.addFunctionAsValueAlternative(
    "scale_callback",
    FunctionTag::Double,
    {},
    "scale");
  inlet.addDouble("scale", "Scale factor");

  // Sidre marks the internal callback schema group without storing the callable.
  const auto* callbackGroup = inlet.getGlobalContainer()
                                .getChildFunctions()
                                .at("scale_callback")
                                ->sidreGroup();
  ASSERT_TRUE(callbackGroup->hasView(axom::inlet::detail::FUNCTION_VALUE_ALTERNATIVE_FLAG));
  const std::int8_t alternativeFlag =
    callbackGroup->getView(axom::inlet::detail::FUNCTION_VALUE_ALTERNATIVE_FLAG)->getScalar();
  EXPECT_EQ(alternativeFlag, 1);

  inlet.write(SphinxWriter(sphinxFile));
  inlet.write(JSONSchemaWriter(jsonFile));

  const auto readFile = [](const std::string& path) {
    std::ifstream stream(path);
    return std::string(std::istreambuf_iterator<char>(stream),
                       std::istreambuf_iterator<char>());
  };
  const std::string sphinx = readFile(sphinxFile);
  const std::string json = readFile(jsonFile);
  std::remove(sphinxFile.c_str());
  std::remove(jsonFile.c_str());

  EXPECT_NE(sphinx.find("scale"), std::string::npos);
  EXPECT_EQ(sphinx.find("scale_callback"), std::string::npos);
  EXPECT_NE(json.find("scale"), std::string::npos);
  EXPECT_EQ(json.find("scale_callback"), std::string::npos);
}

TEST(inlet_function, function_value_alternative_restart_does_not_create_internal_field)
{
  axom::sidre::DataStore datastore;
  {
    auto reader = std::make_unique<LuaReader>();
    reader->parseString("scale = 2.0");
    Inlet inlet(std::move(reader), datastore.getRoot());
    inlet.addFunctionAsValueAlternative(FunctionTag::Double, {}, "scale");
    inlet.addDouble("scale");
  }

  {
    auto reader = std::make_unique<LuaReader>();
    Inlet restart(std::move(reader), datastore.getRoot(), true, true);

    EXPECT_TRUE(restart.contains("scale"));
    EXPECT_DOUBLE_EQ(restart.get<double>("scale"), 2.0);
    EXPECT_EQ(restart.getGlobalContainer().getChildFields().size(), 1u);
    EXPECT_TRUE(restart.getGlobalContainer().getChildFunctions().empty());
    EXPECT_FALSE(restart.getGlobalContainer().containsFunctionValueAlternative("scale"));
  }
}

TEST(inlet_function, returned_function_keeps_lua_state_alive)
{
  // An extracted callback must retain its Lua state after Inlet is destroyed.
  std::function<double(double)> callback;
  {
    auto inlet = createBasicInlet(
      "offset = 3.0; function foo (value) return value + offset end");
    inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::Double});
    callback = inlet["foo"].get<std::function<double(double)>>();
  }

  EXPECT_DOUBLE_EQ(callback(4.0), 7.0);
}

TEST(inlet_function, returned_functions_share_their_lua_state)
{
  std::function<double()> increment;
  std::function<double()> current;
  {
    auto inlet = createBasicInlet(
      "value = 0\n"
      "function increment () value = value + 1; return value end\n"
      "function current () return value end");
    inlet.addFunction("increment", FunctionTag::Double, {});
    inlet.addFunction("current", FunctionTag::Double, {});
    increment = inlet["increment"].get<std::function<double()>>();
    current = inlet["current"].get<std::function<double()>>();
  }

  EXPECT_DOUBLE_EQ(current(), 0.0);
  EXPECT_DOUBLE_EQ(increment(), 1.0);
  EXPECT_DOUBLE_EQ(current(), 1.0);
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

struct FooDictionary
{
  std::unordered_map<std::string, Foo> values;
};

template <>
struct FromInlet<FooDictionary>
{
  FooDictionary operator()(const axom::inlet::Container& base)
  {
    return {base["foo"].get<std::unordered_map<std::string, Foo>>()};
  }
};

struct FooWithValueAlternative
{
  std::function<double()> bar;
};

template <>
struct FromInlet<FooWithValueAlternative>
{
  FooWithValueAlternative operator()(const axom::inlet::Container& base)
  {
    return {base["bar_callback"]};
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

TEST(inlet_function, function_path_override_in_array_of_struct)
{
  std::string testString =
    "foo = { [7] = { bar = true, "
    "                callback = function (v) return 2*v end }, "
    "       [12] = { bar = false, "
    "                 callback = function (v) return 3*v end } "
    "}";
  auto inlet = createBasicInlet(testString);

  auto& arr_container = inlet.addStructArray("foo");
  arr_container.addBool("bar");
  arr_container.addFunction("baz",
                            FunctionTag::Vector,
                            {FunctionTag::Vector},
                            "",
                            "callback");

  auto foos = inlet["foo"].get<std::unordered_map<int, Foo>>();
  EXPECT_FLOAT_EQ(foos[7].baz({1, 2, 3})[0], 2);
  EXPECT_FLOAT_EQ(foos[12].baz({1, 2, 3})[0], 3);
}

TEST(inlet_function, explicit_exact_function_input_path_in_array_of_struct)
{
  auto inlet = createBasicInlet(
    "shared_callback = function (v) return 4*v end; "
    "foo = { [7] = { bar = true }, [12] = { bar = false } }");

  auto& arr_container = inlet.addStructArray("foo");
  arr_container.addBool("bar");
  arr_container.addFunction("baz",
                            FunctionTag::Vector,
                            {FunctionTag::Vector},
                            InputPath::exact("shared_callback"));

  auto foos = inlet["foo"].get<std::unordered_map<int, Foo>>();
  EXPECT_FLOAT_EQ(foos[7].baz({1, 2, 3})[0], 4);
  EXPECT_FLOAT_EQ(foos[12].baz({1, 2, 3})[0], 4);
}

TEST(inlet_function, explicit_relative_function_input_path_in_nested_dictionary_of_struct)
{
  auto inlet = createBasicInlet(
    "groups = { "
    "  [0] = { foo = { first = { bar = true, "
    "                            callback = function (v) return 2*v end }, "
    "                  second = { bar = false, "
    "                             callback = function (v) return 3*v end } } }, "
    "  [1] = { foo = { third = { bar = true, "
    "                            callback = function (v) return 4*v end } } } }");

  auto& group_container = inlet.addStructArray("groups");
  auto& dict_container = group_container.addStructDictionary("foo");
  dict_container.addBool("bar");
  // Resolve "callback" from each dictionary value, not the enclosing schema.
  dict_container.addFunction(
    "baz",
    FunctionTag::Vector,
    {FunctionTag::Vector},
    InputPath::relativeToCollectionElement("callback"));

  auto groups = inlet["groups"].get<std::unordered_map<int, FooDictionary>>();
  EXPECT_FLOAT_EQ(groups[0].values["first"].baz({1, 2, 3})[0], 2);
  EXPECT_FLOAT_EQ(groups[0].values["second"].baz({1, 2, 3})[0], 3);
  EXPECT_FLOAT_EQ(groups[1].values["third"].baz({1, 2, 3})[0], 4);
}

TEST(inlet_function, explicit_relative_function_value_alternative_in_array_of_struct)
{
  auto inlet = createBasicInlet(
    "foo = { [7] = { bar = function () return 2 end }, "
    "        [12] = { bar = function () return 3 end } }");

  auto& arr_container = inlet.addStructArray("foo");
  arr_container.addDouble("bar");
  arr_container.addFunctionAsValueAlternative(
    "bar_callback",
    FunctionTag::Double,
    {},
    InputPath::relativeToCollectionElement("bar"));

  EXPECT_TRUE(inlet.verify());
  auto foos =
    inlet["foo"].get<std::unordered_map<int, FooWithValueAlternative>>();
  EXPECT_DOUBLE_EQ(foos[7].bar(), 2.0);
  EXPECT_DOUBLE_EQ(foos[12].bar(), 3.0);
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

TEST(inlet_function, explicit_relative_function_input_path_in_nested_struct)
{
  std::string testString =
    "quux = { [0] = { foo = { callback = function (x) return x + 1 end } }, "
    "         [1] = { foo = { callback = function (x) return x + 3 end } } }";
  auto inlet = createBasicInlet(testString);

  auto& quux_schema = inlet.addStructArray("quux");
  auto& foo_schema = quux_schema.addStruct("foo");
  foo_schema.addFunction("bar",
                         FunctionTag::Double,
                         {FunctionTag::Double},
                         InputPath::relativeToCollectionElement("callback"));

  auto foos = inlet["quux"].get<std::vector<FooWithScalarFunc>>();
  EXPECT_DOUBLE_EQ(foos[0].bar(4.0), 5.0);
  EXPECT_DOUBLE_EQ(foos[1].bar(4.0), 7.0);
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

  result = RUN_ALL_TESTS();

  return result;
}
