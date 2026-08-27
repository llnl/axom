// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/SphinxWriter.hpp"

#include "gtest/gtest.h"

#include <array>
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

TEST(inlet_function, vector_function_rejects_malformed_table_returns)
{
  // Lua vectors must be dense numeric sequences with a supported dimension.
  const std::array<std::string, 6> results {{
    "2.0",
    "{}",
    "{1, 2, 3, 4}",
    "{[1] = 1, [3] = 3}",
    "{1, 'two'}",
    "{1, 2, label = 3}",
  }};

  for(const auto& result : results)
  {
    auto inlet = createBasicInlet("function foo () return " + result + " end");
    auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {});
    ASSERT_TRUE(func);
    EXPECT_THROW(func.call<FunctionType::Vector>(), axom::inlet::InletError);
  }
}

TEST(inlet_function, lua_callback_failures_are_catchable)
{
  auto inlet = createBasicInlet(R"(
    function runtime_error () error('callback failed') end
    function wrong_type () return 'not a number' end
  )");

  auto wrongType = inlet.reader().getFunction("wrong_type", FunctionTag::Double, {});
  ASSERT_TRUE(wrongType);
  EXPECT_THROW(wrongType.call<FunctionType::Double>(), axom::inlet::InletError);

  auto runtimeError = inlet.reader().getFunction("runtime_error", FunctionTag::Double, {});
  ASSERT_TRUE(runtimeError);

  try
  {
    runtimeError.call<FunctionType::Double>();
    FAIL() << "Expected the Lua callback to throw";
  }
  catch(const axom::inlet::InletError& error)
  {
    EXPECT_NE(std::string(error.what()).find("callback failed"), std::string::npos);
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

TEST(inlet_function, function_value_alternative_selects_supplied_representation)
{
  auto inlet = createBasicInlet(R"(
    function label () return 'computed' end
    scale = 4.0
  )");

  inlet.addFunctionAsValueAlternative("label", FunctionTag::String, {});
  inlet.addString("label");
  inlet.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  inlet.addDouble("scale");

  EXPECT_TRUE(inlet.verify());
  auto& container = inlet.getGlobalContainer();
  EXPECT_FALSE(inlet.contains("label"));
  ASSERT_TRUE(container.containsFunctionValueAlternative("label"));
  EXPECT_EQ(container.getFunctionValueAlternative("label").call<std::string>(), "computed");

  EXPECT_TRUE(inlet.contains("scale"));
  EXPECT_FALSE(container.containsFunctionValueAlternative("scale"));
  EXPECT_DOUBLE_EQ(inlet["scale"].get<double>(), 4.0);
}

TEST(inlet_function, function_value_alternative_rejects_unrelated_wrong_type)
{
  auto inlet = createBasicInlet("foo = 'not a number or function'");
  inlet.addFunctionAsValueAlternative("foo", FunctionTag::Double, {});
  inlet.addDouble("foo");

  EXPECT_FALSE(inlet.verify());
  EXPECT_FALSE(inlet.contains("foo"));
  EXPECT_FALSE(inlet.getGlobalContainer().containsFunctionValueAlternative("foo"));
  // The input exists even though neither schema entry accepts its type.
  EXPECT_TRUE(inlet.isUserProvided("foo"));
  EXPECT_FALSE(inlet.getGlobalContainer().exists());
}

TEST(inlet_function, function_value_alternative_is_container_independent)
{
  for(const bool functionOnRoot : {true, false})
  {
    auto inlet = createBasicInlet("group = { value = function() return 2.0 end }");
    inlet.getGlobalContainer().strict();
    auto& group = inlet.addStruct("group");

    // The alternative and the concrete entry may be declared through different
    // Containers, so long as the alternative comes first.
    if(functionOnRoot)
    {
      inlet.addFunctionAsValueAlternative("group/value", FunctionTag::Double, {});
      group.addDouble("value");
    }
    else
    {
      group.addFunctionAsValueAlternative("value", FunctionTag::Double, {});
      inlet.addDouble("group/value");
    }

    EXPECT_TRUE(inlet.verify());
    EXPECT_TRUE(inlet.unexpectedNames().empty());
    EXPECT_TRUE(inlet.getGlobalContainer().containsFunctionValueAlternative("group/value"));
    EXPECT_TRUE(group.containsFunctionValueAlternative("value"));
    EXPECT_FALSE(inlet.contains("group/value"));
    EXPECT_TRUE(inlet.isUserProvided("group/value"));
    EXPECT_TRUE(group.isUserProvided("value"));
    EXPECT_TRUE(group.exists());
    EXPECT_DOUBLE_EQ(group.getFunctionValueAlternative("value").call<double>(), 2.0);
    EXPECT_EQ(group.getFunctionValueAlternativeNames(), std::vector<std::string> {"value"});
    EXPECT_TRUE(inlet.getGlobalContainer().getFunctionValueAlternativeNames().empty());
  }
}

TEST(inlet_function, function_value_alternative_names_are_sorted)
{
  auto inlet = createBasicInlet(R"(
    zeta = function() return 1.0 end
    alpha = function() return 2.0 end
    plain = 3.0
  )");

  inlet.addFunctionAsValueAlternative("zeta", FunctionTag::Double, {});
  inlet.addFunctionAsValueAlternative("alpha", FunctionTag::Double, {});
  // Declared but not supplied by the input, so it is not reported
  inlet.addFunctionAsValueAlternative("plain", FunctionTag::Double, {});
  inlet.addDouble("plain");

  EXPECT_TRUE(inlet.verify());
  EXPECT_EQ(inlet.getGlobalContainer().getFunctionValueAlternativeNames(),
            (std::vector<std::string> {"alpha", "zeta"}));
}

TEST(inlet_function, function_value_alternative_array_is_container_independent)
{
  for(const bool functionOnRoot : {true, false})
  {
    auto inlet = createBasicInlet("group = { values = function() return {1.0, 2.0, 3.0} end }");
    auto& group = inlet.addStruct("group");

    if(functionOnRoot)
    {
      inlet.addFunctionAsValueAlternative("group/values", FunctionTag::Vector, {});
    }
    else
    {
      group.addFunctionAsValueAlternative("values", FunctionTag::Vector, {});
    }
    inlet.addDoubleArray("group/values");

    EXPECT_TRUE(inlet.verify());
    EXPECT_FALSE(inlet.contains("group/values"));
    EXPECT_TRUE(group.containsFunctionValueAlternative("values"));
    const auto result = group.getFunctionValueAlternative("values").call<FunctionType::Vector>();
    EXPECT_DOUBLE_EQ(result[0], 1.0);
    EXPECT_DOUBLE_EQ(result[1], 2.0);
    EXPECT_DOUBLE_EQ(result[2], 3.0);
  }
}

TEST(inlet_function, function_value_alternative_declaration_is_idempotent)
{
  // Matches addFunction's behavior: redeclaring returns the existing entry
  auto inlet = createBasicInlet("group = { scale = function() return 2.0 end }");
  auto& group = inlet.addStruct("group");

  auto& first = inlet.addFunctionAsValueAlternative("group/scale", FunctionTag::Double, {});
  auto& second = group.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  EXPECT_EQ(&first, &second);

  EXPECT_EQ(1u, group.getChildFunctions().size());
  EXPECT_TRUE(inlet.getGlobalContainer().getChildFunctions().empty());
  EXPECT_DOUBLE_EQ(group.getFunctionValueAlternative("scale").call<double>(), 2.0);
}

TEST(inlet_function, function_value_alternative_is_not_documented_under_its_schema_name)
{
  // The alternative is stored under an internal schema name,
  // which must not leak into generated documentation
  const std::string docFile = "inlet_function_value_alternative_docs.rst";
  {
    auto inlet = createBasicInlet("scale = function() return 2.0 end");
    inlet.addFunctionAsValueAlternative("scale", FunctionTag::Double, {}, "a scale callback");
    inlet.addDouble("scale", "a scale");
    EXPECT_TRUE(inlet.verify());
    // The alternative is a real schema entry, just under an internal name
    EXPECT_EQ(1u, inlet.getGlobalContainer().getChildFunctions().size());
    inlet.write(axom::inlet::SphinxWriter(docFile));
  }

  std::ifstream stream(docFile);
  ASSERT_TRUE(stream.good());
  const std::string contents {std::istreambuf_iterator<char> {stream},
                              std::istreambuf_iterator<char> {}};
  EXPECT_EQ(std::string::npos, contents.find("_inlet_function_alternative"));
  EXPECT_NE(std::string::npos, contents.find("scale"));
}

TEST(inlet_function, function_value_alternative_can_be_required)
{
  // The alternative is an ordinary Function, so schema constraints apply to it
  auto missing = createBasicInlet("scale = 2.0");
  missing.addFunctionAsValueAlternative("scale", FunctionTag::Double, {}).required();
  missing.addDouble("scale");
  EXPECT_FALSE(missing.verify());

  auto supplied = createBasicInlet("scale = function() return 2.0 end");
  supplied.addFunctionAsValueAlternative("scale", FunctionTag::Double, {}).required();
  supplied.addDouble("scale");
  EXPECT_TRUE(supplied.verify());
}

TEST(inlet_function, function_value_alternative_rejects_invalid_declaration)
{
  auto inlet = createBasicInlet("");
  axom::slic::ScopedAbortToThrow abortGuard;

  EXPECT_THROW(inlet.addFunctionAsValueAlternative("", FunctionTag::Double, {}),
               axom::slic::SlicAbortException);
  EXPECT_THROW(inlet.addFunctionAsValueAlternative("value", FunctionTag::Void, {}),
               axom::slic::SlicAbortException);
  EXPECT_TRUE(inlet.getGlobalContainer().getFunctionValueAlternativeNames().empty());
}

TEST(inlet_function, function_value_alternative_rejects_declaration_after_the_value)
{
  // Declaring the concrete entry first cannot work: it has already been read and
  // marked as being of the wrong type. Without this check the schema author sees
  // a verification failure blaming the input instead of the schema.
  axom::slic::ScopedAbortToThrow abortGuard;

  auto scalar = createBasicInlet("scale = function() return 2.0 end");
  scalar.addDouble("scale");
  EXPECT_THROW(scalar.addFunctionAsValueAlternative("scale", FunctionTag::Double, {}),
               axom::slic::SlicAbortException);

  auto collection = createBasicInlet("values = function() return {1.0, 2.0} end");
  collection.addDoubleArray("values");
  EXPECT_THROW(collection.addFunctionAsValueAlternative("values", FunctionTag::Vector, {}),
               axom::slic::SlicAbortException);

  // Declaring through a parent Container is rejected the same way
  auto nested = createBasicInlet("group = { value = function() return 2.0 end }");
  nested.addStruct("group").addDouble("value");
  EXPECT_THROW(nested.addFunctionAsValueAlternative("group/value", FunctionTag::Double, {}),
               axom::slic::SlicAbortException);
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

struct FooWithValueAlternative
{
  double bar;
};

template <>
struct FromInlet<FooWithValueAlternative>
{
  FooWithValueAlternative operator()(const axom::inlet::Container& base)
  {
    if(base.containsFunctionValueAlternative("bar"))
    {
      return {base.getFunctionValueAlternative("bar").call<double>()};
    }
    return {base["bar"].get<double>()};
  }
};

struct FooWithValueAlternativeDictionary
{
  std::unordered_map<std::string, FooWithValueAlternative> values;
};

template <>
struct FromInlet<FooWithValueAlternativeDictionary>
{
  FooWithValueAlternativeDictionary operator()(const axom::inlet::Container& base)
  {
    return {base["foo"].get<std::unordered_map<std::string, FooWithValueAlternative>>()};
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

TEST(inlet_function, function_value_alternative_in_nested_dictionary_of_struct)
{
  // Both representations appear in the same collection, so the alternative must
  // be expanded across the collection's elements exactly as the concrete entry is
  auto inlet = createBasicInlet(R"(
    groups = {
      [0] = {
        foo = {
          first = {bar = 2},
          second = {bar = function () return 3 end}
        }
      }
    }
  )");

  auto& groups = inlet.addStructArray("groups");
  auto& foos = groups.addStructDictionary("foo");
  foos.addFunctionAsValueAlternative("bar", FunctionTag::Double, {});
  foos.addDouble("bar");

  EXPECT_TRUE(inlet.verify());
  const auto values =
    inlet["groups"].get<std::unordered_map<int, FooWithValueAlternativeDictionary>>();
  EXPECT_DOUBLE_EQ(values.at(0).values.at("first").bar, 2.0);
  EXPECT_DOUBLE_EQ(values.at(0).values.at("second").bar, 3.0);
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

  result = RUN_ALL_TESTS();

  return result;
}
