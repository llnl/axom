// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/Inlet.hpp"

#include "gtest/gtest.h"

#include <array>
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
