// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file LuaReader.cpp
 *
 * \brief This file contains the class implementation of the LuaReader.
 *******************************************************************************
 */

#include <fstream>

#include "axom/inlet/LuaReader.hpp"

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/inlet/inlet_utils.hpp"

#include "axom/fmt.hpp"
#include "axom/slic.hpp"

#include "axom/sol.hpp"

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

namespace axom
{
namespace inlet
{
namespace detail
{
/*!
 *******************************************************************************
 * \brief Extracts an object from sol into a concrete type, implemented to support
 * extracting to a VariantKey
 * 
 * \tparam T The type to extract to
 *******************************************************************************
 */
template <typename T>
T extractAs(const axom::sol::object& obj)
{
  // By default, just ask sol to cast it
  return obj.as<T>();
}
/// \overload
template <>
VariantKey extractAs(const axom::sol::object& obj)
{
  // FIXME: Floating-point indices?
  if(obj.get_type() == axom::sol::type::number)
  {
    return obj.as<int>();
  }
  else
  {
    return obj.as<std::string>();
  }
}

/*!
 *******************************************************************************
 * \brief Recursive name retrieval function - adds the names of all descendents
 * of \a node as an Inlet-style path
 * 
 * \param [in] ignores The vector of paths to ignore, used for pre-loaded entries
 * in Lua's global table
 * \param [in] table The Lua table to "visit"
 * \param [in] prefix The Inlet-style path to \a table relative to the "root" of
 * the input file
 * \param [out] names The vector of paths to add to
 *******************************************************************************
 */
void nameRetrievalHelper(const std::vector<std::string>& ignores,
                         const axom::sol::table& table,
                         const std::string& prefix,
                         std::vector<std::string>& names)
{
  auto toString = [](const VariantKey& key) {
    return key.type() == InletType::String ? static_cast<std::string>(key)
                                           : std::to_string(static_cast<int>(key));
  };
  for(const auto& entry : table)
  {
    const auto variantKey = detail::extractAs<VariantKey>(entry.first);
    const std::string fullName = utilities::string::appendPrefix(prefix, toString(variantKey));
    if(std::find(ignores.begin(), ignores.end(), fullName) == ignores.end())
    {
      names.push_back(fullName);
      if(entry.second.get_type() == axom::sol::type::table && (ignores.back() != fullName))
      {
        nameRetrievalHelper(ignores, entry.second, fullName, names);
      }
    }
  }
}

}  // end namespace detail

LuaReader::LuaReader()
{
  m_lua = std::make_shared<axom::sol::state>();
  m_lua->open_libraries(axom::sol::lib::base,
                        axom::sol::lib::math,
                        axom::sol::lib::string,
                        axom::sol::lib::package);
  auto vec_type = m_lua->new_usertype<FunctionType::Vector>(
    "Vector",  // Name of the class in Lua
    // Add make_vector as a constructor to enable "new Vector(x,y,z)"
    // Use lambdas for 2D and "default" cases - default arguments cannot be
    // propagated automatically
    "new",
    axom::sol::factories([](double x, double y, double z) { return FunctionType::Vector {x, y, z}; },
                         [](double x, double y) { return FunctionType::Vector {x, y}; },
                         // Assume three for a default constructor
                         [] { return FunctionType::Vector {}; }),
    // Add vector addition operation
    axom::sol::meta_function::addition,
    [](const FunctionType::Vector& u, const FunctionType::Vector& v) {
      SLIC_ASSERT_MSG(u.dim == v.dim,
                      "[Inlet] Operands to InletVector addition are of different dimension");
      return FunctionType::Vector {u.vec + v.vec, u.dim};
    },
    axom::sol::meta_function::subtraction,
    [](const FunctionType::Vector& u, const FunctionType::Vector& v) {
      SLIC_ASSERT_MSG(u.dim == v.dim,
                      "[Inlet] Operands to InletVector subtraction are of "
                      "different dimension");
      return FunctionType::Vector {u.vec - v.vec, u.dim};
    },
    // Needs to be resolved in the same way as operator+
    axom::sol::meta_function::unary_minus,
    [](const FunctionType::Vector& u) { return FunctionType::Vector {-u.vec, u.dim}; },
    // To allow both "directions" of a scalar multiplication, the overloads
    // have to be manually specified + resolved
    axom::sol::meta_function::multiplication,
    axom::sol::overload([](const FunctionType::Vector& u,
                           const double a) { return FunctionType::Vector {a * u.vec, u.dim}; },
                        [](const double a, const FunctionType::Vector& u) {
                          return FunctionType::Vector {a * u.vec, u.dim};
                        }),
    // Separate functions from get/set via index - subtract 1 as lua is 1-indexed
    axom::sol::meta_function::index,
    [](const FunctionType::Vector& vec, const int key) { return vec[key - 1]; },
    // A lambda is used here as the set-via-returned reference is insufficient
    axom::sol::meta_function::new_index,
    [](FunctionType::Vector& vec, const int key, const double value) { vec[key - 1] = value; },
    // Set up the mathematical operations by name
    "norm",
    [](const FunctionType::Vector& u) { return u.vec.norm(); },
    "squared_norm",
    [](const FunctionType::Vector& u) { return u.vec.squared_norm(); },
    "unitVector",
    [](const FunctionType::Vector& u) { return FunctionType::Vector {u.vec.unitVector(), u.dim}; },
    "dot",
    [](const FunctionType::Vector& u, const FunctionType::Vector& v) {
      SLIC_ASSERT_MSG(u.dim == v.dim,
                      "[Inlet] Operands to InletVector dot product are of "
                      "different dimension");
      return u.vec.dot(v.vec);
    },
    // Implemented as a member function like dot products are, for simplicity
    "cross",
    // Needs to be resolved as it is an overloaded static method
    [](const FunctionType::Vector& u, const FunctionType::Vector& v) {
      SLIC_ASSERT_MSG(u.dim == v.dim,
                      "[Inlet] Operands to InletVector cross product are of "
                      "different dimension");
      return FunctionType::Vector {primal::Vector3D::cross_product(u.vec, v.vec), u.dim};
    },
    "dim",
    axom::sol::property([](const FunctionType::Vector& u) { return u.dim; }),
    "x",
    axom::sol::property([](const FunctionType::Vector& u) { return u.vec[0]; }),
    "y",
    axom::sol::property([](const FunctionType::Vector& u) { return u.vec[1]; }),
    "z",
    axom::sol::property([](const FunctionType::Vector& u) { return u.vec[2]; }));

  // Pass the preloaded globals as both the set to ignore and the set to add
  // to, such that only the top-level preloaded globals are added
  detail::nameRetrievalHelper(m_preloaded_globals, m_lua->globals(), "", m_preloaded_globals);
}

bool LuaReader::parseFile(const std::string& filePath)
{
  if(!axom::utilities::filesystem::pathExists(filePath))
  {
    SLIC_WARNING(fmt::format("Inlet: Given Lua input file does not exist: {0}", filePath));
    return false;
  }

  auto script = m_lua->script_file(filePath);
  if(!script.valid())
  {
    SLIC_WARNING(fmt::format("Inlet: Given Lua input file is invalid: {0}", filePath));
  }
  return script.valid();
}

bool LuaReader::parseString(const std::string& luaString)
{
  if(luaString.empty())
  {
    SLIC_WARNING("Inlet: Given an empty Lua string to parse.");
    return false;
  }
  m_lua->script(luaString);
  return true;
}

// TODO allow alternate delimiter at sidre level
#define SCOPE_DELIMITER '/'

ReaderResult LuaReader::getBool(const std::string& id, bool& value) { return getValue(id, value); }

ReaderResult LuaReader::getDouble(const std::string& id, double& value)
{
  return getValue(id, value);
}

ReaderResult LuaReader::getInt(const std::string& id, int& value) { return getValue(id, value); }

ReaderResult LuaReader::getString(const std::string& id, std::string& value)
{
  return getValue(id, value);
}

ReaderResult LuaReader::getIntMap(const std::string& id, std::unordered_map<int, int>& values)
{
  return getMap(id, values, axom::sol::type::number);
}

ReaderResult LuaReader::getDoubleMap(const std::string& id, std::unordered_map<int, double>& values)
{
  return getMap(id, values, axom::sol::type::number);
}

ReaderResult LuaReader::getBoolMap(const std::string& id, std::unordered_map<int, bool>& values)
{
  return getMap(id, values, axom::sol::type::boolean);
}

ReaderResult LuaReader::getStringMap(const std::string& id,
                                     std::unordered_map<int, std::string>& values)
{
  return getMap(id, values, axom::sol::type::string);
}

ReaderResult LuaReader::getIntMap(const std::string& id, std::unordered_map<VariantKey, int>& values)
{
  return getMap(id, values, axom::sol::type::number);
}

ReaderResult LuaReader::getDoubleMap(const std::string& id,
                                     std::unordered_map<VariantKey, double>& values)
{
  return getMap(id, values, axom::sol::type::number);
}

ReaderResult LuaReader::getBoolMap(const std::string& id, std::unordered_map<VariantKey, bool>& values)
{
  return getMap(id, values, axom::sol::type::boolean);
}

ReaderResult LuaReader::getStringMap(const std::string& id,
                                     std::unordered_map<VariantKey, std::string>& values)
{
  return getMap(id, values, axom::sol::type::string);
}

template <typename Iter>
bool LuaReader::traverseToTable(Iter begin, Iter end, axom::sol::table& table)
{
  // Nothing to traverse
  if(begin == end)
  {
    return true;
  }

  if(!(*m_lua)[*begin].valid())
  {
    return false;
  }

  // Use the first one to index into the global lua state
  table = (*m_lua)[*begin];
  ++begin;

  // Then use the remaining keys to walk down to the requested table
  for(auto curr = begin; curr != end; ++curr)
  {
    auto key = *curr;
    bool is_int = conduit::utils::string_is_integer(key);
    int key_as_int = conduit::utils::string_to_value<int>(key);
    if(is_int && table[key_as_int].valid())
    {
      table = table[key_as_int];
    }
    else if(table[key].valid())
    {
      table = table[key];
    }
    else
    {
      return false;
    }
  }
  return true;
}

ReaderResult LuaReader::getIndices(const std::string& id, std::vector<int>& indices)
{
  return getIndicesInternal(id, indices);
}

ReaderResult LuaReader::getIndices(const std::string& id, std::vector<VariantKey>& indices)
{
  return getIndicesInternal(id, indices);
}

// A set of pure functions for handling the conversion of Lua functions to C++
// callables
namespace detail
{
/*!
 *****************************************************************************
 * \brief Templated function for calling a sol function
 *
 * \param [in] func The sol function of unknown concrete type
 * \tparam Args The argument types of the function
 *
 * \return A checkable version of the function's result
 *****************************************************************************
 */
template <typename... Args>
axom::sol::protected_function_result callWith(const axom::sol::protected_function& func,
                                              Args&&... args)
{
  auto tentative_result = func(std::forward<Args>(args)...);
  SLIC_ERROR_IF(!tentative_result.valid(),
                "[Inlet] Lua function call failed, argument types possibly incorrect");
  return tentative_result;
}

/*!
 *****************************************************************************
 * \brief Templated function for extracting a concrete type from a sol function
 * result, used to allow for returning nonprimitive types, specifically, vectors
 *
 * \param [in] res The sol result of unknown concrete type
 * \tparam Ret The return type of the function
 *
 * \return The function's result
 *****************************************************************************
 */
template <typename Ret>
Ret extractResult(axom::sol::protected_function_result&& res)
{
  axom::sol::optional<Ret> option = res;
  SLIC_ERROR_IF(!option, "[Inlet] Lua function call failed, return types possibly incorrect");
  return option.value();
}

template <>
FunctionType::Void extractResult<FunctionType::Void>(axom::sol::protected_function_result&&)
{ }

/*!
 *****************************************************************************
 * \brief Creates a std::function given a Lua function and template parameters
 * corresponding to the function signature
 *
 * \param [in] func The sol object containing the lua function of unknown signature
 * \tparam Ret The return type of the function
 * \tparam Args... The argument types of the function
 *
 * \return A std::function that wraps the lua function
 * 
 * \note This is needed as a layer of indirection for bindArgType so it can
 * properly deduce the constructor call
 *****************************************************************************
 */
template <typename Ret, typename... Args>
std::function<Ret(typename detail::inlet_function_arg_type<Args>::type...)> buildStdFunction(
  axom::sol::protected_function&& func)
{
  // Generalized lambda capture needed to move into lambda
  return [func(std::move(func))](typename detail::inlet_function_arg_type<Args>::type... args) {
    return extractResult<Ret>(callWith(func, args...));
  };
}

/*!
 *****************************************************************************
 * \brief Adds argument types to a parameter pack based on the contents
 * of a std::vector of type tags
 *
 * \param [in] func The sol object containing the lua function of unknown signature
 * \param [in] arg_types The vector of argument types
 * 
 * \tparam I The number of arguments processed, or "stack size", used to mitigate
 * infinite compile-time recursion
 * \tparam Ret The function's return type
 * \tparam Args... The function's current arguments (already processed), remaining
 * arguments are in the arg_types vector
 *
 * \return A callable wrapper
 *****************************************************************************
 */
template <std::size_t I, typename Ret, typename... Args>
typename std::enable_if<(I > MAX_NUM_ARGS), FunctionVariant>::type bindArgType(
  axom::sol::protected_function&&,
  const std::vector<FunctionTag>&)
{
  SLIC_ERROR("[Inlet] Maximum number of function arguments exceeded: " << I);
  return {};
}

template <std::size_t I, typename Ret, typename... Args>
typename std::enable_if<I <= MAX_NUM_ARGS, FunctionVariant>::type bindArgType(
  axom::sol::protected_function&& func,
  const std::vector<FunctionTag>& arg_types)
{
  if(arg_types.size() == I)
  {
    return buildStdFunction<Ret, Args...>(std::move(func));
  }
  else
  {
    switch(arg_types[I])
    {
    case FunctionTag::Vector:
      return bindArgType<I + 1, Ret, Args..., FunctionType::Vector>(std::move(func), arg_types);
    case FunctionTag::Double:
      return bindArgType<I + 1, Ret, Args..., double>(std::move(func), arg_types);
    case FunctionTag::String:
      return bindArgType<I + 1, Ret, Args..., std::string>(std::move(func), arg_types);
    default:
      SLIC_ERROR("[Inlet] Unexpected function argument type");
    }
  }
  return {};  // Never reached but needed as errors do not imply control flow as with exceptions
}

/*!
 *****************************************************************************
 * \brief Performs a type-checked access to a Lua table
 *
 * \param [in]  proxy The axom::sol::proxy object to retrieve from
 * \param [out] val The value to write to, if it is of the correct type
 *
 * \return ReaderResult::Success if the object was of the correct type,
 * ReaderResult::WrongType otherwise
 *****************************************************************************
 */
template <typename Proxy, typename Value>
ReaderResult checkedGet(const Proxy& proxy, Value& val)
{
  axom::sol::optional<Value> option = proxy;
  if(option)
  {
    val = option.value();
    return ReaderResult::Success;
  }
  return ReaderResult::WrongType;
}

}  // end namespace detail

FunctionVariant LuaReader::getFunction(const std::string& id,
                                       const FunctionTag ret_type,
                                       const std::vector<FunctionTag>& arg_types)
{
  auto lua_func = getFunctionInternal(id);
  if(lua_func)
  {
    switch(ret_type)
    {
    case FunctionTag::Vector:
      return detail::bindArgType<0u, FunctionType::Vector>(std::move(lua_func), arg_types);
    case FunctionTag::Double:
      return detail::bindArgType<0u, double>(std::move(lua_func), arg_types);
    case FunctionTag::Void:
      return detail::bindArgType<0u, void>(std::move(lua_func), arg_types);
    case FunctionTag::String:
      return detail::bindArgType<0u, std::string>(std::move(lua_func), arg_types);
    default:
      SLIC_ERROR("[Inlet] Unexpected function return type");
    }
  }
  return {};  // Return an empty function to indicate that the function was not found
}

template <typename T>
ReaderResult LuaReader::getValue(const std::string& id, T& value)
{
  std::vector<std::string> tokens = axom::utilities::string::split(id, SCOPE_DELIMITER);

  if(tokens.size() == 1)
  {
    if((*m_lua)[tokens[0]].valid())
    {
      return detail::checkedGet((*m_lua)[tokens[0]], value);
    }
    return ReaderResult::NotFound;
  }

  axom::sol::table t;
  // Don't traverse through the last token as it doesn't contain a table
  if(traverseToTable(tokens.begin(), tokens.end() - 1, t))
  {
    if(t[tokens.back()].valid())
    {
      return detail::checkedGet(t[tokens.back()], value);
    }
  }

  return ReaderResult::NotFound;
}

std::vector<std::string> LuaReader::getAllNames()
{
  std::vector<std::string> result;
  detail::nameRetrievalHelper(m_preloaded_globals, m_lua->globals(), "", result);
  return result;
}

template <typename Key, typename Val>
ReaderResult LuaReader::getMap(const std::string& id,
                               std::unordered_map<Key, Val>& values,
                               axom::sol::type type)
{
  values.clear();
  std::vector<std::string> tokens = axom::utilities::string::split(id, SCOPE_DELIMITER);

  axom::sol::table t;
  if(tokens.empty() || !traverseToTable(tokens.begin(), tokens.end(), t))
  {
    return ReaderResult::NotFound;
  }

  // Allows for filtering out keys of incorrect type
  const auto is_correct_key_type = [](const axom::sol::type type) {
    bool is_number = type == axom::sol::type::number;
    // Arrays only
    if(std::is_same<Key, int>::value)
    {
      return is_number;
    }
    // Dictionaries can have both string-valued and numeric keys
    else
    {
      return is_number || (type == axom::sol::type::string);
    }
  };
  bool contains_other_type = false;
  for(const auto& entry : t)
  {
    // Gets only indexed items in the table.
    if(is_correct_key_type(entry.first.get_type()) && entry.second.get_type() == type)
    {
      values[detail::extractAs<Key>(entry.first)] = detail::extractAs<Val>(entry.second);
    }
    else
    {
      contains_other_type = true;
    }
  }
  return collectionRetrievalResult(contains_other_type, !values.empty());
}

template <typename T>
ReaderResult LuaReader::getIndicesInternal(const std::string& id, std::vector<T>& indices)
{
  std::vector<std::string> tokens = axom::utilities::string::split(id, SCOPE_DELIMITER);

  axom::sol::table t;

  if(tokens.empty() || !traverseToTable(tokens.begin(), tokens.end(), t))
  {
    return ReaderResult::NotFound;
  }

  indices.clear();

  // std::transform ends up being messier here
  for(const auto& entry : t)
  {
    indices.push_back(detail::extractAs<T>(entry.first));
  }
  return ReaderResult::Success;
}

axom::sol::protected_function LuaReader::getFunctionInternal(const std::string& id)
{
  std::vector<std::string> tokens = axom::utilities::string::split(id, SCOPE_DELIMITER);
  axom::sol::protected_function lua_func;

  if(tokens.size() == 1)
  {
    if((*m_lua)[tokens[0]].valid())
    {
      lua_func = (*m_lua)[tokens[0]];
      detail::checkedGet((*m_lua)[tokens[0]], lua_func);
    }
  }
  else
  {
    axom::sol::table t;
    // Don't traverse through the last token as it doesn't contain a table
    if(traverseToTable(tokens.begin(), tokens.end() - 1, t) && t[tokens.back()].valid())
    {
      detail::checkedGet(t[tokens.back()], lua_func);
    }
  }
  return lua_func;
}

}  // end namespace inlet
}  // end namespace axom
