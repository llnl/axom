// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file LuaReader.cpp
 *
 * \brief This file contains the class implementation of the LuaReader.
 *******************************************************************************
 */

#include <array>
#include <fstream>
#include <memory>
#include <stdexcept>

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

bool extractVariantValue(const axom::sol::object& obj, VariantValue& value)
{
  switch(obj.get_type())
  {
  case axom::sol::type::boolean:
    value = obj.as<bool>();
    return true;
  case axom::sol::type::number:
  {
    const double double_value = obj.as<double>();
    const int int_value = obj.as<int>();
    value = (static_cast<double>(int_value) == double_value) ? VariantValue {int_value}
                                                             : VariantValue {double_value};
    return true;
  }
  case axom::sol::type::string:
    value = obj.as<std::string>();
    return true;
  default:
    return false;
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

ReaderResult LuaReader::getVariantMap(const std::string& id,
                                      std::unordered_map<int, VariantValue>& values)
{
  return getVariantMapInternal(id, values);
}

ReaderResult LuaReader::getVariantMap(const std::string& id,
                                      std::unordered_map<VariantKey, VariantValue>& values)
{
  return getVariantMapInternal(id, values);
}

template <typename Iter>
bool LuaReader::traverseToTable(Iter begin, Iter end, axom::sol::table& table)
{
  // Nothing to traverse
  if(begin == end)
  {
    return true;
  }

  axom::sol::object object = (*m_lua)[*begin];
  if(!object.valid() || object.get_type() != axom::sol::type::table)
  {
    return false;
  }

  // Use the first one to index into the global lua state
  table = object.as<axom::sol::table>();
  ++begin;

  // Then use the remaining keys to walk down to the requested table
  for(auto curr = begin; curr != end; ++curr)
  {
    auto key = *curr;
    bool is_int = conduit::utils::string_is_integer(key);
    axom::sol::object child;
    if(is_int)
    {
      const int key_as_int = conduit::utils::string_to_value<int>(key);
      if(table[key_as_int].valid())
      {
        child = table[key_as_int];
      }
    }
    if(!child.valid() && table[key].valid())
    {
      child = table[key];
    }
    if(!child.valid())
    {
      return false;
    }

    if(child.get_type() != axom::sol::type::table)
    {
      return false;
    }
    table = child.as<axom::sol::table>();
  }
  return true;
}

axom::sol::object LuaReader::getObject(const std::string& id)
{
  const auto tokens = axom::utilities::string::split(id, SCOPE_DELIMITER);
  if(tokens.empty())
  {
    return {};
  }

  if(tokens.size() == 1)
  {
    return (*m_lua)[tokens.front()];
  }

  axom::sol::table parent;
  if(!traverseToTable(tokens.begin(), tokens.end() - 1, parent))
  {
    return {};
  }

  const auto& key = tokens.back();
  const bool is_int = conduit::utils::string_is_integer(key);
  if(is_int)
  {
    const int key_as_int = conduit::utils::string_to_value<int>(key);
    if(parent[key_as_int].valid())
    {
      return parent[key_as_int];
    }
  }
  return parent[key];
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
 * \param [in] args Arguments forwarded to the Lua function
 * \tparam Args The argument types of the function
 *
 * \return A checkable version of the function's result
 *****************************************************************************
 */
template <typename... Args>
axom::sol::protected_function_result callWith(const axom::sol::protected_function& func,
                                              Args&&... args)
{
  // Lua functions are exposed to clients as std::functions that can be invoked
  // after schema verification. Use a catchable failure here so those clients can
  // add context; SLIC errors may abort or only log and continue.
  auto tentative_result = func(std::forward<Args>(args)...);
  if(!tentative_result.valid())
  {
    axom::sol::error err = tentative_result;
    throw std::runtime_error(fmt::format("[Inlet] Lua function call failed: {0}", err.what()));
  }
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
  if(!option)
  {
    // A failed result conversion is a runtime input error for this function
    // call. Throwing avoids dereferencing an empty optional after a SLIC log.
    throw std::runtime_error("[Inlet] Lua function call failed, return types possibly incorrect");
  }
  return option.value();
}

template <>
FunctionType::Void extractResult<FunctionType::Void>(axom::sol::protected_function_result&&)
{ }

template <>
FunctionType::Vector extractResult<FunctionType::Vector>(axom::sol::protected_function_result&& res)
{
  // Keep Vector.new(...) returns supported, but also accept raw numeric Lua
  // tables so input decks can write idiomatic vector callbacks such as
  // function() return {1.0, 2.0, 3.0} end.
  axom::sol::optional<FunctionType::Vector> vector_option = res;
  if(vector_option)
  {
    return vector_option.value();
  }

  axom::sol::optional<axom::sol::table> table_option = res;
  if(table_option)
  {
    axom::sol::table table = table_option.value();
    std::array<double, 3> values {{0., 0., 0.}};
    std::array<bool, 3> seen {{false, false, false}};
    int count = 0;

    for(const auto& entry : table)
    {
      if(entry.first.get_type() != axom::sol::type::number)
      {
        throw std::runtime_error(
          "[Inlet] Lua vector function return must only contain numeric indices");
      }

      const double numeric_index = entry.first.as<double>();
      const int index = entry.first.as<int>();
      if(static_cast<double>(index) != numeric_index || index < 1 || index > 3)
      {
        throw std::runtime_error(
          "[Inlet] Lua vector function return indices must be integers between 1 and 3");
      }
      if(entry.second.get_type() != axom::sol::type::number)
      {
        throw std::runtime_error(
          "[Inlet] Lua vector function return components must be numeric");
      }

      values[index - 1] = entry.second.as<double>();
      seen[index - 1] = true;
      ++count;
    }

    if(count < 1 || count > 3)
    {
      throw std::runtime_error(fmt::format(
        "[Inlet] Lua vector function returned a table with {0} entries; "
        "expected 1 to 3 numeric entries",
        count));
    }
    for(int i = 0; i < count; ++i)
    {
      if(!seen[i])
      {
        throw std::runtime_error(
          "[Inlet] Lua vector function return indices must be contiguous starting at 1");
      }
    }

    return FunctionType::Vector {values.data(), count};
  }

  throw std::runtime_error("[Inlet] Lua function call failed, return types possibly incorrect");
}

/*!
 *****************************************************************************
 * \brief Creates a std::function given a Lua function and template parameters
 * corresponding to the function signature
 *
 * \param [in] func The sol object containing the lua function of unknown signature
 * \param [in] lua_state Shared ownership of the Lua state used by \a func
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
  axom::sol::protected_function&& func,
  std::shared_ptr<axom::sol::state> lua_state)
{
  // Keep the Lua state alive for the lifetime of callbacks returned to callers.
  return [lua_state(std::move(lua_state)),
          func(std::move(func))](typename detail::inlet_function_arg_type<Args>::type... args) {
    SLIC_ASSERT(lua_state);
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
 * \param [in] lua_state Shared ownership of the Lua state used by \a func
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
  const std::vector<FunctionTag>&,
  std::shared_ptr<axom::sol::state>)
{
  SLIC_ERROR("[Inlet] Maximum number of function arguments exceeded: " << I);
  return {};
}

template <std::size_t I, typename Ret, typename... Args>
typename std::enable_if<I <= MAX_NUM_ARGS, FunctionVariant>::type bindArgType(
  axom::sol::protected_function&& func,
  const std::vector<FunctionTag>& arg_types,
  std::shared_ptr<axom::sol::state> lua_state)
{
  if(arg_types.size() == I)
  {
    return buildStdFunction<Ret, Args...>(std::move(func), std::move(lua_state));
  }
  else
  {
    switch(arg_types[I])
    {
    case FunctionTag::Vector:
      return bindArgType<I + 1, Ret, Args..., FunctionType::Vector>(std::move(func),
                                                                    arg_types,
                                                                    std::move(lua_state));
    case FunctionTag::Double:
      return bindArgType<I + 1, Ret, Args..., double>(std::move(func),
                                                      arg_types,
                                                      std::move(lua_state));
    case FunctionTag::String:
      return bindArgType<I + 1, Ret, Args..., std::string>(std::move(func),
                                                           arg_types,
                                                           std::move(lua_state));
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
  if(proxy.template is<Value>())
  {
    val = proxy.template as<Value>();
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
    FunctionVariant function;
    switch(ret_type)
    {
    case FunctionTag::Vector:
      function =
        detail::bindArgType<0u, FunctionType::Vector>(std::move(lua_func), arg_types, m_lua);
      break;
    case FunctionTag::Double:
      function = detail::bindArgType<0u, double>(std::move(lua_func), arg_types, m_lua);
      break;
    case FunctionTag::Void:
      function = detail::bindArgType<0u, void>(std::move(lua_func), arg_types, m_lua);
      break;
    case FunctionTag::String:
      function = detail::bindArgType<0u, std::string>(std::move(lua_func), arg_types, m_lua);
      break;
    default:
      SLIC_ERROR("[Inlet] Unexpected function return type");
    }
    return function;
  }
  return {};  // Return an empty function to indicate that the function was not found
}

template <typename T>
ReaderResult LuaReader::getValue(const std::string& id, T& value)
{
  const auto object = getObject(id);
  if(!object.valid())
  {
    return ReaderResult::NotFound;
  }

  return detail::checkedGet(object, value);
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
  const auto object = getObject(id);
  if(!object.valid())
  {
    return ReaderResult::NotFound;
  }
  if(object.get_type() != axom::sol::type::table)
  {
    return ReaderResult::WrongType;
  }

  const auto table = object.as<axom::sol::table>();

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
  for(const auto& entry : table)
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

template <typename Key>
ReaderResult LuaReader::getVariantMapInternal(const std::string& id,
                                              std::unordered_map<Key, VariantValue>& values)
{
  values.clear();
  const auto object = getObject(id);
  if(!object.valid())
  {
    return ReaderResult::NotFound;
  }
  if(object.get_type() != axom::sol::type::table)
  {
    return ReaderResult::WrongType;
  }

  const auto table = object.as<axom::sol::table>();

  const auto is_correct_key_type = [](const axom::sol::type type) {
    const bool is_number = type == axom::sol::type::number;
    if(std::is_same<Key, int>::value)
    {
      return is_number;
    }
    else
    {
      return is_number || (type == axom::sol::type::string);
    }
  };

  bool contains_other_type = false;
  for(const auto& entry : table)
  {
    VariantValue value;
    if(is_correct_key_type(entry.first.get_type()) && detail::extractVariantValue(entry.second, value))
    {
      values[detail::extractAs<Key>(entry.first)] = value;
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
  indices.clear();
  const auto object = getObject(id);
  if(!object.valid())
  {
    return ReaderResult::NotFound;
  }
  if(object.get_type() != axom::sol::type::table)
  {
    return ReaderResult::WrongType;
  }

  const auto table = object.as<axom::sol::table>();
  // std::transform ends up being messier here
  for(const auto& entry : table)
  {
    indices.push_back(detail::extractAs<T>(entry.first));
  }
  return ReaderResult::Success;
}

axom::sol::protected_function LuaReader::getFunctionInternal(const std::string& id)
{
  axom::sol::protected_function lua_func;
  const auto object = getObject(id);
  if(object.valid())
  {
    detail::checkedGet(object, lua_func);
  }
  return lua_func;
}

}  // end namespace inlet
}  // end namespace axom
