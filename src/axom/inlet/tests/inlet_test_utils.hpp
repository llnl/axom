// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef INLET_TEST_UTILS_HPP
#define INLET_TEST_UTILS_HPP

#include <vector>
#include <string>

#include "gtest/gtest.h"

#include "axom/inlet/YAMLReader.hpp"
#include "axom/inlet/JSONReader.hpp"

#ifdef AXOM_USE_SOL
  #include "axom/inlet/LuaReader.hpp"
#endif

namespace axom
{
namespace inlet
{
namespace detail
{
/*!
 *******************************************************************************
 * \class LuaTranslator
 * \brief A converter class that translates Lua to YAML or JSON
 * \note This class is not a fully functional Lua parser and should not be treated
 * as such.  It is designed only to parse the subset of Lua that maps to valid
 * YAML/JSON.
 *******************************************************************************
 */
class LuaTranslator
{
public:
  /*!
   *****************************************************************************
   * \brief Converts a Lua string to YAML
   * \param [in] luaString The string to convert
   * \note This function does not check for syntactic validity.  It is the
   * responsibility of the caller to pass a valid Lua string.
   *****************************************************************************
   */
  static std::string convertYAML(const std::string& luaString);

  /*!
   *****************************************************************************
   * \brief Converts a Lua string to JSON
   * \param [in] luaString The string to convert
   * \note This function does not check for syntactic validity.  It is the
   * responsibility of the caller to pass a valid Lua string.
   *****************************************************************************
   */
  static std::string convertJSON(const std::string& luaString);

private:
  /*!
   *****************************************************************************
   * \brief Adds a token to a vector of strings, separating any punctuation
   * that may precede or follow the token.
   * \param [in] token The token to add
   * \param [inout] tokens The list of tokens to append to
   *****************************************************************************
   */
  static void add_token(std::string&& token, std::vector<std::string>& tokens);

  /*!
   *****************************************************************************
   * \brief Splits a Lua string into tokens
   * \param [in] text The Lua code to lex/tokenize
   *****************************************************************************
   */
  static std::vector<std::string> tokenize(const std::string& text);
};

/*!
 *******************************************************************************
 * \brief Converts a Lua string to a string accepted by the reader of specified
 * type
 * \param [in] luaString The Lua code representing the contents of an input file
 * \tparam InletReader The inlet::Reader implementation intended to receive
 * the result of this function
 *******************************************************************************
 */
template <typename InletReader>
inline std::string fromLuaTo(const std::string& luaString)
{
  return luaString;
}
/// \overload
template <>
inline std::string fromLuaTo<axom::inlet::YAMLReader>(const std::string& luaString)
{
  return LuaTranslator::convertYAML(luaString);
}
/// \overload
template <>
inline std::string fromLuaTo<axom::inlet::JSONReader>(const std::string& luaString)
{
  return LuaTranslator::convertJSON(luaString);
}

#ifdef AXOM_USE_SOL
using ReaderTypes =
  ::testing::Types<axom::inlet::LuaReader, axom::inlet::YAMLReader, axom::inlet::JSONReader>;
#else
using ReaderTypes = ::testing::Types<axom::inlet::YAMLReader, axom::inlet::JSONReader>;
#endif

}  // namespace detail

}  // namespace inlet

}  // namespace axom

#endif
