// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/klee/ShapeSet.hpp"

#include <string>
#include <istream>
#include <unordered_map>
#include <variant>

namespace axom
{
namespace klee
{
/// Input file formats supported by Klee.
enum class InputFormat
{
  YAML,
  Lua
};

/// Runtime Lua chunk evaluated before deck parsing; exported bindings become initial Lua globals.
struct LuaBindingsChunk
{
  std::string source;
  std::string label {"<lua bindings>"};
};

/// Primitive value types that may be set as initial Lua globals.
using InputVariableValue = std::variant<bool, int, double, std::string>;

/// Variables to set as ordinary mutable globals before a Lua input deck is evaluated.
using InputVariables = std::unordered_map<std::string, InputVariableValue>;

/**
 * Read a ShapeSet from an input stream.
 *
 * \param stream the stream from which to read the ShapeSet
 * \note This overload reads YAML for backward compatibility.
 * \throws KleeError if parsing, schema verification, or semantic validation fails
 */
ShapeSet readShapeSet(std::istream& stream);

/**
 * Read a ShapeSet from an input stream.
 *
 * \param stream the stream from which to read the ShapeSet
 * \param format the input file format to use
 * \throws KleeError if parsing, schema verification, or semantic validation fails,
 *         or if the requested input format is unsupported by this build
 */
ShapeSet readShapeSet(std::istream& stream, InputFormat format);

/**
 * Read a ShapeSet from an input stream with caller-provided input variables.
 *
 * \param stream the stream from which to read the ShapeSet
 * \param format the input deck format to use
 * \param variables primitive values to set as initial mutable Lua globals
 * \note Input variables are supported only for Lua input decks.
 * \throws runtime_error if the input is invalid
 */
ShapeSet readShapeSet(std::istream &stream, InputFormat format, const InputVariables &variables);

/**
 * Read a ShapeSet from an input stream with caller-provided Lua bindings.
 *
 * \param stream the stream from which to read the ShapeSet
 * \param format the input deck format to use
 * \param bindings Lua chunk evaluated before deck parsing; must return a table
 *        of exported bindings, which become initial mutable Lua globals
 * \note Lua bindings are supported only for Lua input decks.
 * \throws runtime_error if the input is invalid
 */
ShapeSet readShapeSet(std::istream &stream, InputFormat format, const LuaBindingsChunk &bindings);

/**
 * Read a ShapeSet from an input stream with caller-provided input variables and Lua bindings.
 *
 * \param stream the stream from which to read the ShapeSet
 * \param format the input deck format to use
 * \param variables primitive values to set as initial mutable Lua globals
 * \param bindings Lua chunk evaluated before deck parsing; must return a table
 *        of exported bindings, which become initial mutable Lua globals
 * \note Input variables and Lua bindings are supported only for Lua input decks.
 * \throws runtime_error if the input is invalid
 */
ShapeSet readShapeSet(std::istream &stream,
                      InputFormat format,
                      const InputVariables &variables,
                      const LuaBindingsChunk &bindings);

/**
 * Read a ShapeSet from a specified file
 *
 * \param filePath the file from which to read the ShapeSet
 * \note The input format is inferred from the file extension. Files without
 * an extension are read as YAML for backward compatibility.
 * \return the ShapeSet read from the file
 * \throws KleeError if the extension is unsupported or if parsing,
 *         schema verification, or semantic validation fails
 */
ShapeSet readShapeSet(const std::string& filePath);

/**
 * Read a ShapeSet from a specified file using an explicit input format.
 *
 * \param filePath the file from which to read the ShapeSet
 * \param format the input file format to use, regardless of the file extension
 * \return the ShapeSet read from the file
 * \throws KleeError if parsing, schema verification, or semantic validation fails,
 *         or if the requested input format is unsupported by this build
 */
ShapeSet readShapeSet(const std::string& filePath, InputFormat format);

/**
 * Read a ShapeSet from a specified file with caller-provided input variables.
 *
 * \param filePath the file from which to read the ShapeSet
 * \param variables primitive values to set as initial mutable Lua globals
 * \note The input format is inferred from the file extension. Input variables
 *       are supported only for Lua input decks.
 * \return the ShapeSet read from the file
 * \throws runtime_error if the input is invalid
 */
ShapeSet readShapeSet(const std::string &filePath, const InputVariables &variables);

/**
 * Read a ShapeSet from a specified file with caller-provided Lua bindings.
 *
 * \param filePath the file from which to read the ShapeSet
 * \param bindings Lua chunk evaluated before deck parsing; must return a table
 *        of exported bindings, which become initial mutable Lua globals
 * \note The input format is inferred from the file extension. Lua bindings are
 *       supported only for Lua input decks.
 * \return the ShapeSet read from the file
 * \throws runtime_error if the input is invalid
 */
ShapeSet readShapeSet(const std::string &filePath, const LuaBindingsChunk &bindings);

/**
 * Read a ShapeSet from a specified file with caller-provided input variables and Lua bindings.
 *
 * \param filePath the file from which to read the ShapeSet
 * \param variables primitive values to set as initial mutable Lua globals
 * \param bindings Lua chunk evaluated before deck parsing; must return a table
 *        of exported bindings, which become initial mutable Lua globals
 * \note The input format is inferred from the file extension. Input variables
 *       and Lua bindings are supported only for Lua input decks.
 * \return the ShapeSet read from the file
 * \throws runtime_error if the input is invalid
 */
ShapeSet readShapeSet(const std::string &filePath,
                      const InputVariables &variables,
                      const LuaBindingsChunk &bindings);

}  // namespace klee
}  // namespace axom
