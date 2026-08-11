// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/klee/ShapeSet.hpp"

#include <string>
#include <istream>
#include <optional>

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

/// Lua initialization chunk evaluated before deck parsing in an isolated environment.
struct LuaInitializationChunk
{
  /// Lua source to evaluate before the input deck.
  std::string source;

  /// Source label used in diagnostics from this chunk.
  std::string label {"<lua initialization>"};
};

/// Optional caller-provided initialization for a Lua input deck.
struct LuaInputOptions
{
  /// Isolated chunk whose returned table entries become initial mutable Lua globals.
  std::optional<LuaInitializationChunk> initialization;
};

/**
 * Read a ShapeSet from an input stream.
 *
 * \param stream the stream from which to read the ShapeSet
 * \note This overload reads YAML for backward compatibility.
 * \return the ShapeSet read from the stream
 * \throws KleeError if parsing, schema verification, or semantic validation fails
 */
ShapeSet readShapeSet(std::istream& stream);

/**
 * Read a ShapeSet from an input stream.
 *
 * \param stream the stream from which to read the ShapeSet
 * \param format the input file format to use
 * \return the ShapeSet read from the stream
 * \throws KleeError if parsing, schema verification, or semantic validation fails,
 *         or if the requested input format is unsupported by this build
 */
ShapeSet readShapeSet(std::istream& stream, InputFormat format);

/**
 * Read a ShapeSet from an input stream with caller-provided Lua inputs.
 *
 * \param stream the stream from which to read the ShapeSet
 * \param format the input deck format to use
 * \param options optional initialization for a Lua input deck
 * \note Lua input options are supported only for Lua input decks.
 * \return the ShapeSet read from the stream
 * \throws KleeError if the input or Lua input options are invalid
 */
ShapeSet readShapeSet(std::istream &stream,
                      InputFormat format,
                      const LuaInputOptions &options);

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
 * Read a ShapeSet from a specified file with caller-provided Lua inputs.
 *
 * \param filePath the file from which to read the ShapeSet
 * \param options optional initialization for a Lua input deck
 * \note The input format is inferred from the file extension.
 *       Lua input options are supported only for Lua input decks.
 * \return the ShapeSet read from the file
 * \throws KleeError if the input or Lua input options are invalid
 */
ShapeSet readShapeSet(const std::string &filePath, const LuaInputOptions &options);

/**
 * Read a ShapeSet from a specified file using an explicit format and
 * caller-provided Lua inputs.
 *
 * \param filePath the file from which to read the ShapeSet
 * \param format the input file format to use, regardless of the file extension
 * \param options optional initialization for a Lua input deck
 * \note Lua input options are supported only for Lua input decks.
 * \return the ShapeSet read from the file
 * \throws KleeError if the input or Lua input options are invalid
 */
ShapeSet readShapeSet(const std::string &filePath,
                      InputFormat format,
                      const LuaInputOptions &options);

}  // namespace klee
}  // namespace axom
