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

/**
 * \brief Lua source to evaluate before the Klee input
 *
 * The chunk must return a table. Its keys must be ASCII Lua identifiers that are
 * neither keywords nor preloaded globals. Values retain their original Lua types,
 * including userdata such as Vector.new(1, 2). The table entries become mutable
 * deck globals and exported functions retain access to the chunk's locals and environment.
 *
 * The environment shares preloaded objects with the deck.
 */
struct LuaInitializationChunk
{
  /// Nonempty Lua source that returns the table of exported globals.
  std::string source;

  /// Diagnostic label. An empty label uses the default label shown below.
  std::string label {"<lua initialization>"};
};

/// Optional caller-provided initialization for a Lua input deck.
struct LuaInputOptions
{
  /// Chunk evaluated once per readShapeSet() call. Omit to read without initialization.
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
 * Read a ShapeSet from an input stream with optional Lua initialization.
 *
 * \param stream the stream from which to read the ShapeSet
 * \param format the input deck format to use
 * \param options optional initialization for a Lua input deck
 * \return the ShapeSet read from the stream
 * \throws KleeError if initialization, parsing, validation, or callback evaluation fails,
 *         or if the input format is unsupported by this build
 */
ShapeSet readShapeSet(std::istream& stream, InputFormat format, const LuaInputOptions& options);

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
 * Read a ShapeSet from a file with optional Lua initialization.
 *
 * \param filePath the file from which to read the ShapeSet
 * \param options optional extra initialization for a Lua input deck
 * \note The input format is inferred from the file extension; no extension means YAML.
 *       A populated initialization requires Lua input; empty options also accept YAML.
 * \return the ShapeSet read from the file
 * \throws KleeError if initialization, parsing, validation, or callback evaluation
 *         fails, or if the file extension or input format is unsupported
 */
ShapeSet readShapeSet(const std::string& filePath, const LuaInputOptions& options);

/**
 * Read a ShapeSet from a file using an explicit format and optional Lua initialization.
 *
 * \param filePath the file from which to read the ShapeSet
 * \param format the input file format to use, regardless of the file extension
 * \param options optional extra initialization for a Lua input deck
 * \return the ShapeSet read from the file
 * \throws KleeError if initialization, parsing, validation, or callback evaluation fails,
 *         or if the input format is unsupported by this build
 */
ShapeSet readShapeSet(const std::string& filePath, InputFormat format, const LuaInputOptions& options);

}  // namespace klee
}  // namespace axom
