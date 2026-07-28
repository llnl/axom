// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/klee/ShapeSet.hpp"

#include <string>
#include <istream>

namespace axom
{
namespace klee
{
/// Input deck formats supported by Klee.
enum class InputFormat
{
  YAML,
  Lua
};

/**
 * Read a ShapeSet from an input stream.
 *
 * \param stream the stream from which to read the ShapeSet
 * \note This overload reads YAML for backward compatibility.
 * \throws KleeError if parsing, schema verification, or semantic validation fails
 */
ShapeSet readShapeSet(std::istream &stream);

/**
 * Read a ShapeSet from an input stream.
 *
 * \param stream the stream from which to read the ShapeSet
 * \param format the input deck format to use
 * \throws KleeError if parsing, schema verification, or semantic validation fails,
 *         or if the requested input format is unsupported by this build
 */
ShapeSet readShapeSet(std::istream &stream, InputFormat format);

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
ShapeSet readShapeSet(const std::string &filePath);

/**
 * Read a ShapeSet from a specified file using an explicit input format.
 *
 * \param filePath the file from which to read the ShapeSet
 * \param format the input deck format to use, regardless of the file extension
 * \return the ShapeSet read from the file
 * \throws KleeError if parsing, schema verification, or semantic validation fails,
 *         or if the requested input format is unsupported by this build
 */
ShapeSet readShapeSet(const std::string &filePath, InputFormat format);

}  // namespace klee
}  // namespace axom
