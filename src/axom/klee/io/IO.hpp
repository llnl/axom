// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_IO_HPP_
#define AXOM_KLEE_IO_HPP_

#include "axom/klee/ShapeSet.hpp"

#include <string>
#include <istream>

namespace axom
{
namespace klee
{
/**
 * Read a ShapeSet from an input stream.
 *
 * \param stream the stream from which to read the ShapeSet
 * \throws runtime_error if the input is invalid
 */
ShapeSet readShapeSet(std::istream &stream);

/**
 * Read a ShapeSet from a specified file
 *
 * \param filePath the file from which to read the ShapeSet
 * \return the ShapeSet read from the file
 * \throws runtime_error if the input is invalid
 */
ShapeSet readShapeSet(const std::string &filePath);

}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_IO_HPP_
