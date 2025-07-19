// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/FileUtilities.hpp"

#include "axom/klee/ShapeSet.hpp"

#include <utility>
#include <stdexcept>

namespace axom
{
namespace klee
{
void ShapeSet::setShapes(std::vector<Shape> shapes) { m_shapes = std::move(shapes); }

void ShapeSet::setPath(const std::string &path) { m_path = path; }

void ShapeSet::setDimensions(Dimensions dimensions)
{
  if(dimensions != Dimensions::Two && dimensions != Dimensions::Three)
  {
    throw std::logic_error(
      "Invalid ShapeSet dimensions. Must be Dimensions::Two or Dimensions::Three");
  }
  m_dimensions = dimensions;
}

Dimensions ShapeSet::getDimensions() const
{
  if(m_dimensions == Dimensions::Unspecified)
  {
    throw std::logic_error("ShapeSet dimensions have not been set yet");
  }

  return m_dimensions;
}

}  // namespace klee
}  // namespace axom
