// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/ShapeSet.hpp"

#include <utility>
#include <stdexcept>

#include "axom/core/utilities/FileUtilities.hpp"

namespace axom
{
namespace klee
{
void ShapeSet::setShapes(std::vector<Shape> shapes) { m_shapes = std::move(shapes); }

void ShapeSet::setPath(const std::string &path) { m_path = path; }

void ShapeSet::setDimensions(Dimensions dimensions)
{
  m_dimensions = dimensions;
  m_dimensionsHaveBeenSet = true;
}

Dimensions ShapeSet::getDimensions() const
{
  if(!m_dimensionsHaveBeenSet)
  {
    throw std::logic_error("Can only query the ShapeSet dimensions after calling setShapes()");
  }

  return m_dimensions;
}

}  // namespace klee
}  // namespace axom
