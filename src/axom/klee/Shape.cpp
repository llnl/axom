// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Shape.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace axom
{
namespace klee
{
namespace
{
/**
 * Check whether a given container contains a value
 *
 * \tparam Container the type of the container
 * \param container the container to check for the value
 * \param value the value to look for
 * \return whether the container contains the value
 */
template <typename Container>
bool contains(const Container &container, const typename Container::value_type &value)
{
  auto endIter = std::end(container);
  return std::find(std::begin(container), endIter, value) != endIter;
}
}  // unnamed namespace

Shape::Shape(std::string name,
             std::string material,
             std::vector<std::string> materialsReplaced,
             std::vector<std::string> materialsNotReplaced,
             Geometry geometry)
  : m_name(std::move(name))
  , m_material(std::move(material))
  , m_materialsReplaced(std::move(materialsReplaced))
  , m_materialsNotReplaced(std::move(materialsNotReplaced))
  , m_geometry(std::move(geometry))
{
  if(!m_materialsNotReplaced.empty() && !m_materialsReplaced.empty())
  {
    throw std::logic_error(
      "Can't set both the list of materials to replace "
      "and materials to not replace");
  }
}

bool Shape::replaces(const std::string &material) const
{
  if(!m_materialsReplaced.empty())
  {
    return contains(m_materialsReplaced, material);
  }
  else if(!m_materialsNotReplaced.empty())
  {
    return !contains(m_materialsNotReplaced, material);
  }
  return true;
}

}  // namespace klee
}  // namespace axom
