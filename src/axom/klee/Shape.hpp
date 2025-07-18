// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_SHAPE_HPP_
#define AXOM_KLEE_SHAPE_HPP_

#include "axom/klee/Geometry.hpp"

#include <string>
#include <vector>

namespace axom
{
namespace klee
{
/**
 * A Shape describes a region of space, a material it is made of, and a list
 * of materials which it does/doesn't replace.
 */
class Shape
{
public:
  /**
   * Create a new Shape.
   *
   * \param name the name of the shape
   * \param material the shape's material.
   * \param materialsReplaced the materials which can be replaced. If empty,
   * all materials can be replaced unless materialsNotReplaced is set.
   * \param materialsNotReplaced the materials which cannot be replaced. If
   * empty, all materials can be replaced unless materialsReplaced is set.
   * \param geometry the geometric properties of this shape
   * \throws std::logic_error if both materialsReplaced and materialsNotReplaced have entries.
   */
  Shape(std::string name,
        std::string material,
        std::vector<std::string> materialsReplaced,
        std::vector<std::string> materialsNotReplaced,
        Geometry geometry);

  /**
     * Get the name of this shape.
     * \return the shape's name
     */
  const std::string &getName() const { return m_name; }

  /**
   * Get the material this shape is made of.
   * \return the shape's material.
   */
  const std::string &getMaterial() const { return m_material; }

  /**
   * Check whether this shape can replace the given material (within the
   * space defined by this shape).
   *
   * \param material the material to check
   * \return whether this shape replaces the given material
   */
  bool replaces(const std::string &material) const;

  /**
   * Get the description fo the geometry for this shape.
   *
   * \return the shape's geometry
   */
  const Geometry &getGeometry() const { return m_geometry; }

  /**
   * Get the vector of materials that this shape can replace.
   * \return A reference to the material name vector.
   */
  const std::vector<std::string> &getMaterialsReplaced() const { return m_materialsReplaced; }

  /**
   * Get the vector of materials that this shape cannot replace.
   * \return A reference to the material name vector.
   */
  const std::vector<std::string> &getMaterialsNotReplaced() const { return m_materialsNotReplaced; }

private:
  std::string m_name;
  std::string m_material;
  std::vector<std::string> m_materialsReplaced;
  std::vector<std::string> m_materialsNotReplaced;
  Geometry m_geometry;
};

}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_SHAPE_HPP_
