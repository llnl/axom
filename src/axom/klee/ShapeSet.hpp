// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_SHAPESET_HPP_
#define AXOM_KLEE_SHAPESET_HPP_

#include "axom/klee/Dimensions.hpp"
#include "axom/klee/Shape.hpp"

#include <string>
#include <vector>

namespace axom
{
namespace klee
{
/// A ShapeSet represents a document in the common shape format.
class ShapeSet
{
public:
  /**
   * Set the shapes in this set.
   *
   * \param shapes all the shapes in this set
   */
  void setShapes(std::vector<Shape> shapes);

  /**
   * Get the shapes in this set.
   *
   * \return the shapes in this set
   */
  std::vector<Shape> const &getShapes() const { return m_shapes; }

  /**
   * Set the file path from which this ShapeSet was created. This must be
   * set for resolvePath() to work.
   *
   * \param path the ShapeSet's path
   */
  void setPath(const std::string &path);

  /**
   * Get the path of the file from which this ShapeSet was created.
   *
   * \return the path of the file. Can be empty.
   */
  const std::string &getPath() const { return m_path; }

  /**
   * Sets the dimensions for all shapes in the ShapeSet.
   *
   * \param dimensions the dimension for all the shapes
   * \note This function must be called before calling \a getDimensions()
   */
  void setDimensions(Dimensions dimensions);

  /**
   * Returns the dimension of the ShapeSet.
   *
   * \pre Only valid after \a setDimensions() has been called on this instance
   * \sa setDimensions()
   */
  Dimensions getDimensions() const;

private:
  std::vector<Shape> m_shapes;
  std::string m_path;
  Dimensions m_dimensions {Dimensions::Unspecified};
};

}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_SHAPESET_HPP_
