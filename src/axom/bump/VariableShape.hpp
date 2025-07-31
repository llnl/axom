// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_UTILITIES_VARIABLE_SHAPE_HPP_
#define AXOM_BUMP_UTILITIES_VARIABLE_SHAPE_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/slic.hpp"
#include "axom/bump/views/Shapes.hpp"

#include <iostream>

namespace axom
{
namespace bump
{

/**
 * \brief Represent various shapes that consist of points.
 *
 * \tparam T The coordinate precision type.
 * \tparam NDIMS The spatial coordinate dimension.
 * \tparam N The max number of points that the shape can contain.
 */
template <typename T, int NDIMS, int N = 8>
class VariableShape
{
public:
  using PointType = axom::primal::Point<T, NDIMS>;

  /*!
   * \brief Return the shape id type.
   * \return The shape id type.
   */
  AXOM_HOST_DEVICE int id() const { return m_shapeId; }

  /*!
   * \brief Return the number of points in the shape.
   * \return The number of points in the shape.
   */
  AXOM_HOST_DEVICE axom::IndexType size() const { return m_points.size(); }

  /*!
   * \brief Add a point to the shape.
   * \param pt The point to add.
   */
  AXOM_HOST_DEVICE void push_back(const PointType &pt) { m_points.push_back(pt); }

  /*!
   * \brief Return the \a index'th point.
   * \param index The index of the point to return.
   * \return The desired point.
   */
  AXOM_HOST_DEVICE const PointType &operator[](axom::IndexType index) const
  {
    return m_points[index];
  }

  /*!
   * \brief Return unsigned shape volume.
   * \return Unsigned shape volume.
   */
  AXOM_HOST_DEVICE double volume() const
  {
    double retval = 0.;
    if(m_shapeId == axom::bump::views::Tet_ShapeID)
    {
      axom::primal::Tetrahedron<T, 3> tet(m_points[0], m_points[1], m_points[2], m_points[3]);
      retval = tet.volume();
    }
    else if(m_shapeId == axom::bump::views::Pyramid_ShapeID)
    {
      axom::primal::Tetrahedron<T, 3> tets[2];
      splitPyramid(tets);
      for(int i = 0; i < 2; i++)
      {
        retval += tets[i].volume();
      }
    }
    else if(m_shapeId == axom::bump::views::Wedge_ShapeID)
    {
      axom::primal::Tetrahedron<T, 3> tets[3];
      splitWedge(tets);
      for(int i = 0; i < 3; i++)
      {
        retval += tets[i].volume();
      }
    }
    else if(m_shapeId == axom::bump::views::Hex_ShapeID)
    {
      axom::primal::Hexahedron<T, 3> hex(m_points[0],
                                         m_points[1],
                                         m_points[2],
                                         m_points[3],
                                         m_points[4],
                                         m_points[5],
                                         m_points[6],
                                         m_points[7]);
      retval = hex.volume();
    }
    else
    {
      SLIC_ASSERT("Unsupported shape type");
    }
    return retval;
  }

  /*!
   * \brief Split the shape into tets.
   * \param[out] tets The output array of tets that make up the shape.
   */
  AXOM_HOST_DEVICE void splitPyramid(axom::primal::Tetrahedron<T, NDIMS> tets[2]) const
  {
    SLIC_ASSERT(m_shapeId == axom::bump::views::Pyramid_ShapeID);
    tets[0] = axom::primal::Tetrahedron<T, NDIMS>(m_points[0], m_points[1], m_points[3], m_points[4]);
    tets[1] = axom::primal::Tetrahedron<T, NDIMS>(m_points[1], m_points[2], m_points[3], m_points[4]);
  }

  /*!
   * \brief Split the shape into tets.
   * \param[out] tets The output array of tets that make up the shape.
   */
  AXOM_HOST_DEVICE void splitWedge(axom::primal::Tetrahedron<T, NDIMS> tets[3]) const
  {
    SLIC_ASSERT(m_shapeId == axom::bump::views::Wedge_ShapeID);
    tets[0] = axom::primal::Tetrahedron<T, NDIMS>(m_points[0], m_points[1], m_points[2], m_points[3]);
    tets[1] = axom::primal::Tetrahedron<T, NDIMS>(m_points[3], m_points[1], m_points[2], m_points[5]);
    tets[2] = axom::primal::Tetrahedron<T, NDIMS>(m_points[3], m_points[4], m_points[1], m_points[5]);
  }

  int m_shapeId;
  axom::StaticArray<PointType, N> m_points;
};

/// Printing method for VariableShape objects.
template <typename T, int NDIMS, int N = 8>
std::ostream &operator<<(std::ostream &os, const VariableShape<T, NDIMS, N> &obj)
{
  os << "{shapeId=" << obj.m_shapeId << ", points={";
  for(int i = 0; i < obj.m_points.size(); i++)
  {
    if(i > 0)
    {
      os << ", ";
    }
    os << obj.m_points[i];
  }
  os << "}}";
  return os;
}

}  // namespace bump
}  // namespace axom

#endif
