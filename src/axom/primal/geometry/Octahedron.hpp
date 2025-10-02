// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef OCTAHEDRON_HPP_
#define OCTAHEDRON_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "axom/slic/interface/slic.hpp"

#include <ostream>  // for std::ostream

namespace axom
{
namespace primal
{
/*!
 * \class Octahedron
 *
 * \brief Represents an octahedral geometric shape defined by six points.
 *
 * \accelerated
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of spatial dimensions
 *
 * There are six vertices in the octahedron, labelled P through U as the
 * constructor's arguments.  They are accessible using the square-brackets
 * operator, with P being index 0, Q index 1, through U as index 5.
 *
 * Imagine a regular octahedron with two parallel triangles, top and bottom---
 * the end-caps, as it were.  If you look "down", normal to the end-cap triangles,
 * you will see that the vertices of the top triangle protrude beyond the
 * edges of the bottom triangle (and vice versa).  Here's a diagram showing
 * just the end-cap triangles, omitting the other edges for clarity.
 *
 * \verbatim
 *
 *         P
 *         /\
 *    Q -------- U
 *      \      /
 *      /\    /\
 *    R --\  /-- T
 *         \/
 *         S
 *
 * \endverbatim
 *
 * Now imagine looking from the side, edge-on to the end-caps.  If you unroll
 * the "side-wall" triangles, you get a triangle strip.  Here's another diagram
 * showing just the triangle strip, with the same points labeled.  Points P
 * and Q are repeated so we can show all eight faces.
 *
 * \verbatim
 *
 *        Q --- S --- U --- Q
 *       / \   / \   / \   /
 *      /   \ /   \ /   \ /
 *     P --- R --- T --- P
 *
 * \endverbatim
 *
 */
template <typename T, int NDIMS = 3>
class Octahedron
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

  static constexpr int NUM_VERTS = 6;

public:
  /*!
   * \brief Default Octahedron constructor. Creates a degenerate octahedron.
   */
  AXOM_HOST_DEVICE Octahedron() { }

  /*!
   * \brief Creates an octahedron from the 6 points p,q,r,s,t,u.
   * \param [in] Point corresponding to vertex p of the octahedron.
   * \param [in] Point corresponding to vertex q of the octahedron.
   * \param [in] Point corresponding to vertex r of the octahedron.
   * \param [in] Point corresponding to vertex s of the octahedron.
   * \param [in] Point corresponding to vertex t of the octahedron.
   * \param [in] Point corresponding to vertex u of the octahedron.
   *
   * p is opposite s, q is opposite t, r is opposite u.
   */
  AXOM_HOST_DEVICE
  Octahedron(const PointType& p,
             const PointType& q,
             const PointType& r,
             const PointType& s,
             const PointType& t,
             const PointType& u)
  {
    m_points[0] = p;
    m_points[1] = q;
    m_points[2] = r;
    m_points[3] = s;
    m_points[4] = t;
    m_points[5] = u;
  }

  /*!
   * \brief Octahedron constructor from an array of Points
   *
   * \param [in] pts An array containing at least 6 Points.
   *
   * \note It is the responsiblity of the caller to pass
   *       an array with at least 6 Points
   */
  AXOM_HOST_DEVICE
  explicit Octahedron(const PointType* pts)
  {
    for(int i = 0; i < NUM_VERTS; i++)
    {
      m_points[i] = pts[i];
    }
  }

  /*!
   * \brief Octahedron constructor from an Array of Points
   *
   * \param [in] pts An ArrayView containing at 6 Points.
   */
  AXOM_HOST_DEVICE
  explicit Octahedron(const axom::ArrayView<PointType> pts)
  {
    SLIC_ASSERT(pts.size() == NUM_VERTS);

    for(int i = 0; i < NUM_VERTS; i++)
    {
      m_points[i] = pts[i];
    }
  }

  /*!
   * \brief Octahedron constructor from an initializer list of Points
   *
   * \param [in] pts an initializer list containing 6 Points
   */
  AXOM_HOST_DEVICE
  explicit Octahedron(std::initializer_list<PointType> pts)
  {
    SLIC_ASSERT(pts.size() == NUM_VERTS);

    int i = 0;
    for(const auto& pt : pts)
    {
      m_points[i] = pt;
      i++;
    }
  }

  /*!
   * \brief Return the number of vertices in an Octahedron.
   *
   * \return The number of vertices in an Octahedron.
   */
  AXOM_HOST_DEVICE static constexpr int numVertices() { return NUM_VERTS; }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, 3, 4, or 5
   */
  AXOM_HOST_DEVICE PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, 3, 4, or 5
   */
  AXOM_HOST_DEVICE const PointType& operator[](int idx) const
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_VERTS);
    return m_points[idx];
  }

  /*!
    * \brief Test if this Octahedron is equal to another, within a tolerance
    *
    * \note This function can be an expensive operation
    */
  AXOM_HOST_DEVICE
  bool equals(const Octahedron& other, double eps = 1.e-24) const
  {
    // Two octs are equal if each vertex is closer than eps to a vertex of the other.
    int matched[NUM_VERTS];
    for(int i = 0; i < NUM_VERTS; ++i)
    {
      matched[i] = 0;
    }
    for(int ourvert = 0; ourvert < NUM_VERTS; ++ourvert)
    {
      for(int theirvert = 0; theirvert < NUM_VERTS; ++theirvert)
      {
        if(!matched[theirvert] && squared_distance(m_points[ourvert], other[theirvert]) < eps)
        {
          matched[theirvert] = 1;
        }
      }
    }
    int matchedcount = 0;
    for(int i = 0; i < NUM_VERTS; ++i)
    {
      matchedcount += matched[i];
    }

    return (matchedcount == NUM_VERTS);
  }

  /*!
   * \brief Simple formatted print of an octahedron instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{" << m_points[0] << " " << m_points[1] << " " << m_points[2] << " " << m_points[3]
       << " " << m_points[4] << " " << m_points[5] << "}";

    return os;
  }

private:
  PointType m_points[NUM_VERTS];
};

//------------------------------------------------------------------------------
/// Free functions implementing Octahedron's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Octahedron<T, NDIMS>& oct)
{
  oct.print(os);
  return os;
}

} /* namespace primal */
} /* namespace axom */

#endif /* OCTAHEDRON_HPP_ */
