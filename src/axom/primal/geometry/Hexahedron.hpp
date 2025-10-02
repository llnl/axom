// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef HEXAHEDRON_HPP_
#define HEXAHEDRON_HPP_

#include "axom/core/StackArray.hpp"

#include "axom/core/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "axom/slic/interface/slic.hpp"

#include <ostream>  // for std::ostream

namespace axom
{
namespace primal
{
/*!
 * \class Hexahedron
 *
 * \brief Represents an hexahedral geometric shape defined by eight points.
 *
 * \accelerated
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of spatial dimensions
 *
 * There are eight vertices in the hexahedron, labelled P through W as the
 * constructor's arguments.  They are accessible using the square-brackets
 * operator, with P being index 0, Q index 1, through W as index 7.
 *
 * Here's a diagram showing
 * a cube with the labeled vertices.
 *
 * \verbatim
 *
 * S +---------+ R              +y
 *   |\        |\
 *   |  \      |  \              ^
 *   | W + --------+ V           |
 * P +---|-----+ Q |             |
 *   \   |     \   |             -----> +x
 *    \  |      \  |              \
 *     \ |       \ |               \
 *   T  +----------+ U              >
 *                                   +z
 * \endverbatim
 *
 */
template <typename T, int NDIMS = 3>
class Hexahedron
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using TetrahedronType = Tetrahedron<T, NDIMS>;
  using NumArrayType = axom::NumericArray<T, NDIMS>;

  enum
  {
    NUM_HEX_VERTS = 8,
    NUM_TRIANGULATE = 24
  };

public:
  /*!
   * \brief Default constructor. Creates a degenerate hexahedron.
   */
  Hexahedron() = default;

  /*!
   * \brief Creates an hexahedron from the 8 points p,q,r,s,t,u,v,w.
   * \param [in] Point corresponding to vertex p of the hexahedron.
   * \param [in] Point corresponding to vertex q of the hexahedron.
   * \param [in] Point corresponding to vertex r of the hexahedron.
   * \param [in] Point corresponding to vertex s of the hexahedron.
   * \param [in] Point corresponding to vertex t of the hexahedron.
   * \param [in] Point corresponding to vertex u of the hexahedron.
   * \param [in] Point corresponding to vertex v of the hexahedron.
   * \param [in] Point corresponding to vertex w of the hexahedron.
   *
   */
  AXOM_HOST_DEVICE
  Hexahedron(const PointType& p,
             const PointType& q,
             const PointType& r,
             const PointType& s,
             const PointType& t,
             const PointType& u,
             const PointType& v,
             const PointType& w)
  {
    m_points[0] = p;
    m_points[1] = q;
    m_points[2] = r;
    m_points[3] = s;
    m_points[4] = t;
    m_points[5] = u;
    m_points[6] = v;
    m_points[7] = w;
  }

  /*!
   * \brief Hexahedron constructor from an array of Points
   *
   * \param [in] pts An array containing at least 8 Points.
   *
   * \note It is the responsiblity of the caller to pass
   *       an array with at least 8 Points
   */
  AXOM_HOST_DEVICE
  explicit Hexahedron(const PointType* pts)
  {
    for(int i = 0; i < NUM_HEX_VERTS; i++)
    {
      m_points[i] = pts[i];
    }
  }

  /*!
   * \brief Hexahedron constructor from an Array of Points.
   *
   * \param [in] pts An ArrayView containing 8 Points.
   */
  AXOM_HOST_DEVICE
  explicit Hexahedron(const axom::ArrayView<PointType> pts)
  {
    SLIC_ASSERT(pts.size() == NUM_HEX_VERTS);

    for(int i = 0; i < NUM_HEX_VERTS; i++)
    {
      m_points[i] = pts[i];
    }
  }

  /*!
   * \brief Hexahedron constructor from an initializer list of Points
   *
   * \param [in] pts an initializer list containing 8 Points
   */
  AXOM_HOST_DEVICE
  explicit Hexahedron(std::initializer_list<PointType> pts)
  {
    SLIC_ASSERT(pts.size() == NUM_HEX_VERTS);

    int i = 0;
    for(const auto& pt : pts)
    {
      m_points[i] = pt;
      i++;
    }
  }

  /*!
   * \brief Return the number of vertices in a Hexahedron.
   *
   * \return The number of vertices in a Hexahedron.
   */
  AXOM_HOST_DEVICE static constexpr int numVertices() { return NUM_HEX_VERTS; }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, 3, 4, 5, 6, or 7
   */
  AXOM_HOST_DEVICE PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_HEX_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, 3, 4, 5, 6, or 7
   */
  AXOM_HOST_DEVICE const PointType& operator[](int idx) const
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_HEX_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Computes the average of the hexahedron's vertex
   *  positions
   *
   * \return A point at the mean of the hexahedron's vertices
   */
  AXOM_HOST_DEVICE
  PointType vertexMean() const
  {
    PointType sum;

    for(int i = 0; i < NUM_HEX_VERTS; ++i)
    {
      sum.array() += m_points[i].array();
    }
    sum.array() /= NUM_HEX_VERTS;

    return sum;
  }

  /**
   * \brief Finds the signed volume of the hexahedron
   *
   * Uses algorithm in
   * Grandy, J. Efficient computation of volume of hexahedral cells.
   * United States: N. p., 1997. Web. doi:10.2172/632793.
   *
   * \return The volume of the hexahedron
   *
   * \note Indices are different from the original paper
   */
  AXOM_HOST_DEVICE
  double signedVolume() const
  {
    constexpr double scale = 1. / 12.;
    return scale *
      (VectorType::scalar_triple_product(
         VectorType(m_points[6] - m_points[3]) + VectorType(m_points[5] - m_points[0]),
         VectorType(m_points[6] - m_points[4]),
         VectorType(m_points[7] - m_points[0])) +
       VectorType::scalar_triple_product(
         VectorType(m_points[5] - m_points[0]),
         VectorType(m_points[6] - m_points[4]) + VectorType(m_points[2] - m_points[0]),
         VectorType(m_points[6] - m_points[1])) +
       VectorType::scalar_triple_product(
         VectorType(m_points[6] - m_points[3]),
         VectorType(m_points[2] - m_points[0]),
         VectorType(m_points[6] - m_points[1]) + VectorType(m_points[7] - m_points[0])));
  }

  /*!
   * \brief Returns the absolute (unsigned) volume of the hexahedron
   * \sa signedVolume()
   */
  AXOM_HOST_DEVICE
  double volume() const { return axom::utilities::abs(signedVolume()); }

  /**
   * \brief Method to decompose a Hexahedron
   *        into 24 Tetrahedrons.
   *        Each Tetrahedron consists of 4 points:
   *          (1) The mean of all Hexahedron points
   *          (2,3) Two adjacent vertices on a Hexahedron face
   *          (4) The mean of the current Hexahedron face
   *
   * \param tets [out] The tetrahedrons
   *
   * \note Assumes tets is pre-allocated
   *
   * \tparam TetIndexable An indexable container of NUM_TRIANGULATE
   * TetrahedronType objects.
   */
  template <typename TetIndexable>
  AXOM_HOST_DEVICE void triangulate(TetIndexable& tets) const
  {
    // Hex center (hc)
    PointType hc = vertexMean();

    //Face means (fm)
    PointType fm1 = PointType::midpoint(PointType::midpoint(m_points[0], m_points[1]),
                                        PointType::midpoint(m_points[2], m_points[3]));

    PointType fm2 = PointType::midpoint(PointType::midpoint(m_points[0], m_points[1]),
                                        PointType::midpoint(m_points[4], m_points[5]));

    PointType fm3 = PointType::midpoint(PointType::midpoint(m_points[0], m_points[3]),
                                        PointType::midpoint(m_points[4], m_points[7]));

    PointType fm4 = PointType::midpoint(PointType::midpoint(m_points[1], m_points[2]),
                                        PointType::midpoint(m_points[5], m_points[6]));

    PointType fm5 = PointType::midpoint(PointType::midpoint(m_points[2], m_points[3]),
                                        PointType::midpoint(m_points[6], m_points[7]));

    PointType fm6 = PointType::midpoint(PointType::midpoint(m_points[4], m_points[5]),
                                        PointType::midpoint(m_points[6], m_points[7]));

    // Initialize tets
    tets[0] = TetrahedronType(hc, m_points[1], m_points[0], fm1);
    tets[1] = TetrahedronType(hc, m_points[0], m_points[3], fm1);
    tets[2] = TetrahedronType(hc, m_points[3], m_points[2], fm1);
    tets[3] = TetrahedronType(hc, m_points[2], m_points[1], fm1);

    tets[4] = TetrahedronType(hc, m_points[4], m_points[0], fm2);
    tets[5] = TetrahedronType(hc, m_points[0], m_points[1], fm2);
    tets[6] = TetrahedronType(hc, m_points[1], m_points[5], fm2);
    tets[7] = TetrahedronType(hc, m_points[5], m_points[4], fm2);

    tets[8] = TetrahedronType(hc, m_points[3], m_points[0], fm3);
    tets[9] = TetrahedronType(hc, m_points[0], m_points[4], fm3);
    tets[10] = TetrahedronType(hc, m_points[4], m_points[7], fm3);
    tets[11] = TetrahedronType(hc, m_points[7], m_points[3], fm3);

    tets[12] = TetrahedronType(hc, m_points[5], m_points[1], fm4);
    tets[13] = TetrahedronType(hc, m_points[1], m_points[2], fm4);
    tets[14] = TetrahedronType(hc, m_points[2], m_points[6], fm4);
    tets[15] = TetrahedronType(hc, m_points[6], m_points[5], fm4);

    tets[16] = TetrahedronType(hc, m_points[6], m_points[2], fm5);
    tets[17] = TetrahedronType(hc, m_points[2], m_points[3], fm5);
    tets[18] = TetrahedronType(hc, m_points[3], m_points[7], fm5);
    tets[19] = TetrahedronType(hc, m_points[7], m_points[6], fm5);

    tets[20] = TetrahedronType(hc, m_points[7], m_points[4], fm6);
    tets[21] = TetrahedronType(hc, m_points[4], m_points[5], fm6);
    tets[22] = TetrahedronType(hc, m_points[5], m_points[6], fm6);
    tets[23] = TetrahedronType(hc, m_points[6], m_points[7], fm6);
  }

  /*!
    * \brief Test if this Hexahedron is equal to another, within a tolerance
    *
    * \note This function can be an expensive operation
    */
  AXOM_HOST_DEVICE
  bool equals(const Hexahedron& other, double eps = 1.e-24) const
  {
    // Two hexes are equal if each vertex is closer than eps to a vertex of the other.
    int matched[NUM_HEX_VERTS];
    for(int i = 0; i < NUM_HEX_VERTS; ++i)
    {
      matched[i] = 0;
    }
    for(int ourvert = 0; ourvert < NUM_HEX_VERTS; ++ourvert)
    {
      for(int theirvert = 0; theirvert < NUM_HEX_VERTS; ++theirvert)
      {
        if(!matched[theirvert] && squared_distance(m_points[ourvert], other[theirvert]) < eps)
        {
          matched[theirvert] = 1;
        }
      }
    }
    int matchedcount = 0;
    for(int i = 0; i < NUM_HEX_VERTS; ++i)
    {
      matchedcount += matched[i];
    }

    return (matchedcount == NUM_HEX_VERTS);
  }

  /*!
   * \brief Simple formatted print of an hexahedron instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{" << m_points[0] << " " << m_points[1] << " " << m_points[2] << " " << m_points[3] << " "
       << m_points[4] << " " << m_points[5] << " " << m_points[6] << " " << m_points[7] << "}";

    return os;
  }

private:
  PointType m_points[NUM_HEX_VERTS] {PointType {},
                                     PointType {},
                                     PointType {},
                                     PointType {},
                                     PointType {},
                                     PointType {},
                                     PointType {},
                                     PointType {}};
};

//------------------------------------------------------------------------------
/// Free functions implementing Hexahedron's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Hexahedron<T, NDIMS>& hex)
{
  hex.print(os);
  return os;
}

} /* namespace primal */
} /* namespace axom */

#endif /* HEXAHEDRON_HPP_ */
