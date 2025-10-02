// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_TETRAHEDRON_HPP_
#define AXOM_PRIMAL_TETRAHEDRON_HPP_

#include "axom/core.hpp"

#include "axom/primal/constants.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Sphere.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/fmt.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
/*!
 * \class Tetrahedron
 *
 * \brief Represents a tetrahedral geometric shape defined by four points.
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of spatial dimensions
 */
template <typename T, int NDIMS>
class Tetrahedron
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using SphereType = Sphere<T, NDIMS>;

  static constexpr int NUM_VERTS = 4;

public:
  /// \brief Default constructor. Creates a degenerate tetrahedron.
  AXOM_HOST_DEVICE Tetrahedron() : m_points {PointType(), PointType(), PointType(), PointType()} { }

  /*!
   * \brief Tetrahedron constructor from 4 points A,B,C,D.
   * \param [in] A Vertex A of the tetrahedron
   * \param [in] B Vertex B of the tetrahedron
   * \param [in] C Vertex C of the tetrahedron
   * \param [in] D Vertex D of the tetrahedron
   * 
   * \note The orientation of the tetrahedron is determined from the 
   * scalar triple product of the vectors from B, C and D to A, respectively,
   * i.e. \f$ dot(B-A, cross(C-A, D-A)) \f$.
   */
  AXOM_HOST_DEVICE
  Tetrahedron(const PointType& A, const PointType& B, const PointType& C, const PointType& D)
    : m_points {A, B, C, D}
  { }

  /*!
   * \brief Tetrahedron constructor from an array of Points
   *
   * \param [in] pts An array containing at least 4 Points.
   *
   * \note It is the responsiblity of the caller to pass
   *       an array with at least 4 Points
   */
  AXOM_HOST_DEVICE
  explicit Tetrahedron(const PointType* pts)
  {
    for(int i = 0; i < NUM_VERTS; i++)
    {
      m_points[i] = pts[i];
    }
  }

  /*!
   * \brief Tetrahedron constructor from an Array of Points.
   *
   * \param [in] pts An ArrayView containing 4 Points.
   */
  AXOM_HOST_DEVICE
  explicit Tetrahedron(const axom::ArrayView<PointType> pts)
  {
    SLIC_ASSERT(pts.size() == NUM_VERTS);

    for(int i = 0; i < NUM_VERTS; i++)
    {
      m_points[i] = pts[i];
    }
  }

  /*!
   * \brief Tetrahedron constructor from an initializer list of Points
   *
   * \param [in] pts an initializer list containing 4 Points
   */
  AXOM_HOST_DEVICE
  explicit Tetrahedron(std::initializer_list<PointType> pts)
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
   * \brief Return the number of vertices in a Tetrahedron.
   *
   * \return The number of vertices in a Tetrahedron.
   */
  AXOM_HOST_DEVICE static constexpr int numVertices() { return NUM_VERTS; }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, or 3
   */
  AXOM_HOST_DEVICE PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, or 3
   */
  AXOM_HOST_DEVICE const PointType& operator[](int idx) const
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Returns whether the tetrahedron is degenerate
   * \return true iff the tetrahedron is degenerate (has near zero volume)
   */
  AXOM_HOST_DEVICE
  bool degenerate(double eps = 1.0e-12) const
  {
    return axom::utilities::isNearlyEqual(ppedVolume(), 0.0, eps);
  }

  /*!
   * \brief Returns the barycentric coordinates of a point within a tetrahedron
   *
   * \param [in] p The point at which we want to compute barycentric coordinates
   * \param [in] skipNormalization Determines if the result should be normalized
   * by the volume of the tetrahedron. One might want to skip this if they were
   * only interested in the relative weights rather than their actual amounts
   * 
   * \post The barycentric coordinates sum to 1 when \a skipNormalization is false
   * Otherwise, the sum of coordinates will be proportional to volume of the tetrahedron
   * (Specifically, they should sum to the parallelpiped volume, which is 6x the volume).
   */
  AXOM_HOST_DEVICE
  Point<double, 4> physToBarycentric(const PointType& p, bool skipNormalization = false) const
  {
    Point<double, 4> bary;

    const PointType& A = m_points[0];
    const PointType& B = m_points[1];
    const PointType& C = m_points[2];
    const PointType& D = m_points[3];

    const auto pA = A - p;
    const auto pB = B - p;
    const auto pC = C - p;
    const auto pD = D - p;

    const double detA = VectorType::scalar_triple_product(pB, pC, pD);
    const double detB = -VectorType::scalar_triple_product(pC, pD, pA);
    const double detC = VectorType::scalar_triple_product(pD, pA, pB);
    const double detD = -VectorType::scalar_triple_product(pA, pB, pC);

    if(!skipNormalization)
    {
      const double vol = VectorType::scalar_triple_product(B - A, C - A, D - A);
      const double EPS = vol >= 0 ? primal::PRIMAL_TINY : -primal::PRIMAL_TINY;

      // Compute one over denominator; offset by a tiny amount to avoid division by zero
      const double ood = 1. / (vol + EPS);

      bary[0] = detA * ood;
      bary[1] = detB * ood;
      bary[2] = detC * ood;
      bary[3] = detD * ood;

      // Replace the smallest entry with the difference of 1 from the sum of the others
      const int amin = axom::abs(bary.array()).argMin();
      bary[amin] = 0.;
      bary[amin] = 1. - bary.array().sum();
    }
    else
    {
      bary[0] = detA;
      bary[1] = detB;
      bary[2] = detC;
      bary[3] = detD;
    }

    return bary;
  }

  /*!
   * \brief Returns the physical coordinates of a barycentric point
   * \param [in] bary Barycentric coordinates relative to this tetrahedron
   */
  AXOM_HOST_DEVICE
  PointType baryToPhysical(const Point<double, 4>& bary) const
  {
    SLIC_CHECK_MSG(axom::utilities::isNearlyEqual(1., bary[0] + bary[1] + bary[2] + bary[3]),
                   "Barycentric coordinates must sum to (near) one.");

    PointType res;
    for(int i = 0; i < NUM_VERTS; ++i)
    {
      res.array() += bary[i] * m_points[i].array();
    }

    return res;
  }

  /*! 
   * \brief Returns whether point \a p is contained within the tetrahedron (within tolerance \a eps)
   */
  AXOM_HOST_DEVICE
  bool contains(const PointType& p, double eps = 1e-8) const
  {
    const auto bC = physToBarycentric(p);
    return (bC[0] >= (0.0 - eps)) && (bC[0] <= (1.0 + eps))  //
      && (bC[1] >= (0.0 - eps)) && (bC[1] <= (1.0 + eps))    //
      && (bC[2] >= (0.0 - eps)) && (bC[2] <= (1.0 + eps))    //
      && (bC[3] >= (0.0 - eps)) && (bC[3] <= (1.0 + eps));
  }

  /*!
   * \brief Simple formatted print of a tetrahedron instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{" << m_points[0] << " " << m_points[1] << " " << m_points[2] << " " << m_points[3] << "}";

    return os;
  }

  /*!
   * \brief Returns the signed volume of the tetrahedron
   * \sa volume()
   */
  AXOM_HOST_DEVICE
  double signedVolume() const
  {
    constexpr double scale = 1. / 6.;
    return scale * ppedVolume();
  }

  /*!
   * \brief Returns the absolute (unsigned) volume of the tetrahedron
   * \sa signedVolume()
   */
  AXOM_HOST_DEVICE
  double volume() const { return axom::utilities::abs(signedVolume()); }

  /*!
   * \brief Swaps the order of vertices if the signed volume of the
   *        tetrahedron is negative. Signed volume will become positive.
   * \sa signedVolume()
   */
  AXOM_HOST_DEVICE
  void checkAndFixOrientation()
  {
    if(signedVolume() < 0)
    {
      axom::utilities::swap<PointType>(m_points[1], m_points[2]);
    }
  }

  /**
   * \brief Returns the circumsphere (circumscribing sphere) of the tetrahedron
   *
   * Derived from formula on https://mathworld.wolfram.com/Circumsphere.html
   * but uses observation that radius can be derived from distance to a vertex
   * \note This function is only available for 3D tetrahedra in 3D
   */
  template <int TDIM = NDIMS>
  typename std::enable_if<TDIM == 3, SphereType>::type circumsphere() const
  {
    const PointType& p0 = m_points[0];
    const PointType& p1 = m_points[1];
    const PointType& p2 = m_points[2];
    const PointType& p3 = m_points[3];

    // It's useful to separate the x-, y- and z- components of
    // the difference vectors for p1, p2, and p3 from p0
    const VectorType vx {p1[0] - p0[0], p2[0] - p0[0], p3[0] - p0[0]};
    const VectorType vy {p1[1] - p0[1], p2[1] - p0[1], p3[1] - p0[1]};
    const VectorType vz {p1[2] - p0[2], p2[2] - p0[2], p3[2] - p0[2]};

    // We also need their squared norms
    const VectorType sq {vx[0] * vx[0] + vy[0] * vy[0] + vz[0] * vz[0],
                         vx[1] * vx[1] + vy[1] * vy[1] + vz[1] * vz[1],
                         vx[2] * vx[2] + vy[2] * vy[2] + vz[2] * vz[2]};

    // Compute one over denominator using a small offset to avoid division by zero
    const double a = VectorType::scalar_triple_product(vx, vy, vz);
    const double EPS = (a >= 0) ? primal::PRIMAL_TINY : -primal::PRIMAL_TINY;
    const double ood = 1. / (2 * a + EPS);

    // Compute offset from p0 to center
    const auto center_offset = ood *
      VectorType {VectorType::scalar_triple_product(sq, vy, vz),
                  VectorType::scalar_triple_product(sq, vz, vx),
                  VectorType::scalar_triple_product(sq, vx, vy)};

    return SphereType(p0 + center_offset, center_offset.norm());
  }

private:
  /*!
   * \brief Computes the signed volume of a parallelepiped defined by the
   * three edges of the tetrahedron incident to its first vertex
   *
   * \note The ppedVolume is a factor of 6 greater than that of the tetrahedron
   * \return The signed parallelepiped volume
   * \sa signedVolume(), volume()
   */
  AXOM_HOST_DEVICE
  double ppedVolume() const
  {
    return NDIMS != 3 ? 0.
                      : VectorType::scalar_triple_product(m_points[1] - m_points[0],
                                                          m_points[2] - m_points[0],
                                                          m_points[3] - m_points[0]);
  }

private:
  PointType m_points[4];
};

//------------------------------------------------------------------------------
/// Free functions implementing Tetrahedron's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Tetrahedron<T, NDIMS>& tet)
{
  tet.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Tetrahedron using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Tetrahedron<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_TETRAHEDRON_HPP_
