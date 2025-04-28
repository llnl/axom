// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_ORIENTEDBOUNDINGBOX_HPP_
#define AXOM_PRIMAL_ORIENTEDBOUNDINGBOX_HPP_

#include <vector>

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/NumericLimits.hpp"

#include "axom/core/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class OrientedBoundingBox;

/// \name Forward Declared Overloaded Operators
///@{

/*!
 * \brief Equality comparison operator for oriented bounding boxes.
 *
 * Two oriented bounding boxes are equal when they have the same bounds and
 * axes.
 */
template <typename T, int NDIMS>
bool operator==(const OrientedBoundingBox<T, NDIMS>& lhs, const OrientedBoundingBox<T, NDIMS>& rhs);

/*!
 * \brief Inequality comparison operator for oriented bounding boxes.
 *
 * Two bounding boxes are unequal when they have different bounds
 */
template <typename T, int NDIMS>
bool operator!=(const OrientedBoundingBox<T, NDIMS>& lhs, const OrientedBoundingBox<T, NDIMS>& rhs);

/*!
 * \brief Overloaded output operator for oriented bounding boxes.
 */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const OrientedBoundingBox<T, NDIMS>& box);

// Forward declared methods
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS> compute_oriented_bounding_box(const Point<T, NDIMS>* pts, int n);

///@}

/*!
 * \class OrientedBoundingBox
 *
 * \brief OrientedBoundingBox represents an oriented bounding box defined by
 * its centroid, axes, and extents.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */

template <typename T, int NDIMS>
class OrientedBoundingBox
{
public:
  using CoordType = T;
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using OrientedBoxType = OrientedBoundingBox<T, NDIMS>;
  using BoxType = BoundingBox<T, NDIMS>;

public:
  /*!
   * \brief Constructor. Creates an invalid oriented bounding box by setting
   * the first coordinate of the first axis to the smallest number of type T.
   */
  OrientedBoundingBox();

  /*!
   * \brief Constructor. Creates an oriented bounding box containing a single
   * point. Extents are all zero, and axes are set to identity.
   */
  explicit OrientedBoundingBox(const PointType& pt);

  /*!
   * \brief Constructor. Creates an oriented bounding box from a collection of
   * points.
   *
   * Initialize the random number generator before using this constructor
   * (with srand()).  This constructor uses eigen_solve to find a good fit
   * around the passed-in points, which uses the random number generator.
   */
  OrientedBoundingBox(const PointType* pts, int n);

  /*!
   * \brief Constructor. Creates an oriented bounding box with given center,
   * axes, and extents. Normalizes axes and sets any negative extent to 0.
   * \param [in] c center of OBB
   * \param [in] u axes of OBB
   * \param [in] e extents of OBB
   * \note Axes are made orthonormal, but if they don't span, some will be
   * set to 0.
   */
  OrientedBoundingBox(const PointType& c, const VectorType (&u)[NDIMS], const VectorType& e);

  /*!
   * \brief Copy Constructor.
   * \param [in] other The oriented bounding box to copy
   */
  OrientedBoundingBox(const OrientedBoundingBox& other);

  /*!
   * \brief Destructor.
   */
  ~OrientedBoundingBox() { }

  /*!
   * \brief Resets the bounds to those of the default constructor.
   * \note This invalidates the bounding box (i.e. isValid() will be false).
   */
  void clear();

  /*!
   * \brief Returns the centroid (midpoint) of the bounding box.
   * \return Point at the bounding box centroid.
   */
  const PointType& getCentroid() const { return m_c; }

  /*!
   * \brief Return the axes.
   * \return Pointer to the axes.
   */
  const VectorType* getAxes() const { return m_u; };

  /*!
   * \brief Returns the extents of the oriented bounding box.
   * \return Extents in a vector.
   */
  const VectorType& getExtents() const { return m_e; }

  /*!
   * \brief Returns the vertices of the oriented bounding box.
   * \return vertices in a std::vector<VectorType>.
   */
  std::vector<PointType> vertices() const;

  /*!
   * \brief Returns the dimension of the ambient space for this bounding box.
   * \return d the dimension of this bounding box instance.
   * \post d >= 1.
   */
  int dimension() const { return NDIMS; };

  /*!
   * \brief Expands the box so it contains the passed in point.
   * \param pt [in] point to add in
   * \note Simply expands the extents of this box so that it contains pt,
   * this loses optimality of the box. When constructing a box from a
   * collection of points, the constructor which takes an array of points
   * should be preferred.
   * \post this->contains(pt) == true
   */
  void addPoint(PointType pt);

  /*!
   * \brief Expands the box so it contains the passed in box.
   * \param obb [in] OBB to add in
   * \note If obb is invalid, this is a no-op
   * \note Expands the extents so this box contains obb. This loses optimality of the box
   * for a better fit use merge_boxes in primal/compute_bounding_box.
   * \post this->contains(obb) == true
   */
  void addBox(OrientedBoxType obb);

  /*!
   * \brief Expands the extents by the given amount.
   * \param [in] expansionAmount an absolute amount to expand
   * \note This function checks to ensure the bounding box is valid afterwards.
   * \return A reference to the bounding box after it has been expanded.
   */
  OrientedBoundingBox& expand(T expansionAmount);

  /*!
   * \brief Returns the product of the extents.
   * \return Product of the extents.
   */
  T volume() const;

  /*!
   * \brief Scales the bounding box about its center by a given amount.
   * \param [in] scaleFactor the multiplicative factor by which to scale
   * \note Checks to ensure that the bounding box is valid after inflation.
   * \note If scaleFactor is less than 1, the bounding box will shrink.
   * \note If scaleFactor is 0, the bounding box will shrink to its midpoint
   *  towards the center, and we fix the bounds after shrinking
   * \note if scaleFactor is less than 0, we use -scaleFactor.
   *  \return A reference to the bounding box after it has been scaled.
   */
  OrientedBoundingBox& scale(double scaleFactor);

  /*!
   * \brief Shifts the bounding box by a fixed displacement.
   * \param [in] displacement the amount with which to move the bounding box
   * \note This function checks to ensure the bounding box is valid afterwards.
   * \return A reference to the bounding box after it has been shifted.
   */
  OrientedBoundingBox& shift(const VectorType& displacement);

  /*!
   * \brief Overloaded assignment operator.
   * \param [in] rhs oriented bounding box instance on the right-hand side
   * \return
   */
  OrientedBoundingBox& operator=(const OrientedBoundingBox& rhs);

  /*!
   * \brief Checks whether the box contains the point.
   * \param [in] otherPt the point that we are checking
   * \return status true if point inside the box, else false.
   *
   * \note This function assumes all intervals are closed
   * (i.e. contain their boundaries).  We may need to deal with open
   * and half open boundaries in the future.
   */
  template <typename OtherType>
  bool contains(const Point<OtherType, NDIMS>& otherPt, double EPS = 1E-8) const;

  /*!
   * \brief Checks whether the box fully contains another bounding box.
   * \param [in] otherOBB the bounding box that we are checking
   * \return status true if bb is inside the box, else false.
   * \note We are allowing the other bounding box to have a different coordinate
   *  type. This should work as long as the two Ts are comparable with
   *  operator<().
   */
  template <typename OtherType>
  bool contains(const OrientedBoundingBox<OtherType, NDIMS>& otherOBB, double EPS = 1E-8) const;

  /*!
   * \brief Returns the same point but in the local coordinates of this OBB.
   * \param [in] pt the pt in world coords we want to convert to local coords
   * \return Vector which is just pt in local coordinates, i.e. is the
   *  projection of the difference between pt and m_c onto the axes m_u.
   */
  PointType toLocal(const PointType& pt) const;

  /*!
   * \brief Returns the furthest point in this OBB from pt.
   * \param [in] pt the pt we want to find the furthest point from
   * \note If there are multiple furthest points, any one may be returned.
   * \return furthest point in this OBB from pt
   */
  PointType furthestPoint(const PointType& pt) const;

  /*!
   * \brief Checks if this bounding box is valid.
   */
  bool isValid() const;

  /*!
   * \brief Subdivides this oriented bounding box instance into two sub-boxes by
   *  splitting along the longest dimension.
   * \param [in,out] right the right sub-box.
   * \param [in,out] left  the left sub-box.
   */
  void bisect(OrientedBoxType& right, OrientedBoxType& left) const;

  /*!
   * \brief Simple formatted print of a bounding box instance.
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const;

  /// @}

  /*!
   * \brief Returns the list of vertices in an OBB instance.
   * \param  [in] bb the OBB whose vertices will be returned
   * \param [out] pnts the list of points
   * \post pnts.size() == (1 << NDIMS)
   */
  static void getPoints(const OrientedBoundingBox<T, NDIMS>& obb, std::vector<Point<T, NDIMS>>& pnts)
  {
    pnts = obb.vertices();
  }

private:
  /// \name Internal Helper Routines
  /// @{

  /*!
   * \brief Recursively enumerates the vertices of this box, storing them in l.
   */
  static void vertex_enum(std::vector<Point<T, NDIMS>>& l,
                          int i,
                          const Vector<T, NDIMS> u[NDIMS],
                          const Vector<T, NDIMS>& e,
                          Vector<T, NDIMS> curr)
  {
    if(i == NDIMS)  // base case
    {
      l.push_back(Point<T, NDIMS>(curr.array()));
    }
    else
    {
      curr += e[i] * u[i];
      vertex_enum(l, i + 1, u, e, curr);
      curr -= (static_cast<T>(2.) * e[i]) * u[i];
      vertex_enum(l, i + 1, u, e, curr);
    }
  }

  /// @}

  /*!
   * \brief Iterates through the axes, making each orthogonal to the previous
   * ones and then normalizing it. If after a given axis' components in the
   * previous axes' direction is subtracted off it is nearly 0, it is set to 0.
   * This is the case when a box is "flat" in one or more dimensions. If any
   * extent is negative, it is set to zero.
   */
  void checkAndFix();

private:
  PointType m_c;          // centroid
  VectorType m_u[NDIMS];  // axes (all unit length)
  VectorType m_e;         // extents (each coord >= 0 or it's invalid)
};

} /* namespace primal */

} /* namespace axom */

//------------------------------------------------------------------------------
//  OrientedBoundingBox implementation
//------------------------------------------------------------------------------

namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS>::OrientedBoundingBox()
{
  this->clear();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS>::OrientedBoundingBox(const Point<T, NDIMS>& pt)
{
  this->m_c = Point<T, NDIMS>(pt);

  for(int i = 0; i < NDIMS; i++)
  {
    // initialize ith axis to ith standard basis vector
    this->m_u[i] = Vector<T, NDIMS>();
    this->m_u[i][i] = static_cast<T>(1.);

    // initialize ith extent to 0
    this->m_e[i] = T();
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS>::OrientedBoundingBox(const PointType* pts, int n)
{
  if(n <= 0)
  {
    this->clear();
    return;
  }

  numerics::Matrix<T> covar = numerics::Matrix<T>(NDIMS, NDIMS);
  axom::NumericArray<T, NDIMS> c;  // centroid

  for(int i = 0; i < n; i++)
  {
    c += pts[i].array();
  }
  c /= static_cast<T>(n);

  // save space for pts minus the centroid
  axom::NumericArray<T, NDIMS> diff;

  for(int i = 0; i < n; i++)
  {
    diff = pts[i].array() - c;
    for(int j = 0; j < NDIMS; j++)
    {
      for(int k = 0; k < NDIMS; k++)
      {
        covar(j, k) += diff[j] * diff[k];
      }
    }
  }

  // average the covariance matrix
  numerics::matrix_scalar_multiply(covar, static_cast<T>(1. / n));

  // make room for the eigenvectors and eigenvalues
  T u[NDIMS * NDIMS];
  T lambdas[NDIMS];

  int eigen_res = numerics::eigen_solve<T>(covar, NDIMS, u, lambdas);
  AXOM_UNUSED_VAR(eigen_res);
  SLIC_ASSERT(eigen_res);

  // save the axes
  for(int i = 0; i < NDIMS; ++i)
  {
    this->m_u[i] = Vector<T, NDIMS>(u + NDIMS * i);
  }

  // compute the extents
  Vector<T, NDIMS> maxima;
  T dot;
  for(int i = 0; i < n; ++i)
  {
    for(int j = 0; j < NDIMS; ++j)
    {
      diff = pts[i].array() - c;
      dot = utilities::abs<T>(numerics::dot_product<T>(&(m_u[j][0]), &diff[0], NDIMS));
      if(maxima[j] < dot)
      {
        maxima[j] = dot;
      }
    }
  }

  // save the extents
  this->m_e = maxima;

  // save the centroid
  this->m_c = Point<T, NDIMS>(c);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS>::OrientedBoundingBox(const Point<T, NDIMS>& c,
                                                   const Vector<T, NDIMS> (&u)[NDIMS],
                                                   const Vector<T, NDIMS>& e)
{
  this->m_c = Point<T, NDIMS>(c);
  for(int i = 0; i < NDIMS; i++)
  {
    this->m_u[i] = Vector<T, NDIMS>(u[i]);
  }
  this->m_e = Vector<T, NDIMS>(e);
  this->checkAndFix();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS>::OrientedBoundingBox(const OrientedBoundingBox& other)
{
  this->m_c = other.getCentroid();
  this->m_e = other.getExtents();
  const Vector<T, NDIMS>* u = other.getAxes();
  for(int i = 0; i < NDIMS; i++)
  {
    this->m_u[i] = u[i];
  }
  this->checkAndFix();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void OrientedBoundingBox<T, NDIMS>::clear()
{
  (this->m_u[0])[0] = axom::numeric_limits<T>::lowest();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void OrientedBoundingBox<T, NDIMS>::addPoint(Point<T, NDIMS> pt)
{
  Point<T, NDIMS> pt_l = this->toLocal(pt);
  for(int i = 0; i < NDIMS; i++)
  {
    T proj = pt_l[i];
    if(proj < T())
    {
      proj = -proj;
    }

    if(proj > this->m_e[i])
    {
      this->m_e[i] = proj;
    }
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void OrientedBoundingBox<T, NDIMS>::addBox(OrientedBoundingBox<T, NDIMS> obb)
{
  if(!obb.isValid())
  {
    // don't do anything
    return;
  }

  for(const auto& vert : obb.vertices())
  {
    this->addPoint(vert);
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::vector<Point<T, NDIMS>> OrientedBoundingBox<T, NDIMS>::vertices() const
{
  std::vector<Point<T, NDIMS>> res;
  res.reserve(1 << NDIMS);
  Vector<T, NDIMS> curr(this->m_c);

  OrientedBoundingBox<T, NDIMS>::vertex_enum(res, 0, this->m_u, this->m_e, curr);
  return res;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS>& OrientedBoundingBox<T, NDIMS>::expand(T expansionAmount)
{
  for(int i = 0; i < NDIMS; i++)
  {
    this->m_e[i] += expansionAmount;
  }
  this->checkAndFix();

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS>& OrientedBoundingBox<T, NDIMS>::scale(double scaleFactor)
{
  if(scaleFactor < T())
  {
    scaleFactor = -scaleFactor;
  }
  for(int i = 0; i < NDIMS; i++)
  {
    this->m_e[i] *= scaleFactor;
  }
  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS>& OrientedBoundingBox<T, NDIMS>::shift(const Vector<T, NDIMS>& displacement)
{
  axom::NumericArray<T, NDIMS>& coords = (this->m_c).array();
  axom::NumericArray<T, NDIMS> c;
  for(int i = 0; i < NDIMS; i++)
  {
    c[i] = coords[i] + displacement[i];
  }
  this->m_c = Point<T, NDIMS>(c);
  checkAndFix();
  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS>& OrientedBoundingBox<T, NDIMS>::operator=(const OrientedBoundingBox& rhs)
{
  if(this != &rhs)
  {
    this->m_e = rhs.getExtents();
    this->m_c = rhs.getCentroid();
    const Vector<T, NDIMS>* u = rhs.getAxes();
    for(int i = 0; i < NDIMS; i++)
    {
      this->m_u[i] = u[i];
    }
    this->checkAndFix();
  }
  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
template <typename OtherType>
bool OrientedBoundingBox<T, NDIMS>::contains(const Point<OtherType, NDIMS>& otherPt, double EPS) const
{
  if(!(this->isValid()))
  {
    return false;
  }
  Vector<T, NDIMS> pt_l(this->toLocal(otherPt));
  T proj;
  T margin = static_cast<T>(EPS);
  for(int i = 0; i < NDIMS; i++)
  {
    proj = pt_l[i];

    // make proj nonnegative
    if(proj < T())
    {
      proj = -proj;
    }

    if(proj > (m_e[i] + margin))
    {
      return false;
    }
  }
  return true;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
template <typename OtherType>
bool OrientedBoundingBox<T, NDIMS>::contains(const OrientedBoundingBox<OtherType, NDIMS>& otherOBB,
                                             double EPS) const
{
  std::vector<Point<T, NDIMS>> l = otherOBB.vertices();

  int size = static_cast<int>(l.size());

  for(int i = 0; i < size; i++)
  {
    if(!this->contains(l[i], EPS))
    {
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Point<T, NDIMS> OrientedBoundingBox<T, NDIMS>::toLocal(const Point<T, NDIMS>& pt) const
{
  Vector<T, NDIMS> d(m_c, pt);

  Point<T, NDIMS> res;
  for(int i = 0; i < NDIMS; i++)
  {
    res[i] = d.dot(this->m_u[i]);
  }
  return res;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Point<T, NDIMS> OrientedBoundingBox<T, NDIMS>::furthestPoint(const PointType& pt) const
{
  Point<T, NDIMS> pt_l = this->toLocal(pt);
  Vector<T, NDIMS> res;
  for(int i = 0; i < NDIMS; i++)
  {
    // since the local coordinates are individually constrained, we can simply
    // maximize the objective in each direction independently, meaning there
    // is a simple analytical solution

    if(pt_l[i] > T())
    {
      res -= (this->m_e[i]) * (this->m_u[i]);
    }
    else
    {
      res += (this->m_e[i]) * (this->m_u[i]);
    }
  }
  return Point<T, NDIMS>(res.array());
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
bool OrientedBoundingBox<T, NDIMS>::isValid() const
{
  return !((this->m_u[0])[0] == axom::numeric_limits<T>::lowest());
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void OrientedBoundingBox<T, NDIMS>::bisect(OrientedBoxType& right, OrientedBoxType& left) const
{
  T max = this->m_e[0];
  int bi = 0;
  for(int i = 1; i < NDIMS; i++)
  {
    if(max < this->m_e[i])
    {
      max = this->m_e[i];
      bi = i;
    }
  }

  Vector<T, NDIMS> c1 = Vector<T, NDIMS>(this->m_c);
  Vector<T, NDIMS> c2 = Vector<T, NDIMS>(this->m_c);

  c1 += (static_cast<T>(0.5) * this->m_e[bi]) * m_u[bi];
  c2 -= (static_cast<T>(0.5) * this->m_e[bi]) * m_u[bi];

  Vector<T, NDIMS> e(this->m_e);
  e[bi] *= (static_cast<T>(0.5));

  right = OrientedBoundingBox<T, NDIMS>(Point<T, NDIMS>(c1.array()), this->m_u, e);
  left = OrientedBoundingBox<T, NDIMS>(Point<T, NDIMS>(c2.array()), this->m_u, e);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& OrientedBoundingBox<T, NDIMS>::print(std::ostream& os) const
{
  // TODO: see if I can improve this
  os << "{centroid:" << this->m_c << "; extents:" << this->m_e;
  for(int i = 0; i < NDIMS; i++)
  {
    os << " axis " << i << ":" << this->m_u[i];
  }
  os << "}";
  return os;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void OrientedBoundingBox<T, NDIMS>::checkAndFix()
{
  // if this is invalid, don't do anything!
  if(!this->isValid())
  {
    return;
  }

  // do Gram-Schmidt
  for(int i = 0; i < NDIMS; i++)
  {
    // make orthogonal
    for(int j = 0; j < i; j++)
    {
      numerics::make_orthogonal<T>(this->m_u[i].data(), this->m_u[j].data(), NDIMS);
    }

    bool res = numerics::normalize<T>(this->m_u[i].data(), NDIMS);

    if(!res)  // it's nearly 0, so in the span of the others; set it to 0
    {
      this->m_u[i] = Vector<T, NDIMS>();
    }
  }

  // check and fix extents
  for(int i = 0; i < NDIMS; i++)
  {
    if(this->m_e[i] < T())
    {
      this->m_e[i] = T();
    }
  }
}

//------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
//------------------------------------------------------------------------------

template <typename T, int NDIMS>
bool operator==(const OrientedBoundingBox<T, NDIMS>& lhs, const OrientedBoundingBox<T, NDIMS>& rhs)
{
  return lhs.contains(rhs) && rhs.contains(lhs);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
bool operator!=(const OrientedBoundingBox<T, NDIMS>& lhs, const OrientedBoundingBox<T, NDIMS>& rhs)
{
  return !(lhs == rhs);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const OrientedBoundingBox<T, NDIMS>& box)
{
  return box.print(os);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_ORIENTEDBOUNDINGBOX_HPP_
