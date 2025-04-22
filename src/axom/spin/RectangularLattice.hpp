// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_RECTANGULAR_LATTICE_HPP_
#define AXOM_SPIN_RECTANGULAR_LATTICE_HPP_

#include "axom/config.hpp"
#include "axom/core/NumericArray.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/constants.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <cmath>
#include <iostream>

namespace axom
{
namespace spin
{
/*!
 * \class RectangularLattice
 * \brief A rectangular lattice maps all of space (of dimension NDIMS)
 *        into a rectangular grid of cells identified by integer coordinates.
 *
 * \tparam NDIMS the dimension of the space
 * \tparam SpaceCoordType the type of the spatial coordinates. Default: double
 * \tparam CellCoordType the type of the integral lattice coordinates.
 *         Default: int
 *
 * A rectangular lattice maps every point in space (a SpacePoint) to a unique
 * rectangular cell with integral coordinates (a GridCell). Each lattice cell is
 * identified with the integer grid point at its lowest corner in each
 * dimension.
 *
 * RectangularLattice also maps GridCell coordinates back to spatial coordinates
 * or to the bounding box of the GridCell.
 *
 * GridCells follow a half-open boundary convention. Points on their lower
 * boundaries get mapped to the cell, while points on their upper boundaries
 * get mapped to neighboring cells.
 *
 * A RectangularLattice is defined by an origin (a SpacePoint)
 * and a grid spacing (a SpaceVector).
 *
 * \note Grid spacing coordinates that are really small (magnitude less
 * than primal::PRIMAL_TINY := 1e-50) are snapped to zero to avoid division by zero.
 * \sa primal::PRIMAL_TINY
 */
template <int NDIMS, typename SpaceCoordType = double, typename CellCoordType = int>
class RectangularLattice
{
public:
  using GridCell = primal::Point<CellCoordType, NDIMS>;
  using SpacePoint = primal::Point<SpaceCoordType, NDIMS>;
  using SpaceVector = primal::Vector<SpaceCoordType, NDIMS>;
  using SpatialBoundingBox = primal::BoundingBox<SpaceCoordType, NDIMS>;

public:
  /*!
   * \brief Default constructor
   *
   * \note Sets the origin to 0 and spacing to 1 in each dimension
   */
  RectangularLattice()
    : m_origin(SpaceCoordType(0))
    , m_spacing(SpacePoint(SpaceCoordType(1)))
    , m_invSpacing(SpacePoint(SpaceCoordType(1)))
  { }

  /*!
   * \brief Constructor from a given origin.
   *
   * \param origin The lattice's origin
   * \note Spacing will be set to default (vector of ones)
   */
  RectangularLattice(const SpacePoint& origin)
    : m_origin(origin)
    , m_spacing(SpacePoint(SpaceCoordType(1)))
    , m_invSpacing(SpacePoint(SpaceCoordType(1)))
  { }

  /*!
   * \brief Constructor from a given origin and spacing.
   *
   * \param origin The lattice's origin
   * \param spacing The lattice's spacing
   *
   * \note The magnitude of the spacing coordinates should be greater than zero.
   * If they are less than primal::PRIMAL_TINY, the lattice will be degenerate in
   * that dimension.
   */
  RectangularLattice(const SpacePoint& origin, const SpaceVector& spacing)
    : m_origin(origin)
    , m_spacing(spacing)
  {
    // Note: We use an inverted spacing, m_invSpacing, for efficiency.
    // It trades divisions for multiplications and handles 0-sized spacings
    // The following helper function sets this up.
    initializeSpacingAndInvSpacing();
  }

  /*!
   * \brief Constructor from SpaceCoordType arrays
   *
   * \param origin_date An array containing the origin's coordinates
   * \param spacing_data An array containing the spacing coordinates
   *
   * \pre When origin_data is not NULL, it must have at least NDIMS entries
   * \note Origin will be set to the zero vector if origin_data pointer is NULL
   *
   * \pre When spacing_data is not NULL, it must have at least NDIMS entries
   * \note Spacing will be set to vector or ones if pointer is NULL
   *
   * \note The magnitude of the spacing coordinates should be greater than zero.
   * If they are less than EPS = primal::PRIMAL_TINY, the lattice will be 
   * degenerate in that dimension.
   */
  RectangularLattice(SpaceCoordType* origin_data, SpaceCoordType* spacing_data)
  {
    m_origin = (origin_data != nullptr) ? SpacePoint(origin_data) : SpacePoint::zero();

    m_spacing =
      (spacing_data != nullptr) ? SpaceVector(spacing_data) : SpaceVector(SpaceCoordType(1));

    // Note: We use an inverted spacing, m_invSpacing, for efficiency.
    // It trades divisions for multiplications and handles 0-sized spacings
    // The following helper function sets this up.
    initializeSpacingAndInvSpacing();
  }

  /*! Accessor for lattice origin   */
  const SpacePoint& origin() const { return m_origin; }

  /*! Accessor for lattice spacing   */
  const SpaceVector& spacing() const { return m_spacing; }

  /*! Returns the lattice cell associated with the given space point pt */
  AXOM_HOST_DEVICE
  GridCell gridCell(const SpacePoint& pt) const
  {
    GridCell cell;

    for(int i = 0; i < NDIMS; ++i)
    {
      // Note: Always round down to negative infinity
      // uses inverted spacing, which handles zero-sized spacings
      cell[i] =
        static_cast<CellCoordType>(axom::utilities::floor((pt[i] - m_origin[i]) * m_invSpacing[i]));
    }

    return cell;
  }

  /*! Returns the space point associated with lowest corner of GridCell cell */
  SpacePoint spacePoint(const GridCell& cell) const
  {
    SpacePoint pt;
    for(int i = 0; i < NDIMS; ++i)
    {
      pt[i] = m_origin[i] + (m_spacing[i] * cell[i]);
    }

    return pt;
  }

  /*!
   * \brief Find the spatial bounding box of a lattice cell
   *
   * \note The bounding box covers all points in space that map to gridCell
   */
  SpatialBoundingBox cellBounds(const GridCell& cell) const
  {
    return SpatialBoundingBox(spacePoint(cell),
                              spacePoint(GridCell(cell.array() + GridCell(1).array())));
  }

  /*! Simple formatted print of a rectangular lattice */
  std::ostream& print(std::ostream& os) const
  {
    os << "{ origin:" << m_origin << "; spacing:" << m_spacing << " }";
    return os;
  }

private:
  /*!
   * \brief Helper function to snap really small coordinates for the grid
   * spacing to zero and to initialize the inverted spacing.
   *
   * A spacing coordinate is considered really small when its magnitude
   * is less than primal::PRIMAL_TINY.
   *
   * For each coordinate i, the inverted coordinate will be:
   *     m_invSpacing[i] = 1. / m_spacing[i]
   *
   * as long as m_spacing is sufficiently large. Otherwise, it will be zero.
   */
  void initializeSpacingAndInvSpacing()
  {
    constexpr SpaceCoordType EPS = primal::PRIMAL_TINY;
    constexpr SpaceCoordType ZERO = SpaceCoordType(0.);
    constexpr SpaceCoordType ONE = SpaceCoordType(1.);

    for(int i = 0; i < NDIMS; ++i)
    {
      // snap really small values to zero
      if(axom::utilities::isNearlyEqual(m_spacing[i], ZERO, EPS))
      {
        m_spacing[i] = ZERO;
      }

      // compute the inverted spacing coordinate
      // It is OK to compare to zero here due to snapping
      m_invSpacing[i] = (m_spacing[i] != ZERO) ? ONE / m_spacing[i] : ZERO;
    }
  }

private:
  SpacePoint m_origin;       /// Origin of the grid
  SpaceVector m_spacing;     /// Spacing of the grid
  SpaceVector m_invSpacing;  /// Inverse spacing (for efficiency)
};

/*!
 * \brief Helper function to create a Rectangular Lattice from
 * a supplied bounding box and grid resolution
 *
 * \param bbox The bounding box from which to construct a RectangularLattice
 * \param gridRes The resolution of the input grid
 *
 * This is a convenience function to simplify working with bounding boxes.
 * It extracts the lattice spacing from the supplied bounding box
 * and resolution, and sets the lattice origin to the bounding box's
 * minimum corner position.
 *
 * \note If the bounding box range along a dimension is near zero (i.e. smaller
 * than primal::PRIMAL_TINY, the grid resolution in that dimension will be 
 * set to zero in that dimension.
 */
template <int NDIMS, typename SpaceCoordType, typename CellCoordType>
RectangularLattice<NDIMS, SpaceCoordType, CellCoordType> rectangular_lattice_from_bounding_box(
  const primal::BoundingBox<SpaceCoordType, NDIMS>& bbox,
  const NumericArray<CellCoordType, NDIMS>& gridRes)
{
  using LatticeType = RectangularLattice<NDIMS, SpaceCoordType, CellCoordType>;
  using SpaceVector = typename LatticeType::SpaceVector;

  SpaceVector spacing;

  // Use the resolution and the box range to compute the spacing
  for(int i = 0; i < NDIMS; ++i)
  {
    if(gridRes[i] != CellCoordType(0))
    {
      spacing[i] = (bbox.getMax()[i] - bbox.getMin()[i]) / gridRes[i];
    }
  }

  return LatticeType(bbox.getMin(), spacing);
}

/// --------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
/// --------------------------------------------------------------------------------

/*! Equality operator on two RectangularLattices */
template <int NDIMS, typename SpaceCoordType, typename CellCoordType>
bool operator==(const RectangularLattice<NDIMS, SpaceCoordType, CellCoordType>& lhs,
                const RectangularLattice<NDIMS, SpaceCoordType, CellCoordType>& rhs)
{
  return lhs.origin() == rhs.origin() && lhs.spacing() == rhs.spacing();
}

/*! Inequality operator on two RectangularLattices */
template <int NDIMS, typename SpaceCoordType, typename CellCoordType>
bool operator!=(const RectangularLattice<NDIMS, SpaceCoordType, CellCoordType>& lhs,
                const RectangularLattice<NDIMS, SpaceCoordType, CellCoordType>& rhs)
{
  return !(lhs == rhs);
}

/*! Stream output operator on a RectangularLattice */
template <int NDIMS, typename SpaceCoordType, typename CellCoordType>
std::ostream& operator<<(std::ostream& os,
                         const RectangularLattice<NDIMS, SpaceCoordType, CellCoordType>& lattice)
{
  return lattice.print(os);
}

}  // end namespace spin
}  // end namespace axom

#endif  // AXOM_SPIN_RECTANGULAR_LATTICE_HPP_
