// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_LAGRANGE_SHAPEFUNCTION_HPP_
#define MINT_LAGRANGE_SHAPEFUNCTION_HPP_

// Axom includes
#include "axom/core/Macros.hpp"

// Mint includes
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/mint/fem/FEBasisTypes.hpp"
#include "axom/mint/fem/shape_functions/ShapeFunction.hpp"

// Slic includes
#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace mint
{
/*!
 * \brief Defines the Lagrange family of Finite Elements
 *
 * \tparam CELLTYPE the cell type of the element, e.g., mint::QUAD, etc.
 *
 * \note This is the default implementation. Only stubs are defined at this
 *  level.This class is specialized according to cell type.
 */
template <CellType CELLTYPE>
class Lagrange : public ShapeFunction<Lagrange<CELLTYPE>>
{
public:
  /*!
   * \brief Returns the cell type of this instance.
   * \return c the cell type, e.g., MINT_QUAD etc.
   * \see CellTypes.hpp
   *
   * \note This method is implemented in specialized instances.
   */
  static CellType getCellType()
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
    return CELLTYPE;
  }

  /*!
   * \brief Returns the Basis function type.
   * \return b the basis thpe, e.g., MINT_LAGRANGE_BASIS
   * \see FEBasisTypes.hpp
   *
   * \note This method is implemented in specialized instances.
   */
  static int getType()
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
    return MINT_UNDEFINED_BASIS;
  }

  /*!
   * \brief Returns the number of degrees of freedom on this Finite Element.
   * \return nodfs the number of degrees of freedom.
   *
   * \note This method is implemented in specialized instances.
   */
  static int getNumDofs()
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the max number of newton iterations for this Finite Element
   * \return N the max number of newton iterations
   * \post N >= 1
   *
   * \note This method is implemented in specialized instances.
   */
  static int getMaxNewtonIters()
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the dimension of the reference element.
   * \return dim the dimension of the reference element
   * \post dim >= 1 && dim <= 3
   *
   * \note This method is implemented in specialized instances.
   */
  static int getDimension()
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the min coordinate of the element's reference space
   * \return min the min coordinate of the element's reference space, e.g., 0
   *
   * \note This method is implemented in specialized instances.
   */
  static int getMin()
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the max coordinate of the element's reference space
   * \return max the max coordinate of the element's reference space, e.g., 1
   *
   * \note This method is implemented in specialized instances.
   */
  static int getMax()
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the centroid of the reference element.
   *
   * \param [out] center ndims long buffer to store the centroid
   * \pre center != nullptr
   *
   * \note This method is implemented in specialized instances.
   */
  static void getCenter(double* AXOM_UNUSED_PARAM(center))
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
  }

  /*!
   * \brief Returns the coordinates of the reference element.
   *
   * \param [out] coords ndims*ndofs long buffer to store the coordinates
   * \pre coords != nullptr
   *
   * \note The coordinates are arranged in column-major flat array layout.
   * \note This method is implemented in specialized instances.
   */
  static void getCoords(double* AXOM_UNUSED_PARAM(coords))
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
  }

  /*!
   * \brief Computes the shape functions of the Finite Element at the given
   *  natural coordinates, \f$ \xi \f$
   *
   * \param [in] nc natural coordinates at which to compute the shape functions
   * \param [out] phi buffer (ndofs long) to store the shape functions
   *
   * \pre nc != nullptr
   * \pre phi != nullptr
   *
   * \note This method is implemented in specialized instances.
   */
  static void computeShape(const double* AXOM_UNUSED_PARAM(nc), double* AXOM_UNUSED_PARAM(phi))
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
  }

  /*!
   * \brief Computes the shape function first derivatives for the Finite
   *  Element at the given natural coordinates, \f$ \xi \f$
   *
   * \param [in] nc natural coordinates at which to compute the derivatives
   * \param [out] phidot buffer (ndofs*ndims long) for the derivatives
   *
   * \pre nc != nullptr
   * \pre phidot != nullptr
   *
   * \note This method is implemented in specialized instances.
   */
  static void computeDerivatives(const double* AXOM_UNUSED_PARAM(nc),
                                 double* AXOM_UNUSED_PARAM(phidot))
  {
    constexpr int cell_value = mint::cellTypeToInt(CELLTYPE);
    AXOM_STATIC_ASSERT(cell_value >= 0 && cell_value < mint::NUM_CELL_TYPES);
    Lagrange<CELLTYPE>::checkCellType();
  }

private:
  /*!
   * \brief Checks if the CELLTYPE is valid and supported in the Lagrange basis.
   */
  static void checkCellType()
  {
    if(CELLTYPE != UNDEFINED_CELL)
    {
      SLIC_ERROR("Lagrange shape functions for [" << getCellInfo(CELLTYPE).name
                                                  << "] are not defined!");
    }
    else
    {
      SLIC_ERROR("Invalid CellType: " << cellTypeToInt(CELLTYPE));
    }
  }
};

} /* namespace mint */
} /* namespace axom */

// Template Specializations for Lagrange
#include "axom/mint/fem/shape_functions/lagrange/lagrange_hexa_27.hpp"
#include "axom/mint/fem/shape_functions/lagrange/lagrange_hexa_8.hpp"
#include "axom/mint/fem/shape_functions/lagrange/lagrange_prism_6.hpp"
#include "axom/mint/fem/shape_functions/lagrange/lagrange_pyra_5.hpp"
#include "axom/mint/fem/shape_functions/lagrange/lagrange_quad_4.hpp"
#include "axom/mint/fem/shape_functions/lagrange/lagrange_quad_9.hpp"
#include "axom/mint/fem/shape_functions/lagrange/lagrange_tetra_4.hpp"
#include "axom/mint/fem/shape_functions/lagrange/lagrange_tri_3.hpp"

#endif /* MINT_LAGRANGE_SHAPEFUNCTION_HPP_ */
