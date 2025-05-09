// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FEM_BASIS_HPP_
#define MINT_FEM_BASIS_HPP_

// Mint includes
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/mint/fem/shape_functions/Lagrange.hpp"
#include "axom/mint/fem/shape_functions/ShapeFunction.hpp"
#include "axom/mint/fem/FEBasisTypes.hpp"

/*!
 * \def REGISTER_LAGRANGE_BASIS( C )
 *
 * \brief Macro used to register a new Lagrange Finite Element. This macro
 *  expands to a specialization of the FEBasis trait class that binds a
 *  Finite Element basis to a cell type.
 */
#define REGISTER_LAGRANGE_BASIS(C)                                    \
  template <>                                                         \
  struct FEBasis<MINT_LAGRANGE_BASIS, C>                              \
  {                                                                   \
    static const int BasisType = MINT_LAGRANGE_BASIS;                 \
    typedef mint::ShapeFunction<mint::Lagrange<C>> ShapeFunctionType; \
  }

namespace axom
{
namespace mint
{
/*!
 * \brief FEBasis is a traits class that binds a Finite Element basis type,
 *  e.g., MINT_LAGRANGE_BASIS, to a particular cell type, e.g., MINT_QUAD.
 *
 * \note This is an empty class that is intended to be specialized
 *
 * \see ShapeFunction
 * \see FEBasisTypes
 */
template <int BasisType, CellType CELLTYPE>
struct FEBasis
{ };

// Lagrange Basis
REGISTER_LAGRANGE_BASIS(mint::QUAD);
REGISTER_LAGRANGE_BASIS(mint::TRIANGLE);
REGISTER_LAGRANGE_BASIS(mint::TET);
REGISTER_LAGRANGE_BASIS(mint::HEX);
REGISTER_LAGRANGE_BASIS(mint::PRISM);
REGISTER_LAGRANGE_BASIS(mint::PYRAMID);

REGISTER_LAGRANGE_BASIS(mint::QUAD9);
REGISTER_LAGRANGE_BASIS(mint::HEX27);

} /* namespace mint */
} /* namespace axom */

// undef internal macros
#undef REGISTER_LAGRANGE_BASIS

#endif /* MINT_FEM_BASIS_HPP_ */
