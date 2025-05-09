// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FEM_BASIS_TYPES_HPP_
#define MINT_FEM_BASIS_TYPES_HPP_

#include <string>

/*!
 * \enum FEBasisTypes
 *
 * \brief Enumerates the different types of supported Finite Element basis.
 *
 * \see FEBasis
 */
enum FEBasisTypes
{
  MINT_UNDEFINED_BASIS = -1, /*!< Undefined basis type */
  MINT_LAGRANGE_BASIS,       /*!< Lagrange basis type */

  MINT_NUM_BASIS_TYPES /*!< MINT_NUM_BASIS */
};

namespace axom
{
namespace mint
{
/*!
 * \brief Array of strings corresponding to each Finite Element Basis.
 *
 * \note The length of the array is MINT_NUM_BASIS
 *
 * \note It is used to get a string representation for a Finite Element basis,
 *  e.g., for, debugging etc.
 */
static const std::string basis_name[] = {
  "LAGRANGE_BASIS",
};

}  // namespace mint
}  // namespace axom

#endif
