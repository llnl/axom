// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_CONSTANTS_HPP_
#define AXOM_PRIMAL_CONSTANTS_HPP_

/*! 
 * \file constants.hpp
 * \brief This file contains some useful constants to use throughout primal
 */

namespace axom
{
namespace primal
{
/// Small constant that can be useful for e.g. offsetting denominators to avoid division by zero
static constexpr double PRIMAL_TINY = 1e-50;
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CONSTANTS_
