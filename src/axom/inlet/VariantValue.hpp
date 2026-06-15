// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file VariantValue.hpp
 *
 * \brief This file contains Inlet's variant value type for mixed primitive
 * collections.
 *******************************************************************************
 */

#ifndef INLET_VARIANT_VALUE_HPP
#define INLET_VARIANT_VALUE_HPP

#include <string>
#include <variant>

namespace axom
{
namespace inlet
{

using VariantValue = std::variant<bool, int, double, std::string>;

}  // namespace inlet
}  // namespace axom

#endif
