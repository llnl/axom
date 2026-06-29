// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#pragma once

#include "axom/core/utilities/Units.hpp"

#include <string>

namespace axom
{
namespace inlet
{
class Proxy;
}
namespace klee
{

using utilities::LengthUnit;

namespace internal
{

/*!
 * \brief This function parses a string and returns a LengthUnit. It is a compatibility
 *        function that throws a KleeError if the unit is invalid.
 *
 * \param unitsAsProxy The Inlet proxy from which to get the unit string.
 *
 * \return A LengthUnit containing the unit type.
 */
LengthUnit parseLengthUnits(const inlet::Proxy &unitsAsProxy);

}  // namespace internal
}  // namespace klee
}  // namespace axom

