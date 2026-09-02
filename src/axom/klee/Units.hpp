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
 * \brief Parse a unit string, reporting an invalid value at the supplied input path.
 *
 * \param [in] unitsAsString the unit string to parse
 * \param [in] path the input path to report on failure
 *
 * \return A LengthUnit containing the unit type.
 * \throws KleeError if the unit string is invalid
 */
LengthUnit parseLengthUnits(const std::string& unitsAsString, const std::string& path);

/*!
 * \brief This function parses a string and returns a LengthUnit. It is a compatibility
 *        function that throws a KleeError if the unit is invalid.
 *
 * \param unitsAsProxy The Inlet proxy from which to get the unit string.
 *
 * \return A LengthUnit containing the unit type.
 * \throws KleeError if the unit string is invalid
 */
LengthUnit parseLengthUnits(const inlet::Proxy& unitsAsProxy);

}  // namespace internal
}  // namespace klee
}  // namespace axom
