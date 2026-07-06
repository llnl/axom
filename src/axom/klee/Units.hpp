// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_UNITS_HPP
#define AXOM_KLEE_UNITS_HPP

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

using utilities::getLengthUnitName;
using utilities::LengthUnit;

LengthUnit parseLengthUnits(const std::string &unitsAsString, const std::string &path);

LengthUnit parseLengthUnits(const inlet::Proxy &unitsAsProxy);

}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_UNITS_HPP
