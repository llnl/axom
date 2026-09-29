// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Units.hpp"

#include "axom/inlet.hpp"
#include "axom/klee/KleeError.hpp"

#include <stdexcept>

namespace axom
{
namespace klee
{
namespace internal
{
/**
 * Parse length units from a string.
 *
 * \param unitsAsString the unit string to parse
 * \param path the input path to report on failure
 * \return the parsed length unit
 * \throws KleeError if the unit string is invalid
 */
LengthUnit parseLengthUnits(const std::string& unitsAsString, const std::string& path)
{
  try
  {
    return utilities::getLengthUnit(unitsAsString);
  }
  catch(const std::invalid_argument& ex)
  {
    throw KleeError({path, ex.what()});
  }
}

LengthUnit parseLengthUnits(const inlet::Proxy& unitsAsProxy)
{
  return parseLengthUnits(unitsAsProxy.get<std::string>(), unitsAsProxy.name());
}
}  // namespace internal
}  // namespace klee
}  // namespace axom
