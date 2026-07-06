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
LengthUnit parseLengthUnits(const std::string &unitsAsString, const std::string &path)
{
  try
  {
    return utilities::getLengthUnit(unitsAsString);
  }
  catch(const std::invalid_argument &ex)
  {
    throw KleeError({path, ex.what()});
  }
}

LengthUnit parseLengthUnits(const inlet::Proxy &unitsAsProxy)
{
  return parseLengthUnits(unitsAsProxy.get<std::string>(), unitsAsProxy.name());
}

}  // namespace klee
}  // namespace axom
