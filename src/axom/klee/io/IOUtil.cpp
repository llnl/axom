// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "IOUtil.hpp"

#include "axom/inlet.hpp"
#include "axom/klee/KleeError.hpp"

namespace axom
{
namespace klee
{
namespace internal
{
std::tuple<LengthUnit, LengthUnit> getOptionalStartAndEndUnits(const inlet::Container& container)
{
  bool hasStartUnits = container.contains("start_units");
  bool hasEndUnits = container.contains("end_units");
  if(container.contains("units"))
  {
    if(hasStartUnits || hasEndUnits)
    {
      throw KleeError({container.name(), "Can't specify 'units' with 'start_units' or 'end_units'"});
    }
    auto units = internal::parseLengthUnits(container["units"]);
    return std::make_tuple(units, units);
  }
  else if(hasStartUnits || hasEndUnits)
  {
    if(!(hasStartUnits && hasEndUnits))
    {
      throw KleeError({container.name(), "Must specify both 'start_units' and 'end_units'"});
    }
    auto startUnits = internal::parseLengthUnits(container["start_units"]);
    auto endUnits = internal::parseLengthUnits(container["end_units"]);
    return std::make_tuple(startUnits, endUnits);
  }
  return std::make_tuple(LengthUnit::unspecified, LengthUnit::unspecified);
}

std::tuple<LengthUnit, LengthUnit> getStartAndEndUnits(const inlet::Container& container)
{
  auto units = getOptionalStartAndEndUnits(container);
  if(std::get<0>(units) == LengthUnit::unspecified)
  {
    throw KleeError({container.name(), "Did not specify units"});
  }
  return units;
}

void defineUnitsSchema(inlet::Container& container,
                       const char* unitsDescription,
                       const char* startUnitsDescription,
                       const char* endUnitsDescription)
{
  container.addString("start_units", startUnitsDescription);
  container.addString("end_units", endUnitsDescription);
  container.addString("units", unitsDescription);
  // Don't do custom validator here because getOptionalStartAndEndUnits()
  // verifies the right combination is specified. If we were to add a
  // custom validator, we would have to repeat some of the logic when
  // figuring out which fields to use.
}

inlet::VerifiableScalar& defineDimensionsField(inlet::Container& parent,
                                               const char* name,
                                               const char* description)
{
  return parent.addInt(name, description).range(2, 3);
}

Dimensions toDimensions(const inlet::Proxy& dimProxy)
{
  return static_cast<Dimensions>(dimProxy.get<int>());
}

}  // namespace internal
}  // namespace klee
}  // namespace axom
