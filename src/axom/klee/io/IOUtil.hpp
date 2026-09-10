// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/klee/Dimensions.hpp"
#include "axom/klee/Units.hpp"

#include <tuple>

namespace axom
{
namespace inlet
{
class Container;
class Proxy;
class VerifiableScalar;
}  // namespace inlet

namespace klee
{
namespace internal
{
/**
 * Get the start and end units in a Container.
 *
 * The Container may either have a "units" field, or a "start_units" and "end_units".
 * In the first case, "units" will be used for both the start and end. In the
 * second, both must be present. In the case where no units are present at
 * all, both returned units will be LengthUnit::unspecified.
 *
 * \param container the Container from which to get the units
 * \return the start and end units
 * \throws KleeError if an invalid combination of fields is specified
 */
std::tuple<LengthUnit, LengthUnit> getOptionalStartAndEndUnits(const inlet::Container& container);

/**
 * Get the start and end units in a Container.
 *
 * The Container may either have a "units" field, or a "start_units" and "end_units".
 * In the first case, "units" will be used for both the start and end. In the
 * second, both must be present.
 *
 * \param container the Container from which to get the units
 * \return the start and end units
 * \throws KleeError if an invalid combination of fields is
 * specified or if no units are specified.
 */
std::tuple<LengthUnit, LengthUnit> getStartAndEndUnits(const inlet::Container& container);

/**
 * Define the schema for units. This is the schema that will be
 * expected by getOptionalStartAndEndUnits() and getStartAndEndUnits().
 *
 * @param container the container to which to add the expected fields
 * @param unitsDescription the description of the "units" field
 * @param startUnitsDescription the description of the "start_units" field
 * @param endUnitsDescription the description of the "end_units" field
 */
void defineUnitsSchema(inlet::Container& container,
                       const char* unitsDescription = "",
                       const char* startUnitsDescription = "",
                       const char* endUnitsDescription = "");

/**
 * Define a field which can hold a number of dimensions
 *
 * @param parent the parent Container on which to define the field
 * @param name the name of the field
 * @param description and optional description of the field
 * @return the field, which can have additional restrictions set on it
 */
inlet::VerifiableScalar& defineDimensionsField(inlet::Container& parent,
                                               const char* name,
                                               const char* description = "");

/**
 * Convert the given proxy to a Dimensions object. The field should have been
 * created by defineDimensionsField()
 *
 * @param dimProxy the proxy to the dimensions field
 * @return the value of the dimensions
 */
Dimensions toDimensions(const inlet::Proxy& dimProxy);

}  // namespace internal
}  // namespace klee
}  // namespace axom
