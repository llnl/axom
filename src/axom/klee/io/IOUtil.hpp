// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_IO_UTIL_HPP_
#define AXOM_KLEE_IO_UTIL_HPP_

#include "axom/klee/Dimensions.hpp"
#include "axom/klee/Units.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <tuple>
#include <vector>

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
 * Convert the given field to a std::vector<double>, ensuring that it
 * has the expected number of entries.
 *
 * @param field the field to convert
 * @param expectedDims the expected dimensionality of the array
 * @param fieldName the name of the field (used for error reporting)
 * @return the field as a std::vector<double>
 */
std::vector<double> toDoubleVector(inlet::Proxy const &field,
                                   Dimensions expectedDims,
                                   char const *fieldName);

/**
 * Convert the specified field to a Point3D, ensuring that it
 * has the expected number of entries.
 *
 * @param parent the parent of the field
 * @param fieldName the name of the field
 * @param expectedDims the expected dimensionality of the point
 * @return the field as a primal::Point3D
 */
primal::Point3D toPoint(inlet::Container const &parent, char const *fieldName, Dimensions expectedDims);

/**
 * Convert the specified field to a Point3D, ensuring that it
 * has the expected number of entries. If the field is not present, the default value is used.
 *
 * @param parent the parent of the field
 * @param fieldName the name of the field
 * @param expectedDims the expected dimensionality of the point
 * @param defaultValue the default value of the field if it is not present
 * @return the field as a primal::Point3D
 */
primal::Point3D toPoint(inlet::Container const &parent,
                        char const *fieldName,
                        Dimensions expectedDims,
                        const primal::Point3D &defaultValue);

/**
 * Convert the specified field to a Vector3D, ensuring that it
 * has the expected number of entries.
 *
 * @param parent the parent of the field
 * @param fieldName the name of the field
 * @param expectedDims the expected dimensionality of the vector
 * @return the field as a primal::Vector3D
 */
primal::Vector3D toVector(inlet::Container const &parent,
                          char const *fieldName,
                          Dimensions expectedDims);

/**
 * Convert the specified field to a Vector3D, ensuring that it
 * has the expected number of entries. If the field is not present, the
 * default value is used.
 *
 * @param parent the parent of the field
 * @param fieldName the name of the field
 * @param expectedDims the expected dimensionality of the vector
 * @param defaultValue the default value of the field if it is not present
 * @return the field as a primal::Vector3D
 */
primal::Vector3D toVector(inlet::Container const &parent,
                          char const *fieldName,
                          Dimensions expectedDims,
                          const primal::Vector3D &defaultValue);

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
std::tuple<LengthUnit, LengthUnit> getOptionalStartAndEndUnits(const inlet::Container &container);

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
std::tuple<LengthUnit, LengthUnit> getStartAndEndUnits(const inlet::Container &container);

/**
 * Define the schema for units. This is the schema that will be
 * expected by getOptionalStartAndEndUnits() and getStartAndEndUnits().
 *
 * @param container the container to which to add the expected fields
 * @param unitsDescription the description of the "units" field
 * @param startUnitsDescription the description of the "start_units" field
 * @param endUnitsDescription the description of the "end_units" field
 */
void defineUnitsSchema(inlet::Container &container,
                       const char *unitsDescription = "",
                       const char *startUnitsDescription = "",
                       const char *endUnitsDescription = "");

/**
 * Define a field which can hold a number of dimensions
 *
 * @param parent the parent Container on which to define the field
 * @param name the name of the field
 * @param description and optional description of the field
 * @return the field, which can have additional restrictions set on it
 */
inlet::VerifiableScalar &defineDimensionsField(inlet::Container &parent,
                                               const char *name,
                                               const char *description = "");

/**
 * Convert the given proxy to a Dimensions object. The field should have been
 * created by defineDimensionsField()
 *
 * @param dimProxy the proxy to the dimensions field
 * @return the value of the dimensions
 */
Dimensions toDimensions(const inlet::Proxy &dimProxy);

}  // namespace internal
}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_IO_UTIL_HPP_
