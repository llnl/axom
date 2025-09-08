// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_UNITS_HPP
#define AXOM_KLEE_UNITS_HPP

#include <string>

namespace axom
{
namespace inlet
{
class Proxy;
}
namespace klee
{
/**
 * Units of length in which users can express lengths and in which client
 * codes can request them.
 */
enum class LengthUnit
{
  km,
  m,
  dm,
  cm,
  mm,
  um,
  nm,
  angstrom,
  miles,
  feet,
  inches,
  mils,
  unspecified
};

/**
 * Convert a string to a LengthUnit.
 *
 * \param unitsAsString the units as a string
 * \param path the Path where the length units were found in the document.
 * Used for error reporting.
 * \return the parsed units
 * \throws KleeError if the string does not represent known
 * units
 */
LengthUnit parseLengthUnits(const std::string &unitsAsString, const std::string &path);

/**
 * Convert a proxy to a LengthUnit.
 *
 * \param unitsAsProxy the units as a proxy
 * \return the parsed units
 * \throws KleeError if the string does not represent known units
 */
LengthUnit parseLengthUnits(const inlet::Proxy &unitsAsProxy);

/**
 * Get the conversion factor to convert from the given source units to the target units.
 *
 * \param sourceUnits the original units
 * \param targetUnits the units to convert to
 * \return the value by which to multiply lengths in the original units
 * to get the target units
 */
double getConversionFactor(LengthUnit sourceUnits, LengthUnit targetUnits);

/**
 * Convert a value from one set of units to another.
 *
 * \param sourceValue the value of the length in the original units
 * \param sourceUnits the original units
 * \param targetUnits the units to convert to
 * \return the value of the length in the new units
 */
double convert(double sourceValue, LengthUnit sourceUnits, LengthUnit targetUnits);

/**
 * Convert multiple lengths in place.
 *
 * \tparam T the type containing the units. Must be iterable.
 *
 * \param values the value of the length in the original units
 * \param sourceUnits the original units
 * \param targetUnits the units to convert to
 */
template <typename T>
void convertAll(T &values, LengthUnit sourceUnits, LengthUnit targetUnits)
{
  double factor = getConversionFactor(sourceUnits, targetUnits);
  for(double &value : values)
  {
    value *= factor;
  }
}

}  // namespace klee
}  // namespace axom
#endif  // AXOM_KLEE_UNITS_HPP
