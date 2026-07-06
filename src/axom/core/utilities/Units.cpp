// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/Units.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace axom
{
namespace utilities
{
namespace
{
/**
 * A simple functor which can be used for hashing length units. Needed because
 * C++11 does not support hashing enums by default (fixed in 14, this can
 * go away when we commit to that standard).
 */
struct LengthUnitHash
{
  std::size_t operator()(LengthUnit unit) const { return static_cast<std::size_t>(unit); }
};

std::string toLower(std::string str)
{
  std::transform(str.begin(), str.end(), str.begin(), [](unsigned char ch) {
    return static_cast<char>(std::tolower(ch));
  });
  return str;
}

std::string unrecognizedUnitsMessage(const std::string &unitsAsString)
{
  std::string message = "Unrecognized units: ";
  message += unitsAsString;
  return message;
}
}  // namespace

LengthUnit getLengthUnit(const std::string &unit)
{
  const std::string lowerUnit = toLower(unit);

  static const std::unordered_map<std::string, LengthUnit> UNITS_BY_NAME {
    {"inch", LengthUnit::inches},
    {"inches", LengthUnit::inches},
    {"in", LengthUnit::inches},
    {"foot", LengthUnit::feet},
    {"feet", LengthUnit::feet},
    {"ft", LengthUnit::feet},
    {"mile", LengthUnit::miles},
    {"miles", LengthUnit::miles},
    {"mi", LengthUnit::miles},
    {"mil", LengthUnit::mils},
    {"mils", LengthUnit::mils},
    {"a", LengthUnit::angstrom},
    {"angstrom", LengthUnit::angstrom},
    {"angstroms", LengthUnit::angstrom},
    {"am", LengthUnit::am},
    {"attometer", LengthUnit::am},
    {"attometers", LengthUnit::am},
    {"attometre", LengthUnit::am},
    {"attometres", LengthUnit::am},
    {"fm", LengthUnit::fm},
    {"femtometer", LengthUnit::fm},
    {"femtometers", LengthUnit::fm},
    {"femtometre", LengthUnit::fm},
    {"femtometres", LengthUnit::fm},
    {"pm", LengthUnit::pm},
    {"picometer", LengthUnit::pm},
    {"picometers", LengthUnit::pm},
    {"picometre", LengthUnit::pm},
    {"picometres", LengthUnit::pm},
    {"nm", LengthUnit::nm},
    {"nanometer", LengthUnit::nm},
    {"nanometers", LengthUnit::nm},
    {"nanometre", LengthUnit::nm},
    {"nanometres", LengthUnit::nm},
    {"um", LengthUnit::um},
    {"micrometer", LengthUnit::um},
    {"micrometers", LengthUnit::um},
    {"micrometre", LengthUnit::um},
    {"micrometres", LengthUnit::um},
    {"micron", LengthUnit::um},
    {"microns", LengthUnit::um},
    {"mm", LengthUnit::mm},
    {"millimeter", LengthUnit::mm},
    {"millimeters", LengthUnit::mm},
    {"millimetre", LengthUnit::mm},
    {"millimetres", LengthUnit::mm},
    {"cm", LengthUnit::cm},
    {"centimeter", LengthUnit::cm},
    {"centimeters", LengthUnit::cm},
    {"centimetre", LengthUnit::cm},
    {"centimetres", LengthUnit::cm},
    {"dm", LengthUnit::dm},
    {"decimeter", LengthUnit::dm},
    {"decimeters", LengthUnit::dm},
    {"decimetre", LengthUnit::dm},
    {"decimetres", LengthUnit::dm},
    {"m", LengthUnit::m},
    {"meter", LengthUnit::m},
    {"meters", LengthUnit::m},
    {"metre", LengthUnit::m},
    {"metres", LengthUnit::m},
    {"dam", LengthUnit::dam},
    {"decameter", LengthUnit::dam},
    {"decameters", LengthUnit::dam},
    {"decametre", LengthUnit::dam},
    {"decametres", LengthUnit::dam},
    {"hm", LengthUnit::hm},
    {"hectometer", LengthUnit::hm},
    {"hectometers", LengthUnit::hm},
    {"hectometre", LengthUnit::hm},
    {"hectometres", LengthUnit::hm},
    {"km", LengthUnit::km},
    {"kilometer", LengthUnit::km},
    {"kilometers", LengthUnit::km},
    {"kilometre", LengthUnit::km},
    {"kilometres", LengthUnit::km},
  };

  auto iter = UNITS_BY_NAME.find(lowerUnit);
  if(iter == UNITS_BY_NAME.end())
  {
    throw std::invalid_argument(unrecognizedUnitsMessage(unit));
  }

  return iter->second;
}

std::string getLengthUnitName(LengthUnit unit)
{
  static const std::unordered_map<LengthUnit, std::string, LengthUnitHash> UNIT_NAMES {
    {LengthUnit::am, "am"},
    {LengthUnit::fm, "fm"},
    {LengthUnit::pm, "pm"},
    {LengthUnit::nm, "nm"},
    {LengthUnit::um, "um"},
    {LengthUnit::mm, "mm"},
    {LengthUnit::cm, "cm"},
    {LengthUnit::dm, "dm"},
    {LengthUnit::m, "m"},
    {LengthUnit::dam, "dam"},
    {LengthUnit::hm, "hm"},
    {LengthUnit::km, "km"},
    {LengthUnit::miles, "miles"},
    {LengthUnit::feet, "ft"},
    {LengthUnit::inches, "in"},
    {LengthUnit::mils, "mils"},
    {LengthUnit::angstrom, "A"},
    {LengthUnit::unspecified, "unspecified"}};

  auto iter = UNIT_NAMES.find(unit);
  if(iter == UNIT_NAMES.end())
  {
    throw std::invalid_argument("Unknown length unit");
  }
  return iter->second;
}

double getConversionFactor(LengthUnit sourceUnits, LengthUnit targetUnits)
{
  static const std::unordered_map<LengthUnit, double, LengthUnitHash> CONVERSION_TO_CM {
    {LengthUnit::am, 1e-16},
    {LengthUnit::fm, 1e-13},
    {LengthUnit::pm, 1e-10},
    {LengthUnit::km, 1e5},
    {LengthUnit::hm, 1e4},
    {LengthUnit::dam, 1e3},
    {LengthUnit::m, 1e2},
    {LengthUnit::dm, 1e1},
    {LengthUnit::cm, 1.0},
    {LengthUnit::mm, 1e-1},
    {LengthUnit::um, 1e-4},
    {LengthUnit::nm, 1e-7},
    {LengthUnit::angstrom, 1e-8},
    {LengthUnit::miles, 2.54 * 12.0 * 5280},
    {LengthUnit::feet, 2.54 * 12.0},
    {LengthUnit::inches, 2.54},
    {LengthUnit::mils, 2.54e-3}};

  if(sourceUnits == LengthUnit::unspecified || targetUnits == LengthUnit::unspecified)
  {
    throw std::invalid_argument("Cannot convert with unspecified units");
  }

  if(sourceUnits == targetUnits)
  {
    return 1.0;
  }

  double conversionToCm = CONVERSION_TO_CM.find(sourceUnits)->second;
  if(targetUnits == LengthUnit::cm)
  {
    return conversionToCm;
  }
  return conversionToCm / CONVERSION_TO_CM.find(targetUnits)->second;
}

double convert(double sourceValue, LengthUnit sourceUnits, LengthUnit targetUnits)
{
  return sourceValue * getConversionFactor(sourceUnits, targetUnits);
}

}  // namespace utilities
}  // namespace axom
