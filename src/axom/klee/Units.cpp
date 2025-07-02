// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Proxy.hpp"

#include "axom/klee/Units.hpp"
#include "axom/klee/KleeError.hpp"

#include <stdexcept>
#include <unordered_map>

namespace axom
{
namespace klee
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
}  // namespace

LengthUnit parseLengthUnits(const std::string &unitsAsString, const std::string &path)
{
  static const std::unordered_map<std::string, LengthUnit> UNITS_BY_NAME {
    {"km", LengthUnit::km},
    {"m", LengthUnit::m},
    {"dm", LengthUnit::dm},
    {"cm", LengthUnit::cm},
    {"mm", LengthUnit::mm},
    {"um", LengthUnit::um},
    {"nm", LengthUnit::nm},
    {"A", LengthUnit::angstrom},
    {"miles", LengthUnit::miles},
    {"ft", LengthUnit::feet},
    {"feet", LengthUnit::feet},
    {"in", LengthUnit::inches},
    {"inches", LengthUnit::inches},
    {"mils", LengthUnit::mils}};

  auto iter = UNITS_BY_NAME.find(unitsAsString);

  if(iter == UNITS_BY_NAME.end())
  {
    std::string message = "Unrecognized units: ";
    message += unitsAsString;
    throw KleeError({path, message});
  }

  return iter->second;
}

LengthUnit parseLengthUnits(const inlet::Proxy &unitsAsProxy)
{
  return parseLengthUnits(unitsAsProxy.get<std::string>(), unitsAsProxy.name());
}

double getConversionFactor(LengthUnit sourceUnits, LengthUnit targetUnits)
{
  static const std::unordered_map<LengthUnit, double, LengthUnitHash> CONVERSION_TO_CM {
    {LengthUnit::km, 1e5},
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

}  // namespace klee
}  // namespace axom