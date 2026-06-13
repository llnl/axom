// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "shaping_helpers.hpp"

namespace axom
{
namespace quest
{
namespace shaping
{
namespace
{
constexpr const char* SHAPE_INOUT_PREFIX = "inout_";
constexpr const char* MATERIAL_INOUT_PREFIX = "mat_inout_";
constexpr const char* VOLUME_FRACTION_PREFIX = "vol_frac_";
constexpr const char* SHAPE_VOLUME_FRACTION_PREFIX = "shape_vol_frac_";

std::string extractSuffixedName(const std::string& fieldName, const std::string& prefix)
{
  return axom::utilities::string::startsWith(fieldName, prefix) ? fieldName.substr(prefix.size())
                                                                : std::string {};
}
}  // namespace

std::string shapeInOutFieldName(const std::string& shapeName)
{
  return axom::fmt::format("{}{}", SHAPE_INOUT_PREFIX, shapeName);
}

std::string materialInOutFieldName(const std::string& materialName)
{
  return axom::fmt::format("{}{}", MATERIAL_INOUT_PREFIX, materialName);
}

std::string volumeFractionFieldName(const std::string& materialName)
{
  return axom::fmt::format("{}{}", VOLUME_FRACTION_PREFIX, materialName);
}

std::string shapeVolumeFractionFieldName(const std::string& shapeName)
{
  return axom::fmt::format("{}{}", SHAPE_VOLUME_FRACTION_PREFIX, shapeName);
}

std::string materialNameFromMaterialInOutFieldName(const std::string& fieldName)
{
  return extractSuffixedName(fieldName, MATERIAL_INOUT_PREFIX);
}

std::string materialNameFromVolumeFractionFieldName(const std::string& fieldName)
{
  return extractSuffixedName(fieldName, VOLUME_FRACTION_PREFIX);
}

}  // namespace shaping
}  // namespace quest
}  // namespace axom
