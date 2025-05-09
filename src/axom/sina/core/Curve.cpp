// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Curve.cpp
 *
 * \brief   Implementation file for Sina Curve class
 *
 ******************************************************************************
 */

#include "axom/sina/core/Curve.hpp"
#include "axom/sina/core/ConduitUtil.hpp"

#include <memory>

namespace axom
{
namespace sina
{

namespace
{
constexpr auto CURVE_TYPE_NAME = "curve";
constexpr auto VALUES_KEY = "value";
constexpr auto UNITS_KEY = "units";
constexpr auto TAGS_KEY = "tags";
}  // namespace

Curve::Curve(std::string name_, std::vector<double> values_)
  : name {std::move(name_)}
  , values {std::move(values_)}
  , units {}
  , tags {}
{ }

Curve::Curve(std::string name_, double const *values_, std::size_t numValues)
  : name {std::move(name_)}
  , values {values_, values_ + numValues}
  , units {}
  , tags {}
{ }

Curve::Curve(std::string name_, conduit::Node const &curveAsNode)
  : name {std::move(name_)}
  , values {}
  , units {}
  , tags {}
{
  auto &valuesAsNode = getRequiredField(VALUES_KEY, curveAsNode, CURVE_TYPE_NAME);
  values = toDoubleVector(valuesAsNode, VALUES_KEY);

  units = getOptionalString(UNITS_KEY, curveAsNode, CURVE_TYPE_NAME);

  if(curveAsNode.has_child(TAGS_KEY))
  {
    tags = toStringVector(curveAsNode[TAGS_KEY], TAGS_KEY);
  }
}

void Curve::setUnits(std::string units_) { units = std::move(units_); }

void Curve::setTags(std::vector<std::string> tags_) { tags = std::move(tags_); }

conduit::Node Curve::toNode() const
{
  conduit::Node asNode;
  asNode[VALUES_KEY] = values;
  if(!units.empty())
  {
    asNode[UNITS_KEY] = units;
  }
  if(!tags.empty())
  {
    addStringsToNode(asNode, TAGS_KEY, tags);
  }
  return asNode;
}

}  // namespace sina
}  // namespace axom
