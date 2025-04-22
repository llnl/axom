// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file CurveSet.cpp
 *
 * \brief   Implementation file for Sina CurveSet class
 *
 * \sa Curve.hpp
 *
 ******************************************************************************
 */

#include "axom/sina/core/CurveSet.hpp"

#include <utility>

#include "axom/sina/core/ConduitUtil.hpp"

namespace axom
{
namespace sina
{

namespace
{

constexpr auto INDEPENDENT_KEY = "independent";
constexpr auto DEPENDENT_KEY = "dependent";

/**
 * Add a curve to the given curve map.
 *
 * @param curve the curve to add
 * @param curves the CurveMap to which to add the curve
 */
void addCurve(Curve &&curve, CurveSet::CurveMap &curves)
{
  auto &curveName = curve.getName();
  auto existing = curves.find(curveName);
  if(existing == curves.end())
  {
    curves.insert(std::make_pair(curveName, curve));
  }
  else
  {
    existing->second = curve;
  }
}

/**
 * Extract a CurveMap from the given node.
 *
 * @param parent the parent node
 * @param childNodeName the name of the child node
 * @return a CurveMap representing the specified child
 */
CurveSet::CurveMap extractCurveMap(conduit::Node const &parent, std::string const &childNodeName)
{
  CurveSet::CurveMap curveMap;
  if(!parent.has_child(childNodeName))
  {
    return curveMap;
  }

  auto &mapAsNode = parent.child(childNodeName);
  for(auto iter = mapAsNode.children(); iter.has_next();)
  {
    auto &curveAsNode = iter.next();
    std::string curveName = iter.name();
    Curve curve {curveName, curveAsNode};
    curveMap.insert(std::make_pair(std::move(curveName), std::move(curve)));
  }

  return curveMap;
};

/**
 * Create a Conduit node to represent the given CurveMap.
 *
 * @param curveMap the CurveMap to convert
 * @return the map as a node
 */
conduit::Node createCurveMapNode(CurveSet::CurveMap const &curveMap)
{
  conduit::Node mapNode;
  mapNode.set_dtype(conduit::DataType::object());
  for(auto &entry : curveMap)
  {
    mapNode.add_child(entry.first) = entry.second.toNode();
  }
  return mapNode;
}

}  // namespace

CurveSet::CurveSet(std::string name_)
  : name {std::move(name_)}
  , independentCurves {}
  , dependentCurves {}
{ }

CurveSet::CurveSet(std::string name_, conduit::Node const &node)
  : name {std::move(name_)}
  , independentCurves {extractCurveMap(node, INDEPENDENT_KEY)}
  , dependentCurves {extractCurveMap(node, DEPENDENT_KEY)}
{ }

void CurveSet::addIndependentCurve(Curve curve) { addCurve(std::move(curve), independentCurves); }

void CurveSet::addDependentCurve(Curve curve) { addCurve(std::move(curve), dependentCurves); }

conduit::Node CurveSet::toNode() const
{
  conduit::Node asNode;
  asNode[INDEPENDENT_KEY] = createCurveMapNode(independentCurves);
  asNode[DEPENDENT_KEY] = createCurveMapNode(dependentCurves);
  return asNode;
}

}  // namespace sina
}  // namespace axom
