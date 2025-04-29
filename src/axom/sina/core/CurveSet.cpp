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
#include <algorithm>

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
 * @param nameList the vector of curve names to add the curve's name to.
                   Used for tracking insertion order for codes.
 */
void addCurve(Curve &&curve, CurveSet::CurveMap &curves, std::vector<std::string> &nameList)
{
  auto &curveName = curve.getName();
  auto existing = curves.find(curveName);
  if(existing == curves.end())
  {
    curves.insert(std::make_pair(curveName, curve));
    nameList.emplace_back(curveName);
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
 * @return a struct containing the curveMap and curveOrder.
 */
CurveSet::curveNodeInfo extractCurveMap(conduit::Node const &parent, std::string const &childNodeName)
{
  CurveSet::CurveMap curveMap;
  std::vector<std::string> curveNames;

  if(!parent.has_child(childNodeName))
  {
    return CurveSet::curveNodeInfo {curveMap, curveNames};
  }

  auto &mapAsNode = parent.child(childNodeName);
  for(auto iter = mapAsNode.children(); iter.has_next();)
  {
    auto &curveAsNode = iter.next();
    std::string curveName = iter.name();
    curveNames.emplace_back(curveName);
    Curve curve {curveName, curveAsNode};
    curveMap.insert(std::make_pair(std::move(curveName), std::move(curve)));
  }

  return CurveSet::curveNodeInfo {std::move(curveMap), std::move(curveNames)};
};

/**
 * Create a Conduit node to represent the given CurveMap.
 *
 * @param curveMap the CurveMap to convert
 * @param nameList the insertion-order "index" of CurveMap
 * @param curveOrder how nameList should be sorted if not oldest-first, ex: alphabetical.
 * @return the map as a node
 */
conduit::Node createCurveMapNode(CurveSet::CurveMap const &curveMap,
                                 std::vector<std::string> const &nameList,
                                 CurveSet::CurveOrder const curveOrder)
{
  conduit::Node mapNode;
  mapNode.set_dtype(conduit::DataType::object());
  // Copy for sorting
  std::vector<std::string> orderedNameList = nameList;
  switch(curveOrder)
  {
  case CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST:
    break;
  case CurveSet::CurveOrder::REGISTRATION_NEWEST_FIRST:
    std::reverse(orderedNameList.begin(), orderedNameList.end());
    break;
  case CurveSet::CurveOrder::ALPHABETIC:
    std::sort(orderedNameList.begin(), orderedNameList.end());
    break;
  case CurveSet::CurveOrder::REVERSE_ALPHABETIC:
    std::sort(orderedNameList.begin(), orderedNameList.end(), std::greater<std::string>());
    break;
  }
  for(auto &curveName : orderedNameList)
  {
    auto expectedCurve = curveMap.find(curveName);
    // Warn if not found? Should this be allowed to happen?
    if(expectedCurve != curveMap.end())
    {
      mapNode.add_child(curveName) = expectedCurve->second.toNode();
    }
  }
  return mapNode;
}

}  // namespace

CurveSet::CurveSet(std::string name_)
  : name {std::move(name_)}
  , independentCurves {}
  , dependentCurves {}
  , independentCurveNameOrder {}
  , dependentCurveNameOrder {}
{ }

CurveSet::CurveSet(std::string name_, conduit::Node const &node)
{
  name = std::move(name_);
  auto independentCurveInfo = extractCurveMap(node, INDEPENDENT_KEY);
  independentCurves = std::move(independentCurveInfo.curveMap);
  independentCurveNameOrder = std::move(independentCurveInfo.curveOrder);
  auto dependentCurveInfo = extractCurveMap(node, DEPENDENT_KEY);
  dependentCurves = std::move(dependentCurveInfo.curveMap);
  dependentCurveNameOrder = std::move(dependentCurveInfo.curveOrder);
}

void CurveSet::addIndependentCurve(Curve curve)
{
  addCurve(std::move(curve), independentCurves, independentCurveNameOrder);
}

void CurveSet::addDependentCurve(Curve curve)
{
  addCurve(std::move(curve), dependentCurves, dependentCurveNameOrder);
}

conduit::Node CurveSet::toNode(CurveOrder curveOrder) const
{
  conduit::Node asNode;
  asNode[INDEPENDENT_KEY] =
    createCurveMapNode(independentCurves, independentCurveNameOrder, curveOrder);
  asNode[DEPENDENT_KEY] = createCurveMapNode(dependentCurves, dependentCurveNameOrder, curveOrder);
  return asNode;
}

}  // namespace sina
}  // namespace axom
