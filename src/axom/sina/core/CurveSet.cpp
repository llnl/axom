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
#include <set>

#include "axom/sina/core/ConduitUtil.hpp"

namespace axom
{
namespace sina
{

// Definition with initialization
CurveSet::CurveOrder sinaDefaultCurveOrder = CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST;

// Function to set the default curve order
void setDefaultCurveOrder(CurveSet::CurveOrder order) { sinaDefaultCurveOrder = order; }

namespace
{

constexpr auto INDEPENDENT_KEY = "independent";
constexpr auto DEPENDENT_KEY = "dependent";

/**
 * Reset the default sinaCurveOrder, all records created after this will us e that
 */
/**
 * Add a curve to the given curve map.
 *
 * Helper function, users use addIndependentCurve/addDependentCurve.
 *
 * @param curve the curve to add
 * @param curves the CurveMap to which to add the curve
 * @param nameList the vector of curve names to add the curve's name to.
                   Used for tracking insertion order for codes.
 */
void addCurve(Curve &&curve, CurveSet::CurveMap &curves, std::vector<std::string> &nameList)
{
  std::string curveName = curve.getName();  // Make a COPY before moving
  auto existing = curves.find(curveName);
  if(existing == curves.end())
  {
    curves.insert(std::make_pair(curveName, std::move(curve)));  // Explicit move
    nameList.emplace_back(std::move(curveName));                 // Move the copy into the list
  }
  else
  {
    existing->second = std::move(curve);  // Explicit move
  }
}

/**
 * Set a list of curve names as the canonical insertion order.
 *
 * Helper, users use independent/dependent specific ones as above.
 */
bool applyCustomCurveOrder(const std::vector<std::string> &newOrder,
                           std::vector<std::string> &oldOrder)
{
  if(newOrder.size() != oldOrder.size())
  {
    return false;
  }
  std::set<std::string> newOrderSet(newOrder.begin(), newOrder.end());
  std::set<std::string> oldOrderSet(oldOrder.begin(), oldOrder.end());
  if(newOrderSet == oldOrderSet)
  {
    oldOrder = newOrder;
    return true;
  };
  return false;
}

/**
 * Extract a CurveMap from the given node.
 *
 * @param parent the parent node
 * @param childNodeName the name of the child node
 * @return a struct containing the curveMap and ordered curve names.
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
  , orderedIndependentCurveNames {}
  , orderedDependentCurveNames {}
{ }

CurveSet::CurveSet(std::string name_, conduit::Node const &node)
{
  name = std::move(name_);
  auto independentCurveInfo = extractCurveMap(node, INDEPENDENT_KEY);
  independentCurves = std::move(independentCurveInfo.curveMap);
  orderedIndependentCurveNames = std::move(independentCurveInfo.orderedCurveNames);
  auto dependentCurveInfo = extractCurveMap(node, DEPENDENT_KEY);
  dependentCurves = std::move(dependentCurveInfo.curveMap);
  orderedDependentCurveNames = std::move(dependentCurveInfo.orderedCurveNames);
}

void CurveSet::addIndependentCurve(Curve curve)
{
  addCurve(std::move(curve), independentCurves, orderedIndependentCurveNames);
}

void CurveSet::addDependentCurve(Curve curve)
{
  addCurve(std::move(curve), dependentCurves, orderedDependentCurveNames);
}

bool CurveSet::applyCustomIndependentCurveOrder(const std::vector<std::string> newOrderedCurveNames)
{
  return applyCustomCurveOrder(newOrderedCurveNames, orderedIndependentCurveNames);
}

bool CurveSet::applyCustomDependentCurveOrder(const std::vector<std::string> newOrderedCurveNames)
{
  return applyCustomCurveOrder(newOrderedCurveNames, orderedDependentCurveNames);
}

conduit::Node CurveSet::toNode(CurveOrder curveOrder) const
{
  conduit::Node asNode;
  asNode[INDEPENDENT_KEY] =
    createCurveMapNode(independentCurves, orderedIndependentCurveNames, curveOrder);
  asNode[DEPENDENT_KEY] = createCurveMapNode(dependentCurves, orderedDependentCurveNames, curveOrder);
  return asNode;
}

conduit::Node CurveSet::toNode() const { return toNode(sinaDefaultCurveOrder); }

}  // namespace sina
}  // namespace axom
