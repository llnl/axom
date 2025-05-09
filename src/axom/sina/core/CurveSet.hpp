// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_CURVESET_HPP
#define SINA_CURVESET_HPP

/*!
 ******************************************************************************
 *
 * \file CurveSet.hpp
 *
 * \brief   Header file for Sina CurveSet class
 *
 * \sa Curve.hpp
 *
 ******************************************************************************
 */

#include <string>
#include <unordered_map>

#include "conduit.hpp"

#include "axom/sina/core/Curve.hpp"

namespace axom
{
namespace sina
{

/**
 * \brief A CurveSet represents an entry in a record's "curve_set".
 *
 * A CurveSet consist of a set of independent and dependent curves. Each curve
 * is a list of numbers along with optional units and tags.
 *
 * \sa Record
 * \sa Curve
 */
class CurveSet
{
public:
  /**
     * An unordered map of Curve objects.
     */
  using CurveMap = std::unordered_map<std::string, Curve>;

  /**
     * \brief Create a CurveSet with the given name
     *
     * \param name the name of the CurveSet
     */
  explicit CurveSet(std::string name);

  /**
     * \brief Create a CurveSet from the given Conduit node.
     *
     * \param name the name of the CurveSet
     * \param node the Conduit node representing the CurveSet
     */
  CurveSet(std::string name, conduit::Node const &node);

  /**
     * \brief Get the name of the this CurveSet.
     *
     * \return the curve set's name
     */
  std::string const &getName() const { return name; }

  /**
     * \brief Add an independent curve.
     *
     * \param curve the curve to add
     */
  void addIndependentCurve(Curve curve);

  /**
     * \brief Add a dependent curve.
     *
     * \param curve the curve to add
     */
  void addDependentCurve(Curve curve);

  /**
     * \brief Get a map of all the independent curves.
     *
     * \return a map of all the independent curves
     */
  CurveMap const &getIndependentCurves() const { return independentCurves; }

  /**
     * \brief Get a map of all the dependent curves.
     *
     * \return a map of all the dependent curves
     */
  CurveMap const &getDependentCurves() const { return dependentCurves; }

  /**
     * \brief Convert his CurveSet to a Conduit node.
     *
     * \return the Node representation of this CurveSet
     */
  conduit::Node toNode() const;

private:
  std::string name;
  CurveMap independentCurves;
  CurveMap dependentCurves;
};

}  // namespace sina
}  // namespace axom

#endif  //SINA_CURVESET_HPP
