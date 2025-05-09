// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_CURVE_HPP
#define SINA_CURVE_HPP

/*!
 ******************************************************************************
 *
 * \file Curve.hpp
 *
 * \brief   Header file for Sina Curve class
 *
 ******************************************************************************
 */

#include <string>
#include <vector>

#include "conduit.hpp"

namespace axom
{
namespace sina
{

/**
 * A Curve represents a 1-dimensional curve inside a CurveSet.
 */
class Curve
{
public:
  /**
     * \brief Create a Curve with the given name and values
     *
     * \param name the name of the curve
     * \param values the curve's values
     */
  Curve(std::string name, std::vector<double> values);

  /**
     * \brief Create a Curve with the given name and values
     *
     * \param name the name of the curve
     * \param values the curve's values
     * \param numValues the number of values.
     */
  Curve(std::string name, double const *values, std::size_t numValues);

  /**
     * \brief Create a Curve by deserializing a conduit node.
     *
     * \param name the name of the curve
     * \param curveAsNode the serialized version of a curve
     */
  Curve(std::string name, conduit::Node const &curveAsNode);

  /**
     * \brief Get the curve's name.
     *
     * \return the curve's name
     */
  std::string const &getName() const { return name; }

  /**
     * \brief Get the values of the curve.
     *
     * \return the curve's values
     */
  std::vector<double> const &getValues() const { return values; }

  /**
     * \brief Set the units of the values.
     *
     * \param units the value's units
     */
  void setUnits(std::string units);

  /**
     * \brief Get the units of the values.
     *
     * \return the value's units
     */
  std::string const &getUnits() const { return units; }

  /**
     * \brief Set the tags for this curve.
     *
     * \param tags the curve's tags
     */
  void setTags(std::vector<std::string> tags);

  /**
     * \brief Get the tags for this curve.
     *
     * \return the curve's tags
     */
  std::vector<std::string> const &getTags() const { return tags; }

  /**
     * \brief Convert this curve to a Conduit node.
     *
     * \return a Conduit representation of this curve
     */
  conduit::Node toNode() const;

private:
  std::string name;
  std::vector<double> values;
  std::string units;
  std::vector<std::string> tags;
};

}  // namespace sina
}  // namespace axom

#endif  //SINA_CURVE_HPP
