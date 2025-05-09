// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_JSONUTIL_HPP
#define SINA_JSONUTIL_HPP

/*!
 ******************************************************************************
 *
 * \file ConduitUtil.hpp
 *
 * \brief   Header file for Sina Conduit utility functions
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
 * \brief Get a required field from a conduit Node.
 *
 * \param fieldName the name of the field to get
 * \param parent the parent object from which to get the field
 * \param parentType a user-friendly name of the type of the parent to use
 *                   in an error message if the field doesn't exist.
 * \return the requested field as a Node
 * \throws std::invalid_argument if the field does not exist
 */
conduit::Node const &getRequiredField(std::string const &fieldName,
                                      conduit::Node const &parent,
                                      std::string const &parentType);

/**
 * \brief Get the value of a required field from a conduit Node. The field value
 *        must be a string.
 *
 * \param fieldName the name of the field to get
 * \param parent the parent object from which to get the field
 * \param parentType a user-friendly name of the type of the parent to use
 *                   in an error message if the field doesn't exist.
 * \return the value of the requested field
 * \throws std::invalid_argument if the field does not exist or is not a string
 */
std::string getRequiredString(std::string const &fieldName,
                              conduit::Node const &parent,
                              std::string const &parentType);

/**
 * \brief Get the value of a required field from a conduit Node. The field value
 *        must be a double.
 *
 * \param fieldName the name of the field to get
 * \param parent the parent object from which to get the field
 * \param parentType a user-friendly name of the type of the parent to use
 *                   in an error message if the field doesn't exist.
 * \return the value of the requested field
 * \throws std::invalid_argument if the field does not exist or is not a double
 */
double getRequiredDouble(std::string const &fieldName,
                         conduit::Node const &parent,
                         std::string const &parentType);

/**
 * \brief Get the value of an optional field from a conduit Node. The field value
 *        must be a string if it is present.
 *
 * \param fieldName the name of the field to get
 * \param parent the parent object from which to get the field
 * \param parentType a user-friendly name of the type of the parent to use
 *                   in an error message if the field doesn't exist.
 * \return the value of the requested field, or an empty string if it
 *         does not exist
 * \throws std::invalid_argument if the field exists but is not a string
 */
std::string getOptionalString(std::string const &fieldName,
                              conduit::Node const &parent,
                              std::string const &parentType);

/**
 * \brief Convert the given node to a vector of doubles.
 *
 * \param node the node to convert
 * \param name the name of the node, used in error reporting
 * \return the node as a list of doubles
 * \throws std::invalid_argument if the node is not a list of doubles
 */
std::vector<double> toDoubleVector(conduit::Node const &node, std::string const &name);

/**
 * \brief Convert the given node to a vector of strings.
 *
 * \param node the node to convert
 * \param name the name of the node, used in error reporting
 * \return the node as a list of strings
 * \throws std::invalid_argument if the node is not a list of strings
 */
std::vector<std::string> toStringVector(conduit::Node const &node, std::string const &name);

/**
 * \brief Add a vector of strings to a Node. This operation's not natively
 *        part of Conduit.
 *
 * \param parent the node to add the strings to
 * \param child_name the name of the child (aka the name of the field)
 * \param string_values the data values for the field
 */
void addStringsToNode(conduit::Node &parent,
                      const std::string &child_name,
                      std::vector<std::string> const &string_values);

}  // end namespace sina
}  // end namespace axom

#endif  //SINA_JSONUTIL_HPP
