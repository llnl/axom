// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file ConduitUtil.cpp
 *
 * \brief   Implementation file for Sina Conduit utility functions
 *
 ******************************************************************************
 */

#include "axom/sina/core/ConduitUtil.hpp"

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include "conduit.hpp"

namespace axom
{
namespace sina
{

namespace
{
/**
 * Get the given field as string. If it is not a string, an exception is
 * thrown with a user-friendly message.
 *
 * @param field the value of the field
 * @param fieldName the name of the field
 * @param parentType the name of the parent which contained the field
 * @return the avlue of the field
 * @throws std::invalid_argument if the field is not a string
 */
std::string getExpectedString(conduit::Node const &field,
                              std::string const &fieldName,
                              std::string const &parentType)
{
  if(!field.dtype().is_string())
  {
    std::ostringstream message;
    message << "The field '" << fieldName << "' for objects of type '" << parentType
            << "' must be a string, was '" << field.dtype().name() << "'";
    throw std::invalid_argument(message.str());
  }
  return field.as_string();
}
}  // namespace

conduit::Node const &getRequiredField(std::string const &fieldName,
                                      conduit::Node const &parent,
                                      std::string const &parentType)
{
  if(!parent.has_child(fieldName))
  {
    std::ostringstream message;
    message << "The field '" << fieldName << "' is required for objects of type '" << parentType
            << "'";
    throw std::invalid_argument(message.str());
  }
  return parent.child(fieldName);
}

std::string getRequiredString(std::string const &fieldName,
                              conduit::Node const &parent,
                              std::string const &parentType)
{
  conduit::Node const &field = getRequiredField(fieldName, parent, parentType);
  return getExpectedString(field, fieldName, parentType);
}

std::string getOptionalString(std::string const &fieldName,
                              conduit::Node const &parent,
                              std::string const &parentType)
{
  if(!parent.has_child(fieldName) || parent.child(fieldName).dtype().is_empty())
  {
    return "";
  }
  return getExpectedString(parent.child(fieldName), fieldName, parentType);
}

double getRequiredDouble(std::string const &fieldName,
                         conduit::Node const &parent,
                         std::string const &parentType)
{
  auto &ref = getRequiredField(fieldName, parent, parentType);
  if(!ref.dtype().is_number())
  {
    std::ostringstream message;
    message << "The field '" << fieldName << "' for objects of type '" << parentType
            << "' must be a double";
    throw std::invalid_argument(message.str());
  }
  return ref.as_double();
}

void addStringsToNode(conduit::Node &parent,
                      std::string const &child_name,
                      std::vector<std::string> const &string_values)
{
  // If the child already exists, add_child returns it
  conduit::Node &child_node = parent.add_child(child_name);
  for(auto &value : string_values)
  {
    auto &list_entry = child_node.append();
    list_entry.set(value);
  }

  // If there were no children, this will be a null node rather than a list.
  // We prefer empty lists to null values in our serialized JSON, so force
  // this to be a list
  if(child_node.number_of_children() == 0)
  {
    child_node.set_dtype(conduit::DataType::list());
  }
}

std::vector<double> toDoubleVector(conduit::Node const &node, std::string const &name)
{
  if(node.dtype().is_list() && node.dtype().number_of_elements() == 0)
  {
    return std::vector<double> {};
  }
  conduit::Node asDoubles;
  try
  {
    node.to_double_array(asDoubles);
  }
  catch(conduit::Error const &err)
  {
    std::ostringstream errStream;
    errStream << "Error trying to convert node \"" << name << "\" into a list of doubles"
              << err.what();
    throw std::invalid_argument(errStream.str());
  }
  double const *start = asDoubles.as_double_ptr();
  auto count = static_cast<std::vector<double>::size_type>(asDoubles.dtype().number_of_elements());
  return std::vector<double> {start, start + count};
}

std::vector<std::string> toStringVector(conduit::Node const &node, std::string const &name)
{
  std::vector<std::string> converted;
  if(!node.dtype().is_list())
  {
    std::ostringstream errStream;
    errStream << "Error trying to convert node \"" << name
              << "\" into a list of strings. It is not a list. " << node.to_json_default();
    throw std::invalid_argument(errStream.str());
  }
  for(auto iter = node.children(); iter.has_next();)
  {
    auto &child = iter.next();
    if(child.dtype().is_string())
    {
      converted.emplace_back(child.as_string());
    }
    else
    {
      std::ostringstream errStream;
      errStream << "Error trying to convert node \"" << name
                << "\" into a list of strings. A value is not a string. " << node.to_json_default();
      throw std::invalid_argument(errStream.str());
    }
  }
  return converted;
}

}  // namespace sina
}  // namespace axom
