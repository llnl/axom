// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Datum.cpp
 *
 * \brief   Implementation file for Sina Datum class
 *
 ******************************************************************************
 */

#include "axom/sina/core/Datum.hpp"
#include "conduit.hpp"

#include <utility>
#include <sstream>
#include <stdexcept>

namespace
{

char const VALUE_FIELD[] = "value";
char const UNITS_FIELD[] = "units";
char const TAGS_FIELD[] = "tags";
char const DATA_PARENT_TYPE[] = "data";

}  // namespace

namespace axom
{
namespace sina
{

Datum::Datum(const std::string &value_) : stringValue {value_}
{
  //Set type to String, as we know it uses strings
  type = ValueType::String;
}

Datum::Datum(const double &value_) : scalarValue {value_}
{
  //Set type to Scalar, as we know it uses doubles
  type = ValueType::Scalar;
}

Datum::Datum(const std::vector<std::string> &value_) : stringArrayValue {value_}
{
  //Set type to StringArray, as we know it uses an array of strings
  type = ValueType::StringArray;
}

Datum::Datum(const std::vector<double> &value_) : scalarArrayValue {value_}
{
  //Set type to ScalarArray, as we know it uses an array of doubles
  type = ValueType::ScalarArray;
}

Datum::Datum(conduit::Node const &asNode)
{
  //Need to determine what type of Datum we have: Scalar (double), String,
  //or list of one of those two.
  conduit::Node valueNode = getRequiredField(VALUE_FIELD, asNode, DATA_PARENT_TYPE);
  if(valueNode.dtype().is_string())
  {
    stringValue = valueNode.as_string();
    type = ValueType::String;
  }
  else if(valueNode.dtype().is_number() && valueNode.dtype().number_of_elements() == 1)
  {
    scalarValue = valueNode.to_double();
    type = ValueType::Scalar;
  }
  // There are two different ways to end up with a "list" of numbers in conduit, but
  // only one of them tests True for is_list. This handles the other.
  else if(valueNode.dtype().is_number())
  {
    type = ValueType::ScalarArray;
    // What's passed in could be an array of any numeric type
    // We pass a cast copy into captureNode
    conduit::Node captureNode;
    valueNode.to_float64_array(captureNode);
    std::vector<double> array_as_vect(
      captureNode.as_double_ptr(),
      captureNode.as_double_ptr() + captureNode.dtype().number_of_elements());
    scalarArrayValue = array_as_vect;
  }
  else if(valueNode.dtype().is_list())
  {
    //An empty list is assumed to be an empty list of doubles.
    //This only works because this field is immutable!
    //If this ever changes, or if Datum's type is used directly to make
    //decisions (ex: Sina deciding where to store data), this logic
    //should be revisited.
    if(valueNode.number_of_children() == 0 || valueNode[0].dtype().is_number())
    {
      type = ValueType::ScalarArray;
    }
    else if(valueNode[0].dtype().is_string())
    {
      type = ValueType::StringArray;
    }
    else
    {
      std::ostringstream message;
      message << "The only valid types for an array '" << VALUE_FIELD
              << "' are strings and numbers. Got '" << valueNode.to_json() << "'";
      throw std::invalid_argument(message.str());
    }

    auto itr = valueNode.children();
    while(itr.has_next())
    {
      conduit::Node const &entry = itr.next();
      if(entry.dtype().is_string() && type == ValueType::StringArray)
      {
        stringArrayValue.emplace_back(entry.as_string());
      }
      else if(entry.dtype().is_number() && type == ValueType::ScalarArray)
      {
        scalarArrayValue.emplace_back(entry.to_double());
      }
      else
      {
        std::ostringstream message;
        message << "If the required field '" << VALUE_FIELD
                << "' is an array, it must consist of only strings or only numbers, "
                << "but got '" << entry.dtype().name() << "' (" << entry.to_json() << ")";
        throw std::invalid_argument(message.str());
      }
    }
  }
  else
  {
    std::ostringstream message;
    message << "The required field '" << VALUE_FIELD
            << "' must be a string, double, list of strings, or list of doubles.";
    throw std::invalid_argument(message.str());
  }

  //Get the units, if there are any
  units = getOptionalString(UNITS_FIELD, asNode, DATA_PARENT_TYPE);

  //Need to grab the tags and add them to a vector of strings
  if(asNode.has_child(TAGS_FIELD))
  {
    auto tagNodeIter = asNode[TAGS_FIELD].children();
    while(tagNodeIter.has_next())
    {
      auto &tag = tagNodeIter.next();
      if(tag.dtype().is_string())
      {
        tags.emplace_back(std::string(tag.as_string()));
      }
      else
      {
        std::ostringstream message;
        message << "The optional field '" << TAGS_FIELD << "' must be an array of strings. Found '"
                << tag.dtype().name() << "' instead.";
        throw std::invalid_argument(message.str());
      }
    }
  }
}

void Datum::setUnits(const std::string &units_) { units = units_; }

void Datum::setTags(const std::vector<std::string> &tags_) { tags = tags_; }

conduit::Node Datum::toNode() const
{
  conduit::Node asNode;
  switch(type)
  {
  case ValueType::Scalar:
    asNode[VALUE_FIELD] = scalarValue;
    break;
  case ValueType::String:
    asNode[VALUE_FIELD] = stringValue;
    break;
  case ValueType::ScalarArray:
    asNode[VALUE_FIELD] = scalarArrayValue;
    break;
  case ValueType::StringArray:
    addStringsToNode(asNode, VALUE_FIELD, stringArrayValue);
    break;
  }
  if(tags.size() > 0) addStringsToNode(asNode, TAGS_FIELD, tags);
  if(!units.empty()) asNode[UNITS_FIELD] = units;
  return asNode;
};

}  // namespace sina
}  // namespace axom
