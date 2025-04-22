// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file JSONSchemaWriter.cpp
 *
 * \brief This file contains the class implementation of the JSONSchemaWriter.
 *******************************************************************************
 */

#include "axom/inlet/JSONSchemaWriter.hpp"

#include <iostream>

#include "axom/slic.hpp"
#include "axom/inlet/Container.hpp"

namespace axom
{
namespace inlet
{
namespace detail
{
/*!
 *******************************************************************************
 * \brief Writes range information (maximum/minimum) to a node corresponding
 * to the properties of a field
 * 
 * \param [in] view The view containing a size-2 array with max/min data
 * \param [inout] node The node to write to
 * 
 * \return The type of the min/max values
 *******************************************************************************
 */
sidre::TypeID recordRange(const sidre::View* view, conduit::Node& node)
{
  auto type = view->getTypeID();
  if(type == sidre::INT_ID)
  {
    const int* range = view->getData();
    node["minimum"] = range[0];
    node["maximum"] = range[1];
  }
  else
  {
    const double* range = view->getData();
    node["minimum"] = range[0];
    node["maximum"] = range[1];
  }

  return type;
}

/*!
 *******************************************************************************
 * \brief Writes information on the set of valid values to a node corresponding
 * to the properties of a field
 * 
 * \param [in] view The view containing the array of valid numerical values
 * \param [inout] node The node to write to
 * 
 * \return The type of the valid values
 *******************************************************************************
 */
sidre::TypeID recordEnum(const sidre::View* view, conduit::Node& node)
{
  sidre::TypeID type = view->getTypeID();
  auto size = view->getNumElements();
  if(type == sidre::INT_ID)
  {
    const int* values = view->getData();
    for(int i = 0; i < size; i++)
    {
      node["enum"].append() = values[i];
    }
  }
  else
  {
    const double* values = view->getData();
    for(int i = 0; i < size; i++)
    {
      node["enum"].append() = values[i];
    }
  }
  return type;
}
/*!
 *******************************************************************************
 * \overload
 * 
 * \param [in] group The group containing the array of valid string values
 * 
 * \note This version is intended for string enumerations, which must be stored
 * in a Group instead of a view
 *******************************************************************************
 */
sidre::TypeID recordEnum(const sidre::Group* group, conduit::Node& node)
{
  for(auto& view : group->views())
  {
    node["enum"].append() = std::string(view.getString());
  }
  return sidre::CHAR8_STR_ID;
}

/*!
 *******************************************************************************
 * \brief Writes information on a Field object to a node corresponding
 * to the properties in the schema
 * 
 * \param [in] field The field whose JSON schema is to be extracted
 * \param [inout] node The node to write to
 *******************************************************************************
 */
void recordFieldSchema(const Field& field, conduit::Node& node)
{
  const static auto jsonTypeNames = []() {
    std::unordered_map<std::underlying_type<sidre::TypeID>::type, std::string> result;
    result[sidre::TypeID::INT_ID] = "integer";
    result[sidre::TypeID::CHAR8_STR_ID] = "string";
    result[sidre::TypeID::FLOAT64_ID] = "number";
    result[sidre::TypeID::INT8_ID] = "boolean";
    return result;
  }();

  // The type information isn't guaranteed to be anywhere, so we try the
  // value, the valid values, the range, etc
  axom::sidre::TypeID type = sidre::NO_TYPE_ID;

  const auto sidreGroup = field.sidreGroup();

  if(sidreGroup->hasView("description"))
  {
    node["description"] = std::string(sidreGroup->getView("description")->getString());
  }

  if(sidreGroup->hasView("range"))
  {
    type = recordRange(sidreGroup->getView("range"), node);
  }
  else if(sidreGroup->hasView("validValues"))
  {
    type = recordEnum(sidreGroup->getView("validValues"), node);
  }
  else if(sidreGroup->hasGroup("validStringValues"))
  {
    type = recordEnum(sidreGroup->getGroup("validStringValues"), node);
  }

  if(sidreGroup->hasView("defaultValue"))
  {
    auto default_view = sidreGroup->getView("defaultValue");
    type = default_view->getTypeID();
    switch(type)
    {
    case sidre::INT_ID:
      node["default"] = static_cast<int>(default_view->getData());
      break;
    case sidre::FLOAT64_ID:
      node["default"] = static_cast<double>(default_view->getData());
      break;
    case sidre::CHAR8_STR_ID:
      node["default"] = std::string(default_view->getString());
      break;
    case sidre::INT8_ID:
      node["default"] = static_cast<std::int8_t>(default_view->getData());
      break;
    default:
      break;
    }
  }

  if(sidreGroup->hasView("value"))
  {
    type = sidreGroup->getView("value")->getTypeID();
  }

  // If we couldn't figure out what type it is, we can't enforce it in the schema
  if(type != sidre::NO_TYPE_ID)
  {
    node["type"] = jsonTypeNames.at(type);
  }
}

/*!
 *******************************************************************************
 * \brief Augments the paths in the input file corresponding to collections
 * with a label (via appending just after the path)
 * 
 * \param [inout] target The string to modify
 * \param [in] collectionPaths The paths of collections to look for
 * \param [in] label The label to append to collection paths found in the target
 *******************************************************************************
 */
void augmentCollectionPaths(std::string& target,
                            const std::vector<std::string>& collectionPaths,
                            const std::string& label)
{
  for(const auto& path : collectionPaths)
  {
    // Needs to be at the beginning - since the path is absolute
    if(target.find(path) == 0)
    {
      // Remove the name of the selected element - get the slash following the name
      const auto pos = target.find('/', path.length() + 1);
      if(pos != std::string::npos)
      {
        target.erase(path.length() + 1, pos - path.length());
      }
      // Insert the extension just after it
      target.insert(path.length(), "/" + label);
    }
  }
}

/*!
 *******************************************************************************
 * \brief Recursively removes empty nodes from a conduit::Node - used as null
 * values cannot be present in the final JSON schema
 * 
 * \param [inout] parent The node to prune from
 *******************************************************************************
 */
void pruneEmptyNodes(conduit::Node& parent)
{
  // Can't use iterators here due to iterator invalidation by deletion,
  // node.remove() also does not return an iterator like the STL container erase()
  conduit::index_t idx = 0;
  while(idx < parent.number_of_children())
  {
    auto& child = parent[idx];
    if(child.dtype().is_empty())
    {
      parent.remove(idx);
    }
    else
    {
      pruneEmptyNodes(child);
      idx++;
    }
  }
}

}  // namespace detail

JSONSchemaWriter::JSONSchemaWriter(const std::string& fileName) : m_fileName(fileName)
{
  m_schemaRoot["$schema"] = "https://json-schema.org/draft/2020-12/schema#";
  m_schemaRoot["title"] = "Input File Options";
  m_schemaRoot["description"] = "The schema for the input file";
  m_schemaRoot["type"] = "object";
}

void JSONSchemaWriter::documentContainer(const Container& container)
{
  const auto sidreGroup = container.sidreGroup();

  // The "level" of the schema to add properties to - adjusted "up" if we're an element
  // of a collection
  const sidre::Group* propertyGroup = sidreGroup;
  if(isCollectionGroup(sidreGroup->getParent()->getName()))
  {
    propertyGroup = sidreGroup->getParent();
  }
  // FIXME: Can we use in-progress path class for some of this logic?
  // JBE: I think this produces incorrect schema when the Inlet root is not the Sidre root,
  // but using the Path::dirname/parent functions will fix this
  std::string filteredPathName = propertyGroup->getPathName();
  // The location that the schema for collection elements should be defined in
  // Note: additionalItems is **not** analogous to additionalProperties
  static const auto arrayElementSchema = "items";
  static const auto dictionaryElementSchema = "additionalProperties";
  detail::augmentCollectionPaths(filteredPathName, m_ArrayPaths, arrayElementSchema);
  detail::augmentCollectionPaths(filteredPathName, m_DictionaryPaths, dictionaryElementSchema);
  std::vector<std::string> tokens = utilities::string::split(filteredPathName, '/');
  auto iter = std::find(tokens.begin(), tokens.end(), detail::COLLECTION_GROUP_NAME);
  // Replace collection group annotations with a token corresponding to the correct path
  while(iter != tokens.end())
  {
    auto afterRemoved = tokens.erase(iter);
    // If the collection tag wasn't the last element, collapse the thing following it
    // with the thing preceding it
    if(afterRemoved != tokens.end())
    {
      // e.g. for {"foo", COLLECTION_GROUP_NAME, "baz"}, we want just {"foo/baz"}
      // so after removing the GROUP_NAME we prepend to the "baz" element and remove
      // the "foo" element
      // in practice "baz" will always be "items" or "additionalProperties"
      *afterRemoved = utilities::string::appendPrefix(*(afterRemoved - 1), *afterRemoved);
      tokens.erase(afterRemoved - 1);
    }
    iter = std::find(tokens.begin(), tokens.end(), detail::COLLECTION_GROUP_NAME);
  }
  std::string containerPath = fmt::format("properties/{}", fmt::join(tokens, "/properties/"));

  conduit::Node& containerNode =
    container.name().empty() ? m_schemaRoot : m_schemaRoot[containerPath];

  // Annotate with the "object" as the default since almost everything nested
  // (structs and dictionaries of structs) is considered an "object"
  if(!containerNode.has_child("type"))
  {
    containerNode["type"] = "object";
  }

  // We need to record the names of collections so they can be used to properly
  // query the schema when documenting their elements
  if(isCollectionGroup(container.name()))
  {
    auto indices = detail::collectionIndices(container);
    // Check to see if it's an array - have to watch out for mixed index types
    if(std::all_of(indices.begin(), indices.end(), [](const VariantKey& key) {
         return key.type() == InletType::Integer;
       }))
    {
      containerNode["type"] = "array";
      m_ArrayPaths.push_back(filteredPathName);
    }
    else
    {
      m_DictionaryPaths.push_back(filteredPathName);
    }
  }

  if(sidreGroup->getName() != "" && sidreGroup->hasView("description"))
  {
    containerNode["description"] = std::string(sidreGroup->getView("description")->getString());
  }

  // If this is the collection container for a primitive array, we need to treat it specially
  // (by only visiting one element to get the type info), otherwise, visit each sub-field
  if(isCollectionGroup(container.name()) && !sidreGroup->hasView(detail::STRUCT_COLLECTION_FLAG) &&
     !container.getChildFields().empty())
  {
    const auto& location =
      (containerNode["type"].as_string() == "array") ? arrayElementSchema : dictionaryElementSchema;
    detail::recordFieldSchema(*container.getChildFields().begin()->second, containerNode[location]);
  }
  else
  {
    for(const auto& fieldEntry : container.getChildFields())
    {
      const auto name = utilities::string::removeBeforeDelimiter(fieldEntry.first);
      auto& childNode = containerNode["properties"][name];
      detail::recordFieldSchema(*fieldEntry.second, childNode);
      if(fieldEntry.second->isRequired())
      {
        containerNode["required"].append() = name;
      }
    }
  }
}

void JSONSchemaWriter::finalize()
{
  detail::pruneEmptyNodes(m_schemaRoot);
  m_schemaRoot.save(m_fileName, "json");
}

}  // namespace inlet
}  // namespace axom
