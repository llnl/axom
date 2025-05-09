// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file ConduitReader.cpp
 *
 * \brief This file contains the class implementation of the ConduitReader.
 *******************************************************************************
 */

#include "axom/inlet/ConduitReader.hpp"

#include <fstream>
#include <numeric>

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/inlet/inlet_utils.hpp"

#include "axom/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
ConduitReader::ConduitReader(const std::string& protocol) : m_protocol(protocol)
{
  SLIC_ERROR_IF((protocol != "yaml") && (protocol != "json"),
                fmt::format("Inlet: Only 'json' and 'yaml' protocols are "
                            "supported by ConduitReader, got: {0}",
                            protocol));
}

bool ConduitReader::parseFile(const std::string& filePath)
{
  if(!axom::utilities::filesystem::pathExists(filePath))
  {
    SLIC_WARNING(fmt::format("Inlet: Given input file does not exist: {0}", filePath));
    return false;
  }
  bool success = true;
  // Temporarily enable exceptions for Conduit errors so they can be caught and displayed
  sidre::DataStore::setConduitDefaultMessageHandlers();
  try
  {
    m_root.load(filePath, m_protocol);
  }
  catch(const conduit::Error& e)
  {
    SLIC_WARNING(fmt::format("[Inlet]: Failed to parse {0}:\n{1}", m_protocol, e.message()));
    success = false;
  }
  sidre::DataStore::setConduitSLICMessageHandlers();
  return success;
}

bool ConduitReader::parseString(const std::string& stringToRead)
{
  if(stringToRead.empty())
  {
    SLIC_WARNING("Inlet: Given an empty string to parse.");
    return false;
  }
  bool success = true;
  sidre::DataStore::setConduitDefaultMessageHandlers();
  try
  {
    m_root.parse(stringToRead, m_protocol);
  }
  catch(const conduit::Error& e)
  {
    SLIC_WARNING(fmt::format("[Inlet]: Failed to parse {0}:\n{1}", m_protocol, e.message()));
    success = false;
  }
  sidre::DataStore::setConduitSLICMessageHandlers();
  return success;
}

namespace detail
{
// TODO: allow alternate delimiter at sidre level
const static char SCOPE_DELIMITER = '/';

/*!
 *******************************************************************************
 * \brief Traverses a path starting from a Node
 *
 * \param [in] root The node from which to begin traversal
 * \param [in] id   The path to traverse as an Inlet path
 * \note Needed for paths containing integers as a conversion is required
 * \note Inlet paths differ from Conduit paths in that they can contain integers,
 * e.g., "/path/to/7/foo" would correspond to node["path"]["to"][7]["foo"].
 * The traversal prefers the string-valued name but will attempt to convert to
 * integer if no such child exists.
 *******************************************************************************
 */
const conduit::Node* traverseNode(const conduit::Node& root, const std::string& id)
{
  if(root.has_path(id))
  {
    return &root[id];
  }

  const conduit::Node* node = &root;
  std::vector<std::string> tokens = axom::utilities::string::split(id, SCOPE_DELIMITER);
  for(const auto& token : tokens)
  {
    // Prefer the string name, but if it doesn't exist, try converting to int
    if(node->has_child(token))
    {
      node = &((*node)[token]);
    }
    else
    {
      bool is_int = conduit::utils::string_is_integer(token);
      int token_as_int = conduit::utils::string_to_value<int>(token);
      if(is_int && (token_as_int < node->number_of_children()))
      {
        node = &((*node)[token_as_int]);
      }
      else
      {
        // Bail out and return an empty node to indicate failure
        return {};
      }
    }
  }
  return node;
}

/*!
 *******************************************************************************
 * \brief Copies a Conduit array into an unordered map, using zero-based array
 * indices as the keys
 *
 * \param [in] array The array to copy from
 * \param [out] map  The map to copy into
 * \note Implementing to allow for widening/narrowing conversions
 *******************************************************************************
 */
template <typename ConduitType, typename MapValueType>
void arrayToMap(const conduit::DataArray<ConduitType>& array,
                std::unordered_map<int, MapValueType>& map)
{
  map.clear();
  for(conduit::index_t i = 0; i < array.number_of_elements(); i++)
  {
    // No begin/end iterators are provided by DataArray
    map[i] = array[i];
  }
}

/*!
 *******************************************************************************
 * \brief Recursive name retrieval function - adds the names of all descendents
 * of @p node as an Inlet-style path
 * 
 * \param [in] node The Conduit node to "visit"
 * \param [out] names The set of paths to add to
 *******************************************************************************
 */
void nameRetrievalHelper(const conduit::Node& node, std::vector<std::string>& names)
{
  // Conduit paths use [0] for array indices, Inlet does not, so they need
  // to be removed - e.g., foo/[0]/bar vs foo/0/bar
  auto filter_name = [](std::string name) {
    name.erase(std::remove(name.begin(), name.end(), '['), name.end());
    name.erase(std::remove(name.begin(), name.end(), ']'), name.end());
    return name;
  };
  for(const auto& child : node.children())
  {
    names.push_back(filter_name(child.path()));
    nameRetrievalHelper(child, names);
  }
}

}  // namespace detail

ReaderResult ConduitReader::getValue(const conduit::Node* node, int& value)
{
  if(!node)
  {
    return ReaderResult::NotFound;
  }
  // Match LuaReader functionality - narrow from floating-point but exclude bool
  if(node->dtype().is_number() && !node->dtype().is_uint8())
  {
    value = node->to_int();
    return ReaderResult::Success;
  }
  return node->dtype().is_empty() ? ReaderResult::NotFound : ReaderResult::WrongType;
}

ReaderResult ConduitReader::getValue(const conduit::Node* node, std::string& value)
{
  if(!node)
  {
    return ReaderResult::NotFound;
  }
  if(node->dtype().is_string())
  {
    value = node->as_string();
    return ReaderResult::Success;
  }
  return node->dtype().is_empty() ? ReaderResult::NotFound : ReaderResult::WrongType;
}

ReaderResult ConduitReader::getValue(const conduit::Node* node, double& value)
{
  if(!node)
  {
    return ReaderResult::NotFound;
  }
  // Match LuaReader functionality - promote from integer but not bool
  if(node->dtype().is_number() && !node->dtype().is_uint8())
  {
    value = node->to_double();
    return ReaderResult::Success;
  }
  return node->dtype().is_empty() ? ReaderResult::NotFound : ReaderResult::WrongType;
}

ReaderResult ConduitReader::getValue(const conduit::Node* node, bool& value)
{
  if(!node)
  {
    return ReaderResult::NotFound;
  }
  // Boolean literals don't appear to be parsed as such - they are strings
  if((m_protocol == "yaml") && node->dtype().is_string())
  {
    std::string as_str = node->as_string();
    // YAML 1.2 spec, section 10.3.2
    // FIXME: Converting the string to lowercase is not strictly correct, it
    // allows for things like tRue and falsE
    utilities::string::toLower(as_str);
    if(as_str == "true")
    {
      value = true;
      return ReaderResult::Success;
    }
    else if(as_str == "false")
    {
      value = false;
      return ReaderResult::Success;
    }
  }
  else if((m_protocol == "json") && node->dtype().is_uint8())
  {
    value = node->as_uint8();
    return ReaderResult::Success;
  }
  return node->dtype().is_empty() ? ReaderResult::NotFound : ReaderResult::WrongType;
}

ReaderResult ConduitReader::getBool(const std::string& id, bool& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

ReaderResult ConduitReader::getDouble(const std::string& id, double& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

ReaderResult ConduitReader::getInt(const std::string& id, int& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

ReaderResult ConduitReader::getString(const std::string& id, std::string& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

ReaderResult ConduitReader::getIntMap(const std::string& id, std::unordered_map<int, int>& values)
{
  return getArray(id, values);
}

ReaderResult ConduitReader::getDoubleMap(const std::string& id,
                                         std::unordered_map<int, double>& values)
{
  return getArray(id, values);
}

ReaderResult ConduitReader::getBoolMap(const std::string& id, std::unordered_map<int, bool>& values)
{
  return getArray(id, values);
}

ReaderResult ConduitReader::getStringMap(const std::string& id,
                                         std::unordered_map<int, std::string>& values)
{
  return getArray(id, values);
}

ReaderResult ConduitReader::getIntMap(const std::string& id,
                                      std::unordered_map<VariantKey, int>& values)
{
  return getDictionary(id, values);
}

ReaderResult ConduitReader::getDoubleMap(const std::string& id,
                                         std::unordered_map<VariantKey, double>& values)
{
  return getDictionary(id, values);
}

ReaderResult ConduitReader::getBoolMap(const std::string& id,
                                       std::unordered_map<VariantKey, bool>& values)
{
  return getDictionary(id, values);
}

ReaderResult ConduitReader::getStringMap(const std::string& id,
                                         std::unordered_map<VariantKey, std::string>& values)
{
  return getDictionary(id, values);
}

ReaderResult ConduitReader::getIndices(const std::string& id, std::vector<int>& indices)
{
  indices.clear();
  const auto node_ptr = detail::traverseNode(m_root, id);
  if(!node_ptr)
  {
    return ReaderResult::NotFound;
  }
  const auto& node = *node_ptr;
  int num_elements = node.number_of_children();
  // Primitive arrays do not count as lists
  if(!node.dtype().is_list())
  {
    num_elements = node.dtype().number_of_elements();
  }

  indices.resize(num_elements);
  // Arrays in YAML/JSON are contiguous so we don't need to query the input file
  std::iota(indices.begin(), indices.end(), 0);
  return ReaderResult::Success;
}

ReaderResult ConduitReader::getIndices(const std::string& id, std::vector<VariantKey>& indices)
{
  indices.clear();
  const auto node_ptr = detail::traverseNode(m_root, id);
  if(!node_ptr)
  {
    return ReaderResult::NotFound;
  }
  const auto& node = *node_ptr;
  if(!node.dtype().is_object())
  {
    // If it's not an object, try integer indexing
    std::vector<int> int_indices;
    const auto result = getIndices(id, int_indices);
    if(result == ReaderResult::Success)
    {
      for(const int idx : int_indices)
      {
        indices.emplace_back(idx);
      }
    }
    return result;
  }
  for(const auto& child : node.children())
  {
    indices.push_back(child.name());
  }
  return ReaderResult::Success;
}

FunctionVariant ConduitReader::getFunction(const std::string&,
                                           const FunctionTag,
                                           const std::vector<FunctionTag>&)
{
  SLIC_ERROR("[Inlet] Conduit YAML/JSON does not support functions");
  return {};
}

std::vector<std::string> ConduitReader::getAllNames()
{
  std::vector<std::string> result;
  detail::nameRetrievalHelper(m_root, result);
  return result;
}

template <typename T>
ReaderResult ConduitReader::getDictionary(const std::string& id,
                                          std::unordered_map<VariantKey, T>& values)
{
  values.clear();
  const auto node_ptr = detail::traverseNode(m_root, id);
  if(!node_ptr)
  {
    return ReaderResult::NotFound;
  }
  const auto& node = *node_ptr;
  // If it's empty, then the dictionary must have been empty, which counts as successful
  if(node.dtype().is_empty())
  {
    return ReaderResult::Success;
  }
  if(!node.dtype().is_object())
  {
    return ReaderResult::WrongType;
  }

  bool contains_other_type = false;
  for(const auto& child : node.children())
  {
    const auto name = child.name();

    T value;
    // Inlet allows for heterogenous collections, but a failure here must be reported
    const auto result = getValue(&child, value);
    if(result == ReaderResult::Success)
    {
      values[name] = value;
    }
    else
    {
      contains_other_type = true;
    }
  }
  return collectionRetrievalResult(contains_other_type, !values.empty());
}

template <typename T>
ReaderResult ConduitReader::getArray(const std::string& id, std::unordered_map<int, T>& values)
{
  values.clear();
  const auto node_ptr = detail::traverseNode(m_root, id);
  if(!node_ptr)
  {
    return ReaderResult::NotFound;
  }
  const auto& node = *node_ptr;
  // If it's empty, then the array must have been empty, which counts as successful
  if(node.dtype().is_empty())
  {
    return ReaderResult::Success;
  }
  // Truly primitive (i.e., not string) types are contiguous so we grab the array pointer
  else if(node.dtype().number_of_elements() > 1)
  {
    // The template parameter is not enough to know the type of the conduit array
    // as widening/narrowing conversions are supported
    if(node.dtype().is_floating_point())
    {
      detail::arrayToMap(node.as_double_array(), values);
    }
    else if(node.dtype().is_int32())
    {
      detail::arrayToMap(node.as_int32_array(), values);
    }
    else if(node.dtype().is_int64())
    {
      detail::arrayToMap(node.as_int64_array(), values);
    }
    else
    {
      return ReaderResult::WrongType;
    }
  }
  else if(!node.dtype().is_list() && !node.dtype().is_object())
  {
    // Single-element arrays will be just the element itself
    // If it's a single element, we know the index is zero
    T value;
    const auto result = getValue(&node, value);
    if(result == ReaderResult::Success)
    {
      values[0] = value;
    }
    else
    {
      return result;
    }
  }
  // String arrays are not supported natively so the node is directly iterated over
  else
  {
    conduit::index_t index = 0;
    bool contains_other_type = false;
    for(const auto& child : node.children())
    {
      T value;
      // Inlet allows for heterogenous collections, but a failure here must be reported
      const auto result = getValue(&child, value);
      if(result == ReaderResult::Success)
      {
        values[index] = value;
      }
      else
      {
        contains_other_type = true;
      }
      index++;
    }
    return collectionRetrievalResult(contains_other_type, !values.empty());
  }
  return ReaderResult::Success;
}

}  // end namespace inlet
}  // end namespace axom
