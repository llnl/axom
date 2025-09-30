// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Container.hpp"

#include "axom/slic.hpp"
#include "axom/core/Path.hpp"
#include "axom/inlet/inlet_utils.hpp"
#include "axom/inlet/Proxy.hpp"

namespace axom
{
namespace inlet
{
Container::Container(const std::string& name,
                     const std::string& description,
                     Reader& reader,
                     axom::sidre::Group* sidreRootGroup,
                     std::vector<std::string>& unexpectedNames,
                     bool docEnabled,
                     bool reconstruct)
  : m_name(name)
  , m_reader(reader)
  , m_sidreRootGroup(sidreRootGroup)
  , m_unexpectedNames(unexpectedNames)
  , m_docEnabled(docEnabled)
{
  SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr, "Inlet's Sidre Datastore class not set");

  if(m_name == "")
  {
    m_sidreGroup = m_sidreRootGroup;
  }
  else
  {
    if(!m_sidreRootGroup->hasGroup(name))
    {
      m_sidreGroup = m_sidreRootGroup->createGroup(name);
      m_sidreGroup->createViewString("InletType", "Container");
    }
    else
    {
      m_sidreGroup = m_sidreRootGroup->getGroup(name);
    }
  }

  if(reconstruct)
  {
    for(auto& group : m_sidreGroup->groups())
    {
      if(group.isUsingMap() && group.hasView("InletType"))
      {
        const std::string inletType = group.getView("InletType")->getString();
        const std::string childName = utilities::string::appendPrefix(m_name, group.getName());

        if(inletType == "Container")
        {
          m_containerChildren.emplace(childName,
                                      std::make_unique<Container>(childName,
                                                                  "",
                                                                  m_reader,
                                                                  m_sidreRootGroup,
                                                                  m_unexpectedNames,
                                                                  m_docEnabled,
                                                                  true));
        }
        else if(inletType == "Field")
        {
          // FIXME: We probably need to write the type to the datastore
          m_fieldChildren.emplace(childName,
                                  std::make_unique<Field>(&group,
                                                          m_sidreRootGroup,
                                                          sidre::DataTypeId::NO_TYPE_ID,
                                                          m_docEnabled));
        }
      }
    }
  }

  if(!description.empty())
  {
    if(m_sidreGroup->hasView("description"))
    {
      //TODO: warn user?
      m_sidreGroup->destroyViewAndData("description");
    }
    m_sidreGroup->createViewString("description", description);
  }
}

template <typename Func>
void Container::forEachCollectionElement(Func&& func) const
{
  for(const auto& index : detail::collectionIndices(*this))
  {
    func(getContainer(detail::indexToString(index)));
  }
}

template <typename OutputIt, typename Func>
bool Container::transformFromNestedElements(OutputIt output, const std::string& name, Func&& func) const
{
  for(Container& container : m_nested_aggregates)
  {
    *output++ = func(container, {});
  }

  if(isStructCollection())
  {
    for(const auto& indexPath : detail::collectionIndicesWithPaths(*this, name))
    {
      *output++ = func(getContainer(indexPath.first), indexPath.second);
    }
  }
  return isStructCollection() || !m_nested_aggregates.empty();
}

Container& Container::addContainer(const std::string& name, const std::string& description)
{
  auto currContainer = this;
  Path path(name);

  // Create intermediate Containers if they don't already exist
  for(const auto& pathPart : path.parts())
  {
    // Need to watch out for empty paths here
    auto currContainerName = Path::join({currContainer->m_name, pathPart});
    // Leave the description empty for intermediate containers
    const std::string currDescr = (pathPart == name) ? description : "";
    if(!currContainer->hasChild<Container>(pathPart))
    {
      // Will the copy always be elided here with a move ctor
      // or do we need std::piecewise_construct/std::forward_as_tuple?
      const auto& emplaceResult =
        currContainer->m_containerChildren.emplace(currContainerName,
                                                   std::make_unique<Container>(currContainerName,
                                                                               currDescr,
                                                                               m_reader,
                                                                               m_sidreRootGroup,
                                                                               m_unexpectedNames,
                                                                               m_docEnabled));
      // emplace_result is a pair whose first element is an iterator to the inserted element
      currContainer = emplaceResult.first->second.get();
    }
    else
    {
      currContainer = currContainer->m_containerChildren[currContainerName].get();
    }
  }

  return *currContainer;
}

Container& Container::addStruct(const std::string& name, const std::string& description)
{
  auto& base_container = addContainer(name, description);
  transformFromNestedElements(
    std::back_inserter(base_container.m_nested_aggregates),
    name,
    [&name, &description](Container& subcontainer, const std::string&) -> Container& {
      return subcontainer.addStruct(name, description);
    });
  return base_container;
}

Verifiable<Container>& Container::addBoolArray(const std::string& name, const std::string& description)
{
  return addPrimitiveArray<bool>(name, description);
}

Verifiable<Container>& Container::addIntArray(const std::string& name, const std::string& description)
{
  return addPrimitiveArray<int>(name, description);
}

Verifiable<Container>& Container::addDoubleArray(const std::string& name,
                                                 const std::string& description)
{
  return addPrimitiveArray<double>(name, description);
}

Verifiable<Container>& Container::addStringArray(const std::string& name,
                                                 const std::string& description)
{
  return addPrimitiveArray<std::string>(name, description);
}

template <typename Key>
Container& Container::addStructCollection(const std::string& name, const std::string& description)
{
  auto& container =
    addContainer(utilities::string::appendPrefix(name, detail::COLLECTION_GROUP_NAME), description);

  transformFromNestedElements(
    std::back_inserter(container.m_nested_aggregates),
    name,
    [&name, &description](Container& subcontainer, const std::string&) -> Container& {
      return subcontainer.addStructCollection<Key>(name, description);
    });

  if(isStructCollection())
  {
    markAsStructCollection(*container.m_sidreGroup);
  }
  else
  {
    std::vector<Key> indices;
    std::string fullName = utilities::string::appendPrefix(m_name, name);
    fullName = utilities::string::removeAllInstances(fullName, detail::COLLECTION_GROUP_NAME + "/");
    detail::updateUnexpectedNames(fullName, m_unexpectedNames);
    const auto result = m_reader.getIndices(fullName, indices);
    if(result == ReaderResult::Success)
    {
      container.addIndicesGroup(indices, description, true);
    }
    markRetrievalStatus(*container.m_sidreGroup, result);
    markAsStructCollection(*container.m_sidreGroup);
  }
  return container;
}

Container& Container::addStructArray(const std::string& name, const std::string& description)
{
  return addStructCollection<int>(name, description);
}

Verifiable<Container>& Container::addBoolDictionary(const std::string& name,
                                                    const std::string& description)
{
  return addPrimitiveArray<bool>(name, description, true);
}

Verifiable<Container>& Container::addIntDictionary(const std::string& name,
                                                   const std::string& description)
{
  return addPrimitiveArray<int>(name, description, true);
}

Verifiable<Container>& Container::addDoubleDictionary(const std::string& name,
                                                      const std::string& description)
{
  return addPrimitiveArray<double>(name, description, true);
}

Verifiable<Container>& Container::addStringDictionary(const std::string& name,
                                                      const std::string& description)
{
  return addPrimitiveArray<std::string>(name, description, true);
}

Container& Container::addStructDictionary(const std::string& name, const std::string& description)
{
  return addStructCollection<VariantKey>(name, description);
}

axom::sidre::Group* Container::createSidreGroup(const std::string& name,
                                                const std::string& description)
{
  SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr, "[Inlet] Sidre Datastore Group not set");

  if(m_sidreRootGroup->hasGroup(name))
  {
    SLIC_WARNING("[Inlet] Cannot add value that already exists: " + name);
    setWarningFlag(m_sidreRootGroup);
    return nullptr;
  }

  axom::sidre::Group* sidreGroup = m_sidreRootGroup->createGroup(name);
  sidreGroup->createViewString("InletType", "Field");
  SLIC_ASSERT_MSG(sidreGroup != nullptr, "[Inlet] Sidre failed to create group");
  if(description != "")
  {
    sidreGroup->createViewString("description", description);
  }

  return sidreGroup;
}

Field& Container::addField(axom::sidre::Group* sidreGroup,
                           axom::sidre::DataTypeId type,
                           const std::string& fullName,
                           const std::string& name)
{
  const size_t found = name.find_last_of("/");
  auto currContainer = this;
  if(found != std::string::npos)
  {
    // This will add any intermediate Containers (if not present) before adding the field
    currContainer = &addContainer(name.substr(0, found));
  }
  const auto& emplace_result = currContainer->m_fieldChildren.emplace(
    fullName,
    std::make_unique<Field>(sidreGroup, m_sidreRootGroup, type, m_docEnabled));
  // emplace_result is a pair whose first element is an iterator to the inserted element
  return *(emplace_result.first->second);
}

Function& Container::addFunctionInternal(axom::sidre::Group* sidreGroup,
                                         FunctionVariant&& func,
                                         const std::string& fullName,
                                         const std::string& name)
{
  const size_t found = name.find_last_of("/");
  auto currContainer = this;
  if(found != std::string::npos)
  {
    // This will add any intermediate Containers (if not present) before adding the field
    currContainer = &addContainer(name.substr(0, found));
  }
  const auto& emplace_result = currContainer->m_functionChildren.emplace(
    fullName,
    std::make_unique<Function>(sidreGroup, m_sidreRootGroup, std::move(func)));
  // emplace_result is a pair whose first element is an iterator to the inserted element
  return *(emplace_result.first->second);
}

template <typename T, typename SFINAE>
VerifiableScalar& Container::addPrimitive(const std::string& name,
                                          const std::string& description,
                                          bool forArray,
                                          T val,
                                          const std::string& pathOverride)
{
  // If it has indices, we're adding a primitive field to an array
  // of structs, so we need to iterate over the subcontainers
  // corresponding to elements of the array
  std::vector<std::reference_wrapper<VerifiableScalar>> fields;
  const bool is_nested = transformFromNestedElements(
    std::back_inserter(fields),
    name,
    [&name, &description, forArray, &val](Container& subcontainer,
                                          const std::string& path) -> VerifiableScalar& {
      return subcontainer.addPrimitive<T>(name, description, forArray, val, path);
    });

  if(is_nested)
  {
    // Create an aggregate field so requirements can be collectively imposed
    // on all elements of the array
    m_aggregate_fields.emplace_back(std::move(fields));

    // Remove when C++17 is available
    return m_aggregate_fields.back();
  }
  else
  {
    // Otherwise actually add a Field
    std::string fullName = utilities::string::appendPrefix(m_name, name);
    // First check if the field already exists
    auto iter = m_fieldChildren.find(fullName);
    if(iter != m_fieldChildren.end())
    {
      return *iter->second;
    }
    axom::sidre::Group* sidreGroup = createSidreGroup(fullName, description);
    SLIC_ERROR_IF(sidreGroup == nullptr,
                  fmt::format("Failed to create Sidre group with name '{0}'", fullName));
    // If a pathOverride is specified, needed when Inlet-internal groups
    // are part of fullName
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    lookupPath =
      utilities::string::removeAllInstances(lookupPath, detail::COLLECTION_GROUP_NAME + "/");
    detail::updateUnexpectedNames(lookupPath, m_unexpectedNames);
    auto typeId = addPrimitiveHelper(sidreGroup, lookupPath, forArray, val);
    return addField(sidreGroup, typeId, fullName, name);
  }
}

template <>
axom::sidre::DataTypeId Container::addPrimitiveHelper<bool>(axom::sidre::Group* sidreGroup,
                                                            const std::string& lookupPath,
                                                            bool forArray,
                                                            bool val)
{
  const auto result = m_reader.getBool(lookupPath, val);
  if(forArray || result == ReaderResult::Success)
  {
    sidreGroup->createViewScalar("value", val ? std::int8_t(1) : std::int8_t(0));
  }
  if(!forArray)
  {
    markRetrievalStatus(*sidreGroup, result);
  }
  return axom::sidre::DataTypeId::INT8_ID;
}

template <>
axom::sidre::DataTypeId Container::addPrimitiveHelper<int>(axom::sidre::Group* sidreGroup,
                                                           const std::string& lookupPath,
                                                           bool forArray,
                                                           int val)
{
  const auto result = m_reader.getInt(lookupPath, val);
  if(forArray || result == ReaderResult::Success)
  {
    sidreGroup->createViewScalar("value", val);
  }
  if(!forArray)
  {
    markRetrievalStatus(*sidreGroup, result);
  }
  return axom::sidre::DataTypeId::INT_ID;
}

template <>
axom::sidre::DataTypeId Container::addPrimitiveHelper<double>(axom::sidre::Group* sidreGroup,
                                                              const std::string& lookupPath,
                                                              bool forArray,
                                                              double val)
{
  const auto result = m_reader.getDouble(lookupPath, val);
  if(forArray || result == ReaderResult::Success)
  {
    sidreGroup->createViewScalar("value", val);
  }
  if(!forArray)
  {
    markRetrievalStatus(*sidreGroup, result);
  }
  return axom::sidre::DataTypeId::DOUBLE_ID;
}

template <>
axom::sidre::DataTypeId Container::addPrimitiveHelper<std::string>(axom::sidre::Group* sidreGroup,
                                                                   const std::string& lookupPath,
                                                                   bool forArray,
                                                                   std::string val)
{
  const auto result = m_reader.getString(lookupPath, val);
  if(forArray || result == ReaderResult::Success)
  {
    sidreGroup->createViewString("value", val);
  }
  if(!forArray)
  {
    markRetrievalStatus(*sidreGroup, result);
  }
  return axom::sidre::DataTypeId::CHAR8_STR_ID;
}

namespace detail
{
/*!
  *****************************************************************************
  * \brief Adds the contents of an array to the container
  * 
  * \return The keys that were added
  *****************************************************************************
  */
template <typename T>
std::vector<VariantKey> registerCollection(Container& container,
                                           const std::unordered_map<int, T>& collection)
{
  std::vector<VariantKey> result;
  for(const auto& entry : collection)
  {
    result.push_back(entry.first);
    container.addPrimitive(std::to_string(entry.first), "", true, entry.second);
  }
  return result;
}

/*!
  *****************************************************************************
  * \brief Adds the contents of a dict to the container
  * 
  * \return The keys that were added
  *****************************************************************************
  */
template <typename T>
std::vector<VariantKey> registerCollection(Container& container,
                                           const std::unordered_map<VariantKey, T>& collection)
{
  std::vector<VariantKey> result;
  for(const auto& entry : collection)
  {
    result.push_back(entry.first);
    auto string_key = indexToString(entry.first);
    const auto illegal_char_loc = string_key.find_first_of("/[]");
    SLIC_ERROR_IF(illegal_char_loc != std::string::npos,
                  fmt::format("[Inlet] Dictionary key '{0}' contains illegal character '{1}'",
                              string_key,
                              string_key[illegal_char_loc]));
    SLIC_ERROR_IF(string_key.empty(), "[Inlet] Dictionary key cannot be the empty string");
    container.addPrimitive(string_key, "", true, entry.second);
  }
  return result;
}

/*!
 *****************************************************************************
 * \brief Implementation helper for adding primitive arrays
 * 
 * \note Structs are used for partial template specializations
 *****************************************************************************
 */
template <typename Key, typename Primitive>
struct PrimitiveArrayHelper
{ };

template <typename Key>
struct PrimitiveArrayHelper<Key, bool>
{
  /*!
   *****************************************************************************
   * \brief Finalizes the creation of a collection
   * \param [inout] container The container to add the collection to
   * \param [in] reader The Reader object to read the collection from
   * \param [in] lookupPath The path within the input file to the collection
   * 
   * \return The keys from the collection that was added
   *****************************************************************************
   */
  static std::vector<VariantKey> add(Container& container, Reader& reader, const std::string& lookupPath)
  {
    std::unordered_map<Key, bool> map;
    // Failure to retrieve a map is not necessarily an error
    const auto result = reader.getBoolMap(lookupPath, map);
    markRetrievalStatus(*container.sidreGroup(), result);
    return registerCollection(container, map);
  }
};

template <typename Key>
struct PrimitiveArrayHelper<Key, int>
{
  static std::vector<VariantKey> add(Container& container, Reader& reader, const std::string& lookupPath)
  {
    std::unordered_map<Key, int> map;
    const auto result = reader.getIntMap(lookupPath, map);
    markRetrievalStatus(*container.sidreGroup(), result);
    return registerCollection(container, map);
  }
};

template <typename Key>
struct PrimitiveArrayHelper<Key, double>
{
  static std::vector<VariantKey> add(Container& container, Reader& reader, const std::string& lookupPath)
  {
    std::unordered_map<Key, double> map;
    const auto result = reader.getDoubleMap(lookupPath, map);
    markRetrievalStatus(*container.sidreGroup(), result);
    return registerCollection(container, map);
  }
};

template <typename Key>
struct PrimitiveArrayHelper<Key, std::string>
{
  static std::vector<VariantKey> add(Container& container, Reader& reader, const std::string& lookupPath)
  {
    std::unordered_map<Key, std::string> map;
    const auto result = reader.getStringMap(lookupPath, map);
    markRetrievalStatus(*container.sidreGroup(), result);
    return registerCollection(container, map);
  }
};

void addIndexViewToGroup(sidre::Group& group, const int& index)
{
  group.createViewScalar("", index);
}

void addIndexViewToGroup(sidre::Group& group, const std::string& index)
{
  group.createViewString("", index);
}

void addIndexViewToGroup(sidre::Group& group, const VariantKey& index)
{
  if(index.type() == InletType::String)
  {
    addIndexViewToGroup(group, static_cast<std::string>(index));
  }
  else
  {
    addIndexViewToGroup(group, static_cast<int>(index));
  }
}

/*!
 *****************************************************************************
 * \brief Writes function signature information to the sidre Group associated
 * with an inlet::Function
 * 
 * \param [in] ret_type The function's return type
 * \param [in] arg_types The function's argument types
 * \param [inout] group The group to add to
 * 
 * The structure of the function signature information is as follows:
 * <group>
 *  ├─• function_arguments <integer array>
 *  └─• return_type <integer>
 *****************************************************************************
 */
void addSignatureToGroup(const FunctionTag ret_type,
                         const std::vector<FunctionTag>& arg_types,
                         sidre::Group* group)
{
  group->createViewScalar("return_type", static_cast<int>(ret_type));
  auto args_view =
    group->createViewAndAllocate("function_arguments", sidre::INT_ID, arg_types.size());
  int* args_array = args_view->getArray();
  for(const auto arg_type : arg_types)
  {
    *args_array++ = static_cast<int>(arg_type);
  }
}

std::vector<VariantKey> collectionIndices(const Container& container, bool trimAbsolute)
{
  std::vector<VariantKey> indices;
  const auto sidreGroup = container.sidreGroup();
  // Not having indices is not necessarily an error, as the collection
  // could exist but just be empty
  if(sidreGroup->hasGroup(detail::COLLECTION_INDICES_NAME))
  {
    const auto group = sidreGroup->getGroup(detail::COLLECTION_INDICES_NAME);
    indices.reserve(group->getNumViews());
    for(auto& view : group->views())
    {
      if(view.getTypeID() == axom::sidre::CHAR8_STR_ID)
      {
        std::string string_idx = view.getString();
        VariantKey key = string_idx;
        if(trimAbsolute)
        {
          // If the index is full/absolute, we only care about the last segment of it
          string_idx = utilities::string::removeBeforeDelimiter(string_idx);
          // The basename might be an integer, so check and convert accordingly
          if(conduit::utils::string_is_integer(string_idx))
          {
            key = conduit::utils::string_to_value<int>(string_idx);
          }
          else
          {
            key = string_idx;
          }
        }
        indices.push_back(key);
      }
      else
      {
        indices.push_back(static_cast<int>(view.getData()));
      }
    }
  }
  return indices;
}

std::vector<std::pair<std::string, std::string>> collectionIndicesWithPaths(const Container& container,
                                                                            const std::string& name)
{
  std::vector<std::pair<std::string, std::string>> result;
  for(const auto& index : collectionIndices(container, false))
  {
    Path indexPath = detail::indexToString(index);
    // Since the index is absolute, we only care about the last segment of it
    // But since it's an absolute path then it gets used as the fullPath
    // which is used by the Reader to search in the input file
    const auto fullPath = Path::join({indexPath, name});
    // e.g. fullPath could be foo/1/bar for field "bar" at index 1 of array "foo"
    result.push_back({indexPath.baseName(), fullPath});
  }
  return result;
}

void updateUnexpectedNames(const std::string& accessedName, std::vector<std::string>& unexpectedNames)
{
  std::vector<std::string> accessed_tokens = axom::utilities::string::split(accessedName, '/');
  for(auto iter = unexpectedNames.begin(); iter != unexpectedNames.end();)
  {
    // Check if the possibly unexpected name is an "ancestor" of the accessed name,
    // if it is, then it gets marked as expected via removal
    std::vector<std::string> unexpected_tokens = axom::utilities::string::split(*iter, '/');
    // If it's bigger, it can't be an ancestor so we can bail out early
    if(unexpected_tokens.size() > accessed_tokens.size())
    {
      ++iter;
    }
    // Check if the beginning of the accessed tokens sequence is equal to the tokens in the possibly
    // unexpected name
    else if(std::equal(unexpected_tokens.begin(), unexpected_tokens.end(), accessed_tokens.begin()))
    {
      iter = unexpectedNames.erase(iter);
    }
    else
    {
      ++iter;
    }
  }
}

/*!
 *****************************************************************************
 * \brief Filters through the global list of unexpected names to retrieve the
 * ones that are within the Container that corresponds to \a group
 * 
 * \param [in] group The Sidre group for a Container object
 * \param [in] unexpectedNames The global (within the Inlet hierarchy) list
 * of unexpected names
 *****************************************************************************
 */
std::vector<std::string> filterUnexpectedNames(const sidre::Group* group,
                                               const std::vector<std::string>& unexpectedNames)
{
  const std::string callerName = group->getPathName();
  std::vector<std::string> callerTokens = axom::utilities::string::split(callerName, '/');
  callerTokens.erase(
    std::remove(callerTokens.begin(), callerTokens.end(), detail::COLLECTION_GROUP_NAME),
    callerTokens.end());
  // The items of unexpected items below/within the Container corresponding
  // to the "group" argument
  std::vector<std::string> unexpectedNamesWithinGroup;
  for(const auto& name : unexpectedNames)
  {
    std::vector<std::string> unexpectedTokens = axom::utilities::string::split(name, '/');
    // If it's smaller it can't be a descendant
    if(unexpectedTokens.size() >= callerTokens.size())
    {
      // If the beginning of the unexpected tokens sequence is equal to the tokens
      // of the calling Container's name
      if(std::equal(callerTokens.begin(), callerTokens.end(), unexpectedTokens.begin()))
      {
        unexpectedNamesWithinGroup.push_back(name);
      }
    }
  }
  return unexpectedNamesWithinGroup;
}

}  // end namespace detail

template <typename Key>
void Container::addIndicesGroup(const std::vector<Key>& indices,
                                const std::string& description,
                                const bool add_containers)
{
  sidre::Group* indices_group = m_sidreGroup->createGroup(detail::COLLECTION_INDICES_NAME,
                                                          /* list_format = */ true);
  // For each index, add a container whose name is its index
  // Schema for struct is defined using the returned container
  for(const auto& idx : indices)
  {
    const std::string string_idx =
      utilities::string::removeBeforeDelimiter(detail::indexToString(idx));
    if(add_containers)
    {
      addContainer(string_idx, description);
    }
    std::string absolute = utilities::string::appendPrefix(m_name, detail::indexToString(idx));
    absolute = utilities::string::removeAllInstances(absolute, detail::COLLECTION_GROUP_NAME + "/");
    detail::addIndexViewToGroup(*indices_group, absolute);
  }
}

template <typename T, typename SFINAE>
Verifiable<Container>& Container::addPrimitiveArray(const std::string& name,
                                                    const std::string& description,
                                                    const bool isDict,
                                                    const std::string& pathOverride)
{
  // Adding an array of primitive field to an array of structs
  std::vector<std::reference_wrapper<Verifiable>> containers;
  const bool is_nested = transformFromNestedElements(
    std::back_inserter(containers),
    name,
    [&name, &description, isDict](Container& subcontainer,
                                  const std::string& path) -> Verifiable<Container>& {
      return subcontainer.addPrimitiveArray<T>(name, description, isDict, path);
    });

  if(is_nested)
  {
    m_aggregate_containers.emplace_back(std::move(containers));

    // Remove when C++17 is available
    return m_aggregate_containers.back();
  }
  else
  {
    // "base case", create a container for the field and fill it in with the helper
    auto& container =
      addContainer(utilities::string::appendPrefix(name, detail::COLLECTION_GROUP_NAME), description);
    const std::string& fullName = utilities::string::appendPrefix(m_name, name);
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    lookupPath =
      utilities::string::removeAllInstances(lookupPath, detail::COLLECTION_GROUP_NAME + "/");
    detail::updateUnexpectedNames(lookupPath, m_unexpectedNames);
    std::vector<VariantKey> indices;
    if(isDict)
    {
      indices = detail::PrimitiveArrayHelper<VariantKey, T>::add(container, m_reader, lookupPath);
    }
    else
    {
      indices = detail::PrimitiveArrayHelper<int, T>::add(container, m_reader, lookupPath);
    }
    // Copy the indices to the datastore to keep track of integer vs. string indices
    if(!indices.empty())
    {
      container.addIndicesGroup(indices, description, false);
    }
    return container;
  }
}

VerifiableScalar& Container::addBool(const std::string& name, const std::string& description)
{
  return addPrimitive<bool>(name, description);
}

VerifiableScalar& Container::addDouble(const std::string& name, const std::string& description)
{
  return addPrimitive<double>(name, description);
}

VerifiableScalar& Container::addInt(const std::string& name, const std::string& description)
{
  return addPrimitive<int>(name, description);
}

VerifiableScalar& Container::addString(const std::string& name, const std::string& description)
{
  return addPrimitive<std::string>(name, description);
}

Verifiable<Function>& Container::addFunction(const std::string& name,
                                             const FunctionTag ret_type,
                                             const std::vector<FunctionTag>& arg_types,
                                             const std::string& description,
                                             const std::string& pathOverride)
{
  // If it has indices, we're adding a function to an array
  // of structs, so we need to iterate over the subcontainers
  // corresponding to elements of the array
  std::vector<std::reference_wrapper<Verifiable<Function>>> funcs;

  const bool is_nested = transformFromNestedElements(
    std::back_inserter(funcs),
    name,
    [&name, &ret_type, &arg_types, &description](Container& subcontainer,
                                                 const std::string& path) -> Verifiable<Function>& {
      return subcontainer.addFunction(name, ret_type, arg_types, description, path);
    });
  if(is_nested)
  {
    // Create an aggregate function so requirements can be collectively imposed
    // on all elements of the array
    m_aggregate_funcs.emplace_back(std::move(funcs));

    // Remove when C++17 is available
    return m_aggregate_funcs.back();
  }
  else
  {
    // Otherwise actually add a Function
    std::string fullName = utilities::string::appendPrefix(m_name, name);
    // First check if the function already exists
    auto iter = m_functionChildren.find(fullName);
    if(iter != m_functionChildren.end())
    {
      return *iter->second;
    }
    axom::sidre::Group* sidreGroup = createSidreGroup(fullName, description);
    SLIC_ERROR_IF(sidreGroup == nullptr,
                  fmt::format("Failed to create Sidre group with name '{0}'", fullName));
    detail::addSignatureToGroup(ret_type, arg_types, sidreGroup);
    // If a pathOverride is specified, needed when Inlet-internal groups
    // are part of fullName
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    lookupPath =
      utilities::string::removeAllInstances(lookupPath, detail::COLLECTION_GROUP_NAME + "/");
    detail::updateUnexpectedNames(lookupPath, m_unexpectedNames);
    auto func = m_reader.getFunction(lookupPath, ret_type, arg_types);
    return addFunctionInternal(sidreGroup, std::move(func), fullName, name);
  }
}

Proxy Container::operator[](const std::string& name) const
{
  const bool has_container = hasContainer(name);
  const bool has_field = hasField(name);
  const bool has_func = hasFunction(name);

  // Ambiguous case - both a container and field exist with the same name
  if((has_container && has_field) || (has_field && has_func) || (has_container && has_func))
  {
    const std::string msg = fmt::format(
      "[Inlet] Ambiguous lookup - more than one of a container/field/function "
      "with name '{0}' exist",
      name);
    SLIC_ERROR(msg);
    return Proxy();
  }

  else if(has_container)
  {
    return Proxy(getContainer(name));
  }

  else if(has_field)
  {
    return Proxy(getField(name));
  }

  else if(has_func)
  {
    return Proxy(getFunction(name));
  }

  // Neither exists
  else
  {
    std::string msg =
      fmt::format("[Inlet] No container, field, or function with name '{0}' exists", name);
    SLIC_ERROR(msg);
    return Proxy();
  }
}

Container& Container::required(bool isRequired)
{
  // If it's a struct collection we set the individual fields as required,
  // and also the collection container itself, as the user would expect that marking
  // a struct collection as required means that it is non-empty
  if(isStructCollection())
  {
    forEachCollectionElement([isRequired](Container& container) { container.required(isRequired); });
  }

  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Container specific Sidre Datastore Group not set");
  setFlag(*m_sidreGroup, *m_sidreRootGroup, detail::REQUIRED_FLAG, isRequired);
  return *this;
}

bool Container::isRequired() const
{
  if(isStructCollection())
  {
    bool result = false;
    forEachCollectionElement([&result](Container& container) {
      if(container.isRequired())
      {
        result = true;
      }
    });
    return result;
  }

  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Container specific Sidre Datastore Group not set");
  return checkFlag(*m_sidreGroup, *m_sidreRootGroup, detail::REQUIRED_FLAG);
}

Container& Container::strict(bool isStrict)
{
  if(isStructCollection())
  {
    forEachCollectionElement([isStrict](Container& container) { container.strict(isStrict); });
  }
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Container specific Sidre Datastore Group not set");
  setFlag(*m_sidreGroup, *m_sidreRootGroup, detail::STRICT_FLAG, isStrict);
  return *this;
}

Container& Container::registerVerifier(Verifier lambda)
{
  if(isStructCollection())
  {
    forEachCollectionElement([&lambda](Container& container) { container.registerVerifier(lambda); });
  }
  else
  {
    SLIC_WARNING_IF(m_verifier,
                    fmt::format("[Inlet] Verifier for Container "
                                "already set: {0}",
                                m_name));
    m_verifier = lambda;
  }
  return *this;
}

bool Container::verify(std::vector<VerificationError>* errors) const
{
  // Whether the calling container has anything in it
  // If the name is empty then we're the global (root) container, which we always
  // consider to be defined
  const bool this_container_defined = isUserProvided() || m_name.empty();

  // If this container was required, make sure something was defined in it
  bool verified = verifyRequired(*m_sidreGroup, this_container_defined, "Container", errors);

  // Verify this Container if a lambda was configured
  if(this_container_defined && m_verifier && !m_verifier(*this, errors))
  {
    verified = false;
    const std::string msg = fmt::format("[Inlet] Container failed verification: {0}", m_name);
    INLET_VERIFICATION_WARNING(m_name, msg, errors);
  }

  // If the strict flag is set
  if(checkFlag(*m_sidreGroup, *m_sidreRootGroup, detail::STRICT_FLAG))
  {
    // Unexpected names below/within the calling object
    const auto currUnexpectedNames = unexpectedNames();
    verified = verified && currUnexpectedNames.empty();
    for(const auto& name : currUnexpectedNames)
    {
      const std::string msg =
        fmt::format("[Inlet] Container '{0}' contained unexpected child: {1}", m_name, name);
      INLET_VERIFICATION_WARNING(m_name, msg, errors);
    }
  }

  // Checking the child objects is not needed if the container wasn't defined in the input file
  if(this_container_defined)
  {
    // Verify the child Fields of this Container
    for(const auto& field : m_fieldChildren)
    {
      if(!field.second->verify(errors))
      {
        verified = false;
      }
    }
    // Verify the child Containers of this Container
    for(const auto& container : m_containerChildren)
    {
      if(!container.second->verify(errors))
      {
        verified = false;
      }
    }

    // Verify the child Functions of this Container
    for(const auto& function : m_functionChildren)
    {
      if(!function.second->verify(errors))
      {
        verified = false;
      }
    }
  }
  // If this has a collection group, it always needs to be verified, as annotations
  // may have been applied to the collection group and not the calling group
  else if(hasContainer(detail::COLLECTION_GROUP_NAME))
  {
    verified = verified && getContainer(detail::COLLECTION_GROUP_NAME).verify(errors);
  }

  return verified;
}

template <>
std::unordered_map<std::string, std::unique_ptr<Container>> Container::* Container::getChildren()
{
  return &Container::m_containerChildren;
}

template <>
std::unordered_map<std::string, std::unique_ptr<Field>> Container::* Container::getChildren()
{
  return &Container::m_fieldChildren;
}

template <>
std::unordered_map<std::string, std::unique_ptr<Function>> Container::* Container::getChildren()
{
  return &Container::m_functionChildren;
}

template <typename T>
bool Container::hasChild(const std::string& childName) const
{
  const auto& children = this->*getChildren<T>();
  return children.find(utilities::string::appendPrefix(m_name, childName)) != children.end();
}

template <typename T>
T* Container::getChildInternal(const std::string& childName) const
{
  Path path(childName);
  const std::string baseName = path.baseName();
  auto currContainer = this;

  // Traverse the intermediate paths
  const auto parent = path.parent();
  for(const auto& pathPart : parent.parts())
  {
    if(currContainer->hasChild<Container>(pathPart))
    {
      currContainer =
        currContainer->m_containerChildren.at(Path::join({currContainer->m_name, pathPart})).get();
    }
    else
    {
      return nullptr;
    }
  }

  if(currContainer->hasChild<T>(baseName))
  {
    const auto& children = currContainer->*getChildren<T>();
    return children.at(Path::join({currContainer->m_name, baseName})).get();
  }
  return nullptr;
}

bool Container::hasContainer(const std::string& containerName) const
{
  return static_cast<bool>(getChildInternal<Container>(containerName));
}

bool Container::hasField(const std::string& fieldName) const
{
  return static_cast<bool>(getChildInternal<Field>(fieldName));
}

bool Container::hasFunction(const std::string& fieldName) const
{
  return static_cast<bool>(getChildInternal<Function>(fieldName));
}

Container& Container::getContainer(const std::string& containerName) const
{
  auto container = getChildInternal<Container>(containerName);
  if(!container)
  {
    SLIC_ERROR(fmt::format("[Inlet] Container not found: {0}", containerName));
  }
  return *container;
}

Field& Container::getField(const std::string& fieldName) const
{
  auto field = getChildInternal<Field>(fieldName);
  if(!field)
  {
    SLIC_ERROR(fmt::format("[Inlet] Field not found: {0}", fieldName));
  }
  return *field;
}

Function& Container::getFunction(const std::string& funcName) const
{
  auto func = getChildInternal<Function>(funcName);
  if(!func)
  {
    SLIC_ERROR(fmt::format("[Inlet] Function not found: {0}", funcName));
  }
  return *func;
}

std::string Container::name() const { return m_name; }

std::vector<std::string> Container::unexpectedNames() const
{
  return detail::filterUnexpectedNames(m_sidreGroup, m_unexpectedNames);
}

bool Container::contains(const std::string& name) const
{
  if(auto container = getChildInternal<Container>(name))
  {
    // Check if the container itself exists
    return container->exists();
  }
  else if(auto field = getChildInternal<Field>(name))
  {
    // Check if the field itself exists
    return field->exists();
  }
  else if(auto function = getChildInternal<Function>(name))
  {
    // call operator bool on the function itself
    return static_cast<bool>(*function);
  }
  return false;
}

bool Container::exists() const
{
  // Check if any of its child containers exist
  const bool has_containers = std::any_of(
    m_containerChildren.begin(),
    m_containerChildren.end(),
    [](const decltype(m_containerChildren)::value_type& entry) { return entry.second->exists(); });

  const bool has_fields = std::any_of(
    m_fieldChildren.begin(),
    m_fieldChildren.end(),
    [](const decltype(m_fieldChildren)::value_type& entry) { return entry.second->exists(); });

  // Functions cannot be defaulted and thus have an unambiguous operator bool
  const bool has_functions = std::any_of(m_functionChildren.begin(),
                                         m_functionChildren.end(),
                                         [](const decltype(m_functionChildren)::value_type& entry) {
                                           return static_cast<bool>(*entry.second);
                                         });

  return has_containers || has_fields || has_functions;
}

bool Container::isUserProvided() const
{
  // Check if any of its child containers had user-provided fields
  const bool has_containers = std::any_of(m_containerChildren.begin(),
                                          m_containerChildren.end(),
                                          [](const decltype(m_containerChildren)::value_type& entry) {
                                            return entry.second->isUserProvided();
                                          });

  const bool has_fields = std::any_of(m_fieldChildren.begin(),
                                      m_fieldChildren.end(),
                                      [](const decltype(m_fieldChildren)::value_type& entry) {
                                        return entry.second->isUserProvided();
                                      });

  // Functions cannot be defaulted and thus have an unambiguous operator bool
  const bool has_functions = std::any_of(m_functionChildren.begin(),
                                         m_functionChildren.end(),
                                         [](const decltype(m_functionChildren)::value_type& entry) {
                                           return static_cast<bool>(*entry.second);
                                         });

  return has_containers || has_fields || has_functions;
}

bool Container::isUserProvided(const std::string& name) const
{
  if(auto container = getChildInternal<Container>(name))
  {
    // Check if the container itself was provided by the user
    return container->isUserProvided();
  }
  else if(auto field = getChildInternal<Field>(name))
  {
    // Check if the field itself was provided by the user
    return field->isUserProvided();
  }
  else if(auto function = getChildInternal<Function>(name))
  {
    // call operator bool on the function itself
    return static_cast<bool>(*function);
  }
  return false;
}

const std::unordered_map<std::string, std::unique_ptr<Container>>& Container::getChildContainers() const
{
  return m_containerChildren;
}

const std::unordered_map<std::string, std::unique_ptr<Field>>& Container::getChildFields() const
{
  return m_fieldChildren;
}

const std::unordered_map<std::string, std::unique_ptr<Function>>& Container::getChildFunctions() const
{
  return m_functionChildren;
}

}  // namespace inlet
}  // namespace axom
