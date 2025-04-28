// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file SphinxWriter.cpp
 *
 * \brief This file contains the class implementation of the SphinxWriter.
 *******************************************************************************
 */

#include "axom/inlet/SphinxWriter.hpp"

#include <iostream>

#include "axom/slic.hpp"
#include "axom/inlet/Container.hpp"

namespace axom
{
namespace inlet
{
namespace detail
{
/**
 * \brief Determines whether a container is trivial (contains no fields/functions in its subtree)
 * 
 * \param [in] container The container to evaluate
 */
bool isTrivial(const Container& container)
{
  if(!container.getChildFields().empty() || !container.getChildFunctions().empty())
  {
    return false;
  }
  using value_type = std::decay<decltype(container.getChildContainers())>::type::value_type;
  return std::all_of(container.getChildContainers().begin(),
                     container.getChildContainers().end(),
                     [](const value_type& entry) { return isTrivial(*entry.second); });
}

/**
 * \brief Converts an enumeration to its underlying type
 * \param [in] e The enumeration value to convert
 * This function should be removed once C++23 is available
 * \see https://en.cppreference.com/w/cpp/utility/to_underlying
 */
template <typename E>
constexpr typename std::underlying_type<E>::type to_underlying(const E e)
{
  return static_cast<typename std::underlying_type<E>::type>(e);
}

}  // namespace detail

SphinxWriter::SphinxWriter(const std::string& fileName)
  : m_fieldColLabels(
      {"Field Name", "Description", "Default Value", "Range/Valid Values", "Required"})
  , m_functionColLabels({"Function Name", "Description", "Signature", "Required"})
{
  m_fileName = fileName;
  m_oss << ".. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX\n";
  m_oss << ".. |check|      unicode:: U+2611 .. CHECKED BOX\n\n";
  writeTitle("Input file Options");
}

void SphinxWriter::documentContainer(const Container& container)
{
  const auto sidreGroup = container.sidreGroup();
  const std::string pathName = sidreGroup->getPathName();
  std::string containerName = sidreGroup->getName();
  bool isSelectedElement = false;

  // If the container is empty, ignore it
  if(detail::isTrivial(container))
  {
    return;
  }

  // Replace the "implementation-defined" name with something a bit more readable
  if(isCollectionGroup(containerName))
  {
    containerName = "Collection contents:";
  }

  // If we've gotten to this point and are an element of an array/dict,
  // mark it as the selected element
  if(sidreGroup->getParent()->getName() == detail::COLLECTION_GROUP_NAME)
  {
    // The collection that this Container is a part of
    const std::string collectionName = sidreGroup->getParent()->getParent()->getPathName();
    isSelectedElement = true;
  }

  m_inletContainerPathNames.push_back(pathName);
  auto& currContainer =
    m_rstTables.emplace(pathName, ContainerData {m_fieldColLabels, m_functionColLabels}).first->second;
  currContainer.containerName = containerName;
  currContainer.isSelectedElement = isSelectedElement;
  if(containerName != "" && sidreGroup->hasView("description"))
  {
    currContainer.description = sidreGroup->getView("description")->getString();
  }

  for(const auto& field_entry : container.getChildFields())
  {
    extractFieldMetadata(field_entry.second->sidreGroup(), currContainer);
  }

  for(const auto& function_entry : container.getChildFunctions())
  {
    extractFunctionMetadata(function_entry.second->sidreGroup(), currContainer);
  }
}

void SphinxWriter::finalize()
{
  writeAllTables();
  m_outFile.open(m_fileName);
  m_outFile << m_oss.str();
  m_outFile.close();
}

void SphinxWriter::writeTitle(const std::string& title)
{
  if(title != "")
  {
    std::string equals = std::string(title.length(), '=');
    m_oss << equals << "\n" << title << "\n" << equals << "\n";
  }
}

void SphinxWriter::writeSubtitle(const std::string& sub)
{
  if(sub != "")
  {
    std::string dashes = std::string(sub.length(), '-');
    m_oss << "\n" << dashes << "\n" << sub << "\n" << dashes << "\n\n";
  }
}

void SphinxWriter::writeTable(const std::string& title,
                              const std::vector<std::vector<std::string>>& rstTable)
{
  SLIC_WARNING_IF(rstTable.size() <= 1,
                  "[Inlet] Vector for corresponding rst table must be nonempty");
  std::string result = ".. list-table:: " + title;
  std::string widths = ":widths:";
  // This would be easier with an iterator adaptor like back_inserter but for
  // concatenation
  for(std::size_t i = 0u; i < rstTable.front().size(); i++)
  {
    widths += " 25";
  }
  result += "\n   " + widths + "\n";
  result += "   :header-rows: 1\n   :stub-columns: 1\n\n";
  for(unsigned int i = 0; i < rstTable.size(); ++i)
  {
    result += "   * - ";
    for(unsigned int j = 0; j < rstTable[i].size(); ++j)
    {
      if(j != 0)
      {
        result += "     - ";
      }
      result += rstTable[i][j] + "\n";
    }
  }
  m_oss << result;
}

void SphinxWriter::writeAllTables()
{
  for(std::string& pathName : m_inletContainerPathNames)
  {
    auto& currContainer = m_rstTables.at(pathName);
    // If we're displaying a selected element, the title and description
    // will already have been printed
    if(currContainer.isSelectedElement)
    {
      m_oss << "The input schema defines a collection of this container.\n";
      m_oss << "For brevity, only one instance is displayed here.\n\n";
    }
    else
    {
      writeSubtitle(currContainer.containerName);
      if(currContainer.description != "")
      {
        m_oss << "Description: " << currContainer.description << "\n\n";
      }
    }
    if(currContainer.fieldTable.size() > 1)
    {
      writeTable("Fields", currContainer.fieldTable);
    }
    if(currContainer.functionTable.size() > 1)
    {
      writeTable("Functions", currContainer.functionTable);
    }
  }
}

std::string SphinxWriter::getValueAsString(const axom::sidre::View* view)
{
  axom::sidre::TypeID type = view->getTypeID();
  if(type == axom::sidre::TypeID::INT8_ID)
  {
    std::int8_t val = view->getData();
    return val ? "True" : "False";
  }
  else if(type == axom::sidre::TypeID::INT_ID)
  {
    int val = view->getData();
    return std::to_string(val);
  }
  else if(type == axom::sidre::TypeID::DOUBLE_ID)
  {
    double val = view->getData();
    return std::to_string(val);
  }
  return view->getString();
}

std::string SphinxWriter::getRangeAsString(const axom::sidre::View* view)
{
  std::ostringstream oss;
  oss.precision(3);
  oss << std::scientific;

  axom::sidre::TypeID type = view->getTypeID();
  if(type == axom::sidre::INT_ID)
  {
    const int* range = view->getData();
    oss << range[0] << " to " << range[1];
  }
  else
  {
    const double* range = view->getData();
    oss << range[0] << " to " << range[1];
  }
  return oss.str();
}

std::string SphinxWriter::getValidValuesAsString(const axom::sidre::View* view)
{
  const int* range = view->getData();
  size_t size = view->getBuffer()->getNumElements();
  std::string result = "";
  for(size_t i = 0; i < size; ++i)
  {
    if(i == size - 1)
    {
      result += std::to_string(range[i]);
    }
    else
    {
      result += std::to_string(range[i]) + ", ";
    }
  }
  return result;
}

std::string SphinxWriter::getValidStringValues(const axom::sidre::Group* sidreGroup)
{
  auto idx = sidreGroup->getFirstValidViewIndex();
  std::string validValues = "";
  while(axom::sidre::indexIsValid(idx))
  {
    validValues += std::string(sidreGroup->getView(idx)->getString());
    idx = sidreGroup->getNextValidViewIndex(idx);
    if(axom::sidre::indexIsValid(idx))
    {
      validValues += ", ";
    }
  }
  return validValues;
}

void SphinxWriter::extractFieldMetadata(const axom::sidre::Group* sidreGroup,
                                        ContainerData& currentContainer)
{
  std::vector<std::string> fieldAttributes(m_fieldColLabels.size());

  fieldAttributes[0] = sidreGroup->getName();

  if(sidreGroup->hasView("description"))
  {
    fieldAttributes[1] = std::string(sidreGroup->getView("description")->getString());
  }

  if(sidreGroup->hasView("defaultValue"))
  {
    fieldAttributes[2] = getValueAsString(sidreGroup->getView("defaultValue"));
  }

  if(sidreGroup->hasView("range"))
  {
    fieldAttributes[3] = getRangeAsString(sidreGroup->getView("range"));
  }
  else if(sidreGroup->hasView("validValues"))
  {
    fieldAttributes[3] = getValidValuesAsString(sidreGroup->getView("validValues"));
  }
  else if(sidreGroup->hasGroup("validStringValues"))
  {
    fieldAttributes[3] = getValidStringValues(sidreGroup->getGroup("validStringValues"));
  }

  if(sidreGroup->hasView("required"))
  {
    std::int8_t required = sidreGroup->getView("required")->getData();
    fieldAttributes[4] = required ? "|check|" : "|uncheck|";
  }
  else
  {
    fieldAttributes[4] = "|uncheck|";
  }

  currentContainer.fieldTable.push_back(fieldAttributes);
}

std::string SphinxWriter::getSignatureAsString(const axom::sidre::Group* sidreGroup)
{
  using underlying = std::underlying_type<FunctionTag>::type;
  static const auto type_names = []() {
    std::unordered_map<underlying, std::string> result;
    result[detail::to_underlying(FunctionTag::Vector)] = "Vector";
    result[detail::to_underlying(FunctionTag::Double)] = "Double";
    result[detail::to_underlying(FunctionTag::Void)] = "Void";
    result[detail::to_underlying(FunctionTag::String)] = "String";
    return result;
  }();

  // View::getData<T> does not have a const version...
  const auto ret_type = static_cast<underlying>(sidreGroup->getView("return_type")->getData());

  const auto args_view = sidreGroup->getView("function_arguments");
  const underlying* arg_tags = args_view->getData();
  const int num_args = args_view->getNumElements();
  std::vector<std::string> arg_types(num_args);
  for(int i = 0; i < num_args; i++)
  {
    arg_types[i] = type_names.at(arg_tags[i]);
  }
  return fmt::format("{0}({1})", type_names.at(ret_type), fmt::join(arg_types, ", "));
}

void SphinxWriter::extractFunctionMetadata(const axom::sidre::Group* sidreGroup,
                                           ContainerData& currentContainer)
{
  std::vector<std::string> functionAttributes(m_functionColLabels.size());

  functionAttributes[0] = sidreGroup->getName();

  if(sidreGroup->hasView("description"))
  {
    functionAttributes[1] = std::string(sidreGroup->getView("description")->getString());
  }

  functionAttributes[2] = getSignatureAsString(sidreGroup);

  if(sidreGroup->hasView("required"))
  {
    std::int8_t required = sidreGroup->getView("required")->getData();
    functionAttributes[3] = required ? "|check|" : "|uncheck|";
  }
  else
  {
    functionAttributes[3] = "|uncheck|";
  }

  currentContainer.functionTable.push_back(functionAttributes);
}

}  // namespace inlet
}  // namespace axom
