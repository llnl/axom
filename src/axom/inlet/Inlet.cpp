// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Inlet.cpp
 *
 * \brief This file contains the class implementation of Inlet, the main class
 *        for the Inlet component.
 *******************************************************************************
 */

#include "axom/inlet/Inlet.hpp"

#include "axom/fmt.hpp"
#include "axom/slic.hpp"
#include "axom/inlet/inlet_utils.hpp"
#include <algorithm>

namespace axom
{
namespace inlet
{
Container& Inlet::addStruct(const std::string& name, const std::string& description)
{
  return m_globalContainer.addStruct(name, description);
}

VerifiableScalar& Inlet::addBool(const std::string& name, const std::string& description)
{
  return m_globalContainer.addBool(name, description);
}

VerifiableScalar& Inlet::addDouble(const std::string& name, const std::string& description)
{
  return m_globalContainer.addDouble(name, description);
}

VerifiableScalar& Inlet::addInt(const std::string& name, const std::string& description)
{
  return m_globalContainer.addInt(name, description);
}

VerifiableScalar& Inlet::addString(const std::string& name, const std::string& description)
{
  return m_globalContainer.addString(name, description);
}

namespace detail
{
/*!
 *******************************************************************************
 * \brief Recursive helper function for traversing an Inlet tree for documentation
 * generation purposes
 * 
 * \param [inout] writer The Writer to use for documentation
 * \param [in] container The current container to write
 *******************************************************************************
 */
void writerHelper(Writer& writer, const Container& container)
{
  // Use a pre-order traversal for readability
  writer.documentContainer(container);
  // Only visit a single element of a *struct* collection
  if(isCollectionGroup(container.name()) &&
     container.sidreGroup()->hasView(detail::STRUCT_COLLECTION_FLAG))
  {
    auto indices = detail::collectionIndices(container);
    // Just use the first index
    if(!indices.empty())
    {
      writerHelper(
        writer,
        *container.getChildContainers().at(
          utilities::string::appendPrefix(container.name(), detail::indexToString(indices[0]))));
    }
  }
  else
  {
    for(const auto& sub_container_entry : container.getChildContainers())
    {
      writerHelper(writer, *sub_container_entry.second);
    }
  }
}

}  // end namespace detail

void Inlet::write(Writer&& writer)
{
  if(m_docEnabled)
  {
    detail::writerHelper(writer, m_globalContainer);
    writer.finalize();
  }
}

bool Inlet::verify(std::vector<VerificationError>* errors) const
{
  return m_globalContainer.verify(errors);
}

}  // end namespace inlet
}  // end namespace axom
