// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file AttrValues.cpp
 *
 * \brief   Implementation file for AttrValues class.
 *
 ******************************************************************************
 */

// Associated header file
#include "AttrValues.hpp"

// Other axom headers
#include "axom/core/Types.hpp"
#include "axom/slic/interface/slic.hpp"

// Sidre component headers

namespace axom
{
namespace sidre
{
/*
 *************************************************************************
 *
 * Return true if the attribute has been explicitly set; else false.
 *
 *************************************************************************
 */
bool AttrValues::hasValue(const Attribute* attr) const
{
  // check for empty values
  if(attr == nullptr || m_values == nullptr)
  {
    return false;
  }

  const auto iattr = static_cast<std::size_t>(attr->getIndex());
  if(iattr >= m_values->size())
  {
    // This attribute has not been set for this View.
    return false;
  }

  Node& value = (*m_values)[iattr];
  return !isEmpty(value);
}

/*
 *************************************************************************
 *
 * Set Attribute to its default value. If the Attribute has not been
 * set yet, then it is already the default value.  Otherwise, reset
 * the Node to empty.  Future calls to hasAttributeValue will return
 * false for the attribute.
 *
 *************************************************************************
 */
bool AttrValues::setToDefault(const Attribute* attr)
{
  if(attr == nullptr)
  {
    return false;
  }

  if(m_values == nullptr)
  {
    // No attributes have been set in this View, already default.
    return true;
  }

  const auto iattr = static_cast<std::size_t>(attr->getIndex());
  if(iattr >= m_values->size())
  {
    // This attribute has not been set for this View, already default.
    return true;
  }

  Node& value = (*m_values)[iattr];
  value.reset();

  return true;
}

/*
 *************************************************************************
 *
 * PRIVATE Create a Node for the Attribute.
 *
 * Create vector of Nodes and push empty nodes up to attr's index.
 *
 *************************************************************************
 */
bool AttrValues::createNode(IndexType iattr)
{
  if(m_values == nullptr)
  {
    m_values = new(std::nothrow) Values();
  }

  if(static_cast<std::size_t>(iattr) >= m_values->size())
  {
    // Create all attributes up to iattr, push back empty Nodes
    m_values->reserve(iattr + 1);
    for(int n = m_values->size(); n < iattr + 1; ++n)
    {
      m_values->emplace_back(Node());
    }
  }

  return true;
}

/*
 *************************************************************************
 *
 * Return a scalar attribute value.
 *
 *************************************************************************
 */
Node::ConstValue AttrValues::getScalar(const Attribute* attr) const
{
  const Node& node = getValueNodeRef(attr);
  return node.value();
}

/*
 *************************************************************************
 *
 * Return a string attribute value.
 *
 * The caller must ensure that attr is not NULL.
 *
 *************************************************************************
 */
const char* AttrValues::getString(const Attribute* attr) const
{
  if(attr == nullptr)
  {
    return nullptr;
  }

  if(attr->getTypeID() != CHAR8_STR_ID)
  {
    SLIC_CHECK_MSG(attr->getTypeID() == CHAR8_STR_ID,
                   "getString: Called on attribute '" << attr->getName() << "' which is type "
                                                      << DataType::id_to_name(attr->getTypeID())
                                                      << ".");
    return nullptr;
  }

  const Node& node = getValueNodeRef(attr);
  return node.as_char8_str();
}

/*
 *************************************************************************
 *
 * Return Reference to Node.
 *
 *************************************************************************
 */
const Node& AttrValues::getValueNodeRef(const Attribute* attr) const
{
  if(attr == nullptr)
  {
    return getEmptyNodeRef();
  }

  if(m_values == nullptr)
  {
    // No attributes have been set in this View;
    return attr->getDefaultNodeRef();
  }

  const auto iattr = static_cast<std::size_t>(attr->getIndex());
  if(iattr >= m_values->size())
  {
    // This attribute has not been set for this View
    return attr->getDefaultNodeRef();
  }

  Node& value = (*m_values)[iattr];

  if(isEmpty(value))
  {
    return attr->getDefaultNodeRef();
  }

  return value;
}

/*
 *************************************************************************
 * Return first valid Attribute index for a set Attribute.
 * (i.e., smallest index over all Attributes).
 *
 * sidre::InvalidIndex is returned if DataStore has no Attributes.
 *************************************************************************
 */
IndexType AttrValues::getFirstValidAttrValueIndex() const
{
  if(m_values == nullptr)
  {
    // No attributes have been set in this View.
    return InvalidIndex;
  }

  for(size_t iattr = 0; iattr < m_values->size(); ++iattr)
  {
    // Find first, non-empty attribute.
    Node& value = (*m_values)[iattr];
    if(!isEmpty(value))
    {
      return iattr;
    }
  }

  return InvalidIndex;
}

/*
 *************************************************************************
 * Return next valid Attribute index for a set Attribute after given
 * index (i.e., smallest index over all Attribute indices larger than
 * given one).
 *
 * sidre::InvalidIndex is returned if there is no valid index greater
 * than given one.
 * getNextAttrValueIndex(InvalidIndex) returns InvalidIndex.
 *************************************************************************
 */
IndexType AttrValues::getNextValidAttrValueIndex(IndexType idx) const
{
  if(idx == InvalidIndex)
  {
    return InvalidIndex;
  }

  idx++;
  while(static_cast<unsigned>(idx) < m_values->size() && isEmpty((*m_values)[idx]))
  {
    idx++;
  }
  return ((static_cast<unsigned>(idx) < m_values->size()) ? idx : InvalidIndex);
}

/*
 *************************************************************************
 *
 * PRIVATE ctor for AttrValues
 *
 *************************************************************************
 */
AttrValues::AttrValues() : m_values(nullptr) { }

/*
 *************************************************************************
 *
 * PRIVATE dtor.
 *
 *************************************************************************
 */
AttrValues::~AttrValues()
{
  if(m_values != nullptr)
  {
    delete m_values;
    m_values = nullptr;
  }
}

} /* end namespace sidre */
} /* end namespace axom */
