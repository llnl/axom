// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file AttrValue.hpp
 *
 * \brief   Header file containing definition of AttrValue class.
 *
 ******************************************************************************
 */

// Standard C++ headers
#include <string>
#include <vector>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/slic/interface/slic.hpp"

// Sidre project headers
#include "axom/sidre/core/Attribute.hpp"
#include "axom/sidre/core/SidreTypes.hpp"

#ifndef SIDRE_ATTRVALUES_HPP_
  #define SIDRE_ATTRVALUES_HPP_

namespace axom
{
namespace sidre
{
class View;

/*!
 * \class AttrValue
 *
 * \brief Store Attribute values.
 *
 * Each attribute is defined by an instance of Attribute in the DataStore.
 * The attribute has an associated type and index.
 *
 * The assumption is made that attributes will be looked up more often than
 * set.  So getAttribute should be optimized over setAttribute.
 *
 * Another assumption is that the View should not pay for attributes if they
 * are not being used.
 *
 * Space can be minimized by creating more common Attribute first so that they
 * will have a lower index.
 *
 * Methods which accept a Attribute pointer will check for
 * nullptr but will not print a message.  The calling routines in
 * the View class are expected to print any error message.  This will
 * avoid multiple messages for the same error.  For example, if the
 * index cannot be converted to an Attribute pointer in the View
 * class, an error message will be printing and then a NULL pointer
 * passed to the AttrValues class which will not print another message.
 */
class AttrValues
{
public:
  /// Friend declarations to constrain usage via controlled access to private members.
  friend class View;

private:
  //DISABLE_DEFAULT_CTOR(AttrValues);
  DISABLE_MOVE_AND_ASSIGNMENT(AttrValues);

  /*!
   * \brief Return true if attribute is set in value.
   */
  static bool isEmpty(Node& value) { return value.schema().dtype().is_empty(); }

  /*!
   * \brief Return true if the attribute has been explicitly set; else false.
   */
  bool hasValue(const Attribute* attr) const;

  /*!
   * \brief Create a Conduit Node to store an attribute.
   *
   * Create vector of Nodes and push empty nodes up to attr's index.
   * Called as part of View::createAttributeScalar and
   * View::createAttributeString.
   *
   * \return true
   */
  bool createNode(IndexType idx);

  /*!
   * \brief Set attribute to its default value.
   */
  bool setToDefault(const Attribute* attr);

  /*!
   * \brief Set attribute value from a scalar.
   */
  template <typename ScalarType>
  bool setScalar(const Attribute* attr, ScalarType value)
  {
    const DataTypeId arg_id = detail::SidreTT<ScalarType>::id;
    if(arg_id != attr->getTypeID())
    {
      SLIC_CHECK_MSG(arg_id == attr->getTypeID(),
                     "setScalar: Incorrect type for attribute '"
                       << attr->getName() << "' of type " << attr->getDefaultNodeRef().dtype().name()
                       << ": " << DataType::id_to_name(arg_id) << ".");
      return false;
    }

    const IndexType iattr = attr->getIndex();
    const bool ok = createNode(iattr);
    if(ok)
    {
      (*m_values)[iattr] = value;
    }
    return ok;
  }

  /*!
   * \brief Set attribute value from a string.
   */
  bool setString(const Attribute* attr, const std::string& value)
  {
    const DataTypeId arg_id = CHAR8_STR_ID;
    if(arg_id != attr->getTypeID())
    {
      SLIC_CHECK_MSG(arg_id == attr->getTypeID(),
                     "setString: Incorrect type for attribute '"
                       << attr->getName() << "' of type " << attr->getDefaultNodeRef().dtype().name()
                       << ": " << DataType::id_to_name(arg_id) << ".");
      return false;
    }

    const IndexType iattr = attr->getIndex();
    const bool ok = createNode(iattr);
    if(ok)
    {
      (*m_values)[iattr] = value;
    }
    return ok;
  }

  /*!
   * \brief Set attribute value from a Node.
   *
   * Used when restoring attributes from a file.
   * The type of node is not check.
   */
  bool setNode(const Attribute* attr, const Node& node)
  {
    const IndexType iattr = attr->getIndex();
    const bool ok = createNode(iattr);
    if(ok)
    {
      (*m_values)[iattr] = node;
    }
    return ok;
  }

  /*!
   * \brief Remove all attribute values.
   */
  void clear()
  {
    if(m_values != nullptr)
    {
      delete m_values;
      m_values = nullptr;
    }
  }

  /*!
   * \brief Return a scalar attribute value.
   */
  Node::ConstValue getScalar(const Attribute* attr) const;

  /*!
   * \brief Return a string attribute value.
   */
  const char* getString(const Attribute* attr) const;

  /*!
   * \brief Return reference to value Node.
   */
  const Node& getValueNodeRef(const Attribute* attr) const;

  /*!
   * \brief Return a reference to an empty Node.
   *
   * Used as error return value from getValueNodeRef.
   */
  const Node& getEmptyNodeRef() const
  {
    static const Node empty;
    return empty;
  }

  /*!
   * \brief Return first valid Attribute index for a set Attribute.
   *        (i.e., smallest index over all Attributes).
   *
   * sidre::InvalidIndex is returned if AttrValue has no Attributes.
   */
  IndexType getFirstValidAttrValueIndex() const;

  /*!
   * \brief Return next valid Attribute index for a set Attribute
   *        after given index (i.e., smallest index over all Attribute
   *        indices larger than given one).
   *
   * sidre::InvalidIndex is returned if there is no valid index greater
   * than given one.
   * getNextAttrValueIndex(InvalidIndex) returns InvalidIndex.
   */
  IndexType getNextValidAttrValueIndex(IndexType idx) const;

  ///@{
  //!  @name Private AttrValues ctor and dtor
  //!        (callable only by DataStore methods).

  /*!
   *  \brief Private ctor.
   */
  AttrValues();

  /// \brief Private copy ctor.
  AttrValues(const AttrValues& source);

  /// \brief Private dtor.
  ~AttrValues();

  ///@}

  ///////////////////////////////////////////////////////////////////
  using Values = std::vector<Node>;
  ///////////////////////////////////////////////////////////////////

  /// Attributes values.
  Values* m_values;
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_ATTRVALUES_HPP_ */
