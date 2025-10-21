// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Attribute.hpp
 *
 * \brief   Header file containing definition of Attribute class.
 *
 ******************************************************************************
 */

#ifndef SIDRE_ATTRIBUTE_HPP_
#define SIDRE_ATTRIBUTE_HPP_

// Standard C++ headers
#include <string>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/slic/interface/slic.hpp"

// Sidre project headers
#include "axom/sidre/core/SidreTypes.hpp"

namespace axom
{
namespace sidre
{
/*!
 * \class Attribute
 *
 * \brief An Attribute holds scalar metadata describing a View.
 *
 */
class Attribute
{
public:
  /*!
   * Friend declarations to constrain usage via controlled access to
   * private members.
   */
  friend class DataStore;

  /*!
   * \brief Return const reference to name of Attribute object.
   */
  const std::string& getName() const { return m_name; }

  /*!
   * \brief Return the unique index of this Attribute object.
   */
  IndexType getIndex() const { return m_index; }

  /*!
   * \brief Set default value of Attribute. Return true if successfully changed.
   *
   * The type of the default cannot be changed after the attribute is created.
   */
  template <typename ScalarType>
  bool setDefaultScalar(ScalarType value)
  {
    DataTypeId arg_id = detail::SidreTT<ScalarType>::id;
    if(m_default_value.dtype().is_empty() || arg_id == m_default_value.dtype().id())
    {
      m_default_value = value;
      return true;
    }
    else
    {
      SLIC_CHECK_MSG(arg_id == m_default_value.dtype().id(),
                     "setDefaultScalar: Cannot change type of attribute '"
                       << m_name << "' from " << m_default_value.dtype().name() << " to "
                       << DataType::id_to_name(arg_id) << ".");
      return false;
    }
  }

  /*!
   * \brief Set default value of Attribute. Return true if successfully changed.
   *
   * The type of the default cannot be changed after the attribute is created.
   */
  bool setDefaultString(const std::string& value)
  {
    DataTypeId arg_id = CHAR8_STR_ID;
    if(m_default_value.dtype().is_empty() || arg_id == m_default_value.dtype().id())
    {
      m_default_value = value;
      return true;
    }
    else
    {
      SLIC_CHECK_MSG(arg_id == m_default_value.dtype().id(),
                     "setDefaultString: Cannot change type of attribute '"
                       << m_name << "' from " << m_default_value.dtype().name() << " to "
                       << DataType::id_to_name(arg_id) << ".");
      return false;
    }
  }

  /*!
   * \brief Set default value of Attribute to a Node.
   */
  bool setDefaultNodeRef(Node& node)
  {
    m_default_value = node;
    return true;
  }

  /*!
   * \brief Return default value of Attribute.
   */
  const Node& getDefaultNodeRef() const { return m_default_value; }

  /*!
   * \brief Return type of Attribute.
   */
  TypeID getTypeID() const { return static_cast<TypeID>(m_default_value.dtype().id()); }

private:
  DISABLE_DEFAULT_CTOR(Attribute);
  DISABLE_MOVE_AND_ASSIGNMENT(Attribute);

  ///@{
  //!  @name Private Attribute ctor and dtor
  //!        (callable only by DataStore methods).

  /*!
   *  \brief Private ctor that creates an Attribute with given name
   *         which has no data associated with it.
   */
  Attribute(const std::string& name);

  /*!
   * \brief Private copy ctor.
   */
  Attribute(const Attribute& source);

  /*!
   * \brief Private dtor.
   */
  ~Attribute();

  ///@}

  /// Name of this Attribute object.
  std::string m_name;

  /// Attribute's unique index within DataStore object that created it.
  IndexType m_index;

  /// Default value of attribute.
  Node m_default_value;
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_ATTRIBUTE_HPP_ */
