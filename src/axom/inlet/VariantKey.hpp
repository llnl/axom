// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file VariantKey.hpp
 *
 * \brief This file contains the class definition of VariantKey, Inlet's generic
 * associative array key type
 *******************************************************************************
 */

#ifndef INLET_KEY_HPP
#define INLET_KEY_HPP

#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \enum InletType
 *
 * \brief Enumeration of basic types for things in inlet
 *******************************************************************************
 */
enum class InletType
{
  Nothing,
  Bool,
  String,
  Integer,
  // TODO: Unsigned integer
  Double,
  Object,
  Collection,
  Function
};

/*!
 *******************************************************************************
 * \class VariantKey
 *
 * \brief Provides a generic key type for mixed-key associative arrays/dictionaries
 * 
 * \note Can be thought of as limited implementation of std::variant<int, std::string>
 *******************************************************************************
 */
class VariantKey
{
public:
  /*!
   *****************************************************************************
   * \brief Parameterized constructors for initializing the variant
   * \param [in] key The key to initialize with 
   *****************************************************************************
   */
  VariantKey(const int key) : m_int(key), m_type(VariantKeyType::Integer) { }
  /// \overload
  VariantKey(const std::string& key) : m_string(key), m_type(VariantKeyType::String) { }
  /// \overload
  VariantKey(std::string&& key) : m_string(std::move(key)), m_type(VariantKeyType::String) { }
  /// \overload
  // Explicit const char[] overload required to convert to std::string instead of ptr -> int
  VariantKey(const char key[]) : m_string(key), m_type(VariantKeyType::String) { }

  /*!
   *****************************************************************************
   * \brief Parameterized assignment operators for re-initializing the variant
   * \param [in] key The key to initialize with 
   *****************************************************************************
   */
  VariantKey& operator=(const int key);
  /// \overload
  VariantKey& operator=(const std::string& key);
  /// \overload
  VariantKey& operator=(std::string&& key);
  /// \overload
  VariantKey& operator=(const char key[]);

  /*!
   *****************************************************************************
   * \brief Implicit conversion operators to make usage more convenient
   *****************************************************************************
   */
  operator int() const;
  operator const std::string&() const;

  /*!
   *****************************************************************************
   * \brief Returns the type of the active member
   *****************************************************************************
   */
  InletType type() const;

  /*!
   *****************************************************************************
   * \brief Comparison operator, returns true iff the active types are the same
   * and the corresponding active members compare equal
   *****************************************************************************
   */
  bool operator==(const VariantKey& other) const;

private:
  /*!
   *****************************************************************************
   * \brief Subset of InletType containing only the types present in the union
   *****************************************************************************
   */
  enum class VariantKeyType
  {
    Integer,
    String
  };

  // Integer and string keys
  // With only two possible types a union is overkill
  int m_int {};
  std::string m_string;

  // Active key type
  VariantKeyType m_type;
};

/*!
 *******************************************************************************
 * \brief Inserts the key into a stream
 *
 * \param [inout] out The stream to insert into
 * \param [in] key The key to print
 * 
 * Prints the active member to the stream
 *******************************************************************************
 */
std::ostream& operator<<(std::ostream& out, const VariantKey& key);

}  // end namespace inlet
}  // end namespace axom

namespace std
{
/*!
 *******************************************************************************
 * \brief Specialization of std::hash so VariantKey is compatible with
 * std::unordered_map
 *******************************************************************************
 */
template <>
struct hash<axom::inlet::VariantKey>
{
  size_t operator()(const axom::inlet::VariantKey& key) const
  {
    if(key.type() == axom::inlet::InletType::Integer)
    {
      return std::hash<int> {}(key);
    }
    else
    {
      return std::hash<std::string> {}(key);
    }
  }
};
}  // end namespace std

#endif  // INLET_KEY_HPP
