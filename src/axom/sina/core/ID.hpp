// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_ID_HPP
#define SINA_ID_HPP

/*!
 ******************************************************************************
 *
 * \file ID.hpp
 *
 * \brief   Header file for Sina ID class
 *
 * The Sina schema allows records to have either a local ID or a global ID.
 * When a global ID is specified, that will be used in the database. When a
 * local ID is specified, an ID will be automatically generated when inserting
 * the record into the database.
 *
 ******************************************************************************
 */

#include <string>

#include "conduit.hpp"

namespace axom
{
namespace sina
{

/**
 * Represents whether an ID is local or global.
 */
enum class IDType
{
  Local,
  Global
};

/**
 * \brief The ID of a Record
 *
 * An ID is used to represent the ID of a Record. This class holds both the
 * actual ID and whether it is a local or global ID.
 */
class ID
{
public:
  /**
     * \brief Create a new ID.
     *
     * \param id the actual value of the ID
     * \param type whether the ID is local or global
     */
  ID(std::string id, IDType type);

  /**
     * \brief Get the value of the ID.
     *
     * \return the actual ID
     */
  std::string const &getId() const noexcept { return id; }

  /**
     * \brief Get the type of the ID.
     *
     * \return whether the ID is local or global
     */
  IDType getType() const noexcept { return type; }

private:
  std::string id;
  IDType type;
};

namespace internal
{

/**
 * \brief An object representing a pair of ID fields
 *
 * IDField instances are used to describe a pair of ID fields in a schema
 * object which correspond to global and local names for the field. For
 * example, the "id" and "local_id" fields in records, or the
 * "subject"/"local_subject" and "object"/"local_object" pairs in
 * relationships.
 */
class IDField
{
public:
  /**
     * \brief Construct a new IDField.
     *
     * \param value the value of the ID
     * \param localName the name of the local variant of the field
     * \param globalName the name of the global variant of the field
     */
  IDField(ID value, std::string localName, std::string globalName);

  /**
     * \brief Construct an IDField by looking for its values in a conduit Node.
     *
     * \param parentObject the conduit Node containing the ID field
     * \param localName the local name of the field
     * \param globalName the global name of the field
     */
  IDField(conduit::Node const &parentObject, std::string localName, std::string globalName);

  /**
     * \brief Get the value of this field.
     *
     * \return the ID describing the field's value
     */
  ID const &getID() const noexcept { return value; }

  /**
     * \brief Get the name to use for this field when the ID is local.
     *
     * \return the name of the local ID field
     */
  std::string const &getLocalName() const noexcept { return localName; }

  /**
     * \brief Get the name to use for this field when the ID is global.
     *
     * \return the name of the global ID field
     */
  std::string const &getGlobalName() const noexcept { return globalName; }

  /**
     * \brief Add this field to the given Node.
     *
     * \param object the Node to which to add the field
     */
  void addTo(conduit::Node &object) const;

private:
  ID value;
  std::string localName;
  std::string globalName;
};

}  // namespace internal
}  // namespace sina
}  // namespace axom

#endif  //SINA_ID_HPP
