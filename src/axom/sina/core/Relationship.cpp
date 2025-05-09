// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Relationship.cpp
 *
 * \brief   Implementation file for Sina Relationship class
 *
 ******************************************************************************
 */

#include "axom/sina/core/Relationship.hpp"

#include <utility>

#include "axom/sina/core/ConduitUtil.hpp"

namespace axom
{
namespace sina
{

namespace
{
char const GLOBAL_SUBJECT_KEY[] = "subject";
char const LOCAL_SUBJECT_KEY[] = "local_subject";
char const GLOBAL_OBJECT_KEY[] = "object";
char const LOCAL_OBJECT_KEY[] = "local_object";
char const PREDICATE_KEY[] = "predicate";
}  // namespace

Relationship::Relationship(ID subject_, std::string predicate_, ID object_)
  : subject {std::move(subject_), LOCAL_SUBJECT_KEY, GLOBAL_SUBJECT_KEY}
  , object {std::move(object_), LOCAL_OBJECT_KEY, GLOBAL_OBJECT_KEY}
  , predicate {std::move(predicate_)}
{ }

Relationship::Relationship(conduit::Node const &asNode)
  : subject {asNode, LOCAL_SUBJECT_KEY, GLOBAL_SUBJECT_KEY}
  , object {asNode, LOCAL_OBJECT_KEY, GLOBAL_OBJECT_KEY}
  , predicate {getRequiredString(PREDICATE_KEY, asNode, "Relationship")}
{ }

conduit::Node Relationship::toNode() const
{
  conduit::Node relationshipNode;
  relationshipNode[PREDICATE_KEY] = predicate;
  subject.addTo(relationshipNode);
  object.addTo(relationshipNode);
  return relationshipNode;
}

}  // namespace sina
}  // namespace axom
