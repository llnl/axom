// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Run.cpp
 *
 * \brief   Implementation file for Sina Run class
 *
 * \sa Record.hpp
 *
 ******************************************************************************
 */

#include "axom/sina/core/Run.hpp"

#include <utility>

#include "axom/sina/core/ConduitUtil.hpp"

namespace axom
{
namespace sina
{

namespace
{
char const RUN_TYPE[] = "run";
char const APPLICATION_FIELD[] = "application";
char const VERSION_FIELD[] = "version";
char const USER_FIELD[] = "user";
}  // namespace

Run::Run(sina::ID id, std::string application_, std::string version_, std::string user_)
  : Record {std::move(id), RUN_TYPE}
  , application {std::move(application_)}
  , version {std::move(version_)}
  , user {std::move(user_)}
{ }

Run::Run(conduit::Node const &asNode)
  : Record(asNode)
  , application {getRequiredString(APPLICATION_FIELD, asNode, RUN_TYPE)}
  , version {getOptionalString(VERSION_FIELD, asNode, RUN_TYPE)}
  , user {getOptionalString(USER_FIELD, asNode, RUN_TYPE)}
{ }

conduit::Node Run::toNode(CurveSet::CurveOrder curveOrder) const
{
  auto asNode = Record::toNode(curveOrder);
  asNode[APPLICATION_FIELD] = application;
  asNode[VERSION_FIELD] = version;
  asNode[USER_FIELD] = user;
  return asNode;
}

conduit::Node Run::toNode() const { return toNode(getDefaultCurveOrder()); }

void addRunLoader(RecordLoader &loader)
{
  loader.addTypeLoader(RUN_TYPE,
                       [](conduit::Node const &value) { return std::make_unique<Run>(value); });
}

}  // namespace sina
}  // namespace axom
