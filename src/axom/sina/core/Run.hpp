// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_RUN_HPP
#define SINA_RUN_HPP

/*!
 ******************************************************************************
 *
 * \file Run.hpp
 *
 * \brief   Header file for Sina Run class
 *
 * \sa Record.hpp
 *
 ******************************************************************************
 */

#include "axom/sina/core/Record.hpp"

namespace axom
{
namespace sina
{

/**
 * \brief A sub-type of Record representing a single run of an applicaiton
 *
 * A Run is a subtype of Record corresponding to a single run of an application, as
 * specified in the Sina schema. A Run has a few additional fields required in addition
 * to the id required by a Record (type is automatically set to "run"):
 *
 *  - application: the application/code used to create the Run
 *  - version: the version of the application used to create the Run
 *  - user: the username of the person who ran the application that generated this Run
 *
 * To create a Run:
 * \code
 *   axom::sina::ID run1ID{"run1", axom::sina::IDType::Local};
 *   std::unique_ptr<axom::sina::Record> run1{new axom::sina::Run{run1ID, "My Sim Code", "1.2.3", "jdoe"}};
 * \endcode
 *
 */
class Run : public Record
{
public:
  /**
     * \brief Create a new Run.
     *
     * \param id the run's ID
     * \param application the application that was run
     * \param version (optional) the version of the application
     * \param user (optional) the user who executed the run
     */
  Run(ID id, std::string application, std::string version = "", std::string user = "");

  /**
     * \brief Create a Run from its representation as a conduit Node
     *
     * \param asNode the run as a Node
     */
  explicit Run(conduit::Node const &asNode);

  /**
     * \brief Get the application that was run.
     *
     * \return the application's name
     */
  std::string const &getApplication() const { return application; }

  /**
     * \brief Get the version of the application that was run.
     *
     * \return the application's version
     */
  std::string const &getVersion() const { return version; }

  /**
     * \brief Get the name of the user who ran the application.
     *
     * \return the user's name
     */
  std::string const &getUser() const { return user; }

  conduit::Node toNode(CurveSet::CurveOrder curveOrder) const override;
  conduit::Node toNode() const;

private:
  std::string application;
  std::string version;
  std::string user;
};

/**
 * \brief Add a type loader to the given RecordLoader for loading Run instances.
 *
 * \param loader the RecordLoader to which to add the function for loading
 *               Run instances.
 */
void addRunLoader(RecordLoader &loader);

}  // namespace sina
}  // namespace axom

#endif  //SINA_RUN_HPP
