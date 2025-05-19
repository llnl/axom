// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EXAMPLES_MIR_APPLICATION_HPP
#define AXOM_MIR_EXAMPLES_MIR_APPLICATION_HPP
#include "axom/config.hpp"
#include "axom/core.hpp"  // for axom macros

#include <conduit.hpp>

#include <string>

/*!
 * \brief The concentric circles application.
 */
class MIRApplication
{
public:
  /// Constructor
  MIRApplication();

  /*!
   * \brief Initialize the application from command line args.
   * \return 0 on success; less than zero otherwise.
   */
  int initialize(int argc, char **argv);

  /*!
   * \brief Execute the main application logic.
   * \return 0 on success; less than zero otherwise.
   */
  int execute();

protected:
  /*!
   * \brief Returns whether a structured mesh is needed.
   * \return True if structured mesh is needed; false otherwise.
   */
  bool requiresStructuredMesh(const std::string &method) const;

  /*!
   * \brief Invoke the MIR appropriate for the selected runtime policy.
   * \return 0 on success; less than zero otherwise.
   */
  int runMIR();

  /*!
   * \brief Make any adjustments to the mesh.
   */
  virtual void adjustMesh(conduit::Node &);

  /*!
   * \brief Save the mesh to a file.
   *
   * \param path The filepath where the file will be saved.
   * \param n_mesh The mesh to be saved.
   */
  virtual void saveMesh(const conduit::Node &n_mesh, const std::string &path);

  /*!
   * \brief A static error handler for Conduit.
   */
  static void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1);

  bool handler;
  int gridSize;
  int numCircles;
  int dimension;
  int numTrials;
  bool writeFiles;
  std::string outputFilePath;
  std::string method;
  axom::runtime_policy::Policy policy;
  std::string annotationMode;
  std::string protocol;
};

#endif
