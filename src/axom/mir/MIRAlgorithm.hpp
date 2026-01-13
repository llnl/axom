// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_ALGORITHM_HPP_
#define AXOM_MIR_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/ArrayView.hpp"

#include <conduit/conduit.hpp>

#include <vector>
#include <string>

namespace axom
{
namespace mir
{
/*!
 \brief Base class for Material Interface Reconstruction (MIR) algorithms.
 */
class MIRAlgorithm
{
public:
  MIRAlgorithm() = default;
  virtual ~MIRAlgorithm() = default;

  /*!
    \brief Perform material interface reconstruction on the mesh supplied by \a n_input.

    \param[in] n_input The node that contains a mesh domain with the topology and matset
                       to be used for MIR.
    \param[in] n_options A node that contains options that help govern MIR execution.

\code{.yaml}
options:
  topology: main
  matset: matset
  topologyName: mirtopo
  coordsetName: mircoords
  matsetName: cleanmat
  fields:
    - temperature
    - pressure
  selectedZones: [0,1,6,9]
  mapping: 0
\endcode
    The "topology" option specifies which topology we'll reconstruct. It must have an associated matset.
    "new_topology" is the name of the topology that will be created in the output node.
    "new_coordset" is the name of the new coordset that will be created in the output node. If it is not provided then the name of the topology's coordset will be used.  
    "new_matset" is the name of the new matset that will be created in the output node. If it is not provided then the name of the topology's matset will be used.  
    "fields" is the name of the fields to map to the new topology. If fields is specified but empty, no fields will be mapped. If fields is not present then all fields will be mapped.
    "zones" is a list of zone indices from the topology that need to be reconstructed. If not present then all zones will be considered.
    "mapping" indicates whether we should include an originalElements field on the new topology to indicate where each new zone came from in the original topology.

    \param[out] n_output A node that will contain the new entities.

    */
  virtual void execute(const conduit::Node &n_input,
                       const conduit::Node &n_options,
                       conduit::Node &n_output);

protected:
  /*!
   * \brief Set up the new domain from the old one and invoke executeDomain.
   *
   * \param n_domain The input domain.
   * \param n_options The MIR options.
   * \param n_newDomain The output domain.
   */
  void executeSetup(const conduit::Node &n_domain,
                    const conduit::Node &n_options,
                    conduit::Node &n_newDomain);

  /*!
   * \brief Perform material interface reconstruction on a single domain. Derived classes
   *        must implement this method and any device-specific coding gets handled under it.
   *
   * \param[in] n_topo The Conduit node containing the topology that will be used for MIR.
   * \param[in] n_coordset The Conduit node containing the coordset.
   * \param[in] n_fields The Conduit node containing the fields.
   * \param[in] n_matset The Conduit node containing the matset.
   * \param[in] n_options The Conduit node containing the options that help govern MIR execution.
   *
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   * \param[out] n_newMatset A Conduit node that will contain the new matset.
   * 
   */
  virtual void executeDomain(const conduit::Node &n_topo,
                             const conduit::Node &n_coordset,
                             const conduit::Node &n_fields,
                             const conduit::Node &n_matset,
                             const conduit::Node &n_options,
                             conduit::Node &n_newTopo,
                             conduit::Node &n_newCoordset,
                             conduit::Node &n_newFields,
                             conduit::Node &n_newMatset) = 0;

  /*!
   * \brief Update names in some of the objects when we can tell they have been renamed.
   *
   * \param origTopoName The topology name of the input mesh.
   * \param newTopoName The topology name of the output mesh.
   * \param origCoordsetName The coordset name of the input mesh.
   * \param newCoordsetName The coordset name of the output mesh.
   * \param origMatsetName The matset name of the input mesh.
   * \param newMatsetName The matset name of the output mesh.
   * \param n_newTopo The node that contains the output topology.
   * \param n_newCoordset The node that contains the output coordset.
   * \param n_newFields The node that contains the output fields.
   * \param n_newMatset The node that contains the output matset.
   *
   * \note This method is used internally mainly when MIR copies the input mesh to the output when MIR is no-op.
   */
  void updateNames(const std::string &origTopoName,
                   const std::string &newTopoName,
                   const std::string &origCoordsetName,
                   const std::string &newCoordsetName,
                   const std::string &origMatsetName,
                   const std::string &newMatsetName,
                   conduit::Node &n_newTopo,
                   conduit::Node &n_newCoordset,
                   conduit::Node &n_newFields,
                   conduit::Node &n_newMatset);

  /*!
   * \brief Copy state from the src domain to the destination domain.
   * \param srcState The node that contains the state in the source domain.
   * \param destState The node that contains the state in the destination domain.
   */
  void copyState(const conduit::Node &srcState, conduit::Node &destState) const;

  /*!
   * \brief This is a utility method for printing a Conduit node with large limits
   *        for lines and element counts.
   *
   * \param n The Conduit node to print.
   */
  void printNode(const conduit::Node &n) const;

  /*!
   * \brief Save a Blueprint mesh to disk (YAML and HDF5, if available).
   *
   * \param n_mesh The mesh to save.
   * \param filebase The base filename to use when writing files. Extensions may be added.
   */
  void saveMesh(const conduit::Node &n_mesh, const std::string &filebase) const;

  /*!
   * \brief Return the local path name, stripping off a domain path prefix. Blueprint domains
   *        are typically written under top level nodes with names like "domain_{:05}" or "domain_{:07}".
   *        This method will strip off any path prefix beginning with "domain".
   *
   * \param[in] obj The object whose local path we want.
   *
   * \return The path without the domain prefix.
   */
  std::string localPath(const conduit::Node &obj) const;
};

}  // end namespace mir
}  // end namespace axom

#endif
