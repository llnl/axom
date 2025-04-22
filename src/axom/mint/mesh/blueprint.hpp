// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_MESH_BLUEPRINT_HPP_
#define MINT_MESH_BLUEPRINT_HPP_

// mint includes
#include "axom/mint/config.hpp"          // for compile-time definitions
#include "axom/mint/mesh/CellTypes.hpp"  // for Topology

// C/C++ includes
#include <string>  // for std::string

/*!
 * \file
 *
 * \brief Provides convenience methods for navigating and querying a
 *  Sidre group hierarchy according to the conventions outlined in the
 *  <a href="http://llnl-conduit.readthedocs.io/en/latest/"> mesh blueprint</a>.
 *
 */

namespace axom
{
// Sidre Forward Declarations
namespace sidre
{
class Group;
}

namespace mint
{
namespace blueprint
{
#ifdef AXOM_MINT_USE_SIDRE

/*!
 * \brief Checks if the given root group conforms to the mesh blueprint.
 *
 * \param [in] group pointer to the blueprint root group.
 * \return status true if the group conforms to the blueprint, else, false.
 *
 * \note A root group conforms to the computational mesh blueprint if it
 *  satisfies the following conditions:
 *  (a) it has a "coordsets" child group
 *  (b) it has a "topologies" child group
 *  (c) it has a "fields" child group
 *
 * \note This method emits warnings for each missing group it encounters.
 *
 * \pre group != nullptr
 */
bool isValidRootGroup(const sidre::Group* group);

/*!
 * \brief Checks if the given topology group conforms to the blueprint.
 *
 * \param [in] topo pointer to a topology group.
 * \return status true if the group conforms to the blueprint, else, false.
 *
 * \note This routine performs cursory mesh-agnostic checks, i.e., checks that
 *  the necessary child views exists etc., but, it does not validate the content
 *  of those view at this level.
 *
 * \note A valid topology group conforms to the computational mesh blueprint if
 *  it satisfies the following conditions:
 *  (a) it has a "type" view
 *  (b) it has a "coordset" view
 *  (c) the "type" and "coordset" views are both string views.
 *
 * \pre topo != nullptr
 */
bool isValidTopologyGroup(const sidre::Group* topo);

/*!
 * \brief Checks if the given coordset group conforms to the blueprint.
 *
 * \param [in] coordset pointer to a coordset group.
 * \return status true if the group conforms to the blueprint, else, false.
 *
 * \note This routine performs cursory mesh-agnostic checks, i.e., checks that
 *  the necessary child views exists etc., but, it does not validate the content
 *  of those view at this level.
 *
 * \note A valid coordset group conforms to the computational mesh blueprint
 *  if  it satisfies the following conditions:
 *  (a) it has a "type" view
 *  (b) the "type" view is a string view
 *
 * \pre coordset != nullptr
 */
bool isValidCoordsetGroup(const sidre::Group* coordset);

/*!
 * \brief Returns the requested mesh topology group.
 *
 * \param [in] group the blueprint mesh root group.
 * \param [in] topo optional argument to specify the topology name.
 *
 * \note If a topology is not specified with the second argument, the
 *  code returns the first group under the topologies group.
 *
 * \return topo pointer to the requested topology group or nullptr, if
 *  the requested topology group does not exists or is invalid.
 *
 * \pre group != nullptr
 * \post topo != nullptr
 * \post blueprint::isValidTopologyGroup( topo ) == true
 *
 */
const sidre::Group* getTopologyGroup(const sidre::Group* group, const std::string& topo = "");

/*!
 * \brief Initialize the topology group.
 *
 * \param [in] group the blueprint mesh root group.
 * \param [in] topo the name of the topology to initialize.
 * \param [in] coordset the coordset to associate with the given topology.
 * \param [in] type the type of the given topology.
 *
 * \pre group != nullptr
 */
void initializeTopologyGroup(sidre::Group* group,
                             const std::string& topo,
                             const std::string& coordset,
                             const std::string& type);

/*!
 * \brief Returns the coordset group associated with the given topology group.
 *
 * \param [in] group the mesh blueprint root group.
 * \param [in] topology the topology group associated with the requested
 *  coordset.
 *
 * \return coordset pointer to the coordset associated with given topology,
 *  else false.
 *
 * \pre blueprint::isValidRootGroup( group )
 * \pre blueprint::isValidTopologyGroup( topology )
 * \post blueprint::isValidCoordsetGroup( coordset )
 */
const sidre::Group* getCoordsetGroup(const sidre::Group* group, const sidre::Group* topology);

/*!
 * \brief Returns the coordset group associated with the given topology group.
 *
 * \param [in] group the mesh blueprint root group.
 * \param [in] coords option argument to specify the coordset name.
 *
 * \return coordset pointer to the coordset associated with given topology,
 *  else false.
 *
 * \pre blueprint::isValidRootGroup( group )
 * \pre blueprint::isValidTopologyGroup( topology )
 * \post blueprint::isValidCoordsetGroup( coordset )
 */
const sidre::Group* getCoordsetGroup(const sidre::Group* group, const std::string& coords = "");

/*!
 * \brief Returns the mesh type and dimension given a root group that conforms
 *  to the computational mesh blueprint.
 *
 * \param [out] mesh_type the mint mesh type
 * \param [out] dimension the mesh dimension
 * \param [in] group pointer to the root group
 * \param [in] topo optional argument corresponding to the mesh topology.
 *
 * \note If a topology is not specified, the code assumes that the first group
 *  under the 'topologies' group corresponds to the mesh topology.
 *
 * \pre group != nullptr
 *
 * \see MeshTypes
 */
void getMeshTypeAndDimension(int& mesh_type,
                             int& dimension,
                             const sidre::Group* group,
                             const std::string& topology = "");

/*
 * \brief Return the whether the mesh has mixed cell types.
 *
 * \param [in] group pointer to the root group.
 * \param [in] topo optional argument corresponding to the mesh topology.
 *
 * \note If a topology is not specified, the code assumes that the first group
 *  under the 'topologies' group corresponds to the mesh topology.
 */
bool hasMixedCellTypes(const sidre::Group* group, const std::string& topo = "");

/*!
 * \brief Get the nodal dimensions and the global node extent of a
 *  StructuredMesh from sidre.
 *
 * \param [in] dim the expected dimension of the mesh.
 * \param [out] node_dims buffer to store the nodal dimensions.
 * \param [out] node_ext buffer to store the nodal extent.
 * \param [in] coordset pointer to the associated coordset group.
 *
 * \pre 1 <= dim <= 3
 * \pre node_dims != nullptr
 * \pre node_ext != nullptr
 * \pre coordset != nullptr
 * \pre isValidCoordsetGroup( coordset )
 *
 * \see setStructuredMeshProperties()
 */
void getStructuredMeshProperties(int dimension,
                                 IndexType node_dims[3],
                                 int64 node_ext[6],
                                 const sidre::Group* coordset);

/*!
 * \brief Put the nodal dimensions and the global node extent of a
 *  StructuredMesh into sidre.
 *
 * \param [in] dim the dimension of the mesh.
 * \param [in] node_dims buffer holding the nodal dimensions.
 * \param [in] node_ext buffer holding the nodal extent.
 * \param [out] coordset pointer to the associated coordset group.
 *
 * \pre 1 <= dim <= 3
 * \pre node_dims != nullptr
 * \pre node_ext != nullptr
 * \pre coordset != nullptr
 *
 * \see getStructuredMeshProperties()
 */
void setStructuredMeshProperties(int dimension,
                                 const IndexType node_dims[3],
                                 const int64 node_ext[6],
                                 sidre::Group* coordset);

/*!
 * \brief Put the global node extent of a StructuredMesh into sidre.
 *
 * \param [out] coordset pointer to the associated coordset group.
 * \param [in] node_ext buffer holding the nodal extent.
 *
 * \pre coordset != nullptr
 * \pre node_ext != nullptr
 *
 * \see setStructuredMeshProperties()
 */
void setExtent(sidre::Group* coordset, const int64 node_ext[6]);

/*!
 * \brief Returns the origin, spacing and extent of a uniform mesh from the
 *  corresponding coordset and topology groups associated with the mesh.
 *
 * \param [in] dim the expected dimension of the mesh.
 * \param [out] origin  buffer to store the origin of the uniform mesh.
 * \param [out] spacing buffer to store the spacing of the uniform mesh.
 * \param [in] coordset pointer to the associated coordset group.
 *
 * \note The `origin` and `spacing` pointers must point to buffers that have
 *  at least `dim` entries.
 *
 * \pre 1 <= dim <= 3
 * \pre origin != nullptr
 * \pre spacing != nullptr
 * \pre coordset != nullptr
 * \pre isValidCoordsetGroup( coordset )
 *
 * \see setUniformMeshProperties()
 */
void getUniformMeshProperties(int dim, double* origin, double* spacing, const sidre::Group* coordset);

/*!
 * \brief Populates the specified Coordset & Topology groups with the metadata
 *  for a uniform mesh
 *
 * \param [in] dim the mesh dimension.
 * \param [in] origin pointer to array with the origin's coordinates.
 * \param [in] spacing pointer to array with the spacing along each dimension.
 * \param [out] coordset pointer to the coordset group.
 *
 * \note The `origin` and `spacing` pointers must point to buffers that have
 *  at least `dim` entries.
 *
 * \pre 1 <= dim <= 3
 * \pre origin != nullptr
 * \pre spacing != nullptr
 * \pre coordset != nullptr
 *
 * \post isValidCoordsetGroup( coordset )
 *
 * \see getUniformMeshProperties()
 */
void setUniformMeshProperties(int dim,
                              const double* origin,
                              const double* spacing,
                              sidre::Group* coordset);

#endif /* AXOM_MINT_USE_SIDRE */

} /* namespace blueprint */

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_MESH_BLUEPRINT_HPP_ */
