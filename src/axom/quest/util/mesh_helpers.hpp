// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_MESH_HELPERS__HPP_
#define AXOM_QUEST_MESH_HELPERS__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#if defined(AXOM_USE_SIDRE)
  #include "axom/sidre.hpp"
#endif
#if defined(AXOM_USE_CONDUIT)
  #include <conduit/conduit_node.hpp>
#endif

#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

/*!
 * \brief Some helper functions to create and process mesh data
 */

namespace axom
{
namespace quest
{
namespace util
{
#ifdef AXOM_USE_MFEM
/*!
 * /brief Creates a 2D Cartesian mfem mesh lying within a given bounding box
 *
 * \param bbox The bounding box for the mesh
 * \param res The resolution of the mesh
 * \param polynomial_order The polynomial order for the mesh curvature
 * \param reorder_space_filling Reorder the mesh elements according to a space filling curve
 *
 * \return A pointer to the generated mesh
 * \note The user will be responsible for the mesh's memory
 */
mfem::Mesh* make_cartesian_mfem_mesh_2D(const primal::BoundingBox<double, 2>& bbox,
                                        const NumericArray<int, 2>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling = true);

/*!
 * /brief Creates a 3D Cartesian mfem mesh lying within a given bounding box
 *
 * \param bbox The bounding box for the mesh
 * \param res The resolution of the mesh
 * \param polynomial_order The polynomial order for the mesh curvature
 * \param reorder_space_filling Reorder the mesh elements according to a space filling curve
 *
 * \return A pointer to the generated mesh
 * \note The user will be responsible for the mesh's memory
 */
mfem::Mesh* make_cartesian_mfem_mesh_3D(const primal::BoundingBox<double, 3>& bbox,
                                        const NumericArray<int, 3>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling = true);
#endif

#if defined(AXOM_USE_SIDRE)
/*!
 * /brief Creates a 3D Cartesian blueprint mesh lying within a given bounding box
 *
 * \param meshGrp Put the mesh in this Group
 * \param bbox The bounding box for the mesh
 * \param res The resolution of the mesh
 * \param topologyName Name of the blueprint topoloyy to use.
 * \param coordsetName Name of the blueprint coordset to use.
 * \param runtimePolicy Runtime policy, see axom::runtime_policy.
 *        Memory in \c meshGrp must be compatible with the
 *        specified policy.
 *
 * All data that grows with mesh size will be allocated by the
 * default allocator id.  All other data are allocated with
 * axom::execution_space<axom::SEQ_EXEC>::allocatorID().
 *
 * \return The meshGrp pointer
 */
axom::sidre::Group* make_structured_blueprint_box_mesh_3d(
  axom::sidre::Group* meshGrp,
  const primal::BoundingBox<double, 3>& bbox,
  const NumericArray<int, 3>& res,
  const std::string& topologyName = "mesh",
  const std::string& coordsetName = "coords",
  axom::runtime_policy::Policy runtimePolicy = axom::runtime_policy::Policy::seq);

/*!
 * \brief Make an unstructured blueprint box mesh.
 *
 * All data that grows with mesh size will be allocated by the
 * default allocator id.  All other data are allocated with
 * axom::execution_space<axom::SEQ_EXEC>::allocatorID().
 */
axom::sidre::Group* make_unstructured_blueprint_box_mesh_3d(
  axom::sidre::Group* meshGrp,
  const primal::BoundingBox<double, 3>& bbox,
  const NumericArray<int, 3>& res,
  const std::string& topologyName = "mesh",
  const std::string& coordsetName = "coords",
  axom::runtime_policy::Policy runtimePolicy = axom::runtime_policy::Policy::seq);

/*!
 * /brief Creates a 2D Cartesian blueprint mesh lying within a given bounding box
 *
 * \param meshGrp Put the mesh in this Group
 * \param bbox The bounding box for the mesh
 * \param res The resolution of the mesh
 * \param topologyName Name of the blueprint topoloyy to use.
 * \param coordsetName Name of the blueprint coordset to use.
 * \param runtimePolicy Runtime policy, see axom::runtime_policy.
 *        Memory in \c meshGrp must be compatible with the
 *        specified policy.
 *
 * All data that grows with mesh size will be allocated by the
 * default allocator id.  All other data are allocated with
 * axom::execution_space<axom::SEQ_EXEC>::allocatorID().
 *
 * \return The meshGrp pointer
 */
axom::sidre::Group* make_structured_blueprint_box_mesh_2d(
  axom::sidre::Group* meshGrp,
  const primal::BoundingBox<double, 2>& bbox,
  const NumericArray<int, 2>& res,
  const std::string& topologyName = "mesh",
  const std::string& coordsetName = "coords",
  axom::runtime_policy::Policy runtimePolicy = axom::runtime_policy::Policy::seq);

/*!
 * All data that grows with mesh size will be allocated by the
 * default allocator id.  All other data are allocated with
 * axom::execution_space<axom::SEQ_EXEC>::allocatorID().
 */
axom::sidre::Group* make_unstructured_blueprint_box_mesh_2d(
  axom::sidre::Group* meshGrp,
  const primal::BoundingBox<double, 2>& bbox,
  const NumericArray<int, 2>& res,
  const std::string& topologyName = "mesh",
  const std::string& coordsetName = "coords",
  axom::runtime_policy::Policy runtimePolicy = axom::runtime_policy::Policy::seq);

/*!
 * \brief Convert a structured explicit blueprint mesh to unstructured.
 * \param meshGrp Put the mesh in this Group
 * \param topologyName Name of the blueprint topoloyy to use.
 * \param runtimePolicy Runtime policy, see axom::runtime_policy.
 *        Memory in \c meshGrp must be compatible with the
 *        specified policy.
 *
 * All input mesh data are expected to have the allocator id of
 * meshGrp->getDefaultAllocatorID().  On output, they will also have
 * the same allocator id, even if the intermediate steps have to
 * transfers memory to and from another space.
 */
void convert_blueprint_structured_explicit_to_unstructured_3d(axom::sidre::Group* meshGrp,
                                                              const std::string& topoName,
                                                              axom::runtime_policy::Policy runtimePolicy);

template <typename ExecSpace>
void convert_blueprint_structured_explicit_to_unstructured_impl_3d(axom::sidre::Group* meshGrp,
                                                                   const std::string& topoName);

void convert_blueprint_structured_explicit_to_unstructured_2d(axom::sidre::Group* meshGrp,
                                                              const std::string& topoName,
                                                              axom::runtime_policy::Policy runtimePolicy);

template <typename ExecSpace>
void convert_blueprint_structured_explicit_to_unstructured_impl_2d(axom::sidre::Group* meshGrp,
                                                                   const std::string& topoName);

  #if defined(AXOM_USE_CONDUIT)
/*!
 * \brief Check if blueprint mesh is valid.
 */
bool verifyBlueprintMesh(const axom::sidre::Group* meshGrp, conduit::Node info);
  #endif

#endif

/*!
 * @brief Fill in 3D structured mesh cartesian coordinates.
 * \param runtimePolicy Runtime execution space selector
 * \param xView Vertex x-values array
 * \param yView Vertex y-values array
 * \param zView Vertex z-values array
 * \param domainBox Physical domain
 *
 * Array data must be in a memory space accessible by the
 * selected \c runtimePolicy.
 */
void fill_cartesian_coords_3d(axom::runtime_policy::Policy runtimePolicy,
                              const primal::BoundingBox<double, 3>& domainBox,
                              axom::ArrayView<double, 3>& xView,
                              axom::ArrayView<double, 3>& yView,
                              axom::ArrayView<double, 3>& zView);

/*!
 * @brief Fill in 2D structured mesh cartesian coordinates.
 * \param runtimePolicy Runtime execution space selector
 * \param xView Vertex x-values array
 * \param yView Vertex y-values array
 * \param domainBox Physical domain
 *
 * Array data must be in a memory space accessible by the
 * selected \c runtimePolicy.
 */
void fill_cartesian_coords_2d(axom::runtime_policy::Policy runtimePolicy,
                              const primal::BoundingBox<double, 2>& domainBox,
                              axom::ArrayView<double, 2>& xView,
                              axom::ArrayView<double, 2>& yView);

/*!
 * @brief Fill in 3D structured mesh cartesian coordinates.
 * \tparam ExecSpace Execution space, e.g. \c axom::SEQ_EXEC.
 *         @see axom::execution_space
 * \param xView Vertex x-values array
 * \param yView Vertex y-values array
 * \param zView Vertex z-values array
 * \param domainBox Physical domain
 */
template <typename ExecSpace>
void fill_cartesian_coords_3d_impl(const primal::BoundingBox<double, 3>& domainBox,
                                   axom::ArrayView<double, 3>& xView,
                                   axom::ArrayView<double, 3>& yView,
                                   axom::ArrayView<double, 3>& zView);

/*!
 * @brief Fill in 2D structured mesh cartesian coordinates.
 * \tparam ExecSpace Execution space, e.g. \c axom::SEQ_EXEC.
 *         @see axom::execution_space
 * \param xView Vertex x-values array
 * \param yView Vertex y-values array
 * \param domainBox Physical domain
 */
template <typename ExecSpace>
void fill_cartesian_coords_2d_impl(const primal::BoundingBox<double, 2>& domainBox,
                                   axom::ArrayView<double, 2>& xView,
                                   axom::ArrayView<double, 2>& yView);

}  // namespace util
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_MESH_HELPERS__HPP_
