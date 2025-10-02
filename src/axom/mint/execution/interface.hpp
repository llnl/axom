// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_EXECUTION_INTERFACE_HPP_
#define MINT_EXECUTION_INTERFACE_HPP_

#include "axom/config.hpp"                          // compile-time definitions
#include "axom/core/Macros.hpp"                     // for AXOM_STATIC_ASSERT
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits

#include "axom/mint/execution/xargs.hpp"                   // for xargs
#include "axom/mint/execution/internal/for_all_cells.hpp"  // for_all_cells()
#include "axom/mint/execution/internal/for_all_nodes.hpp"  // for_all_nodes()
#include "axom/mint/execution/internal/for_all_faces.hpp"  // for_all_faces()

#include "axom/mint/mesh/Mesh.hpp"  // for  mint::Mesh

// Slic includes
#include "axom/slic/interface/slic.hpp"  // for SLIC macros

// C/C++ includes
#include <type_traits>  // for std::is_same

/*!
 * \file
 *
 * \brief Defines an execution interface for generic loop traversals and various
 *  mesh traversals, e.g., loop over nodes, cells, etc. The traversals may be
 *  executed in parallel, on the CPU, GPU, or other accelerator device, or
 *  in serial on the CPU according to the execution policy, which is supplied as
 *  a template argument to each of the traversal functions.
 *
 *  The general interface for the mesh traversal functions takes in a mint Mesh
 *  object and the loop body encapsulated in a lambda expression, conforming to
 *  the following template:
 *  \code
 *     for_all_[entity]< [exec_policy], [xargs] >( meshPtr, AXOM_LAMBDA(...) {
 *          // loop body
 *     } );
 *  \endcode
 *
 *  where:
 *
 *   * <b> [entity] </b> <br/>
 *     The entity suffix specifies the mesh entity being traversed, e.g., nodes,
 *     cells, faces, etc,.
 *
 *   * <b> [exec_policy] </b> <br />
 *     The execution policy indicates <em> how </em> and <em> where </em>
 *     the data corresponding to the specified mesh entity will be traversed.
 *     For example, an execution policy may indicate that the mesh traversal
 *     will be executed in parallel on the GPU, or CPU. A list of the
 *     currently supported execution policies and a brief description is given
 *     in execution_space.hpp
 *
 *   * <b> [AXOM_LAMBDA()] </b> <br />
 *     The AXOM_LAMBDA argument encapsulates the loop body, i.e., the kernel
 *     that is to be executed at each mesh entity. By default, the lambda
 *     expression takes the ID of the corresponding mesh entity as an argument.
 *     Additional lambda arguments may be specified by supplying the
 *     <em> xargs </em> template argument to the mesh traversal function.
 *
 *   * <b> [xargs] </b> <br />
 *     The xargs is an optional template argument to the mesh traversal function
 *     that specifies any additional arguments that the supplied lambda
 *     expression takes.A list of the available xargs types and a brief
 *     description is given in xargs.hpp
 *
 *  \note For parallel execution, the implementation relies heavily on the
 *   <a href="https://github.com/LLNL/RAJA"> RAJA </a> programming model
 *   abstraction layer.
 *
 *  \see execution_space.hpp
 *  \see xargs.hpp
 */

namespace axom
{
namespace mint
{
/// \name Mesh Node Traversal Functions
/// @{

/*!
 * \brief Loops over the nodes of the given mesh.
 *
 * \param [in] m pointer to the mesh object
 * \param [in] kernel user-supplied lambda consisting of the kernel
 *
 * \pre m != nullptr
 *
 * \tparam ExecPolicy the execution policy, e.g., serial or parallel
 * \tparam ArgType object indicating the arguments to the kernel
 * \tparam KernelType
 *
 * Usage Example:
 * \code
 *
 *  for_all_nodes< exec >( m,
 *    AXOM_LAMBDA( IndexType nodeID )
 *    { ... }
 *  );
 *
 *  for_all_nodes< exec, xargs::ij >( m,
 *    AXOM_LAMBDA( IndexType nodeID, IndexType i, IndexType j )
 *    { ... }
 *  );
 *
 *  for_all_nodes< exec, xargs::ijk >( m,
 *    AXOM_LAMDA( IndexType nodeIdx, IndexType i, IndexType j, IndexType k)
 *    { ... }
 *  );
 *
 *  for_all_nodes< exec, xargs::xy >( m,
 *    AXOM_LAMBA( IndexType nodeIdx, double x, double y)
 *    { ... }
 *  );
 *
 *  for_all_nodes< exec, xargs::xyz >( m,
 *    AXOM_LAMBA( IndexType nodeIdx, double x, double y, double z)
 *    { ... }
 *  );
 *
 * \endcode
 *
 * \see execution_space.hpp
 * \see xargs.hpp
 */
/// @{
template <typename ExecPolicy, typename ArgType = xargs::index, typename MeshType, typename KernelType>
inline void for_all_nodes(const MeshType* m, KernelType&& kernel)
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT(execution_space<ExecPolicy>::valid());
  AXOM_STATIC_ASSERT(xargs_traits<ArgType>::valid());

  constexpr bool valid_mesh_type = std::is_base_of<Mesh, MeshType>::value;
  AXOM_STATIC_ASSERT(valid_mesh_type);

  // run-time sanity checks
  SLIC_ASSERT(m != nullptr);

  // dispatch
  internal::for_all_nodes_impl<ExecPolicy>(ArgType(), *m, std::forward<KernelType>(kernel));
}

template <typename ExecPolicy, typename ArgType = xargs::index, typename KernelType>
inline void for_all_nodes(const Mesh* m, KernelType&& kernel)
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT(execution_space<ExecPolicy>::valid());
  AXOM_STATIC_ASSERT(xargs_traits<ArgType>::valid());

  // run-time sanity checks
  SLIC_ASSERT(m != nullptr);

  //dispatch
  internal::for_all_nodes<ExecPolicy>(ArgType(), *m, std::forward<KernelType>(kernel));
}

/// @}

/// @}

/// \name Mesh Cell Traversal Functions
/// @{

/*!
 * \brief Loops over all the cells of a given mesh.
 *
 * \param [in] m pointer to the mesh object.
 * \param [in] kernel user-supplied kernel to execute on each cell.
 *
 * \pre m != nullptr
 *
 * \tparam ExecPolicy the execution policy, e.g., serial or parallel
 * \tparam ArgType object indicating the arguments to the kernel
 *
 * Usage Example:
 * \code
 *
 *  for_all_cells< exec >( m,
 *    AXOM_LAMBDA( IndexType cellID )
 *    { ... }
 *  );
 *
 *  for_all_cells< exec, xargs::ij >( m,
 *    AXOM_LAMBDA( IndexType cellID, IndexType i, IndexType j )
 *    { ... }
 *  );
 *
 *  for_all_cells< exec, xargs::ijk >( m,
 *    AXOM_LAMBDA( IndexType cellID, IndexType i, IndexType j, IndexType k )
 *    { ... }
 *  );
 *
 *  for_all_cells< exec, xargs::nodeids >( m,
 *    AXOM_LAMBDA( IndexType cellID, const IndexType* nodeIDs, IndexType N )
 *    { ... }
 *  );
 *
 *  for_all_cells< exec, xargs::coords >( m,
 *    AXOM_LAMBDA( IndexType cellID, numerics::Matrix<double>& coords,
 *                 const IndexType * nodeIDs )
 *    { ... }
 *  );
 *
 *  for_all_cells< exec, xargs::faceids >( m,
 *    AXOM_LAMBDA( IndexType cellID, const IndexType* faceIDs, IndexType N )
 *    { ... }
 *  );
 *
 * \endcode
 *
 * \see execution_space.hpp
 * \see xargs.hpp
 */
/// @{

template <typename ExecPolicy, typename ArgType = xargs::index, typename MeshType, typename KernelType>
inline void for_all_cells(const MeshType* m, KernelType&& kernel)
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT(execution_space<ExecPolicy>::valid());
  AXOM_STATIC_ASSERT(xargs_traits<ArgType>::valid());

  constexpr bool valid_mesh_type = std::is_base_of<Mesh, MeshType>::value;
  AXOM_STATIC_ASSERT(valid_mesh_type);

  // run-time sanity checks
  SLIC_ASSERT(m != nullptr);

  // dispatch
  internal::for_all_cells_impl<ExecPolicy>(ArgType(), *m, std::forward<KernelType>(kernel));
}

template <typename ExecPolicy, typename ArgType = xargs::index, typename KernelType>
inline void for_all_cells(const Mesh* m, KernelType&& kernel)
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT(execution_space<ExecPolicy>::valid());
  AXOM_STATIC_ASSERT(xargs_traits<ArgType>::valid());

  // run-time sanity checks
  SLIC_ASSERT(m != nullptr);

  //dispatch
  internal::for_all_cells<ExecPolicy>(ArgType(), *m, std::forward<KernelType>(kernel));
}

/// @}

/// @}

/// \name Mesh Face Traversal Functions
/// @{

/*!
 * \brief Loops over all the faces of a given mesh.
 *
 * \param [in] m pointer to the mesh object.
 * \param [in] kernel user-supplied kernel to execute on each face.
 *
 * \pre m != nullptr
 *
 * \tparam ExecPolicy the execution policy, e.g., serial or parallel
 * \tparam ArgType object indicating the arguments to the kernel
 *
 * Usage Example:
 * \code
 *
 *   for_all_faces< exec >( m,
 *    AXOM_LAMBDA( IndexType faceID )
 *    { ... }
 *  );
 *
 *  for_all_faces< exec, xargs::nodeids >( m,
 *    AXOM_LAMBDA( IndexType faceID, const IndexType* nodeIDs, IndexType N )
 *    { ... }
 *  );
 *
 *  for_all_faces< exec, xargs::coords >( m,
 *    AXOM_LAMBDA( IndexType faceID, numerics::Matrix<double>& coords,
 *                 const IndexType * nodeIDs )
 *    { ... }
 *  );
 *
 *  for_all_faces< exec, xargs::cellids >( m,
 *    AXOM_LAMBDA( IndexType faceID, IndexType cellIDOne, IndexType cellIDTwo )
 *    { ... }
 *  );
 *
 * \endcode
 *
 * \note A face can be associated with one or two cells, depending on whether
 *  it is an external boundary face or interior face. By convention, if a face
 *  is an external boundary face, then only cellIDOne exists and cellIDTwo
 *  will be set to -1.
 *
 * \see execution_space.hpp
 * \see xargs.hpp
 */
/// @{

template <typename ExecPolicy, typename ArgType = xargs::index, typename MeshType, typename KernelType>
inline void for_all_faces(const MeshType* m, KernelType&& kernel)
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT(execution_space<ExecPolicy>::valid());
  AXOM_STATIC_ASSERT(xargs_traits<ArgType>::valid());

  constexpr bool valid_mesh_type = std::is_base_of<Mesh, MeshType>::value;
  AXOM_STATIC_ASSERT(valid_mesh_type);

  // run-time sanity checks
  SLIC_ASSERT(m != nullptr);

  // dispatch
  internal::for_all_faces_impl<ExecPolicy>(ArgType(), *m, std::forward<KernelType>(kernel));
}

template <typename ExecPolicy, typename ArgType = xargs::index, typename KernelType>
inline void for_all_faces(const Mesh* m, KernelType&& kernel)
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT(execution_space<ExecPolicy>::valid());
  AXOM_STATIC_ASSERT(xargs_traits<ArgType>::valid());

  // run-time sanity checks
  SLIC_ASSERT(m != nullptr);

  //dispatch
  internal::for_all_faces<ExecPolicy>(ArgType(), *m, std::forward<KernelType>(kernel));
}

/// @}

/// @}

}  // namespace mint
}  // namespace axom
#endif /* MINT_EXECUTION_INTERFACE_HPP_ */
