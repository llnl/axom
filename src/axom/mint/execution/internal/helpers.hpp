// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_EXECUTION_HELPERS_HPP_
#define MINT_EXECUTION_HELPERS_HPP_

// mint includes
#include "axom/mint/config.hpp"     // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"  // for Mesh

#include "axom/core/Macros.hpp"
#include "axom/core/StackArray.hpp"       // for axom::StackArray
#include "axom/core/numerics/Matrix.hpp"  // for Matrix

namespace axom
{
namespace mint
{
namespace internal
{
/*!
 * \brief Iterate over the objects (cells or faces) in a mesh and for each
 *  object construct a NDIM x NNODES matrix of the nodal coordinates of the
 *  object and call a kernel with the object ID, the coordinate matrix, and the
 *  node IDs. 
 *
 * \tparam ExecPolicy the execution policy
 * \tparam NDIM the number of coordinate dimensions.
 * \tparam NNODES the number of nodes per object.
 * \tparam FOR_ALL_FUNCTOR the type of the for_all_nodes functor.
 * \tparam KernelType the type of the supplied lambda expression.
 *
 * \param [in] for_all_nodes functor that iterates over the objects in the mesh
 *  and for each object calls a function with the object ID, the node IDs and
 *  the number of nodes.
 * \param [in] m the Mesh to iterate over.
 * \param [in] kernel the kernel to call on each object.
 *
 */

template <typename ExecPolicy, int NDIM, int NNODES, typename FOR_ALL_FUNCTOR, typename MeshType, typename KernelType>
inline void for_all_coords(const FOR_ALL_FUNCTOR& for_all_nodes, const MeshType& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_UNIFORM_MESH, "Not valid for UniformMesh.");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_RECTILINEAR_MESH, "Not valid for RectilinearMesh.");

  AXOM_STATIC_ASSERT_MSG(NDIM >= 1 && NDIM <= 3, "NDIM must be a valid dimension.");
  AXOM_STATIC_ASSERT_MSG(NNODES > 0, "NNODES must be greater than zero.");

  constexpr bool valid_mesh_type = std::is_base_of<Mesh, MeshType>::value;
  AXOM_STATIC_ASSERT(valid_mesh_type);

  SLIC_ERROR_IF(m.getDimension() != NDIM, "Dimension mismatch!");

  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  constexpr bool NO_COPY = true;

  IndexType coordinate_size = m.getNumberOfNodes();

  // Extract coordinate values into an axom::ArrayView
  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE), coordinate_size);
  auto y_vals_h = (NDIM > 1)
    ? axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE), coordinate_size)
    : axom::ArrayView<const double>();
  auto z_vals_h = (NDIM > 2)
    ? axom::ArrayView<const double>(m.getCoordinateArray(Z_COORDINATE), coordinate_size)
    : axom::ArrayView<const double>();

  // Move xyz values onto device
  axom::Array<double> x_vals_d = axom::Array<double>(x_vals_h, device_allocator);
  auto x_vals_view = x_vals_d.view();

  axom::Array<double> y_vals_d =
    (NDIM > 1) ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();
  auto y_vals_view = (NDIM > 1) ? y_vals_d.view() : y_vals_h;

  axom::Array<double> z_vals_d =
    (NDIM > 2) ? axom::Array<double>(z_vals_h, device_allocator) : axom::Array<double>();
  auto z_vals_view = (NDIM > 2) ? z_vals_d.view() : z_vals_h;

  for_all_nodes(
    ExecPolicy(),
    m,
    AXOM_LAMBDA(IndexType objectID, const IndexType* nodeIDs, IndexType numNodes) {
      AXOM_UNUSED_VAR(numNodes);
      assert(numNodes == NNODES);

      double localCoords[NDIM * NNODES];
      for(int i = 0; i < NNODES; ++i)
      {
        const int i_offset = NDIM * i;

        localCoords[i_offset] = x_vals_view[nodeIDs[i]];
        if(NDIM > 1)
        {
          localCoords[i_offset + 1] = y_vals_view[nodeIDs[i]];
        }
        if(NDIM > 2)
        {
          localCoords[i_offset + 2] = z_vals_view[nodeIDs[i]];
        }
      }

      numerics::Matrix<double> coordsMatrix(NDIM, NNODES, localCoords, NO_COPY);
      kernel(objectID, coordsMatrix, nodeIDs);
    });
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_EXECUTION_HELPERS_HPP_ */