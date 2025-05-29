// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FOR_ALL_NODES_HPP_
#define MINT_FOR_ALL_NODES_HPP_

// Axom core includes
#include "axom/config.hpp"                          // compile time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits
#include "axom/core/execution/for_all.hpp"          // for axom::for_all

// mint includes
#include "axom/mint/execution/xargs.hpp"       // for xargs
#include "axom/mint/config.hpp"                // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"             // for mint::Mesh
#include "axom/mint/mesh/RectilinearMesh.hpp"  // for mint::RectilinearMesh
#include "axom/mint/mesh/StructuredMesh.hpp"   // for mint::StructuredMesh
#include "axom/mint/mesh/UniformMesh.hpp"      // for mint::UniformMesh
#include "axom/core/execution/nested_for_exec.hpp"

#include "axom/core/StackArray.hpp"  // for axom::StackArray

namespace axom
{
namespace mint
{
namespace internal
{
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::index, const Mesh& m, KernelType&& kernel)
{
  const IndexType numNodes = m.getNumberOfNodes();
  axom::for_all<ExecPolicy>(numNodes, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::index, const Mesh& m, KernelType&& kernel)
{
  for_all_nodes_impl<ExecPolicy>(xargs::index(), m, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::ij, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::ij is only valid for 2D structured meshes!");

  const IndexType jp = m.nodeJp();
  const IndexType Ni = m.getNodeResolution(I_DIRECTION);
  const IndexType Nj = m.getNodeResolution(J_DIRECTION);

  axom::StackArray<IndexType, 2> i_range {{0, Ni}}, j_range {{0, Nj}};
  axom::for_all<ExecPolicy>(
    i_range,
    j_range,
    AXOM_LAMBDA(IndexType i, IndexType j) {
      const IndexType nodeIdx = i + j * jp;
      kernel(nodeIdx, i, j);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::ij, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(!m.isStructured(), "xargs::ijk is only valid for 2D structured meshes!");
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::ijk is only valid for 2D structured meshes!");

  const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
  for_all_nodes_impl<ExecPolicy>(xargs::ij(), sm, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::ijk, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::ijk is only valid for 3D structured meshes!");

  const IndexType jp = m.nodeJp();
  const IndexType kp = m.nodeKp();
  const IndexType Ni = m.getNodeResolution(I_DIRECTION);
  const IndexType Nj = m.getNodeResolution(J_DIRECTION);
  const IndexType Nk = m.getNodeResolution(K_DIRECTION);

  axom::StackArray<IndexType, 2> i_range {{0, Ni}}, j_range {{0, Nj}}, k_range {{0, Nk}};
  axom::for_all<ExecPolicy>(
    i_range,
    j_range,
    k_range,
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType nodeIdx = i + j * jp + k * kp;
      kernel(nodeIdx, i, j, k);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::ijk, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(!m.isStructured(), "xargs::ijk is only valid for 3D structured meshes!");
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::ijk is only valid for 3D structured meshes!");

  const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
  for_all_nodes_impl<ExecPolicy>(xargs::ijk(), sm, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::x, const UniformMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 1, "xargs::x is only valid for 1D meshes");

  const double x0 = m.getOrigin()[0];
  const double dx = m.getSpacing()[0];

  for_all_nodes_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType nodeID) {
      const double x = x0 + nodeID * dx;
      kernel(nodeID, x);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::x, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 1, "xargs::x is only valid for 1D meshes");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_UNIFORM_MESH, "Not valid for UniformMesh.");
  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();
  IndexType coordinate_size = m.getNumberOfNodes();

  // extract coordinate values into an axom::Array
  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE), coordinate_size);

  SLIC_ASSERT(x_vals_h.data() != nullptr);

  // Move x values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, device_allocator) : axom::Array<double>();
  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;

  for_all_nodes_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType nodeID) { kernel(nodeID, x_vals_view[nodeID]); });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::x, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 1, "xargs::x is only valid for 1D meshes");

  const int mesh_type = m.getMeshType();
  if(mesh_type == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::x(), um, std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_nodes_impl<ExecPolicy>(xargs::x(), m, std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xy, const UniformMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::xy is only valid for 2D meshes");
  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  // extract origin and spacing into an axom::Array
  auto origin_h = axom::ArrayView<const double>(m.getOrigin().begin(), 3);
  auto spacing_h = axom::ArrayView<const double>(m.getSpacing().begin(), 3);

  // Move origin and spacing values onto device
  axom::Array<double> origin_d =
    on_device ? axom::Array<double>(origin_h, device_allocator) : axom::Array<double>();
  auto origin_view = on_device ? origin_d.view() : origin_h;

  axom::Array<double> spacing_d =
    on_device ? axom::Array<double>(spacing_h, device_allocator) : axom::Array<double>();
  auto spacing_view = on_device ? spacing_d.view() : spacing_h;

  for_all_nodes_impl<ExecPolicy>(
    xargs::ij(),
    m,
    AXOM_LAMBDA(IndexType nodeID, IndexType i, IndexType j) {
      const double x = origin_view[0] + i * spacing_view[0];
      const double y = origin_view[1] + j * spacing_view[1];
      kernel(nodeID, x, y);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xy, const RectilinearMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::xy is only valid for 2D meshes");

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  // extract coordinate values into an axom::Array
  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE),
                                                m.getNodeResolution(X_COORDINATE));
  auto y_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE),
                                                m.getNodeResolution(Y_COORDINATE));

  SLIC_ASSERT(x_vals_h.data() != nullptr);
  SLIC_ASSERT(y_vals_h.data() != nullptr);

  // Move xy values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, device_allocator) : axom::Array<double>();
  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;

  axom::Array<double> y_vals_d =
    on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();
  auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

  for_all_nodes_impl<ExecPolicy>(
    xargs::ij(),
    m,
    AXOM_LAMBDA(IndexType nodeID, IndexType i, IndexType j) {
      kernel(nodeID, x_vals_view[i], y_vals_view[j]);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xy, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::xy is only valid for 2D meshes");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_UNIFORM_MESH, "Not valid for UniformMesh.");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_RECTILINEAR_MESH, "Not valid for RectilinearMesh.");

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();
  IndexType coordinate_size = m.getNumberOfNodes();

  // extract coordinate values into an axom::Array
  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE), coordinate_size);
  auto y_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE), coordinate_size);

  SLIC_ASSERT(x_vals_h.data() != nullptr);
  SLIC_ASSERT(y_vals_h.data() != nullptr);

  // Move xy values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, device_allocator) : axom::Array<double>();
  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;

  axom::Array<double> y_vals_d =
    on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();
  auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

  for_all_nodes_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType nodeID) { kernel(nodeID, x_vals_view[nodeID], y_vals_view[nodeID]); });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::xy, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::xy is only valid for 3D meshes");

  const int mesh_type = m.getMeshType();
  if(mesh_type == STRUCTURED_RECTILINEAR_MESH)
  {
    const RectilinearMesh& rm = static_cast<const RectilinearMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::xy(), rm, std::forward<KernelType>(kernel));
  }
  else if(mesh_type == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::xy(), um, std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_nodes_impl<ExecPolicy>(xargs::xy(), m, std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xyz, const UniformMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::xyz is only valid for 3D meshes");
  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  // extract origin and spacing into an axom::Array
  auto origin_h = axom::ArrayView<const double>(m.getOrigin().begin(), 3);
  auto spacing_h = axom::ArrayView<const double>(m.getSpacing().begin(), 3);

  // Move origin and spacing values onto device
  axom::Array<double> origin_d =
    on_device ? axom::Array<double>(origin_h, device_allocator) : axom::Array<double>();
  auto origin_view = on_device ? origin_d.view() : origin_h;

  axom::Array<double> spacing_d =
    on_device ? axom::Array<double>(spacing_h, device_allocator) : axom::Array<double>();
  auto spacing_view = on_device ? spacing_d.view() : spacing_h;

  for_all_nodes_impl<ExecPolicy>(
    xargs::ijk(),
    m,
    AXOM_LAMBDA(IndexType nodeID, IndexType i, IndexType j, IndexType k) {
      const double x = origin_view[0] + i * spacing_view[0];
      const double y = origin_view[1] + j * spacing_view[1];
      const double z = origin_view[2] + k * spacing_view[2];
      kernel(nodeID, x, y, z);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xyz, const RectilinearMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::xyz is only valid for 3D meshes");
  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  // extract coordinate values into an axom::Array
  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE),
                                                m.getNodeResolution(X_COORDINATE));
  auto y_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE),
                                                m.getNodeResolution(Y_COORDINATE));
  auto z_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Z_COORDINATE),
                                                m.getNodeResolution(Z_COORDINATE));

  // Move xyz values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, device_allocator) : axom::Array<double>();
  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;

  axom::Array<double> y_vals_d =
    on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();
  auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

  axom::Array<double> z_vals_d =
    on_device ? axom::Array<double>(z_vals_h, device_allocator) : axom::Array<double>();
  auto z_vals_view = on_device ? z_vals_d.view() : z_vals_h;

  for_all_nodes_impl<ExecPolicy>(
    xargs::ijk(),
    m,
    AXOM_LAMBDA(IndexType nodeID, IndexType i, IndexType j, IndexType k) {
      kernel(nodeID, x_vals_view[i], y_vals_view[j], z_vals_view[k]);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes_impl(xargs::xyz, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::xyz is only valid for 3D meshes");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_UNIFORM_MESH, "Not valid for UniformMesh.");
  SLIC_ERROR_IF(m.getMeshType() == STRUCTURED_RECTILINEAR_MESH, "Not valid for RectilinearMesh.");

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();
  IndexType coordinate_size = m.getNumberOfNodes();

  // extract coordinate values into an axom::Array
  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE), coordinate_size);
  auto y_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE), coordinate_size);
  auto z_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Z_COORDINATE), coordinate_size);

  SLIC_ASSERT(x_vals_h.data() != nullptr);
  SLIC_ASSERT(y_vals_h.data() != nullptr);
  SLIC_ASSERT(z_vals_h.data() != nullptr);

  // Move xyz values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, device_allocator) : axom::Array<double>();
  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;

  axom::Array<double> y_vals_d =
    on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();
  auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

  axom::Array<double> z_vals_d =
    on_device ? axom::Array<double>(z_vals_h, device_allocator) : axom::Array<double>();
  auto z_vals_view = on_device ? z_vals_d.view() : z_vals_h;
#if !defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA)
  // If we have CUDA but not RAJA then we are doing serial execution and cannot
  // mark this intermediate lambda as AXOM_LAMBDA to comply with CUDA's lambda
  // rules.
  for_all_nodes_impl<ExecPolicy>(xargs::index(), m, [&](IndexType nodeID) {
    kernel(nodeID, x_vals_view[nodeID], y_vals_view[nodeID], z_vals_view[nodeID]);
  });
#else
  for_all_nodes_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType nodeID) {
      kernel(nodeID, x_vals_view[nodeID], y_vals_view[nodeID], z_vals_view[nodeID]);
    });
#endif
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_nodes(xargs::xyz, const mint::Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::xyz is only valid for 3D meshes");

  const int mesh_type = m.getMeshType();
  if(mesh_type == STRUCTURED_RECTILINEAR_MESH)
  {
    const RectilinearMesh& rm = static_cast<const RectilinearMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::xyz(), rm, std::forward<KernelType>(kernel));
  }
  else if(mesh_type == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_nodes_impl<ExecPolicy>(xargs::xyz(), um, std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_nodes_impl<ExecPolicy>(xargs::xyz(), m, std::forward<KernelType>(kernel));
  }
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINTFOR_ALL_NODES_STRUCTURED_HPP_ */
