// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FOR_ALL_CELLS_HPP_
#define MINT_FOR_ALL_CELLS_HPP_

// Axom core includes
#include "axom/config.hpp"                          // compile time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits
#include "axom/core/execution/for_all.hpp"          // for axom::for_all

// mint includes
#include "axom/mint/execution/xargs.hpp"  // for xargs

#include "axom/mint/config.hpp"                 // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"              // for Mesh
#include "axom/mint/mesh/StructuredMesh.hpp"    // for StructuredMesh
#include "axom/mint/mesh/UniformMesh.hpp"       // for UniformMesh
#include "axom/mint/mesh/RectilinearMesh.hpp"   // for RectilinearMesh
#include "axom/mint/mesh/CurvilinearMesh.hpp"   // for CurvilinearMesh
#include "axom/mint/mesh/UnstructuredMesh.hpp"  // for UnstructuredMesh
#include "axom/mint/execution/internal/helpers.hpp"
#include "axom/core/execution/nested_for_exec.hpp"

#include "axom/core/StackArray.hpp"       // for axom::StackArray
#include "axom/core/numerics/Matrix.hpp"  // for Matrix

namespace axom
{
namespace mint
{
namespace internal
{
//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::index, const Mesh& m, KernelType&& kernel)
{
  const IndexType numCells = m.getNumberOfCells();
  axom::for_all<ExecPolicy>(numCells, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::index, const Mesh& m, KernelType&& kernel)
{
  for_all_cells_impl<ExecPolicy>(xargs::index(), m, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::ij, const StructuredMesh& m, KernelType&& kernel)
{
  // run-time checks
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::ij is only valid for 2D meshes!");

  const IndexType jp = m.cellJp();
  const IndexType Ni = m.getCellResolution(I_DIRECTION);
  const IndexType Nj = m.getCellResolution(J_DIRECTION);

  axom::StackArray<IndexType, 2> i_range {{0, Ni}}, j_range {{0, Nj}};
  axom::for_all<ExecPolicy>(
    i_range,
    j_range,
    AXOM_LAMBDA(IndexType i, IndexType j) {
      const IndexType cellID = i + j * jp;
      kernel(cellID, i, j);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::ij, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(!m.isStructured(), "xargs::ij is only valid on structured meshes!");
  SLIC_ERROR_IF(m.getDimension() != 2, "xargs::ij is only valid on 2D meshes!");

  const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
  for_all_cells_impl<ExecPolicy>(xargs::ij(), sm, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::ijk, const StructuredMesh& m, KernelType&& kernel)
{
  const IndexType Ni = m.getCellResolution(I_DIRECTION);
  const IndexType Nj = m.getCellResolution(J_DIRECTION);
  const IndexType Nk = m.getCellResolution(K_DIRECTION);

  const IndexType jp = m.cellJp();
  const IndexType kp = m.cellKp();

  axom::StackArray<IndexType, 2> i_range {{0, Ni}}, j_range {{0, Nj}}, k_range {{0, Nk}};
  axom::for_all<ExecPolicy>(
    i_range,
    j_range,
    k_range,
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType cellID = i + j * jp + k * kp;
      kernel(cellID, i, j, k);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::ijk, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(!m.isStructured(), "xargs::ijk is only valid on structured meshes!");
  SLIC_ERROR_IF(m.getDimension() != 3, "xargs::ijk is only valid for 3D structured meshes!");

  const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
  for_all_cells_impl<ExecPolicy>(xargs::ijk(), sm, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::nodeids, const StructuredMesh& m, KernelType&& kernel)
{
  const int dimension = m.getDimension();
  const IndexType nodeJp = m.nodeJp();
  const IndexType nodeKp = m.nodeKp();
  const StackArray<IndexType, 8>& offsets = m.getCellNodeOffsetsArray();

  // Note: gcc@10.3.1 emits a '-Warray-bounds' warning in callers of this function
  // about the sizes of nodeIds and coords not matching due to the runtime switch on dimension

  if(dimension == 1)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::index(),
      m,
      AXOM_LAMBDA(IndexType cellID) {
        IndexType cell_connectivity[2] = {cellID, cellID + 1};
        kernel(cellID, cell_connectivity, 2);
      });
  }
  else if(dimension == 2)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        IndexType cell_connectivity[4];

        for(int ii = 0; ii < 4; ++ii)
        {
          cell_connectivity[ii] = n0 + offsets[ii];
        }

        kernel(cellID, cell_connectivity, 4);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);

    for_all_cells_impl<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        IndexType cell_connectivity[8];

        for(int ii = 0; ii < 8; ++ii)
        {
          cell_connectivity[ii] = n0 + offsets[ii];
        }

        kernel(cellID, cell_connectivity, 8);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::nodeids,
                               const UnstructuredMesh<MIXED_SHAPE>& m,
                               KernelType&& kernel)
{
  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  auto cell_connectivity_h =
    axom::ArrayView<const IndexType>(m.getCellNodesArray(), m.getCellNodesSize());
  auto cell_offsets_h =
    axom::ArrayView<const IndexType>(m.getCellNodesOffsetsArray(), m.getNumberOfCells() + 1);

  // Move cell connectivity and cell offsets onto device
  axom::Array<IndexType> cell_connectivity_d = on_device
    ? axom::Array<IndexType>(cell_connectivity_h, device_allocator)
    : axom::Array<IndexType>();
  auto cell_connectivity_view = on_device ? cell_connectivity_d.view() : cell_connectivity_h;

  axom::Array<IndexType> cell_offsets_d =
    on_device ? axom::Array<IndexType>(cell_offsets_h, device_allocator) : axom::Array<IndexType>();
  auto cell_offsets_view = on_device ? cell_offsets_d.view() : cell_offsets_h;

  for_all_cells_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType cellID) {
      const IndexType N = cell_offsets_view[cellID + 1] - cell_offsets_view[cellID];
      kernel(cellID, &cell_connectivity_view[cell_offsets_view[cellID]], N);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::nodeids,
                               const UnstructuredMesh<SINGLE_SHAPE>& m,
                               KernelType&& kernel)
{
  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  auto cell_connectivity_h =
    axom::ArrayView<const IndexType>(m.getCellNodesArray(), m.getCellNodesSize());

  // Move cell connectivity onto device
  axom::Array<IndexType> cell_connectivity_d = on_device
    ? axom::Array<IndexType>(cell_connectivity_h, device_allocator)
    : axom::Array<IndexType>();
  auto cell_connectivity_view = on_device ? cell_connectivity_d.view() : cell_connectivity_h;

  const IndexType stride = m.getNumberOfCellNodes();

  for_all_cells_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType cellID) {
      kernel(cellID, &cell_connectivity_view[cellID * stride], stride);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::nodeids, const Mesh& m, KernelType&& kernel)
{
  if(m.isStructured())
  {
    const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::nodeids(), sm, std::forward<KernelType>(kernel));
  }
  else if(m.hasMixedCellTypes())
  {
    const UnstructuredMesh<MIXED_SHAPE>& um = static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::nodeids(), um, std::forward<KernelType>(kernel));
  }
  else
  {
    const UnstructuredMesh<SINGLE_SHAPE>& um = static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::nodeids(), um, std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::faceids, const StructuredMesh& m, KernelType&& kernel)
{
  const int dimension = m.getDimension();
  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);

  if(dimension == 2)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType AXOM_UNUSED_PARAM(i), IndexType j) {
        IndexType faces[4];

        /* The I_DIRECTION faces */
        faces[0] = cellID + j;
        faces[1] = faces[0] + 1;

        /* The J_DIRECTION faces */
        faces[2] = cellID + numIFaces;
        faces[3] = faces[2] + ICellResolution;

        kernel(cellID, faces, 4);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);

    const IndexType JCellResolution = m.getCellResolution(J_DIRECTION);
    const IndexType totalIJfaces = numIFaces + m.getTotalNumFaces(J_DIRECTION);
    const IndexType cellKp = m.cellKp();

    for_all_cells_impl<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType AXOM_UNUSED_PARAM(i), IndexType j, IndexType k) {
        IndexType faces[6];

        /* The I_DIRECTION faces */
        faces[0] = cellID + j + JCellResolution * k;
        faces[1] = faces[0] + 1;

        /* The J_DIRECTION faces */
        faces[2] = cellID + numIFaces + ICellResolution * k;
        faces[3] = faces[2] + ICellResolution;

        /* The K_DIRECTION faces */
        faces[4] = cellID + totalIJfaces;
        faces[5] = faces[4] + cellKp;

        kernel(cellID, faces, 6);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::faceids,
                               const UnstructuredMesh<SINGLE_SHAPE>& m,
                               KernelType&& kernel)
{
  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  const IndexType num_faces = m.getNumberOfCellFaces();

  // extract cells_to_faces into an axom::ArrayView
  auto cells_to_faces_h =
    axom::ArrayView<const IndexType>(m.getCellFacesArray(), m.getCellNodesSize());

  // Move cells_to_faces values onto device
  axom::Array<IndexType> cells_to_faces_d = on_device
    ? axom::Array<IndexType>(cells_to_faces_h, device_allocator)
    : axom::Array<IndexType>();
  auto cells_to_faces_v = on_device ? cells_to_faces_d.view() : cells_to_faces_h;

  for_all_cells_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType cellID) {
      kernel(cellID, cells_to_faces_v.data() + cellID * num_faces, num_faces);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::faceids,
                               const UnstructuredMesh<MIXED_SHAPE>& m,
                               KernelType&& kernel)
{
  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  // extract cells_to_faces and offsets into an axom::ArrayView
  auto cells_to_faces_h =
    axom::ArrayView<const IndexType>(m.getCellFacesArray(), m.getCellNodesSize());
  auto offsets_h =
    axom::ArrayView<const IndexType>(m.getCellFacesOffsetsArray(), m.getNumberOfCells() + 1);

  // Move cells_to_faces and offsets values onto device
  axom::Array<IndexType> cells_to_faces_d = on_device
    ? axom::Array<IndexType>(cells_to_faces_h, device_allocator)
    : axom::Array<IndexType>();
  auto cells_to_faces_v = on_device ? cells_to_faces_d.view() : cells_to_faces_h;

  axom::Array<IndexType> offsets_d =
    on_device ? axom::Array<IndexType>(offsets_h, device_allocator) : axom::Array<IndexType>();
  auto offsets_v = on_device ? offsets_d.view() : offsets_h;

  for_all_cells_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType cellID) {
      const IndexType num_faces = offsets_v[cellID + 1] - offsets_v[cellID];
      kernel(cellID, cells_to_faces_v.data() + offsets_v[cellID], num_faces);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::faceids, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() == 1,
                "For all cells with face IDs only supported for 2D and 3D meshes");

  if(m.isStructured())
  {
    const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::faceids(), sm, std::forward<KernelType>(kernel));
  }
  else if(m.hasMixedCellTypes())
  {
    const UnstructuredMesh<MIXED_SHAPE>& um = static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::faceids(), um, std::forward<KernelType>(kernel));
  }
  else
  {
    const UnstructuredMesh<SINGLE_SHAPE>& um = static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::faceids(), um, std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::coords, const UniformMesh& m, KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  const int dimension = m.getDimension();
  const double* origin = m.getOrigin();
  const double* spacing = m.getSpacing();
  const IndexType nodeJp = m.nodeJp();
  const IndexType nodeKp = m.nodeKp();

  const double x0 = origin[0];
  const double dx = spacing[0];

  const double y0 = origin[1];
  const double dy = spacing[1];

  const double z0 = origin[2];
  const double dz = spacing[2];

  // Note: gcc@10.3.1 emits a '-Warray-bounds' warning in callers of this function
  // about the sizes of nodeIds and coords not matching due to the runtime switch on dimension

  if(dimension == 1)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::index(),
      m,
      AXOM_LAMBDA(IndexType cellID) {
        const IndexType nodeIDs[2] = {cellID, cellID + 1};
        double coords[2] = {x0 + nodeIDs[0] * dx, x0 + nodeIDs[1] * dx};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else if(dimension == 2)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeJp, n0 + nodeJp};

        double coords[8] = {x0 + i * dx,
                            y0 + j * dy,
                            x0 + (i + 1) * dx,
                            y0 + j * dy,
                            x0 + (i + 1) * dx,
                            y0 + (j + 1) * dy,
                            x0 + i * dx,
                            y0 + (j + 1) * dy};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);
    for_all_cells_impl<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[8] = {n0,
                                      n0 + 1,
                                      n0 + 1 + nodeJp,
                                      n0 + nodeJp,
                                      n0 + nodeKp,
                                      n0 + 1 + nodeKp,
                                      n0 + 1 + nodeJp + nodeKp,
                                      n0 + nodeJp + nodeKp};

        double coords[24] = {
          x0 + i * dx,       y0 + j * dy,       z0 + k * dz,       x0 + (i + 1) * dx,
          y0 + j * dy,       z0 + k * dz,       x0 + (i + 1) * dx, y0 + (j + 1) * dy,
          z0 + k * dz,       x0 + i * dx,       y0 + (j + 1) * dy, z0 + k * dz,
          x0 + i * dx,       y0 + j * dy,       z0 + (k + 1) * dz, x0 + (i + 1) * dx,
          y0 + j * dy,       z0 + (k + 1) * dz, x0 + (i + 1) * dx, y0 + (j + 1) * dy,
          z0 + (k + 1) * dz, x0 + i * dx,       y0 + (j + 1) * dy, z0 + (k + 1) * dz};

        numerics::Matrix<double> coordsMatrix(dimension, 8, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::coords, const RectilinearMesh& m, KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  const int dimension = m.getDimension();
  const IndexType nodeJp = m.nodeJp();
  const IndexType nodeKp = m.nodeKp();

  // Extract x coordinate values into an axom::ArrayView
  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE),
                                                m.getNodeResolution(X_COORDINATE));

  // Move x values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, device_allocator) : axom::Array<double>();
  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;

  if(dimension == 1)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::index(),
      m,
      AXOM_LAMBDA(IndexType cellID) {
        const IndexType nodeIDs[2] = {cellID, cellID + 1};
        double coords[2] = {x_vals_view[nodeIDs[0]], x_vals_view[nodeIDs[1]]};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else if(dimension == 2)
  {
    // Extract y coordinate values into an axom::ArrayView
    auto y_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE),
                                                  m.getNodeResolution(Y_COORDINATE));

    // Move y values onto device
    axom::Array<double> y_vals_d =
      on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();
    auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

    for_all_cells_impl<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeJp, n0 + nodeJp};

        double coords[8] = {x_vals_view[i],
                            y_vals_view[j],
                            x_vals_view[i + 1],
                            y_vals_view[j],
                            x_vals_view[i + 1],
                            y_vals_view[j + 1],
                            x_vals_view[i],
                            y_vals_view[j + 1]};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);

    // Extract yz coordinate values into an axom::ArrayView
    auto y_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE),
                                                  m.getNodeResolution(Y_COORDINATE));
    auto z_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Z_COORDINATE),
                                                  m.getNodeResolution(Z_COORDINATE));

    // Move yz values onto device
    axom::Array<double> y_vals_d =
      on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();
    axom::Array<double> z_vals_d =
      on_device ? axom::Array<double>(z_vals_h, device_allocator) : axom::Array<double>();

    auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;
    auto z_vals_view = on_device ? z_vals_d.view() : z_vals_h;

    for_all_cells_impl<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType cellID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[8] = {n0,
                                      n0 + 1,
                                      n0 + 1 + nodeJp,
                                      n0 + nodeJp,
                                      n0 + nodeKp,
                                      n0 + 1 + nodeKp,
                                      n0 + 1 + nodeJp + nodeKp,
                                      n0 + nodeJp + nodeKp};

        double coords[24] = {
          x_vals_view[i],     y_vals_view[j],     z_vals_view[k],     x_vals_view[i + 1],
          y_vals_view[j],     z_vals_view[k],     x_vals_view[i + 1], y_vals_view[j + 1],
          z_vals_view[k],     x_vals_view[i],     y_vals_view[j + 1], z_vals_view[k],
          x_vals_view[i],     y_vals_view[j],     z_vals_view[k + 1], x_vals_view[i + 1],
          y_vals_view[j],     z_vals_view[k + 1], x_vals_view[i + 1], y_vals_view[j + 1],
          z_vals_view[k + 1], x_vals_view[i],     y_vals_view[j + 1], z_vals_view[k + 1]};

        numerics::Matrix<double> coordsMatrix(dimension, 8, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
struct for_all_cell_nodes_functor
{
  template <typename ExecPolicy, typename MeshType, typename KernelType>
  inline void operator()(ExecPolicy AXOM_UNUSED_PARAM(policy),
                         const MeshType& m,
                         KernelType&& kernel) const
  {
    constexpr bool valid_mesh_type = std::is_base_of<Mesh, MeshType>::value;
    AXOM_STATIC_ASSERT(valid_mesh_type);

    for_all_cells_impl<ExecPolicy>(xargs::nodeids(), m, std::forward<KernelType>(kernel));
  }
};

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells_impl(xargs::coords, const CurvilinearMesh& m, KernelType&& kernel)
{
  const int dimension = m.getDimension();
  if(dimension == 1)
  {
    for_all_coords<ExecPolicy, 1, 2>(for_all_cell_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
  else if(dimension == 2)
  {
    for_all_coords<ExecPolicy, 2, 4>(for_all_cell_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_coords<ExecPolicy, 3, 8>(for_all_cell_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType, Topology TOPO>
inline void for_all_cells_impl(xargs::coords, const UnstructuredMesh<TOPO>& m, KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();
  IndexType coordinate_size = m.getNumberOfNodes();

  const int dimension = m.getDimension();

  // Extract x coordinate values into an axom::ArrayView
  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE), coordinate_size);

  // Move x values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, device_allocator) : axom::Array<double>();
  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;

  if(dimension == 1)
  {
    for_all_cells_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType cellID, const IndexType* nodeIDs, IndexType AXOM_UNUSED_PARAM(numNodes)) {
        double coords[2] = {x_vals_view[nodeIDs[0]], x_vals_view[nodeIDs[1]]};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else if(dimension == 2)
  {
    // Extract y coordinate values into an axom::ArrayView
    auto y_vals_h =
      axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE), coordinate_size);

    // Move y values onto device
    axom::Array<double> y_vals_d =
      on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();
    auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

    for_all_cells_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType cellID, const IndexType* nodeIDs, IndexType numNodes) {
        double coords[2 * MAX_CELL_NODES];
        for(int i = 0; i < numNodes; ++i)
        {
          const IndexType nodeID = nodeIDs[i];
          coords[2 * i] = x_vals_view[nodeID];
          coords[2 * i + 1] = y_vals_view[nodeID];
        }

        numerics::Matrix<double> coordsMatrix(dimension, numNodes, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    SLIC_ASSERT(dimension == 3);

    // Extract coordinate values into an axom::ArrayView
    auto y_vals_h =
      axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE), coordinate_size);
    auto z_vals_h =
      axom::ArrayView<const double>(m.getCoordinateArray(Z_COORDINATE), coordinate_size);

    // Move yz values onto device
    axom::Array<double> y_vals_d =
      on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();
    auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

    axom::Array<double> z_vals_d =
      on_device ? axom::Array<double>(z_vals_h, device_allocator) : axom::Array<double>();
    auto z_vals_view = on_device ? z_vals_d.view() : z_vals_h;

    for_all_cells_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType cellID, const IndexType* nodeIDs, IndexType numNodes) {
        double coords[3 * MAX_CELL_NODES];
        for(int i = 0; i < numNodes; ++i)
        {
          const IndexType nodeID = nodeIDs[i];
          coords[3 * i] = x_vals_view[nodeID];
          coords[3 * i + 1] = y_vals_view[nodeID];
          coords[3 * i + 2] = z_vals_view[nodeID];
        }

        numerics::Matrix<double> coordsMatrix(dimension, numNodes, coords, NO_COPY);
        kernel(cellID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_cells(xargs::coords, const Mesh& m, KernelType&& kernel)
{
  if(m.getMeshType() == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::coords(), um, std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == STRUCTURED_RECTILINEAR_MESH)
  {
    const RectilinearMesh& rm = static_cast<const RectilinearMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::coords(), rm, std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == STRUCTURED_CURVILINEAR_MESH)
  {
    const CurvilinearMesh& cm = static_cast<const CurvilinearMesh&>(m);
    for_all_cells_impl<ExecPolicy>(xargs::coords(), cm, std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == UNSTRUCTURED_MESH)
  {
    if(m.hasMixedCellTypes())
    {
      const UnstructuredMesh<MIXED_SHAPE>& um = static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);

      for_all_cells_impl<ExecPolicy>(xargs::coords(), um, std::forward<KernelType>(kernel));
    }
    else
    {
      const UnstructuredMesh<SINGLE_SHAPE>& um =
        static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);

      for_all_cells_impl<ExecPolicy>(xargs::coords(), um, std::forward<KernelType>(kernel));
    }
  }
  else
  {
    SLIC_ERROR("Unsupported mesh type: " << m.getMeshType());
  }
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_FOR_ALL_CELLS_HPP_ */
