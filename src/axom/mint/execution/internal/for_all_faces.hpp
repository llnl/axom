// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FOR_ALL_FACES_HPP_
#define MINT_FOR_ALL_FACES_HPP_

// Axom core includes
#include "axom/config.hpp"                          // compile time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits
#include "axom/core/execution/for_all.hpp"          // for axom::for_all

// mint includes
#include "axom/mint/execution/xargs.hpp"             // for xargs
#include "axom/mint/config.hpp"                      // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"                   // for Mesh
#include "axom/mint/mesh/StructuredMesh.hpp"         // for StructuredMesh
#include "axom/mint/mesh/UniformMesh.hpp"            // for UniformMesh
#include "axom/mint/mesh/RectilinearMesh.hpp"        // for RectilinearMesh
#include "axom/mint/mesh/CurvilinearMesh.hpp"        // for CurvilinearMesh
#include "axom/mint/execution/internal/helpers.hpp"  // for for_all_coords
#include "axom/core/execution/nested_for_exec.hpp"

#include "axom/core/numerics/Matrix.hpp"  // for Matrix

namespace axom
{
namespace mint
{
namespace internal
{
namespace helpers
{
//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_I_faces(xargs::ij, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "Mesh must be 2D.");

  const IndexType INodeResolution = m.getNodeResolution(I_DIRECTION);
  const IndexType Ni = INodeResolution;
  const IndexType Nj = m.getCellResolution(J_DIRECTION);

  axom::StackArray<IndexType, 2> i_range {{0, Ni}}, j_range {{0, Nj}};
  axom::for_all<ExecPolicy>(
    i_range,
    j_range,
    AXOM_LAMBDA(IndexType i, IndexType j) {
      const IndexType faceID = i + j * INodeResolution;
      kernel(faceID, i, j);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_I_faces(xargs::ijk, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "Mesh must be a 3D.");

  const IndexType INodeResolution = m.getNodeResolution(I_DIRECTION);
  const IndexType numIFacesInKSlice = INodeResolution * m.getCellResolution(J_DIRECTION);
  const IndexType Ni = INodeResolution;
  const IndexType Nj = m.getCellResolution(J_DIRECTION);
  const IndexType Nk = m.getCellResolution(K_DIRECTION);

  axom::StackArray<IndexType, 2> i_range {{0, Ni}}, j_range {{0, Nj}}, k_range {{0, Nk}};
  axom::for_all<ExecPolicy>(
    i_range,
    j_range,
    k_range,
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType faceID = i + j * INodeResolution + k * numIFacesInKSlice;
      kernel(faceID, i, j, k);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_J_faces(xargs::ij, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 2, "Mesh must be 2D.");

  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);
  const IndexType Ni = ICellResolution;
  const IndexType Nj = m.getNodeResolution(J_DIRECTION);

  axom::StackArray<IndexType, 2> i_range {{0, Ni}}, j_range {{0, Nj}};
  axom::for_all<ExecPolicy>(
    i_range,
    j_range,
    AXOM_LAMBDA(IndexType i, IndexType j) {
      const IndexType faceID = numIFaces + i + j * ICellResolution;
      kernel(faceID, i, j);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_J_faces(xargs::ijk, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "Mesh must be 3D.");

  const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);
  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType numJFacesInKSlice = ICellResolution * m.getNodeResolution(J_DIRECTION);
  const IndexType Ni = ICellResolution;
  const IndexType Nj = m.getNodeResolution(J_DIRECTION);
  const IndexType Nk = m.getCellResolution(K_DIRECTION);

  axom::StackArray<IndexType, 2> i_range {{0, Ni}}, j_range {{0, Nj}}, k_range {{0, Nk}};
  axom::for_all<ExecPolicy>(
    i_range,
    j_range,
    k_range,
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType jp = j * ICellResolution;
      const IndexType kp = k * numJFacesInKSlice;
      const IndexType faceID = numIFaces + i + jp + kp;
      kernel(faceID, i, j, k);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_K_faces(xargs::ijk, const StructuredMesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() != 3, "Mesh must be 3D.");

  const IndexType numIJFaces = m.getTotalNumFaces(I_DIRECTION) + m.getTotalNumFaces(J_DIRECTION);
  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType cellKp = m.cellKp();
  const IndexType Ni = ICellResolution;
  const IndexType Nj = m.getCellResolution(J_DIRECTION);
  const IndexType Nk = m.getNodeResolution(K_DIRECTION);

  axom::StackArray<IndexType, 2> i_range {{0, Ni}}, j_range {{0, Nj}}, k_range {{0, Nk}};
  axom::for_all<ExecPolicy>(
    i_range,
    j_range,
    k_range,
    AXOM_LAMBDA(IndexType i, IndexType j, IndexType k) {
      const IndexType jp = j * ICellResolution;
      const IndexType kp = k * cellKp;
      const IndexType faceID = numIJFaces + i + jp + kp;
      kernel(faceID, i, j, k);
    });
}

} /* namespace helpers */

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::index, const Mesh& m, KernelType&& kernel)
{
  const IndexType numFaces = m.getNumberOfFaces();
  axom::for_all<ExecPolicy>(numFaces, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces(xargs::index, const Mesh& m, KernelType&& kernel)
{
  return for_all_faces_impl<ExecPolicy>(xargs::index(), m, std::forward<KernelType>(kernel));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::nodeids, const StructuredMesh& m, KernelType&& kernel)
{
  const IndexType dimension = m.getDimension();
  const IndexType* offsets = m.getCellNodeOffsetsArray();
  const IndexType cellNodeOffset3 = offsets[3];

  if(dimension == 2)
  {
    const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);

    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType AXOM_UNUSED_PARAM(i), IndexType AXOM_UNUSED_PARAM(j)) {
        IndexType nodes[2];
        nodes[0] = faceID;
        nodes[1] = nodes[0] + cellNodeOffset3;
        kernel(faceID, nodes, 2);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType AXOM_UNUSED_PARAM(i), IndexType j) {
        const IndexType shiftedID = faceID - numIFaces;
        IndexType nodes[2];
        nodes[0] = shiftedID + j;
        nodes[1] = nodes[0] + 1;
        kernel(faceID, nodes, 2);
      });
  }
  else
  {
    SLIC_ERROR_IF(dimension != 3, "for_all_faces is only valid for 2 or 3D meshes.");

    const IndexType numIFaces = m.getTotalNumFaces(I_DIRECTION);
    const IndexType numIJFaces = numIFaces + m.getTotalNumFaces(J_DIRECTION);
    const IndexType INodeResolution = m.getNodeResolution(I_DIRECTION);
    const IndexType JNodeResolution = m.getNodeResolution(J_DIRECTION);
    const IndexType KFaceNodeStride =
      m.getCellResolution(I_DIRECTION) + m.getCellResolution(J_DIRECTION) + 1;

    const IndexType cellNodeOffset2 = offsets[2];
    const IndexType cellNodeOffset4 = offsets[4];
    const IndexType cellNodeOffset5 = offsets[5];
    const IndexType cellNodeOffset7 = offsets[7];

    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID,
                  IndexType AXOM_UNUSED_PARAM(i),
                  IndexType AXOM_UNUSED_PARAM(j),
                  IndexType k) {
        IndexType nodes[4];
        nodes[0] = faceID + k * INodeResolution;
        nodes[1] = nodes[0] + cellNodeOffset4;
        nodes[2] = nodes[0] + cellNodeOffset7;
        nodes[3] = nodes[0] + cellNodeOffset3;
        kernel(faceID, nodes, 4);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType AXOM_UNUSED_PARAM(i), IndexType j, IndexType k) {
        const IndexType shiftedID = faceID - numIFaces;
        IndexType nodes[4];
        nodes[0] = shiftedID + j + k * JNodeResolution;
        nodes[1] = nodes[0] + 1;
        nodes[2] = nodes[0] + cellNodeOffset5;
        nodes[3] = nodes[0] + cellNodeOffset4;
        kernel(faceID, nodes, 4);
      });

    helpers::for_all_K_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType AXOM_UNUSED_PARAM(i), IndexType j, IndexType k) {
        const IndexType shiftedID = faceID - numIJFaces;
        IndexType nodes[4];
        nodes[0] = shiftedID + j + k * KFaceNodeStride;
        nodes[1] = nodes[0] + 1;
        nodes[2] = nodes[0] + cellNodeOffset2;
        nodes[3] = nodes[0] + cellNodeOffset3;
        kernel(faceID, nodes, 4);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::nodeids,
                               const UnstructuredMesh<SINGLE_SHAPE>& m,
                               KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getNumberOfFaces() <= 0,
                "No faces in the mesh, perhaps you meant to call "
                  << "UnstructuredMesh::initializeFaceConnectivity first.");

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  auto faces_to_nodes_h =
    axom::ArrayView<const IndexType>(m.getFaceNodesArray(), m.getFaceNodesSize());

  // Move faces to nodes onto device
  axom::Array<IndexType> faces_to_nodes_d = on_device
    ? axom::Array<IndexType>(faces_to_nodes_h, device_allocator)
    : axom::Array<IndexType>();

  auto faces_to_nodes_view = on_device ? faces_to_nodes_d.view() : faces_to_nodes_h;

  const IndexType num_nodes = m.getNumberOfFaceNodes();

  for_all_faces_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType faceID) {
      kernel(faceID, faces_to_nodes_view.data() + faceID * num_nodes, num_nodes);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::nodeids,
                               const UnstructuredMesh<MIXED_SHAPE>& m,
                               KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getNumberOfFaces() <= 0,
                "No faces in the mesh, perhaps you meant to call "
                  << "UnstructuredMesh::initializeFaceConnectivity first.");

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  auto faces_to_nodes_h =
    axom::ArrayView<const IndexType>(m.getFaceNodesArray(), m.getFaceNodesSize());
  auto offsets_h =
    axom::ArrayView<const IndexType>(m.getFaceNodesOffsetsArray(), m.getNumberOfFaces() + 1);

  // Move faces to nodes and offsets onto device
  axom::Array<IndexType> faces_to_nodes_d = on_device
    ? axom::Array<IndexType>(faces_to_nodes_h, device_allocator)
    : axom::Array<IndexType>();
  axom::Array<IndexType> offsets_d =
    on_device ? axom::Array<IndexType>(offsets_h, device_allocator) : axom::Array<IndexType>();

  auto faces_to_nodes_view = on_device ? faces_to_nodes_d.view() : faces_to_nodes_h;
  auto offsets_view = on_device ? offsets_d.view() : offsets_h;

  for_all_faces_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType faceID) {
      const IndexType num_nodes = offsets_view[faceID + 1] - offsets_view[faceID];
      kernel(faceID, faces_to_nodes_view.data() + offsets_view[faceID], num_nodes);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces(xargs::nodeids, const Mesh& m, KernelType&& kernel)
{
  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  if(m.isStructured())
  {
    const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::nodeids(), sm, std::forward<KernelType>(kernel));
  }
  else if(m.hasMixedCellTypes())
  {
    const UnstructuredMesh<MIXED_SHAPE>& um = static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::nodeids(), um, std::forward<KernelType>(kernel));
  }
  else
  {
    const UnstructuredMesh<SINGLE_SHAPE>& um = static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::nodeids(), um, std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::cellids, const StructuredMesh& m, KernelType&& kernel)
{
  const IndexType ICellResolution = m.getCellResolution(I_DIRECTION);
  const IndexType JCellResolution = m.getCellResolution(J_DIRECTION);
  const IndexType cellJp = m.cellJp();
  const int dimension = m.getDimension();

  if(dimension == 2)
  {
    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        IndexType cellIDTwo = i + j * cellJp;
        IndexType cellIDOne = cellIDTwo - 1;
        if(i == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(i == ICellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        IndexType cellIDTwo = i + j * cellJp;
        IndexType cellIDOne = cellIDTwo - cellJp;
        if(j == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(j == JCellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });
  }
  else
  {
    SLIC_ERROR_IF(dimension != 3, "for_all_faces only valid for 2 or 3D meshes.");

    const IndexType KCellResolution = m.getCellResolution(K_DIRECTION);
    const IndexType cellKp = m.cellKp();

    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        IndexType cellIDTwo = i + j * cellJp + k * cellKp;
        IndexType cellIDOne = cellIDTwo - 1;
        if(i == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(i == ICellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        IndexType cellIDTwo = i + j * cellJp + k * cellKp;
        IndexType cellIDOne = cellIDTwo - cellJp;
        if(j == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(j == JCellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });

    helpers::for_all_K_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        IndexType cellIDTwo = i + j * cellJp + k * cellKp;
        IndexType cellIDOne = cellIDTwo - cellKp;
        if(k == 0)
        {
          cellIDOne = cellIDTwo;
          cellIDTwo = -1;
        }
        else if(k == KCellResolution)
        {
          cellIDTwo = -1;
        }

        kernel(faceID, cellIDOne, cellIDTwo);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, Topology TOPO, typename KernelType>
inline void for_all_faces_impl(xargs::cellids, const UnstructuredMesh<TOPO>& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getNumberOfFaces() <= 0,
                "No faces in the mesh, perhaps you meant to call "
                  << "UnstructuredMesh::initializeFaceConnectivity first.");

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  auto faces_to_cells_h =
    axom::ArrayView<const IndexType>(m.getFaceCellsArray(), 2 * m.getNumberOfFaces());

  // Move faces to cells onto device
  axom::Array<IndexType> faces_to_cells_d = on_device
    ? axom::Array<IndexType>(faces_to_cells_h, device_allocator)
    : axom::Array<IndexType>();

  auto faces_to_cells_view = on_device ? faces_to_cells_d.view() : faces_to_cells_h;

  for_all_faces_impl<ExecPolicy>(
    xargs::index(),
    m,
    AXOM_LAMBDA(IndexType faceID) {
      const IndexType offset = 2 * faceID;
      kernel(faceID, faces_to_cells_view[offset], faces_to_cells_view[offset + 1]);
    });
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces(xargs::cellids, const Mesh& m, KernelType&& kernel)
{
  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  if(m.isStructured())
  {
    const StructuredMesh& sm = static_cast<const StructuredMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::cellids(), sm, std::forward<KernelType>(kernel));
  }
  else if(m.hasMixedCellTypes())
  {
    const UnstructuredMesh<MIXED_SHAPE>& um = static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::cellids(), um, std::forward<KernelType>(kernel));
  }
  else
  {
    const UnstructuredMesh<SINGLE_SHAPE>& um = static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::cellids(), um, std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::coords, const UniformMesh& m, KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

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

  if(dimension == 2)
  {
    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = {n0, n0 + nodeJp};

        double coords[4] = {x0 + i * dx, y0 + j * dy, x0 + i * dx, y0 + (j + 1) * dy};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = {n0, n0 + 1};

        double coords[4] = {x0 + i * dx, y0 + j * dy, x0 + (i + 1) * dx, y0 + j * dy};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + nodeKp, n0 + nodeJp + nodeKp, n0 + nodeJp};

        double coords[12] = {x0 + i * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + i * dx,
                             y0 + j * dy,
                             z0 + (k + 1) * dz,
                             x0 + i * dx,
                             y0 + (j + 1) * dy,
                             z0 + (k + 1) * dz,
                             x0 + i * dx,
                             y0 + (j + 1) * dy,
                             z0 + k * dz};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeKp, n0 + nodeKp};

        double coords[12] = {x0 + i * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + (i + 1) * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + (i + 1) * dx,
                             y0 + j * dy,
                             z0 + (k + 1) * dz,
                             x0 + i * dx,
                             y0 + j * dy,
                             z0 + (k + 1) * dz};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_K_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeJp, n0 + nodeJp};

        double coords[12] = {x0 + i * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + (i + 1) * dx,
                             y0 + j * dy,
                             z0 + k * dz,
                             x0 + (i + 1) * dx,
                             y0 + (j + 1) * dy,
                             z0 + k * dz,
                             x0 + i * dx,
                             y0 + (j + 1) * dy,
                             z0 + k * dz};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::coords, const RectilinearMesh& m, KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  const int dimension = m.getDimension();
  const IndexType nodeJp = m.nodeJp();
  const IndexType nodeKp = m.nodeKp();

  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE),
                                                m.getNodeResolution(X_COORDINATE));
  auto y_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE),
                                                m.getNodeResolution(Y_COORDINATE));

  // Move xy values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, device_allocator) : axom::Array<double>();
  axom::Array<double> y_vals_d =
    on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();

  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;
  auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

  if(dimension == 2)
  {
    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = {n0, n0 + nodeJp};

        double coords[4] = {x_vals_view[i], y_vals_view[j], x_vals_view[i], y_vals_view[j + 1]};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ij(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j) {
        const IndexType n0 = i + j * nodeJp;
        const IndexType nodeIDs[2] = {n0, n0 + 1};

        double coords[4] = {x_vals_view[i], y_vals_view[j], x_vals_view[i + 1], y_vals_view[j]};

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    auto z_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Z_COORDINATE),
                                                  m.getNodeResolution(Z_COORDINATE));

    // Move z values onto device
    axom::Array<double> z_vals_d =
      on_device ? axom::Array<double>(z_vals_h, device_allocator) : axom::Array<double>();

    auto z_vals_view = on_device ? z_vals_d.view() : z_vals_h;

    helpers::for_all_I_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + nodeKp, n0 + nodeJp + nodeKp, n0 + nodeJp};

        double coords[12] = {x_vals_view[i],
                             y_vals_view[j],
                             z_vals_view[k],
                             x_vals_view[i],
                             y_vals_view[j],
                             z_vals_view[k + 1],
                             x_vals_view[i],
                             y_vals_view[j + 1],
                             z_vals_view[k + 1],
                             x_vals_view[i],
                             y_vals_view[j + 1],
                             z_vals_view[k]};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_J_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeKp, n0 + nodeKp};

        double coords[12] = {x_vals_view[i],
                             y_vals_view[j],
                             z_vals_view[k],
                             x_vals_view[i + 1],
                             y_vals_view[j],
                             z_vals_view[k],
                             x_vals_view[i + 1],
                             y_vals_view[j],
                             z_vals_view[k + 1],
                             x_vals_view[i],
                             y_vals_view[j],
                             z_vals_view[k + 1]};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });

    helpers::for_all_K_faces<ExecPolicy>(
      xargs::ijk(),
      m,
      AXOM_LAMBDA(IndexType faceID, IndexType i, IndexType j, IndexType k) {
        const IndexType n0 = i + j * nodeJp + k * nodeKp;
        const IndexType nodeIDs[4] = {n0, n0 + 1, n0 + 1 + nodeJp, n0 + nodeJp};

        double coords[12] = {x_vals_view[i],
                             y_vals_view[j],
                             z_vals_view[k],
                             x_vals_view[i + 1],
                             y_vals_view[j],
                             z_vals_view[k],
                             x_vals_view[i + 1],
                             y_vals_view[j + 1],
                             z_vals_view[k],
                             x_vals_view[i],
                             y_vals_view[j + 1],
                             z_vals_view[k]};

        numerics::Matrix<double> coordsMatrix(dimension, 4, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
struct for_all_face_nodes_functor
{
  template <typename ExecPolicy, typename MeshType, typename KernelType>
  inline void operator()(ExecPolicy AXOM_UNUSED_PARAM(policy),
                         const MeshType& m,
                         KernelType&& kernel) const
  {
    constexpr bool valid_mesh_type = std::is_base_of<Mesh, MeshType>::value;
    AXOM_STATIC_ASSERT(valid_mesh_type);

    for_all_faces_impl<ExecPolicy>(xargs::nodeids(), m, std::forward<KernelType>(kernel));
  }
};

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces_impl(xargs::coords, const CurvilinearMesh& m, KernelType&& kernel)
{
  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  const int dimension = m.getDimension();
  if(dimension == 2)
  {
    for_all_coords<ExecPolicy, 2, 2>(for_all_face_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
  else
  {
    for_all_coords<ExecPolicy, 3, 4>(for_all_face_nodes_functor(),
                                     m,
                                     std::forward<KernelType>(kernel));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType, Topology TOPO>
inline void for_all_faces_impl(xargs::coords, const UnstructuredMesh<TOPO>& m, KernelType&& kernel)
{
  constexpr bool NO_COPY = true;

  SLIC_ASSERT(m.getDimension() > 1 && m.getDimension() <= 3);

  constexpr bool on_device = axom::execution_space<ExecPolicy>::onDevice();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  const int dimension = m.getDimension();

  IndexType coordinate_size = m.getNumberOfNodes();

  auto x_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(X_COORDINATE), coordinate_size);
  auto y_vals_h = axom::ArrayView<const double>(m.getCoordinateArray(Y_COORDINATE), coordinate_size);

  // Move xy values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, device_allocator) : axom::Array<double>();
  axom::Array<double> y_vals_d =
    on_device ? axom::Array<double>(y_vals_h, device_allocator) : axom::Array<double>();

  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;
  auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

  if(dimension == 2)
  {
    for_all_faces_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType faceID, const IndexType* nodeIDs, IndexType numNodes) {
        double coords[2 * MAX_FACE_NODES];
        for(int i = 0; i < numNodes; ++i)
        {
          const IndexType nodeID = nodeIDs[i];
          coords[2 * i] = x_vals_view[nodeID];
          coords[2 * i + 1] = y_vals_view[nodeID];
        }

        numerics::Matrix<double> coordsMatrix(dimension, 2, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
  else
  {
    auto z_vals_h =
      axom::ArrayView<const double>(m.getCoordinateArray(Z_COORDINATE), coordinate_size);

    // Move z values onto device
    axom::Array<double> z_vals_d =
      on_device ? axom::Array<double>(z_vals_h, device_allocator) : axom::Array<double>();

    auto z_vals_view = on_device ? z_vals_d.view() : z_vals_h;

    for_all_faces_impl<ExecPolicy>(
      xargs::nodeids(),
      m,
      AXOM_LAMBDA(IndexType faceID, const IndexType* nodeIDs, IndexType numNodes) {
        double coords[3 * MAX_FACE_NODES];
        for(int i = 0; i < numNodes; ++i)
        {
          const IndexType nodeID = nodeIDs[i];
          coords[3 * i] = x_vals_view[nodeID];
          coords[3 * i + 1] = y_vals_view[nodeID];
          coords[3 * i + 2] = z_vals_view[nodeID];
        }

        numerics::Matrix<double> coordsMatrix(dimension, numNodes, coords, NO_COPY);
        kernel(faceID, coordsMatrix, nodeIDs);
      });
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, typename KernelType>
inline void for_all_faces(xargs::coords, const Mesh& m, KernelType&& kernel)
{
  SLIC_ERROR_IF(m.getDimension() <= 1 || m.getDimension() > 3, "Invalid dimension");

  if(m.getMeshType() == STRUCTURED_UNIFORM_MESH)
  {
    const UniformMesh& um = static_cast<const UniformMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::coords(), um, std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == STRUCTURED_RECTILINEAR_MESH)
  {
    const RectilinearMesh& rm = static_cast<const RectilinearMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::coords(), rm, std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == STRUCTURED_CURVILINEAR_MESH)
  {
    const CurvilinearMesh& cm = static_cast<const CurvilinearMesh&>(m);
    for_all_faces_impl<ExecPolicy>(xargs::coords(), cm, std::forward<KernelType>(kernel));
  }
  else if(m.getMeshType() == UNSTRUCTURED_MESH)
  {
    if(m.hasMixedCellTypes())
    {
      const UnstructuredMesh<MIXED_SHAPE>& um = static_cast<const UnstructuredMesh<MIXED_SHAPE>&>(m);

      for_all_faces_impl<ExecPolicy>(xargs::coords(), um, std::forward<KernelType>(kernel));
    }
    else
    {
      const UnstructuredMesh<SINGLE_SHAPE>& um =
        static_cast<const UnstructuredMesh<SINGLE_SHAPE>&>(m);

      for_all_faces_impl<ExecPolicy>(xargs::coords(), um, std::forward<KernelType>(kernel));
    }
  }
  else
  {
    SLIC_ERROR("Unknown mesh type.");
  }
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_FOR_ALL_FACES_HPP_ */
