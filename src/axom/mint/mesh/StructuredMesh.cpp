// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/mesh/StructuredMesh.hpp"
#include "axom/mint/mesh/MeshTypes.hpp"
#include "axom/mint/mesh/blueprint.hpp"

#include <cstring> /* for memcpy() */
#include "axom/core/NumericLimits.hpp"

namespace axom
{
namespace mint
{
//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
bool validStructuredMeshType(int type)
{
  return ((type == STRUCTURED_CURVILINEAR_MESH) || (type == STRUCTURED_RECTILINEAR_MESH) ||
          (type == STRUCTURED_UNIFORM_MESH));
}

inline int dim(const IndexType& AXOM_UNUSED_PARAM(Ni), const IndexType& Nj, const IndexType& Nk)
{
  return ((Nk >= 1) ? 3 : ((Nj >= 1) ? 2 : 1));
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------

void StructuredMesh::setExtent(int ndims, const int64* extent)
{
  SLIC_ASSERT(0 < ndims && ndims <= 3);
  SLIC_ASSERT(extent != nullptr);

  std::memset(m_node_extent, 0, sizeof(m_node_extent));

  int64 numNodes = 1;
  for(int dim = 0; dim < ndims; ++dim)
  {
    const int64 min = extent[2 * dim];
    const int64 max = extent[2 * dim + 1];
    SLIC_ASSERT(min <= max);
    m_node_extent[2 * dim] = min;
    m_node_extent[2 * dim + 1] = max;
    numNodes *= max - min + 1;
  }

  SLIC_ASSERT(numNodes == getNumberOfNodes());
  AXOM_UNUSED_VAR(numNodes);

#ifdef AXOM_MINT_USE_SIDRE
  if(hasSidreGroup())
  {
    blueprint::setExtent(getCoordsetGroup(), m_node_extent);
  }
#endif
}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh(int meshType, IndexType Ni, IndexType Nj, IndexType Nk)
  : Mesh(dim(Ni, Nj, Nk), meshType)
{
  SLIC_ERROR_IF(!validStructuredMeshType(m_type), "invalid structured mesh type!");

  SLIC_ERROR_IF(Ni <= 1, "Ni must be greater or equal to 2");
  m_node_dims[0] = Ni;
  if(m_ndims > 1)
  {
    SLIC_ERROR_IF(Nj <= 1, "Nj must be greater or equal to 2");
    m_node_dims[1] = Nj;
  }
  if(m_ndims > 2)
  {
    SLIC_ERROR_IF(Nk <= 1, "Nk must be greater or equal to 2");
    m_node_dims[2] = Nk;
  }

  structuredInit();
}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh(sidre::Group* group, const std::string& topo) : Mesh(group, topo)
{
  SLIC_ERROR_IF(!validStructuredMeshType(m_type), "invalid structured mesh type!");

  blueprint::getStructuredMeshProperties(m_ndims, m_node_dims, m_node_extent, getCoordsetGroup());
  structuredInit();
}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh(int meshType,
                               IndexType Ni,
                               IndexType Nj,
                               IndexType Nk,
                               sidre::Group* group,
                               const std::string& topo,
                               const std::string& coordset)
  : Mesh(dim(Ni, Nj, Nk), meshType, group, topo, coordset)
{
  SLIC_ERROR_IF(!validStructuredMeshType(m_type), "invalid structured mesh type!");

  SLIC_ERROR_IF(Ni <= 1, "Ni must be greater or equal to 2");
  m_node_dims[0] = Ni;
  if(m_ndims > 1)
  {
    SLIC_ERROR_IF(Nj <= 1, "Nj must be greater or equal to 2");
    m_node_dims[1] = Nj;
  }
  if(m_ndims > 2)
  {
    SLIC_ERROR_IF(Nk <= 1, "Nk must be greater or equal to 2");
    m_node_dims[2] = Nk;
  }

  std::string topo_type;
  if(meshType == STRUCTURED_UNIFORM_MESH)
  {
    topo_type = "uniform";
  }
  else if(meshType == STRUCTURED_RECTILINEAR_MESH)
  {
    topo_type = "rectilinear";
  }
  else
  {
    topo_type = "structured";
  }

  blueprint::initializeTopologyGroup(m_group, m_topology, m_coordset, topo_type);
  SLIC_ERROR_IF(!blueprint::isValidTopologyGroup(getTopologyGroup()), "invalid topology group!");

  blueprint::setStructuredMeshProperties(m_ndims, m_node_dims, m_node_extent, getCoordsetGroup());
  structuredInit();
}

#endif

//------------------------------------------------------------------------------
void StructuredMesh::structuredInit()
{
  for(int dim = 0; dim < m_ndims; ++dim)
  {
    SLIC_ERROR_IF(getNodeResolution(dim) < 2, "invalid extent");
  }

  /* Initialize the node meta data. */
  m_node_jp = (m_ndims > 1) ? getNodeResolution(0) : axom::numeric_limits<IndexType>::max();
  m_node_kp =
    (m_ndims > 2) ? m_node_jp * getNodeResolution(1) : axom::numeric_limits<IndexType>::max();

  /* Initialize the cell meta data */
  for(int dim = 0; dim < m_ndims; ++dim)
  {
    m_cell_dims[dim] = getNodeResolution(dim) - 1;
  }

  m_cell_jp = (m_ndims > 1) ? getCellResolution(0) : axom::numeric_limits<IndexType>::max();
  m_cell_kp =
    (m_ndims > 2) ? m_cell_jp * getCellResolution(1) : axom::numeric_limits<IndexType>::max();

  /* Build the cell to node offsets. */
  m_cell_node_offsets[0] = 0;
  m_cell_node_offsets[1] = 1;
  m_cell_node_offsets[2] = (m_ndims == 1) ? axom::numeric_limits<IndexType>::min() : (1 + nodeJp());
  m_cell_node_offsets[3] = nodeJp();

  m_cell_node_offsets[4] = nodeKp();
  m_cell_node_offsets[5] = (m_ndims < 3) ? axom::numeric_limits<IndexType>::min() : (1 + nodeKp());
  m_cell_node_offsets[6] =
    (m_ndims < 3) ? axom::numeric_limits<IndexType>::min() : (1 + nodeJp() + nodeKp());
  m_cell_node_offsets[7] =
    (m_ndims < 3) ? axom::numeric_limits<IndexType>::min() : (nodeJp() + nodeKp());

  /* Initialize the face meta data */
  if(m_ndims == 2)
  {
    m_total_faces[0] = getNodeResolution(0) * getCellResolution(1);
    m_total_faces[1] = getCellResolution(0) * getNodeResolution(1);
  }
  else if(m_ndims == 3)
  {
    m_total_faces[0] = getNodeResolution(0) * getCellResolution(1) * getCellResolution(2);
    m_total_faces[1] = getCellResolution(0) * getNodeResolution(1) * getCellResolution(2);
    m_total_faces[2] = getCellResolution(0) * getCellResolution(1) * getNodeResolution(2);
  }

  m_total_IJ_faces = m_total_faces[0] + m_total_faces[1];
  m_num_I_faces_in_k_slice = getNodeResolution(0) * getCellResolution(1);
  m_num_J_faces_in_k_slice = getCellResolution(0) * getNodeResolution(1);

  /* Initialize the edge meta data */
  if(m_ndims == 3)
  {
    m_num_edges = getCellResolution(0) * getNodeResolution(1) * getNodeResolution(2) +
      getCellResolution(0) * getNodeResolution(2) * getNodeResolution(1) +
      getCellResolution(1) * getNodeResolution(2) * getNodeResolution(0);
  }

  /* Initialize the fields */
  m_mesh_fields[NODE_CENTERED]->setResizeRatio(0.0);
  m_mesh_fields[CELL_CENTERED]->setResizeRatio(0.0);
  m_mesh_fields[FACE_CENTERED]->setResizeRatio(0.0);
  m_mesh_fields[EDGE_CENTERED]->setResizeRatio(0.0);

  m_mesh_fields[NODE_CENTERED]->resize(getNumberOfNodes());
  m_mesh_fields[CELL_CENTERED]->resize(getNumberOfCells());
  m_mesh_fields[FACE_CENTERED]->resize(getNumberOfFaces());
  m_mesh_fields[EDGE_CENTERED]->resize(getNumberOfEdges());
}

} /* end namespace mint */
} /* end namespace axom */
