// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mint/mesh/Mesh.hpp"

// axom includes
#include "axom/core/Types.hpp"

// mint includes
#include "axom/mint/mesh/CurvilinearMesh.hpp"
#include "axom/mint/mesh/FieldData.hpp"
#include "axom/mint/mesh/ParticleMesh.hpp"
#include "axom/mint/mesh/RectilinearMesh.hpp"
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
#endif

#include "axom/mint/mesh/blueprint.hpp"  // for blueprint functions

// C/C++ includes
#include <cstring>  // for strcmp()

namespace axom
{
namespace mint
{
#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
// IMPLEMENTATION OF FREE METHODS
//------------------------------------------------------------------------------
Mesh* getMesh(sidre::Group* group, const std::string& topo)
{
  SLIC_ERROR_IF(group == nullptr, "supplied group is null");

  int mesh_type = UNDEFINED_MESH;
  int dimension = -1;
  blueprint::getMeshTypeAndDimension(mesh_type, dimension, group, topo);

  Mesh* m = nullptr;
  switch(mesh_type)
  {
  case STRUCTURED_CURVILINEAR_MESH:
    m = new CurvilinearMesh(group, topo);
    break;
  case STRUCTURED_RECTILINEAR_MESH:
    m = new RectilinearMesh(group, topo);
    break;
  case STRUCTURED_UNIFORM_MESH:
    m = new UniformMesh(group, topo);
    break;
  case UNSTRUCTURED_MESH:
    if(blueprint::hasMixedCellTypes(group, topo))
    {
      m = new UnstructuredMesh<MIXED_SHAPE>(group, topo);
    }
    else
    {
      m = new UnstructuredMesh<SINGLE_SHAPE>(group, topo);
    }
    break;
  case PARTICLE_MESH:
    m = new ParticleMesh(group, topo);
    break;
  default:
    SLIC_ERROR("undefined mesh_type [" << mesh_type << "]\n");
  }  // END switch

  SLIC_ASSERT(m != nullptr);
  return m;
}

#endif /* AXOM_MINT_USE_SIDRE */

//------------------------------------------------------------------------------
// MESH IMPLEMENTATION
//------------------------------------------------------------------------------
Mesh::Mesh(int ndims, int type)
  : m_ndims(ndims)
  , m_type(type)
  , m_block_idx(-1)
  , m_part_idx(-1)
  , m_explicit_coords(false)
  , m_explicit_connectivity(false)
  , m_has_mixed_topology(false)
#ifdef AXOM_MINT_USE_SIDRE
  , m_group(nullptr)
  , m_topology()
  , m_coordset()
#endif
{
  SLIC_ERROR_IF(!validMeshType(), "invalid mesh type=" << m_type);
  SLIC_ERROR_IF(!validDimension(), "invalid mesh dimension=" << m_ndims);
  allocateFieldData();
}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
Mesh::Mesh(sidre::Group* group, const std::string& topo)
  : m_ndims(-1)
  , m_type(UNDEFINED_MESH)
  , m_block_idx(-1)
  , m_part_idx(-1)
  , m_explicit_coords(false)
  , m_explicit_connectivity(false)
  , m_has_mixed_topology(false)
  , m_group(group)
  , m_topology(topo)
  , m_coordset()
{
  SLIC_ERROR_IF(m_group == nullptr, "NULL sidre group");
  SLIC_ERROR_IF(!blueprint::isValidRootGroup(m_group), "root group does not conform to blueprint");

  blueprint::getMeshTypeAndDimension(m_type, m_ndims, m_group, m_topology);
  m_topology = getTopologyGroup()->getName();
  m_coordset = blueprint::getCoordsetGroup(m_group, getTopologyGroup())->getName();

  sidre::Group* state_group = m_group->getGroup("state");
  if(state_group != nullptr && state_group->hasChildGroup(m_topology))
  {
    state_group = state_group->getGroup(m_topology);
    if(state_group->hasChildView("block_id"))
    {
      m_block_idx = state_group->getView("block_id")->getScalar();
    }

    if(state_group->hasChildView("partition_id"))
    {
      m_part_idx = state_group->getView("partition_id")->getScalar();
    }
  }

  SLIC_ERROR_IF(!validMeshType(), "invalid mesh type=" << m_type);
  SLIC_ERROR_IF(!validDimension(), "invalid mesh dimension=" << m_ndims);

  allocateFieldData();
}

//------------------------------------------------------------------------------
Mesh::Mesh(int ndims, int type, sidre::Group* group, const std::string& topo, const std::string& coordset)
  : m_ndims(ndims)
  , m_type(type)
  , m_block_idx(-1)
  , m_part_idx(-1)
  , m_explicit_coords(false)
  , m_explicit_connectivity(false)
  , m_has_mixed_topology(false)
  , m_group(group)
  , m_topology(topo)
  , m_coordset()
{
  SLIC_ERROR_IF(!validMeshType(), "invalid mesh type=" << m_type);
  SLIC_ERROR_IF(!validDimension(), "invalid mesh dimension=" << m_ndims);
  SLIC_ERROR_IF(m_group == nullptr, "NULL sidre group");

  // provide default names if not specified
  m_topology = (topo.empty()) ? "t1" : topo;
  m_coordset = (coordset.empty()) ? "c1" : coordset;

  if(!m_group->hasChildGroup("state"))
  {
    m_group->createGroup("state");
  }

  sidre::Group* state_group = m_group->getGroup("state")->createGroup(m_topology);
  state_group->createView("block_id")->setScalar(m_block_idx);
  state_group->createView("partition_id")->setScalar(m_part_idx);

  // create the coordset group
  if(!m_group->hasChildGroup("coordsets"))
  {
    m_group->createGroup("coordsets");
  }
  m_group->getGroup("coordsets")->createGroup(m_coordset);

  // create the topology group for this mesh
  if(!m_group->hasChildGroup("topologies"))
  {
    m_group->createGroup("topologies");
  }
  m_group->getGroup("topologies")->createGroup(m_topology);

  // create the fields group  for this mesh
  if(!m_group->hasChildGroup("fields"))
  {
    m_group->createGroup("fields");
  }

  allocateFieldData();
}

//------------------------------------------------------------------------------
sidre::Group* Mesh::getCoordsetGroup()
{
  const sidre::Group* g = blueprint::getCoordsetGroup(m_group, m_coordset);
  return (const_cast<sidre::Group*>(g));
}

//------------------------------------------------------------------------------
sidre::Group* Mesh::getTopologyGroup()
{
  const sidre::Group* g = blueprint::getTopologyGroup(m_group, m_topology);
  return (const_cast<sidre::Group*>(g));
}

#endif

//------------------------------------------------------------------------------
Mesh::~Mesh() { deallocateFieldData(); }

//------------------------------------------------------------------------------
void Mesh::setBlockId(int ID)
{
  m_block_idx = ID;
#ifdef AXOM_MINT_USE_SIDRE
  if(hasSidreGroup())
  {
    sidre::Group* state_group = m_group->getGroup("state")->getGroup(m_topology);
    SLIC_ASSERT(state_group != nullptr);
    sidre::View* block_view = state_group->getView("block_id");
    SLIC_ASSERT(block_view != nullptr);
    block_view->setScalar(m_block_idx);
  }
#endif
}

//------------------------------------------------------------------------------
void Mesh::setPartitionId(int ID)
{
  m_part_idx = ID;
#ifdef AXOM_MINT_USE_SIDRE
  if(hasSidreGroup())
  {
    sidre::Group* state_group = m_group->getGroup("state")->getGroup(m_topology);
    SLIC_ASSERT(state_group != nullptr);
    sidre::View* partition_view = state_group->getView("partition_id");
    SLIC_ASSERT(partition_view != nullptr);
    partition_view->setScalar(m_part_idx);
  }
#endif
}

//------------------------------------------------------------------------------
void Mesh::allocateFieldData()
{
#ifdef AXOM_MINT_USE_SIDRE
  if(hasSidreGroup())
  {
    sidre::Group* fields_group = m_group->getGroup("fields");
    SLIC_ASSERT(fields_group != nullptr);

    for(int assoc = 0; assoc < NUM_FIELD_ASSOCIATIONS; ++assoc)
    {
      m_mesh_fields[assoc] = new FieldData(assoc, fields_group, m_topology);
    }

    return;
  }
#endif

  for(int assoc = 0; assoc < NUM_FIELD_ASSOCIATIONS; ++assoc)
  {
    m_mesh_fields[assoc] = new FieldData(assoc);
  }
}

//------------------------------------------------------------------------------
void Mesh::deallocateFieldData()
{
  for(int i = 0; i < NUM_FIELD_ASSOCIATIONS; ++i)
  {
    SLIC_ASSERT(m_mesh_fields[i] != nullptr);

    delete m_mesh_fields[i];
    m_mesh_fields[i] = nullptr;
  }
}

} /* namespace mint */
} /* namespace axom */
