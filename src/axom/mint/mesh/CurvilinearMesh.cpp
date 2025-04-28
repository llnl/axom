// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/mesh/CurvilinearMesh.hpp"

// mint includes
#include "axom/mint/mesh/blueprint.hpp"        // for blueprint functions
#include "axom/mint/mesh/MeshTypes.hpp"        // STRUCTURED_CURVILINEAR_MESH
#include "axom/mint/mesh/MeshCoordinates.hpp"  // for MeshCoordinates class

#include "axom/mint/mesh/internal/MeshHelpers.hpp"  // for internal helpers

// slic includes
#include "axom/slic/interface/slic.hpp"  // for SLIC macros

namespace axom
{
namespace mint
{
//------------------------------------------------------------------------------
// CURVILINEAR MESH IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh(IndexType Ni, IndexType Nj, IndexType Nk)
  : StructuredMesh(STRUCTURED_CURVILINEAR_MESH, Ni, Nj, Nk)
  , m_coordinates(new mint::MeshCoordinates(m_ndims, getNumberOfNodes()))
{
  initialize();

  // sanity checks
  SLIC_ASSERT(m_coordinates != nullptr);
  SLIC_ASSERT(getNumberOfNodes() == m_coordinates->numNodes());
  SLIC_ASSERT(m_coordinates->dimension() == m_ndims);
}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh(IndexType Ni, double* x, IndexType Nj, double* y, IndexType Nk, double* z)
  : StructuredMesh(STRUCTURED_CURVILINEAR_MESH, Ni, Nj, Nk)
  , m_coordinates(new mint::MeshCoordinates(getNumberOfNodes(), x, y, z))
{
  initialize();

  // sanity checks
  SLIC_ASSERT(m_coordinates != nullptr);
  SLIC_ASSERT(getNumberOfNodes() == m_coordinates->numNodes());
  SLIC_ASSERT(m_coordinates->dimension() == m_ndims);
}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh(sidre::Group* group, const std::string& topo)
  : StructuredMesh(group, topo)
  , m_coordinates(new MeshCoordinates(getCoordsetGroup()))
{
  SLIC_ERROR_IF(m_type != STRUCTURED_CURVILINEAR_MESH,
                "supplied Sidre group does not correspond to a CurvilinearMesh");

  initialize();

  // sanity checks
  SLIC_ASSERT(m_coordinates != nullptr);
  SLIC_ASSERT(getNumberOfNodes() == m_coordinates->numNodes());
  SLIC_ASSERT(m_coordinates->dimension() == m_ndims);
}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh(sidre::Group* group,
                                 const std::string& topo,
                                 const std::string& coordset,
                                 IndexType Ni,
                                 IndexType Nj,
                                 IndexType Nk)
  : StructuredMesh(STRUCTURED_CURVILINEAR_MESH, Ni, Nj, Nk, group, topo, coordset)
{
  m_coordinates =
    new mint::MeshCoordinates(getCoordsetGroup(), m_ndims, getNumberOfNodes(), getNumberOfNodes());

  initialize();

  // sanity checks
  SLIC_ASSERT(m_coordinates != nullptr);
  SLIC_ASSERT(getNumberOfNodes() == m_coordinates->numNodes());
  SLIC_ASSERT(m_coordinates->dimension() == m_ndims);
}

#endif

//------------------------------------------------------------------------------
CurvilinearMesh::~CurvilinearMesh()
{
  delete m_coordinates;
  m_coordinates = nullptr;
}

//------------------------------------------------------------------------------
void CurvilinearMesh::initialize()
{
  m_explicit_coords = true;
  m_explicit_connectivity = false;
  m_has_mixed_topology = false;
}

} /* namespace mint */
} /* namespace axom */
