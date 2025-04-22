// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/mesh/UniformMesh.hpp"

// mint includes
#include "axom/mint/mesh/blueprint.hpp"  // for blueprint functions
#include "axom/mint/config.hpp"          // for compile-time definitions
#include "axom/mint/mesh/MeshTypes.hpp"  // for STRUCTURED_UNIFORM_MESH

#include "axom/mint/mesh/internal/MeshHelpers.hpp"  // for internal helper

// core includes
#include "axom/core/utilities/Utilities.hpp"  // for utilities::isNearlyEqual

namespace axom
{
namespace mint
{
//------------------------------------------------------------------------------
// UNIFORM MESH IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
UniformMesh::UniformMesh(const double* lower_bound,
                         const double* upper_bound,
                         IndexType Ni,
                         IndexType Nj,
                         IndexType Nk)
  : StructuredMesh(STRUCTURED_UNIFORM_MESH, Ni, Nj, Nk)
{
  SLIC_ERROR_IF(lower_bound == nullptr, "supplied null for lower_bound");
  SLIC_ERROR_IF(upper_bound == nullptr, "supplied null for upper_bound");

  setSpacingAndOrigin(lower_bound, upper_bound);
}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
UniformMesh::UniformMesh(sidre::Group* group, const std::string& topo) : StructuredMesh(group, topo)
{
  SLIC_ERROR_IF(m_type != STRUCTURED_UNIFORM_MESH,
                "supplied Sidre group does not correspond to a UniformMesh!");

  blueprint::getUniformMeshProperties(m_ndims, m_origin, m_h, getCoordsetGroup());
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh(sidre::Group* group,
                         const std::string& topo,
                         const std::string& coordset,
                         const double* lower_bound,
                         const double* upper_bound,
                         IndexType Ni,
                         IndexType Nj,
                         IndexType Nk)
  : StructuredMesh(STRUCTURED_UNIFORM_MESH, Ni, Nj, Nk, group, topo, coordset)
{
  SLIC_ERROR_IF(lower_bound == nullptr, "supplied null for lower_bound");
  SLIC_ERROR_IF(upper_bound == nullptr, "supplied null for upper_bound");

  // STEP 0: initialize mesh
  setSpacingAndOrigin(lower_bound, upper_bound);

  // STEP 1: populate sidre
  blueprint::setUniformMeshProperties(m_ndims, m_origin, m_h, getCoordsetGroup());
}

#endif

//------------------------------------------------------------------------------
void UniformMesh::getNode(IndexType nodeID, double* node) const
{
  SLIC_ASSERT(0 <= nodeID && nodeID < getNumberOfNodes());
  SLIC_ASSERT(node != nullptr);

  IndexType i = -1;
  IndexType j = -1;
  IndexType k = -1;
  switch(m_ndims)
  {
  case 1:
    node[0] = evaluateCoordinate(nodeID, I_DIRECTION);
    break;
  case 2:
    getNodeGridIndex(nodeID, i, j);
    node[0] = evaluateCoordinate(i, I_DIRECTION);
    node[1] = evaluateCoordinate(j, J_DIRECTION);
    break;
  default:
    SLIC_ASSERT(m_ndims == 3);
    getNodeGridIndex(nodeID, i, j, k);
    node[0] = evaluateCoordinate(i, I_DIRECTION);
    node[1] = evaluateCoordinate(j, J_DIRECTION);
    node[2] = evaluateCoordinate(k, K_DIRECTION);
  }  // END switch()
}

void UniformMesh::setSpacingAndOrigin(const double* lo, const double* hi)
{
  SLIC_ASSERT(lo != nullptr);
  SLIC_ASSERT(hi != nullptr);

  for(int dim = 0; dim < m_ndims; ++dim)
  {
    m_origin[dim] = lo[dim];
    double dx = hi[dim] - lo[dim];
    SLIC_ERROR_IF(utilities::isNearlyEqual(dx, 0.0) || dx < 0.0, "supplied invalid bounds!");
    m_h[dim] = dx / getCellResolution(dim);
  }
}

} /* namespace mint */
} /* namespace axom */
