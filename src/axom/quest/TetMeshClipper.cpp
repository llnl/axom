// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT
  #include "conduit_blueprint.hpp"
#endif

#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/TetMeshClipper.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{

TetMeshClipper::TetMeshClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("TetMesh") : name)
  , m_topoName(kGeom.getBlueprintTopology())
  , m_cellCount(0)
{
  SLIC_ASSERT(!m_topoName.empty());

  extractClipperInfo();
}

bool TetMeshClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  return false;  // TODO: implement TetMeshClipper::labelInOut
}

bool TetMeshClipper::getShapeAsTets(quest::ShapeeMesh& shapeeMesh, axom::Array<TetrahedronType>& tets)
{
  const int hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int allocId = shapeeMesh.getAllocatorId();

  if(tets.getAllocatorID() != allocId || tets.size() != m_cellCount)
  {
    tets = axom::Array<TetrahedronType>(m_cellCount, m_cellCount, allocId);
  }

  /*
    1. Initialize a mint Mesh intermediary from the blueprint mesh.
       mint::getMesh() utility for this makes saves some coding.
       Note that blueprint mesh is on host.
    2. Populate a tet array on the host (because mint data works only for host).
       Host array could be tets or a temporary array.
    3. If needed, copy temporary array to tets.
  */

  axom::Array<TetrahedronType> tmpTets(0, 0, hostAllocId);

  // Host loop writes directly to tets if possible, else a temporary array.
  axom::Array<TetrahedronType>& tetsOnHost =
    axom::execution_space<axom::SEQ_EXEC>::usesAllocId(tets.getAllocatorID()) ? tets : tmpTets;

  if(&tetsOnHost == &tmpTets)
  {
    tmpTets.resize(m_cellCount, TetrahedronType());
  }

  axom::sidre::DataStore ds;
  auto* bpMeshGrp = ds.getRoot()->createGroup("blueprintMesh");
  bpMeshGrp->importConduitTree(*m_bpMesh);
  std::shared_ptr<mint::UnstructuredMesh<axom::mint::Topology::SINGLE_SHAPE>> mintMesh {
    (mint::UnstructuredMesh<axom::mint::Topology::SINGLE_SHAPE>*)axom::mint::getMesh(bpMeshGrp,
                                                                                     m_topoName)};

  // Initialize tetrahedra
  axom::Array<IndexType> nodeIds(4);
  axom::Array<Point3D> pts(4);
  for(int i = 0; i < m_cellCount; i++)
  {
    mintMesh->getCellNodeIDs(i, nodeIds.data());

    mintMesh->getNode(nodeIds[0], pts[0].data());
    mintMesh->getNode(nodeIds[1], pts[1].data());
    mintMesh->getNode(nodeIds[2], pts[2].data());
    mintMesh->getNode(nodeIds[3], pts[3].data());

    tetsOnHost[i] = TetrahedronType({pts[0], pts[1], pts[2], pts[3]});
  }

  if(&tets != &tetsOnHost)
  {
    axom::copy(tets.data(), tetsOnHost.data(), sizeof(TetrahedronType) * tets.size());
  }

  return true;
}

void TetMeshClipper::extractClipperInfo()
{
  m_topoName = m_info.fetch_existing("topologyName").to_string();

  m_bpMesh = &m_info.fetch_existing("unstructuredBlueprintTetMesh");

  SLIC_ASSERT(
    m_bpMesh->fetch_existing("topologies").fetch_existing(m_topoName).fetch_existing("type").as_string() ==
    "unstructured");

  conduit::Node& topoNode = m_bpMesh->fetch_existing("topologies").fetch_existing(m_topoName);

  bool isMultiDomain = conduit::blueprint::mesh::is_multi_domain(*m_bpMesh);
  SLIC_ERROR_IF(isMultiDomain, "TetMeshClipper does not support multi-domain tet meshes yet.");

  SLIC_ASSERT(conduit::blueprint::mesh::topology::dims(topoNode) == 3);

  m_cellCount = conduit::blueprint::mesh::topology::length(topoNode);
}

}  // end namespace quest
}  // end namespace axom
