// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubes.cpp requires conduit"
#endif
#include "conduit_blueprint.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/detail/MarchingCubesSingleDomain.hpp"
#include "axom/quest/detail/MarchingCubesImpl.hpp"
#include "axom/fmt.hpp"

namespace axom::quest
{
const axom::StackArray<axom::IndexType, 2> twoZeros {0, 0};

MarchingCubes::MarchingCubes(RuntimePolicy runtimePolicy,
                             int allocatorID,
                             MarchingCubesDataParallelism dataParallelism)
  : m_runtimePolicy(runtimePolicy)
  , m_allocatorID(allocatorID)
  , m_dataParallelism(dataParallelism)
  , m_singles()
  , m_topologyName()
  , m_fcnFieldName()
  , m_fcnPath()
  , m_maskFieldName()
  , m_maskPath()
  , m_facetIndexOffsets(0, 0)
  , m_facetCount(0)
  , m_caseIdsFlat(0, 0, m_allocatorID)
  , m_crossingFlags(0, 0, m_allocatorID)
  , m_scannedFlags(0, 0, m_allocatorID)
  , m_facetIncrs(0, 0, m_allocatorID)
  , m_nodeCount(0)
  , m_facetNodeIds(twoZeros, m_allocatorID)
  , m_facetNodeCoords(twoZeros, m_allocatorID)
  , m_facetParentIds(0, 0, m_allocatorID)
  , m_facetDomainIds(0, 0, m_allocatorID)
{ }

// Set the object up for a blueprint mesh state.
void MarchingCubes::setMesh(const conduit::Node& bpMesh,
                            const std::string& topologyName,
                            const std::string& maskField)
{
  const conduit::Node* mdMesh = &bpMesh;
  if(bpMesh.has_path("topologies/" + topologyName))
  {
    m_singleDomainMesh.reset();
    m_singleDomainMesh.append().set_external(bpMesh);
    mdMesh = &m_singleDomainMesh;
  }
  else if(conduit::blueprint::mesh::is_multi_domain(bpMesh))
  {
    m_singleDomainMesh.reset();
  }
  else
  {
    // Neither a single domain carrying the requested topology nor a valid multi-domain mesh.
    // Error out here since wrapping it would defers the failure into an opaque fetch_existing() below.
    SLIC_ERROR(
      axom::fmt::format("MarchingCubes::setMesh: the input mesh is neither a multi-domain "
                        "Blueprint mesh nor a single domain containing topology '{}'.",
                        topologyName));
    return;
  }

  m_topologyName = topologyName;
  m_maskFieldName = maskField;

  /*
    To avoid slow memory allocations (especially on GPUs) keep the
    single-domain objects around and just re-initialize them.  Arrays
    will be cleared, but not deallocated.  The actual number of
    domains is m_domainCount, not m_singles.size().  To *really*
    deallocate memory, deallocate the MarchingCubes object.
  */
  auto newDomainCount = conduit::blueprint::mesh::number_of_domains(*mdMesh);

  if(m_singles.size() < newDomainCount)
  {
    auto tmpSize = m_singles.size();
    m_singles.resize(newDomainCount);
    for(int d = tmpSize; d < newDomainCount; ++d)
    {
      m_singles[d].reset(new detail::marching_cubes::MarchingCubesSingleDomain(*this));
    }
  }

  for(int d = 0; d < newDomainCount; ++d)
  {
    const auto& dom = mdMesh->child(d);
    m_singles[d]->setDomain(dom, m_topologyName, maskField);
  }
  for(int d = newDomainCount; d < m_singles.size(); ++d)
  {
    m_singles[d]->getImpl().clearDomain();
  }

  m_domainCount = newDomainCount;
}

void MarchingCubes::setFunctionField(const std::string& fcnField)
{
  m_fcnFieldName = fcnField;
  m_fcnPath = "fields/" + fcnField;
  for(axom::IndexType d = 0; d < m_domainCount; ++d)
  {
    m_singles[d]->setFunctionField(fcnField);
  }
}

void MarchingCubes::setUseBumpBackend(bool useBump)
{
#if !defined(AXOM_USE_BUMP)
  SLIC_ERROR_IF(useBump,
                "MarchingCubes bump backend requires Axom to be configured "
                "with the bump component.");
#endif
  m_useBumpBackend = useBump;
}

void MarchingCubes::computeIsocontour(double contourVal)
{
  AXOM_ANNOTATE_SCOPE("MarchingCubes::computeIsoContour");

  // Mark and scan domains while adding up their
  // facet counts to get the total facet counts.
  m_facetIndexOffsets.resize(m_singles.size());
  m_nodeIndexOffsets.resize(m_singles.size());
  for(axom::IndexType d = 0; d < m_domainCount; ++d)
  {
    auto& single = *m_singles[d];
    single.setContourValue(contourVal);
    single.setMaskValue(m_maskVal);
    single.setParentCellIdMode(m_parentCellIdMode);
    single.setRobustnessPolicy(m_robustnessPolicy);
    single.markCrossings();
    single.scanCrossings();
    m_facetIndexOffsets[d] = m_facetCount;
    m_nodeIndexOffsets[d] = m_nodeCount;
    m_facetCount += single.getContourCellCount();
    m_nodeCount += single.getContourNodeCount();
  }

  allocateOutputBuffers();

  // Tell singles where to put contour data.
  auto facetNodeIdsView = m_facetNodeIds.view();
  auto facetNodeCoordsView = m_facetNodeCoords.view();
  auto facetParentIdsView = m_facetParentIds.view();
  for(axom::IndexType d = 0; d < m_domainCount; ++d)
  {
    m_singles[d]->getImpl().setOutputBuffers(facetNodeIdsView,
                                             facetNodeCoordsView,
                                             facetParentIdsView,
                                             m_facetIndexOffsets[d],
                                             m_nodeIndexOffsets[d]);
  }

  for(axom::IndexType d = 0; d < m_domainCount; ++d)
  {
    m_singles[d]->computeFacets();
  }

  for(axom::IndexType d = 0; d < m_domainCount; ++d)
  {
    const auto domainId = m_singles[d]->getDomainId(d);
    const auto domainFacetCount =
      (d < m_domainCount - 1 ? m_facetIndexOffsets[d + 1] : m_facetCount) - m_facetIndexOffsets[d];
    m_facetDomainIds.fill(domainId, domainFacetCount, m_facetIndexOffsets[d]);
  }
}

axom::IndexType MarchingCubes::getContourNodeCount() const { return m_nodeCount; }

void MarchingCubes::clearOutput()
{
  m_facetCount = 0;
  m_nodeCount = 0;
  m_facetNodeIds.clear();
  m_facetNodeCoords.clear();
  m_facetParentIds.clear();
  m_facetDomainIds.clear();
}

void MarchingCubes::populateContourMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
                                        const std::string& cellIdField,
                                        const std::string& domainIdField) const
{
  AXOM_ANNOTATE_SCOPE("MarchingCubes::populateContourMesh");
  if(!cellIdField.empty() && !mesh.hasField(cellIdField, axom::mint::CELL_CENTERED))
  {
    // Create cellId field, currently the multidimensional index of the parent cell.
    mesh.createField<axom::IndexType>(cellIdField, axom::mint::CELL_CENTERED);
  }

  if(!domainIdField.empty() && !mesh.hasField(domainIdField, axom::mint::CELL_CENTERED))
  {
    mesh.createField<DomainIdType>(domainIdField, axom::mint::CELL_CENTERED);
  }

  // Reserve space once for all local domains.
  const axom::IndexType contourCellCount = getContourCellCount();
  const axom::IndexType contourNodeCount = getContourNodeCount();
  mesh.reserveCells(contourCellCount);
  mesh.reserveNodes(contourNodeCount);

  if(m_facetCount)
  {
    // Put nodes and cells into the mesh.
    // If data is not in host memory, copy to temporary host memory first.
    axom::MemorySpace internalMemorySpace = axom::detail::getAllocatorSpace(m_allocatorID);
    const bool hostAndInternalMemoriesAreSeparate = internalMemorySpace != axom::MemorySpace::Dynamic
#ifdef AXOM_USE_UMPIRE
      && internalMemorySpace != axom::MemorySpace::Host
#endif
      ;
    const int hostAllocatorId = hostAndInternalMemoriesAreSeparate
      ?
#ifdef AXOM_USE_UMPIRE
      axom::detail::getAllocatorID<axom::MemorySpace::Host>()
#else
      axom::detail::getAllocatorID<axom::MemorySpace::Dynamic>()
#endif
      : m_allocatorID;

    if(hostAndInternalMemoriesAreSeparate)
    {
      axom::Array<double, 2> tmpfacetNodeCoords(m_facetNodeCoords, hostAllocatorId);
      axom::Array<axom::IndexType, 2> tmpfacetNodeIds(m_facetNodeIds, hostAllocatorId);
      mesh.appendNodes(tmpfacetNodeCoords.data(), contourNodeCount);
      mesh.appendCells(tmpfacetNodeIds.data(), contourCellCount);
    }
    else
    {
      mesh.appendNodes(m_facetNodeCoords.data(), contourNodeCount);
      mesh.appendCells(m_facetNodeIds.data(), contourCellCount);
    }

    if(!cellIdField.empty())
    {
      // Put parent cell ids into the mesh.
      axom::IndexType* cellIdPtr =
        mesh.getFieldPtr<axom::IndexType>(cellIdField, axom::mint::CELL_CENTERED);
      axom::copy(cellIdPtr, m_facetParentIds.data(), m_facetCount * sizeof(axom::IndexType));
    }

    if(!domainIdField.empty())
    {
      // Put parent domain ids into the mesh.
      auto* domainIdPtr = mesh.getFieldPtr<DomainIdType>(domainIdField, axom::mint::CELL_CENTERED);
      axom::copy(domainIdPtr, m_facetDomainIds.data(), m_facetCount * sizeof(axom::IndexType));
    }
  }
}

void MarchingCubes::populateContourMeshBlueprint(conduit::Node& bpMesh, bool triangulate) const
{
  AXOM_ANNOTATE_SCOPE("MarchingCubes::populateContourMeshBlueprint");
  bpMesh.reset();

  SLIC_ERROR_IF(!m_useBumpBackend,
                "MarchingCubes Blueprint contour output is available only when "
                "setUseBumpBackend(true) was used.");

  for(axom::IndexType d = 0; d < m_domainCount; ++d)
  {
    const auto& single = *m_singles[d];
    const auto& impl = single.getImpl();
    SLIC_ERROR_IF(!impl.hasContourMeshBlueprint(),
                  "MarchingCubes has no Blueprint contour output. "
                  "Call computeIsocontour() before requesting it.");

    conduit::Node& outDom = bpMesh.append();
    impl.copyContourMeshBlueprint(outDom, triangulate);
    if(!outDom.has_path("state/domain_id"))
    {
      outDom["state/domain_id"] = single.getDomainId(static_cast<int32_t>(d));
    }
  }
}

void MarchingCubes::relinquishContourDataBlueprint(conduit::Node& bpMesh)
{
  AXOM_ANNOTATE_SCOPE("MarchingCubes::relinquishContourDataBlueprint");
  bpMesh.reset();

  SLIC_ERROR_IF(!m_useBumpBackend,
                "MarchingCubes Blueprint contour output is available only when "
                "setUseBumpBackend(true) was used.");

  for(axom::IndexType d = 0; d < m_domainCount; ++d)
  {
    auto& single = *m_singles[d];
    auto& impl = single.getImpl();
    SLIC_ERROR_IF(!impl.hasContourMeshBlueprint(),
                  "MarchingCubes has no Blueprint contour output. "
                  "Call computeIsocontour() before requesting it.");

    conduit::Node& outDom = bpMesh.append();
    impl.relinquishContourMeshBlueprint(outDom);
    if(!outDom.has_path("state/domain_id"))
    {
      outDom["state/domain_id"] = single.getDomainId(static_cast<int32_t>(d));
    }
  }

  clearOutput();
}

void MarchingCubes::allocateOutputBuffers()
{
  AXOM_ANNOTATE_SCOPE("MarchingCubes::allocateOutputBuffers");
  if(!m_singles.empty())
  {
    int ndim = m_singles[0]->spatialDimension();
    m_facetNodeIds.resize(axom::StackArray<axom::IndexType, 2> {m_facetCount, ndim}, 0);
    m_facetNodeCoords.resize(axom::StackArray<axom::IndexType, 2> {m_nodeCount, ndim}, 0.0);
    m_facetParentIds.resize(axom::StackArray<axom::IndexType, 1> {m_facetCount}, 0);
    m_facetDomainIds.resize(axom::StackArray<axom::IndexType, 1> {m_facetCount}, 0);
  }
}

}  // end namespace axom::quest
