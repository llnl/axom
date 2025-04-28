// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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

namespace axom
{
namespace quest
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
  SLIC_ASSERT_MSG(conduit::blueprint::mesh::is_multi_domain(bpMesh),
                  "MarchingCubes class input mesh must be in multidomain format.");

  m_topologyName = topologyName;
  m_maskFieldName = maskField;

  /*
    To avoid slow memory allocations (especially on GPUs) keep the
    single-domain objects around and just re-initialize them.  Arrays
    will be cleared, but not deallocated.  The actual number of
    domains is m_domainCount, not m_singles.size().  To *really*
    deallocate memory, deallocate the MarchingCubes object.
  */
  auto newDomainCount = conduit::blueprint::mesh::number_of_domains(bpMesh);

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
    const auto& dom = bpMesh.child(d);
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

void MarchingCubes::computeIsocontour(double contourVal)
{
  AXOM_ANNOTATE_SCOPE("MarchingCubes::computeIsoContour");

  // Mark and scan domains while adding up their
  // facet counts to get the total facet counts.
  m_facetIndexOffsets.resize(m_singles.size());
  for(axom::IndexType d = 0; d < m_domainCount; ++d)
  {
    auto& single = *m_singles[d];
    single.setContourValue(contourVal);
    single.setMaskValue(m_maskVal);
    single.markCrossings();
    single.scanCrossings();
    m_facetIndexOffsets[d] = m_facetCount;
    m_facetCount += single.getContourCellCount();
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
                                             m_facetIndexOffsets[d]);
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

axom::IndexType MarchingCubes::getContourNodeCount() const
{
  axom::IndexType contourNodeCount =
    (m_domainCount > 0) ? m_facetCount * m_singles[0]->spatialDimension() : 0;
  return contourNodeCount;
}

void MarchingCubes::clearOutput()
{
  m_facetCount = 0;
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

void MarchingCubes::allocateOutputBuffers()
{
  AXOM_ANNOTATE_SCOPE("MarchingCubes::allocateOutputBuffers");
  if(!m_singles.empty())
  {
    int ndim = m_singles[0]->spatialDimension();
    const auto nodeCount = m_facetCount * ndim;
    m_facetNodeIds.resize(axom::StackArray<axom::IndexType, 2> {m_facetCount, ndim}, 0);
    m_facetNodeCoords.resize(axom::StackArray<axom::IndexType, 2> {nodeCount, ndim}, 0.0);
    m_facetParentIds.resize(axom::StackArray<axom::IndexType, 1> {m_facetCount}, 0);
    m_facetDomainIds.resize(axom::StackArray<axom::IndexType, 1> {m_facetCount}, 0);
  }
}

}  // end namespace quest
}  // end namespace axom
