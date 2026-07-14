// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubesImpl.hpp requires conduit"
#endif
#include "conduit_blueprint.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/core/MDMapping.hpp"
#include "axom/quest/MeshViewUtil.hpp"
#include "axom/quest/detail/MarchingCubesSingleDomain.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/constants.hpp"
#include "axom/core/execution/nested_for_exec.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
namespace detail
{
namespace marching_cubes
{
/*!
  @brief Computations for MarchingCubesSingleDomain

  Spatial dimension and execution space are here as template
  parameters, to keep out of higher level classes MarchingCubes and
  MarchingCubesSingleDomain.

  ExecSpace is the general execution space, like axom::SEQ_EXEC and
  axom::CUDA_EXEC<256>.
*/
template <int DIM, typename ExecSpace, typename SequentialExecSpace>
class MarchingCubesImpl : public MarchingCubesSingleDomain::ImplBase
{
public:
  using Point = axom::primal::Point<double, DIM>;
  using MIdx = axom::StackArray<axom::IndexType, DIM>;
  using MDMapper = axom::MDMapping<DIM>;
  using FacetIdType = int;
  using CrossingFlagType = axom::quest::MarchingCubes::CrossingFlagType;

  static constexpr auto MemorySpace = execution_space<ExecSpace>::memory_space;

  AXOM_HOST MarchingCubesImpl(int allocatorID,
                              axom::Array<std::uint16_t>& caseIdsFlat,
                              axom::Array<CrossingFlagType>& crossingFlags,
                              axom::Array<axom::IndexType>& scannedFlags,
                              axom::Array<axom::IndexType>& facetIncrs)
    : m_allocatorID(allocatorID)
    , m_caseIds()
    , m_caseIdsMDMapper()
    , m_caseIdsFlat(caseIdsFlat)
    , m_crossingFlags(crossingFlags)
    , m_scannedFlags(scannedFlags)
    , m_facetIncrs(facetIncrs)
    , m_crossingCases(0, 0, m_allocatorID)
    , m_crossingParentIds(0, 0, m_allocatorID)
    , m_firstFacetIds(0, 0, m_allocatorID)
  {
    SLIC_ASSERT(caseIdsFlat.getAllocatorID() == allocatorID);
    SLIC_ASSERT(crossingFlags.getAllocatorID() == allocatorID);
    SLIC_ASSERT(scannedFlags.getAllocatorID() == allocatorID);
    SLIC_ASSERT(facetIncrs.getAllocatorID() == allocatorID);
  }

  /*!
    @brief Initialize data to a blueprint domain.
    @param dom Blueprint structured mesh domain
    @param topologyName Name of mesh topology (see blueprint
           mesh documentation)
    @param maskFieldName Name of integer cell mask function is in dom

    Set up views to domain data and allocate other data to work on the
    given domain.

    The above data from the domain MUST be in a memory space
    compatible with ExecSpace.
  */
  AXOM_HOST void setDomain(const conduit::Node& dom,
                           const std::string& topologyName,
                           const std::string& maskFieldName) override
  {
    // Time this due to potentially slow memory allocation
    AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::initialize");
    clearDomain();

    SLIC_ASSERT(conduit::blueprint::mesh::topology::dims(
                  dom.fetch_existing(axom::fmt::format("topologies/{}", topologyName))) == DIM);

    m_mvu = axom::quest::MeshViewUtil<DIM, MemorySpace>(dom, topologyName);

    m_bShape = m_mvu.getCellShape();
    m_coordsViews = m_mvu.getConstCoordsViews(false);
    if(!maskFieldName.empty())
    {
      m_maskView = m_mvu.template getConstFieldView<int>(maskFieldName, false);
    }
  }

  AXOM_HOST void setDataParallelism(MarchingCubesDataParallelism dataPar) override
  {
    constexpr MarchingCubesDataParallelism autoPolicy =
      std::is_same<ExecSpace, axom::SEQ_EXEC>::value ? MarchingCubesDataParallelism::hybridParallel
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
      : std::is_same<ExecSpace, axom::OMP_EXEC>::value ? MarchingCubesDataParallelism::hybridParallel
#endif
                                                       : MarchingCubesDataParallelism::fullParallel;

    m_dataParallelism = dataPar;

    if(m_dataParallelism == axom::quest::MarchingCubesDataParallelism::byPolicy)
    {
      m_dataParallelism = autoPolicy;
    }
  }

  /*!
    @brief Set the scale field name
    @param fcnFieldName Name of nodal function is in dom
  */
  void setFunctionField(const std::string& fcnFieldName) override
  {
    m_fcnView = m_mvu.template getConstFieldView<double>(fcnFieldName, false);
  }

  void setContourValue(double contourVal) override { m_contourVal = contourVal; }

  void setMaskValue(int maskVal) override { m_maskVal = maskVal; }

  /*!
    @brief Implementation of virtual markCrossings.

    Virtual methods cannot be templated, so this implementation
    delegates to an implementation templated on DIM.
  */
  void markCrossings() override
  {
    AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::markCrossings");

    m_caseIdsFlat.resize(m_mvu.getCellCount(), 0);
    m_caseIdsFlat.fill(0);

    // Choose caseIds stride order to match function stride order.
    MDMapper fcnMDMapper(m_fcnView.strides());
    m_caseIdsMDMapper.initializeShape(m_bShape, fcnMDMapper.slowestDirs());
    m_caseIds = axom::ArrayView<std::uint16_t, DIM, MemorySpace>(m_caseIdsFlat.data(),
                                                                 m_bShape,
                                                                 m_caseIdsMDMapper.strides());
    SLIC_ASSERT_MSG(MDMapper(m_caseIds.strides()).getStrideOrder() == fcnMDMapper.getStrideOrder(),
                    "Mismatched order is inefficient.");

    markCrossings_dim();
  }

  //!@brief Populate m_caseIds with crossing indices.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type markCrossings_dim()
  {
    MarkCrossings_Util mcu(m_caseIds, m_fcnView, m_maskView, m_contourVal, m_maskVal);

    const auto order = m_caseIdsMDMapper.getStrideOrder();
    if(int(order) & int(axom::ArrayStrideOrder::COLUMN))
    {
      axom::for_all<ExecSpace>(
        m_bShape,
        AXOM_LAMBDA(axom::IndexType i, axom::IndexType j) { mcu.computeCaseId(i, j); });
    }
    else
    {
      axom::StackArray<axom::IndexType, 2> shapeJI {{m_bShape[1], m_bShape[0]}};
      axom::for_all<ExecSpace>(
        shapeJI,
        AXOM_LAMBDA(axom::IndexType j, axom::IndexType i) { mcu.computeCaseId(i, j); });
    }
  }

  //!@brief Populate m_caseIds with crossing indices.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type markCrossings_dim()
  {
    MarkCrossings_Util mcu(m_caseIds, m_fcnView, m_maskView, m_contourVal, m_maskVal);

    auto order = m_caseIdsMDMapper.getStrideOrder();
    // order ^= axom::ArrayStrideOrder::BOTH; // Pick wrong ordering to test behavior.
    if(int(order) & int(axom::ArrayStrideOrder::COLUMN))
    {
      axom::for_all<ExecSpace>(
        m_bShape,
        AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k) {
          mcu.computeCaseId(i, j, k);
        });
    }
    else
    {
      axom::StackArray<axom::IndexType, 3> shapeKJI {{m_bShape[2], m_bShape[1], m_bShape[0]}};
      axom::for_all<ExecSpace>(
        shapeKJI,
        AXOM_LAMBDA(axom::IndexType k, axom::IndexType j, axom::IndexType i) {
          mcu.computeCaseId(i, j, k);
        });
    }
  }

  /*!
    @brief Implementation used by MarchingCubesImpl::markCrossings_dim()
    containing just the objects needed for that part, to be made available
    on devices.
  */
  struct MarkCrossings_Util
  {
    axom::ArrayView<std::uint16_t, DIM, MemorySpace> caseIdsView;
    axom::ArrayView<const double, DIM, MemorySpace> fcnView;
    axom::ArrayView<const int, DIM, MemorySpace> maskView;
    double contourVal;
    int maskVal;
    MarkCrossings_Util(axom::ArrayView<std::uint16_t, DIM, MemorySpace>& caseIds,
                       axom::ArrayView<const double, DIM, MemorySpace>& fcnView_,
                       axom::ArrayView<const int, DIM, MemorySpace>& maskView_,
                       double contourVal_,
                       int maskVal_)
      : caseIdsView(caseIds)
      , fcnView(fcnView_)
      , maskView(maskView_)
      , contourVal(contourVal_)
      , maskVal(maskVal_)
    { }

    //!@brief Compute the case index into cases2D or cases3D.
    AXOM_HOST_DEVICE inline int computeCrossingCase(const double* f) const
    {
      int index = 0;
      for(int n = 0; n < CELL_CORNER_COUNT; ++n)
      {
        if(f[n] >= contourVal)
        {
          const int bit = (1 << n);
          index |= bit;
        }
      }
      return index;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 2>::type computeCaseId(axom::IndexType i,
                                                                                   axom::IndexType j) const
    {
      const bool useZone = maskView.empty() || (maskView(i, j) == maskVal);
      if(useZone)
      {
        // clang-format off
        double nodalValues[CELL_CORNER_COUNT] =
          {fcnView(i    , j    ),
           fcnView(i + 1, j    ),
           fcnView(i + 1, j + 1),
           fcnView(i    , j + 1)};
        // clang-format on
        caseIdsView(i, j) = computeCrossingCase(nodalValues);
      }
    }

    //!@brief Populate m_caseIds with crossing indices.
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 3>::type computeCaseId(axom::IndexType i,
                                                                                   axom::IndexType j,
                                                                                   axom::IndexType k) const
    {
      const bool useZone = maskView.empty() || (maskView(i, j, k) == maskVal);
      if(useZone)
      {
        // clang-format off
        double nodalValues[CELL_CORNER_COUNT] =
          {fcnView(i + 1, j    , k    ),
           fcnView(i + 1, j + 1, k    ),
           fcnView(i    , j + 1, k    ),
           fcnView(i    , j    , k    ),
           fcnView(i + 1, j    , k + 1),
           fcnView(i + 1, j + 1, k + 1),
           fcnView(i    , j + 1, k + 1),
           fcnView(i    , j    , k + 1)};
        // clang-format on
        caseIdsView(i, j, k) = computeCrossingCase(nodalValues);
      }
    }
  };  // MarkCrossings_Util

  void scanCrossings() override
  {
    AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::scanCrossings");
    if(m_dataParallelism == axom::quest::MarchingCubesDataParallelism::hybridParallel)
    {
      AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::scanCrossings:hybridParallel");
      scanCrossings_hybridParallel();
    }
    else if(m_dataParallelism == axom::quest::MarchingCubesDataParallelism::fullParallel)
    {
      AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::scanCrossings:fullParallel");
      scanCrossings_fullParallel();
    }
  }

  void allocateIndexLists()
  {
    AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::allocateIndexLists");
    m_crossingParentIds.resize(m_crossingCount, 0);
    m_crossingCases.resize(m_crossingCount, 0);
    m_facetIncrs.resize(m_crossingCount, 0);
    m_firstFacetIds.resize(1 + m_crossingCount, 0);
  }

  void scanCrossings_fullParallel()
  {
    const axom::IndexType parentCellCount = m_caseIds.size();
    auto caseIdsView = m_caseIds;

    //
    // Initialize m_crossingFlags, prefix-sum it, and count the
    // crossings.
    //

    m_crossingFlags.resize(m_mvu.getCellCount(), 0);
    m_scannedFlags.resize(1 + m_mvu.getCellCount(), 0);

    auto crossingFlagsView = m_crossingFlags.view();
    // Facet counts can change once ambiguous and iso-degenerate cells are
    // resolved, so the scan must use the same logic as facet generation.
    ComputeFacets_Util cfu(m_contourVal, m_caseIdsMDMapper, m_fcnView, m_coordsViews);
    {
      AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::scanCrossings:set_flags");
      axom::for_all<ExecSpace>(
        0,
        parentCellCount,
        AXOM_LAMBDA(axom::IndexType parentCellId) {
          auto numContourCells =
            cfu.contour_cell_count(caseIdsView.flatIndex(parentCellId), parentCellId);
          crossingFlagsView[parentCellId] = bool(numContourCells);
        });
    }

    m_scannedFlags.fill(0, 1, 0);

    {
      AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::scanCrossings:scan_flags");
      axom::inclusive_scan<ExecSpace>(
        axom::ArrayView<CrossingFlagType>(m_crossingFlags.data(), parentCellCount),
        axom::ArrayView<axom::IndexType>(m_scannedFlags.data() + 1, parentCellCount));
    }

    axom::copy(&m_crossingCount,
               m_scannedFlags.data() + m_scannedFlags.size() - 1,
               sizeof(axom::IndexType));

    //
    // Generate crossing info in compact arrays.
    //
    allocateIndexLists();
    auto scannedFlagsView = m_scannedFlags.view();
    auto crossingParentIdsView = m_crossingParentIds.view();
    auto crossingCasesView = m_crossingCases.view();
    auto facetIncrsView = m_facetIncrs.view();

    {
      AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::scanCrossings:set_incrs");
      auto loopBody = AXOM_LAMBDA(axom::IndexType parentCellId)
      {
        if(scannedFlagsView[parentCellId] != scannedFlagsView[1 + parentCellId])
        {
          auto crossingId = scannedFlagsView[parentCellId];
          auto caseId = caseIdsView.flatIndex(parentCellId);
          auto facetIncr = cfu.contour_cell_count(caseId, parentCellId);
          crossingParentIdsView[crossingId] = parentCellId;
          crossingCasesView[crossingId] = caseId;
          facetIncrsView[crossingId] = facetIncr;
        }
      };
      axom::for_all<ExecSpace>(0, parentCellCount, loopBody);
    }

    //
    // Prefix-sum the facets counts to get first facet id for each crossing
    // and the total number of facets.
    //

    m_firstFacetIds.fill(0, 1, 0);

    {
      AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::scanCrossings:scan_incrs");
      axom::inclusive_scan<ExecSpace>(
        axom::ArrayView<FacetIncrsType>(m_facetIncrs.data(), m_crossingCount),
        axom::ArrayView<axom::IndexType>(m_firstFacetIds.data() + 1, m_crossingCount));
    }

    axom::copy(&m_facetCount,
               m_firstFacetIds.data() + m_firstFacetIds.size() - 1,
               sizeof(axom::IndexType));
    // m_firstFacetIds.resize(m_crossingCount);
  }

  void scanCrossings_hybridParallel()
  {
    //
    // Compute number of crossings in m_caseIds
    //
    const axom::IndexType parentCellCount = m_caseIds.size();
    auto caseIdsView = m_caseIds;
    ComputeFacets_Util cfu(m_contourVal, m_caseIdsMDMapper, m_fcnView, m_coordsViews);
    axom::ReduceSum<ExecSpace, axom::IndexType> vsum(0);
    axom::for_all<ExecSpace>(
      parentCellCount,
      AXOM_LAMBDA(axom::IndexType n) {
        vsum += bool(cfu.contour_cell_count(caseIdsView.flatIndex(n), n));
      });
    m_crossingCount = static_cast<axom::IndexType>(vsum.get());

    //
    // Allocate space for crossing info
    //
    allocateIndexLists();
    auto crossingParentIdsView = m_crossingParentIds.view();
    auto crossingCasesView = m_crossingCases.view();
    auto facetIncrsView = m_facetIncrs.view();

    axom::IndexType* crossingId =
      axom::allocate<axom::IndexType>(1, axom::detail::getAllocatorID<MemorySpace>());

    auto loopBody = AXOM_LAMBDA(axom::IndexType n)
    {
      auto caseId = caseIdsView.flatIndex(n);
      auto ccc = cfu.contour_cell_count(caseId, n);
      if(ccc != 0)
      {
        crossingParentIdsView[*crossingId] = n;
        crossingCasesView[*crossingId] = caseId;
        facetIncrsView[*crossingId] = ccc;
        ++(*crossingId);
      }
    };

    /*
      loopBody isn't data-parallel and shouldn't be parallelized.
      This contrived for_all forces it to run sequentially in the memory
      space associated with the implementation.
    */
    axom::for_all<SequentialExecSpace>(1, [=] AXOM_HOST_DEVICE(axom::IndexType /* i */) {
      *crossingId = 0;
      for(axom::IndexType n = 0; n < parentCellCount; ++n)
      {
        loopBody(n);
      }
    });
    axom::IndexType finalCrossingId = 0;
    axom::copy(&finalCrossingId, crossingId, sizeof(axom::IndexType));
    SLIC_ASSERT(finalCrossingId == m_crossingCount);

    axom::deallocate(crossingId);

    m_firstFacetIds.fill(0, 1, 0);

    const auto firstFacetIdsView = m_firstFacetIds.view();
    axom::inclusive_scan<ExecSpace>(
      axom::ArrayView<FacetIncrsType>(facetIncrsView.data(), m_crossingCount),
      axom::ArrayView<axom::IndexType>(firstFacetIdsView.data() + 1, m_crossingCount));
    axom::copy(&m_facetCount,
               m_firstFacetIds.data() + m_firstFacetIds.size() - 1,
               sizeof(axom::IndexType));
    // m_firstFacetIds.resize(m_crossingCount);
  }

  void computeFacets() override
  {
    AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::computeFacets");
    const auto firstFacetIdsView = m_firstFacetIds.view();
    const auto crossingParentIdsView = m_crossingParentIds.view();
    const auto crossingCasesView = m_crossingCases.view();

    // Internal contour mesh data to populate
    axom::ArrayView<axom::IndexType, 2> facetNodeIdsView = m_facetNodeIds;
    axom::ArrayView<double, 2> facetNodeCoordsView = m_facetNodeCoords;
    axom::ArrayView<axom::IndexType> facetParentIdsView = m_facetParentIds;
    const axom::IndexType facetIndexOffset = m_facetIndexOffset;

    ComputeFacets_Util cfu(m_contourVal, m_caseIdsMDMapper, m_fcnView, m_coordsViews);

    auto gen_for_parent_cell = AXOM_LAMBDA(axom::IndexType crossingId)
    {
      auto parentCellId = crossingParentIdsView[crossingId];
      auto caseId = crossingCasesView[crossingId];
      Point cornerCoords[CELL_CORNER_COUNT];
      double cornerValues[CELL_CORNER_COUNT];
      cfu.get_corner_coords_and_values(parentCellId, cornerCoords, cornerValues);

      auto additionalFacets = firstFacetIdsView[crossingId + 1] - firstFacetIdsView[crossingId];
      auto firstFacetId = facetIndexOffset + firstFacetIdsView[crossingId];
      // Keep canonical cube cases on the legacy LUT path so the exhaustive
      // single-cube regression stays stable, but switch known problematic
      // non-canonical cells to the resolved path.
      if(cfu.should_use_resolved_path(cornerValues))
      {
        Point resolvedTriangles[ComputeFacets_Util::MAX_RESOLVED_TRIANGLES][DIM];
        const int resolvedCount =
          cfu.generate_resolved_triangles(cornerCoords, cornerValues, resolvedTriangles);

        for(int triIdx = 0; triIdx < resolvedCount; ++triIdx)
        {
          const axom::IndexType newFacetId = firstFacetId + triIdx;
          const axom::IndexType firstCornerId = newFacetId * DIM;
          facetParentIdsView[newFacetId] = parentCellId;

          for(axom::IndexType d = 0; d < DIM; ++d)
          {
            const axom::IndexType newCornerId = firstCornerId + d;
            facetNodeIdsView[newFacetId][d] = newCornerId;
            for(int comp = 0; comp < DIM; ++comp)
            {
              facetNodeCoordsView(newCornerId, comp) = resolvedTriangles[triIdx][d][comp];
            }
          }
        }
      }
      else
      {
        if(cfu.has_iso_vertex(cornerValues))
        {
          Point classicTriangles[ComputeFacets_Util::MAX_RESOLVED_TRIANGLES][DIM];
          const int classicCount =
            cfu.generate_classic_triangles(caseId, cornerCoords, cornerValues, classicTriangles);

          for(int triIdx = 0; triIdx < classicCount; ++triIdx)
          {
            const axom::IndexType newFacetId = firstFacetId + triIdx;
            const axom::IndexType firstCornerId = newFacetId * DIM;
            facetParentIdsView[newFacetId] = parentCellId;

            for(axom::IndexType d = 0; d < DIM; ++d)
            {
              const axom::IndexType newCornerId = firstCornerId + d;
              facetNodeIdsView[newFacetId][d] = newCornerId;
              for(int comp = 0; comp < DIM; ++comp)
              {
                facetNodeCoordsView(newCornerId, comp) = classicTriangles[triIdx][d][comp];
              }
            }
          }
        }
        else
        {
          for(axom::IndexType fId = 0; fId < additionalFacets; ++fId)
          {
            axom::IndexType newFacetId = firstFacetId + fId;
            axom::IndexType firstCornerId = newFacetId * DIM;

            facetParentIdsView[newFacetId] = parentCellId;

            for(axom::IndexType d = 0; d < DIM; ++d)
            {
              axom::IndexType newCornerId = firstCornerId + d;
              facetNodeIdsView[newFacetId][d] = newCornerId;

              int edge = cases_table(caseId, fId * DIM + d);
              cfu.linear_interp(edge, cornerCoords, cornerValues, &facetNodeCoordsView(newCornerId, 0));
            }
          }
        }
      }
    };

    axom::for_all<ExecSpace>(0, m_crossingCount, gen_for_parent_cell);
  }

  /*!
    @brief Implementation used by MarchingCubesImpl::computeFacets().
    containing just the objects needed for that part, to be made available
    on devices.
  */
  struct ComputeFacets_Util
  {
    double contourVal;
    axom::MDMapping<DIM> mapping;
    axom::ArrayView<const double, DIM, MemorySpace> fcnView;
    axom::StackArray<axom::ArrayView<const double, DIM, MemorySpace>, DIM> coordsViews;
    ComputeFacets_Util(
      double contourVal_,
      const axom::MDMapping<DIM>& parentMDMapper,
      const axom::ArrayView<const double, DIM, MemorySpace>& fcnView_,
      const axom::StackArray<axom::ArrayView<const double, DIM, MemorySpace>, DIM> coordsViews_)
      : contourVal(contourVal_)
      , mapping(parentMDMapper)
      , fcnView(fcnView_)
      , coordsViews(coordsViews_)
    { }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2>::type get_corner_coords_and_values(
      IndexType parentCellId,
      Point cornerCoords[],
      double cornerValues[]) const
    {
      const auto& x = coordsViews[0];
      const auto& y = coordsViews[1];

      const auto c = mapping.toMultiIndex(parentCellId);
      const auto& i = c[0];
      const auto& j = c[1];

      // clang-format off
      cornerCoords[0] = { x(i  , j  ), y(i  , j  ) };
      cornerCoords[1] = { x(i+1, j  ), y(i+1, j  ) };
      cornerCoords[2] = { x(i+1, j+1), y(i+1, j+1) };
      cornerCoords[3] = { x(i  , j+1), y(i  , j+1) };

      cornerValues[0] = fcnView(i  , j  );
      cornerValues[1] = fcnView(i+1, j  );
      cornerValues[2] = fcnView(i+1, j+1);
      cornerValues[3] = fcnView(i  , j+1);
      // clang-format on
    }
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_corner_coords_and_values(
      IndexType parentCellId,
      Point cornerCoords[],
      double cornerValues[]) const
    {
      const auto& x = coordsViews[0];
      const auto& y = coordsViews[1];
      const auto& z = coordsViews[2];

      const auto c = mapping.toMultiIndex(parentCellId);
      const auto& i = c[0];
      const auto& j = c[1];
      const auto& k = c[2];

      // clang-format off
      cornerCoords[0] = { x(i+1, j  , k  ), y(i+1, j  , k  ), z(i+1, j  , k  ) };
      cornerCoords[1] = { x(i+1, j+1, k  ), y(i+1, j+1, k  ), z(i+1, j+1, k  ) };
      cornerCoords[2] = { x(i  , j+1, k  ), y(i  , j+1, k  ), z(i  , j+1, k  ) };
      cornerCoords[3] = { x(i  , j  , k  ), y(i  , j  , k  ), z(i  , j  , k  ) };
      cornerCoords[4] = { x(i+1, j  , k+1), y(i+1, j  , k+1), z(i+1, j  , k+1) };
      cornerCoords[5] = { x(i+1, j+1, k+1), y(i+1, j+1, k+1), z(i+1, j+1, k+1) };
      cornerCoords[6] = { x(i  , j+1, k+1), y(i  , j+1, k+1), z(i  , j+1, k+1) };
      cornerCoords[7] = { x(i  , j  , k+1), y(i  , j  , k+1), z(i  , j  , k+1) };

      cornerValues[0] = fcnView(i+1, j  , k  );
      cornerValues[1] = fcnView(i+1, j+1, k  );
      cornerValues[2] = fcnView(i  , j+1, k  );
      cornerValues[3] = fcnView(i  , j  , k  );
      cornerValues[4] = fcnView(i+1, j  , k+1);
      cornerValues[5] = fcnView(i+1, j+1, k+1);
      cornerValues[6] = fcnView(i  , j+1, k+1);
      cornerValues[7] = fcnView(i  , j  , k+1);
      // clang-format on
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2>::type get_corner_values(
      IndexType parentCellId,
      double cornerValues[]) const
    {
      Point unusedCornerCoords[4];
      get_corner_coords_and_values(parentCellId, unusedCornerCoords, cornerValues);
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_corner_values(
      IndexType parentCellId,
      double cornerValues[]) const
    {
      Point unusedCornerCoords[8];
      get_corner_coords_and_values(parentCellId, unusedCornerCoords, cornerValues);
    }

    static constexpr int MAX_RESOLVED_TRIANGLES = 12;
    static constexpr double ISO_REL_TOL = 1.0e-12;
    static constexpr double ISO_ABS_TOL = 1.0e-14;

    struct IsoPoint
    {
      Point pt;
      std::uint8_t kind = 0;
      std::uint8_t v0 = 0;
      std::uint8_t v1 = 0;
    };

    AXOM_HOST_DEVICE bool is_on_iso(double value) const
    {
      return axom::utilities::isNearlyEqualRelative(value, contourVal, ISO_REL_TOL, ISO_ABS_TOL);
    }

    AXOM_HOST_DEVICE bool same_scalar_value(double lhs, double rhs) const
    {
      return axom::utilities::isNearlyEqualRelative(lhs, rhs, ISO_REL_TOL, ISO_ABS_TOL);
    }

    AXOM_HOST_DEVICE bool is_near_zero(double value) const { return same_scalar_value(value, 0.0); }

    AXOM_HOST_DEVICE double relative_value(double value) const { return value - contourVal; }

    AXOM_HOST_DEVICE bool is_canonical_classic_cell(const double cornerValues[]) const
    {
      const double referenceMagnitude = axom::utilities::abs(cornerValues[0] - contourVal);
      if(is_on_iso(cornerValues[0]))
      {
        return false;
      }

      constexpr double magnitudeTol = 1.0e-12;
      for(int i = 0; i < 8; ++i)
      {
        const double magnitude = axom::utilities::abs(cornerValues[i] - contourVal);
        if(is_on_iso(cornerValues[i]) ||
           !axom::utilities::isNearlyEqual(magnitude, referenceMagnitude, magnitudeTol))
        {
          return false;
        }
      }

      return true;
    }

    AXOM_HOST_DEVICE int sign_class(double value) const
    {
      if(is_on_iso(value))
      {
        return 0;
      }

      return value > contourVal ? 1 : -1;
    }

    AXOM_HOST_DEVICE int compute_classic_case_id(const double cornerValues[]) const
    {
      int caseId = 0;
      for(int n = 0; n < 8; ++n)
      {
        if(cornerValues[n] >= contourVal)
        {
          caseId |= (1 << n);
        }
      }
      return caseId;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, int>::type contour_cell_count(
      int caseId,
      IndexType /* parentCellId */) const
    {
      return MarchingCubesImpl<DIM, ExecSpace, SequentialExecSpace>::template num_contour_cells<2>(
        caseId);
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, int>::type contour_cell_count(
      int caseId,
      IndexType parentCellId) const
    {
      double cornerValues[8];
      get_corner_values(parentCellId, cornerValues);
      if(!should_use_resolved_path(cornerValues))
      {
        if(has_iso_vertex(cornerValues))
        {
          // Exact-iso vertices can collapse multiple LUT triangles onto the
          // same geometry, so re-count after filtering duplicates.
          Point unitCubeCoords[8];
          get_unit_cube_corner_coords(unitCubeCoords);
          Point triangles[MAX_RESOLVED_TRIANGLES][DIM];
          return generate_classic_triangles(caseId, unitCubeCoords, cornerValues, triangles);
        }
        return MarchingCubesImpl<DIM, ExecSpace, SequentialExecSpace>::template num_contour_cells<3>(
          caseId);
      }

      Point unitCubeCoords[8];
      get_unit_cube_corner_coords(unitCubeCoords);
      Point triangles[MAX_RESOLVED_TRIANGLES][DIM];
      return generate_resolved_triangles(unitCubeCoords, cornerValues, triangles);
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, bool>::type should_use_resolved_path(
      const double /* cornerValues */[]) const
    {
      return false;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, bool>::type should_use_resolved_path(
      const double cornerValues[]) const
    {
      const int caseId = compute_classic_case_id(cornerValues);
      const bool canonicalClassicCell = is_canonical_classic_cell(cornerValues);

      if(is_legacy_all_iso_face_case(cornerValues))
      {
        return false;
      }

      if(is_single_iso_edge_cap_case(cornerValues))
      {
        return false;
      }

      // These interior tunnel cases are topologically ambiguous in practice
      // once non-canonical sampled values are introduced.
      if(is_interior_ambiguous_tunnel_case(caseId) && !canonicalClassicCell)
      {
        return true;
      }

      return has_iso_edge(cornerValues) ||
        (has_resolved_ambiguous_face(cornerValues) && !canonicalClassicCell);
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_unit_cube_corner_coords(
      Point cornerCoords[]) const
    {
      cornerCoords[0] = {1.0, 0.0, 0.0};
      cornerCoords[1] = {1.0, 1.0, 0.0};
      cornerCoords[2] = {0.0, 1.0, 0.0};
      cornerCoords[3] = {0.0, 0.0, 0.0};
      cornerCoords[4] = {1.0, 0.0, 1.0};
      cornerCoords[5] = {1.0, 1.0, 1.0};
      cornerCoords[6] = {0.0, 1.0, 1.0};
      cornerCoords[7] = {0.0, 0.0, 1.0};
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, bool>::type has_iso_vertex(const double[]) const
    {
      return false;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, bool>::type has_iso_vertex(
      const double cornerValues[]) const
    {
      for(int vertex = 0; vertex < 8; ++vertex)
      {
        if(is_on_iso(cornerValues[vertex]))
        {
          return true;
        }
      }

      return false;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, bool>::type has_iso_edge(
      const double cornerValues[]) const
    {
      for(int edge = 0; edge < 12; ++edge)
      {
        int v0 = 0;
        int v1 = 0;
        get_hex_edge_vertices(edge, v0, v1);
        if(is_on_iso(cornerValues[v0]) && is_on_iso(cornerValues[v1]))
        {
          return true;
        }
      }

      return false;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, bool>::type is_single_iso_edge_cap_case(
      const double cornerValues[]) const
    {
      int isoVertices[8];
      int isoVertexCount = 0;
      for(int vertex = 0; vertex < 8; ++vertex)
      {
        if(is_on_iso(cornerValues[vertex]))
        {
          isoVertices[isoVertexCount++] = vertex;
        }
      }

      if(isoVertexCount != 2)
      {
        return false;
      }

      int isoEdge = -1;
      int isoEdgeCount = 0;
      for(int edge = 0; edge < 12; ++edge)
      {
        int v0 = 0;
        int v1 = 0;
        get_hex_edge_vertices(edge, v0, v1);
        if(is_on_iso(cornerValues[v0]) && is_on_iso(cornerValues[v1]))
        {
          isoEdge = edge;
          ++isoEdgeCount;
        }
      }

      if(isoEdgeCount != 1)
      {
        return false;
      }

      int edgeV0 = 0;
      int edgeV1 = 0;
      get_hex_edge_vertices(isoEdge, edgeV0, edgeV1);
      const bool endpointsMatch = (isoVertices[0] == edgeV0 || isoVertices[0] == edgeV1) &&
        (isoVertices[1] == edgeV0 || isoVertices[1] == edgeV1);
      if(!endpointsMatch)
      {
        return false;
      }

      const int caseId = compute_classic_case_id(cornerValues);
      return caseId == 243 || caseId == 252 || caseId == 63 || caseId == 207;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, bool>::type has_resolved_ambiguous_face(
      const double cornerValues[]) const
    {
      for(int face = 0; face < 6; ++face)
      {
        int faceVertices[4];
        get_hex_face_vertices(face, faceVertices);

        int pattern = 0;
        bool hasZero = false;
        for(int i = 0; i < 4; ++i)
        {
          const int sign = sign_class(cornerValues[faceVertices[i]]);
          if(sign == 0)
          {
            hasZero = true;
            break;
          }

          if(sign > 0)
          {
            pattern |= (1 << i);
          }
        }

        if(hasZero || (pattern != 0x5 && pattern != 0xA))
        {
          continue;
        }

        const double decider = ambiguous_face_decider(cornerValues, faceVertices);
        const double centerValue = ambiguous_face_center_value(cornerValues, faceVertices);
        if(!is_near_zero(decider) || !is_near_zero(centerValue))
        {
          return true;
        }
      }

      return false;
    }

    AXOM_HOST_DEVICE bool is_interior_ambiguous_tunnel_case(int caseId) const
    {
      return caseId == 60 || caseId == 90 || caseId == 102 || caseId == 153 || caseId == 165 ||
        caseId == 195;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, bool>::type is_legacy_all_iso_face_case(
      const double cornerValues[]) const
    {
      for(int face = 0; face < 6; ++face)
      {
        int faceVertices[4];
        int oppositeVertices[4];
        get_hex_face_vertices(face, faceVertices);
        get_hex_opposite_face_vertices(face, oppositeVertices);

        bool allFaceIso = true;
        for(int i = 0; i < 4; ++i)
        {
          if(!is_on_iso(cornerValues[faceVertices[i]]))
          {
            allFaceIso = false;
            break;
          }
        }

        if(!allFaceIso)
        {
          continue;
        }

        int oppositeSign = 0;
        bool uniformOppositeSide = true;
        for(int i = 0; i < 4; ++i)
        {
          const int sign = sign_class(cornerValues[oppositeVertices[i]]);
          if(sign == 0)
          {
            uniformOppositeSide = false;
            break;
          }

          if(i == 0)
          {
            oppositeSign = sign;
          }
          else if(sign != oppositeSign)
          {
            uniformOppositeSide = false;
            break;
          }
        }

        if(uniformOppositeSide)
        {
          return true;
        }
      }

      return false;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, int>::type
    generate_resolved_triangles(const Point[4], const double[4], Point triangles[][DIM]) const
    {
      AXOM_UNUSED_VAR(triangles);
      return 0;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, int>::type generate_resolved_triangles(
      const Point cornerCoords[8],
      const double cornerValues[8],
      Point triangles[][DIM]) const
    {
      if(!has_iso_vertex(cornerValues) && !has_iso_edge(cornerValues))
      {
        return generate_face_resolved_triangles(cornerCoords, cornerValues, triangles);
      }

      return generate_tet_resolved_triangles(cornerCoords, cornerValues, triangles);
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, int>::type
    generate_classic_triangles(int, const Point[4], const double[4], Point triangles[][DIM]) const
    {
      AXOM_UNUSED_VAR(triangles);
      return 0;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, int>::type generate_classic_triangles(
      int caseId,
      const Point cornerCoords[8],
      const double cornerValues[8],
      Point triangles[][DIM]) const
    {
      const int lutTriangleCount =
        MarchingCubesImpl<DIM, ExecSpace, SequentialExecSpace>::template num_contour_cells<3>(caseId);
      int triangleCount = 0;

      for(int triIdx = 0; triIdx < lutTriangleCount; ++triIdx)
      {
        Point triangle[DIM];
        for(int d = 0; d < DIM; ++d)
        {
          const int edge = cases_table(caseId, triIdx * DIM + d);
          linear_interp(edge, cornerCoords, cornerValues, &triangle[d][0]);
        }
        // Iso-vertex snapping can make distinct LUT entries geometrically
        // identical, so filter duplicates before writing facets.
        append_triangle_if_unique(triangles, triangleCount, MAX_RESOLVED_TRIANGLES, triangle);
      }

      return triangleCount;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, int>::type generate_tet_resolved_triangles(
      const Point cornerCoords[8],
      const double cornerValues[8],
      Point triangles[][DIM]) const
    {
      Point generated[MAX_RESOLVED_TRIANGLES][DIM];
      int triangleCount = 0;

      for(int tet = 0; tet < 6; ++tet)
      {
        int tetVertices[4];
        get_cube_tet_vertices(tet, tetVertices);
        append_tet_triangles(cornerCoords,
                             cornerValues,
                             tetVertices,
                             generated,
                             triangleCount,
                             MAX_RESOLVED_TRIANGLES);
      }

      for(int triIdx = 0; triIdx < triangleCount; ++triIdx)
      {
        for(int v = 0; v < DIM; ++v)
        {
          triangles[triIdx][v] = generated[triIdx][v];
        }
      }

      return triangleCount;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, int>::type generate_face_resolved_triangles(
      const Point cornerCoords[8],
      const double cornerValues[8],
      Point triangles[][DIM]) const
    {
      bool edgeHasCrossing[12];
      Point edgeCrossings[12];
      for(int edge = 0; edge < 12; ++edge)
      {
        edgeHasCrossing[edge] = false;
        int v0 = 0;
        int v1 = 0;
        get_hex_edge_vertices(edge, v0, v1);
        if(sign_class(cornerValues[v0]) != sign_class(cornerValues[v1]))
        {
          edgeHasCrossing[edge] = true;
          linear_interp_vertices(v0, v1, cornerCoords, cornerValues, &edgeCrossings[edge][0]);
        }
      }

      int neighbors[12][2];
      int degree[12];
      for(int edge = 0; edge < 12; ++edge)
      {
        degree[edge] = 0;
        neighbors[edge][0] = -1;
        neighbors[edge][1] = -1;
      }

      for(int face = 0; face < 6; ++face)
      {
        int faceVertices[4];
        int faceEdges[4];
        get_hex_face_vertices(face, faceVertices);
        get_hex_face_edges(face, faceEdges);

        int crossingEdges[4];
        int crossingCount = 0;
        for(int i = 0; i < 4; ++i)
        {
          if(edgeHasCrossing[faceEdges[i]])
          {
            crossingEdges[crossingCount++] = i;
          }
        }

        if(crossingCount == 2)
        {
          connect_face_edges(faceEdges[crossingEdges[0]],
                             faceEdges[crossingEdges[1]],
                             neighbors,
                             degree);
        }
        else if(crossingCount == 4)
        {
          if(connect_ambiguous_face_diagonal_02(cornerValues, faceVertices))
          {
            connect_face_edges(faceEdges[0], faceEdges[1], neighbors, degree);
            connect_face_edges(faceEdges[2], faceEdges[3], neighbors, degree);
          }
          else
          {
            connect_face_edges(faceEdges[1], faceEdges[2], neighbors, degree);
            connect_face_edges(faceEdges[3], faceEdges[0], neighbors, degree);
          }
        }
      }

      bool visited[12];
      for(int edge = 0; edge < 12; ++edge)
      {
        visited[edge] = false;
      }

      int triangleCount = 0;
      for(int start = 0; start < 12; ++start)
      {
        if(!edgeHasCrossing[start] || degree[start] == 0 || visited[start])
        {
          continue;
        }

        int polygonEdges[12];
        int polygonSize = 0;
        int prev = -1;
        int curr = start;

        while(curr >= 0 && !visited[curr] && polygonSize < 12)
        {
          polygonEdges[polygonSize++] = curr;
          visited[curr] = true;

          int next = neighbors[curr][0];
          if(next == prev)
          {
            next = neighbors[curr][1];
          }

          prev = curr;
          curr = next;
          if(curr == start)
          {
            break;
          }
        }

        if(polygonSize < 3)
        {
          continue;
        }

        for(int i = 1; i + 1 < polygonSize && triangleCount < MAX_RESOLVED_TRIANGLES; ++i)
        {
          triangles[triangleCount][0] = edgeCrossings[polygonEdges[0]];
          triangles[triangleCount][1] = edgeCrossings[polygonEdges[i]];
          triangles[triangleCount][2] = edgeCrossings[polygonEdges[i + 1]];
          ++triangleCount;
        }
      }

      return triangleCount;
    }

    AXOM_HOST_DEVICE double ambiguous_face_decider(const double cornerValues[],
                                                   const int faceVertices[4]) const
    {
      // Asymptotic decider for a bilinear quad face, using values relative to
      // the isovalue.
      const double f0 = relative_value(cornerValues[faceVertices[0]]);
      const double f1 = relative_value(cornerValues[faceVertices[1]]);
      const double f2 = relative_value(cornerValues[faceVertices[2]]);
      const double f3 = relative_value(cornerValues[faceVertices[3]]);
      return f0 * f2 - f1 * f3;
    }

    AXOM_HOST_DEVICE double ambiguous_face_center_value(const double cornerValues[],
                                                        const int faceVertices[4]) const
    {
      return 0.25 *
        (relative_value(cornerValues[faceVertices[0]]) +
         relative_value(cornerValues[faceVertices[1]]) +
         relative_value(cornerValues[faceVertices[2]]) +
         relative_value(cornerValues[faceVertices[3]]));
    }

    AXOM_HOST_DEVICE bool connect_ambiguous_face_diagonal_02(const double cornerValues[],
                                                             const int faceVertices[4]) const
    {
      const double decider = ambiguous_face_decider(cornerValues, faceVertices);
      if(!is_near_zero(decider))
      {
        return decider > 0.0;
      }

      // Degenerate saddle: preserve the previous center-sign tie-breaker.
      const bool positiveDiag02 =
        cornerValues[faceVertices[0]] > contourVal && cornerValues[faceVertices[2]] > contourVal;
      const bool centerPositive = ambiguous_face_center_value(cornerValues, faceVertices) > 0.0;
      return (positiveDiag02 && centerPositive) || (!positiveDiag02 && !centerPositive);
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_hex_edge_vertices(int edgeIdx,
                                                                                    int& v0,
                                                                                    int& v1) const
    {
      const int hex_edge_table[] = {
        0, 1, 1, 2, 2, 3, 3, 0,  // base
        4, 5, 5, 6, 6, 7, 7, 4,  // top
        0, 4, 1, 5, 2, 6, 3, 7   // vertical
      };

      v0 = hex_edge_table[2 * edgeIdx];
      v1 = hex_edge_table[2 * edgeIdx + 1];
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_hex_face_vertices(
      int faceIdx,
      int faceVertices[4]) const
    {
      const int face_vertex_table[6][4] =
        {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 4, 5, 1}, {3, 2, 6, 7}, {1, 5, 6, 2}, {0, 3, 7, 4}};

      for(int i = 0; i < 4; ++i)
      {
        faceVertices[i] = face_vertex_table[faceIdx][i];
      }
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_hex_face_edges(int faceIdx,
                                                                                 int faceEdges[4]) const
    {
      const int face_edge_table[6][4] =
        {{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 4, 9, 0}, {2, 10, 6, 11}, {9, 5, 10, 1}, {3, 11, 7, 8}};

      for(int i = 0; i < 4; ++i)
      {
        faceEdges[i] = face_edge_table[faceIdx][i];
      }
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type connect_face_edges(int edgeA,
                                                                                 int edgeB,
                                                                                 int neighbors[12][2],
                                                                                 int degree[12]) const
    {
      if(degree[edgeA] < 2)
      {
        neighbors[edgeA][degree[edgeA]++] = edgeB;
      }

      if(degree[edgeB] < 2)
      {
        neighbors[edgeB][degree[edgeB]++] = edgeA;
      }
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_hex_opposite_face_vertices(
      int faceIdx,
      int faceVertices[4]) const
    {
      const int opposite_face_table[6][4] =
        {{4, 5, 6, 7}, {0, 1, 2, 3}, {3, 7, 6, 2}, {0, 1, 5, 4}, {0, 4, 7, 3}, {1, 2, 6, 5}};

      for(int i = 0; i < 4; ++i)
      {
        faceVertices[i] = opposite_face_table[faceIdx][i];
      }
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_cube_tet_vertices(int tetIdx,
                                                                                    int tetVertices[4]) const
    {
      const int tet_vertex_table[6][4] =
        {{3, 0, 1, 5}, {3, 0, 4, 5}, {3, 7, 4, 5}, {3, 7, 6, 5}, {3, 2, 6, 5}, {3, 2, 1, 5}};

      for(int i = 0; i < 4; ++i)
      {
        tetVertices[i] = tet_vertex_table[tetIdx][i];
      }
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_tet_edge_vertices(int edgeIdx,
                                                                                    int& v0,
                                                                                    int& v1) const
    {
      const int tet_edge_table[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
      v0 = tet_edge_table[edgeIdx][0];
      v1 = tet_edge_table[edgeIdx][1];
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type get_tet_face_vertices(
      int faceIdx,
      int faceVertices[3]) const
    {
      const int tet_face_table[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
      for(int i = 0; i < 3; ++i)
      {
        faceVertices[i] = tet_face_table[faceIdx][i];
      }
    }

    AXOM_HOST_DEVICE bool same_point(const Point& lhs, const Point& rhs) const
    {
      for(int d = 0; d < DIM; ++d)
      {
        if(!axom::utilities::isNearlyEqual(lhs[d], rhs[d]))
        {
          return false;
        }
      }
      return true;
    }

    AXOM_HOST_DEVICE bool point_precedes_lexicographically(const Point& lhs, const Point& rhs) const
    {
      for(int d = 0; d < DIM; ++d)
      {
        if(lhs[d] < rhs[d])
        {
          return true;
        }
        if(rhs[d] < lhs[d])
        {
          return false;
        }
      }

      return false;
    }

    AXOM_HOST_DEVICE bool same_triangle(const Point lhs[DIM], const Point rhs[DIM]) const
    {
      bool matched[DIM];
      for(int i = 0; i < DIM; ++i)
      {
        matched[i] = false;
      }
      for(int i = 0; i < DIM; ++i)
      {
        bool found = false;
        for(int j = 0; j < DIM; ++j)
        {
          if(!matched[j] && same_point(lhs[i], rhs[j]))
          {
            matched[j] = true;
            found = true;
            break;
          }
        }

        if(!found)
        {
          return false;
        }
      }
      return true;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, double>::type triangle_measure(
      const Point triangle[DIM]) const
    {
      const Point& a = triangle[0];
      const Point& b = triangle[1];
      const Point& c = triangle[2];

      const double abx = b[0] - a[0];
      const double aby = b[1] - a[1];
      const double abz = b[2] - a[2];
      const double acx = c[0] - a[0];
      const double acy = c[1] - a[1];
      const double acz = c[2] - a[2];

      const double crossX = aby * acz - abz * acy;
      const double crossY = abz * acx - abx * acz;
      const double crossZ = abx * acy - aby * acx;
      return crossX * crossX + crossY * crossY + crossZ * crossZ;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, bool>::type
    add_unique_iso_point(IsoPoint points[], int& pointCount, const IsoPoint& candidate) const
    {
      for(int i = 0; i < pointCount; ++i)
      {
        if(same_point(points[i].pt, candidate.pt))
        {
          return false;
        }
      }

      points[pointCount++] = candidate;
      return true;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, bool>::type point_lies_on_tet_face(
      const IsoPoint& point,
      const int tetVertices[4],
      const int faceVertices[3]) const
    {
      if(point.kind == 0)
      {
        for(int i = 0; i < 3; ++i)
        {
          if(point.v0 == tetVertices[faceVertices[i]])
          {
            return true;
          }
        }
        return false;
      }

      bool hasV0 = false;
      bool hasV1 = false;
      for(int i = 0; i < 3; ++i)
      {
        const int faceVertex = tetVertices[faceVertices[i]];
        hasV0 |= point.v0 == faceVertex;
        hasV1 |= point.v1 == faceVertex;
      }
      return hasV0 && hasV1;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type append_triangle_if_unique(
      Point triangles[][DIM],
      int& triangleCount,
      int triangleCapacity,
      const Point triangle[DIM]) const
    {
      if(triangle_measure(triangle) <= axom::primal::PRIMAL_TINY)
      {
        return;
      }

      for(int triIdx = 0; triIdx < triangleCount; ++triIdx)
      {
        if(same_triangle(triangles[triIdx], triangle))
        {
          return;
        }
      }

      if(triangleCount < triangleCapacity)
      {
        for(int i = 0; i < DIM; ++i)
        {
          triangles[triangleCount][i] = triangle[i];
        }
        ++triangleCount;
      }
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type append_tet_triangles(
      const Point cornerCoords[8],
      const double cornerValues[8],
      const int tetVertices[4],
      Point triangles[][DIM],
      int& triangleCount,
      int triangleCapacity) const
    {
      int zeroVertexCount = 0;
      for(int i = 0; i < 4; ++i)
      {
        zeroVertexCount += is_on_iso(cornerValues[tetVertices[i]]) ? 1 : 0;
      }

      if(zeroVertexCount == 4)
      {
        return;
      }

      IsoPoint points[4];
      int pointCount = 0;

      for(int edgeIdx = 0; edgeIdx < 6; ++edgeIdx)
      {
        int localV0 = 0;
        int localV1 = 0;
        get_tet_edge_vertices(edgeIdx, localV0, localV1);

        const int v0 = tetVertices[localV0];
        const int v1 = tetVertices[localV1];
        const int sign0 = sign_class(cornerValues[v0]);
        const int sign1 = sign_class(cornerValues[v1]);

        if(sign0 == 0 && sign1 == 0)
        {
          IsoPoint p0;
          p0.pt = cornerCoords[v0];
          p0.kind = 0;
          p0.v0 = static_cast<std::uint8_t>(v0);
          add_unique_iso_point(points, pointCount, p0);

          IsoPoint p1;
          p1.pt = cornerCoords[v1];
          p1.kind = 0;
          p1.v0 = static_cast<std::uint8_t>(v1);
          add_unique_iso_point(points, pointCount, p1);
          continue;
        }

        if(sign0 == 0 || sign1 == 0 || sign0 != sign1)
        {
          IsoPoint point;
          point.kind = (sign0 == 0 || sign1 == 0) ? 0 : 1;
          if(point.kind == 0)
          {
            point.v0 = static_cast<std::uint8_t>(sign0 == 0 ? v0 : v1);
            point.pt = cornerCoords[point.v0];
          }
          else
          {
            point.v0 = static_cast<std::uint8_t>((v0 < v1) ? v0 : v1);
            point.v1 = static_cast<std::uint8_t>((v0 < v1) ? v1 : v0);
            linear_interp_vertices(v0, v1, cornerCoords, cornerValues, &point.pt[0]);
          }

          add_unique_iso_point(points, pointCount, point);
        }
      }

      if(pointCount < 3)
      {
        return;
      }

      if(pointCount == 3)
      {
        Point triangle[DIM] = {points[0].pt, points[1].pt, points[2].pt};
        append_triangle_if_unique(triangles, triangleCount, triangleCapacity, triangle);
        return;
      }

      int neighbors[4][2];
      int degree[4] = {0, 0, 0, 0};
      for(int i = 0; i < 4; ++i)
      {
        neighbors[i][0] = -1;
        neighbors[i][1] = -1;
      }

      for(int faceIdx = 0; faceIdx < 4; ++faceIdx)
      {
        int faceVertices[3];
        get_tet_face_vertices(faceIdx, faceVertices);

        int facePointIds[4];
        int facePointCount = 0;
        for(int pointIdx = 0; pointIdx < pointCount; ++pointIdx)
        {
          if(point_lies_on_tet_face(points[pointIdx], tetVertices, faceVertices))
          {
            facePointIds[facePointCount++] = pointIdx;
          }
        }

        if(facePointCount == 2)
        {
          const int a = facePointIds[0];
          const int b = facePointIds[1];
          if(degree[a] < 2 && (neighbors[a][0] != b && neighbors[a][1] != b))
          {
            neighbors[a][degree[a]++] = b;
          }
          if(degree[b] < 2 && (neighbors[b][0] != a && neighbors[b][1] != a))
          {
            neighbors[b][degree[b]++] = a;
          }
        }
      }

      int ordered[4] = {0, -1, -1, -1};
      int orderedCount = 1;
      int prev = -1;
      int curr = 0;
      while(orderedCount < pointCount)
      {
        int next = -1;
        for(int i = 0; i < degree[curr]; ++i)
        {
          if(neighbors[curr][i] != prev)
          {
            next = neighbors[curr][i];
            break;
          }
        }

        if(next < 0)
        {
          break;
        }

        ordered[orderedCount++] = next;
        prev = curr;
        curr = next;
      }

      if(orderedCount != pointCount)
      {
        Point triangleA[DIM] = {points[0].pt, points[1].pt, points[2].pt};
        Point triangleB[DIM] = {points[0].pt, points[2].pt, points[3].pt};
        append_triangle_if_unique(triangles, triangleCount, triangleCapacity, triangleA);
        append_triangle_if_unique(triangles, triangleCount, triangleCapacity, triangleB);
        return;
      }

      for(int i = 1; i + 1 < orderedCount; ++i)
      {
        Point triangle[DIM] = {points[ordered[0]].pt, points[ordered[i]].pt, points[ordered[i + 1]].pt};
        append_triangle_if_unique(triangles, triangleCount, triangleCapacity, triangle);
      }
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type linear_interp_vertices(
      int n1,
      int n2,
      const Point cornerCoords[8],
      const double nodeValues[8],
      double* crossingPt) const
    {
      int a = n1;
      int b = n2;
      // Normalize edge endpoint order so neighboring cells interpolate shared
      // edges identically even when they visit the edge in opposite directions.
      if(point_precedes_lexicographically(cornerCoords[b], cornerCoords[a]) ||  //
         (same_point(cornerCoords[a], cornerCoords[b]) && b < a))
      {
        const int tmp = a;
        a = b;
        b = tmp;
      }

      const double f1 = nodeValues[a];
      const double f2 = nodeValues[b];
      const Point& p1 = cornerCoords[a];
      const Point& p2 = cornerCoords[b];

      if(is_on_iso(f1))
      {
        for(int d = 0; d < DIM; ++d)
        {
          crossingPt[d] = p1[d];
        }
        return;
      }

      if(is_on_iso(f2))
      {
        for(int d = 0; d < DIM; ++d)
        {
          crossingPt[d] = p2[d];
        }
        return;
      }

      if(same_scalar_value(f1, f2))
      {
        for(int d = 0; d < DIM; ++d)
        {
          crossingPt[d] = 0.5 * (p1[d] + p2[d]);
        }
        return;
      }

      double weight = (contourVal - f1) / (f2 - f1);
      if(weight < 0.0)
      {
        weight = 0.0;
      }
      else if(weight > 1.0)
      {
        weight = 1.0;
      }

      for(int d = 0; d < DIM; ++d)
      {
        crossingPt[d] = p1[d] + weight * (p2[d] - p1[d]);
      }
    }

    //!@brief Interpolate for the contour location crossing a parent edge.
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2>::type linear_interp(
      int edgeIdx,
      const Point cornerCoords[4],
      const double nodeValues[4],
      double* /* Point& */ crossingPt) const
    {
      // STEP 0: get the edge node indices
      // 2 nodes define the edge.  n1 and n2 are the indices of
      // the nodes w.r.t. the square or cubic zone.  There is a
      // agreed-on ordering of these indices in the arrays xx, yy,
      // zz, nodeValues, crossingPt.
      int n1 = edgeIdx;
      int n2 = (edgeIdx == 3) ? 0 : edgeIdx + 1;

      // STEP 1: get the fields and coordinates from the two points
      const double f1 = nodeValues[n1];
      const double f2 = nodeValues[n2];

      const Point& p1 = cornerCoords[n1];
      const Point& p2 = cornerCoords[n2];

      // STEP 2: check whether the interpolated point is at one of the two corners.
      if(is_on_iso(f1))
      {
        crossingPt[0] = p1[0];
        crossingPt[1] = p1[1];  // crossingPt = p1;
        return;
      }

      if(is_on_iso(f2))
      {
        crossingPt[0] = p2[0];
        crossingPt[1] = p2[1];  // crossingPt = p2;
        return;
      }

      if(same_scalar_value(f1, f2))
      {
        crossingPt[0] = 0.5 * (p1[0] + p2[0]);
        crossingPt[1] = 0.5 * (p1[1] + p2[1]);
        return;
      }

      // STEP 3: point is in between the edge points, interpolate its position
      double w = (contourVal - f1) / (f2 - f1);
      if(w < 0.0)
      {
        w = 0.0;
      }
      else if(w > 1.0)
      {
        w = 1.0;
      }
      for(int d = 0; d < DIM; ++d)
      {
        crossingPt[d] = p1[d] + w * (p2[d] - p1[d]);
      }
    }

    //!@brief Interpolate for the contour location crossing a parent edge.
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type linear_interp(
      int edgeIdx,
      const Point cornerCoords[8],
      const double nodeValues[8],
      double* /* Point& */ crossingPt) const
    {
      // STEP 0: get the edge node indices
      // 2 nodes define the edge.  n1 and n2 are the indices of
      // the nodes w.r.t. the square or cubic zone.  There is a
      // agreed-on ordering of these indices in the arrays
      // cornerCoords, nodeValues, hex_edge_table.
      const int hex_edge_table[] = {
        0, 1, 1, 2, 2, 3, 3, 0,  // base
        4, 5, 5, 6, 6, 7, 7, 4,  // top
        0, 4, 1, 5, 2, 6, 3, 7   // vertical
      };

      int n1 = hex_edge_table[edgeIdx * 2];
      int n2 = hex_edge_table[edgeIdx * 2 + 1];

      linear_interp_vertices(n1, n2, cornerCoords, nodeValues, crossingPt);
    }
  };  // ComputeFacets_Util

  // These 4 functions provide access to the look-up table
  // whether on host or device.  Is there a more elegant way
  // to put static 1D and 2D arrays on both host and device?  BTNG.

  template <int TDIM = DIM>
  static AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 2, int>::type num_contour_cells(int iCase)
  {
#define _MC_LOOKUP_NUM_SEGMENTS
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_NUM_SEGMENTS
    SLIC_ASSERT(iCase >= 0 && iCase < 16);
    return num_segments[iCase];
  }

  template <int TDIM = DIM>
  static AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 2, int>::type cases_table(int iCase,
                                                                                           int iEdge)
  {
#define _MC_LOOKUP_CASES2D
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_CASES2D
    SLIC_ASSERT(iCase >= 0 && iCase < 16);
    return cases2D[iCase][iEdge];
  }

  template <int TDIM = DIM>
  static AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 3, int>::type num_contour_cells(int iCase)
  {
#define _MC_LOOKUP_NUM_TRIANGLES
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_NUM_TRIANGLES
    SLIC_ASSERT(iCase >= 0 && iCase < 256);
    return num_triangles[iCase];
  }

  template <int TDIM = DIM>
  static AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 3, int>::type cases_table(int iCase,
                                                                                           int iEdge)
  {
#define _MC_LOOKUP_CASES3D
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_CASES3D
    SLIC_ASSERT(iCase >= 0 && iCase < 256);
    return cases3D[iCase][iEdge];
  }

  //!@brief Compute the case index into cases2D or cases3D.
  AXOM_HOST_DEVICE inline int compute_crossing_case(const double* f) const
  {
    int index = 0;
    for(int n = 0; n < CELL_CORNER_COUNT; ++n)
    {
      if(f[n] >= m_contourVal)
      {
        const int bit = (1 << n);
        index |= bit;
      }
    }
    return index;
  }

  /*!
    @brief Constructor.
  */
  MarchingCubesImpl() { }

  /*!
    @brief Clear computed data (without deallocating memory).

    After clearing, you can change the field, contour value
    and recompute the contour.
  */
  void clearDomain() override
  {
    m_caseIdsFlat.clear();
    m_crossingFlags.clear();
    m_scannedFlags.clear();
    m_crossingParentIds.clear();
    m_facetIncrs.clear();
    m_firstFacetIds.clear();
    m_crossingCount = 0;
    m_facetCount = 0;
  }

private:
  const int m_allocatorID;

  axom::quest::MeshViewUtil<DIM, MemorySpace> m_mvu;
  MIdx m_bShape;  //!< @brief Blueprint cell data shape.

  // Views of parent domain data.
  // DIM coordinate components, each on a DIM-dimensional mesh.
  using CoordViews = axom::StackArray<axom::ArrayView<const double, DIM, MemorySpace>, DIM>;
  CoordViews m_coordsViews;
  axom::ArrayView<const double, DIM, MemorySpace> m_fcnView;
  axom::ArrayView<const int, DIM, MemorySpace> m_maskView;

  /*!
    @brief Crossing case for each computational mesh cell.

    This is a multidim view into 1D data from m_caseIdsFlat,
    set up with help from m_caseIdsMDMapper.
  */
  axom::ArrayView<std::uint16_t, DIM, MemorySpace> m_caseIds;
  /*!
    @brief Multidim mapping to handle data ordering in
    m_caseIdsFlat.

    We want caseIds ordering to match m_fcnView, but Array
    only supports column-major ordering currently.  To control
    the order, we put caseIds in a 1D array and construct a
    multidim view with the ordering we want.
  */
  axom::MDMapping<DIM> m_caseIdsMDMapper;

  // Array references refer to shared Arrays in MarchingCubes.

  //!@brief Crossing case for each computational mesh cell.
  axom::Array<std::uint16_t>& m_caseIdsFlat;

  //!@brief Whether a parent cell crosses the contour.
  axom::Array<CrossingFlagType>& m_crossingFlags;

  //!@brief Prefix sum of m_crossingFlags
  axom::Array<axom::IndexType>& m_scannedFlags;

  /*!
    TODO: Facet increment values lie in [0, 5], but are wastefully
    stored in 32 bits because of a ROCM scan implementation that adds
    them in the input type without promoting them to the bigger output
    type.  When ROCM supports the promotion and RAJA uses it, we can
    change this type to something more efficient.
  */
  using FacetIncrsType = axom::IndexType;
  //!@brief Number of surface mesh facets added by each crossing.
  axom::Array<FacetIncrsType>& m_facetIncrs;

  //!@brief Number of parent cells crossing the contour surface.
  axom::IndexType m_crossingCount = 0;

  //!@brief Number of contour surface cells from all crossings.
  axom::IndexType m_facetCount = 0;
  axom::IndexType getContourCellCount() const override { return m_facetCount; }

  //!@brief Case ids for found crossings.
  axom::Array<std::int16_t> m_crossingCases;

  //!@brief Parent cell id (flat index into m_caseIds) for each crossing.
  axom::Array<axom::IndexType, 1, MemorySpace> m_crossingParentIds;

  //!@brief First index of facets for each crossing.
  axom::Array<axom::IndexType, 1, MemorySpace> m_firstFacetIds;

  //!@brief Number of corners (nodes) on each parent cell.
  static constexpr std::uint8_t CELL_CORNER_COUNT = (DIM == 3) ? 8 : 4;

  double m_contourVal = 0.0;
  int m_maskVal = 1;

  axom::StackArray<axom::IndexType, DIM> emptyShape()
  {
    axom::StackArray<axom::IndexType, DIM> rval;
    for(int d = 0; d < DIM; ++d)
    {
      rval[d] = 0;
    }
    return rval;
  }
};

}  // namespace marching_cubes
}  // namespace detail
}  // namespace quest
}  // namespace axom
