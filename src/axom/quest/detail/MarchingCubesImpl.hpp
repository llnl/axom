// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

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
    {
      AXOM_ANNOTATE_SCOPE("MarchingCubesImpl::scanCrossings:set_flags");
      axom::for_all<ExecSpace>(
        0,
        parentCellCount,
        AXOM_LAMBDA(axom::IndexType parentCellId) {
          auto numContourCells = num_contour_cells(caseIdsView.flatIndex(parentCellId));
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
          auto facetIncr = num_contour_cells(caseId);
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
    axom::ReduceSum<ExecSpace, axom::IndexType> vsum(0);
    axom::for_all<ExecSpace>(
      parentCellCount,
      AXOM_LAMBDA(axom::IndexType n) { vsum += bool(num_contour_cells(caseIdsView.flatIndex(n))); });
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
      auto ccc = num_contour_cells(caseId);
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
      This contrived for_all forces it to run sequentially.
    */
    axom::for_all<axom::SEQ_EXEC>(1, [=] AXOM_HOST_DEVICE(axom::IndexType /* i */) {
      *crossingId = 0;
      for(axom::IndexType n = 0; n < parentCellCount; ++n)
      {
        loopBody(n);
      }
    });
    SLIC_ASSERT(*crossingId == m_crossingCount);

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
      if(axom::utilities::isNearlyEqual(contourVal, f1) || axom::utilities::isNearlyEqual(f1, f2))
      {
        crossingPt[0] = p1[0];
        crossingPt[1] = p1[1];  // crossingPt = p1;
        return;
      }

      if(axom::utilities::isNearlyEqual(contourVal, f2))
      {
        crossingPt[0] = p2[0];
        crossingPt[1] = p2[1];  // crossingPt = p2;
        return;
      }

      // STEP 3: point is in between the edge points, interpolate its position
      constexpr double ptiny = axom::primal::PRIMAL_TINY;
      const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
      const double w = (contourVal - f1) / df;
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

      // STEP 1: get the fields and coordinates from the two points
      const double f1 = nodeValues[n1];
      const double f2 = nodeValues[n2];

      const Point& p1 = cornerCoords[n1];
      const Point& p2 = cornerCoords[n2];

      // STEP 2: check whether the interpolated point is at one of the two corners.
      if(axom::utilities::isNearlyEqual(contourVal, f1) || axom::utilities::isNearlyEqual(f1, f2))
      {
        crossingPt[0] = p1[0];
        crossingPt[1] = p1[1];
        crossingPt[2] = p1[2];  // crossingPt = p1;
        return;
      }

      if(axom::utilities::isNearlyEqual(contourVal, f2))
      {
        crossingPt[0] = p2[0];
        crossingPt[1] = p2[1];
        crossingPt[2] = p2[2];  // crossingPt = p2;
        return;
      }

      // STEP 3: point is not at corner; interpolate its position
      constexpr double ptiny = axom::primal::PRIMAL_TINY;
      const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
      const double w = (contourVal - f1) / df;
      for(int d = 0; d < DIM; ++d)
      {
        crossingPt[d] = p1[d] + w * (p2[d] - p1[d]);
      }
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
