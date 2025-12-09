// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MESHCLIPPERIMPL_HPP_
#define AXOM_MESHCLIPPERIMPL_HPP_

#include "axom/config.hpp"

#ifndef AXOM_USE_RAJA
  #error "quest::MeshClipper requires RAJA."
#endif

#include "axom/core/numerics/matvecops.hpp"
#include "axom/quest/MeshClipperStrategy.hpp"
#include "axom/quest/MeshClipper.hpp"
#include "axom/spin/BVH.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"
#include "RAJA/RAJA.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{
namespace detail
{

/*!
 * @brief Implementation of MeshClipper::Impl
 *
 * This class should be thought of as a part of the MeshClipper code,
 * even though it's in a different file.  Abstract base class
 * MeshClipper::Impl defines interfaces for MeshClipper methods that
 * should be implemented in the same execution space.  This class
 * implements those methods with the execution space as a template
 * parameter.
 */
template <typename ExecSpace>
class MeshClipperImpl : public MeshClipper::Impl
{
public:
  using LabelType = MeshClipper::LabelType;
  using Point3DType = primal::Point<double, 3>;
  using BoundingBoxType = primal::BoundingBox<double, 3>;
  using TetrahedronType = primal::Tetrahedron<double, 3>;
  using OctahedronType = primal::Octahedron<double, 3>;

  MeshClipperImpl(MeshClipper& clipper) : MeshClipper::Impl(clipper) { }

  void initVolumeOverlaps(const axom::ArrayView<MeshClipperStrategy::LabelType>& labels,
                          axom::ArrayView<double> ovlap) override
  {
    const axom::IndexType cellCount = getShapeMesh().getCellCount();
    SLIC_ASSERT(labels.size() == cellCount);
    SLIC_ASSERT(ovlap.size() == cellCount);

    auto cellVolumes = getShapeMesh().getCellVolumes();

    /*
     * Overlap volumes is cell volume for cells inside geometry.
     * and zero for cells outside geometry.
     * Cells on boundary are zeroed for accumulating by clipping process.
    */
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        auto& l = labels[i];
        ovlap[i] = l == LabelType::LABEL_IN ? cellVolumes[i] : 0.0;
      });

    return;
  }

  void zeroVolumeOverlaps(axom::ArrayView<double> ovlap) override
  {
    SLIC_ASSERT(ovlap.size() == getShapeMesh().getCellCount());
    ovlap.fill(0.0);
    return;
  }

  void addVolumesOfInteriorTets(axom::ArrayView<const axom::IndexType> cellsOnBdry,
                                axom::ArrayView<const LabelType> tetLabels,
                                axom::ArrayView<double> ovlap) override
  {
    auto meshTets = getShapeMesh().getCellsAsTets();

    const axom::IndexType hexCount = cellsOnBdry.size();

    SLIC_ASSERT(tetLabels.size() == NUM_TETS_PER_HEX * hexCount);

    axom::for_all<ExecSpace>(
      hexCount,
      AXOM_LAMBDA(axom::IndexType ih) {
        const axom::IndexType hexId = cellsOnBdry[ih];
        const LabelType* tetLabelsForHex = &tetLabels[NUM_TETS_PER_HEX * ih];
        for(int it = 0; it < NUM_TETS_PER_HEX; ++it)
        {
          if(tetLabelsForHex[it] == LabelType::LABEL_IN)
          {
            const axom::IndexType tetId = hexId * NUM_TETS_PER_HEX + it;
            const auto& tet = meshTets[tetId];
            ovlap[hexId] += tet.volume();
          }
        }
      });
  }

  //! @brief Make a list of indices where labels have value LABEL_ON.
  void collectOnIndices(const axom::ArrayView<LabelType>& labels,
                        axom::Array<axom::IndexType>& onIndices) override
  {
    if(labels.empty())
    {
      return;
    };

    AXOM_ANNOTATE_SCOPE("MeshClipper:collect_unlabeleds");
    /*!
     * 1. Generate tmpLabels, having a value of 1 where labels is LABEL_ON and zero elsewhere.
     * 2. Inclusive scan on tmpLabels to generate values that step up at LABEL_ON cells.
     * 3. Find unlabeled cells by seeing where tmpLabels changes values.
     *    (Handle first cell separately, then loop from second cell on.)
     *    Note that tmpLabels holds non-decreasing values.  By populating
     *    onIndices based on where tmpLabels changes, we never write to
     *    the same index more than once.  Write conflicts are thus avoided.
     *    Thanks to Jason Burmark for recommending this approach.
     */
    using ScanPolicy = typename axom::execution_space<ExecSpace>::loop_policy;

    const axom::IndexType labelCount = labels.size();

    axom::Array<axom::IndexType> tmpLabels(labels.shape(), labels.getAllocatorID());
    auto tmpLabelsView = tmpLabels.view();
    axom::for_all<ExecSpace>(
      labelCount,
      AXOM_LAMBDA(axom::IndexType ci) { tmpLabelsView[ci] = labels[ci] == LabelType::LABEL_ON; });

    RAJA::inclusive_scan_inplace<ScanPolicy>(RAJA::make_span(tmpLabels.data(), tmpLabels.size()),
                                             RAJA::operators::plus<axom::IndexType> {});

    axom::IndexType onCount;  // Count of tets labeled ON.
    axom::copy(&onCount, &tmpLabels.back(), sizeof(onCount));

    if(onIndices.size() < onCount || onIndices.getAllocatorID() != labels.getAllocatorID())
    {
      onIndices = axom::Array<axom::IndexType> {axom::ArrayOptions::Uninitialized(),
                                                onCount,
                                                onCount,
                                                labels.getAllocatorID()};
    }
    auto onIndicesView = onIndices.view();

    LabelType firstLabel = LabelType::LABEL_IN;
    axom::copy(&firstLabel, &labels[0], sizeof(firstLabel));
    if(firstLabel == LabelType::LABEL_ON)
    {
      axom::IndexType zero = 0;
      axom::copy(&onIndices[0], &zero, sizeof(zero));
    }

    axom::for_all<ExecSpace>(
      1,
      labelCount,
      AXOM_LAMBDA(axom::IndexType i) {
        if(tmpLabelsView[i] != tmpLabelsView[i - 1])
        {
          onIndicesView[tmpLabelsView[i - 1]] = i;
        }
      });
  }

  void remapTetIndices(axom::ArrayView<const axom::IndexType> cellIndices,
                       axom::ArrayView<axom::IndexType> tetIndices) override
  {
    if(tetIndices.empty())
    {
      return;
    }

    axom::for_all<ExecSpace>(
      tetIndices.size(),
      AXOM_LAMBDA(axom::IndexType i) {
        auto tetIdIn = tetIndices[i];
        auto cellIdFake = tetIdIn / NUM_TETS_PER_HEX;
        auto cellIdTrue = cellIndices[cellIdFake];
        auto tetIdInCell = tetIdIn % NUM_TETS_PER_HEX;
        auto tetIdOut = cellIdTrue * NUM_TETS_PER_HEX + tetIdInCell;
        tetIndices[i] = tetIdOut;
      });
  }

  /*
   * Clip tets from the mesh with tets or octs from the clipping
   * geometry.  This implementation was lifted from IntersectionShaper
   * and modified to work both tet and oct representations of the
   * geometry.
   */
  void computeClipVolumes3D(axom::ArrayView<double> ovlap, conduit::Node& statistics) override
  {
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

    ShapeMesh& shapeMesh = getShapeMesh();
    const int allocId = shapeMesh.getAllocatorID();
    const IndexType cellCount = shapeMesh.getCellCount();

    /*
     * Geometry as discrete tets or octs, and their bounding boxes.
     */
    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    axom::Array<BoundingBoxType> pieceBbs;
    spin::BVH<3, ExecSpace, double> bvh;
    bool useTets = getDiscreteGeometry(geomAsTets, geomAsOcts, pieceBbs, bvh);
    auto geomTetsView = geomAsTets.view();
    auto geomOctsView = geomAsOcts.view();

    /*
     * Find which shape bounding boxes intersect hexahedron bounding boxes
     */

    AXOM_ANNOTATE_BEGIN("MeshClipper:find_candidates");
    axom::Array<IndexType> counts(cellCount, cellCount, allocId);
    axom::Array<IndexType> offsets(cellCount, cellCount, allocId);
    axom::Array<IndexType> candidates;
    bvh.findBoundingBoxes(offsets, counts, candidates, cellCount, shapeMesh.getCellBoundingBoxes());
    AXOM_ANNOTATE_END("MeshClipper:find_candidates");

    const auto countsView = counts.view();
    const int candidateCount = candidates.size();

    AXOM_ANNOTATE_BEGIN("allocate scratch space");
    // Initialize hexahedron indices and shape candidates
    axom::Array<IndexType> hexIndices(candidateCount * NUM_TETS_PER_HEX,
                                      candidateCount * NUM_TETS_PER_HEX,
                                      allocId);
    auto hexIndicesView = hexIndices.view();

    axom::Array<IndexType> shapeCandidates(candidateCount * NUM_TETS_PER_HEX,
                                           candidateCount * NUM_TETS_PER_HEX,
                                           allocId);
    auto shapeCandidatesView = shapeCandidates.view();

    // Tetrahedra from hexes
    auto cellsAsTets = shapeMesh.getCellsAsTets();
    auto meshTetVolumes = getShapeMesh().getTetVolumes();

    // Index into 'tets'
    axom::Array<IndexType> tetIndices(candidateCount * NUM_TETS_PER_HEX,
                                      candidateCount * NUM_TETS_PER_HEX,
                                      allocId);
    auto tetIndicesView = tetIndices.view();
    AXOM_ANNOTATE_END("allocate scratch space");

    // New total number of candidates after omitting degenerate shapes
    AXOM_ANNOTATE_BEGIN("newTotalCandidates memory");
    IndexType tetCandidatesCount = 0;
    IndexType* tetCandidatesCountPtr = &tetCandidatesCount;
    if(!axom::execution_space<ExecSpace>::usesMemorySpace(MemorySpace::Dynamic))
    {
      // Use temporary space compatible with runtime policy.
      tetCandidatesCountPtr = axom::allocate<IndexType>(1, allocId);
      axom::copy(tetCandidatesCountPtr, &tetCandidatesCount, sizeof(tetCandidatesCount));
    }
    AXOM_ANNOTATE_END("newTotalCandidates memory");

    const auto offsetsView = offsets.view();
    const auto candidatesView = candidates.view();
    {
      AXOM_ANNOTATE_SCOPE("init_candidates");
      axom::for_all<ExecSpace>(
        cellCount,
        AXOM_LAMBDA(axom::IndexType i) {
          for(int j = 0; j < countsView[i]; j++)
          {
            int shapeIdx = candidatesView[offsetsView[i] + j];

            for(int k = 0; k < NUM_TETS_PER_HEX; k++)
            {
              IndexType idx = RAJA::atomicAdd<ATOMIC_POL>(tetCandidatesCountPtr, IndexType {1});
              hexIndicesView[idx] = i;
              shapeCandidatesView[idx] = shapeIdx;
              tetIndicesView[idx] = i * NUM_TETS_PER_HEX + k;
            }
          }
        });
    }

    tetCandidatesCount = NUM_TETS_PER_HEX * candidates.size();
#if defined(AXOM_DEBUG)
    // Verifying: this should always pass.
    if(tetCandidatesCountPtr != &tetCandidatesCount)
    {
      axom::copy(&tetCandidatesCount, tetCandidatesCountPtr, sizeof(IndexType));
    }
    SLIC_ASSERT(tetCandidatesCount == candidateCount * NUM_TETS_PER_HEX);
#endif

    SLIC_DEBUG(
      axom::fmt::format("Running clip loop on {} candidate pieces for of all {} hexes in the mesh",
                        tetCandidatesCount,
                        cellCount));

    /*
     * Initialize statistics and copy some to allocId space for computing.
     */

    const std::int64_t zero = 0;
    std::int64_t& clipCount = *(statistics["clips"] = zero).as_int64_ptr();
    std::int64_t& contribCount = *(statistics["contribs"] = zero).as_int64_ptr();
    statistics["candidates"].set_int64(candidateCount);

    std::int64_t* clipCountPtr = axom::allocate<std::int64_t>(1, allocId);
    std::int64_t* contribCountPtr = axom::allocate<std::int64_t>(1, allocId);
    axom::copy(clipCountPtr, &clipCount, sizeof(zero));
    axom::copy(contribCountPtr, &contribCount, sizeof(zero));

    const auto screenLevel = m_myClipper.getScreenLevel();
    constexpr double EPS1 = EPS;

    AXOM_ANNOTATE_BEGIN("MeshClipper:clipLoop_notScreened");
    if(useTets)
    {
      axom::for_all<ExecSpace>(
        tetCandidatesCount,
        AXOM_LAMBDA(axom::IndexType i) {
          const int tetIndex = tetIndicesView[i];

          // Skip degenerate mesh tets.
          // Tet screening can filter out degenerate tets, but this method
          // assumes no tet screening.
          if(axom::utilities::isNearlyEqual(meshTetVolumes[tetIndex], 0.0, 1e-10))
          {
            return;
          }

          const int cellId = hexIndicesView[i];
          const int pieceId = shapeCandidatesView[i];
          const auto& cellTet = cellsAsTets[tetIndex];
          const TetrahedronType& geomPiece = geomTetsView[pieceId];

          double volume = 0.0;
          LabelType tmpLabel =
            computeMeshTetGeomPieceOverlap(cellTet, geomPiece, volume, screenLevel);
          if(tmpLabel == LabelType::LABEL_IN || tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
            RAJA::atomicAdd<ATOMIC_POL>(contribCountPtr, std::int64_t(volume >= EPS1));
          }
          if(tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(clipCountPtr, std::int64_t(1));
          }
        });
    }
    else  // useOcts
    {
      axom::for_all<ExecSpace>(
        tetCandidatesCount,
        AXOM_LAMBDA(axom::IndexType i) {
          const int tetIndex = tetIndicesView[i];

          // Skip degenerate mesh tets.
          if(axom::utilities::isNearlyEqual(meshTetVolumes[tetIndex], 0.0, 1e-10))
          {
            return;
          }

          const int cellId = hexIndicesView[i];
          const auto& cellTet = cellsAsTets[tetIndex];
          const int pieceId = shapeCandidatesView[i];
          const OctahedronType& geomPiece = geomOctsView[pieceId];

          double volume = 0.0;
          LabelType tmpLabel =
            computeMeshTetGeomPieceOverlap(cellTet, geomPiece, volume, screenLevel);
          if(tmpLabel == LabelType::LABEL_IN || tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
            RAJA::atomicAdd<ATOMIC_POL>(contribCountPtr, std::int64_t(volume >= EPS1));
          }
          if(tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(clipCountPtr, std::int64_t(1));
          }
        });
    }
    AXOM_ANNOTATE_END("MeshClipper:clipLoop_notScreened");
    clipCount = tetCandidatesCount;
    axom::copy(&contribCount, contribCountPtr, sizeof(contribCount));
    axom::deallocate(contribCountPtr);

    if(tetCandidatesCountPtr != &tetCandidatesCount)
    {
      axom::deallocate(tetCandidatesCountPtr);
    }
  }  // end of computeClipVolumes3D() function

  /*
   * Clip tets from the mesh with tets or octs from the clipping
   * geometry.  This implementation is like the above except that it
   * limits clipping to a subset of mesh cells labeled as potentially
   * on the boundary.
   */
  void computeClipVolumes3D(const axom::ArrayView<axom::IndexType>& cellIndices,
                            axom::ArrayView<double> ovlap,
                            conduit::Node& statistics) override

  {
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

    ShapeMesh& shapeMesh = getShapeMesh();
    const int allocId = shapeMesh.getAllocatorID();
    const IndexType cellCount = shapeMesh.getCellCount();

    /*
     * Geometry as discrete tets or octs, and their bounding boxes.
     */
    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    axom::Array<BoundingBoxType> pieceBbs;
    spin::BVH<3, ExecSpace, double> bvh;
    bool useTets = getDiscreteGeometry(geomAsTets, geomAsOcts, pieceBbs, bvh);
    auto geomTetsView = geomAsTets.view();
    auto geomOctsView = geomAsOcts.view();

    /*
     * Find which shape bounding boxes intersect hexahedron bounding boxes
     */

    AXOM_ANNOTATE_BEGIN("MeshClipper:find_candidates");
    // Find which shape bounding boxes intersect hexahedron bounding boxes
    // Create a temporary subset of cell bounding boxes,
    // containing only those listed in cellIndices.
    axom::ArrayView<const BoundingBoxType> cellBbsView = shapeMesh.getCellBoundingBoxes();
    const axom::IndexType limitedCellCount = cellIndices.size();
    axom::Array<BoundingBoxType> limitedCellBbs(limitedCellCount, limitedCellCount, allocId);
    axom::ArrayView<BoundingBoxType> limitedCellBbsView = limitedCellBbs.view();
    axom::for_all<ExecSpace>(
      limitedCellBbsView.size(),
      AXOM_LAMBDA(axom::IndexType i) { limitedCellBbsView[i] = cellBbsView[cellIndices[i]]; });

    axom::Array<IndexType> counts(limitedCellCount, limitedCellCount, allocId);
    axom::Array<IndexType> offsets(limitedCellCount, limitedCellCount, allocId);
    axom::Array<IndexType> candidates;
    bvh.findBoundingBoxes(offsets, counts, candidates, limitedCellCount, limitedCellBbsView);
    AXOM_ANNOTATE_END("MeshClipper:find_candidates");

    const auto countsView = counts.view();
    const int candidateCount = candidates.size();

    AXOM_ANNOTATE_BEGIN("allocate scratch space");
    // Initialize hexahedron indices and shape candidates
    axom::Array<IndexType> hexIndices(candidateCount * NUM_TETS_PER_HEX,
                                      candidateCount * NUM_TETS_PER_HEX,
                                      allocId);
    auto hexIndicesView = hexIndices.view();

    axom::Array<IndexType> shapeCandidates(candidateCount * NUM_TETS_PER_HEX,
                                           candidateCount * NUM_TETS_PER_HEX,
                                           allocId);
    auto shapeCandidatesView = shapeCandidates.view();

    // Tetrahedrons from hexes
    auto cellsAsTets = shapeMesh.getCellsAsTets();
    auto meshTetVolumes = getShapeMesh().getTetVolumes();

    // Index into 'tets'
    axom::Array<IndexType> tetIndices(candidateCount * NUM_TETS_PER_HEX,
                                      candidateCount * NUM_TETS_PER_HEX,
                                      allocId);
    auto tetIndicesView = tetIndices.view();
    AXOM_ANNOTATE_END("allocate scratch space");

    // New total number of candidates after omitting degenerate shapes
    AXOM_ANNOTATE_BEGIN("newTotalCandidates memory");
    IndexType tetCandidatesCount = 0;
    IndexType* tetCandidatesCountPtr = &tetCandidatesCount;
    if(!axom::execution_space<ExecSpace>::usesMemorySpace(MemorySpace::Dynamic))
    {
      // Use temporary space compatible with runtime policy.
      tetCandidatesCountPtr = axom::allocate<IndexType>(1, allocId);
      axom::copy(tetCandidatesCountPtr, &tetCandidatesCount, sizeof(IndexType));
    }
    AXOM_ANNOTATE_END("newTotalCandidates memory");

    const auto offsetsView = offsets.view();
    const auto candidatesView = candidates.view();
    {
      AXOM_ANNOTATE_SCOPE("init_candidates");
      axom::for_all<ExecSpace>(
        limitedCellCount,
        AXOM_LAMBDA(axom::IndexType i) {
          for(int j = 0; j < countsView[i]; j++)
          {
            int shapeIdx = candidatesView[offsetsView[i] + j];

            for(int k = 0; k < NUM_TETS_PER_HEX; k++)
            {
              IndexType idx = RAJA::atomicAdd<ATOMIC_POL>(tetCandidatesCountPtr, IndexType {1});
              hexIndicesView[idx] = i;
              shapeCandidatesView[idx] = shapeIdx;
              tetIndicesView[idx] = i * NUM_TETS_PER_HEX + k;
            }
          }
        });
    }

    SLIC_DEBUG(axom::fmt::format(
      "Running clip loop on {} candidate pieces for the select {} hexes of the full {} mesh cells",
      tetCandidatesCount,
      cellIndices.size(),
      cellCount));

    tetCandidatesCount = NUM_TETS_PER_HEX * candidates.size();
#if defined(AXOM_DEBUG)
    // Verifying: this should always pass.
    if(tetCandidatesCountPtr != &tetCandidatesCount)
    {
      axom::copy(&tetCandidatesCount, tetCandidatesCountPtr, sizeof(IndexType));
    }
    SLIC_ASSERT(tetCandidatesCount == candidateCount * NUM_TETS_PER_HEX);
#endif

    /*
     * Initialize statistics and copy some to allocId space for computing.
     */

    const std::int64_t zero = 0;
    std::int64_t& clipCount = *(statistics["clips"] = zero).as_int64_ptr();
    std::int64_t& contribCount = *(statistics["contribs"] = zero).as_int64_ptr();
    statistics["candidates"].set_int64(candidateCount);

    std::int64_t* clipCountPtr = axom::allocate<std::int64_t>(1, allocId);
    std::int64_t* contribCountPtr = axom::allocate<std::int64_t>(1, allocId);
    axom::copy(clipCountPtr, &clipCount, sizeof(zero));
    axom::copy(contribCountPtr, &contribCount, sizeof(zero));

    const auto screenLevel = m_myClipper.getScreenLevel();

    AXOM_ANNOTATE_BEGIN("MeshClipper:clipLoop_hexScreened");
    if(useTets)
    {
      axom::for_all<ExecSpace>(
        tetCandidatesCount,
        AXOM_LAMBDA(axom::IndexType i) {
          int tetIndex =
            tetIndicesView[i];  // index into BVH results, implicit because BVH results specify hexes, not tets.
          int tetIndex1 = tetIndex / NUM_TETS_PER_HEX;
          int tetIndex2 = tetIndex % NUM_TETS_PER_HEX;
          tetIndex = cellIndices[tetIndex1] * NUM_TETS_PER_HEX +
            tetIndex2;  // Now it indexes into the full tets-from-hexes array.

          // Skip degenerate mesh tets.
          // Tet screening can filter out degenerate tets, but this method
          // assumes no tet screening.
          if(axom::utilities::isNearlyEqual(meshTetVolumes[tetIndex], 0.0, 1e-10))
          {
            return;
          }

          int cellId = hexIndicesView[i];  // index into limited mesh hex array
          cellId = cellIndices[cellId];    // Now, it indexes into the full hex array.

          const int pieceId = shapeCandidatesView[i];  // index into pieces array
          const auto& cellTet = cellsAsTets[tetIndex];
          const TetrahedronType& geomPiece = geomTetsView[pieceId];

          double volume = 0.0;
          LabelType tmpLabel =
            computeMeshTetGeomPieceOverlap(cellTet, geomPiece, volume, screenLevel);
          if(tmpLabel == LabelType::LABEL_IN || tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
            RAJA::atomicAdd<ATOMIC_POL>(contribCountPtr, std::int64_t(volume >= EPS));
          }
          if(tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(clipCountPtr, std::int64_t(1));
          }
        });
    }
    else  // useOcts
    {
      axom::for_all<ExecSpace>(
        tetCandidatesCount,
        AXOM_LAMBDA(axom::IndexType i) {
          int tetIndex =
            tetIndicesView[i];  // index into BVH results, implicit because BVH results specify hexes, not tets.
          int tetIndex1 = tetIndex / NUM_TETS_PER_HEX;
          int tetIndex2 = tetIndex % NUM_TETS_PER_HEX;
          tetIndex = cellIndices[tetIndex1] * NUM_TETS_PER_HEX +
            tetIndex2;  // Now it indexes into the full tets-from-hexes array.

          // Skip degenerate mesh tets.
          if(axom::utilities::isNearlyEqual(meshTetVolumes[tetIndex], 0.0, 1e-10))
          {
            return;
          }

          int cellId = hexIndicesView[i];  // index into limited mesh hex array
          cellId = cellIndices[cellId];    // Now, it indexes into the full hex array.

          const int pieceId = shapeCandidatesView[i];  // index into pieces array
          const OctahedronType& geomPiece = geomOctsView[pieceId];
          const auto& cellTet = cellsAsTets[tetIndex];

          double volume = 0.0;
          LabelType tmpLabel =
            computeMeshTetGeomPieceOverlap(cellTet, geomPiece, volume, screenLevel);
          if(tmpLabel == LabelType::LABEL_IN || tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
            RAJA::atomicAdd<ATOMIC_POL>(contribCountPtr, std::int64_t(volume >= EPS));
          }
          if(tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(clipCountPtr, std::int64_t(1));
          }
        });
    }
    AXOM_ANNOTATE_END("MeshClipper:clipLoop_hexScreened");

    clipCount = tetCandidatesCount;
    axom::copy(&contribCount, contribCountPtr, sizeof(contribCount));
    axom::deallocate(contribCountPtr);

    if(tetCandidatesCountPtr != &tetCandidatesCount)
    {
      axom::deallocate(tetCandidatesCountPtr);
    }
  }  // end of computeClipVolumes3D() function

  /*
   * Clip tets of from the mesh with tets or octs from the clipping
   * geometry.  This implementation is like the two above except that
   * it limits clipping to a subset of mesh tets labeled as
   * potentially on the boundary.
   */
  void computeClipVolumes3DTets(const axom::ArrayView<axom::IndexType>& tetIndices,
                                axom::ArrayView<double> ovlap,
                                conduit::Node& statistics) override

  {
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

    ShapeMesh& shapeMesh = getShapeMesh();
    auto meshTets = getShapeMesh().getCellsAsTets();

    const int allocId = shapeMesh.getAllocatorID();

    /*
     * Geometry as discrete tets or octs, and their bounding boxes.
     */
    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    axom::Array<BoundingBoxType> pieceBbs;
    spin::BVH<3, ExecSpace, double> bvh;
    bool useTets = getDiscreteGeometry(geomAsTets, geomAsOcts, pieceBbs, bvh);
    auto geomTetsView = geomAsTets.view();
    auto geomOctsView = geomAsOcts.view();

    /*
     * Find which shape bounding boxes intersect hexahedron bounding boxes
     */

    AXOM_ANNOTATE_BEGIN("MeshClipper:find_candidates");
    // Create a temporary subset of tet bounding boxes,
    // containing only those listed in tetIndices.
    // The BVH searches on this array.
    const axom::IndexType tetCount = tetIndices.size();
    axom::Array<BoundingBoxType> tetBbs(tetCount, tetCount, allocId);
    axom::ArrayView<BoundingBoxType> tetBbsView = tetBbs.view();
    axom::for_all<ExecSpace>(
      tetCount,
      AXOM_LAMBDA(axom::IndexType i) {
        auto& tetBb = tetBbsView[i];
        axom::IndexType tetId = tetIndices[i];
        const auto& tet = meshTets[tetId];
        for(int j = 0; j < 4; ++j) tetBb.addPoint(tet[j]);
      });

    axom::Array<IndexType> counts(tetCount, tetCount, allocId);
    axom::Array<IndexType> offsets(tetCount, tetCount, allocId);
    axom::Array<IndexType> candidates;
    auto countsView = counts.view();
    auto offsetsView = offsets.view();
#if 1
    // Get the BVH traverser for doing the 2-pass search manually.
    const auto bvhTraverser = bvh.getTraverser();
    /*
     * Predicate for traversing the BVH.  We enter BVH nodes
     * whose bounding boxes intersect the query bounding box.
     */
    auto traversePredicate = [] AXOM_HOST_DEVICE(const BoundingBoxType& queryBbox,
                                                 const BoundingBoxType& bvhBbox) -> bool {
                               // Could replace this with mesh tet as first arg.
                               // Should instantiate by template.  called at non-leaves.
                               return queryBbox.intersectsWith(bvhBbox);
                             };

    /*
     * First pass: count number of collisions each of the tetBbs makes
     * with the BVH leaves.  Populate the counts array.
     */
    axom::ReduceSum<ExecSpace, IndexType> totalCountReduce(0);
    axom::for_all<ExecSpace>(
      tetCount,
      AXOM_LAMBDA(axom::IndexType iTet) {
        axom::IndexType count = 0;
        auto countCollisions = [&](std::int32_t currentNode, const std::int32_t* leafNodes) {
                                 // countCollisions is only called at the leaves.
                                 AXOM_UNUSED_VAR(currentNode);
                                 AXOM_UNUSED_VAR(leafNodes);
                                 ++count;
                                 // Eventually, bypass if tet and candidate can be shown not to collide.
                                 };
        bvhTraverser.traverse_tree(tetBbsView[iTet], countCollisions, traversePredicate);
        countsView[iTet] = count;
        totalCountReduce += count;
      });

    // Compute the offsets array using a prefix scan of counts.
    axom::exclusive_scan<ExecSpace>(counts, offsets);
    const IndexType nCollisions = totalCountReduce.get();

    /*
     * Allocate 2 arrays to hold info about the meshTet/geometry piece collisions.
     * - candidates: geometry pieces in a potential collision, actually their indices.
     * - candToTetIdId: indicates the meshTets in the collision,
     *   where candToTetIdId[i] corresponds to meshTets[tetIndices[i]].
     */
    candidates = axom::Array<IndexType>(nCollisions, nCollisions, allocId);
    axom::Array<IndexType> candToTetIdId(candidates.size(), candidates.size(), allocId);
    auto candidatesView = candidates.view();
    auto candToTetIdIdView = candToTetIdId.view();

    /*
     * Second pass: Populate tet-candidate piece collision arrays.
     */
    axom::for_all<ExecSpace>(
      tetCount,
      AXOM_LAMBDA(axom::IndexType iTet) {
        auto offset = offsetsView[iTet];

        // PrimitiveType cellTet = tetBbsView[iTet]; // Eventually, use the tet, not its BB.
        // Record indices of the tet and the candidate that collided.
        // Eventually, bypass if tet and candidate can be shown not to collide.
        auto recordCollision = [&](std::int32_t currentNode, const std::int32_t* leafs) {
                                 candToTetIdIdView[offset] = iTet;
                                 candidatesView[offset] = leafs[currentNode];
                                 ++offset;
                               };

        bvhTraverser.traverse_tree(tetBbsView[iTet], recordCollision, traversePredicate);
      });
#else
    bvh.findBoundingBoxes(offsets, counts, candidates, tetBbsView.size(), tetBbsView);
    auto candidatesView = candidates.view();

    axom::Array<IndexType> candToTetIdId(candidates.size(), candidates.size(), allocId);
    auto candToTetIdIdView = candToTetIdId.view();
    axom::for_all<ExecSpace>(
      tetCount,
      AXOM_LAMBDA(axom::IndexType i) {
        auto count = countsView[i];
        auto offset = offsetsView[i];
        for(int j = 0; j < count; ++j) candToTetIdIdView[offset + j] = i;
      });
#endif
    AXOM_ANNOTATE_END("MeshClipper:find_candidates");

    SLIC_DEBUG(axom::fmt::format(
      "Running clip loop on {} candidate pieces for the select {} tets of the full {} mesh cells",
      candidates.size(),
      tetCount,
      shapeMesh.getCellCount()));

    /*
     * Initialize statistics and copy some to allocId space for computing.
     */

    const std::int64_t zero = 0;
    std::int64_t& clipCount = *(statistics["clips"] = zero).as_int64_ptr();
    std::int64_t& contribCount = *(statistics["contribs"] = zero).as_int64_ptr();
    std::int64_t& candidateCount = *(statistics["candidates"] = zero).as_int64_ptr();
    candidateCount = candidates.size();

    std::int64_t* clipCountPtr = axom::allocate<std::int64_t>(1, allocId);
    std::int64_t* contribCountPtr = axom::allocate<std::int64_t>(1, allocId);
    axom::copy(clipCountPtr, &clipCount, sizeof(zero));
    axom::copy(contribCountPtr, &contribCount, sizeof(zero));

    const auto screenLevel = m_myClipper.getScreenLevel();

    AXOM_ANNOTATE_BEGIN("MeshClipper:clipLoop_tetScreened");
    if(useTets)
    {
      axom::for_all<ExecSpace>(
        candidates.size(),
        AXOM_LAMBDA(axom::IndexType iCand) {
          auto tetIdId = candToTetIdIdView[iCand];
          auto tetId = tetIndices[tetIdId];

          auto cellId = tetId / NUM_TETS_PER_HEX;
          auto pieceId = candidatesView[iCand];
          const auto& meshTet = meshTets[tetId];
          const TetrahedronType& geomPiece = geomTetsView[pieceId];

          double volume = 0.0;
          LabelType tmpLabel =
            computeMeshTetGeomPieceOverlap(meshTet, geomPiece, volume, screenLevel);
          if(tmpLabel == LabelType::LABEL_IN || tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
            RAJA::atomicAdd<ATOMIC_POL>(contribCountPtr, std::int64_t(volume >= EPS));
          }
          if(tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(clipCountPtr, std::int64_t(1));
          }
        });
    }
    else  // useOcts
    {
      axom::for_all<ExecSpace>(
        candidates.size(),
        AXOM_LAMBDA(axom::IndexType iCand) {
          auto tetIdId = candToTetIdIdView[iCand];
          auto tetId = tetIndices[tetIdId];

          auto cellId = tetId / NUM_TETS_PER_HEX;
          auto pieceId = candidatesView[iCand];
          const auto& meshTet = meshTets[tetId];
          const OctahedronType& geomPiece = geomOctsView[pieceId];

          double volume = 0.0;
          LabelType tmpLabel =
            computeMeshTetGeomPieceOverlap(meshTet, geomPiece, volume, screenLevel);
          if(tmpLabel == LabelType::LABEL_IN || tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
            RAJA::atomicAdd<ATOMIC_POL>(contribCountPtr, std::int64_t(volume >= EPS));
          }
          if(tmpLabel == LabelType::LABEL_ON)
          {
            RAJA::atomicAdd<ATOMIC_POL>(clipCountPtr, std::int64_t(1));
          }
        });
    }
    AXOM_ANNOTATE_END("MeshClipper:clipLoop_tetScreened");

    axom::copy(&contribCount, contribCountPtr, sizeof(contribCount));
    axom::copy(&clipCount, clipCountPtr, sizeof(clipCount));
    axom::deallocate(contribCountPtr);
    axom::deallocate(clipCountPtr);

    SLIC_DEBUG(axom::fmt::format(""));
  }  // end of computeClipVolumes3DTets() function

  /*!
   * @brief Get the geometry in discrete pieces,
   *   which can be tets or octs, and place them in a search tree.
   * @return whether geometry are tetrahedra instead of octahedra.
   */
  bool getDiscreteGeometry(axom::Array<axom::primal::Tetrahedron<double, 3>>& geomAsTets,
                           axom::Array<axom::primal::Octahedron<double, 3>>& geomAsOcts,
                           axom::Array<BoundingBoxType>& pieceBbs,
                           spin::BVH<3, ExecSpace, double>& bvh)
  {
    auto& strategy = getStrategy();
    ShapeMesh& shapeMesh = getShapeMesh();
    int allocId = shapeMesh.getAllocatorID();

    AXOM_ANNOTATE_BEGIN("MeshClipper:get_geometry");
    const bool useOcts = strategy.getGeometryAsOcts(shapeMesh, geomAsOcts);
    const bool useTets = strategy.getGeometryAsTets(shapeMesh, geomAsTets);
    AXOM_ANNOTATE_END("MeshClipper:get_geometry");

    if(useTets)
    {
      SLIC_ASSERT(geomAsTets.getAllocatorID() == allocId);
    }
    else
    {
      SLIC_ASSERT(geomAsOcts.getAllocatorID() == allocId);
    }
    if(useTets == useOcts)
    {
      SLIC_ERROR(
        axom::fmt::format("Problem with MeshClipperStrategy implementation '{}'."
                          "  Implementations that don't provide a specializedClip function"
                          " must provide exactly one of either getGeometryAsOcts() or"
                          " getGeometryAsTets().   This implementation provides {}.",
                          strategy.name(),
                          int(useOcts) + int(useTets)));
    }

    SLIC_DEBUG(axom::fmt::format("Geometry {} has {} discrete {}s",
                                 strategy.name(),
                                 useTets ? geomAsTets.size() : geomAsOcts.size(),
                                 useTets ? "tet" : "oct"));

    /*
     * Get the bounding boxes for the discrete geometry pieces.
     * If debug build, check for degenerate pieces.
     */
    const axom::IndexType bbCount = useTets ? geomAsTets.size() : geomAsOcts.size();
    pieceBbs = axom::Array<BoundingBoxType>(bbCount, bbCount, allocId);
    axom::ArrayView<BoundingBoxType> pieceBbsView = pieceBbs.view();

    if(useTets)
    {
      auto geomTetsView = geomAsTets.view();
      axom::for_all<ExecSpace>(
        pieceBbsView.size(),
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomTetsView[i]);
#if defined(AXOM_DEBUG)
          SLIC_ASSERT(!geomTetsView[i].degenerate());
#endif
        });
    }
    else
    {
      auto geomOctsView = geomAsOcts.view();
      axom::for_all<ExecSpace>(
        pieceBbsView.size(),
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomOctsView[i]);
        });
    }

    bvh.setAllocatorID(allocId);
    bvh.setTolerance(EPS);
    bvh.setScaleFactor(BVH_SCALE_FACTOR);
    bvh.initialize(pieceBbsView, pieceBbsView.size());

    return useTets;
  }

  /*!
   * @brief Volume of a tetrahedron from discretized geometry.
   */
  AXOM_HOST_DEVICE inline double geomPieceVolume(const TetrahedronType& tet)
  {
    return tet.volume();
  }

  /*!
   * @brief Volume of a octahedron from discretized geometry.
   *
   * Assumes octahedron is convex.
   */
  AXOM_HOST_DEVICE inline double geomPieceVolume(const OctahedronType& oct)
  {
    TetrahedronType tets[] = {TetrahedronType(oct[0], oct[3], oct[1], oct[2]),
                              TetrahedronType(oct[0], oct[3], oct[2], oct[4]),
                              TetrahedronType(oct[0], oct[3], oct[4], oct[5]),
                              TetrahedronType(oct[0], oct[3], oct[5], oct[1])};
    double octVol = 0.0;
    for(int i = 0; i < 4; ++i)
    {
      double tetVol = tets[i].volume();
      SLIC_ASSERT(tetVol >= -EPS);  // Tet may be degenerate but not inverted.
      octVol += tetVol;
    }
    return octVol;
  }

  /*!
   * @brief Compute overlap volume between a tet (from the shape mesh)
   * and a piece (tet or oct) of the discretized geometry.
   *
   * Becaue primal::clip is so expensive, we do a conservative
   * overlap check on @c meshTet and @c geomPiece to avoid clipping.
   *
   * @return results of check whether the piece is IN/ON/OUT of the tet.
   *
   * @tparam TetOrOctType either a TetrahedronType or OctahedronType,
   * the two types a geometry can be discretized into.
   */
  template <typename TetOrOctType>
  AXOM_HOST_DEVICE inline LabelType computeMeshTetGeomPieceOverlap(const TetrahedronType& meshTet,
                                                                   const TetOrOctType& geomPiece,
                                                                   double& overlapVolume,
                                                                   int screenLevel)
  {
    constexpr bool tryFixOrientation = false;
    if(screenLevel >= 3)
    {
      LabelType geomLabel = labelPieceInOutOfTet(meshTet, geomPiece);
      if(geomLabel == LabelType::LABEL_OUT)
      {
        overlapVolume = 0.0;
        return geomLabel;
      }
      if(geomLabel == LabelType::LABEL_IN)
      {
        overlapVolume = geomPieceVolume(geomPiece);
        return geomLabel;
      }
    }

    const auto poly = primal::clip<double>(meshTet, geomPiece, EPS, tryFixOrientation);
    if(poly.numVertices() >= 4)
    {
      // Poly is valid
      overlapVolume = poly.volume();
      SLIC_ASSERT(overlapVolume >= 0);
    }
    else
    {
      overlapVolume = 0.0;
    }

    return LabelType::LABEL_ON;
  }

  /*!
   * @brief Compute whether a tetrahedron or octhedron is inside,
   * outside or on the boundary of a reference tetrahedron,
   * and conservatively label it as on, if not known.
   *
   * @internal To reduce repeatedly computing toUnitTet for
   * the same tet, precompute that in the calling function
   * and use it instead of tet.
   */
  template <typename TetOrOctType>
  AXOM_HOST_DEVICE inline LabelType labelPieceInOutOfTet(const TetrahedronType& tet,
                                                         const TetOrOctType& piece)
  {
    Point3DType unitTet[] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    axom::primal::experimental::CoordinateTransformer toUnitTet(&tet[0], unitTet);

    /*
     * Count (transformed) piece vertices above/below unitTet as unitTet
     * rests on its 4 sides.  Sides 0-2 are perpendicular to the axes.
     * Side 3 is the diagonal side.
     */
    int vsAbove[4] = {0, 0, 0, 0};
    int vsBelow[4] = {0, 0, 0, 0};
    int vsTetSide[4] = {0, 0, 0, 0};
    for(int i = 0; i < TetOrOctType::NUM_VERTS; ++i)
    {
      auto pVert = toUnitTet.getTransformed(piece[i]);
      // h4 is height of pVert above the diagonal face, scaled by sqrt(3).
      // h4 of 1 coresponds to the unitTet's height of sqrt(3).
      double h4 = 1 - (pVert[0] + pVert[1] + pVert[2]);
      vsAbove[0] += pVert[0] >= 1;
      vsAbove[1] += pVert[1] >= 1;
      vsAbove[2] += pVert[2] >= 1;
      vsAbove[3] += h4 >= 1;
      vsBelow[0] += pVert[0] <= 0;
      vsBelow[1] += pVert[1] <= 0;
      vsBelow[2] += pVert[2] <= 0;
      vsBelow[3] += h4 <= 0;
      vsTetSide[0] += pVert[0] >= 0;
      vsTetSide[1] += pVert[1] >= 0;
      vsTetSide[2] += pVert[2] >= 0;
      vsTetSide[3] += h4 >= 0;
    }
    if(vsAbove[0] == TetOrOctType::NUM_VERTS || vsAbove[1] == TetOrOctType::NUM_VERTS ||
       vsAbove[2] == TetOrOctType::NUM_VERTS || vsAbove[3] == TetOrOctType::NUM_VERTS ||
       vsBelow[0] == TetOrOctType::NUM_VERTS || vsBelow[1] == TetOrOctType::NUM_VERTS ||
       vsBelow[2] == TetOrOctType::NUM_VERTS || vsBelow[3] == TetOrOctType::NUM_VERTS)
    {
      return LabelType::LABEL_OUT;
    }
    if(vsTetSide[0] == TetOrOctType::NUM_VERTS && vsTetSide[1] == TetOrOctType::NUM_VERTS &&
       vsTetSide[2] == TetOrOctType::NUM_VERTS && vsTetSide[3] == TetOrOctType::NUM_VERTS)
    {
      return LabelType::LABEL_IN;
    }
    return LabelType::LABEL_ON;
  }

  void getLabelCounts(axom::ArrayView<const LabelType> labels,
                      std::int64_t& inCount,
                      std::int64_t& onCount,
                      std::int64_t& outCount) override
  {
    AXOM_ANNOTATE_SCOPE("MeshClipper:count_labels");
    using ReducePolicy = typename axom::execution_space<ExecSpace>::reduce_policy;
    using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
    RAJA::ReduceSum<ReducePolicy, std::int64_t> inSum(0);
    RAJA::ReduceSum<ReducePolicy, std::int64_t> onSum(0);
    RAJA::ReduceSum<ReducePolicy, std::int64_t> outSum(0);
    RAJA::forall<LoopPolicy>(
      RAJA::RangeSegment(0, labels.size()),
      AXOM_LAMBDA(axom::IndexType cellId) {
        const auto& label = labels[cellId];
        if(label == LabelType::LABEL_OUT)
        {
          outSum += 1;
        }
        else if(label == LabelType::LABEL_IN)
        {
          inSum += 1;
        }
        else
        {
          onSum += 1;
        }
      });
    inCount = static_cast<std::int64_t>(inSum.get());
    onCount = static_cast<std::int64_t>(onSum.get());
    outCount = static_cast<std::int64_t>(outSum.get());
  }

private:
  static constexpr double EPS = 1e-10;
  static constexpr double BVH_SCALE_FACTOR = 1.0;
  static constexpr int MAX_VERTS_FOR_TET_CLIPPING = 32;
  static constexpr int MAX_NBRS_PER_VERT_FOR_TET_CLIPPING = 8;
  static constexpr int MAX_VERTS_FOR_OCT_CLIPPING = 32;
  static constexpr int MAX_NBRS_PER_VERT_FOR_OCT_CLIPPING = 8;
};

}  // end namespace detail
}  // namespace experimental
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_MESHCLIPPERIMPL_HPP_
