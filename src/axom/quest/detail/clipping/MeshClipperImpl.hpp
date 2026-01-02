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

#include "axom/quest/MeshClipperStrategy.hpp"
#include "axom/quest/MeshClipper.hpp"
#include "axom/spin/BVH.hpp"
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

    AXOM_ANNOTATE_SCOPE("MeshClipper::collect_unlabeleds");
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

  // Work space for clip counters.
  struct ClippingStats
  {
    axom::ReduceSum<ExecSpace, IndexType> inSum {0};
    axom::ReduceSum<ExecSpace, IndexType> onSum {0};
    axom::ReduceSum<ExecSpace, IndexType> outSum {0};
    axom::ReduceSum<ExecSpace, IndexType> missSum {0};
    ClippingStats() : inSum(0), onSum(0), outSum(0), missSum(0) { }
    void copyTo(conduit::Node& stats)
    {
      // Place clip counts in statistics container.
      std::int64_t clipsInCount = inSum.get();
      std::int64_t clipsOnCount = onSum.get();
      std::int64_t clipsOutCount = outSum.get();
      std::int64_t clipsMissCount = missSum.get();
      stats["clipsIn"].set_int64(clipsInCount);
      stats["clipsOn"].set_int64(clipsOnCount);
      stats["clipsOut"].set_int64(clipsOutCount);
      stats["clipsMiss"].set_int64(clipsMissCount);
      stats["clipsSum"] = clipsInCount + clipsOnCount + clipsOutCount;
    }
  };

  /*
   * Clip tets from the mesh with tets or octs from the clipping
   * geometry.  This implementation was lifted from IntersectionShaper
   * and modified to work both tet and oct representations of the
   * geometry.
   */
  void computeClipVolumes3D(axom::ArrayView<double> ovlap, conduit::Node& statistics) override
  {
    AXOM_ANNOTATE_SCOPE("MeshClipper::computeClipVolumes3D");

    using BoundingBoxType = primal::BoundingBox<double, 3>;

    ShapeMesh& shapeMesh = getShapeMesh();

    const int allocId = shapeMesh.getAllocatorID();

    const IndexType cellCount = shapeMesh.getCellCount();

    SLIC_INFO(axom::fmt::format(
      "MeshClipper::computeClipVolumes3D: Getting discrete geometry for shape '{}'",
      getStrategy().name()));

    //
    // Get the geometry in discrete pieces, which can be tets or octs.
    //
    auto& strategy = getStrategy();
    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    const bool useOcts = strategy.getGeometryAsOcts(shapeMesh, geomAsOcts);
    const bool useTets = strategy.getGeometryAsTets(shapeMesh, geomAsTets);
    SLIC_ASSERT(useOcts || geomAsOcts.empty());
    SLIC_ASSERT(useTets || geomAsTets.empty());
    if(useTets == useOcts)
    {
      SLIC_ERROR(
        axom::fmt::format("Problem with MeshClipperStrategy implementation '{}'."
                          "  Implementations that don't provide a specializedClip function"
                          " must provide exactly one getGeometryAsOcts() or getGeometryAsTets()."
                          "  This implementation provides {}.",
                          strategy.name(),
                          int(useOcts) + int(useTets)));
    }

    auto geomTetsView = geomAsTets.view();
    auto geomOctsView = geomAsOcts.view();

    SLIC_INFO(axom::fmt::format("{:-^80}", " Inserting shapes' bounding boxes into BVH "));

    //
    // Generate the BVH over the (bounding boxes of the) discretized geometry
    //
    const axom::IndexType bbCount = useTets ? geomAsTets.size() : geomAsOcts.size();
    axom::Array<BoundingBoxType> pieceBbs(bbCount, bbCount, allocId);
    axom::ArrayView<BoundingBoxType> pieceBbsView = pieceBbs.view();

    // Get the bounding boxes for the shapes
    if(useTets)
    {
      axom::for_all<ExecSpace>(
        pieceBbsView.size(),
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomTetsView[i]);
        });
    }
    else
    {
      axom::for_all<ExecSpace>(
        pieceBbsView.size(),
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomOctsView[i]);
        });
    }

    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(pieceBbsView, pieceBbsView.size());

    SLIC_INFO(axom::fmt::format("{:-^80}", " Querying the BVH tree "));

    axom::ArrayView<const BoundingBoxType> cellBbsView = shapeMesh.getCellBoundingBoxes();

    // Find which shape bounding boxes intersect hexahedron bounding boxes
    SLIC_INFO(
      axom::fmt::format("{:-^80}", " Finding shape candidates for each hexahedral element "));

    axom::Array<IndexType> offsets(cellCount, cellCount, allocId);
    axom::Array<IndexType> counts(cellCount, cellCount, allocId);
    axom::Array<IndexType> candidates;
    AXOM_ANNOTATE_BEGIN("bvh.findBoundingBoxes");
    bvh.findBoundingBoxes(offsets, counts, candidates, cellCount, cellBbsView);
    AXOM_ANNOTATE_END("bvh.findBoundingBoxes");

    // Get the total number of candidates
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

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

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = false;

    /*
     * Statistics from the clip loop.
     * Count number of times the piece was found inside/outside/on a mesh tet boundary.
     * Be sure to use kernel-compatible memory.
     */
    ClippingStats clipStats;

    {
      tetCandidatesCount = NUM_TETS_PER_HEX * candidates.size();
      AXOM_ANNOTATE_SCOPE("MeshClipper::clipLoop");
#if defined(AXOM_DEBUG)
      // Verifying: this should always pass.
      if(tetCandidatesCountPtr != &tetCandidatesCount)
      {
        axom::copy(&tetCandidatesCount, tetCandidatesCountPtr, sizeof(IndexType));
      }
      SLIC_ASSERT(tetCandidatesCount == candidateCount * NUM_TETS_PER_HEX);
#endif

      SLIC_INFO(
        axom::fmt::format("Running clip loop on {} candidate tets for of all {} hexes in the mesh",
                          tetCandidatesCount,
                          cellCount));

      if(useTets)
      {
        axom::for_all<ExecSpace>(
          tetCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            const int index = hexIndicesView[i];
            const int shapeIndex = shapeCandidatesView[i];
            const int tetIndex = tetIndicesView[i];
            if(cellsAsTets[tetIndex].degenerate())
            {
              return;
            }

            const auto poly = primal::clip<double>(geomTetsView[shapeIndex],
                                                   cellsAsTets[tetIndex],
                                                   EPS,
                                                   tryFixOrientation);

            // Poly is valid
            if(poly.numVertices() >= 4)
            {
              // Workaround - intermediate volume variable needed for
              // CUDA Pro/E test case correctness
              double volume = poly.volume();
              SLIC_ASSERT(volume >= 0);
              RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + index, volume);
            }
          });
      }
      else  // useOcts
      {
        axom::for_all<ExecSpace>(
          tetCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            const int index = hexIndicesView[i];
            const int shapeIndex = shapeCandidatesView[i];
            const int tetIndex = tetIndicesView[i];
            if(cellsAsTets[tetIndex].degenerate())
            {
              return;
            }

            const auto poly = primal::clip<double>(geomOctsView[shapeIndex],
                                                   cellsAsTets[tetIndex],
                                                   EPS,
                                                   tryFixOrientation);

            // Poly is valid
            if(poly.numVertices() >= 4)
            {
              // Workaround - intermediate volume variable needed for
              // CUDA Pro/E test case correctness
              double volume = poly.volume();
              SLIC_ASSERT(volume >= 0);
              RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + index, volume);
            }
          });
      }
    }

    AXOM_ANNOTATE_END("MeshClipper:clipLoop_notScreened");

    clipStats.copyTo(statistics);
    statistics["clipsCandidates"].set_int64(tetCandidatesCount);

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
    AXOM_ANNOTATE_SCOPE("MeshClipper::computeClipVolumes3D:limited");

    using BoundingBoxType = primal::BoundingBox<double, 3>;

    ShapeMesh& shapeMesh = getShapeMesh();

    const int allocId = shapeMesh.getAllocatorID();

    const IndexType cellCount = shapeMesh.getCellCount();

    SLIC_INFO(axom::fmt::format(
      "MeshClipper::computeClipVolumes3D: Getting discrete geometry for shape '{}'",
      getStrategy().name()));

    auto& strategy = getStrategy();
    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    const bool useOcts = strategy.getGeometryAsOcts(shapeMesh, geomAsOcts);
    const bool useTets = strategy.getGeometryAsTets(shapeMesh, geomAsTets);
    SLIC_ASSERT(useOcts || geomAsOcts.empty());
    SLIC_ASSERT(useTets || geomAsTets.empty());
    if(useTets == useOcts)
    {
      SLIC_ERROR(
        axom::fmt::format("Problem with MeshClipperStrategy implementation '{}'."
                          "  Implementations that don't provide a specializedClip function"
                          " must provide exactly one getGeometryAsOcts() or getGeometryAsTets()."
                          "  This implementation provides {}.",
                          strategy.name(),
                          int(useOcts) + int(useTets)));
    }

    auto geomTetsView = geomAsTets.view();
    auto geomOctsView = geomAsOcts.view();

    SLIC_INFO(axom::fmt::format("{:-^80}", " Inserting shapes' bounding boxes into BVH "));

    // Generate the BVH tree over the shape's discretized geometry
    // axis-aligned bounding boxes.  "pieces" refers to tets or octs.
    const axom::IndexType bbCount = useTets ? geomAsTets.size() : geomAsOcts.size();
    axom::Array<BoundingBoxType> pieceBbs(bbCount, bbCount, allocId);
    axom::ArrayView<BoundingBoxType> pieceBbsView = pieceBbs.view();

    // Get the bounding boxes for the shapes
    if(useTets)
    {
      axom::for_all<ExecSpace>(
        pieceBbsView.size(),
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomTetsView[i]);
        });
    }
    else
    {
      axom::for_all<ExecSpace>(
        pieceBbsView.size(),
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomOctsView[i]);
        });
    }

    // Insert shapes' Bounding Boxes into BVH.
    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(pieceBbsView, pieceBbsView.size());

    SLIC_INFO(axom::fmt::format("{:-^80}", " Querying the BVH tree "));

    // Create a temporary subset of cell bounding boxes,
    // containing only those listed in cellIndices.
    const axom::IndexType limitedCellCount = cellIndices.size();
    axom::ArrayView<const BoundingBoxType> cellBbsView = shapeMesh.getCellBoundingBoxes();
    axom::Array<BoundingBoxType> limitedCellBbs(limitedCellCount, limitedCellCount, allocId);
    axom::ArrayView<BoundingBoxType> limitedCellBbsView = limitedCellBbs.view();
    axom::for_all<ExecSpace>(
      limitedCellBbsView.size(),
      AXOM_LAMBDA(axom::IndexType i) { limitedCellBbsView[i] = cellBbsView[cellIndices[i]]; });

    // Find which shape bounding boxes intersect hexahedron bounding boxes
    SLIC_INFO(
      axom::fmt::format("{:-^80}", " Finding shape candidates for each hexahedral element "));

    axom::Array<IndexType> offsets(limitedCellCount, limitedCellCount, allocId);
    axom::Array<IndexType> counts(limitedCellCount, limitedCellCount, allocId);
    axom::Array<IndexType> candidates;
    AXOM_ANNOTATE_BEGIN("bvh.findBoundingBoxes");
    bvh.findBoundingBoxes(offsets, counts, candidates, limitedCellCount, limitedCellBbsView);
    AXOM_ANNOTATE_END("bvh.findBoundingBoxes");

    // Get the total number of candidates
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

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

    SLIC_INFO(axom::fmt::format(
      "Running clip loop on {} candidate tets for the select {} hexes of the full {} cells",
      tetCandidatesCount,
      cellIndices.size(),
      cellCount));

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = false;

    {
      tetCandidatesCount = NUM_TETS_PER_HEX * candidates.size();
      AXOM_ANNOTATE_SCOPE("MeshClipper::clipLoop_limited");
#if defined(AXOM_DEBUG)
      // Verifying: this should always pass.
      if(tetCandidatesCountPtr != &tetCandidatesCount)
      {
        axom::copy(&tetCandidatesCount, tetCandidatesCountPtr, sizeof(IndexType));
      }
      SLIC_ASSERT(tetCandidatesCount == candidateCount * NUM_TETS_PER_HEX);
#endif

      ClippingStats clipStats;

      if(useTets)
      {
        axom::for_all<ExecSpace>(
          tetCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            int index = hexIndicesView[i];  // index into limited mesh hex array
            index = cellIndices[index];     // Now, it indexes into the full hex array.

            const int shapeIndex = shapeCandidatesView[i];  // index into pieces array
            int tetIndex =
              tetIndicesView[i];  // index into BVH results, implicit because BVH results specify hexes, not tets.
            int tetIndex1 = tetIndex / NUM_TETS_PER_HEX;
            int tetIndex2 = tetIndex % NUM_TETS_PER_HEX;
            tetIndex = cellIndices[tetIndex1] * NUM_TETS_PER_HEX +
              tetIndex2;  // Now it indexes into the full tets-from-hexes array.
            if(cellsAsTets[tetIndex].degenerate())
            {
              return;
            }

            const auto poly = primal::clip<double>(geomTetsView[shapeIndex],
                                                   cellsAsTets[tetIndex],
                                                   EPS,
                                                   tryFixOrientation);

            // Poly is valid
            if(poly.numVertices() >= 4)
            {
              // Workaround - intermediate volume variable needed for
              // CUDA Pro/E test case correctness
              double volume = poly.volume();
              SLIC_ASSERT(volume >= 0);
              RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + index, volume);
            }
          });
      }
      else  // useOcts
      {
        axom::for_all<ExecSpace>(
          tetCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            int index = hexIndicesView[i];  // index into limited mesh hex array
            index = cellIndices[index];     // Now, it indexes into the full hex array.

            const int shapeIndex = shapeCandidatesView[i];  // index into pieces array
            int tetIndex =
              tetIndicesView[i];  // index into BVH results, implicit because BVH results specify hexes, not tets.
            int tetIndex1 = tetIndex / NUM_TETS_PER_HEX;
            int tetIndex2 = tetIndex % NUM_TETS_PER_HEX;
            tetIndex = cellIndices[tetIndex1] * NUM_TETS_PER_HEX +
              tetIndex2;  // Now it indexes into the full tets-from-hexes array.
            if(cellsAsTets[tetIndex].degenerate())
            {
              return;
            }

            const auto poly = primal::clip<double>(geomOctsView[shapeIndex],
                                                   cellsAsTets[tetIndex],
                                                   EPS,
                                                   tryFixOrientation);

            // Poly is valid
            if(poly.numVertices() >= 4)
            {
              // Workaround - intermediate volume variable needed for
              // CUDA Pro/E test case correctness
              double volume = poly.volume();
              SLIC_ASSERT(volume >= 0);
              RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + index, volume);
            }
          });
      }

      clipStats.copyTo(statistics);
      statistics["clipsCandidates"].set_int64(tetCandidatesCount);
    }

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
    AXOM_ANNOTATE_SCOPE("MeshClipper::computeClipVolumes3D:limited");

    using BoundingBoxType = primal::BoundingBox<double, 3>;
    using TetrahedronType = primal::Tetrahedron<double, 3>;
    using OctahedronType = primal::Octahedron<double, 3>;

    ShapeMesh& shapeMesh = getShapeMesh();
    auto meshTets = getShapeMesh().getCellsAsTets();

    const int allocId = shapeMesh.getAllocatorID();

    SLIC_INFO(axom::fmt::format(
      "MeshClipper::computeClipVolumes3D: Getting discrete geometry for shape '{}'",
      getStrategy().name()));

    auto& strategy = getStrategy();
    axom::Array<TetrahedronType> geomAsTets;
    axom::Array<OctahedronType> geomAsOcts;
    const bool useOcts = strategy.getGeometryAsOcts(shapeMesh, geomAsOcts);
    const bool useTets = strategy.getGeometryAsTets(shapeMesh, geomAsTets);
    SLIC_ASSERT(useOcts || geomAsOcts.empty());
    SLIC_ASSERT(useTets || geomAsTets.empty());
    if(useTets == useOcts)
    {
      SLIC_ERROR(
        axom::fmt::format("Problem with MeshClipperStrategy implementation '{}'."
                          "  Implementations that don't provide a specializedClip function"
                          " must provide exactly one getGeometryAsOcts() or getGeometryAsTets()."
                          "  This implementation provides {}.",
                          strategy.name(),
                          int(useOcts) + int(useTets)));
    }

    auto geomTetsView = geomAsTets.view();
    auto geomOctsView = geomAsOcts.view();

    SLIC_INFO(axom::fmt::format("{:-^80}", " Inserting shapes' bounding boxes into BVH "));

    // Generate the BVH tree over the shape's discretized geometry
    // axis-aligned bounding boxes.  "pieces" refers to tets or octs.
    const axom::IndexType bbCount = useTets ? geomAsTets.size() : geomAsOcts.size();
    axom::Array<BoundingBoxType> pieceBbs(bbCount, bbCount, allocId);
    axom::ArrayView<BoundingBoxType> pieceBbsView = pieceBbs.view();

    // Get the bounding boxes for the shapes
    if(useTets)
    {
      axom::for_all<ExecSpace>(
        pieceBbsView.size(),
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomTetsView[i]);
        });
    }
    else
    {
      axom::for_all<ExecSpace>(
        pieceBbsView.size(),
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomOctsView[i]);
        });
    }

    // Insert shapes' Bounding Boxes into BVH.
    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(pieceBbsView, pieceBbsView.size());

    SLIC_INFO(axom::fmt::format("{:-^80}", " Querying the BVH tree "));

    // Create a temporary subset of tet bounding boxes,
    // containing only those listed in tetIndices.
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

    auto countsView = counts.view();
    auto offsetsView = offsets.view();

    AXOM_ANNOTATE_BEGIN("bvh.findBoundingBoxes");
    axom::Array<IndexType> candidates;
    bvh.findBoundingBoxes(offsets, counts, candidates, tetBbsView.size(), tetBbsView);
    auto candidatesView = candidates.view();
    AXOM_ANNOTATE_END("bvh.findBoundingBoxes");

    axom::Array<IndexType> candToTetIdId(candidates.size(), candidates.size(), allocId);
    auto candToTetIdIdView = candToTetIdId.view();
    axom::for_all<ExecSpace>(
      tetCount,
      AXOM_LAMBDA(axom::IndexType i) {
        auto count = countsView[i];
        auto offset = offsetsView[i];
        for(int j = 0; j < count; ++j) candToTetIdIdView[offset + j] = i;
      });

    // Find which shape bounding boxes intersect hexahedron bounding boxes
    SLIC_INFO(axom::fmt::format("Finding shape candidates for {} tet elements ", tetIndices.size()));

    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;
    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = false;

    ClippingStats clipStats;

    if(useTets)
    {
      axom::for_all<ExecSpace>(
        candidates.size(),
        AXOM_LAMBDA(axom::IndexType iCand) {
          auto pieceId = candidatesView[iCand];
          const axom::primal::Tetrahedron<double, 3>& geomPiece = geomTetsView[pieceId];

          auto tetIdId = candToTetIdIdView[iCand];
          auto tetId = tetIndices[tetIdId];
          const auto& tet = meshTets[tetId];
          if(tet.degenerate())
          {
            return;
          }

          const auto poly = primal::clip<double>(tet, geomPiece, EPS, tryFixOrientation);

          if(poly.numVertices() >= 4)
          {
            // Poly is valid
            double volume = poly.volume();
            SLIC_ASSERT(volume >= 0);
            auto cellId = tetId / NUM_TETS_PER_HEX;
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
          }
        });
    }
    else  // useOcts
    {
      axom::for_all<ExecSpace>(
        candidates.size(),
        AXOM_LAMBDA(axom::IndexType iCand) {
          auto pieceId = candidatesView[iCand];
          const axom::primal::Octahedron<double, 3>& geomPiece = geomOctsView[pieceId];

          auto tetIdId = candToTetIdIdView[iCand];
          auto tetId = tetIndices[tetIdId];
          const auto& tet = meshTets[tetId];
          if(tet.degenerate())
          {
            return;
          }

          const auto poly = primal::clip<double>(tet, geomPiece, EPS, tryFixOrientation);

          if(poly.numVertices() >= 4)
          {
            // Poly is valid
            double volume = poly.volume();
            SLIC_ASSERT(volume >= 0);
            auto cellId = tetId / NUM_TETS_PER_HEX;
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
          }
        });
    }

    clipStats.copyTo(statistics);
    statistics["clipsCandidates"].set_int64(candidates.size());

  }  // end of computeClipVolumes3DTets() function

  void getLabelCounts(axom::ArrayView<const LabelType> labels,
                      std::int64_t& inCount,
                      std::int64_t& onCount,
                      std::int64_t& outCount) override
  {
    AXOM_ANNOTATE_SCOPE("MeshClipper::getLabelCounts");
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
};

}  // end namespace detail
}  // namespace experimental
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_MESHCLIPPERIMPL_HPP_
