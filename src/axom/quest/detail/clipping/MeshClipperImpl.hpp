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
  using Plane3DType = primal::Plane<double, 3>;
  using BoundingBoxType = primal::BoundingBox<double, 3>;
  using TetrahedronType = primal::Tetrahedron<double, 3>;
  using OctahedronType = primal::Octahedron<double, 3>;
  using CoordTransformer = primal::experimental::CoordinateTransformer<double>;

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

    AXOM_ANNOTATE_SCOPE("MeshClipper:collect_indices");
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

    axom::Array<axom::IndexType> tmpLabels(ArrayOptions::Uninitialized(),
                                           1 + labels.size(),
                                           0,
                                           labels.getAllocatorID());
    tmpLabels.fill(0, 1, 0);
    auto tmpLabelsView = tmpLabels.view();
    axom::ReduceSum<ExecSpace, IndexType> onCountReduce {0};
    axom::for_all<ExecSpace>(
      labelCount,
      AXOM_LAMBDA(axom::IndexType ci) {
        bool isOn = labels[ci] == LabelType::LABEL_ON;
        tmpLabelsView[1 + ci] = isOn;
        onCountReduce += isOn;
      });

    RAJA::inclusive_scan_inplace<ScanPolicy>(RAJA::make_span(tmpLabels.data(), tmpLabels.size()),
                                             RAJA::operators::plus<axom::IndexType> {});

    // Space for output index list
    axom::IndexType onCount = onCountReduce.get();
    if(onIndices.size() < onCount || onIndices.getAllocatorID() != labels.getAllocatorID())
    {
      onIndices = axom::Array<axom::IndexType> {axom::ArrayOptions::Uninitialized(),
                                                onCount,
                                                0,
                                                labels.getAllocatorID()};
    }

    auto onIndicesView = onIndices.view();
    axom::for_all<ExecSpace>(
      1,
      1 + labelCount,
      AXOM_LAMBDA(axom::IndexType i) {
        if(tmpLabelsView[i] != tmpLabelsView[i - 1])
        {
          onIndicesView[tmpLabelsView[i] - 1] = i - 1;
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
    // Get the BVH traverser for doing the 2-pass search manually.
    const auto bvhTraverser = bvh.getTraverser();
    /*
     * Predicate for traversing the BVH.  We enter BVH nodes
     * whose bounding boxes intersect the query bounding box.
     */
    auto traversePredTetId = [=] AXOM_HOST_DEVICE(const IndexType& queryTetId,
                                                  const BoundingBoxType& bvhBbox) -> bool {
      const auto& queryTet = meshTets[tetIndices[queryTetId]];
      return tetBoxCollision(queryTet, bvhBbox);
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
          auto& tetId = tetIndices[iTet];
          const auto& meshTet = meshTets[tetId];

          auto pieceId = leafNodes[currentNode];
          if(useTets)
          {
            const auto& piece = geomTetsView[pieceId];
            if(tetTetCollision(meshTet, piece))
            {
              ++count;
            }
          }
          else
          {
            const auto& piece = geomOctsView[pieceId];
            if(tetOctCollision(meshTet, piece))
            {
              ++count;
            }
          }
        };
        bvhTraverser.traverse_tree(iTet, countCollisions, traversePredTetId);
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

        /*
         * Record indices of the tet and the candidate that collided.
         * Unless tet and candidate can be shown not to collide.
         */
        auto recordCollision = [&](std::int32_t currentNode, const std::int32_t* leafs) {
          auto& tetId = tetIndices[iTet];
          const auto& meshTet = meshTets[tetId];
          auto pieceId = leafs[currentNode];
          bool record = false;
          if(useTets)
          {
            const auto& piece = geomTetsView[pieceId];
            if(tetTetCollision(meshTet, piece))
            {
              record = true;
            }
          }
          else
          {
            const auto& piece = geomOctsView[pieceId];
            if(tetOctCollision(meshTet, piece))
            {
              record = true;
            }
          }
          if(record)
          {
            candToTetIdIdView[offset] = iTet;
            candidatesView[offset] = pieceId;
            ++offset;
          }
        };

        bvhTraverser.traverse_tree(iTet, recordCollision, traversePredTetId);
      });
    AXOM_ANNOTATE_END("MeshClipper:find_candidates");

    SLIC_DEBUG(axom::fmt::format(
      "Running clip loop on {} candidate pieces for the select {} tets of the full {} mesh cells",
      candidates.size(),
      tetCount,
      shapeMesh.getCellCount()));

    ClippingStats clipStats;

    const auto screenLevel = myClipper().getScreenLevel();

    /*
     * Now we have the candidates.  Do the clip loop.
     */
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
          computeMeshTetGeomPieceOverlap(meshTet,
                                         geomPiece,
                                         ovlap.data() + cellId,
                                         clipStats,
                                         screenLevel);
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
          computeMeshTetGeomPieceOverlap(meshTet,
                                         geomPiece,
                                         ovlap.data() + cellId,
                                         clipStats,
                                         screenLevel);
        });
    }
    AXOM_ANNOTATE_END("MeshClipper:clipLoop_tetScreened");

    clipStats.copyTo(statistics);
    statistics["clipsCandidates"].set_int64(candidates.size());
  }  // end of computeClipVolumes3DTets() function

  /*!
   * @brief Get the geometry in discrete pieces,
   *   which can be tets or octs, and place them in a search tree.
   * @return true if geometry is composed of tetrahedra, false if octahedra.
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
  AXOM_HOST_DEVICE static inline double geomPieceVolume(const TetrahedronType& tet)
  {
    return tet.volume();
  }

  /*!
   * @brief Volume of a octahedron from discretized geometry.
   *
   * Assumes octahedron is convex.
   */
  AXOM_HOST_DEVICE static inline double geomPieceVolume(const OctahedronType& oct)
  {
    // Oct vertex indices of the four tets that the oct decomposes into.
    IndexType tIds[4][4] = {{0, 3, 1, 2}, {0, 3, 2, 4}, {0, 3, 4, 5}, {0, 3, 5, 1}};
    double octVol = 0.0;
    for(int i = 0; i < 4; ++i)
    {
      TetrahedronType tet(oct[tIds[i][0]], oct[tIds[i][1]], oct[tIds[i][2]], oct[tIds[i][3]]);
      double tetVol = tet.signedVolume();
      octVol -= tetVol;  // Octs from the discretized geometries are inverted w.r.t. tIDs.
    }
    SLIC_ASSERT(octVol > 0.0);
    return octVol;
  }

  /*!
   * @brief Compute overlap volume between a reference tet (from the shape mesh)
   * and a piece (tet or oct) of the discretized geometry.
   *
   * Because primal::clip is so expensive, we do a conservative
   * overlap check on @c meshTet and @c geomPiece to avoid clipping.
   *
   * @return results of check whether the piece is IN/ON/OUT of the tet.
   *
   * @tparam TetOrOctType either a TetrahedronType or OctahedronType,
   * the two types a geometry can be discretized into.
   */
  template <typename TetOrOctType>
  AXOM_HOST_DEVICE static inline LabelType computeMeshTetGeomPieceOverlap(
    const TetrahedronType& meshTet,
    const TetOrOctType& geomPiece,
    double* overlapVolume,
    const ClippingStats& clipStats,
    int screenLevel)
  {
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;
    constexpr bool tryFixOrientation = false;
    if(screenLevel >= 3)
    {
      LabelType geomLabel = labelPieceInOutOfTet(meshTet, geomPiece);
      if(geomLabel == LabelType::LABEL_OUT)
      {
        clipStats.outSum += 1;
        return geomLabel;
      }
      if(geomLabel == LabelType::LABEL_IN)
      {
        auto contribVol = geomPieceVolume(geomPiece);
        RAJA::atomicAdd<ATOMIC_POL>(overlapVolume, contribVol);
        clipStats.inSum += 1;
        return geomLabel;
      }
    }

    clipStats.onSum += 1;
    const auto poly = primal::clip<double>(meshTet, geomPiece, EPS, tryFixOrientation);
    if(poly.numVertices() >= 4)
    {
      // Poly is valid
      auto contribVol = poly.volume();
      SLIC_ASSERT(contribVol >= 0);
      RAJA::atomicAdd<ATOMIC_POL>(overlapVolume, contribVol);
    }
    else
    {
      clipStats.missSum += 1;
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
  AXOM_HOST_DEVICE static inline LabelType labelPieceInOutOfTet(const TetrahedronType& tet,
                                                                const TetOrOctType& piece)
  {
    Point3DType unitTet[] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    CoordTransformer toUnitTet(&tet[0], unitTet);

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
      // h4 of 1 is right at the unitTet's height of sqrt(3).
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

  /*!
   * @brief Whether a tet and a bounding box (possibly) intersect.
   * Answer may be a false positive but never a false negative
   * (which is why this code lives here instead of in a primal::intersect method).
   */
  AXOM_HOST_DEVICE static inline bool tetBoxCollision(const TetrahedronType& tet,
                                                      const BoundingBoxType& box)
  {
    if(box.contains(tet[0]) || box.contains(tet[1]) || box.contains(tet[2]) || box.contains(tet[3]))
    {
      return true;
    }

    Point3DType unitTet[] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    CoordTransformer toUnitTet(&tet[0], unitTet);

    int vsAbove[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int vsBelow[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for(int i = 0; i < 8; ++i)
    {
      Point3DType boxVert {(i & 1) == 0 ? box.getMin()[0] : box.getMax()[0],
                           (i & 2) == 0 ? box.getMin()[1] : box.getMax()[1],
                           (i & 4) == 0 ? box.getMin()[2] : box.getMax()[2]};
      toUnitTet.transform(boxVert.array());
      // h4 is height of boxVert above the diagonal face, scaled by sqrt(3).
      // h4 of 1 is right at the unitTet's height of sqrt(3).
      double h4 = 1 - (boxVert[0] + boxVert[1] + boxVert[2]);
      vsAbove[0] += boxVert[0] >= 1;
      vsAbove[1] += boxVert[1] >= 1;
      vsAbove[2] += boxVert[2] >= 1;
      vsAbove[3] += h4 >= 1;
      vsBelow[0] += boxVert[0] <= 0;
      vsBelow[1] += boxVert[1] <= 0;
      vsBelow[2] += boxVert[2] <= 0;
      vsBelow[3] += h4 <= 0;
    }
    if(vsAbove[0] == 8 || vsAbove[1] == 8 || vsAbove[2] == 8 || vsAbove[3] == 8 ||
       vsBelow[0] == 8 || vsBelow[1] == 8 || vsBelow[2] == 8 || vsBelow[3] == 8)
    {
      return false;
    }
    return true;
  }

  /*!
   * @brief Whether a tet and the convex hull of an octahedron (possibly) intersect.
   * Answer may be a false positive but never a false negative
   * (which is why this code lives here instead of in a primal::intersect method).
   */
  AXOM_HOST_DEVICE static inline bool tetOctCollision(const TetrahedronType& tet,
                                                      const OctahedronType& oct)
  {
    Point3DType unitTet[] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    CoordTransformer toUnitTet(&tet[0], unitTet);
    int octVertsAbove[OctahedronType::NUM_VERTS] = {0, 0, 0, 0, 0, 0};
    int octVertsBelow[OctahedronType::NUM_VERTS] = {0, 0, 0, 0, 0, 0};
    for(int i = 0; i < OctahedronType::NUM_VERTS; ++i)
    {
      auto octVert = toUnitTet.getTransformed(oct[i]);
      // h4 is height of octVert above the diagonal face, scaled by sqrt(3).
      // h4 of 1 is right at the unitTet's height of sqrt(3).
      double h4 = 1 - (octVert[0] + octVert[1] + octVert[2]);
      octVertsAbove[0] += octVert[0] >= 1;
      octVertsAbove[1] += octVert[1] >= 1;
      octVertsAbove[2] += octVert[2] >= 1;
      octVertsAbove[3] += h4 >= 1;
      octVertsBelow[0] += octVert[0] <= 0;
      octVertsBelow[1] += octVert[1] <= 0;
      octVertsBelow[2] += octVert[2] <= 0;
      octVertsBelow[3] += h4 <= 0;
    }
    if(octVertsAbove[0] == OctahedronType::NUM_VERTS ||
       octVertsAbove[1] == OctahedronType::NUM_VERTS ||
       octVertsAbove[2] == OctahedronType::NUM_VERTS ||
       octVertsAbove[3] == OctahedronType::NUM_VERTS ||
       octVertsBelow[0] == OctahedronType::NUM_VERTS ||
       octVertsBelow[1] == OctahedronType::NUM_VERTS ||
       octVertsBelow[2] == OctahedronType::NUM_VERTS || octVertsBelow[3] == OctahedronType::NUM_VERTS)
    {
      return false;
    }

    // Indices of the vertices of each of the 8 faces of the octagon, oriented inside.
    using ThreeIds = int[3];
    ThreeIds octFIds[8] =
      {{0, 2, 4}, {0, 5, 1}, {2, 1, 3}, {4, 3, 5}, {1, 5, 3}, {3, 4, 2}, {5, 0, 4}, {1, 2, 0}};
    for(int fi = 0; fi < 8; ++fi)
    {
      // Construct  plane of face fi of the oct.
      const ThreeIds& fIds = octFIds[fi];
      auto& v0 = oct[fIds[0]];
      auto& v1 = oct[fIds[1]];
      auto& v2 = oct[fIds[2]];
      axom::primal::Vector<double, 3> r1(v0, v1);
      axom::primal::Vector<double, 3> r2(v0, v2);
      axom::primal::Vector<double, 3> normal = axom::primal::Vector<double, 3>::cross_product(r1, r2);
      if(normal.squared_norm() < EPS)
      {
        continue;
      }  // Skip degenerate face
      axom::primal::Plane<double, 3> octBase(normal, v0);

      // Compute height range of vertices not in face fi.
      double maxHeight = 0.0;
      double minHeight = 0.0;
      const ThreeIds& nonFids = octFIds[(fi + 4) % 8];  // 3 oct vertices not part of face fi.
      for(int vi = 0; vi < 3; ++vi)
      {
        const auto& vert = oct[nonFids[vi]];
        double vertHeight = octBase.signedDistance(vert);
        if(maxHeight < vertHeight)
        {
          maxHeight = vertHeight;
        }
        if(minHeight > vertHeight)
        {
          minHeight = vertHeight;
        }
      }

      // Number of tet vertices outside [minHeight,maxHeight]..
      int tetVertsAbove = 0;
      int tetVertsBelow = 0;
      for(int ti = 0; ti < 4; ++ti)
      {
        const auto& tetVert = tet[ti];
        double tetVertHeight = octBase.signedDistance(tetVert);
        tetVertsAbove += tetVertHeight >= maxHeight;
        tetVertsBelow += tetVertHeight <= minHeight;
      }

      if(tetVertsAbove == TetrahedronType::NUM_VERTS || tetVertsBelow == TetrahedronType::NUM_VERTS)
      {
        return false;
      }
    }
    return true;
  }

  /*!
   * @brief Whether a tet and another tet (possibly) intersect.
   * Answer may be a false positive but never a false negative
   * (which is why this code lives here instead of in a primal::intersect method).
   */
  AXOM_HOST_DEVICE static inline bool tetTetCollision(const TetrahedronType& tetA,
                                                      const TetrahedronType& tetB,
                                                      bool flip = true)
  {
    Point3DType unitTet[] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    CoordTransformer toUnitTet(&tetA[0], unitTet);

    int vsAbove[TetrahedronType::NUM_VERTS] = {0, 0, 0, 0};
    int vsBelow[TetrahedronType::NUM_VERTS] = {0, 0, 0, 0};
    for(int i = 0; i < TetrahedronType::NUM_VERTS; ++i)
    {
      const auto bVert = toUnitTet.getTransformed(tetB[i].array());
      // h4 is height of bVert above the diagonal face, scaled by sqrt(3).
      // h4 of 1 is right at the unitTet's height of sqrt(3).
      double h4 = 1 - (bVert[0] + bVert[1] + bVert[2]);
      vsAbove[0] += bVert[0] >= 1;
      vsAbove[1] += bVert[1] >= 1;
      vsAbove[2] += bVert[2] >= 1;
      vsAbove[3] += h4 >= 1;
      vsBelow[0] += bVert[0] <= 0;
      vsBelow[1] += bVert[1] <= 0;
      vsBelow[2] += bVert[2] <= 0;
      vsBelow[3] += h4 <= 0;
    }
    if(vsAbove[0] == TetrahedronType::NUM_VERTS || vsAbove[1] == TetrahedronType::NUM_VERTS ||
       vsAbove[2] == TetrahedronType::NUM_VERTS || vsAbove[3] == TetrahedronType::NUM_VERTS)
    {
      return false;
    }

    if(flip)
    {
      // Cannot claim no-intersection checking whether tetB above or below tetA.
      // So try checking whether tetA is above or below tetB.
      return tetTetCollision(tetB, tetA, false);
    }

    return true;
  }

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
  static constexpr double EPS = 1e-10;
  static constexpr double BVH_SCALE_FACTOR = 1.0;
};

}  // end namespace detail
}  // namespace experimental
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_MESHCLIPPERIMPL_HPP_
