// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_GEOMETRYCLIPPERDELEGATEEXEC_HPP_
#define AXOM_GEOMETRYCLIPPERDELEGATEEXEC_HPP_

#ifndef AXOM_USE_RAJA
  #error "quest::GeometryClipper requires RAJA."
#endif

#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/quest/GeometryClipper.hpp"
#include "axom/spin/BVH.hpp"
#include "RAJA/RAJA.hpp"

namespace axom
{
namespace quest
{
namespace detail
{

template <typename ExecSpace>
class GeometryClipperDelegateExec : public GeometryClipper::Delegate
{
public:
  using LabelType = GeometryClipper::LabelType;

  GeometryClipperDelegateExec(GeometryClipper& delegator) : GeometryClipper::Delegate(delegator) { }

  void initVolumeOverlaps(const axom::ArrayView<GeometryClipperStrategy::LabelType>& labels,
                          axom::ArrayView<double> ovlap) override
  {
    const axom::IndexType cellCount = getShapeeMesh().getCellCount();
    SLIC_ASSERT(labels.size() == cellCount);
    SLIC_ASSERT(ovlap.size() == cellCount);

    auto cellVolumes = getShapeeMesh().getCellVolumes();

    /*
      Overlap volumes is cell volume for cells inside geometry.
      Cells outside the geometry should be zero as are cells.
      Cells on boundary to be zeroed for accumulating by clipping remainder.
    */
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        auto& l = labels[i];
        ovlap[i] = l == GeometryClipperStrategy::LABEL_IN ? cellVolumes[i] : 0.0;
      });

    return;
  }

  void initVolumeOverlaps(axom::ArrayView<double> ovlap) override
  {
    SLIC_ASSERT(ovlap.size() == getShapeeMesh().getCellCount());
    ovlap.fill(0.0);
    return;
  }

  void addVolumesOfInteriorTets(axom::ArrayView<const axom::IndexType> cellsOnBdry,
                                axom::ArrayView<const LabelType> tetLabels,
                                axom::ArrayView<double> ovlap) override
  {
    auto meshTets = getShapeeMesh().getCellsAsTets();

    const axom::IndexType hexCount = cellsOnBdry.size();

    SLIC_ASSERT(tetLabels.size() == TETS_PER_HEXAHEDRON * hexCount);

    axom::for_all<ExecSpace>(
      hexCount,
      AXOM_LAMBDA(axom::IndexType ih) {
        const axom::IndexType hexId = cellsOnBdry[ih];
        const LabelType* tetLabelsAtHex = &tetLabels[TETS_PER_HEXAHEDRON * ih];
        for(int it = 0; it < TETS_PER_HEXAHEDRON; ++it)
        {
          if(tetLabelsAtHex[it] == GeometryClipperStrategy::LABEL_IN)
          {
            const axom::IndexType tetId = hexId * TETS_PER_HEXAHEDRON + it;
            const auto& tet = meshTets[tetId];
            ovlap[hexId] += tet.volume();
          }
        }
      });
  }

  //! @brief Make an list of indices where labels have value LABEL_ON.
  void collectOnIndices(const axom::ArrayView<LabelType>& labels,
                        axom::Array<axom::IndexType>& onIndices) override
  {
    if(labels.empty()) {return;};

    AXOM_ANNOTATE_SCOPE("GeometryClipper::collect_unlabeleds");
    /*!
      1. Generate tmpLabels, having a value of 1 where labels is LABEL_ON and zero elsewhere.
      2. Inclusive scan on tmpLabels to generate values that step up at unlabeled cells.
      3. Find unlabeled cells by seeing where tmpLabels changes values.
         (Handle first cell separately, then loop from second cell on.)
    */
    using ScanPolicy = typename axom::execution_space<ExecSpace>::loop_policy;

    const axom::IndexType labelCount = labels.size();

    axom::Array<axom::IndexType> tmpLabels(labels.shape(), labels.getAllocatorID());
    auto tmpLabelsView = tmpLabels.view();
    axom::for_all<ExecSpace>(
      labelCount,
      AXOM_LAMBDA(axom::IndexType ci) { tmpLabelsView[ci] = labels[ci] == GeometryClipperStrategy::LABEL_ON; });

    RAJA::inclusive_scan_inplace<ScanPolicy>(RAJA::make_span(tmpLabels.data(), tmpLabels.size()),
                                             RAJA::operators::plus<axom::IndexType> {});

    axom::IndexType unlabeledCount;
    axom::copy(&unlabeledCount, &tmpLabels.back(), sizeof(unlabeledCount));

    // Re-allocate output if it's too small or has wrong allocator.
    if(onIndices.size() < unlabeledCount ||
       onIndices.getAllocatorID() != labels.getAllocatorID())
    {
      onIndices = axom::Array<axom::IndexType> {axom::ArrayOptions::Uninitialized(),
                                                     unlabeledCount,
                                                     unlabeledCount,
                                                     labels.getAllocatorID()};
    }
    auto onIndicesView = onIndices.view();

    LabelType firstLabel = '\0';
    axom::copy(&firstLabel, &labels[0], sizeof(firstLabel));
    if(firstLabel == 1)
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

  void remapTetIndices(axom::ArrayView<axom::IndexType> tetsOnBdry,
                       axom::ArrayView<const axom::IndexType> cellsOnBdry) override
  {
    /*
      cellsOnBdry is a list of cell indices.

      Each cell has N = 24 tets.

      tetsOnBdry are a list of indices referring to the tets in those
      cells.  N tets for each cell.  So the values in tetsOnBDry are
      in [0, cellsOnBdry.size()*24).

      Use the indices in cellsOnBdry to map tetsOnBdry values to the
      indices of the full set of tets.
    */
    if(tetsOnBdry.empty()) { return; }

    axom::for_all<ExecSpace>(
      tetsOnBdry.size(),
      AXOM_LAMBDA(axom::IndexType i) {
        auto tetIdId = tetsOnBdry[i];
        auto cellIdId = tetIdId/24;
        auto cellId = cellsOnBdry[cellIdId];
        auto tetIdInCell = tetIdId % 24;
        auto tetId = cellId * 24 + tetIdInCell;
        tetsOnBdry[i] = tetId;
      });
  }

  void computeClipVolumes3D(axom::ArrayView<double> ovlap) override
  {
    AXOM_ANNOTATE_SCOPE("GeometryClipper::computeClipVolumes3D");

    using BoundingBoxType = primal::BoundingBox<double, 3>;

    ShapeeMesh& shapeeMesh = getDelegator().getShapeeMesh();

    const int allocId = shapeeMesh.getAllocatorID();

    const IndexType cellCount = shapeeMesh.getCellCount();

    SLIC_INFO(axom::fmt::format(
      "GeometryClipper::computeClipVolumes3D: Getting discrete geometry for shape '{}'",
      getStrategy().name()));

    //
    // Get the geometry in discrete pieces, which can be tets or octs.
    //
    auto& strategy = getStrategy();
    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    const bool useOcts = strategy.getGeometryAsOcts(shapeeMesh, geomAsOcts);
    const bool useTets = strategy.getGeometryAsTets(shapeeMesh, geomAsTets);
    SLIC_ASSERT(useOcts || geomAsOcts.empty());
    SLIC_ASSERT(useTets || geomAsTets.empty());
    if(useTets == useOcts)
    {
      SLIC_ERROR(axom::fmt::format("Problem with GeometryClipperStrategy implementation '{}'."
                                   "  Implementations that don't provide a specializedClip function"
                                   " must provide exactly one getGeometryAsOcts() or getGeometryAsTets()."
                                   "  This implementation provides {}.", strategy.name(),
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

    axom::ArrayView<const BoundingBoxType> cellBbsView = shapeeMesh.getCellBoundingBoxes();

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
    axom::Array<IndexType> hexIndices(candidateCount * TETS_PER_HEXAHEDRON,
                                              candidateCount * TETS_PER_HEXAHEDRON,
                                              allocId);
    auto hexIndicesView = hexIndices.view();

    axom::Array<IndexType> shapeCandidates(candidateCount * TETS_PER_HEXAHEDRON,
                                                   candidateCount * TETS_PER_HEXAHEDRON,
                                                   allocId);
    auto shapeCandidatesView = shapeCandidates.view();

    // Tetrahedrons from hexes (24 for each hex)
    auto cellsAsTets = shapeeMesh.getCellsAsTets();

    // Index into 'tets'
    axom::Array<IndexType> tetIndices(candidateCount * TETS_PER_HEXAHEDRON,
                                              candidateCount * TETS_PER_HEXAHEDRON,
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

            for(int k = 0; k < TETS_PER_HEXAHEDRON; k++)
            {
              IndexType idx = RAJA::atomicAdd<ATOMIC_POL>(tetCandidatesCountPtr, IndexType {1});
              hexIndicesView[idx] = i;
              shapeCandidatesView[idx] = shapeIdx;
              tetIndicesView[idx] = i * TETS_PER_HEXAHEDRON + k;
            }
          }
        });
    }

    SLIC_INFO(axom::fmt::format(
      "Running clip loop on {} candidate tets for of all {} hexes in the mesh",
      tetCandidatesCount,
      cellCount));

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = false;

    {
      tetCandidatesCount = TETS_PER_HEXAHEDRON*candidates.size();
      AXOM_ANNOTATE_SCOPE("GeometryClipper::clipLoop");
#if defined(AXOM_DEBUG)
      // Verifying: this should always pass.
      if(tetCandidatesCountPtr != &tetCandidatesCount)
      {
        axom::copy(&tetCandidatesCount, tetCandidatesCountPtr, sizeof(IndexType));
      }
      SLIC_ASSERT(tetCandidatesCount == candidateCount*TETS_PER_HEXAHEDRON);
#endif

      if(useTets)
      {
        axom::for_all<ExecSpace>(
          tetCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            const int index = hexIndicesView[i];
            const int shapeIndex = shapeCandidatesView[i];
            const int tetIndex = tetIndicesView[i];

            const auto poly = primal::clip<double, MAX_VERTS_FOR_TET_CLIPPING, MAX_NBRS_PER_VERT_FOR_TET_CLIPPING>(geomTetsView[shapeIndex],
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
      else // useOcts
      {
        axom::for_all<ExecSpace>(
          tetCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            const int index = hexIndicesView[i];
            const int shapeIndex = shapeCandidatesView[i];
            const int tetIndex = tetIndicesView[i];

            const auto poly = primal::clip<double, MAX_VERTS_FOR_OCT_CLIPPING, MAX_NBRS_PER_VERT_FOR_OCT_CLIPPING>(geomOctsView[shapeIndex],
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

    if(tetCandidatesCountPtr != &tetCandidatesCount)
    {
      axom::deallocate(tetCandidatesCountPtr);
    }
  }  // end of computeClipVolumes3D() function

  void computeClipVolumes3D(const axom::ArrayView<axom::IndexType>& cellIndices,
                            axom::ArrayView<double> ovlap) override

  {
    AXOM_ANNOTATE_SCOPE("GeometryClipper::computeClipVolumes3D:limited");

    using BoundingBoxType = primal::BoundingBox<double, 3>;

    ShapeeMesh& shapeeMesh = getDelegator().getShapeeMesh();

    const int allocId = shapeeMesh.getAllocatorID();

    const IndexType cellCount = shapeeMesh.getCellCount();

    SLIC_INFO(axom::fmt::format(
      "GeometryClipper::computeClipVolumes3D: Getting discrete geometry for shape '{}'",
      getStrategy().name()));

    auto& strategy = getStrategy();
    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    const bool useOcts = strategy.getGeometryAsOcts(shapeeMesh, geomAsOcts);
    const bool useTets = strategy.getGeometryAsTets(shapeeMesh, geomAsTets);
    SLIC_ASSERT(useOcts || geomAsOcts.empty());
    SLIC_ASSERT(useTets || geomAsTets.empty());
    if(useTets == useOcts)
    {
      SLIC_ERROR(axom::fmt::format("Problem with GeometryClipperStrategy implementation '{}'."
                                   "  Implementations that don't provide a specializedClip function"
                                   " must provide exactly one getGeometryAsOcts() or getGeometryAsTets()."
                                   "  This implementation provides {}.", strategy.name(),
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
    axom::ArrayView<const BoundingBoxType> cellBbsView = shapeeMesh.getCellBoundingBoxes();
    axom::Array<BoundingBoxType> limitedCellBbs(limitedCellCount, limitedCellCount, allocId);
    axom::ArrayView<BoundingBoxType> limitedCellBbsView = limitedCellBbs.view();
    axom::for_all<ExecSpace>(limitedCellBbsView.size(),
                             AXOM_LAMBDA(axom::IndexType i) {
                               limitedCellBbsView[i] = cellBbsView[cellIndices[i]];
                               });

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
    axom::Array<IndexType> hexIndices(candidateCount * TETS_PER_HEXAHEDRON,
                                              candidateCount * TETS_PER_HEXAHEDRON,
                                              allocId);
    auto hexIndicesView = hexIndices.view();

    axom::Array<IndexType> shapeCandidates(candidateCount * TETS_PER_HEXAHEDRON,
                                                   candidateCount * TETS_PER_HEXAHEDRON,
                                                   allocId);
    auto shapeCandidatesView = shapeCandidates.view();

    // Tetrahedrons from hexes (24 for each hex)
    auto cellsAsTets = shapeeMesh.getCellsAsTets();

    // Index into 'tets'
    axom::Array<IndexType> tetIndices(candidateCount * TETS_PER_HEXAHEDRON,
                                              candidateCount * TETS_PER_HEXAHEDRON,
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

            for(int k = 0; k < TETS_PER_HEXAHEDRON; k++)
            {
              IndexType idx = RAJA::atomicAdd<ATOMIC_POL>(tetCandidatesCountPtr, IndexType {1});
              hexIndicesView[idx] = i;
              shapeCandidatesView[idx] = shapeIdx;
              tetIndicesView[idx] = i * TETS_PER_HEXAHEDRON + k;
            }
          }
        });
    }

    SLIC_INFO(axom::fmt::format(
      "Running clip loop on {} candidate tets for the select {} hexes of the full {} cells",
      tetCandidatesCount,
      cellIndices.size(), cellCount));

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = false;

    {
      tetCandidatesCount = TETS_PER_HEXAHEDRON*candidates.size();
      AXOM_ANNOTATE_SCOPE("GeometryClipper::clipLoop_limited");
#if defined(AXOM_DEBUG)
      // Verifying: this should always pass.
      if(tetCandidatesCountPtr != &tetCandidatesCount)
      {
        axom::copy(&tetCandidatesCount, tetCandidatesCountPtr, sizeof(IndexType));
      }
      SLIC_ASSERT(tetCandidatesCount == candidateCount*TETS_PER_HEXAHEDRON);
#endif

      if(useTets)
      {
        axom::for_all<ExecSpace>(
          tetCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            int index = hexIndicesView[i]; // index into limited mesh hex array
            index = cellIndices[index]; // Now, it indexes into the full hex array.

            const int shapeIndex = shapeCandidatesView[i]; // index into pieces array
            int tetIndex = tetIndicesView[i]; // index into BVH results, implicit because BVH results specify hexes, not tets.
            int tetIndex1 = tetIndex/TETS_PER_HEXAHEDRON;
            int tetIndex2 = tetIndex%TETS_PER_HEXAHEDRON;
            tetIndex = cellIndices[tetIndex1]*TETS_PER_HEXAHEDRON + tetIndex2; // Now it indexes into the full tets-from-hexes array.

            const auto poly = primal::clip<double, MAX_VERTS_FOR_TET_CLIPPING, MAX_NBRS_PER_VERT_FOR_TET_CLIPPING>(geomTetsView[shapeIndex],
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
      else // useOcts
      {
        axom::for_all<ExecSpace>(
          tetCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            int index = hexIndicesView[i]; // index into limited mesh hex array
            index = cellIndices[index]; // Now, it indexes into the full hex array.

            const int shapeIndex = shapeCandidatesView[i]; // index into pieces array
            int tetIndex = tetIndicesView[i]; // index into BVH results, implicit because BVH results specify hexes, not tets.
            int tetIndex1 = tetIndex/TETS_PER_HEXAHEDRON;
            int tetIndex2 = tetIndex%TETS_PER_HEXAHEDRON;
            tetIndex = cellIndices[tetIndex1]*TETS_PER_HEXAHEDRON + tetIndex2; // Now it indexes into the full tets-from-hexes array.

            const auto poly = primal::clip<double, MAX_VERTS_FOR_OCT_CLIPPING, MAX_NBRS_PER_VERT_FOR_OCT_CLIPPING>(geomOctsView[shapeIndex],
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

    if(tetCandidatesCountPtr != &tetCandidatesCount)
    {
      axom::deallocate(tetCandidatesCountPtr);
    }
  }  // end of computeClipVolumes3D() function

  void computeClipVolumes3DTets(const axom::ArrayView<axom::IndexType>& tetIndices,
                                axom::ArrayView<double> ovlap) override

  {
    AXOM_ANNOTATE_SCOPE("GeometryClipper::computeClipVolumes3D:limited");

    using Point3DType = primal::Point<double, 3>;
    using BoundingBoxType = primal::BoundingBox<double, 3>;
    using TetrahedronType = primal::Tetrahedron<double, 3>;
    using OctahedronType = primal::Octahedron<double, 3>;

    ShapeeMesh& shapeeMesh = getDelegator().getShapeeMesh();
    auto meshTets = getShapeeMesh().getCellsAsTets();

    const int allocId = shapeeMesh.getAllocatorID();

    SLIC_INFO(axom::fmt::format(
      "GeometryClipper::computeClipVolumes3D: Getting discrete geometry for shape '{}'",
      getStrategy().name()));

    auto& strategy = getStrategy();
    axom::Array<TetrahedronType> geomAsTets;
    axom::Array<OctahedronType> geomAsOcts;
    const bool useOcts = strategy.getGeometryAsOcts(shapeeMesh, geomAsOcts);
    const bool useTets = strategy.getGeometryAsTets(shapeeMesh, geomAsTets);
    SLIC_ASSERT(useOcts || geomAsOcts.empty());
    SLIC_ASSERT(useTets || geomAsTets.empty());
    if(useTets == useOcts)
    {
      SLIC_ERROR(axom::fmt::format("Problem with GeometryClipperStrategy implementation '{}'."
                                   "  Implementations that don't provide a specializedClip function"
                                   " must provide exactly one getGeometryAsOcts() or getGeometryAsTets()."
                                   "  This implementation provides {}.", strategy.name(),
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
    axom::for_all<ExecSpace>(tetCount,
                             AXOM_LAMBDA(axom::IndexType i) {
                               auto& tetBb = tetBbsView[i];
                               axom::IndexType tetId = tetIndices[i];
                               const auto& tet = meshTets[tetId];
                               for(int j = 0; j < 4; ++j)
                                 tetBb.addPoint(tet[j]);
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
    axom::for_all<ExecSpace>(tetCount,
      AXOM_LAMBDA(axom::IndexType i) {
        auto count = countsView[i];
        auto offset = offsetsView[i];
        for(int j = 0; j < count; ++j) candToTetIdIdView[offset + j] = i;
    });


    // Find which shape bounding boxes intersect hexahedron bounding boxes
    SLIC_INFO(
      axom::fmt::format("Finding shape candidates for {} tet elements ", tetIndices.size()));

    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;
    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = false;

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

          const auto poly = primal::clip<double,
                                         MAX_VERTS_FOR_TET_CLIPPING,
                                         MAX_NBRS_PER_VERT_FOR_TET_CLIPPING>(
                                           tet,
                                           geomPiece,
                                           EPS,
                                           tryFixOrientation);

          if(poly.numVertices() >= 4)
          {
            // Poly is valid
            double volume = poly.volume();
            SLIC_ASSERT(volume >= 0);
            auto cellId = tetId / TETS_PER_HEXAHEDRON;
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
          }
        });
    }
    else // useOcts
    {
      axom::for_all<ExecSpace>(
        candidates.size(),
        AXOM_LAMBDA(axom::IndexType iCand) {
          auto pieceId = candidatesView[iCand];
          const axom::primal::Octahedron<double, 3>& geomPiece = geomOctsView[pieceId];

          auto tetIdId = candToTetIdIdView[iCand];
          auto tetId = tetIndices[tetIdId];
          const auto& tet = meshTets[tetId];

          const auto poly = primal::clip<double,
                                         MAX_VERTS_FOR_TET_CLIPPING,
                                         MAX_NBRS_PER_VERT_FOR_TET_CLIPPING>(
                                           tet,
                                           geomPiece,
                                           EPS,
                                           tryFixOrientation);

          if(poly.numVertices() >= 4)
          {
            // Poly is valid
            double volume = poly.volume();
            SLIC_ASSERT(volume >= 0);
            auto cellId = tetId / TETS_PER_HEXAHEDRON;
            RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, volume);
          }
        });
    }

  }  // end of computeClipVolumes3DTets() function

  void getLabelCounts(axom::ArrayView<const LabelType> labels,
                      axom::IndexType& inCount,
                      axom::IndexType& onCount,
                      axom::IndexType& outCount) override
  {
    AXOM_ANNOTATE_SCOPE("GeometryClipper::getLabelCounts");
    using ReducePolicy = typename axom::execution_space<ExecSpace>::reduce_policy;
    using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
    RAJA::ReduceSum<ReducePolicy, axom::IndexType> inSum(0);
    RAJA::ReduceSum<ReducePolicy, axom::IndexType> onSum(0);
    RAJA::ReduceSum<ReducePolicy, axom::IndexType> outSum(0);
    RAJA::forall<LoopPolicy>(
      RAJA::RangeSegment(0, labels.size()),
      AXOM_LAMBDA(axom::IndexType cellId) {
        const auto& label = labels[cellId];
        if(label == GeometryClipperStrategy::LABEL_OUT)
        {
          outSum += 1;
        }
        else if(label == GeometryClipperStrategy::LABEL_IN)
        {
          inSum += 1;
        }
        else
        {
          onSum += 1;
        }
      });
    inCount = static_cast<axom::IndexType>(inSum.get());
    onCount = static_cast<axom::IndexType>(onSum.get());
    outCount = static_cast<axom::IndexType>(outSum.get());
  }
private:
  static constexpr int MAX_VERTS_FOR_TET_CLIPPING = 32;
  static constexpr int MAX_NBRS_PER_VERT_FOR_TET_CLIPPING = 8;
  static constexpr int MAX_VERTS_FOR_OCT_CLIPPING = 32;
  static constexpr int MAX_NBRS_PER_VERT_FOR_OCT_CLIPPING = 8;
};

}  // end namespace detail
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_GEOMETRYCLIPPERDELEGATEEXEC_HPP_
