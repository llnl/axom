// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_GEOMETRYCLIPPERDELEGATEEXEC_HPP_
#define AXOM_GEOMETRYCLIPPERDELEGATEEXEC_HPP_

#ifndef AXOM_USE_RAJA
  #error "quest::GeometryClipper requires RAJA."
#endif

#include "axom/quest/GeometryClipper.hpp"
#include "axom/spin/BVH.hpp"
#include "axom/core/WhereMacro.hpp"
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

    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        auto& l = labels[i];
        ovlap[i] = l == 0 ? cellVolumes[i] : 0.0;
      });

    return;
  }

  /*!
    @brief Make an list of indices where labels have value 1.

    1. Generate tmpLabels, having a value of 1 for unlabeled cells and zero elsewhere.
    2. Inclusive scan on tmpLabels to generate values that step up at unlabeled cells.
    3. Find unlabeled cells by seeing where tmpLabels changes values.
       (Handle first cell separately, then loop from second cell on.)
  */
  void collectUnlabeledCellIndices(const axom::ArrayView<LabelType>& labels,
                                   axom::Array<axom::IndexType>& unlabeledCells) override
  {
    using ScanPolicy = typename axom::execution_space<ExecSpace>::loop_policy;

    const axom::IndexType labelCount = labels.size();

    axom::Array<axom::IndexType> tmpLabels(labels.shape(), labels.getAllocatorID());
    auto tmpLabelsView = tmpLabels.view();
    axom::for_all<ExecSpace>(
      labelCount,
      AXOM_LAMBDA(axom::IndexType ci) { tmpLabelsView[ci] = labels[ci] == 1; });

    RAJA::inclusive_scan_inplace<ScanPolicy>(RAJA::make_span(tmpLabels.data(), tmpLabels.size()),
                                             RAJA::operators::plus<axom::IndexType> {});

    axom::IndexType unlabeledCount;
    axom::copy(&unlabeledCount, &tmpLabels.back(), sizeof(unlabeledCount));

    // Re-allocate output if it's too small or has wrong allocator.
    if(unlabeledCells.size() < unlabeledCount ||
       unlabeledCells.getAllocatorID() != labels.getAllocatorID())
    {
      unlabeledCells = axom::Array<axom::IndexType> {axom::ArrayOptions::Uninitialized(),
                                                     unlabeledCount,
                                                     unlabeledCount,
                                                     labels.getAllocatorID()};
    }
    auto unlabeledCellsView = unlabeledCells.view();

    LabelType firstCellLabel = '\0';
    axom::copy(&firstCellLabel, &labels[0], sizeof(firstCellLabel));
    if(firstCellLabel == 1)
    {
      axom::IndexType zero = 0;
      axom::copy(&unlabeledCells[0], &zero, sizeof(zero));
    }

    axom::for_all<ExecSpace>(
      1,
      labelCount,
      AXOM_LAMBDA(axom::IndexType i) {
        if(tmpLabelsView[i] != tmpLabelsView[i - 1])
        {
          unlabeledCellsView[tmpLabelsView[i - 1]] = i;
        }
      });
  }

  void computeClipVolumes3D(axom::ArrayView<double> ovlap) override
  {
    AXOM_ANNOTATE_SCOPE("GeometryClipper::computeClipVolumes3D");

    using BoundingBoxType = primal::BoundingBox<double, 3>;

    ShapeeMesh& shapeeMesh = getDelegator().getShapeeMesh();

    const int allocId = shapeeMesh.getAllocatorId();

    const IndexType cellCount = shapeeMesh.getCellCount();

    constexpr int NUM_TETS_PER_HEX = 24;

    SLIC_INFO(axom::fmt::format(
      "GeometryClipper::computeClipVolumes3D: Getting discrete geometry for shape '{}'",
      getStrategy().name()));

    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    const bool useOcts = getStrategy().getGeometryAsOcts(shapeeMesh, geomAsOcts);
    const bool useTets = getStrategy().getGeometryAsTets(shapeeMesh, geomAsTets);
    SLIC_ASSERT(useTets != useOcts);
    SLIC_ASSERT(useOcts || geomAsOcts.empty());
    SLIC_ASSERT(useTets || geomAsTets.empty());

    auto geomTetsView = geomAsTets.view();
    auto geomOctsView = geomAsOcts.view();

    SLIC_INFO(axom::fmt::format("{:-^80}", " Inserting shapes' bounding boxes into BVH "));

    // Generate the BVH tree over the shape's discretized geometry
    // Axis-aligned bounding boxes
    const axom::IndexType bbCount = useTets ? geomAsTets.size() : geomAsOcts.size();
    axom::Array<BoundingBoxType> pieceBbs(bbCount, bbCount, allocId);
    axom::ArrayView<BoundingBoxType> pieceBbsView = pieceBbs.view();

    // Get the bounding boxes for the shapes
    if(useTets)
    {
      axom::for_all<ExecSpace>(
        bbCount,
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomTetsView[i]);
        });
    }
    else
    {
      axom::for_all<ExecSpace>(
        bbCount,
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomOctsView[i]);
        });
    }

    // Insert shapes' Bounding Boxes into BVH.
    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(pieceBbsView, bbCount);

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
    using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

    const auto counts_device_view = counts.view();
    AXOM_ANNOTATE_BEGIN("populate totalCandidates");
    RAJA::ReduceSum<REDUCE_POL, int> totalCandidates(0);
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) { totalCandidates += counts_device_view[i]; });
    AXOM_ANNOTATE_END("populate totalCandidates");

    AXOM_ANNOTATE_BEGIN("allocate scratch space");
    // Initialize hexahedron indices and shape candidates
    axom::Array<IndexType> hex_indices_device(totalCandidates.get() * NUM_TETS_PER_HEX,
                                              totalCandidates.get() * NUM_TETS_PER_HEX,
                                              allocId);
    auto hex_indices_device_view = hex_indices_device.view();

    axom::Array<IndexType> shape_candidates_device(totalCandidates.get() * NUM_TETS_PER_HEX,
                                                   totalCandidates.get() * NUM_TETS_PER_HEX,
                                                   allocId);
    auto shape_candidates_device_view = shape_candidates_device.view();

    // Tetrahedrons from hexes (24 for each hex)
    auto tets_from_hexes_device_view = shapeeMesh.getCellsAsTets();

    // Index into 'tets'
    axom::Array<IndexType> tet_indices_device(totalCandidates.get() * NUM_TETS_PER_HEX,
                                              totalCandidates.get() * NUM_TETS_PER_HEX,
                                              allocId);
    auto tet_indices_device_view = tet_indices_device.view();
    AXOM_ANNOTATE_END("allocate scratch space");

    // New total number of candidates after omitting degenerate shapes
    AXOM_ANNOTATE_BEGIN("newTotalCandidates memory");
    IndexType totalCandidatesCount = 0;
    IndexType* totalCandidatesCountPtr = &totalCandidatesCount;
    if(!axom::execution_space<ExecSpace>::usesMemorySpace(MemorySpace::Dynamic))
    {
      // Use temporary space compatible with runtime policy.
      totalCandidatesCountPtr = axom::allocate<IndexType>(1, allocId);
      axom::copy(totalCandidatesCountPtr, &totalCandidatesCount, sizeof(totalCandidatesCount));
    }
    AXOM_ANNOTATE_END("newTotalCandidates memory");

    const auto offsets_device_view = offsets.view();
    const auto candidates_device_view = candidates.view();
    {
      AXOM_ANNOTATE_SCOPE("init_candidates");
      axom::for_all<ExecSpace>(
        cellCount,
        AXOM_LAMBDA(axom::IndexType i) {
          for(int j = 0; j < counts_device_view[i]; j++)
          {
            int shapeIdx = candidates_device_view[offsets_device_view[i] + j];

            for(int k = 0; k < NUM_TETS_PER_HEX; k++)
            {
              IndexType idx = RAJA::atomicAdd<ATOMIC_POL>(totalCandidatesCountPtr, IndexType {1});
              hex_indices_device_view[idx] = i;
              shape_candidates_device_view[idx] = shapeIdx;
              tet_indices_device_view[idx] = i * NUM_TETS_PER_HEX + k;
            }
          }
        });
    }

    axom::for_all<ExecSpace>(
      ovlap.size(),
      AXOM_LAMBDA(axom::IndexType i) { ovlap[i] = 0.0; });

    SLIC_INFO(axom::fmt::format(
      "Running clip loop on {} candidate tets for of all {} hexes in the mesh",
      totalCandidatesCount,
      cellCount));

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = false;

    {
      AXOM_ANNOTATE_SCOPE("GeometryClipper::clipLoop");
      // Copy calculated total back to host if needed
      if(totalCandidatesCountPtr != &totalCandidatesCount)
      {
        axom::copy(&totalCandidatesCount, totalCandidatesCountPtr, sizeof(IndexType));
      }

      using PolyhedronType = primal::Polyhedron<double, 3>;

      if(useTets)
      {
        axom::for_all<ExecSpace>(
          totalCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            const int index = hex_indices_device_view[i];
            const int shapeIndex = shape_candidates_device_view[i];
            const int tetIndex = tet_indices_device_view[i];

            const PolyhedronType poly = primal::clip(geomTetsView[shapeIndex],
                                                     tets_from_hexes_device_view[tetIndex],
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
          totalCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            const int index = hex_indices_device_view[i];
            const int shapeIndex = shape_candidates_device_view[i];
            const int tetIndex = tet_indices_device_view[i];

            const PolyhedronType poly = primal::clip(geomOctsView[shapeIndex],
                                                     tets_from_hexes_device_view[tetIndex],
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

    if(totalCandidatesCountPtr != &totalCandidatesCount)
    {
      axom::deallocate(totalCandidatesCountPtr);
    }
  }  // end of computeClipVolumes3D() function

  void computeClipVolumes3D(const axom::ArrayView<axom::IndexType>& cellIndices,
                            axom::ArrayView<double> ovlap) override

  {
    AXOM_UNUSED_VAR(ovlap);
    AXOM_ANNOTATE_SCOPE("GeometryClipper::computeClipVolumes3D:limited");

    using BoundingBoxType = primal::BoundingBox<double, 3>;

    ShapeeMesh& shapeeMesh = getDelegator().getShapeeMesh();

    const int allocId = shapeeMesh.getAllocatorId();

    const IndexType cellCount = shapeeMesh.getCellCount();

    constexpr int NUM_TETS_PER_HEX = 24;

    SLIC_INFO(axom::fmt::format(
      "GeometryClipper::computeClipVolumes3D: Getting discrete geometry for shape '{}'",
      getStrategy().name()));

    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    const bool useOcts = getStrategy().getGeometryAsOcts(shapeeMesh, geomAsOcts);
    const bool useTets = getStrategy().getGeometryAsTets(shapeeMesh, geomAsTets);
    SLIC_ASSERT(useTets != useOcts);
    SLIC_ASSERT(useOcts || geomAsOcts.empty());
    SLIC_ASSERT(useTets || geomAsTets.empty());

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
        bbCount,
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomTetsView[i]);
        });
    }
    else
    {
      axom::for_all<ExecSpace>(
        bbCount,
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] = primal::compute_bounding_box<double, 3>(geomOctsView[i]);
        });
    }

    // Insert shapes' Bounding Boxes into BVH.
    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(pieceBbsView, bbCount);

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
    using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

    const auto counts_device_view = counts.view();
    AXOM_ANNOTATE_BEGIN("populate totalCandidates");
    RAJA::ReduceSum<REDUCE_POL, int> totalCandidates(0);
    axom::for_all<ExecSpace>(
      limitedCellCount,
      AXOM_LAMBDA(axom::IndexType i) { totalCandidates += counts_device_view[i]; });
    AXOM_ANNOTATE_END("populate totalCandidates");
assert(totalCandidates.get() == candidates.size());

    AXOM_ANNOTATE_BEGIN("allocate scratch space");
    // Initialize hexahedron indices and shape candidates
    axom::Array<IndexType> hex_indices_device(totalCandidates.get() * NUM_TETS_PER_HEX,
                                              totalCandidates.get() * NUM_TETS_PER_HEX,
                                              allocId);
    auto hex_indices_device_view = hex_indices_device.view();

    axom::Array<IndexType> shape_candidates_device(totalCandidates.get() * NUM_TETS_PER_HEX,
                                                   totalCandidates.get() * NUM_TETS_PER_HEX,
                                                   allocId);
    auto shape_candidates_device_view = shape_candidates_device.view();

    // Tetrahedrons from hexes (24 for each hex)
    auto tets_from_hexes_device_view = shapeeMesh.getCellsAsTets();

    // Index into 'tets'
    axom::Array<IndexType> tet_indices_device(totalCandidates.get() * NUM_TETS_PER_HEX,
                                              totalCandidates.get() * NUM_TETS_PER_HEX,
                                              allocId);
    auto tet_indices_device_view = tet_indices_device.view();
    AXOM_ANNOTATE_END("allocate scratch space");

    // New total number of candidates after omitting degenerate shapes
    AXOM_ANNOTATE_BEGIN("newTotalCandidates memory");
    IndexType totalCandidatesCount = 0;
    IndexType* totalCandidatesCountPtr = &totalCandidatesCount;
    if(!axom::execution_space<ExecSpace>::usesMemorySpace(MemorySpace::Dynamic))
    {
      // Use temporary space compatible with runtime policy.
      totalCandidatesCountPtr = axom::allocate<IndexType>(1, allocId);
      axom::copy(totalCandidatesCountPtr, &totalCandidatesCount, sizeof(IndexType));
    }
    AXOM_ANNOTATE_END("newTotalCandidates memory");

    const auto offsets_device_view = offsets.view();
    const auto candidates_device_view = candidates.view();
    {
      AXOM_ANNOTATE_SCOPE("init_candidates");
      axom::for_all<ExecSpace>(
        limitedCellCount,
        AXOM_LAMBDA(axom::IndexType i) {
          for(int j = 0; j < counts_device_view[i]; j++)
          {
            int shapeIdx = candidates_device_view[offsets_device_view[i] + j];

            for(int k = 0; k < NUM_TETS_PER_HEX; k++)
            {
              IndexType idx = RAJA::atomicAdd<ATOMIC_POL>(totalCandidatesCountPtr, IndexType {1});
              hex_indices_device_view[idx] = i;
              shape_candidates_device_view[idx] = shapeIdx;
              tet_indices_device_view[idx] = i * NUM_TETS_PER_HEX + k;
            }
          }
        });
    }

    SLIC_INFO(axom::fmt::format(
      "Running clip loop on {} candidate tets for the select {} hexes of the full {} cells",
      totalCandidatesCount,
      cellIndices.size(),
      cellCount));

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = false;

    {
      AXOM_ANNOTATE_SCOPE("GeometryClipper::clipLoop_limited");
      totalCandidatesCount = NUM_TETS_PER_HEX*candidates.size();
#if 0
      // Verifying: this should always pass.
      if(totalCandidatesCountPtr != &totalCandidatesCount)
      {
        axom::copy(&totalCandidatesCount, totalCandidatesCountPtr, sizeof(IndexType));
      }
      SLIC_ASSERT(totalCandidatesCount == NUM_TETS_PER_HEX*candidates.size());
#endif

      using PolyhedronType = primal::Polyhedron<double, 3>;

      if(useTets)
      {
        axom::for_all<ExecSpace>(
          totalCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            int index = hex_indices_device_view[i]; // index into limited mesh hex array
            index = cellIndices[index]; // Now, it indexes into the full hex array.

            const int shapeIndex = shape_candidates_device_view[i]; // index into pieces array
            int tetIndex = tet_indices_device_view[i]; // index into BVH results, implicit because BVH results specify hexes, not tets.
            int tetIndex1 = tetIndex/NUM_TETS_PER_HEX;
            int tetIndex2 = tetIndex%NUM_TETS_PER_HEX;
            tetIndex = cellIndices[tetIndex1]*NUM_TETS_PER_HEX + tetIndex2; // Now it indexes into the full tets-from-hexes array.

            const PolyhedronType poly = primal::clip(geomTetsView[shapeIndex],
                                                     tets_from_hexes_device_view[tetIndex],
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
          totalCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            int index = hex_indices_device_view[i]; // index into limited mesh hex array
            index = cellIndices[index]; // Now, it indexes into the full hex array.

            const int shapeIndex = shape_candidates_device_view[i]; // index into pieces array
            int tetIndex = tet_indices_device_view[i]; // index into BVH results, implicit because BVH results specify hexes, not tets.
            int tetIndex1 = tetIndex/NUM_TETS_PER_HEX;
            int tetIndex2 = tetIndex%NUM_TETS_PER_HEX;
            tetIndex = cellIndices[tetIndex1]*NUM_TETS_PER_HEX + tetIndex2; // Now it indexes into the full tets-from-hexes array.

            const PolyhedronType poly = primal::clip(geomOctsView[shapeIndex],
                                                     tets_from_hexes_device_view[tetIndex],
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

    if(totalCandidatesCountPtr != &totalCandidatesCount)
    {
      axom::deallocate(totalCandidatesCountPtr);
    }
  }  // end of computeClipVolumes3D() function

  void getLabelCounts(axom::ArrayView<const LabelType> labels,
                      axom::IndexType& inCount,
                      axom::IndexType& onCount,
                      axom::IndexType& outCount) override
  {
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
};

}  // end namespace detail
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_GEOMETRYCLIPPERDELEGATEEXEC_HPP_
