// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_GEOMETRYCLIPPERDELEGATEEXEC_HPP_
#define AXOM_GEOMETRYCLIPPERDELEGATEEXEC_HPP_

#include "axom/quest/GeometryClipper.hpp"
#include "axom/spin/BVH.hpp"

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

  GeometryClipperDelegateExec(GeometryClipper& delegator)
    : GeometryClipper::Delegate(delegator)
  { }

  void setCleanVolumeOverlaps(
    const axom::ArrayView<GeometryClipperStrategy::LabelType>& labels,
    axom::ArrayView<double> ovlap)
  {
    const axom::IndexType cellCount = getShapeeMesh().getCellCount();
    SLIC_ASSERT(labels.size() == cellCount);
    SLIC_ASSERT(ovlap.size() == cellCount);

    auto cellVolumes = getShapeeMesh().getCellVolumes();

    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        auto& l = labels[i];
        if(l == 0)
        {
          ovlap[i] = cellVolumes[i];
        }
        else if(l == 2)
        {
          ovlap[i] = 0;
        }
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
                                   axom::Array<axom::IndexType>& unlabeledCells)
  {
    using ScanPolicy = typename axom::execution_space<ExecSpace>::loop_policy;

    const axom::IndexType cellCount = labels.size();

    axom::Array<axom::IndexType> tmpLabels(labels.shape(),
                                           labels.getAllocatorID());
    auto tmpLabelsView = tmpLabels.view();
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType ci) { tmpLabelsView[ci] = labels[ci] == 1; });

    RAJA::inclusive_scan_inplace<ScanPolicy>(
      RAJA::make_span(tmpLabels.data(), tmpLabels.size()),
      RAJA::operators::plus<axom::IndexType> {});

    axom::IndexType unlabeledCount;
    axom::copy(&unlabeledCount, &tmpLabels.back(), sizeof(unlabeledCount));

    // Re-allocate output if it's too small or has wrong allocator.
    if(unlabeledCells.size() < unlabeledCount ||
       unlabeledCells.getAllocatorID() != labels.getAllocatorID())
    {
      unlabeledCells =
        axom::Array<axom::IndexType> {axom::ArrayOptions::Uninitialized(),
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
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        if(tmpLabelsView[i] != tmpLabelsView[i - 1])
        {
          unlabeledCellsView[tmpLabelsView[i - 1]] = i;
        }
      });
  }

  void computeClipVolumes3D(axom::ArrayView<double> ovlap)
  {
    AXOM_ANNOTATE_SCOPE("IntersectionShaper::runShapeQueryImpl");

    using BoundingBoxType = primal::BoundingBox<double, 3>;

    ShapeeMesh& shapeeMesh = m_delegator.getShapeeMesh();

    const int allocId = shapeeMesh.getAllocatorId();

    const IndexType cellCount = shapeeMesh.getCellCount();

    constexpr int NUM_TETS_PER_HEX = 24;

    SLIC_INFO(axom::fmt::format("{:-^80}",
                                " Inserting shapes' bounding boxes into BVH "));

#if 0
    axom::ArrayView<ShapeType> discretizedGeometryView = discretizedGeometry.view();
#endif

    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    const bool useOcts = getGeometryClipperStrategy().getShapeAsOcts(shapeeMesh, geomAsOcts);
    const bool useTets = getGeometryClipperStrategy().getShapeAsTets(shapeeMesh, geomAsTets);
    SLIC_ASSERT(useTets != useOcts);
    SLIC_ASSERT(useOcts || geomAsOcts.empty());
    SLIC_ASSERT(useTets || geomAsTets.empty());

    auto geomTetsView = geomAsTets.view();
    auto geomOctsView = geomAsOcts.view();

    // Generate the BVH tree over the shape's discretized geometry
    // Axis-aligned bounding boxes
    const axom::IndexType bbCount =
      useTets ? geomAsTets.size() : geomAsOcts.size();
    axom::Array<BoundingBoxType> pieceBbs(bbCount, bbCount, allocId);
    axom::ArrayView<BoundingBoxType> pieceBbsView = pieceBbs.view();

    // Get the bounding boxes for the shapes
    if(useTets)
    {
      axom::for_all<ExecSpace>(
        bbCount,
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] =
            primal::compute_bounding_box<double, 3>(geomTetsView[i]);
        });
    }
    else
    {
      axom::for_all<ExecSpace>(
        bbCount,
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] =
            primal::compute_bounding_box<double, 3>(geomOctsView[i]);
        });
    }

    // Insert shapes' Bounding Boxes into BVH.
    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(pieceBbsView, bbCount);

    SLIC_INFO(axom::fmt::format("{:-^80}", " Querying the BVH tree "));

    axom::ArrayView<const BoundingBoxType> hex_bbs_device_view =
      shapeeMesh.getCellBoundingBoxes();

    // Find which shape bounding boxes intersect hexahedron bounding boxes
    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      " Finding shape candidates for each hexahedral element "));

    axom::Array<IndexType> offsets(cellCount, cellCount, allocId);
    axom::Array<IndexType> counts(cellCount, cellCount, allocId);
    axom::Array<IndexType> candidates;
    AXOM_ANNOTATE_BEGIN("bvh.findBoundingBoxes");
    bvh.findBoundingBoxes(offsets,
                          counts,
                          candidates,
                          cellCount,
                          hex_bbs_device_view);
    AXOM_ANNOTATE_END("bvh.findBoundingBoxes");

    // Get the total number of candidates
    using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

    const auto counts_device_view = counts.view();
    AXOM_ANNOTATE_BEGIN("populate totalCandidates");
    RAJA::ReduceSum<REDUCE_POL, int> totalCandidates(0);
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        totalCandidates += counts_device_view[i];
      });
    AXOM_ANNOTATE_END("populate totalCandidates");

    AXOM_ANNOTATE_BEGIN("allocate scratch space");
    // Initialize hexahedron indices and shape candidates
    AXOM_ANNOTATE_BEGIN("allocate hex_indices");
    axom::Array<IndexType> hex_indices_device(
      totalCandidates.get() * NUM_TETS_PER_HEX,
      totalCandidates.get() * NUM_TETS_PER_HEX,
      allocId);
    AXOM_ANNOTATE_END("allocate hex_indices");
    auto hex_indices_device_view = hex_indices_device.view();

    AXOM_ANNOTATE_BEGIN("allocate shape_candidates");
    axom::Array<IndexType> shape_candidates_device(
      totalCandidates.get() * NUM_TETS_PER_HEX,
      totalCandidates.get() * NUM_TETS_PER_HEX,
      allocId);
    AXOM_ANNOTATE_END("allocate shape_candidates");
    auto shape_candidates_device_view = shape_candidates_device.view();

    // Tetrahedrons from hexes (24 for each hex)
    auto tets_from_hexes_device_view = shapeeMesh.getCellsAsTets();

    // Index into 'tets'
    AXOM_ANNOTATE_BEGIN("allocate tet_indices_device");
    axom::Array<IndexType> tet_indices_device(
      totalCandidates.get() * NUM_TETS_PER_HEX,
      totalCandidates.get() * NUM_TETS_PER_HEX,
      allocId);
    AXOM_ANNOTATE_END("allocate tet_indices_device");
    auto tet_indices_device_view = tet_indices_device.view();
    AXOM_ANNOTATE_END("allocate scratch space");

    // New total number of candidates after omitting degenerate shapes
    AXOM_ANNOTATE_BEGIN("newTotalCandidates memory");
    IndexType totalCandidatesCount = 0;
    IndexType* totalCandidatesCountPtr = &totalCandidatesCount;
#if 1
    if(!axom::execution_space<ExecSpace>::usesMemorySpace(MemorySpace::Dynamic))
    {
      // Use temporary space compatible with runtime policy.
      totalCandidatesCountPtr = axom::allocate<IndexType>(1, allocId);
    }
#else
    axom::Array<IndexType> newTotalCandidates_host(1, 1, host_allocator);
    newTotalCandidates_host[0] = 0;
    axom::Array<IndexType> newTotalCandidates_device =
      axom::Array<IndexType>(newTotalCandidates_host, allocId);
    auto newTotalCandidates_device_view = newTotalCandidates_device.view();
#endif
    AXOM_ANNOTATE_END("newTotalCandidates memory");

    SLIC_INFO(
      axom::fmt::format("{:-^80}",
                        " Creating an array of candidate pairs for shaping "));

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
              IndexType idx = RAJA::atomicAdd<ATOMIC_POL>(totalCandidatesCountPtr,
                                                          IndexType {1});
              hex_indices_device_view[idx] = i;
              shape_candidates_device_view[idx] = shapeIdx;
              tet_indices_device_view[idx] = i * NUM_TETS_PER_HEX + k;
            }
          }
        });
    }

    // ovlap.fill(0.0);
    axom::for_all<ExecSpace>(
      ovlap.size(),
      AXOM_LAMBDA(axom::IndexType i) { ovlap[i] = 0.0; });

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      " Calculating element overlap volume from each tet-shape pair "));

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = true;

    {
      AXOM_ANNOTATE_SCOPE("clipLoop");
      // Copy calculated total back to host if needed
#if 1
      if(totalCandidatesCountPtr != &totalCandidatesCount)
      {
        axom::copy(&totalCandidatesCount,
                   totalCandidatesCountPtr,
                   sizeof(IndexType));
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

            const PolyhedronType poly =
              primal::clip(geomTetsView[shapeIndex],
                           tets_from_hexes_device_view[tetIndex],
                           EPS,
                           tryFixOrientation);

            // Poly is valid
            if(poly.numVertices() >= 4)
            {
              // Workaround - intermediate volume variable needed for
              // CUDA Pro/E test case correctness
              double volume = poly.volume();
              RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + index, volume);
            }
          });
      }
      else
      {
        axom::for_all<ExecSpace>(
          totalCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            const int index = hex_indices_device_view[i];
            const int shapeIndex = shape_candidates_device_view[i];
            const int tetIndex = tet_indices_device_view[i];

            const PolyhedronType poly =
              primal::clip(geomOctsView[shapeIndex],
                           tets_from_hexes_device_view[tetIndex],
                           EPS,
                           tryFixOrientation);

            // Poly is valid
            if(poly.numVertices() >= 4)
            {
              // Workaround - intermediate volume variable needed for
              // CUDA Pro/E test case correctness
              double volume = poly.volume();
              RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + index, volume);
            }
          });
      }
#else
      axom::Array<IndexType> newTotalCandidates_calc_host =
        axom::Array<IndexType>(newTotalCandidates_device, host_allocator);

      axom::for_all<ExecSpace>(
        newTotalCandidates_calc_host[0],  // Number of candidates found.
        AXOM_LAMBDA(axom::IndexType i) {
          const int index = hex_indices_device_view[i];
          const int shapeIndex = shape_candidates_device_view[i];
          const int tetIndex = tet_indices_device_view[i];

          const PolyhedronType poly =
            primal::clip(discretizedGeometryView[shapeIndex],
                         tets_from_hexes_device_view[tetIndex],
                         EPS,
                         tryFixOrientation);

          // Poly is valid
          if(poly.numVertices() >= 4)
          {
            // Workaround - intermediate volume variable needed for
            // CUDA Pro/E test case correctness
            double volume = poly.volume();
            RAJA::atomicAdd<ATOMIC_POL>(ovlapView.data() + index, volume);
          }
        });
#endif
    }

    if(totalCandidatesCountPtr != &totalCandidatesCount)
    {
      axom::deallocate(totalCandidatesCountPtr);
    }
  }  // end of computeClipVolumes3D() function

  void computeClipVolumes3D(const axom::ArrayView<axom::IndexType>& cellIndices,
                            axom::ArrayView<double> ovlap)

  {
    AXOM_UNUSED_VAR(ovlap);
    AXOM_UNUSED_VAR(cellIndices);
    AXOM_ANNOTATE_SCOPE("GeometryClipper::clip_volumes_3D");

    using BoundingBoxType = primal::BoundingBox<double, 3>;

    ShapeeMesh& shapeeMesh = m_delegator.getShapeeMesh();

    const int allocId = shapeeMesh.getAllocatorId();

    const IndexType cellCount = shapeeMesh.getCellCount();

    constexpr int NUM_TETS_PER_HEX = 24;

    SLIC_INFO(axom::fmt::format("{:-^80}",
                                " Inserting shapes' bounding boxes into BVH "));

    axom::Array<axom::primal::Tetrahedron<double, 3>> geomAsTets;
    axom::Array<axom::primal::Octahedron<double, 3>> geomAsOcts;
    const bool useOcts = getGeometryClipperStrategy().getShapeAsOcts(shapeeMesh, geomAsOcts);
    const bool useTets = getGeometryClipperStrategy().getShapeAsTets(shapeeMesh, geomAsTets);
    SLIC_ASSERT(useTets != useOcts);
    SLIC_ASSERT(useTets + useOcts == 1);

    auto geomTetsView = geomAsTets.view();
    auto geomOctsView = geomAsOcts.view();

    // Generate the BVH tree over the shape's discretized geometry
    // axis-aligned bounding boxes.  "pieces" refers to tets or octs.
    const axom::IndexType bbCount =
      useTets ? geomAsTets.size() : geomAsOcts.size();
    axom::Array<BoundingBoxType> pieceBbs(bbCount, bbCount, allocId);
    axom::ArrayView<BoundingBoxType> pieceBbsView = pieceBbs.view();

    // Get the bounding boxes for the shapes
    if(useTets)
    {
      axom::for_all<ExecSpace>(
        bbCount,
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] =
            primal::compute_bounding_box<double, 3>(geomTetsView[i]);
        });
    }
    else
    {
      axom::for_all<ExecSpace>(
        bbCount,
        AXOM_LAMBDA(axom::IndexType i) {
          pieceBbsView[i] =
            primal::compute_bounding_box<double, 3>(geomOctsView[i]);
        });
    }

    // Insert shapes' Bounding Boxes into BVH.
    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(pieceBbsView, bbCount);

    SLIC_INFO(axom::fmt::format("{:-^80}", " Querying the BVH tree "));

    axom::ArrayView<const BoundingBoxType> hex_bbs_device_view =
      shapeeMesh.getCellBoundingBoxes();

    // Find which shape bounding boxes intersect hexahedron bounding boxes
    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      " Finding shape candidates for each hexahedral element "));

    axom::Array<IndexType> offsets(cellCount, cellCount, allocId);
    axom::Array<IndexType> counts(cellCount, cellCount, allocId);
    axom::Array<IndexType> candidates;
    AXOM_ANNOTATE_BEGIN("bvh.findBoundingBoxes");
    bvh.findBoundingBoxes(offsets,
                          counts,
                          candidates,
                          cellCount,
                          hex_bbs_device_view);
    AXOM_ANNOTATE_END("bvh.findBoundingBoxes");

    // Get the total number of candidates
    using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

    const auto counts_device_view = counts.view();
    AXOM_ANNOTATE_BEGIN("populate totalCandidates");
    RAJA::ReduceSum<REDUCE_POL, int> totalCandidates(0);
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        totalCandidates += counts_device_view[i];
      });
    AXOM_ANNOTATE_END("populate totalCandidates");

    AXOM_ANNOTATE_BEGIN("allocate scratch space");
    // Initialize hexahedron indices and shape candidates
    AXOM_ANNOTATE_BEGIN("allocate hex_indices");
    axom::Array<IndexType> hex_indices_device(
      totalCandidates.get() * NUM_TETS_PER_HEX,
      totalCandidates.get() * NUM_TETS_PER_HEX,
      allocId);
    AXOM_ANNOTATE_END("allocate hex_indices");
    auto hex_indices_device_view = hex_indices_device.view();

    AXOM_ANNOTATE_BEGIN("allocate shape_candidates");
    axom::Array<IndexType> shape_candidates_device(
      totalCandidates.get() * NUM_TETS_PER_HEX,
      totalCandidates.get() * NUM_TETS_PER_HEX,
      allocId);
    AXOM_ANNOTATE_END("allocate shape_candidates");
    auto shape_candidates_device_view = shape_candidates_device.view();

    // Tetrahedrons from hexes (24 for each hex)
    auto tets_from_hexes_device_view = shapeeMesh.getCellsAsTets();

    // Index into 'tets'
    AXOM_ANNOTATE_BEGIN("allocate tet_indices_device");
    axom::Array<IndexType> tet_indices_device(
      totalCandidates.get() * NUM_TETS_PER_HEX,
      totalCandidates.get() * NUM_TETS_PER_HEX,
      allocId);
    AXOM_ANNOTATE_END("allocate tet_indices_device");
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
    }
    AXOM_ANNOTATE_END("newTotalCandidates memory");

    SLIC_INFO(
      axom::fmt::format("{:-^80}",
                        " Creating an array of candidate pairs for shaping "));

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
              IndexType idx = RAJA::atomicAdd<ATOMIC_POL>(totalCandidatesCountPtr,
                                                          IndexType {1});
              hex_indices_device_view[idx] = i;
              shape_candidates_device_view[idx] = shapeIdx;
              tet_indices_device_view[idx] = i * NUM_TETS_PER_HEX + k;
            }
          }
        });
    }

    // ovlap.fill(0.0);
    axom::for_all<ExecSpace>(
      ovlap.size(),
      AXOM_LAMBDA(axom::IndexType i) { ovlap[i] = 0.0; });

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      " Calculating element overlap volume from each tet-shape pair "));

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = true;

    {
      AXOM_ANNOTATE_SCOPE("clipLoop");
      // Copy calculated total back to host if needed
      if(totalCandidatesCountPtr != &totalCandidatesCount)
      {
        axom::copy(&totalCandidatesCount,
                   totalCandidatesCountPtr,
                   sizeof(IndexType));
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

            const PolyhedronType poly =
              primal::clip(geomTetsView[shapeIndex],
                           tets_from_hexes_device_view[tetIndex],
                           EPS,
                           tryFixOrientation);

            // Poly is valid
            if(poly.numVertices() >= 4)
            {
              // Workaround - intermediate volume variable needed for
              // CUDA Pro/E test case correctness
              double volume = poly.volume();
              RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + index, volume);
            }
          });
      }
      else
      {
        axom::for_all<ExecSpace>(
          totalCandidatesCount,
          AXOM_LAMBDA(axom::IndexType i) {
            const int index = hex_indices_device_view[i];
            const int shapeIndex = shape_candidates_device_view[i];
            const int tetIndex = tet_indices_device_view[i];

            const PolyhedronType poly =
              primal::clip(geomOctsView[shapeIndex],
                           tets_from_hexes_device_view[tetIndex],
                           EPS,
                           tryFixOrientation);

            // Poly is valid
            if(poly.numVertices() >= 4)
            {
              // Workaround - intermediate volume variable needed for
              // CUDA Pro/E test case correctness
              double volume = poly.volume();
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

};

}  // end namespace detail
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_GEOMETRYCLIPPERDELEGATEEXEC_HPP_
