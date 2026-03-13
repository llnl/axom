// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_MATSET_SLICER_HPP
#define AXOM_BUMP_MATSET_SLICER_HPP

#include "axom/core.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/FieldSlicer.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"

#include <conduit.hpp>

namespace axom
{
namespace bump
{

/*!
 * \brief Slices the input matset view and outputs a new matset (unibuffer flavor).
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam MatsetView The matset view type that wraps the Blueprint matset.
 */
template <typename ExecSpace, typename MatsetView>
class MatsetSlicer
{
public:
  using SelectedZonesView = axom::ArrayView<axom::IndexType>;

  /*!
   * \brief Constructor.
   */
  MatsetSlicer(const MatsetView &matsetView)
    : m_matsetView(matsetView)
    , m_allocator_id(axom::execution_space<ExecSpace>::allocatorID())
  { }

  /*!
   * \brief Set the allocator id to use when allocating memory.
   *
   * \param allocator_id The allocator id to use when allocating memory.
   */
  void setAllocatorID(int allocator_id) { m_allocator_id = allocator_id; }

  /*!
   * \brief Get the allocator id to use when allocating memory.
   *
   * \return The allocator id to use when allocating memory.
   */
  int getAllocatorID() const { return m_allocator_id; }

  /*!
   * \brief Slice the input matset and output a new matset.
   *
   * \param matsetView A view that wraps the input matset.
   * \param slice Slice data that contains the zone ids that we're extracting from the matset.
   * \param n_matset The input matset.
   * \param[out] n_newMatset The output matset.
   */
  void execute(const SliceData &slice, const conduit::Node &n_matset, conduit::Node &n_newMatset)
  {
    using MatsetIndex = typename MatsetView::IndexType;
    using MatsetFloat = typename MatsetView::FloatType;
    namespace utils = axom::bump::utilities;
    const axom::ArrayView<axom::IndexType> &selectedZonesView = slice.m_indicesView;
    SLIC_ASSERT(selectedZonesView.size() > 0);

    // Copy the material_map if it exists.
    const char *keys[] = {"topology", "material_map"};
    for(int i = 0; i < 2; i++)
    {
      if(n_matset.has_child(keys[i])) n_newMatset[keys[i]] = n_matset.fetch_existing(keys[i]);
    }

    const auto conduitAllocatorId =
      axom::sidre::ConduitMemory::axomAllocIdToConduit(getAllocatorID());

    // Allocate sizes/offsets.
    AXOM_ANNOTATE_BEGIN("alloc");
    conduit::Node &n_sizes = n_newMatset["sizes"];
    n_sizes.set_allocator(conduitAllocatorId);
    n_sizes.set(conduit::DataType(utils::cpp2conduit<MatsetIndex>::id, selectedZonesView.size()));
    auto sizesView = utils::make_array_view<MatsetIndex>(n_sizes);

    conduit::Node &n_offsets = n_newMatset["offsets"];
    n_offsets.set_allocator(conduitAllocatorId);
    n_offsets.set(conduit::DataType(utils::cpp2conduit<MatsetIndex>::id, selectedZonesView.size()));
    auto offsetsView = utils::make_array_view<MatsetIndex>(n_offsets);
    AXOM_ANNOTATE_END("alloc");

    // Figure out overall size of the matset zones we're keeping.
    AXOM_ANNOTATE_BEGIN("size");
    MatsetView deviceMatsetView(m_matsetView);
    const axom::ArrayView<axom::IndexType> deviceSelectedZonesView(selectedZonesView);
    axom::IndexType totalSize {};
    if constexpr(axom::execution_space<ExecSpace>::onDevice())
    {
      axom::ReduceSum<ExecSpace, MatsetIndex> size_reduce(0);
      axom::for_all<ExecSpace>(
        selectedZonesView.size(),
        AXOM_LAMBDA(axom::IndexType index) {
          const auto nmats = deviceMatsetView.numberOfMaterials(deviceSelectedZonesView[index]);
          sizesView[index] = nmats;
          size_reduce += nmats;
        });
      totalSize = size_reduce.get();
    }
    else
    {
      axom::for_all<ExecSpace>(
        selectedZonesView.size(),
        AXOM_LAMBDA(axom::IndexType index) {
          const auto nmats = deviceMatsetView.numberOfMaterials(deviceSelectedZonesView[index]);
          sizesView[index] = nmats;
        });
    }
    AXOM_ANNOTATE_END("size");

    {
      AXOM_ANNOTATE_SCOPE("scan");
      axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);
    }

    AXOM_ANNOTATE_BEGIN("alloc2");
    if constexpr(!axom::execution_space<ExecSpace>::onDevice())
    {
      if(sizesView.size() > 0)
      {
        const auto lastIndex = sizesView.size() - 1;
        totalSize = offsetsView[lastIndex] + sizesView[lastIndex];
      }
    }

    // Allocate data for the rest of the matset.
    if(sizesView.size() > 0)
    {
      SLIC_ERROR_IF(totalSize == 0, "ReduceSum returned 0 for totalSize.");
    }
    conduit::Node &n_indices = n_newMatset["indices"];
    n_indices.set_allocator(conduitAllocatorId);
    n_indices.set(conduit::DataType(utils::cpp2conduit<MatsetIndex>::id, totalSize));
    auto indicesView = utils::make_array_view<MatsetIndex>(n_indices);

    conduit::Node &n_material_ids = n_newMatset["material_ids"];
    n_material_ids.set_allocator(conduitAllocatorId);
    n_material_ids.set(conduit::DataType(utils::cpp2conduit<MatsetIndex>::id, totalSize));
    auto materialIdsView = utils::make_array_view<MatsetIndex>(n_material_ids);

    conduit::Node &n_volume_fractions = n_newMatset["volume_fractions"];
    n_volume_fractions.set_allocator(conduitAllocatorId);
    n_volume_fractions.set(conduit::DataType(utils::cpp2conduit<MatsetFloat>::id, totalSize));
    auto volumeFractionsView = utils::make_array_view<MatsetFloat>(n_volume_fractions);
    AXOM_ANNOTATE_END("alloc2");

    // Fill in the matset data with the zones we're keeping.
    AXOM_ANNOTATE_BEGIN("copy");
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        const auto size = static_cast<int>(sizesView[index]);
        const auto offset = offsetsView[index];

        auto zoneMat = deviceMatsetView.beginZone(deviceSelectedZonesView[index]);
        for(int i = 0; i < size; i++, zoneMat++)
        {
          const auto destIndex = offset + i;
          materialIdsView[destIndex] = zoneMat.material_id();
          volumeFractionsView[destIndex] = zoneMat.volume_fraction();
          indicesView[destIndex] = destIndex;
        }
      });
    AXOM_ANNOTATE_END("copy");
  }

private:
  MatsetView m_matsetView;
  int m_allocator_id;
};

}  // end namespace bump
}  // end namespace axom

#endif
