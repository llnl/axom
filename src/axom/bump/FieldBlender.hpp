// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_FIELD_BLENDER_HPP_
#define AXOM_BUMP_FIELD_BLENDER_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/utilities/utilities.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/BlendData.hpp"
#include "axom/bump/IndexingPolicies.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
/*!
 * \brief This policy can be used with FieldBlender to select all blend groups.
 */
struct SelectAllPolicy
{
  AXOM_HOST_DEVICE
  static inline IndexType size(const BlendData &blend)
  {
    return blend.m_blendGroupSizesView.size();
  }

  AXOM_HOST_DEVICE
  static inline IndexType selectedIndex(const BlendData & /*blend*/, IndexType index)
  {
    return index;
  }
};

/*!
 * \brief This policy can be used with FieldBlender to select a subset of blend groups, according to m_selectedIndicesView.
 */
struct SelectSubsetPolicy
{
  AXOM_HOST_DEVICE
  static inline IndexType size(const BlendData &blend)
  {
    return blend.m_selectedIndicesView.size();
  }

  AXOM_HOST_DEVICE
  static inline IndexType selectedIndex(const BlendData &blend, IndexType index)
  {
    return blend.m_selectedIndicesView[index];
  }
};

/*!
 * \accelerated
 * \class FieldBlender
 *
 * \brief This class uses BlendData to generate a new blended field from an input field.
 *
 * \tparam ExecSpace The execution space where the work will occur.
 * \tparam SelectionPolicy The selection policy to use.
 * \tparam IndexingPolicy A class that provides operator[] that can transform node indices.
 */
template <typename ExecSpace, typename SelectionPolicy, typename IndexingPolicy = DirectIndexing>
class FieldBlender
{
public:
  /// Constructor
  FieldBlender() : m_indexing(), m_allocator_id(axom::execution_space<ExecSpace>::allocatorID()) { }

  /*!
   * \brief Constructor
   * \param indexing An object used to transform node indices.
   */
  FieldBlender(const IndexingPolicy &indexing)
    : m_indexing(indexing)
    , m_allocator_id(axom::execution_space<ExecSpace>::allocatorID())
  { }

  /*!
   * \brief Set the allocator id to use when allocating memory.
   *
   * \param allocator_id The allocator id to use when allocating memory.
   */
  void setAllocatorID(int allocator_id)
  {
    SLIC_ERROR_IF(!axom::isValidAllocatorID(allocator_id), "Invalid allocator id.");
    SLIC_ERROR_IF(!axom::execution_space<ExecSpace>::usesAllocId(allocator_id),
                  "Allocator id is not compatible with execution space.");
    m_allocator_id = allocator_id;
  }

  /*!
   * \brief Get the allocator id to use when allocating memory.
   *
   * \return The allocator id to use when allocating memory.
   */
  int getAllocatorID() const { return m_allocator_id; }

  /*!
   * \brief Create a new blended field from the \a n_input field and place it in \a n_output.
   *
   * \param blend The BlendData that will be used to make the new field.
   * \param n_input The input field that we're blending.
   * \param n_output The output node that will contain the new field.
   */
  void execute(const BlendData &blend, const conduit::Node &n_input, conduit::Node &n_output) const
  {
    n_output.reset();
    n_output["association"] = n_input["association"];
    n_output["topology"] = n_input["topology"];

    const conduit::Node &n_input_values = n_input["values"];
    conduit::Node &n_output_values = n_output["values"];
    const conduit::index_t nc = n_input_values.number_of_children();
    if(nc > 0)
    {
      for(conduit::index_t i = 0; i < nc; i++)
      {
        const conduit::Node &n_comp = n_input_values[i];
        conduit::Node &n_out_comp = n_output_values[n_comp.name()];
        blendSingleComponent(blend, n_comp, n_out_comp);
      }
    }
    else
    {
      blendSingleComponent(blend, n_input_values, n_output_values);
    }
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  /*!
   * \brief Blend data for a single field component.
   *
   * \param blend The BlendData that will be used to make the new field.
   * \param n_values The input values that we're blending.
   * \param n_output_values The output node that will contain the new field.
   */
  void blendSingleComponent(const BlendData &blend,
                            const conduit::Node &n_values,
                            conduit::Node &n_output_values) const
  {
    // We're allowing selectedIndicesView to be used to select specific blend
    // groups. If the user did not provide that, use all blend groups.
    const auto orig_size = blend.m_originalIdsView.size();
    const auto blend_size = SelectionPolicy::size(blend);
    const auto output_size = orig_size + blend_size;

    const auto conduit_allocator_id =
      axom::sidre::ConduitMemory::axomAllocIdToConduit(getAllocatorID());
    n_output_values.set_allocator(conduit_allocator_id);
    n_output_values.set(conduit::DataType(n_values.dtype().id(), output_size));

    views::nodeToArrayViewSame(n_values, n_output_values, [&](auto comp_view, auto out_view) {
      blendSingleComponentImpl(blend, comp_view, out_view);
    });
  }

  /*!
   * \brief Slice the source view and copy values into the output view.
   *
   * \param comp_view The source values view.
   * \param out_view The output values view.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename SrcView, typename OutputView>
  void blendSingleComponentImpl(const BlendData &blend, SrcView comp_view, OutputView out_view) const
  {
    using value_type = typename decltype(comp_view)::value_type;
    using accum_type = typename utilities::accumulation_traits<value_type>::type;

    // We're allowing selectedIndicesView to be used to select specific blend
    // groups. If the user did not provide that, use all blend groups.
    const auto orig_size = blend.m_originalIdsView.size();
    const auto blend_size = SelectionPolicy::size(blend);
    //    const auto output_size = orig_size + blend_size;

    const IndexingPolicy device_indexing(m_indexing);
    const BlendData device_blend(blend);

    // Copy over some original values to the start of the array.
    axom::for_all<ExecSpace>(
      orig_size,
      AXOM_LAMBDA(axom::IndexType index) {
        const auto src_index = device_blend.m_originalIdsView[index];
        out_view[index] = comp_view[src_index];
      });

    // Append blended values to the end of the array.
    axom::for_all<ExecSpace>(
      blend_size,
      AXOM_LAMBDA(axom::IndexType bgid) {
        // Get the blend group index we want.
        const auto selected_index = SelectionPolicy::selectedIndex(device_blend, bgid);
        const auto start = device_blend.m_blendGroupStartView[selected_index];
        const auto n_values = device_blend.m_blendGroupSizesView[selected_index];
        const auto dest_index = orig_size + bgid;
        if(n_values == 1)
        {
          const auto index = device_blend.m_blendIdsView[start];
          const auto src_index = device_indexing[index];
          out_view[dest_index] = comp_view[src_index];
        }
        else
        {
          const auto end = start + n_values;
          accum_type blended = 0;
          for(IndexType i = start; i < end; i++)
          {
            const auto index = device_blend.m_blendIdsView[i];
            const auto weight = device_blend.m_blendCoeffView[i];
            const auto src_index = device_indexing[index];
            blended += static_cast<accum_type>(comp_view[src_index]) * weight;
          }
          out_view[dest_index] = static_cast<value_type>(blended);
        }
      });
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  IndexingPolicy m_indexing {};
  int m_allocator_id;
};

}  // end namespace bump
}  // end namespace axom

#endif
