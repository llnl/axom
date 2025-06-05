// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_FIELD_BLENDER_HPP_
#define AXOM_BUMP_FIELD_BLENDER_HPP_

#include "axom/core.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/utilities/utilities.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/BlendData.hpp"
#include "axom/bump/IndexingPolicies.hpp"

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
  FieldBlender() : m_indexing() { }

  /*!
   * \brief Constructor
   * \param indexing An object used to transform node indices.
   */
  FieldBlender(const IndexingPolicy &indexing) : m_indexing(indexing) { }

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
    if(n_input_values.number_of_children() > 0)
    {
      for(conduit::index_t i = 0; i < n_input_values.number_of_children(); i++)
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
    const auto origSize = blend.m_originalIdsView.size();
    const auto blendSize = SelectionPolicy::size(blend);
    const auto outputSize = origSize + blendSize;

    // Allocate Conduit data through Axom.
    utilities::ConduitAllocateThroughAxom<ExecSpace> c2a;
    n_output_values.set_allocator(c2a.getConduitAllocatorID());
    n_output_values.set(conduit::DataType(n_values.dtype().id(), outputSize));

    views::Node_to_ArrayView_same(n_values, n_output_values, [&](auto compView, auto outView) {
      blendSingleComponentImpl(blend, compView, outView);
    });
  }

  /*!
   * \brief Slice the source view and copy values into the output view.
   *
   * \param valuesView The source values view.
   * \param outputView The output values view.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename SrcView, typename OutputView>
  void blendSingleComponentImpl(const BlendData &blend, SrcView compView, OutputView outView) const
  {
    using value_type = typename decltype(compView)::value_type;
    using accum_type = typename utilities::accumulation_traits<value_type>::type;

    // We're allowing selectedIndicesView to be used to select specific blend
    // groups. If the user did not provide that, use all blend groups.
    const auto origSize = blend.m_originalIdsView.size();
    const auto blendSize = SelectionPolicy::size(blend);
    //    const auto outputSize = origSize + blendSize;

    const IndexingPolicy deviceIndexing(m_indexing);
    const BlendData deviceBlend(blend);

    // Copy over some original values to the start of the array.
    axom::for_all<ExecSpace>(
      origSize,
      AXOM_LAMBDA(axom::IndexType index) {
        const auto srcIndex = deviceBlend.m_originalIdsView[index];
        outView[index] = compView[srcIndex];
      });

    // Append blended values to the end of the array.
    axom::for_all<ExecSpace>(
      blendSize,
      AXOM_LAMBDA(axom::IndexType bgid) {
        // Get the blend group index we want.
        const auto selectedIndex = SelectionPolicy::selectedIndex(deviceBlend, bgid);
        const auto start = deviceBlend.m_blendGroupStartView[selectedIndex];
        const auto nValues = deviceBlend.m_blendGroupSizesView[selectedIndex];
        const auto destIndex = origSize + bgid;
        if(nValues == 1)
        {
          const auto index = deviceBlend.m_blendIdsView[start];
          const auto srcIndex = deviceIndexing[index];
          outView[destIndex] = compView[srcIndex];
        }
        else
        {
          const auto end = start + nValues;
          accum_type blended = 0;
          for(IndexType i = start; i < end; i++)
          {
            const auto index = deviceBlend.m_blendIdsView[i];
            const auto weight = deviceBlend.m_blendCoeffView[i];
            const auto srcIndex = deviceIndexing[index];
            blended += static_cast<accum_type>(compView[srcIndex]) * weight;
          }
          outView[destIndex] = static_cast<value_type>(blended);
        }
      });
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  IndexingPolicy m_indexing {};
};

}  // end namespace bump
}  // end namespace axom

#endif
