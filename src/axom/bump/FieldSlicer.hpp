// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_FIELD_SLICER_HPP_
#define AXOM_BUMP_FIELD_SLICER_HPP_

#include "axom/core.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/IndexingPolicies.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
/*!
 * \brief Contains the indices to be sliced out of a Blueprint field.
 */
struct SliceData
{
  axom::ArrayView<IndexType> m_indicesView;
};

/*!
 * \accelerated
 * \class FieldSlicer
 *
 * \brief This class uses SliceData to generate a new sliced field from an input field.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam IndexingPolicy A class that provides operator[] that can transform zone indices.
 *
 */
template <typename ExecSpace, typename IndexingPolicy = DirectIndexing>
class FieldSlicer
{
public:
  /// Constructor
  FieldSlicer() : m_indexing() { }

  /*!
   * \brief Constructor
   * \param indexing An object used to transform node indices.
   */
  FieldSlicer(const IndexingPolicy &indexing) : m_indexing(indexing) { }

  /*!
   * \brief Execute the slice on the \a n_input field and store the new sliced field in \a n_output.
   *
   * \param slice    The slice data that indicates how the field will be sliced.
   * \param n_input  A Conduit node containing the field to be sliced.
   * \param n_output A node that will contain the sliced field.
   *
   * \note We assume for now that n_input != n_output.
   */
  void execute(const SliceData &slice, const conduit::Node &n_input, conduit::Node &n_output)
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
        sliceSingleComponent(slice, n_comp, n_out_comp);
      }
    }
    else
    {
      sliceSingleComponent(slice, n_input_values, n_output_values);
    }
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  /*!
   * \brief Slice data for a single field component.
   *
   * \param slice The SliceData that will be used to make the new field.
   * \param n_values The input values that we're slicing.
   * \param n_output_values The output node that will contain the new field.
   */
  void sliceSingleComponent(const SliceData &slice,
                            const conduit::Node &n_values,
                            conduit::Node &n_output_values) const
  {
    namespace utils = axom::bump::utilities;
    const auto outputSize = slice.m_indicesView.size();

    // Allocate Conduit data through Axom.
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    n_output_values.set_allocator(c2a.getConduitAllocatorID());
    n_output_values.set(conduit::DataType(n_values.dtype().id(), outputSize));

    views::Node_to_ArrayView_same(n_values, n_output_values, [&](auto valuesView, auto outputView) {
      sliceSingleComponentImpl(slice, valuesView, outputView);
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
  template <typename ValuesView, typename OutputView>
  void sliceSingleComponentImpl(const SliceData &slice, ValuesView valuesView, OutputView outputView) const
  {
    IndexingPolicy deviceIndexing(m_indexing);
    SliceData deviceSlice(slice);
    axom::for_all<ExecSpace>(
      outputView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        const auto zoneIndex = deviceSlice.m_indicesView[index];
        const auto transformedIndex = deviceIndexing[zoneIndex];
        outputView[index] = valuesView[transformedIndex];
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
