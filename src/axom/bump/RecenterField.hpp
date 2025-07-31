// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_RECENTER_FIELD_HPP_
#define AXOM_BUMP_RECENTER_FIELD_HPP_

#include "axom/core.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/utilities/utilities.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

namespace axom
{
namespace bump
{

/*!
 * \brief Convert a field with one association type to a field of another association type using an o2mrelation.
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 */
template <typename ExecSpace>
class RecenterField
{
public:
  /*!
   * \brief Convert the input field to a different association type using the o2mrelation and store the new field in the output field.
   *
   * \param field       The input field.
   * \param relation    The node that contains an o2mrelation with nodes to zones.
   * \param outField[out] The node that will contain the new field.
   */
  void execute(const conduit::Node &field, const conduit::Node &relation, conduit::Node &outField) const
  {
    const std::string association = field.fetch_existing("association").as_string();

    // Assume that we're flipping the association.
    outField["association"] = (association == "element") ? "vertex" : "element";
    outField["topology"] = field["topology"];

    // Make output values.
    const conduit::Node &n_values = field["values"];
    if(n_values.number_of_children() > 0)
    {
      for(conduit::index_t c = 0; c < n_values.number_of_children(); c++)
      {
        const conduit::Node &n_comp = n_values[c];
        recenterSingleComponent(n_comp, relation, outField["values"][n_comp.name()]);
      }
    }
    else
    {
      recenterSingleComponent(n_values, relation, outField["values"]);
    }
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  /*!
   * \brief Recenter a single field component.
   *
   * \param relation    The node that contains an o2mrelation with nodes to zones.
   * \param n_comp      The input component.
   * \param n_out[out] The node that will contain the new field.
   */
  void recenterSingleComponent(const conduit::Node &n_comp,
                               const conduit::Node &relation,
                               conduit::Node &n_out) const
  {
    namespace utils = axom::bump::utilities;
    // Get the data field for the o2m relation.
    const auto data_paths = conduit::blueprint::o2mrelation::data_paths(relation);

    // Use the o2mrelation to average data from n_comp to the n_out.
    const conduit::Node &n_relvalues = relation[data_paths[0]];
    const conduit::Node &n_sizes = relation["sizes"];
    const conduit::Node &n_offsets = relation["offsets"];
    views::IndexNode_to_ArrayView_same(
      n_relvalues,
      n_sizes,
      n_offsets,
      [&](auto relView, auto sizesView, auto offsetsView) {
        // Allocate Conduit data through Axom.
        const auto relSize = sizesView.size();
        utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
        n_out.set_allocator(c2a.getConduitAllocatorID());
        n_out.set(conduit::DataType(n_comp.dtype().id(), relSize));

        views::Node_to_ArrayView_same(n_out, n_comp, [&](auto outView, auto compView) {
          recenterSingleComponentImpl(relView, sizesView, offsetsView, outView, compView);
        });
      });
  }

  /*!
   * \brief Recenter a single field component.
   *
   * \param relView The view that contains the ids for the relation.
   * \param sizesView The view that contains the sizes for the relation.
   * \param offsetsView The view that contains the offsets for the relation.
   * \param outView The view that contains the out data.
   * \param compView The view that contains the source data.
   */
  template <typename IndexView, typename DataView>
  void recenterSingleComponentImpl(IndexView relView,
                                   IndexView sizesView,
                                   IndexView offsetsView,
                                   DataView outView,
                                   DataView compView) const
  {
    using Precision = typename DataView::value_type;
    using AccumType = typename axom::bump::utilities::accumulation_traits<Precision>::type;
    const auto relSize = sizesView.size();
    axom::for_all<ExecSpace>(
      relSize,
      AXOM_LAMBDA(axom::IndexType relIndex) {
        const auto n = static_cast<axom::IndexType>(sizesView[relIndex]);
        const auto offset = offsetsView[relIndex];

        AccumType sum {};
        for(axom::IndexType i = 0; i < n; i++)
        {
          const auto id = relView[offset + i];
          sum += static_cast<AccumType>(compView[id]);
        }

        outView[relIndex] = static_cast<Precision>(sum / n);
      });
  }
};

}  // end namespace bump
}  // end namespace axom

#endif
