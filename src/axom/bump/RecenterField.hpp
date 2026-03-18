// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_RECENTER_FIELD_HPP_
#define AXOM_BUMP_RECENTER_FIELD_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/utilities/utilities.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"

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
  RecenterField() : m_allocator_id(axom::execution_space<ExecSpace>::allocatorID()) { }

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
   * \brief Convert the input field to a different association type using the o2mrelation and store the new field in the output field.
   *
   * \param field       The input field.
   * \param relation    The node that contains an o2mrelation with nodes to zones.
   * \param out_field[out] The node that will contain the new field.
   */
  void execute(const conduit::Node &field, const conduit::Node &relation, conduit::Node &out_field) const
  {
    const std::string association = field.fetch_existing("association").as_string();

    // Assume that we're flipping the association.
    out_field["association"] = (association == "element") ? "vertex" : "element";
    out_field["topology"] = field["topology"];

    // Make output values.
    const conduit::Node &n_values = field["values"];
    if(n_values.number_of_children() > 0)
    {
      for(conduit::index_t c = 0; c < n_values.number_of_children(); c++)
      {
        const conduit::Node &n_comp = n_values[c];
        recenterSingleComponent(n_comp, relation, out_field["values"][n_comp.name()]);
      }
    }
    else
    {
      recenterSingleComponent(n_values, relation, out_field["values"]);
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
    views::indexNodeToArrayViewSame(
      n_relvalues,
      n_sizes,
      n_offsets,
      [&](auto rel_view, auto sizes_view, auto offsets_view) {
        const auto conduit_allocator_id =
          axom::sidre::ConduitMemory::axomAllocIdToConduit(getAllocatorID());
        const auto rel_size = sizes_view.size();
        n_out.set_allocator(conduit_allocator_id);
        n_out.set(conduit::DataType(n_comp.dtype().id(), rel_size));

        views::nodeToArrayViewSame(n_out, n_comp, [&](auto out_view, auto comp_view) {
          recenterSingleComponentImpl(rel_view, sizes_view, offsets_view, out_view, comp_view);
        });
      });
  }

  /*!
   * \brief Recenter a single field component.
   *
   * \param rel_view The view that contains the ids for the relation.
   * \param sizes_view The view that contains the sizes for the relation.
   * \param offsets_view The view that contains the offsets for the relation.
   * \param out_view The view that contains the out data.
   * \param comp_view The view that contains the source data.
   */
  template <typename IndexView, typename DataView>
  void recenterSingleComponentImpl(IndexView rel_view,
                                   IndexView sizes_view,
                                   IndexView offsets_view,
                                   DataView out_view,
                                   DataView comp_view) const
  {
    using Precision = typename DataView::value_type;
    using AccumType = typename axom::bump::utilities::accumulation_traits<Precision>::type;
    const auto rel_size = sizes_view.size();
    axom::for_all<ExecSpace>(
      rel_size,
      AXOM_LAMBDA(axom::IndexType rel_index) {
        const auto n = static_cast<axom::IndexType>(sizes_view[rel_index]);
        const auto offset = offsets_view[rel_index];

        AccumType sum {};
        for(axom::IndexType i = 0; i < n; i++)
        {
          const auto id = rel_view[offset + i];
          sum += static_cast<AccumType>(comp_view[id]);
        }

        out_view[rel_index] = static_cast<Precision>(sum / n);
      });
  }

private:
  int m_allocator_id;
};

}  // end namespace bump
}  // end namespace axom

#endif
