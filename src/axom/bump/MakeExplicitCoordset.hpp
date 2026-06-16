// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_MAKE_EXPLICIT_COORDSET_HPP_
#define AXOM_BUMP_MAKE_EXPLICIT_COORDSET_HPP_

#include "axom/core.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/utilities/utilities.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/dispatch_coordset.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
/*!
 * \brief Convert a coordset to an explicit coordset.
 *
 * \tparam ExecSpace The execution space where the conversion will execute.
 */
template <typename ExecSpace>
class MakeExplicitCoordset
{
public:
  /*!
   * \brief Converts the supplied coordset into an explicit coordset if it is not already one.
   *
   * \param[inout] n_coordset The coordset to convert.
   * \param allocator_id The allocator id to use when allocating new coordinate memory.
   */
  static void execute(conduit::Node &n_coordset,
                      int allocator_id = axom::execution_space<ExecSpace>::allocatorID())
  {
    const std::string cstype = n_coordset["type"].as_string();
    if(cstype != "explicit")
    {
      conduit::Node n_dest_coordset;
      if(cstype == "uniform")
      {
        axom::bump::views::dispatch_uniform_coordset(n_coordset, [&](auto coordsetView)
        {
          convert(coordsetView, n_dest_coordset, allocator_id);
        });
      }
      else if(cstype == "rectilinear")
      {
        axom::bump::views::dispatch_rectilinear_coordset(n_coordset, [&](auto coordsetView)
        {
          convert(coordsetView, n_dest_coordset, allocator_id);
        });
      }
      else
      {
        SLIC_ERROR(axom::fmt::format("Unsupported coordset type {}.", cstype));
      }

      n_coordset["type"] = "explicit";
      n_coordset["values"].swap(n_dest_coordset["values"]);
    }
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  /*!
   * \brief Copy the coordinates from the input coordset view into values
   *        components in the supplied output coordset node.
   *
   * \param coordsetView A view for retrieving points from the source coordset.
   * \param n_dest_coordset A Conduit node in which to construct the new coordset.
   * \param allocator_id The allocator id to use when allocating new coordinate memory.
   */
  template <typename CoordsetView>
  static void convert(CoordsetView coordsetView, conduit::Node &n_dest_coordset, int allocator_id)
  {
    const auto conduitAllocatorId = axom::sidre::ConduitMemory::axomAllocIdToConduit(allocator_id);

    // Make new coordinate arrays
    const char *names[] = {"x", "y", "z"};
    using value_type = typename CoordsetView::value_type;
    namespace utils = axom::bump::utilities;
    axom::ArrayView<value_type> comps[3];
    conduit::Node &n_values = n_dest_coordset["values"];
    for(int c = 0; c < coordsetView.dimension(); c++)
    {
      conduit::Node &n_comp = n_values[names[c]];
      n_comp.set_allocator(conduitAllocatorId);
      n_comp.set(conduit::DataType(utils::cpp2conduit<value_type>::id, coordsetView.size()));
      comps[c] = utils::make_array_view<value_type>(n_comp);
    }

    // Copy data from the view into the new coordinate array views.
    axom::for_all<ExecSpace>(coordsetView.size(), AXOM_LAMBDA(axom::IndexType i)
    {
      const auto pt = coordsetView[i];
      for(int c = 0; c < coordsetView.dimension(); c++)
      {
        comps[c][i] = pt[c];
      }
    });
  }
};

}  // end namespace bump
}  // end namespace axom
#endif
