// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_MINMAX_HPP_
#define AXOM_BUMP_MINMAX_HPP_

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/reductions.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/export/bump.h"

#include <conduit/conduit.hpp>

#include <utility>

namespace axom
{
namespace bump
{

//------------------------------------------------------------------------------
/*!
 * \brief Get the min/max values for the data in a Conduit node or ArrayView.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam ReturnType The data type of the returned min/max values.
 */
template <typename ExecSpace, typename ReturnType>
struct MinMax
{
  /*!
   * \brief Get the min/max values for the data in a Conduit node.
   *
   * \param[in] n The Conduit node whose data we're checking.
   *
   * \return A pair containing the min,max values in the node.
   */
  static std::pair<ReturnType, ReturnType> execute(const conduit::Node &n)
  {
    SLIC_ASSERT(n.dtype().number_of_elements() > 0);
    std::pair<ReturnType, ReturnType> retval;

    axom::bump::views::Node_to_ArrayView(n, [&](auto nview) { retval = execute(nview); });
    return retval;
  }

  /*!
   * \brief Get the min/max values for the data in an ArrayView.
   *
   * \param[in] n The Conduit node whose data we're checking.
   *
   * \return A pair containing the min,max values in the node.
   */
  template <typename T>
  static std::pair<ReturnType, ReturnType> execute(const axom::ArrayView<T> nview)
  {
    axom::ReduceMin<ExecSpace, T> vmin(axom::numeric_limits<T>::max());
    axom::ReduceMax<ExecSpace, T> vmax(axom::numeric_limits<T>::min());

    axom::for_all<ExecSpace>(
      nview.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        vmin.min(nview[index]);
        vmax.max(nview[index]);
      });

    return std::pair<ReturnType, ReturnType> {static_cast<ReturnType>(vmin.get()),
                                              static_cast<ReturnType>(vmax.get())};
  }
};

}  // end namespace bump
}  // end namespace axom

#endif
