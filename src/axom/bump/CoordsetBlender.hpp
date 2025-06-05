// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_COORDSET_BLENDER_HPP_
#define AXOM_BUMP_COORDSET_BLENDER_HPP_

#include "axom/core.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/FieldBlender.hpp"
#include "axom/bump/BlendData.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
/*!
 * \accelerated
 * \class CoordsetBlender
 *
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam CSVType The coordset view type.
 * \tparam SelectionPolicy The selection policy to use.
 *
 * \brief This class uses BlendData to generate a new blended field from an input coordset.
 *        The output coordset will be explicit.
 */
template <typename ExecSpace, typename CSVType, typename SelectionPolicy>
class CoordsetBlender
{
public:
  using CoordsetViewType = CSVType;

  /*!
   * \brief Create a new blended field from the \a n_input field and place it in \a n_output.
   *
   * \param blend Blend data that describes how to create the new coordset.
   * \param blend The BlendData that will be used to make the new coordset.
   * \param n_input The input coordset that we're blending.
   * \param n_output The output node that will contain the new coordset.
   *
   * \note The coordset view must agree with the coordset in n_input. We pass both
   *       a view and the coordset node since the view may not be able to contain
   *       some coordset metadata and remain trivially copyable.
   */
  void execute(const BlendData &blend,
               const CoordsetViewType &view,
               const conduit::Node &n_input,
               conduit::Node &n_output) const
  {
    using value_type = typename CoordsetViewType::value_type;
    using PointType = typename CoordsetViewType::PointType;
    using VectorType = axom::primal::Vector<value_type, PointType::DIMENSION>;
    namespace utils = axom::bump::utilities;

    // Get the axis names for the output coordset. For uniform, prefer x,y,z
    // instead of i,j,k since we're making an explicit coordset.
    std::vector<std::string> axes(utils::coordsetAxes(n_input));

    const auto nComponents = axes.size();
    SLIC_ASSERT(PointType::DIMENSION == nComponents);

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    n_output.reset();
    n_output["type"] = "explicit";
    conduit::Node &n_values = n_output["values"];

    // Determine output size.
    const auto origSize = blend.m_originalIdsView.size();
    const auto blendSize = SelectionPolicy::size(blend);
    const auto outputSize = origSize + blendSize;

    // Make output nodes using axis names from the input coordset. Make array views too.
    axom::StackArray<axom::ArrayView<value_type>, PointType::DIMENSION> compViews;
    for(size_t i = 0; i < nComponents; i++)
    {
      // Allocate data in the Conduit node and make a view.
      conduit::Node &comp = n_values[axes[i]];
      comp.set_allocator(c2a.getConduitAllocatorID());
      comp.set(conduit::DataType(utils::cpp2conduit<value_type>::id, outputSize));
      compViews[i] = utils::make_array_view<value_type>(comp);
    }

    const CoordsetViewType deviceView(view);
    const BlendData deviceBlend(blend);

    // Copy over some original values to the start of the array.
    axom::for_all<ExecSpace>(
      origSize,
      AXOM_LAMBDA(axom::IndexType index) {
        const auto srcIndex = deviceBlend.m_originalIdsView[index];
        const auto pt = deviceView[srcIndex];

        // Store the point into the Conduit component arrays.
        for(int comp = 0; comp < PointType::DIMENSION; comp++)
        {
          compViews[comp][index] = pt[comp];
        }
      });

    // Append blended values to the end of the array.
    axom::for_all<ExecSpace>(
      blendSize,
      AXOM_LAMBDA(axom::IndexType bgid) {
        // Get the blend group index we want.
        const auto selectedIndex = SelectionPolicy::selectedIndex(deviceBlend, bgid);
        const auto start = deviceBlend.m_blendGroupStartView[selectedIndex];
        const auto nValues = deviceBlend.m_blendGroupSizesView[selectedIndex];
        const auto destIndex = bgid + origSize;

        VectorType blended {};
        if(nValues == 1)
        {
          const auto index = deviceBlend.m_blendIdsView[start];
          blended = VectorType(deviceView[index]);
        }
        else
        {
          const auto end = start + deviceBlend.m_blendGroupSizesView[selectedIndex];
          // Blend points for this blend group.
          for(IndexType i = start; i < end; i++)
          {
            const auto index = deviceBlend.m_blendIdsView[i];
            const auto weight = deviceBlend.m_blendCoeffView[i];
            blended += (VectorType(deviceView[index]) * static_cast<value_type>(weight));
          }
        }

        // Store the point into the Conduit component arrays.
        for(int comp = 0; comp < PointType::DIMENSION; comp++)
        {
          compViews[comp][destIndex] = blended[comp];
        }
      });
  }
};

}  // end namespace bump
}  // end namespace axom

#endif
