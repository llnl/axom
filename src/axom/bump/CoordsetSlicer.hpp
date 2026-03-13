// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_COORDSET_SLICER_HPP_
#define AXOM_BUMP_COORDSET_SLICER_HPP_

#include "axom/core.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/FieldSlicer.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{

/*!
 * \accelerated
 * \class CoordsetSlicer
 *
 * \brief This class uses SliceData to generate a new sliced coordset (pulling out specific points from input coordset).
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam CoordsetView The coordset view type to be operated on.
 *
 */
template <typename ExecSpace, typename CoordsetView>
class CoordsetSlicer
{
public:
  /// Constructor
  CoordsetSlicer(const CoordsetView &coordsetView)
    : m_coordsetView(coordsetView)
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
   * \brief Execute the slice on the \a n_input coordset and store the new sliced coordset in \a n_output.
   *
   * \param slice    The slice data that indicates how the coordset will be sliced.
   * \param n_input  A Conduit node containing the coordset to be sliced.
   * \param n_output A node that will contain the sliced coordset.
   *
   * \note We assume for now that n_input != n_output.
   */
  void execute(const SliceData &slice, const conduit::Node &n_input, conduit::Node &n_output)
  {
    AXOM_ANNOTATE_SCOPE("CoordsetSlicer");
    using value_type = typename CoordsetView::value_type;
    using PointType = typename CoordsetView::PointType;
    namespace utils = axom::bump::utilities;

    // Get the axis names for the output coordset.
    std::vector<std::string> axes(utils::coordsetAxes(n_input));

    const auto nComponents = axes.size();
    SLIC_ASSERT(PointType::DIMENSION == nComponents);
    SLIC_ASSERT(slice.m_indicesView.size() > 0);

    const auto conduitAllocatorId =
      axom::sidre::ConduitMemory::axomAllocIdToConduit(getAllocatorID());

    n_output.reset();
    n_output["type"] = "explicit";
    conduit::Node &n_values = n_output["values"];

    // Determine output size.
    const auto outputSize = slice.m_indicesView.size();

    // Make output nodes using axis names from the input coordset. Make array views too.
    axom::StackArray<axom::ArrayView<value_type>, PointType::DIMENSION> compViews;
    for(size_t i = 0; i < nComponents; i++)
    {
      // Allocate data in the Conduit node and make a view.
      conduit::Node &comp = n_values[axes[i]];
      comp.set_allocator(conduitAllocatorId);
      comp.set(conduit::DataType(utils::cpp2conduit<value_type>::id, outputSize));
      compViews[i] = utils::make_array_view<value_type>(comp);
    }

    // Select the nodes we want in the output.
    const CoordsetView deviceView(m_coordsetView);
    const auto deviceIndicesView = slice.m_indicesView;
    axom::for_all<ExecSpace>(
      outputSize,
      AXOM_LAMBDA(axom::IndexType index) {
        const auto srcIndex = deviceIndicesView[index];
        const auto pt = deviceView[srcIndex];

        // Store the point into the Conduit component arrays.
        for(int comp = 0; comp < PointType::DIMENSION; comp++)
        {
          compViews[comp][index] = pt[comp];
        }
      });
  }

private:
  CoordsetView m_coordsetView;
  int m_allocator_id;
};

}  // end namespace bump
}  // end namespace axom

#endif
