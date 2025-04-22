// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_SELECTED_ZONES_HPP_
#define AXOM_MIR_SELECTED_ZONES_HPP_

#include "axom/core.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{

/*!
 * \brief This class creates a view containing sorted selected zones, given
 *        an optional list of selected zones.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 */
template <typename ExecSpace>
class SelectedZones
{
public:
  /*!
   * \brief Constructor
   *
   * \param nzones The total number of zones in the associated topology.
   * \param n_options The node that contains the clipping options.
   * \param selectionKey The name of the node with the selection data in the options.
   *
   * The n_options node contains options that influence how the class runs.
   * The options can contain a "selectedZones" node that contains an array of
   * zone ids that will be processed. The array should exist in the memory space
   * that is appropriate for the execution space. If this node is not present
   * then all zones will be selected.
   *
   * \code{.yaml}
   *  selectedZones: [0,1,2,3...]
   * \endcode
   */
  SelectedZones(axom::IndexType nzones,
                const conduit::Node &n_options,
                const std::string &selectionKey = std::string("selectedZones"))
    : m_selectedZones()
    , m_selectedZonesView()
    , m_sorted(true)
  {
    buildSelectedZones(nzones, n_options, selectionKey);
  }

  /*!
   * \brief Set whether we need to sort the selected zone ids.
   *
   * \param sorted Whether the ids need to be sorted.
   *
   */
  void setSorted(bool sorted) { m_sorted = sorted; }

  /*!
   * \brief Return a view that contains the list of selected zone ids for the mesh.
   * \return A view that contains the list of selected zone ids for the mesh.
   */
  const axom::ArrayView<axom::IndexType> &view() const { return m_selectedZonesView; }

// The following members are protected (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif

  /*!
   * \brief The options may contain a "selectedZones" (or other provided name) member
   *        that is a list of zones on which to operate. If such an array is present,
   *        copy and sort it. If the zone list is not present, make an array that
   *        selects every zone.
   *
   * \param nzones The total number of zones that are possible.
   * \param n_options A Conduit node that contains the selection.
   * \param selectionKey The name of the node with the selection data in the options.
   *
   * \note selectedZones should contain local zone numbers, which in the case of
   *       strided-structured indexing are the [0..n) zone numbers that exist only
   *       within the selected window.
   */
  void buildSelectedZones(axom::IndexType nzones,
                          const conduit::Node &n_options,
                          const std::string &selectionKey)
  {
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    if(n_options.has_path(selectionKey))
    {
      // Store the zone list in m_selectedZones.
      int badValueCount = 0;
      views::IndexNode_to_ArrayView(n_options[selectionKey], [&](auto zonesView) {
        // It probably does not make sense to request more zones than we have in the mesh.
        SLIC_ASSERT(zonesView.size() <= nzones);

        badValueCount = buildSelectedZones(zonesView, nzones);
      });

      if(badValueCount > 0)
      {
        SLIC_ERROR(axom::fmt::format("Out of range {} values.", selectionKey));
      }
    }
    else
    {
      // Select all zones.
      m_selectedZones = axom::Array<axom::IndexType>(nzones, nzones, allocatorID);
      auto szView = m_selectedZonesView = m_selectedZones.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType zoneIndex) { szView[zoneIndex] = zoneIndex; });
    }
  }

  /*!
   * \brief Help build the selected zones, converting them to axom::IndexType and sorting them.
   *
   * \param zonesView The view that contains the source zone ids.
   * \param nzones The number of zones in the mesh.
   *
   * \return The number of invalid zone ids.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename ZonesViewType>
  int buildSelectedZones(ZonesViewType zonesView, axom::IndexType nzones)
  {
    using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    m_selectedZones = axom::Array<axom::IndexType>(zonesView.size(), zonesView.size(), allocatorID);
    auto szView = m_selectedZonesView = m_selectedZones.view();
    axom::for_all<ExecSpace>(
      szView.size(),
      AXOM_LAMBDA(axom::IndexType index) { szView[index] = zonesView[index]; });

    // Check that the selected zone values are in range.
    RAJA::ReduceSum<reduce_policy, int> errReduce(0);
    axom::for_all<ExecSpace>(
      szView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        const int err = (szView[index] < 0 || szView[index] >= nzones) ? 1 : 0;
        errReduce += err;
      });

    if(m_sorted)
    {
      // Make sure the selectedZones are sorted.
      RAJA::sort<loop_policy>(RAJA::make_span(szView.data(), szView.size()));
    }

    return errReduce.get();
  }

// The following members are protected (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif

  axom::Array<axom::IndexType> m_selectedZones;  // Storage for a list of selected zone ids.
  axom::ArrayView<axom::IndexType> m_selectedZonesView;
  bool m_sorted;
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
