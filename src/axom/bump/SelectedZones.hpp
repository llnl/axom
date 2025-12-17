// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_SELECTED_ZONES_HPP_
#define AXOM_BUMP_SELECTED_ZONES_HPP_

#include "axom/core.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
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
   * \param n_options The node that contains the options.
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
    : m_selectionKey(selectionKey)
    , m_selectedZones()
    , m_selectedZonesView()
    , m_sorted(true)
  {
    buildSelectedZones(nzones, n_options);
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

  /*!
   * \brief Return the selection key for the options.
   *
   * \return The name of the key in the options that this class looks for.
   */
  const std::string &selectionKey() const { return m_selectionKey; }

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
   *
   * \note selectedZones should contain local zone numbers, which in the case of
   *       strided-structured indexing are the [0..n) zone numbers that exist only
   *       within the selected window.
   */
  void buildSelectedZones(axom::IndexType nzones, const conduit::Node &n_options)
  {
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    if(n_options.has_path(m_selectionKey))
    {
      // Store the zone list in m_selectedZones.
      int badValueCount = 0;
      views::IndexNode_to_ArrayView(n_options[m_selectionKey], [&](auto zonesView) {
        // It probably does not make sense to request more zones than we have in the mesh.
        SLIC_ASSERT(zonesView.size() <= nzones);

        badValueCount = buildSelectedZones(zonesView, nzones);
      });

      if(badValueCount > 0)
      {
        SLIC_ERROR(axom::fmt::format("Out of range {} values.", m_selectionKey));
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
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    m_selectedZones = axom::Array<axom::IndexType>(zonesView.size(), zonesView.size(), allocatorID);
    auto szView = m_selectedZonesView = m_selectedZones.view();
    axom::for_all<ExecSpace>(
      szView.size(),
      AXOM_LAMBDA(axom::IndexType index) { szView[index] = zonesView[index]; });

    // Check that the selected zone values are in range.
    axom::ReduceSum<ExecSpace, int> errReduce(0);
    axom::for_all<ExecSpace>(
      szView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        const int err = (szView[index] < 0 || szView[index] >= nzones) ? 1 : 0;
        errReduce += err;
      });

    if(m_sorted)
    {
      // Make sure the selectedZones are sorted.
      axom::sort<ExecSpace>(szView);
    }

    return errReduce.get();
  }

// The following members are protected (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif

  std::string m_selectionKey;
  axom::Array<axom::IndexType> m_selectedZones;  // Storage for a list of selected zone ids.
  axom::ArrayView<axom::IndexType> m_selectedZonesView;
  bool m_sorted;
};

}  // end namespace bump
}  // end namespace axom

#endif
