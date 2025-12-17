// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EQUIZ_ALGORITHM_DETAIL_HPP_
#define AXOM_MIR_EQUIZ_ALGORITHM_DETAIL_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

#include <algorithm>
#include <string>

namespace axom
{
namespace mir
{
using MaterialID = int;
using MaterialIDArray = axom::Array<MaterialID>;
using MaterialIDView = axom::ArrayView<MaterialID>;
using MaterialVF = float;
using MaterialVFArray = axom::Array<MaterialVF>;
using MaterialVFView = axom::ArrayView<MaterialVF>;

constexpr static int NULL_MATERIAL = -1;
constexpr static MaterialVF NULL_MATERIAL_VF = -1.f;

namespace detail
{
/*!
 * \brief This class is an intersection policy compatible with ClipField. It
 *        helps determine clip cases and weights using material-aware logic.
 *
 * \tparam TopologyView The type of topology view being used in the code that
 *                      uses this intersector.
 * \tparam CoordsetView The type of coordset view being used in the code that
 *                      uses this intersector.
 * \tparam MAXMATERIALS The max number of materials to handle.
 */
template <typename TopologyView, typename CoordsetView, int MAXMATERIALS = 10>
class MaterialIntersector
{
public:
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /*!
   * \brief This is a view class for MatsetIntersector that can be used in device code.
   */
  struct View
  {
    static constexpr int INVALID_INDEX = -1;

    /*!
     * \brief Determine the clipping case, taking into account the zone's material
     *        and the current material being added.
     *
     * \param zoneIndex The zone index in the zoneMaterialView field.
     * \param nodeIds A view containing node ids for the current zone.
     *
     * \return The clip case number for the zone.
     */
    AXOM_HOST_DEVICE
    axom::IndexType determineTableCase(axom::IndexType zoneIndex,
                                       const ConnectivityView &nodeIdsView) const
    {
      // Determine the matvf view index for the material that owns the zone.
      int backgroundIndex = INVALID_INDEX;
      int zoneMatID = m_zoneMatNumberView[zoneIndex];
      if(zoneMatID != NULL_MATERIAL) backgroundIndex = matNumberToIndex(zoneMatID);

      axom::IndexType clipcase = 0;
      const auto n = nodeIdsView.size();
      for(IndexType i = 0; i < n; i++)
      {
        const auto nid = nodeIdsView[i];
        SLIC_ASSERT_MSG(
          nid >= 0 && nid < static_cast<ConnectivityType>(m_matvfViews[0].size()),
          axom::fmt::format("Node id {} is not in range [0, {}).", nid, m_matvfViews[0].size()));

        // clang-format off
        MaterialVF vf1 = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][nid] : NULL_MATERIAL_VF;
        MaterialVF vf2 = (m_currentMaterialIndex != INVALID_INDEX) ? m_matvfViews[m_currentMaterialIndex][nid] : 0;
        // clang-format on

        clipcase |= (vf2 > vf1) ? (1 << i) : 0;
      }
      return clipcase;
    }

    /*!
     * \brief Compute the weight of a clip value along an edge (id0, id1) using the clip field and value.
     *
     * \param id0 The mesh node at the start of the edge.
     * \param id1 The mesh node at the end of the edge.
     */
    AXOM_HOST_DEVICE
    float computeWeight(axom::IndexType zoneIndex, ConnectivityType id0, ConnectivityType id1) const
    {
      // Determine the matvf view index for the material that owns the zone.
      int backgroundIndex = INVALID_INDEX;
      int zoneMatID = m_zoneMatNumberView[zoneIndex];
      if(zoneMatID != NULL_MATERIAL) backgroundIndex = matNumberToIndex(zoneMatID);
      // Determine the matvf view index for the current material.
      SLIC_ASSERT_MSG(
        id0 >= 0 && id0 < static_cast<ConnectivityType>(m_matvfViews[0].size()),
        axom::fmt::format("Node id {} is not in range [0, {}).", id0, m_matvfViews[0].size()));
      SLIC_ASSERT_MSG(
        id1 >= 0 && id1 < static_cast<ConnectivityType>(m_matvfViews[0].size()),
        axom::fmt::format("Node id {} is not in range [0, {}).", id1, m_matvfViews[0].size()));

      // Get the volume fractions for mat1, mat2 at the edge endpoints id0, id1.
      MaterialVF vf1[2], vf2[2];
      // clang-format off
      vf1[0] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id0] : NULL_MATERIAL_VF;
      vf1[1] = (backgroundIndex != INVALID_INDEX) ? m_matvfViews[backgroundIndex][id1] : NULL_MATERIAL_VF;
      vf2[0] = (m_currentMaterialIndex != INVALID_INDEX) ? m_matvfViews[m_currentMaterialIndex][id0] : 0;
      vf2[1] = (m_currentMaterialIndex != INVALID_INDEX) ? m_matvfViews[m_currentMaterialIndex][id1] : 0;
      // clang-format on

      float numerator = vf2[0] - vf1[0];
      float denominator = -vf1[0] + vf1[1] + vf2[0] - vf2[1];

      float t = 0.f;
      if(denominator != 0.f)
      {
        t = numerator / denominator;
      }
      t = axom::utilities::clampVal(t, 0.f, 1.f);

      return t;
    }

    /*!
     * \brief Return the volume fraction array index in m_matIndicesView for the
     *        given material number \a matNumber.
     *
     * \param matNumber A material number that occurs in the matset material ids.
     *
     * \return The m_matNumbersView index on success; INVALID_INDEX on failure.
     */
    AXOM_HOST_DEVICE
    inline int matNumberToIndex(int matNumber) const
    {
      auto index = axom::utilities::binary_search(m_matNumbersView, matNumber);
      return (index != -1) ? m_matIndicesView[index] : INVALID_INDEX;
    }

    /// Helper initialization methods for the host.

    void addMaterial(const MaterialVFView &matvf) { m_matvfViews.push_back(matvf); }

    void setMaterialNumbers(const axom::ArrayView<int> &matNumbersView)
    {
      m_matNumbersView = matNumbersView;
    }

    void setMaterialIndices(const axom::ArrayView<int> &matIndicesView)
    {
      m_matIndicesView = matIndicesView;
    }

    void setZoneMaterialID(const axom::ArrayView<int> &zoneMatsView)
    {
      m_zoneMatNumberView = zoneMatsView;
    }

    void setCurrentMaterial(int matNumber, int matNumberIndex)
    {
      m_currentMaterial = matNumber;
      m_currentMaterialIndex = matNumberIndex;
    }

    axom::StaticArray<MaterialVFView, MAXMATERIALS> m_matvfViews {};  //!< Array of volume fraction views
    axom::ArrayView<int> m_matNumbersView {};  //!< Sorted array of material numbers.
    axom::ArrayView<int> m_matIndicesView {};  //!< Array of indices into m_matvfViews for the material numbers.
    axom::ArrayView<int> m_zoneMatNumberView {};  //!< Contains the current material number that owns each zone.
    int m_currentMaterial {};                     //!< The current material.
    int m_currentMaterialIndex {};  //!< The current material's index in the m_matvfViews.
  };

  /*!
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const TopologyView &AXOM_UNUSED_PARAM(topologyView),
                  const CoordsetView &AXOM_UNUSED_PARAM(coordsetView),
                  const conduit::Node &AXOM_UNUSED_PARAM(n_options),
                  const conduit::Node &AXOM_UNUSED_PARAM(n_topology),
                  const conduit::Node &AXOM_UNUSED_PARAM(n_coordset),
                  const conduit::Node &AXOM_UNUSED_PARAM(n_fields))
  { }

  /*!
   * \brief Determine the name of the topology on which to operate.
   * \param n_input The input mesh node.
   * \param n_options The options.
   * \return The name of the toplogy on which to operate.
   */
  std::string getTopologyName(const conduit::Node &AXOM_UNUSED_PARAM(n_input),
                              const conduit::Node &n_options) const
  {
    return n_options["topology"].as_string();
  }

  /// Set various attributes.

  void addMaterial(const MaterialVFView &matvf) { m_view.addMaterial(matvf); }

  void setMaterialNumbers(const axom::ArrayView<int> &matNumbers)
  {
    m_view.setMaterialNumbers(matNumbers);
  }

  void setMaterialIndices(const axom::ArrayView<int> &matIndices)
  {
    m_view.setMaterialIndices(matIndices);
  }

  void setZoneMaterialID(const axom::ArrayView<int> &zoneMatsView)
  {
    m_view.setZoneMaterialID(zoneMatsView);
  }

  void setCurrentMaterial(int matNumber, int matNumberIndex)
  {
    m_view.setCurrentMaterial(matNumber, matNumberIndex);
  }

  /*!
   * \brief Return a new instance of the view.
   * \return A new instance of the view.
   * \note Call this after all values are set.
   */
  View view() const { return m_view; }

private:
  View m_view {};
};

}  // end namespace detail
}  // end namespace mir
}  // end namespace axom

#endif
