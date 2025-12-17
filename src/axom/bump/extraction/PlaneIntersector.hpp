// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internal.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_PLANE_INTERSECTOR_HPP_
#define AXOM_BUMP_PLANE_INTERSECTOR_HPP_

#include "axom/core.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"
#include "axom/bump/utilities/utilities.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
namespace extraction
{

/*!
 * \brief This class helps TableBasedExtractor determine intersection cases and
 *        weights using a plane designated by the options.
 *
 * \tparam TopologyView The type of topology view being used in the code that
 *                      uses this intersector.
 * \tparam CoordsetView The type of coordset view being used in the code that
 *                      uses this intersector.
 */
template <typename TopologyView, typename CoordsetView>
class PlaneIntersector
{
public:
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;
  static constexpr int NDIMS = CoordsetView::dimension();
  using value_type = typename CoordsetView::value_type;
  using PlaneType = axom::primal::Plane<value_type, NDIMS>;

  /*!
   * \brief This is a view class for FieldIntersector that can be used in device code.
   */
  struct View
  {
    /*!
     * \brief Determine the signed distance to the plane.
     *
     * \param nodeId The index of the node to query for distance.
     *
     * \return The signed distance to the plane.
     */
    AXOM_HOST_DEVICE
    value_type distance(axom::IndexType nodeId) const
    {
      return m_plane.signedDistance(m_coordsetView[nodeId]);
    }

    /*!
     * \brief Given a zone index and the node ids that comprise the zone, return
     *        the appropriate table case.
     *
     * \param zoneIndex The zone index.
     * \param nodeIds A view containing node ids for the zone.
     */
    AXOM_HOST_DEVICE
    axom::IndexType determineTableCase(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex),
                                       const ConnectivityView &nodeIds) const
    {
      axom::IndexType caseNumber = 0, numIds = nodeIds.size();
      for(IndexType i = 0; i < numIds; i++)
      {
        const auto dist = distance(nodeIds[i]);
        caseNumber |= (dist > value_type {0}) ? (1 << i) : 0;
      }
      return caseNumber;
    }

    /*!
     * \brief Compute the weight of a clip value along an edge (id0, id1) using the clip field and value.
     *
     * \param id0 The mesh node at the start of the edge.
     * \param id1 The mesh node at the end of the edge.
     *
     * \return A parametric position t [0,1] where we locate \a clipValues in [d0,d1].
     */
    AXOM_HOST_DEVICE
    value_type computeWeight(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex),
                             ConnectivityType id0,
                             ConnectivityType id1) const
    {
      const value_type d0 = distance(id0);
      const value_type d1 = distance(id1);
      constexpr value_type tiny = 1.e-09;
      return axom::utilities::clampVal(
        axom::utilities::abs(d0) / (axom::utilities::abs(d1 - d0) + tiny),
        value_type {0},
        value_type {1});
    }

    CoordsetView m_coordsetView;
    PlaneType m_plane;
  };

  /*!
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const TopologyView &AXOM_UNUSED_PARAM(topologyView),
                  const CoordsetView &coordsetView,
                  const conduit::Node &n_options,
                  const conduit::Node &AXOM_UNUSED_PARAM(n_topology),
                  const conduit::Node &AXOM_UNUSED_PARAM(n_coordset),
                  const conduit::Node &AXOM_UNUSED_PARAM(n_fields))
  {
    // Make a plane from the options.
    SLIC_ASSERT(n_options.has_child("origin"));
    SLIC_ASSERT(n_options.has_child("normal"));
    const auto origin = n_options["origin"].as_double_accessor();
    const auto normal = n_options["normal"].as_double_accessor();
    value_type planeOrigin[NDIMS], planeNormal[NDIMS];
    for(int i = 0; i < NDIMS; i++)
    {
      planeOrigin[i] = static_cast<value_type>(origin[i]);
      // Reverse the normal so generated surface's normals go in same direction as the normal.
      planeNormal[i] = -static_cast<value_type>(normal[i]);
    }

    // Set the plane in the view.
    m_view.m_plane = PlaneType(typename PlaneType::VectorType(planeNormal, NDIMS),
                               typename PlaneType::PointType(planeOrigin, NDIMS));

    // Save the coordset view.
    m_view.m_coordsetView = coordsetView;
  }

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

  /*!
   * \brief Return a new instance of the view.
   * \return A new instance of the view.
   */
  View view() const { return m_view; }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  View m_view {};
};

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom

#endif
