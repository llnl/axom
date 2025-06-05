// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_MAKE_UNSTRUCTURED_HPP_
#define AXOM_BUMP_MAKE_UNSTRUCTURED_HPP_

#include "axom/core.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/utilities/utilities.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/dispatch_structured_topology.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
/**
 * \brief Make an unstructured representation of a structured topology.
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 */
template <typename ExecSpace>
class MakeUnstructured
{
public:
  /**
   * \brief Make an unstructured representation of a structured topology.
   *
   * \tparam ExecSpace The execution space where the work will be done.
   *
   * \param topo      The input topology to be turned into unstructured.
   * \param coordset  The topology's coordset. It will be referenced as an external node in the output \a mesh.
   * \param topoName  The name of the new topology to create.
   * \param[out] mesh The node that will contain the new topology and coordset.
   *
   * \note There are blueprint methods for this sort of thing but this one runs on device.
   */
  static void execute(const conduit::Node &topo,
                      const conduit::Node &coordset,
                      const std::string &topoName,
                      conduit::Node &mesh)
  {
    const std::string type = topo.fetch_existing("type").as_string();
    namespace utils = axom::bump::utilities;
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    mesh["coordsets"][coordset.name()].set_external(coordset);
    conduit::Node &n_newtopo = mesh["topologies"][topoName];
    n_newtopo["coordset"] = coordset.name();

    if(type == "unstructured")
    {
      n_newtopo.set_external(topo);
    }
    else
    {
      n_newtopo["type"] = "unstructured";
      conduit::Node &n_newconn = n_newtopo["elements/connectivity"];
      conduit::Node &n_newsizes = n_newtopo["elements/sizes"];
      conduit::Node &n_newoffsets = n_newtopo["elements/offsets"];
      n_newconn.set_allocator(c2a.getConduitAllocatorID());
      n_newsizes.set_allocator(c2a.getConduitAllocatorID());
      n_newoffsets.set_allocator(c2a.getConduitAllocatorID());

      axom::bump::views::dispatch_structured_topologies(
        topo,
        [&](const std::string &shape, auto &topoView) {
          n_newtopo["elements/shape"] = shape;

          int ptsPerZone = 2;
          if(shape == "quad")
            ptsPerZone = 4;
          else if(shape == "hex")
            ptsPerZone = 8;

          // Allocate new mesh data.
          const auto nzones = topoView.numberOfZones();
          const auto connSize = nzones * ptsPerZone;
          n_newconn.set(conduit::DataType::index_t(connSize));
          n_newsizes.set(conduit::DataType::index_t(nzones));
          n_newoffsets.set(conduit::DataType::index_t(nzones));

          // Make views for the mesh data.
          auto connView = utils::make_array_view<conduit::index_t>(n_newconn);
          auto sizesView = utils::make_array_view<conduit::index_t>(n_newsizes);
          auto offsetsView = utils::make_array_view<conduit::index_t>(n_newoffsets);

          makeUnstructured(ptsPerZone, topoView, connView, sizesView, offsetsView);
        });
    }
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif

  /*!
   * \brief Iterate over the input topology's zones and store their connectivity in
   *        unstructured connectivity nodes.
   *
   * \param ptsPerZone The number of points per zone.
   * \param topoView The input topology view.
   * \param connView The view that will contain the new connectivity.
   * \param sizesView The view that will contain the new sizes.
   * \param offsetsView The view that will contain the new offsets.
   */
  template <typename TopologyView>
  static void makeUnstructured(int ptsPerZone,
                               TopologyView topoView,
                               axom::ArrayView<conduit::index_t> connView,
                               axom::ArrayView<conduit::index_t> sizesView,
                               axom::ArrayView<conduit::index_t> offsetsView)
  {
    // Fill in the new connectivity.
    axom::for_all<ExecSpace>(
      topoView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto zone = topoView.zone(zoneIndex);

        const auto start = zoneIndex * ptsPerZone;
        for(int i = 0; i < ptsPerZone; i++)
        {
          connView[start + i] = static_cast<conduit::index_t>(zone.getId(i));
        }
        sizesView[zoneIndex] = ptsPerZone;
        offsetsView[zoneIndex] = start;
      });
  }
};

}  // end namespace bump
}  // end namespace axom

#endif
