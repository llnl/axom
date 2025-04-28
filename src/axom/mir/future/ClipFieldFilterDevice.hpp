// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_CLIP_FIELD_FILTER_DEVICE_HPP_
#define AXOM_MIR_CLIP_FIELD_FILTER_DEVICE_HPP_

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/dispatch_topology.hpp"

namespace axom
{
namespace mir
{
namespace clipping
{
/**
 * \brief This class instantiates ClipField so it runs on specific device but
 *        can operate on a variety of Blueprint topologies/coordsets/types.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 */
template <typename ExecSpace>
class ClipFieldFilterDevice
{
public:
  /// Constructor
  ClipFieldFilterDevice() { }

  /**
   * \brief Execute the clipping operation using the specified options.
   *
   * \param[in] n_input The Conduit node that contains the topology, coordsets, and fields.
   * \param[in] n_options A Conduit node that contains clipping options.
   * \param[out] n_output A Conduit node that will hold the clipped output mesh. This should be a different node from \a n_input.
   *
   * \note The clipField field must currently be vertex-associated.
   */
  void execute(const conduit::Node &n_input, const conduit::Node &n_options, conduit::Node &n_output)
  {
    ClipOptions opts(n_options);
    const std::string clipFieldName = opts.clipField();

    const conduit::Node &n_fields = n_input.fetch_existing("fields");
    const conduit::Node &n_clipField = n_fields.fetch_existing(clipFieldName);
    const std::string &topoName = n_clipField["topology"].as_string();
    const conduit::Node &n_topo = n_input.fetch_existing("topologies/" + topoName);
    const std::string &coordsetName = n_topo["coordset"].as_string();
    const conduit::Node &n_coordset = n_input.fetch_existing("coordsets/" + coordsetName);

    execute(n_topo,
            n_coordset,
            n_fields,
            n_options,
            n_output["topologies/" + opts.topologyName(topoName)],
            n_output["coordsets/" + opts.coordsetName(coordsetName)],
            n_output["fields"]);
  }

  /**
   * \brief Returns a bitset that records which shapes are supported.
   * \return A bitset that encodes the supported shapes.
   */
  static constexpr int supported_shapes()
  {
    int bitset = 0;
    axom::utilities::setBitOn(bitset, views::Tri_ShapeID);
    axom::utilities::setBitOn(bitset, views::Quad_ShapeID);
    axom::utilities::setBitOn(bitset, views::Pyramid_ShapeID);
    axom::utilities::setBitOn(bitset, views::Wedge_ShapeID);
    axom::utilities::setBitOn(bitset, views::Hex_ShapeID);
    axom::utilities::setBitOn(bitset, views::Mixed_ShapeID);
    return bitset;
  }

  /**
   * \brief Execute the clipping operation using the specified options.
   *
   * \param[in] n_topo The node that contains the input mesh topology.
   * \param[in] n_coordset The node that contains the input mesh coordset.
   * \param[in] n_fields The node that contains the input fields.
   * \param[in] n_options A Conduit node that contains clipping options.
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   *
   * \note The clipField field must currently be vertex-associated. Also, the output topology will be an unstructured topology with mixed shape types.
   */
  void execute(const conduit::Node &n_topo,
               const conduit::Node &n_coordset,
               const conduit::Node &n_fields,
               const conduit::Node &n_options,
               conduit::Node &n_newTopo,
               conduit::Node &n_newCoordset,
               conduit::Node &n_newFields)
  {
#if 0
    // NOTE - there are 2 dispatches here so we can get coordset and topology views.
    //        This instantiates the lambda that creates the ClipField object for
    //        the various view types so we can handle most Blueprint topologies.
    //        However, this expansion makes this file take a long time to build.
    //
    //        Perhaps we could split the topology dispatch somehow to split the
    //        compilation responsibility among more cores for a faster build.
    //
    axom::mir::views::dispatch_coordset(n_coordset, [&](auto coordsetView) {
      // We dispatch to 2D/3D topologies with supported shapes.
      constexpr int dims = axom::mir::views::select_dimensions(2, 3);
      constexpr int shapes = supported_shapes();
      axom::mir::views::dispatch_topology<dims, shapes>(
        n_topo,
        [&](const std::string & /*shape*/, auto topologyView) {
          using TopologyView = decltype(topologyView);
          using CoordsetView = decltype(coordsetView);

          // Don't allow 3D topologies with 2D coordsets.
          constexpr int RuntimeDimension = -1;
          if constexpr((TopologyView::dimension() == 2) ||
                       (TopologyView::dimension() == 3 &&
                        CoordsetView::dimension() == 3) ||
                       (TopologyView::dimension() == RuntimeDimension))
          {
            ClipField<ExecSpace, TopologyView, CoordsetView> clipper(
              topologyView,
              coordsetView);
            clipper.execute(n_topo,
                            n_coordset,
                            n_fields,
                            n_options,
                            n_newTopo,
                            n_newCoordset,
                            n_newFields);
          }
        });
    });
#endif
  }
};

}  // end namespace clipping
}  // end namespace mir
}  // end namespace axom

#endif
