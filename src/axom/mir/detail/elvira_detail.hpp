// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_ELVIRA_ALGORITHM_DETAIL_HPP_
#define AXOM_MIR_ELVIRA_ALGORITHM_DETAIL_HPP_

// Includes happen in the ElviraAlgorithm.hpp header file that includes this file.

namespace axom
{
namespace mir
{
namespace detail
{

/*!
 * \brief Compute a range given by 2 points given a shape and a normal. The
 *        range represents a range such that clipping the shape at range[0]
 *        would yield 0% intersection and clipping the shape at range[1]
 *        would yield 100% intersection.
 *
 * \param shape The shape whose range is being computed.
 * \param normal The normal in the shape along which a clip plane would move.
 * \param[out] range Start and end points of the range that can be used to
 *                   make a clipping plane for the shape.
 */
template <typename ShapeType, typename T, int NDIMS>
AXOM_HOST_DEVICE inline void computeRange(const ShapeType &shape,
                                          const axom::primal::Vector<T, NDIMS> &normal,
                                          axom::primal::Point<T, NDIMS> range[2])
{
  // Compute the shape bounding box.
  const auto bbox = axom::primal::compute_bounding_box(shape);

  const axom::primal::Plane<T, NDIMS> P(normal, bbox.getCentroid(), false);

  // Compute distances from all points in shape to plane.
  range[0] = range[1] = bbox.getCentroid();
  double dist[2] = {0., 0.};
  for(axom::IndexType ip = 0; ip < shape.numVertices(); ip++)
  {
    const auto d = static_cast<T>(P.signedDistance(shape[ip]));
    if(d < dist[0])
    {
      dist[0] = d;
      range[0] = bbox.getCentroid() + (normal * d);
    }
    if(d > dist[1])
    {
      dist[1] = d;
      range[1] = bbox.getCentroid() + (normal * d);
    }
  }
}

/*!
 * \brief Clip a shape until a specified volume is achieved, up to a tolerance
 *        or max number of iterations. The clipped shape is returned. The
 *        origin of the clipping plane is returned as an out argument.
 *
 * \param shape The shape to be clipped.
 * \param normal The normal along which the clipping plane will move.
 * \param _range The start, end points of a range that contains the plane origin.
 * \param matVolume The target shape volume.
 * \param max_iterations The max number of iterations to get the answer.
 * \param tolerance The allowable volume tolerance.
 * \param[out] pt The origin of the clipping plane that was used.
 */
template <typename ClipResultType, typename ShapeType, typename T, int NDIMS>
AXOM_HOST_DEVICE inline ClipResultType clipToVolume(const ShapeType &shape,
                                                    const axom::primal::Vector<T, NDIMS> &normal,
                                                    const axom::primal::Point<T, NDIMS> _range[2],
                                                    double matVolume,
                                                    int max_iterations,
                                                    double tolerance,
                                                    axom::primal::Point<T, NDIMS> &pt)
{
  namespace bputils = axom::mir::utilities::blueprint;
  // The range for the interval
  axom::primal::Point<T, NDIMS> range[2] = {_range[0], _range[1]};
  // This array holds the volumes for the interval.
  double f_t[2] = {0., matVolume};
  // The blend value within the current interval.
  double t_blend = 0.5;

  // The ELVIRA normals point away from the material. Axom's clipping
  // keeps the shape where the normal points into the shape. Reverse
  // the ELVIRA normal when forming the plane so we keep the right piece.
  const auto clipNormal = -normal;

  ClipResultType clippedShape {};
  for(int iterations = 0; iterations < max_iterations; iterations++)
  {
    // Pick the middle of the range and position the plane there.
    pt = axom::primal::Point<T, NDIMS>::lerp(range[0], range[1], t_blend);
    const auto P = axom::primal::Plane<T, NDIMS>(clipNormal, pt, false);

    // Clip the shape at the current plane.
    clippedShape = axom::primal::clip(shape, P);

    // Find the volume of the clipped shape.
    const double fragmentVolume = bputils::ComputeShapeAmount<NDIMS>::execute(clippedShape);
    const double volumeError = axom::utilities::abs(matVolume - fragmentVolume);

    //std::cout << "\titerations=" << iterations << ", t_blend=" << t_blend << ", P=" << P << ", clippedShape=" << clippedShape << ", matVolume=" << matVolume << ", fragmentVolume=" << fragmentVolume << ", volumeError=" << volumeError << std::endl;
    if((volumeError <= tolerance) || (iterations >= max_iterations))
    {
      break;
    }
    else if(fragmentVolume < matVolume)
    {
      range[0] = pt;
      f_t[0] = fragmentVolume;
    }
    else
    {
      range[1] = pt;
      f_t[1] = fragmentVolume;
    }

    // Now, try and find a new value for t_blend, the value within this interval
    // where we think we'll find targetVolume
    //
    //  |           f_t[0]
    //  |           *\___
    //  |           |    \____*targetVolume
    //  |                     |\_______* f_t[1]
    //  |           |
    //  |                     |        |
    //  ------------+---------+--------+--------------
    //              0         t_blend  1
    constexpr double offset = 0.02;
    t_blend = axom::utilities::clampVal((matVolume - f_t[0]) / (f_t[1] - f_t[0]), 0., 1.);
    t_blend = (1. - 2. * offset) * t_blend + offset;
  }

  return clippedShape;
}

/*!
 * \brief Base template for building a new Conduit mesh.
 */
template <typename ExecSpace, typename CoordsetView, typename TopologyView, typename MatsetView, typename PolygonShape, int NDIMS>
struct TopologyBuilder
{ };

/*!
 * \brief Template specialization for making a new 2D Conduit polygon mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 * \tparam CoordsetView The view that wraps the coordset data.
 * \tparam TopologyView The view that wraps the topology data.
 * \tparam MatsetView The view that wraps the matset data.
 * \tparam PolygonShape An axom::primal::Polygon (with various template args filled in).
 */
template <typename ExecSpace, typename CoordsetView, typename TopologyView, typename MatsetView, typename PolygonShape>
class TopologyBuilder<ExecSpace, CoordsetView, TopologyView, MatsetView, PolygonShape, 2>
{
  // The way we build fragments from quads, we should not get any with more than 5 sides.
  static constexpr int MAX_POINTS_PER_FRAGMENT = 5;

  using CoordType = typename CoordsetView::value_type;
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using MaterialID = typename MatsetView::MaterialID;
  using MaterialVF = typename MatsetView::FloatType;

public:
  /*!
   * \brief Create the Blueprint nodes for the output coordset, topology, fields, and matset.
   *
   * \note This is a host-only function. We allocate fixed size arrays for the coordset and
   *       connectivity that will have gaps that consumers will need to skip over using
   *       the provided sizes/offsets.
   */
  void allocate(axom::IndexType numFragments,
                conduit::Node &n_coordset,
                conduit::Node &n_topology,
                conduit::Node &n_fields,
                conduit::Node &n_matset)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const auto numCoordValues = numFragments * MAX_POINTS_PER_FRAGMENT;

    // Set up coordset and allocate data arrays.
    // Note that we overallocate the number of nodes to numCoordValues.
    n_coordset["type"] = "explicit";
    n_coordset["values/x"].set_allocator(c2a.getConduitAllocatorID());
    n_coordset["values/x"].set(conduit::DataType(bputils::cpp2conduit<CoordType>::id, numCoordValues));
    m_view.m_x = bputils::make_array_view<CoordType>(n_coordset["values/x"]);
    n_coordset["values/y"].set_allocator(c2a.getConduitAllocatorID());
    n_coordset["values/y"].set(conduit::DataType(bputils::cpp2conduit<CoordType>::id, numCoordValues));
    m_view.m_y = bputils::make_array_view<CoordType>(n_coordset["values/y"]);

    axom::mir::utilities::fill<ExecSpace>(m_view.m_x, CoordType(0));
    axom::mir::utilities::fill<ExecSpace>(m_view.m_y, CoordType(0));

    // Set up connectivity and allocate data arrays.
    n_topology["type"] = "unstructured";
    n_topology["elements/shape"] = "polygonal";
    conduit::Node &n_conn = n_topology["elements/connectivity"];
    n_conn.set_allocator(c2a.getConduitAllocatorID());
    n_conn.set(conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id,
                                 numFragments * MAX_POINTS_PER_FRAGMENT));
    m_view.m_connectivity = bputils::make_array_view<ConnectivityType>(n_conn);
    axom::mir::utilities::fill<ExecSpace>(m_view.m_connectivity, ConnectivityType(0));

    conduit::Node &n_sizes = n_topology["elements/sizes"];
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set(conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, numFragments));
    m_view.m_sizes = bputils::make_array_view<ConnectivityType>(n_sizes);

    conduit::Node &n_offsets = n_topology["elements/offsets"];
    n_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set(conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, numFragments));
    m_view.m_offsets = bputils::make_array_view<ConnectivityType>(n_offsets);

    // Make new fields.
    n_fields["originalElements/topology"] = n_topology.name();
    n_fields["originalElements/association"] = "element";
    conduit::Node &n_orig_zones = n_fields["originalElements/values"];
    n_orig_zones.set_allocator(c2a.getConduitAllocatorID());
    n_orig_zones.set(conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, numFragments));
    m_view.m_original_zones = bputils::make_array_view<ConnectivityType>(n_orig_zones);

    conduit::Node &n_normal = n_fields["normal"];
    n_normal["topology"] = n_topology.name();
    n_normal["association"] = "element";
    conduit::Node &n_x = n_normal["values/x"];
    conduit::Node &n_y = n_normal["values/y"];
    n_x.set_allocator(c2a.getConduitAllocatorID());
    n_x.set(conduit::DataType(bputils::cpp2conduit<double>::id, numFragments));
    m_view.m_norm_x = bputils::make_array_view<double>(n_x);
    n_y.set_allocator(c2a.getConduitAllocatorID());
    n_y.set(conduit::DataType(bputils::cpp2conduit<double>::id, numFragments));
    m_view.m_norm_y = bputils::make_array_view<double>(n_y);

    // Set up new matset. All of the sizes are numFragments because we're making clean zones.
    n_matset["topology"] = n_topology.name();
    conduit::Node &n_volume_fractions = n_matset["volume_fractions"];
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_volume_fractions.set(conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, numFragments));
    m_view.m_volume_fractions = bputils::make_array_view<MaterialVF>(n_volume_fractions);

    conduit::Node &n_material_ids = n_matset["material_ids"];
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_material_ids.set(conduit::DataType(bputils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_material_ids = bputils::make_array_view<MaterialID>(n_material_ids);

    conduit::Node &n_indices = n_matset["indices"];
    n_indices.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set(conduit::DataType(bputils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_indices = bputils::make_array_view<MaterialID>(n_indices);

    conduit::Node &n_mat_sizes = n_matset["sizes"];
    n_mat_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_mat_sizes.set(conduit::DataType(bputils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_sizes = bputils::make_array_view<MaterialID>(n_mat_sizes);
    axom::mir::utilities::fill<ExecSpace>(m_view.m_mat_sizes, MaterialID(0));

    conduit::Node &n_mat_offsets = n_matset["offsets"];
    n_mat_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_mat_offsets.set(conduit::DataType(bputils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_offsets = bputils::make_array_view<MaterialID>(n_mat_offsets);
  }

  /*!
   * \brief A view that can be used on device to encode a shape into Blueprint arrays.
   */
  struct View
  {
    /*!
     * \brief Encode a polygon from MIR into the views that represent the Conduit data.
     *
     * \param zoneIndex The original zone that yielded the fragment.
     * \param fragmentOffset The offset into the per-fragment data.
     * \param shape The 2D primal polygon we're encoding.
     * \param matId The material id for the fragment.
     * \param normal The normal for the fragment.
     */
    AXOM_HOST_DEVICE
    void addShape(axom::IndexType zoneIndex,
                  axom::IndexType fragmentOffset,
                  const PolygonShape &shape,
                  int matId,
                  const double *normal) const
    {
      const int nverts = shape.numVertices();
      SLIC_ASSERT(nverts <= MAX_POINTS_PER_FRAGMENT);

      // Copy coordinates into coordinate arrays. We might end up with unreferenced coordinates for now.
      auto coordOffset = fragmentOffset * MAX_POINTS_PER_FRAGMENT;
      for(int i = 0; i < nverts; i++)
      {
        m_x[coordOffset + i] = shape[i][0];
        m_y[coordOffset + i] = shape[i][1];
      }

      // Make connectivity.
      auto connOffset = fragmentOffset * MAX_POINTS_PER_FRAGMENT;
      for(int i = 0; i < nverts; i++)
      {
        auto idx = static_cast<ConnectivityType>(connOffset + i);
        m_connectivity[idx] = idx;
      }
      m_sizes[fragmentOffset] = static_cast<ConnectivityType>(nverts);
      m_offsets[fragmentOffset] = connOffset;

      // Save fields.
      m_original_zones[fragmentOffset] = zoneIndex;
      m_norm_x[fragmentOffset] = normal[0];
      m_norm_y[fragmentOffset] = normal[1];

      // Save material data.
      m_volume_fractions[fragmentOffset] = static_cast<MaterialVF>(1.);
      m_material_ids[fragmentOffset] = static_cast<MaterialID>(matId);
      m_mat_sizes[fragmentOffset] = static_cast<MaterialID>(1);
      m_mat_offsets[fragmentOffset] = fragmentOffset;
      m_mat_indices[fragmentOffset] = fragmentOffset;
    }

    // These ArrayView objects expose data that we've allocated in Conduit nodes to represent the new mesh.

    // New coordset data.
    axom::ArrayView<CoordType> m_x {}, m_y {};  //!< X,Y coordinates

    // New topology data.
    axom::ArrayView<ConnectivityType> m_connectivity {}, m_sizes {},
      m_offsets {};  //!< Connectivity data.

    // New field data.
    axom::ArrayView<axom::IndexType> m_original_zones {};  //!< View for originalZone field data.
    axom::ArrayView<double> m_norm_x {}, m_norm_y {};      //!< Fragment normals

    // New matset data
    axom::ArrayView<MaterialVF> m_volume_fractions {};
    axom::ArrayView<MaterialID> m_material_ids {}, m_mat_indices, m_mat_sizes {}, m_mat_offsets {};
  };

  /*!
   * \brief Return a view that can be copied to device.
   *
   * \return A view that be used to store data.
   */
  View view() { return m_view; }

private:
  View m_view;
};

}  // end namespace detail
}  // end namespace mir
}  // end namespace axom

#endif
