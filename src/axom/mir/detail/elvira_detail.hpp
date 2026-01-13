// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_ELVIRA_ALGORITHM_DETAIL_HPP_
#define AXOM_MIR_ELVIRA_ALGORITHM_DETAIL_HPP_

// Most includes happen in the ElviraAlgorithm.hpp header file that includes this file.

#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/ExtractZones.hpp"
#include "axom/bump/ExtractZonesAndMatsetPolyhedral.hpp"
#include "axom/bump/MergeCoordsetPoints.hpp"
#include "axom/bump/MergePolyhedralFaces.hpp"
#include "axom/bump/views/dispatch_coordset.hpp"

#include <sstream>
#include <iomanip>

namespace axom
{
namespace mir
{

//------------------------------------------------------------------------------
/**
 * \brief This class provides a kind of schema over the MIR options, as well
 *        as default values, and some utilities functions.
 */
class ELVIRAOptions : public axom::bump::Options
{
public:
  /**
   * \brief Constructor
   *
   * \param options The node that contains the clipping options.
   */
  ELVIRAOptions(const conduit::Node &options) : axom::bump::Options(options) { }

  /**
   * \brief Get whether the plane equation fields should appear in the output.
   * \return True if plane equation fields should appear, false otherwise.
   */
  bool plane() const { return flagValue("plane"); }

  /**
   * \brief Get whether the output should be a point mesh.
   * \return True if the output should be a point mesh, false otherwise.
   */
  bool pointmesh() const { return flagValue("pointmesh"); }

protected:
  /**
   * \brief Get whether the flag is set in the options.
   *
   * \param key The name of the key that contains the flag.
   *
   * \return True if key is present and set to non-zero, false otherwise.
   */
  bool flagValue(const std::string &key) const
  {
    bool retval = false;
    if(options().has_path(key))
    {
      retval = options().fetch_existing(key).to_int() != 0;
    }
    return retval;
  }
};

namespace detail
{

//------------------------------------------------------------------------------
template <typename T>
struct clip_precision
{
  // The default precision used in axom::primal::clip
  static constexpr double eps = 1.e-10;
};

template <>
struct clip_precision<double>
{
  // Make the axom::primal::clip more precise for double.
  static constexpr double eps = 1.e-15;
};

//------------------------------------------------------------------------------
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
  const auto centroid = bbox.getCentroid();
  range[0] = range[1] = centroid;
  double dist[2] = {0., 0.};
  for(axom::IndexType ip = 0; ip < shape.numVertices(); ip++)
  {
    const auto d = static_cast<T>(P.signedDistance(shape[ip]));
    if(d < dist[0])
    {
      dist[0] = d;
      range[0] = centroid + (normal * d);
    }
    if(d > dist[1])
    {
      dist[1] = d;
      range[1] = centroid + (normal * d);
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
  namespace utils = axom::bump::utilities;
  // The range for the interval
  axom::primal::Point<T, NDIMS> range[2] = {_range[0], _range[1]};
  // This array holds the volumes for the interval.
  double f_t[2] = {0., matVolume};
  // The blend value within the current interval.
  double t_blend = 0.5;

#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
  SLIC_DEBUG("\tclipToVolume: max_iterations=" << max_iterations << ", tolerance=" << tolerance
                                               << ", shape=" << shape);
#endif

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
    clippedShape = axom::primal::clip(shape, P, clip_precision<T>::eps);

    // Find the volume of the clipped shape.
    const double fragmentVolume = utils::ComputeShapeAmount<NDIMS>::execute(clippedShape);
    const double volumeError = axom::utilities::abs(matVolume - fragmentVolume);
#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
    if(ShapeType::PointType::DIMENSION < 3)
    {
      SLIC_DEBUG("\t\titerations=" << iterations << ", t_blend=" << t_blend << ", P=" << P
                                   << ", clippedShape=" << clippedShape
                                   << ", matVolume=" << std::setprecision(16) << matVolume
                                   << ", fragmentVolume=" << std::setprecision(16) << fragmentVolume
                                   << ", volumeError=" << std::setprecision(16) << volumeError);
    }
    else
    {
      SLIC_DEBUG("\t\titerations=" << iterations << ", t_blend=" << t_blend << ", P=" << P
                                   << ", matVolume=" << std::setprecision(16) << matVolume
                                   << ", fragmentVolume=" << std::setprecision(16) << fragmentVolume
                                   << ", volumeError=" << std::setprecision(16) << volumeError);
    }
#endif
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
    constexpr double offset = (NDIMS == 3) ? 0.1 : 0.02;
    t_blend = axom::utilities::clampVal((matVolume - f_t[0]) / (f_t[1] - f_t[0]), 0., 1.);
    t_blend = (1. - 2. * offset) * t_blend + offset;
  }

  return clippedShape;
}

/*!
 * \brief Base template for building a new Conduit mesh.
 */
template <typename ExecSpace, typename CoordsetView, typename TopologyView, typename MatsetView, typename Shape, int NDIMS>
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
  using CoordType = typename CoordsetView::value_type;
  using PointType = typename axom::primal::Point<CoordType, 2>;
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using MaterialID = typename MatsetView::MaterialID;
  using MaterialVF = typename MatsetView::FloatType;

public:
  /*!
   * \brief Create the Blueprint nodes for the output coordset, topology, fields, and matset.
   *
   * \param numFragments The total number of fragments that will be created for MIR.
   *                     This is the total number of fragments in all zones.
   * \param maxCuts The max number of cuts that will happen in any zone.
   * \param[out] n_coordset The Conduit node in which to build the new coordset.
   * \param[out] n_topology The Conduit node in which to build the new topology.
   * \param[out] n_fields The Conduit node in which to build the new fields.
   * \param[out] n_matset The Conduit node in which to build the new matset.
   * \param n_options A Conduit node that contains MIR options.
   *
   * \note This is a host-only function. We allocate fixed size arrays for the coordset and
   *       connectivity that will have gaps that consumers will need to skip over using
   *       the provided sizes/offsets.
   */
  void allocate(axom::IndexType numFragments,
                axom::IndexType maxCuts,
                conduit::Node &n_coordset,
                conduit::Node &n_topology,
                conduit::Node &n_fields,
                conduit::Node &n_matset,
                const conduit::Node &n_options)
  {
    namespace utils = axom::bump::utilities;
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Handle options
    ELVIRAOptions opts(n_options);
    const std::string originalElementsField(opts.originalElementsField());
    m_view.m_makePlane = opts.plane();
    m_view.m_makePointMesh = opts.pointmesh();

    // Figure out the max fragment size.
    m_view.m_MAX_POINTS_PER_FRAGMENT = 4 + maxCuts;

    // Vary the number of coord values depending on whether or not we're making
    // a point mesh.
    const auto numCoordValues =
      m_view.m_makePointMesh ? numFragments : (numFragments * m_view.m_MAX_POINTS_PER_FRAGMENT);

    // Set up coordset and allocate data arrays.
    // Note that we overallocate the number of nodes to numCoordValues.
    n_coordset["type"] = "explicit";
    n_coordset["values/x"].set_allocator(c2a.getConduitAllocatorID());
    n_coordset["values/x"].set(conduit::DataType(utils::cpp2conduit<CoordType>::id, numCoordValues));
    m_view.m_x = utils::make_array_view<CoordType>(n_coordset["values/x"]);
    n_coordset["values/y"].set_allocator(c2a.getConduitAllocatorID());
    n_coordset["values/y"].set(conduit::DataType(utils::cpp2conduit<CoordType>::id, numCoordValues));
    m_view.m_y = utils::make_array_view<CoordType>(n_coordset["values/y"]);

    m_view.m_x.fill(CoordType(0));
    m_view.m_y.fill(CoordType(0));

    // Set up connectivity and allocate data arrays.
    n_topology["type"] = "unstructured";
    n_topology["elements/shape"] = m_view.m_makePointMesh ? "point" : "polygonal";
    conduit::Node &n_conn = n_topology["elements/connectivity"];
    n_conn.set_allocator(c2a.getConduitAllocatorID());
    n_conn.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numCoordValues));
    m_view.m_connectivity = utils::make_array_view<ConnectivityType>(n_conn);
    m_view.m_connectivity.fill(ConnectivityType(0));

    conduit::Node &n_sizes = n_topology["elements/sizes"];
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numFragments));
    m_view.m_sizes = utils::make_array_view<ConnectivityType>(n_sizes);

    conduit::Node &n_offsets = n_topology["elements/offsets"];
    n_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numFragments));
    m_view.m_offsets = utils::make_array_view<ConnectivityType>(n_offsets);

    // Make new fields.
    conduit::Node &n_origElem = n_fields[originalElementsField];
    n_origElem["topology"] = n_topology.name();
    n_origElem["association"] = "element";
    conduit::Node &n_orig_elem_values = n_origElem["values"];
    n_orig_elem_values.set_allocator(c2a.getConduitAllocatorID());
    n_orig_elem_values.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numFragments));
    m_view.m_original_zones = utils::make_array_view<ConnectivityType>(n_orig_elem_values);

    if(m_view.m_makePlane)
    {
      conduit::Node &n_normal = n_fields["normal"];
      n_normal["topology"] = n_topology.name();
      n_normal["association"] = "element";
      conduit::Node &n_x = n_normal["values/x"];
      conduit::Node &n_y = n_normal["values/y"];
      n_x.set_allocator(c2a.getConduitAllocatorID());
      n_x.set(conduit::DataType(utils::cpp2conduit<double>::id, numFragments));
      m_view.m_norm_x = utils::make_array_view<double>(n_x);
      n_y.set_allocator(c2a.getConduitAllocatorID());
      n_y.set(conduit::DataType(utils::cpp2conduit<double>::id, numFragments));
      m_view.m_norm_y = utils::make_array_view<double>(n_y);

      conduit::Node &n_planeOffset = n_fields["offset"];
      n_planeOffset["topology"] = n_topology.name();
      n_planeOffset["association"] = "element";
      conduit::Node &n_values = n_planeOffset["values"];
      n_values.set_allocator(c2a.getConduitAllocatorID());
      n_values.set(conduit::DataType(utils::cpp2conduit<double>::id, numFragments));
      m_view.m_plane_offset = utils::make_array_view<double>(n_values);
    }

    // Set up new matset. All of the sizes are numFragments because we're making clean zones.
    n_matset["topology"] = n_topology.name();
    conduit::Node &n_volume_fractions = n_matset["volume_fractions"];
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_volume_fractions.set(conduit::DataType(utils::cpp2conduit<MaterialVF>::id, numFragments));
    m_view.m_volume_fractions = utils::make_array_view<MaterialVF>(n_volume_fractions);

    conduit::Node &n_material_ids = n_matset["material_ids"];
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_material_ids.set(conduit::DataType(utils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_material_ids = utils::make_array_view<MaterialID>(n_material_ids);

    conduit::Node &n_indices = n_matset["indices"];
    n_indices.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set(conduit::DataType(utils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_indices = utils::make_array_view<MaterialID>(n_indices);

    conduit::Node &n_mat_sizes = n_matset["sizes"];
    n_mat_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_mat_sizes.set(conduit::DataType(utils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_sizes = utils::make_array_view<MaterialID>(n_mat_sizes);
    m_view.m_mat_sizes.fill(MaterialID(0));

    conduit::Node &n_mat_offsets = n_matset["offsets"];
    n_mat_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_mat_offsets.set(conduit::DataType(utils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_offsets = utils::make_array_view<MaterialID>(n_mat_offsets);
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
     * \param pt The origin for the fragment's cutting plane.
     * \param plane The offset for the fragment's cutting plane.
     * \param planeNormal The normal for the fragment's cutting plane.
     */
    AXOM_HOST_DEVICE
    void addShape(axom::IndexType zoneIndex,
                  axom::IndexType fragmentOffset,
                  const PolygonShape &shape,
                  int matId,
                  const PointType &pt,
                  double planeOffset,
                  const double *planeNormal) const
    {
      const int nverts = shape.numVertices();
      SLIC_ASSERT(nverts <= m_MAX_POINTS_PER_FRAGMENT);

      if(m_makePointMesh)
      {
        // Prefer the actual plane origin.
        m_x[fragmentOffset] = pt[0];
        m_y[fragmentOffset] = pt[1];

        // Make connectivity.
        m_connectivity[fragmentOffset] = fragmentOffset;
        m_sizes[fragmentOffset] = 1;
        m_offsets[fragmentOffset] = fragmentOffset;
      }
      else
      {
        // Copy coordinates into coordinate arrays. We might end up with unreferenced coordinates for now.
        auto coordOffset = fragmentOffset * m_MAX_POINTS_PER_FRAGMENT;
        for(int i = 0; i < nverts; i++)
        {
          m_x[coordOffset + i] = shape[i][0];
          m_y[coordOffset + i] = shape[i][1];
        }

        // Make connectivity.
        auto connOffset = fragmentOffset * m_MAX_POINTS_PER_FRAGMENT;
        for(int i = 0; i < nverts; i++)
        {
          auto idx = static_cast<ConnectivityType>(connOffset + i);
          m_connectivity[idx] = idx;
        }
        m_sizes[fragmentOffset] = static_cast<ConnectivityType>(nverts);
        m_offsets[fragmentOffset] = connOffset;
      }
      // Save fields.
      m_original_zones[fragmentOffset] = zoneIndex;
      if(m_makePlane)
      {
        m_norm_x[fragmentOffset] = planeNormal[0];
        m_norm_y[fragmentOffset] = planeNormal[1];

        m_plane_offset[fragmentOffset] = planeOffset;
      }

      // Save material data.
      m_volume_fractions[fragmentOffset] = static_cast<MaterialVF>(1.);
      m_material_ids[fragmentOffset] = static_cast<MaterialID>(matId);
      m_mat_sizes[fragmentOffset] = static_cast<MaterialID>(1);
      m_mat_offsets[fragmentOffset] = fragmentOffset;
      m_mat_indices[fragmentOffset] = fragmentOffset;
#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
      SLIC_DEBUG("\taddShape: zone=" << zoneIndex << ", fragmentOffset=" << fragmentOffset
                                     << ", mat=" << matId << ", shape=" << shape);
#endif
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
    axom::ArrayView<double> m_plane_offset;                //!< Fragment plane offset

    // New matset data
    axom::ArrayView<MaterialVF> m_volume_fractions {};
    axom::ArrayView<MaterialID> m_material_ids {}, m_mat_indices, m_mat_sizes {}, m_mat_offsets {};

    bool m_makePlane {false};
    bool m_makePointMesh {false};

    axom::IndexType m_MAX_POINTS_PER_FRAGMENT {0};
  };

  /*!
   * \brief Return a view that can be copied to device.
   *
   * \return A view that be used to store data.
   */
  View view() { return m_view; }

  /*!
   * \brief Clean the mesh (no-op for 2D)
   */
  void cleanMesh(conduit::Node &AXOM_UNUSED_PARAM(n_coordset),
                 double AXOM_UNUSED_PARAM(point_tolerance),
                 conduit::Node &AXOM_UNUSED_PARAM(n_topology),
                 axom::Array<axom::IndexType> &AXOM_UNUSED_PARAM(selectedIds)) const
  { }

private:
  View m_view;
};

/*!
 * \brief Template specialization for making a new 3D Conduit polyhedral mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 * \tparam CoordsetView The view that wraps the coordset data.
 * \tparam TopologyView The view that wraps the topology data.
 * \tparam MatsetView The view that wraps the matset data.
 * \tparam PHShape An axom::primal::Polyhedron (with various template args filled in).
 */
template <typename ExecSpace, typename CoordsetView, typename TopologyView, typename MatsetView, typename PHShape>
class TopologyBuilder<ExecSpace, CoordsetView, TopologyView, MatsetView, PHShape, 3>
{
  using CoordType = typename CoordsetView::value_type;
  using PointType = typename axom::primal::Point<CoordType, 3>;
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using MaterialID = typename MatsetView::MaterialID;
  using MaterialVF = typename MatsetView::FloatType;

public:
  /*!
   * \brief Create the Blueprint nodes for the output coordset, topology, fields, and matset.
   *
   * \param numFragments The total number of fragments that will be created for MIR.
   *                     This is the total number of fragments in all zones.
   * \param maxCuts The max number of cuts that will happen in any zone.
   * \param[out] n_coordset The Conduit node in which to build the new coordset.
   * \param[out] n_topology The Conduit node in which to build the new topology.
   * \param[out] n_fields The Conduit node in which to build the new fields.
   * \param[out] n_matset The Conduit node in which to build the new matset.
   * \param n_options A Conduit node that contains MIR options.
   *
   * \note This is a host-only function. We allocate fixed size arrays for the coordset and
   *       connectivity that will have gaps that consumers will need to skip over using
   *       the provided sizes/offsets.
   */
  void allocate(axom::IndexType numFragments,
                axom::IndexType maxCuts,
                conduit::Node &n_coordset,
                conduit::Node &n_topology,
                conduit::Node &n_fields,
                conduit::Node &n_matset,
                const conduit::Node &n_options)
  {
    namespace utils = axom::bump::utilities;
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Handle options
    ELVIRAOptions opts(n_options);
    const std::string originalElementsField(opts.originalElementsField());
    m_view.m_makePlane = opts.plane();
    m_view.m_makePointMesh = opts.pointmesh();

    // Figure out some fragment size information given maxCuts, the max number of
    // times a zone will be cut.
    m_view.m_MAX_POINTS_PER_FACE = 6 + maxCuts;
    m_view.m_MAX_FACES_PER_FRAGMENT = 6 + maxCuts;
    m_view.m_MAX_POINTS_PER_FRAGMENT = 8 + maxCuts * 2;

    // Set up coordset and allocate data arrays.
    // Note that we overallocate the number of nodes to numCoordValues.
    const auto numCoordValues =
      m_view.m_makePointMesh ? numFragments : (numFragments * m_view.m_MAX_POINTS_PER_FRAGMENT);
    n_coordset["type"] = "explicit";
    n_coordset["values/x"].set_allocator(c2a.getConduitAllocatorID());
    n_coordset["values/x"].set(conduit::DataType(utils::cpp2conduit<CoordType>::id, numCoordValues));
    m_view.m_x = utils::make_array_view<CoordType>(n_coordset["values/x"]);
    n_coordset["values/y"].set_allocator(c2a.getConduitAllocatorID());
    n_coordset["values/y"].set(conduit::DataType(utils::cpp2conduit<CoordType>::id, numCoordValues));
    m_view.m_y = utils::make_array_view<CoordType>(n_coordset["values/y"]);
    n_coordset["values/z"].set_allocator(c2a.getConduitAllocatorID());
    n_coordset["values/z"].set(conduit::DataType(utils::cpp2conduit<CoordType>::id, numCoordValues));
    m_view.m_z = utils::make_array_view<CoordType>(n_coordset["values/z"]);

    m_view.m_x.fill(CoordType(0));
    m_view.m_y.fill(CoordType(0));
    m_view.m_z.fill(CoordType(0));

    // elements (zone definitions)
    constexpr ConnectivityType UnusedValue = axom::numeric_limits<ConnectivityType>::is_signed
      ? -1
      : axom::numeric_limits<ConnectivityType>::max();
    {
      const auto numConnValues =
        m_view.m_makePointMesh ? numFragments : (numFragments * m_view.m_MAX_FACES_PER_FRAGMENT);
      n_topology["type"] = "unstructured";
      n_topology["elements/shape"] = m_view.m_makePointMesh ? "point" : "polyhedral";
      conduit::Node &n_conn = n_topology["elements/connectivity"];
      n_conn.set_allocator(c2a.getConduitAllocatorID());
      n_conn.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numConnValues));
      m_view.m_connectivity = utils::make_array_view<ConnectivityType>(n_conn);
      m_view.m_connectivity.fill(UnusedValue);

      conduit::Node &n_sizes = n_topology["elements/sizes"];
      n_sizes.set_allocator(c2a.getConduitAllocatorID());
      n_sizes.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numFragments));
      m_view.m_sizes = utils::make_array_view<ConnectivityType>(n_sizes);

      conduit::Node &n_offsets = n_topology["elements/offsets"];
      n_offsets.set_allocator(c2a.getConduitAllocatorID());
      n_offsets.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numFragments));
      m_view.m_offsets = utils::make_array_view<ConnectivityType>(n_offsets);
    }

    // subelements (face definitions)
    if(!m_view.m_makePointMesh)
    {
      n_topology["subelements/shape"] = "polygonal";
      conduit::Node &n_se_conn = n_topology["subelements/connectivity"];
      n_se_conn.set_allocator(c2a.getConduitAllocatorID());
      n_se_conn.set(conduit::DataType(
        utils::cpp2conduit<ConnectivityType>::id,
        numFragments * m_view.m_MAX_FACES_PER_FRAGMENT * m_view.m_MAX_POINTS_PER_FACE));
      m_view.m_subelement_connectivity = utils::make_array_view<ConnectivityType>(n_se_conn);
      m_view.m_subelement_connectivity.fill(UnusedValue);

      conduit::Node &n_se_sizes = n_topology["subelements/sizes"];
      n_se_sizes.set_allocator(c2a.getConduitAllocatorID());
      n_se_sizes.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id,
                                       numFragments * m_view.m_MAX_FACES_PER_FRAGMENT));
      m_view.m_subelement_sizes = utils::make_array_view<ConnectivityType>(n_se_sizes);
      m_view.m_subelement_sizes.fill(ConnectivityType {0});

      conduit::Node &n_se_offsets = n_topology["subelements/offsets"];
      n_se_offsets.set_allocator(c2a.getConduitAllocatorID());
      n_se_offsets.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id,
                                         numFragments * m_view.m_MAX_FACES_PER_FRAGMENT));
      m_view.m_subelement_offsets = utils::make_array_view<ConnectivityType>(n_se_offsets);
      m_view.m_subelement_offsets.fill(UnusedValue);
    }

    // Make new fields.
    conduit::Node &n_origElem = n_fields[originalElementsField];
    n_origElem["topology"] = n_topology.name();
    n_origElem["association"] = "element";
    conduit::Node &n_orig_elem_values = n_origElem["values"];
    n_orig_elem_values.set_allocator(c2a.getConduitAllocatorID());
    n_orig_elem_values.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, numFragments));
    m_view.m_original_zones = utils::make_array_view<ConnectivityType>(n_orig_elem_values);

    if(m_view.m_makePlane)
    {
      conduit::Node &n_normal = n_fields["normal"];
      n_normal["topology"] = n_topology.name();
      n_normal["association"] = "element";
      conduit::Node &n_x = n_normal["values/x"];
      conduit::Node &n_y = n_normal["values/y"];
      conduit::Node &n_z = n_normal["values/z"];
      n_x.set_allocator(c2a.getConduitAllocatorID());
      n_x.set(conduit::DataType(utils::cpp2conduit<double>::id, numFragments));
      m_view.m_norm_x = utils::make_array_view<double>(n_x);
      n_y.set_allocator(c2a.getConduitAllocatorID());
      n_y.set(conduit::DataType(utils::cpp2conduit<double>::id, numFragments));
      m_view.m_norm_y = utils::make_array_view<double>(n_y);
      n_z.set_allocator(c2a.getConduitAllocatorID());
      n_z.set(conduit::DataType(utils::cpp2conduit<double>::id, numFragments));
      m_view.m_norm_z = utils::make_array_view<double>(n_z);

      conduit::Node &n_planeOffset = n_fields["offset"];
      n_planeOffset["topology"] = n_topology.name();
      n_planeOffset["association"] = "element";
      conduit::Node &n_values = n_planeOffset["values"];
      n_values.set_allocator(c2a.getConduitAllocatorID());
      n_values.set(conduit::DataType(utils::cpp2conduit<double>::id, numFragments));
      m_view.m_plane_offset = utils::make_array_view<double>(n_values);
    }

    // Set up new matset. All of the sizes are numFragments because we're making clean zones.
    n_matset["topology"] = n_topology.name();
    conduit::Node &n_volume_fractions = n_matset["volume_fractions"];
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_volume_fractions.set(conduit::DataType(utils::cpp2conduit<MaterialVF>::id, numFragments));
    m_view.m_volume_fractions = utils::make_array_view<MaterialVF>(n_volume_fractions);

    conduit::Node &n_material_ids = n_matset["material_ids"];
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_material_ids.set(conduit::DataType(utils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_material_ids = utils::make_array_view<MaterialID>(n_material_ids);

    conduit::Node &n_indices = n_matset["indices"];
    n_indices.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set(conduit::DataType(utils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_indices = utils::make_array_view<MaterialID>(n_indices);

    conduit::Node &n_mat_sizes = n_matset["sizes"];
    n_mat_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_mat_sizes.set(conduit::DataType(utils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_sizes = utils::make_array_view<MaterialID>(n_mat_sizes);
    m_view.m_mat_sizes.fill(MaterialID(0));

    conduit::Node &n_mat_offsets = n_matset["offsets"];
    n_mat_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_mat_offsets.set(conduit::DataType(utils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_offsets = utils::make_array_view<MaterialID>(n_mat_offsets);
  }

  /*!
   * \brief A view that can be used on device to encode a shape into Blueprint arrays.
   */
  struct View
  {
    /*!
     * \brief Encode a polyhedron from MIR into the views that represent the Conduit data.
     *
     * \param zoneIndex The original zone that yielded the fragment.
     * \param fragmentOffset The offset into the per-fragment data.
     * \param shape The 3D primal polyhedron we're encoding.
     * \param matId The material id for the fragment.
     * \param planeOrigin The origin point for the fragment's cutting plane.
     * \param planeOffset The offset for the fragment's cutting plane.
     * \param planeNormal The normal for the fragment's cutting plane.
     */
    AXOM_HOST_DEVICE
    void addShape(axom::IndexType zoneIndex,
                  axom::IndexType fragmentOffset,
                  const PHShape &shape,
                  int matId,
                  const PointType &planeOrigin,
                  double planeOffset,
                  const double *planeNormal) const
    {
      const int nverts = shape.numVertices();
      SLIC_ASSERT(nverts <= m_MAX_POINTS_PER_FRAGMENT);

      if(m_makePointMesh)
      {
        m_x[fragmentOffset] = planeOrigin[0];
        m_y[fragmentOffset] = planeOrigin[1];
        m_z[fragmentOffset] = planeOrigin[2];

        m_connectivity[fragmentOffset] = fragmentOffset;
        m_sizes[fragmentOffset] = 1;
        m_offsets[fragmentOffset] = fragmentOffset;
      }
      else
      {
        // Copy coordinates into coordinate arrays. We might end up with unreferenced
        // coordinates for now.
        const auto coordOffset = fragmentOffset * m_MAX_POINTS_PER_FRAGMENT;
        for(int i = 0; i < nverts; i++)
        {
          const auto destIndex = coordOffset + i;
          const auto &pt = shape[i];
          m_x[destIndex] = pt[0];
          m_y[destIndex] = pt[1];
          m_z[destIndex] = pt[2];
        }

        // Get pointers to where this shape's faces should be stored in the subelement data.
        const auto faceOffset = fragmentOffset * m_MAX_FACES_PER_FRAGMENT;
        ConnectivityType *subelement_connectivity = m_subelement_connectivity.data() +
          fragmentOffset * (m_MAX_FACES_PER_FRAGMENT * m_MAX_POINTS_PER_FACE);
        ConnectivityType *subelement_sizes = m_subelement_sizes.data() + faceOffset;
        ConnectivityType *subelement_offsets = m_subelement_offsets.data() + faceOffset;
        axom::IndexType numFaces;

        // Get the faces from the actual shape, directly into the subelement data
        // that defines the faces.
        shape.getFaces(subelement_connectivity, subelement_sizes, subelement_offsets, numFaces);

#if defined(AXOM_ELVIRA_DEBUG_MAKE_FRAGMENTS) && !defined(AXOM_DEVICE_CODE)
        std::stringstream ss;
        ss << "addShape: zoneIndex=" << zoneIndex << ", fragmentOffset=" << fragmentOffset
           << ", nverts=" << nverts << ", numFaces=" << numFaces << ", subelement_connectivity={";
        for(int i = 0; i < m_MAX_FACES_PER_FRAGMENT * m_MAX_POINTS_PER_FACE; i++)
        {
          ss << subelement_connectivity[i] << ", ";
        }
        ss << "}, subelement_sizes={";
        for(int i = 0; i < m_MAX_FACES_PER_FRAGMENT; i++)
        {
          ss << subelement_sizes[i] << ", ";
        }
        ss << "}, subelement_offsets={";
        for(int i = 0; i < m_MAX_FACES_PER_FRAGMENT; i++)
        {
          ss << subelement_offsets[i] << ", ";
        }
        ss << "}";
        SLIC_DEBUG(ss.str());
#endif
        // Make the face offsets relative to the whole subelement array.
        axom::IndexType index = 0;
        for(axom::IndexType f = 0; f < numFaces; f++)
        {
          subelement_offsets[f] += fragmentOffset * m_MAX_POINTS_PER_FACE * m_MAX_FACES_PER_FRAGMENT;

#if !defined(AXOM_DEVICE_CODE)
          SLIC_ASSERT_MSG(
            subelement_sizes[f] <= m_MAX_POINTS_PER_FACE,
            axom::fmt::format(
              "Zone {} has {} points in face {} but should have no more than {} points. shape={}",
              zoneIndex,
              subelement_sizes[f],
              f,
              m_MAX_POINTS_PER_FACE,
              shape));
#endif

          for(ConnectivityType i = 0; i < subelement_sizes[f]; i++)
          {
            subelement_connectivity[index++] += fragmentOffset * m_MAX_POINTS_PER_FRAGMENT;
          }
        }

        // Make connectivity.
        auto connOffset = fragmentOffset * m_MAX_FACES_PER_FRAGMENT;
        for(axom::IndexType i = 0; i < numFaces; i++)
        {
          const auto idx = static_cast<ConnectivityType>(connOffset + i);
          m_connectivity[idx] = idx;
        }
        m_sizes[fragmentOffset] = static_cast<ConnectivityType>(numFaces);
        m_offsets[fragmentOffset] = connOffset;
      }

      // Save fields.
      m_original_zones[fragmentOffset] = zoneIndex;
      if(m_makePlane)
      {
        m_norm_x[fragmentOffset] = planeNormal[0];
        m_norm_y[fragmentOffset] = planeNormal[1];
        m_norm_z[fragmentOffset] = planeNormal[2];

        m_plane_offset[fragmentOffset] = planeOffset;
      }

      // Save material data.
      m_volume_fractions[fragmentOffset] = static_cast<MaterialVF>(1.);
      m_material_ids[fragmentOffset] = static_cast<MaterialID>(matId);
      m_mat_sizes[fragmentOffset] = static_cast<MaterialID>(1);
      m_mat_offsets[fragmentOffset] = fragmentOffset;
      m_mat_indices[fragmentOffset] = fragmentOffset;
    }

    // These ArrayView objects expose data that we've allocated in Conduit nodes to represent the new mesh.

    // New coordset data.
    axom::ArrayView<CoordType> m_x {}, m_y {}, m_z {};  //!< X,Y,Z coordinates

    // New topology data.
    axom::ArrayView<ConnectivityType> m_connectivity {}, m_sizes {}, m_offsets {};

    axom::ArrayView<ConnectivityType> m_subelement_connectivity {}, m_subelement_sizes {},
      m_subelement_offsets {};

    // New field data.
    axom::ArrayView<axom::IndexType> m_original_zones {};  //!< View for originalZone field data.
    axom::ArrayView<double> m_norm_x {}, m_norm_y {}, m_norm_z {};  //!< Fragment normals
    axom::ArrayView<double> m_plane_offset;                         //!< Fragment plane offset

    // New matset data
    axom::ArrayView<MaterialVF> m_volume_fractions {};
    axom::ArrayView<MaterialID> m_material_ids {}, m_mat_indices, m_mat_sizes {}, m_mat_offsets {};

    bool m_makePlane {false};
    bool m_makePointMesh {false};

    // Fragment sizing information
    axom::IndexType m_MAX_POINTS_PER_FACE {0};
    axom::IndexType m_MAX_FACES_PER_FRAGMENT {0};
    axom::IndexType m_MAX_POINTS_PER_FRAGMENT {0};
  };

  /*!
   * \brief Return a view that can be copied to device.
   *
   * \return A view that be used to store data.
   */
  View view() { return m_view; }

  /*!
   * \brief Clean the mesh, merging coordinates and faces.
   *
   * \param n_coordset The coordset to clean up.
   * \param point_tolerance The point tolerance used to merge points.
   * \param n_topology The topology to clean up.
   * \param[out] selectedIds An array that indicates the points that were selected during coordset point merging.
   *
   * \note This method invalidates the views in m_view by causing some of their backing arrays to be replaced.
   */
  void cleanMesh(conduit::Node &n_coordset,
                 double point_tolerance,
                 conduit::Node &n_topology,
                 axom::Array<axom::IndexType> &selectedIds) const
  {
    AXOM_ANNOTATE_SCOPE("cleanMesh");

    // If we're making a point mesh then there is nothing to do.
    if(m_view.m_makePointMesh)
    {
      return;
    }

    // _bump_utilities_mergecoordsetpoints_begin
    namespace utils = axom::bump::utilities;
    axom::Array<axom::IndexType> old2new;

    auto newCoordsetView =
      axom::bump::views::make_explicit_coordset<CoordType, CoordsetView::dimension()>::view(
        n_coordset);
    using NewCoordsetView = decltype(newCoordsetView);
    axom::bump::MergeCoordsetPoints<ExecSpace, NewCoordsetView> mcp(newCoordsetView);
    conduit::Node n_mcp_options;
    n_mcp_options["tolerance"] = point_tolerance;
    const bool merged = mcp.execute(n_coordset, n_mcp_options, selectedIds, old2new);
    // _bump_utilities_mergecoordsetpoints_end

    // Changing the coordset changed the nodes so we need to change the subelement/connectivity.
    // Traverse it using sizes/offsets in case some of the connectivity values are being skipped.
    if(merged)
    {
      AXOM_ANNOTATE_SCOPE("rewriting_subelements");
      conduit::Node &n_se_conn = n_topology["subelements/connectivity"];
      const conduit::Node &n_se_sizes = n_topology["subelements/sizes"];
      const conduit::Node &n_se_offsets = n_topology["subelements/offsets"];

      auto se_conn = utils::make_array_view<ConnectivityType>(n_se_conn);
      const auto se_sizes = utils::make_array_view<ConnectivityType>(n_se_sizes);
      const auto se_offsets = utils::make_array_view<ConnectivityType>(n_se_offsets);
      auto old2newView = old2new.view();
      axom::for_all<ExecSpace>(
        se_sizes.size(),
        AXOM_LAMBDA(axom::IndexType index) {
          const auto size = se_sizes[index];
          const auto offset = se_offsets[index];
          for(ConnectivityType i = 0; i < size; i++)
          {
            const auto nodeId = se_conn[offset + i];
            se_conn[offset + i] = old2newView[nodeId];
          }
        });
    }

    // Now merge any faces that can be merged.
    axom::bump::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(n_topology);
  }

private:
  View m_view;
};

//------------------------------------------------------------------------------
/*!
 * \brief Base template for letting us pick different behaviors for making clean zones.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam TopologyView The view type that describes the Blueprint input mesh.
 * \tparam CoordsetView The view type that describes the coordinates.
 * \tparam MatsetView The view type that describes matset.
 * \tparam NDIMS The dimension of the topology.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView, int NDIMS>
struct MakeCleanZones
{
  /*!
   * \brief Make a new mesh containing the selected clean zones.
   *
   * \param cleanZones The clean zone numbers to extract.
   * \param n_mesh The input mesh.
   * \param n_options The algorithm options.
   * \param topologyView The view that wraps the mesh topology.
   * \param coordsetView The view that wraps the coordinates.
   * \param matsetView The view that wraps the matset.
   * \param[out] n_cleanOutput The node that will contain the new mesh.
   */
  static void execute(const axom::ArrayView<axom::IndexType> &cleanZones,
                      const conduit::Node &n_mesh,
                      const conduit::Node &n_options,
                      const TopologyView &topologyView,
                      const CoordsetView &coordsetView,
                      const MatsetView &matsetView,
                      conduit::Node &n_cleanOutput)
  {
    // Make the clean mesh.
    axom::bump::ExtractZonesAndMatset<ExecSpace, TopologyView, CoordsetView, MatsetView> ez(
      topologyView,
      coordsetView,
      matsetView);
    ez.execute(cleanZones, n_mesh, n_options, n_cleanOutput);
  }
};

/*!
 * \brief Partial specialization for making clean zones for 3D zones. We do this to
 *        make better polyhedra upfront, rather than converting to polyhedral at
 *        greater expense later when we merge meshes.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam TopologyView The view type that describes the Blueprint input mesh.
 * \tparam CoordsetView The view type that describes the coordinates.
 * \tparam MatsetView The view type that describes matset.
 * \tparam NDIMS The dimension of the topology.
 */
template <typename ExecSpace, typename TopologyView, typename CoordsetView, typename MatsetView>
struct MakeCleanZones<ExecSpace, TopologyView, CoordsetView, MatsetView, 3>
{
  /*!
   * \brief Make a new mesh containing the selected clean zones.
   *
   * \param cleanZones The clean zone numbers to extract.
   * \param n_mesh The input mesh.
   * \param n_options The algorithm options.
   * \param topologyView The view that wraps the mesh topology.
   * \param coordsetView The view that wraps the coordinates.
   * \param matsetView The view that wraps the matset.
   * \param[out] n_cleanOutput The node that will contain the new mesh.
   */
  static void execute(const axom::ArrayView<axom::IndexType> &cleanZones,
                      const conduit::Node &n_mesh,
                      const conduit::Node &n_options,
                      const TopologyView &topologyView,
                      const CoordsetView &coordsetView,
                      const MatsetView &matsetView,
                      conduit::Node &n_cleanOutput)
  {
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    // Make the clean mesh. We directly make polyhedral output so we can skip
    // mesh improvement stages later when we merge the output with the mixed
    // zone polyhedra.
    axom::bump::ExtractZonesAndMatsetPolyhedral<ExecSpace, IndexingPolicy, CoordsetView, MatsetView> ez(
      topologyView,
      coordsetView,
      matsetView);
    ez.execute(cleanZones, n_mesh, n_options, n_cleanOutput);
  }
};

}  // end namespace detail
}  // end namespace mir
}  // end namespace axom

#endif
