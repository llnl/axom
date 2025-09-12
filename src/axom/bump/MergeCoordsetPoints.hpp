// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_MERGE_COORDSET_POINTS_HPP_
#define AXOM_BUMP_MERGE_COORDSET_POINTS_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/bump/CoordsetSlicer.hpp"
#include "axom/bump/CoordsetExtents.hpp"

#include <conduit/conduit.hpp>

// Uncomment to write files for the coordset inputs/outputs.
//#define AXOM_DEBUG_MERGE_COORDSET_POINTS

#if defined(AXOM_DEBUG_MERGE_COORDSET_POINTS)
  #include <conduit/conduit_relay_io_blueprint.hpp>
#endif

// Includes for appropriate round function.
#if defined(AXOM_DEVICE_CODE)
  #if defined(AXOM_USE_CUDA)
    #include <cuda_runtime.h>
  #endif
  #if defined(AXOM_USE_HIP)
    #include <hip/hip_runtime.h>
  #endif
#else
  #include <cmath>
#endif

namespace axom
{
namespace bump
{
namespace detail
{

/*!
 * \brief Round a number.
 */
template <typename value_type>
struct Rounder
{
  static AXOM_HOST_DEVICE inline value_type execute(value_type value)
  {
#if defined(AXOM_DEVICE_CODE)
    return round(value);
#else
    return std::round(value);
#endif
  }
};

/*!
 * \brief Round a float number.
 */
template <>
struct Rounder<float>
{
  static AXOM_HOST_DEVICE inline float execute(float value)
  {
#if defined(AXOM_DEVICE_CODE)
    return roundf(value);
#else
    return std::round(value);
#endif
  }
};

}  // end namespace detail

/**
 * \brief Merge coordset points that are within a tolerance into the same point
 *        and provide information on how it was done.
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 * \tparam CoordsetView The coordset view type.
 */
template <typename ExecSpace, typename CoordsetView>
class MergeCoordsetPoints
{
public:
  using value_type = typename CoordsetView::value_type;

  /*!
   * \brief Constructor
   *
   * \param coordsetView The coordset view that wraps the coordset to be modified.
   */
  MergeCoordsetPoints(const CoordsetView &coordsetView) : m_coordsetView(coordsetView) { }

  /*!
   * \brief Merge the coordset points using a tolerance and pass out an array of the
   *        points that were selected from the original coordset since they can be
   *        used to slice nodal fields.
   *
   * \param[inout] n_coordset The Conduit node that contains the coordset. The coordset data
   *                          will be replaced.
   * \param n_options A Conduit node that contains the supported options. At present, there
   *                  is just a "tolerance" option that specifies how close points need to
   *                  be in order to be matched.
   * \param[out] selectedIds An array containing the ids of the nodes that are selected from
   *                         the old nodes.
   * \param[out] old2new An array with nnodes elements that contains the new node id for each
   *                     node in the old coordset. This can be useful in rewriting connectivity
   *                     or fields.
   *
   * \note The algorithm eliminates a little precision on the input coordinates to make a
   *       hashed name for the point. Points that are close enough should have the same
   *       hashed name. The names are passed through a Unique filter to select a set of
   *       unique names. Multiple points could result in the same hash name so one of the
   *       points is selected as the new definition for the unique point (the selection
   *       arises from sorting the hash names). This means that among a bunch of points
   *       that are close enough to be grouped, one of them will be picked and it might
   *       not be the best one "esthetically". The points are output in sorted hash name
   *       order too, which is a bit random.
   *
   * \return True if point merging happened; False if no point merging was needed.
   */
  bool execute(conduit::Node &n_coordset,
               const conduit::Node &n_options,
               axom::Array<axom::IndexType> &selectedIds,
               axom::Array<axom::IndexType> &old2new) const
  {
    namespace utils = axom::bump::utilities;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // If the coordset is not explicit then there is nothing to do.
    if(n_coordset["type"].as_string() != "explicit")
    {
      SLIC_ERROR("The input coordset was not explicit.");
      return false;
    }

    // Adjust the name of the timer so we can see the coordinate type.
    AXOM_ANNOTATE_SCOPE(
      axom::fmt::format("MergeCoordsetPoints<{}>", utils::cpp2conduit<value_type>::name));

    // Get any options.
    constexpr double DEFAULT_TOLERANCE = std::is_same<value_type, float>::value
      ? (axom::numeric_limits<float>::epsilon() * 4.f)
      : 1.e-10;
    double tolerance = DEFAULT_TOLERANCE;
    if(n_options.has_child("tolerance"))
    {
      tolerance = axom::utilities::abs(n_options["tolerance"].to_double());
      SLIC_ASSERT(tolerance > 0.);
    }

    //--------------------------------------------------------------------------
    // Create names for the coordinates.
    using KeyType = std::uint64_t;
    axom::Array<KeyType> coordNames;
    createNames<KeyType>(coordNames, tolerance);
    auto coordNamesView = coordNames.view();

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("unique");

    // Make points unique.
    axom::Array<KeyType> uniqueNames;
    axom::bump::Unique<ExecSpace, KeyType>::execute(coordNamesView, uniqueNames, selectedIds);
    const auto uniqueNamesView = uniqueNames.view();
    const auto selectedIdsView = selectedIds.view();
    AXOM_ANNOTATE_END("unique");

    //--------------------------------------------------------------------------
    const auto nnodes = m_coordsetView.numberOfNodes();
    const bool merged = selectedIds.size() < nnodes;
    old2new = axom::Array<axom::IndexType>(nnodes, nnodes, allocatorID);
    auto old2newView = old2new.view();
    if(merged)
    {
      // There are fewer nodes in the selectedIds so we are able to combine nodes.
      SLIC_INFO(axom::fmt::format("Merged {} nodes into {} nodes.", nnodes, selectedIds.size()));

      AXOM_ANNOTATE_BEGIN("old2new");
      // Make a map of nodes in the old coordset to nodes in the new coordset. We
      // do it by looking up the old node name in the new coordset unique names.
      axom::for_all<ExecSpace>(
        nnodes,
        AXOM_LAMBDA(axom::IndexType index) {
          const auto newNodeId =
            axom::utilities::binary_search(uniqueNamesView, coordNamesView[index]);
          // RelWithDebInfo workaround - "old2newView.size() substitutes lambda capture device failure for "nnodes"
          SLIC_ASSERT(newNodeId >= 0 && newNodeId < old2newView.size());
          old2newView[index] = newNodeId;
        });
      AXOM_ANNOTATE_END("old2new");

      //--------------------------------------------------------------------------
      // Use the selectedIds to slice the coordset to make a new coordset that
      // replaces the old one.
      CoordsetSlicer<ExecSpace, CoordsetView> css(m_coordsetView);
      SliceData slice;
      slice.m_indicesView = selectedIdsView;
      conduit::Node n_sliced;
      css.execute(slice, n_coordset, n_sliced);

#if defined(AXOM_DEBUG_MERGE_COORDSET_POINTS)
      // Add a mesh for the input point mesh.
      conduit::Node n_mesh;
      n_mesh["coordsets"][n_coordset.name()].set_external(n_coordset);
      n_mesh["topologies/mesh/type"] = "unstructured";
      n_mesh["topologies/mesh/coordset"] = n_coordset.name();
      n_mesh["topologies/mesh/elements/shape"] = "point";
      std::vector<int> conn(coordNamesView.size());
      std::iota(conn.begin(), conn.end(), 0);
      n_mesh["topologies/mesh/elements/connectivity"].set(conn);
      n_mesh["topologies/mesh/type"] = "unstructured";
      n_mesh["fields/coordNames/association"] = "vertex";
      n_mesh["fields/coordNames/topology"] = "mesh";
      n_mesh["fields/coordNames/values"].set_external(coordNamesView.data(), coordNamesView.size());

      // Add a mesh for the output point mesh.
      n_mesh["coordsets/output"].set_external(n_sliced);
      n_mesh["topologies/merged/type"] = "unstructured";
      n_mesh["topologies/merged/coordset"] = "output";
      n_mesh["topologies/merged/elements/shape"] = "point";
      std::vector<int> merged_conn(uniqueNamesView.size());
      std::iota(merged_conn.begin(), merged_conn.end(), 0);
      n_mesh["topologies/merged/elements/connectivity"].set(merged_conn);
      n_mesh["topologies/merged/type"] = "unstructured";
      n_mesh["fields/mergedCoordNames/association"] = "vertex";
      n_mesh["fields/mergedCoordNames/topology"] = "merged";
      n_mesh["fields/mergedCoordNames/values"].set_external(uniqueNamesView.data(),
                                                            uniqueNamesView.size());

      // Write files to examine the input and output coordsets.
      conduit::Node n_host;
      utils::copy<SEQ_EXEC>(n_host, n_mesh);
      std::string name(axom::execution_space<ExecSpace>::name());
      name = name.substr(1, 3);
      conduit::relay::io::blueprint::save_mesh(n_host, "merge_coordset_points_" + name, "hdf5");
      conduit::relay::io::save(n_host, "merge_coordset_points_" + name + ".yaml", "yaml");
#endif

      // Move the sliced coordset to the input coordset.
      n_coordset.move(n_sliced);
    }
    else
    {
      // It does not look like there is a need for point merging.
      AXOM_ANNOTATE_BEGIN("old2new");

      auto selectedIdsView = selectedIds.view();
      axom::for_all<ExecSpace>(
        nnodes,
        AXOM_LAMBDA(axom::IndexType index) {
          selectedIdsView[index] = index;
          old2newView[index] = index;
        });
    }

    return merged;
  }

  /*!
   * \brief Create names for the coordinates so we can merge them.
   *
   * \param[out] coordNames The "names" of the coordinates that will be used when merging.
   * \param tolerance The tolerance used to merge points.
   */
  template <typename KeyType>
  void createNames(axom::Array<KeyType> &coordNames, double tolerance) const
  {
    constexpr double smallTolerance = 1.e-6;

    // When tolerances that are large enough, and the data are already float,
    // try making names using float.
    if(tolerance >= smallTolerance && std::is_same<value_type, float>::value)
    {
      createNamesInner<KeyType, float>(coordNames, tolerance);
    }
    else
    {
      // The tolerance is small enough to require more precision, or we have
      // double coordinates.
      createNamesInner<KeyType, double>(coordNames, tolerance);
    }
  }

  /*!
   * \brief Create names for the coordinates so we can merge them.
   *
   * \param[out] coordNames The "names" of the coordinates that will be used when merging.
   * \param tolerance The tolerance used to merge points.
   */
  template <typename KeyType, typename Precision>
  void createNamesInner(axom::Array<KeyType> &coordNames, double tolerance) const
  {
    namespace utils = axom::bump::utilities;
    AXOM_ANNOTATE_SCOPE(axom::fmt::format("createNames<{}>", utils::cpp2conduit<Precision>::name));

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    const auto nnodes = m_coordsetView.numberOfNodes();
    const double neg_tolerance = -tolerance;

    coordNames = axom::Array<KeyType>(nnodes, nnodes, allocatorID);
    auto coordNamesView = coordNames.view();
    const auto deviceCoordsetView = m_coordsetView;
    axom::for_all<ExecSpace>(
      nnodes,
      AXOM_LAMBDA(axom::IndexType index) {
        // Get the current point.
        const auto pt = deviceCoordsetView[index];

        // Truncate the point components using the tolerance so we can eliminate
        // precision beyond what we want with the tolerance. The idea is that points
        // that are close enough will hash to the same name.
        Precision truncated[CoordsetView::dimension()];
        for(int d = 0; d < CoordsetView::dimension(); d++)
        {
          const Precision pointValue = static_cast<Precision>(pt[d]);
          Precision value = detail::Rounder<Precision>::execute(pointValue / tolerance) * tolerance;
          if(value > neg_tolerance && value < tolerance)
          {
            value = Precision {0};
          }
          truncated[d] = value;
        }

        // Make a name for this point
        const void *tptr = static_cast<const void *>(truncated);
        coordNamesView[index] =
          axom::utilities::hash_bytes(static_cast<const std::uint8_t *>(tptr),
                                      sizeof(Precision) * CoordsetView::dimension());
      });
  }

  CoordsetView m_coordsetView;
};

}  // end namespace bump
}  // end namespace axom

#endif
