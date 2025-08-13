// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_MERGE_POLYHEDRAL_FACES_HPP_
#define AXOM_BUMP_MERGE_POLYHEDRAL_FACES_HPP_

#include "axom/core.hpp"
#include "axom/bump/utilities/blueprint_utilities.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/HashNaming.hpp"
#include "axom/bump/Unique.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{

/*!
 * \brief Take a Blueprint polyhedral mesh and merge faces and rewrite the
 *        element and subelement connectivity.
 * 
 *        Faces are lists of nodes that define a polygonal face of a polyhedron
 *        and they are defined in the subelement connectivity. Blueprint's
 *        polyhedral topologies define zones using a list of face ids where each
 *        face id refers to a polygonal face in the subelement connectivity.
 *
 *        Faces are considered duplicates when they contain the same list of node
 *        ids, disregarding any difference in node ordering. When a face is found
 *        to have a duplicate, a single face definition is kept and adjacent
 *        polyhedral zones are updated to use the single shared face definition.
 *
 * \tparam ExecSpace The execution space where the algorithm will execute.
 * \tparam ConnectivityType The type used in the topology connectivity arrays.
 */
template <typename ExecSpace, typename ConnectivityType>
class MergePolyhedralFaces
{
public:
  /*!
   * \brief Merge faces using the face definitions and rewrite the element
   *        and subelement connectivity. The old definitions are replaced with
   *        new data.
   *
   * \param n_topology The topology to modify.
   */
  static void execute(conduit::Node &n_topo)
  {
    SLIC_ASSERT(n_topo.fetch_existing("elements/shape").as_string() == "polyhedral");
    namespace utils = axom::bump::utilities;

    AXOM_ANNOTATE_SCOPE("MergePolyhedralFaces");
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Get the data from the topology and make views.
    conduit::Node &n_elem_conn = n_topo["elements/connectivity"];
    conduit::Node &n_elem_sizes = n_topo["elements/sizes"];
    conduit::Node &n_elem_offsets = n_topo["elements/offsets"];
    conduit::Node &n_se_conn = n_topo["subelements/connectivity"];
    conduit::Node &n_se_sizes = n_topo["subelements/sizes"];
    conduit::Node &n_se_offsets = n_topo["subelements/offsets"];
    auto elem_conn = utils::make_array_view<ConnectivityType>(n_elem_conn);
    auto elem_sizes = utils::make_array_view<ConnectivityType>(n_elem_sizes);
    auto elem_offsets = utils::make_array_view<ConnectivityType>(n_elem_offsets);
    const auto se_conn = utils::make_array_view<ConnectivityType>(n_se_conn);
    const auto se_sizes = utils::make_array_view<ConnectivityType>(n_se_sizes);
    const auto se_offsets = utils::make_array_view<ConnectivityType>(n_se_offsets);

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("maxnode");
    // NOTE: the se_conn array may have gap values (typically -1).
    axom::ReduceMax<ExecSpace, ConnectivityType> reduceMaxNodeId(0);
    axom::for_all<ExecSpace>(
      se_conn.size(),
      AXOM_LAMBDA(axom::IndexType index) { reduceMaxNodeId.max(se_conn[index]); });
    const auto maxNodeId = reduceMaxNodeId.get();
    AXOM_ANNOTATE_END("maxnode");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("naming");
    constexpr int MaxPointsPerFace = 16;
    using NamingType = HashNaming<ConnectivityType, MaxPointsPerFace>;
    using KeyType = typename NamingType::KeyType;

    NamingType naming;
    naming.setMaxId(maxNodeId);
    auto namingView = naming.view();

    const axom::IndexType totalFaces = se_sizes.size();
    axom::Array<KeyType> faceNames(totalFaces, totalFaces, allocatorID);
    auto faceNamesView = faceNames.view();
    axom::for_all<ExecSpace>(
      totalFaces,
      AXOM_LAMBDA(axom::IndexType faceIndex) {
        // Get size for current face.
        const auto faceSize = static_cast<int>(se_sizes[faceIndex]);
        // RelWithDebInfo workaround - "16" substitutes lambda capture device failure for "MaxPointsPerFace"
        SLIC_ASSERT(faceSize <= 16);

        KeyType faceName {};
        if(faceSize > 0)
        {
          // Get offset for current face.
          const auto faceOffset = static_cast<axom::IndexType>(se_offsets[faceIndex]);
          SLIC_ASSERT(faceOffset < se_conn.size());

          // Get the ids for the current face.
          const auto faceIds = se_conn.data() + faceOffset;

          // Make a name for the current face.
          faceName = namingView.makeName(faceIds, faceSize);
        }

        faceNamesView[faceIndex] = faceName;
      });
    AXOM_ANNOTATE_END("naming");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("unique");
    // The unique keys.
    axom::Array<KeyType> uniqueKeys;
    // The index of the sorted key in the original keys. We can use this to
    // determine which face to use.
    axom::Array<axom::IndexType> selectedFaces;

    // Make faces unique.
    axom::bump::Unique<ExecSpace, KeyType>::execute(faceNamesView, uniqueKeys, selectedFaces);
    const auto uniqueKeysView = uniqueKeys.view();
    const auto selectedFacesView = selectedFaces.view();
    AXOM_ANNOTATE_END("unique");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("rewriting_subelements");
    conduit::Node n_new_se_sizes;
    n_new_se_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_new_se_sizes.set(
      conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, selectedFaces.size()));
    auto new_se_sizes = utils::make_array_view<ConnectivityType>(n_new_se_sizes);

    conduit::Node n_new_se_offsets;
    n_new_se_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_new_se_offsets.set(
      conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, selectedFaces.size()));
    auto new_se_offsets = utils::make_array_view<ConnectivityType>(n_new_se_offsets);

    // Copy the sizes of the selected faces into new_se_sizes and make new_se_offsets
    axom::ReduceSum<ExecSpace, axom::IndexType> reduceNewSizes(0);
    axom::for_all<ExecSpace>(
      selectedFaces.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        const auto size = se_sizes[selectedFacesView[index]];
        new_se_sizes[index] = size;
        reduceNewSizes += size;
      });
    const axom::IndexType newSEConnSize = reduceNewSizes.get();
    axom::exclusive_scan<ExecSpace>(new_se_sizes, new_se_offsets);

    // Allocate new_se_conn to contain the new face definitions.
    conduit::Node n_new_se_conn;
    n_new_se_conn.set_allocator(c2a.getConduitAllocatorID());
    n_new_se_conn.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, newSEConnSize));
    auto new_se_conn = utils::make_array_view<ConnectivityType>(n_new_se_conn);

    // Copy the selected faces into new_se_conn.
    axom::for_all<ExecSpace>(
      selectedFaces.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        const auto numFaceIds = new_se_sizes[index];
        const auto destOffset = new_se_offsets[index];
        const auto srcOffset = se_offsets[selectedFacesView[index]];
        for(ConnectivityType i = 0; i < numFaceIds; i++)
        {
          new_se_conn[destOffset + i] = se_conn[srcOffset + i];
        }
      });

    // Move the "new" nodes into the Blueprint hierarchy.
    n_se_conn.move(n_new_se_conn);
    n_se_sizes.move(n_new_se_sizes);
    n_se_offsets.move(n_new_se_offsets);
    AXOM_ANNOTATE_END("rewriting_subelements");

    //--------------------------------------------------------------------------
    // Now, rewrite the connectivity with the new face ids.
    AXOM_ANNOTATE_BEGIN("rewriting_elements");

    // Sum the element sizes so we can check for gaps to eliminate.
    axom::ReduceSum<ExecSpace, axom::IndexType> reduceConnSize(0);
    axom::for_all<ExecSpace>(
      elem_sizes.size(),
      AXOM_LAMBDA(axom::IndexType index) { reduceConnSize += elem_sizes[index]; });
    const auto totalConnSize = reduceConnSize.get();

    if(totalConnSize == elem_conn.size())
    {
      // The connectivity has no gaps.

      axom::for_all<ExecSpace>(
        elem_conn.size(),
        AXOM_LAMBDA(axom::IndexType index) {
          const auto originalFaceId = elem_conn[index];

          // Get the "name" of the old face.
          const auto originalFaceKey = faceNamesView[originalFaceId];

          // Look for the index of the "name" in the new uniqueKeys.
          // That will be its face index in the new faces.
          const auto newId = axom::utilities::binary_search(uniqueKeysView, originalFaceKey);
          SLIC_ASSERT(newId != -1);
          elem_conn[index] = static_cast<ConnectivityType>(newId);
        });
    }
    else
    {
      // The connectivity size does not match the expected sizes. It probably
      // contains gaps. We can rewrite it and the offsets.

      conduit::Node n_new_elem_conn;
      n_new_elem_conn.set_allocator(c2a.getConduitAllocatorID());
      n_new_elem_conn.set(conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, totalConnSize));
      auto new_elem_conn = utils::make_array_view<ConnectivityType>(n_new_elem_conn);

      conduit::Node n_new_elem_offsets;
      n_new_elem_offsets.set_allocator(c2a.getConduitAllocatorID());
      n_new_elem_offsets.set(
        conduit::DataType(utils::cpp2conduit<ConnectivityType>::id, elem_sizes.size()));
      auto new_elem_offsets = utils::make_array_view<ConnectivityType>(n_new_elem_offsets);
      axom::exclusive_scan<ExecSpace>(elem_sizes, new_elem_offsets);

      axom::for_all<ExecSpace>(
        elem_sizes.size(),
        AXOM_LAMBDA(axom::IndexType index) {
          const auto srcOffset = elem_offsets[index];
          const auto destOffset = new_elem_offsets[index];
          const auto size = elem_sizes[index];
          for(ConnectivityType i = 0; i < size; i++)
          {
            const auto originalFaceId = elem_conn[srcOffset + i];

            // Get the "name" of the old face.
            const auto originalFaceKey = faceNamesView[originalFaceId];

            // Look for the index of the "name" in the new uniqueKeys.
            // That will be its face index in the new faces.
            const auto newId = axom::utilities::binary_search(uniqueKeysView, originalFaceKey);
            SLIC_ASSERT(newId != -1);
            new_elem_conn[destOffset + i] = static_cast<ConnectivityType>(newId);
          }
        });

      n_elem_conn.move(n_new_elem_conn);
      n_elem_offsets.move(n_new_elem_offsets);
    }

    AXOM_ANNOTATE_END("rewriting_elements");
  }
};

}  // end namespace bump
}  // end namespace axom

#endif
