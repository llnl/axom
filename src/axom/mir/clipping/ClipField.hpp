// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internal.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_CLIP_FIELD_HPP_
#define AXOM_MIR_CLIP_FIELD_HPP_

#include "axom/core.hpp"
#include "axom/mir/clipping/BlendGroupBuilder.hpp"
#include "axom/mir/clipping/ClipCases.h"
#include "axom/mir/clipping/ClipOptions.hpp"
#include "axom/mir/clipping/ClipTableManager.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"
#include "axom/mir/utilities/CoordsetBlender.hpp"
#include "axom/mir/utilities/FieldBlender.hpp"
#include "axom/mir/utilities/FieldSlicer.hpp"
#include "axom/mir/utilities/SelectedZones.hpp"
#include "axom/mir/utilities/utilities.hpp"
#include "axom/mir/views/Shapes.hpp"
#include "axom/mir/views/view_traits.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>
#include <conduit/conduit_relay_io.hpp>
#include <conduit/conduit_relay_io_blueprint.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#include <map>
#include <string>

// Enable code to save some debugging output files.
// #define AXOM_DEBUG_CLIP_FIELD

// Filter out degenerate zones in 2D.
#define AXOM_CLIP_FILTER_DEGENERATES

// Enable code to reduce the number of blend groups by NOT emitting 1-node
// blend groups and instead using a node lookup field for those.
#define AXOM_REDUCE_BLEND_GROUPS

namespace axom
{
namespace mir
{
namespace clipping
{
namespace internal
{
/*!
 * \brief Given an "ST_index" (e.g. ST_TET from clipping definitions), return an appropriate ShapeID value.
 *
 * \param st_index The value we want to translate into a ShapeID value.
 *
 * \return The ShapeID value that matches the st_index, or 0 if there is no match.
 */
template <typename IntegerType>
inline AXOM_HOST_DEVICE IntegerType ST_Index_to_ShapeID(IntegerType st_index)
{
  IntegerType shapeID = 0;
  switch(st_index)
  {
  case ST_LIN:
    shapeID = views::Line_ShapeID;
    break;
  case ST_TRI:
    shapeID = views::Tri_ShapeID;
    break;
  case ST_QUA:
    shapeID = views::Quad_ShapeID;
    break;
  case ST_TET:
    shapeID = views::Tet_ShapeID;
    break;
  case ST_PYR:
    shapeID = views::Pyramid_ShapeID;
    break;
  case ST_WDG:
    shapeID = views::Wedge_ShapeID;
    break;
  case ST_HEX:
    shapeID = views::Hex_ShapeID;
    break;
  }
  return shapeID;
}

/*!
 * \brief Returns a clip table index for the input shapeId.
 * \param shapeId A shapeID (e.g. Tet_ShapeID)
 * \return The clip table index for the shape.
 */
AXOM_HOST_DEVICE
inline int getClipTableIndex(int shapeId)
{
  int index = 0;
  switch(shapeId)
  {
  case views::Tri_ShapeID:
    index = 0;
    break;
  case views::Quad_ShapeID:
    index = 1;
    break;
  case views::Tet_ShapeID:
    index = 2;
    break;
  case views::Pyramid_ShapeID:
    index = 3;
    break;
  case views::Wedge_ShapeID:
    index = 4;
    break;
  case views::Hex_ShapeID:
    index = 5;
    break;
  }
  return index;
}

AXOM_HOST_DEVICE
inline bool color0Selected(int selection) { return axom::utilities::bitIsSet(selection, 0); }

AXOM_HOST_DEVICE
inline bool color1Selected(int selection) { return axom::utilities::bitIsSet(selection, 1); }

AXOM_HOST_DEVICE
inline bool generatedPointIsSelected(unsigned char color, int selection)
{
  return color == NOCOLOR || (color0Selected(selection) && color == COLOR0) ||
    (color1Selected(selection) && color == COLOR1);
}

AXOM_HOST_DEVICE
inline bool shapeIsSelected(unsigned char color, int selection)
{
  return (color0Selected(selection) && color == COLOR0) ||
    (color1Selected(selection) && color == COLOR1);
}

template <typename IdType, int MAXSIZE>
inline AXOM_HOST_DEVICE int unique_count(const IdType *values, int n)
{
  IdType v[MAXSIZE];
  // Start off with one unique element
  int nv = 1;
  v[0] = values[0];
  // Scan the rest
  for(int j = 1; j < n; j++)
  {
    int fi = -1;
    for(int i = 0; i < nv; i++)
    {
      if(values[j] == v[i])
      {
        fi = i;
        break;
      }
    }
    if(fi == -1)
    {
      v[nv++] = values[j];
    }
  }
  return nv;
}

//------------------------------------------------------------------------------
// NOTE - These types were pulled out of ClipField so they could be used in
//        some code that was moved out to handle degeneracies using partial
//        specialization rather than "if constexpr". Put it all back when
//        "if constexpr" is allowed. One nice side-effect is shorter symbol
//        names in the debugger.

using BitSet = std::uint32_t;

/*!
 * \brief Contains data that describes the number and size of zone fragments in the output.
 */
struct FragmentData
{
  IndexType m_finalNumZones {0};
  IndexType m_finalConnSize {0};
  axom::ArrayView<IndexType> m_fragmentsView {};
  axom::ArrayView<IndexType> m_fragmentsSizeView {};
  axom::ArrayView<IndexType> m_fragmentOffsetsView {};
  axom::ArrayView<IndexType> m_fragmentSizeOffsetsView {};
};
//------------------------------------------------------------------------------

#if defined(AXOM_CLIP_FILTER_DEGENERATES)
/*!
 * \brief Replace data in the input Conduit node with a denser version using the mask.
 *
 * \tparam ExecSpace The execution space.
 * \tparam DataView The type of data view that is operated on.
 *
 * \param n_src The Conduit node that contains the data.
 * \param srcView A view that wraps the input Conduit data.
 * \param newSize The new array size.
 * \param maskView The mask for valid data elements.
 * \param maskOffsetsView The offsets view to indicate where to write the new data.
 */
template <typename ExecSpace, typename DataView>
DataView filter(conduit::Node &n_src,
                DataView srcView,
                axom::IndexType newSize,
                axom::ArrayView<int> maskView,
                axom::ArrayView<int> maskOffsetsView)
{
  using value_type = typename DataView::value_type;
  namespace bputils = axom::mir::utilities::blueprint;

  // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
  utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
  const int conduitAllocatorID = c2a.getConduitAllocatorID();

  conduit::Node n_values;
  n_values.set_allocator(conduitAllocatorID);
  n_values.set(conduit::DataType(bputils::cpp2conduit<value_type>::id, newSize));
  auto valuesView = bputils::make_array_view<value_type>(n_values);
  const auto nValues = maskView.size();
  axom::for_all<ExecSpace>(
    nValues,
    AXOM_LAMBDA(axom::IndexType index) {
      if(maskView[index] > 0)
      {
        const auto destIndex = maskOffsetsView[index];
        valuesView[destIndex] = srcView[index];
      }
    });

  n_src.swap(n_values);
  return bputils::make_array_view<value_type>(n_src);
}

/// NOTE - Use partial specialization (instead of the cleaner "if constexpr")
///        for now to implement some 2D-specific behavior.

/*!
 * \brief Base template for degenerate removal.
 */
template <int NDIMS, typename ExecSpace, typename ConnectivityType>
struct DegenerateHandler
{
  /*!
   * \brief Set the size for the current fragment.
   *
   * \param fragmentsView The number of fragments for the szIndex zone.
   * \param connView The new connectivity.
   * \param sizesView The new mesh sizes.
   * \param szIndex The zone index currently being processed.
   * \param nidsThisFragment The number of node ids in the current fragment.
   * \param sizeIndex The write index for the sizes to output.
   * \param outputIndex The write index for the offsets to output.
   *
   * \return True if the fragment is degenerate, false otherwise.
   */
  AXOM_HOST_DEVICE
  static bool setSize(axom::ArrayView<IndexType> AXOM_UNUSED_PARAM(fragmentsView),
                      axom::ArrayView<ConnectivityType> AXOM_UNUSED_PARAM(connView),
                      axom::ArrayView<ConnectivityType> sizesView,
                      axom::IndexType AXOM_UNUSED_PARAM(szIndex),
                      int nIdsThisFragment,
                      int sizeIndex,
                      int &AXOM_UNUSED_PARAM(outputIndex))
  {
    sizesView[sizeIndex] = nIdsThisFragment;
    return false;
  }

  /*!
   * \brief In a previous stage, degenerate shapes were marked as having zero size.
   *        This method filters them out from the auxiliary arrays.
   *
   * \param fragmentData The fragments.
   * \param[inout] n_sizes The node that contains the sizes.
   * \param[inout] n_offsets The node that contains the offsets.
   * \param[inout] n_shapes The node that contains the shapes.
   * \param[inout] n_color The node that contains the color.
   * \param[inout] sizesView The view that wraps sizes (can change on output).
   * \param[inout] offsetsView The view that wraps offsets (can change on output).
   * \param[inout] shapesView The view that wraps shapes (can change on output).
   * \param[inout] colorView The view that wraps colors (can change on output).
   */
  static void filterZeroSizes(FragmentData &AXOM_UNUSED_PARAM(fragmentData),
                              conduit::Node &AXOM_UNUSED_PARAM(n_sizes),
                              conduit::Node &AXOM_UNUSED_PARAM(n_offsets),
                              conduit::Node &AXOM_UNUSED_PARAM(n_shapes),
                              conduit::Node &AXOM_UNUSED_PARAM(n_color),
                              axom::ArrayView<ConnectivityType> &AXOM_UNUSED_PARAM(sizesView),
                              axom::ArrayView<ConnectivityType> &AXOM_UNUSED_PARAM(offsetsView),
                              axom::ArrayView<ConnectivityType> &AXOM_UNUSED_PARAM(shapesView),
                              axom::ArrayView<int> &AXOM_UNUSED_PARAM(colorView))
  { }

  /*!
   * \brief Turns degenerate quads into triangles in-place.
   *
   * \param shapesUsed A BitSet that indicates which shapes are present in the mesh.
   * \param connView A view that contains the connectivity.
   * \param sizesView A view that contains the sizes.
   * \param offsetsView A view that contains the offsets.
   * \param shapesView A view that contains the shapes.
   */
  static BitSet quadtri(BitSet shapesUsed,
                        axom::ArrayView<ConnectivityType> AXOM_UNUSED_PARAM(connView),
                        axom::ArrayView<ConnectivityType> AXOM_UNUSED_PARAM(sizesView),
                        axom::ArrayView<ConnectivityType> AXOM_UNUSED_PARAM(offsetsView),
                        axom::ArrayView<ConnectivityType> AXOM_UNUSED_PARAM(shapesView))
  {
    return shapesUsed;
  }
};

/*!
 * \brief Partial specialization that implements some degeneracy handling for 2D meshes.
 */
template <typename ExecSpace, typename ConnectivityType>
struct DegenerateHandler<2, ExecSpace, ConnectivityType>
{
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

  /*!
   * \brief Set the size for the current fragment.
   *
   * \param fragmentsView The number of fragments for the szIndex zone.
   * \param connView The new connectivity.
   * \param sizesView The new mesh sizes.
   * \param szIndex The zone index currently being processed.
   * \param nidsThisFragment The number of node ids in the current fragment.
   * \param sizeIndex The write index for the sizes to output.
   * \param outputIndex The write index for the offsets to output.
   *
   * \return True if the fragment is degenerate, false otherwise.
   */
  AXOM_HOST_DEVICE
  static bool setSize(axom::ArrayView<IndexType> fragmentsView,
                      axom::ArrayView<ConnectivityType> connView,
                      axom::ArrayView<ConnectivityType> sizesView,
                      axom::IndexType szIndex,
                      int nIdsThisFragment,
                      int sizeIndex,
                      int &outputIndex)
  {
    const int connStart = outputIndex - nIdsThisFragment;

    // Check for degenerate
    const int nUniqueIds =
      internal::unique_count<ConnectivityType, 8>(connView.data() + connStart, nIdsThisFragment);
    const bool thisFragmentDegenerate = nUniqueIds < (nIdsThisFragment - 1);

    // Rewind the outputIndex so we don't emit it in the connectivity.
    if(thisFragmentDegenerate)
    {
      outputIndex = connStart;

      // There is one less fragment than we're expecting in the output.
      fragmentsView[szIndex] -= 1;
    }

    // Mark empty size.
    sizesView[sizeIndex] = thisFragmentDegenerate ? 0 : nIdsThisFragment;

    return thisFragmentDegenerate;
  }

  /*!
   * \brief In a previous stage, degenerate shapes were marked as having zero size.
   *        This method filters them out from the auxiliary arrays.
   *
   * \param fragmentData The fragments.
   * \param[inout] n_sizes The node that contains the sizes.
   * \param[inout] n_offsets The node that contains the offsets.
   * \param[inout] n_shapes The node that contains the shapes.
   * \param[inout] n_color The node that contains the color.
   * \param[inout] sizesView The view that wraps sizes (can change on output).
   * \param[inout] offsetsView The view that wraps offsets (can change on output).
   * \param[inout] shapesView The view that wraps shapes (can change on output).
   * \param[inout] colorView The view that wraps colors (can change on output).
   */
  static void filterZeroSizes(FragmentData &fragmentData,
                              conduit::Node &n_sizes,
                              conduit::Node &n_offsets,
                              conduit::Node &n_shapes,
                              conduit::Node &n_color,
                              axom::ArrayView<ConnectivityType> &sizesView,
                              axom::ArrayView<ConnectivityType> &offsetsView,
                              axom::ArrayView<ConnectivityType> &shapesView,
                              axom::ArrayView<int> &colorView)
  {
    AXOM_ANNOTATE_SCOPE("filterZeroSizes");

    // There were degenerates so the expected number of fragments per zone (m_fragmentsView)
    // was adjusted down. That means redoing the offsets. These need to be up
    // to date to handle zonal fields later.
    axom::exclusive_scan<ExecSpace>(fragmentData.m_fragmentsView, fragmentData.m_fragmentOffsetsView);

    // Use sizesView to make a mask that has 1's where size > 0.
    axom::IndexType nz = fragmentData.m_finalNumZones;
    axom::Array<int> mask(nz, nz, axom::execution_space<ExecSpace>::allocatorID());
    axom::Array<int> maskOffsets(nz, nz, axom::execution_space<ExecSpace>::allocatorID());
    auto maskView = mask.view();
    auto maskOffsetsView = maskOffsets.view();
    RAJA::ReduceSum<reduce_policy, axom::IndexType> mask_reduce(0);
    const axom::ArrayView<ConnectivityType> deviceSizesView = sizesView;
    axom::for_all<ExecSpace>(
      nz,
      AXOM_LAMBDA(axom::IndexType index) {
        const int ival = (deviceSizesView[index] > 0) ? 1 : 0;
        maskView[index] = ival;
        mask_reduce += ival;
      });
    const axom::IndexType filteredZoneCount = mask_reduce.get();

    // Make offsets
    axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

    // Filter sizes, shapes, color using the mask
    sizesView = filter<ExecSpace, axom::ArrayView<ConnectivityType>>(n_sizes,
                                                                     sizesView,
                                                                     filteredZoneCount,
                                                                     maskView,
                                                                     maskOffsetsView);
    offsetsView = filter<ExecSpace, axom::ArrayView<ConnectivityType>>(n_offsets,
                                                                       offsetsView,
                                                                       filteredZoneCount,
                                                                       maskView,
                                                                       maskOffsetsView);
    shapesView = filter<ExecSpace, axom::ArrayView<ConnectivityType>>(n_shapes,
                                                                      shapesView,
                                                                      filteredZoneCount,
                                                                      maskView,
                                                                      maskOffsetsView);
    colorView = filter<ExecSpace, axom::ArrayView<int>>(n_color,
                                                        colorView,
                                                        filteredZoneCount,
                                                        maskView,
                                                        maskOffsetsView);

    // Record the filtered size.
    fragmentData.m_finalNumZones = filteredZoneCount;
  }

  /*!
   * \brief Turns degenerate quads into triangles in-place.
   *
   * \param shapesUsed A BitSet that indicates which shapes are present in the mesh.
   * \param connView A view that contains the connectivity.
   * \param sizesView A view that contains the sizes.
   * \param offsetsView A view that contains the offsets.
   * \param shapesView A view that contains the shapes.
   */
  static BitSet quadtri(BitSet shapesUsed,
                        axom::ArrayView<ConnectivityType> connView,
                        axom::ArrayView<ConnectivityType> sizesView,
                        axom::ArrayView<ConnectivityType> offsetsView,
                        axom::ArrayView<ConnectivityType> shapesView)
  {
    if(axom::utilities::bitIsSet(shapesUsed, views::Quad_ShapeID))
    {
      AXOM_ANNOTATE_SCOPE("quadtri");
      const axom::IndexType numOutputZones = shapesView.size();
      RAJA::ReduceBitOr<reduce_policy, BitSet> shapesUsed_reduce(0);
      axom::for_all<ExecSpace>(
        numOutputZones,
        AXOM_LAMBDA(axom::IndexType index) {
          if(shapesView[index] == views::Quad_ShapeID)
          {
            const auto offset = offsetsView[index];
            ConnectivityType pts[4];
            int npts = 0;
            for(int current = 0; current < 4; current++)
            {
              int next = (current + 1) % 4;
              ConnectivityType curNode = connView[offset + current];
              ConnectivityType nextNode = connView[offset + next];
              if(curNode != nextNode)
              {
                pts[npts++] = curNode;
              }
            }

            if(npts == 3)
            {
              shapesView[index] = views::Tri_ShapeID;
              sizesView[index] = 3;
              connView[offset] = pts[0];
              connView[offset + 1] = pts[1];
              connView[offset + 2] = pts[2];
              // Repeat the last point (it won't be used though).
              connView[offset + 3] = pts[2];
            }
          }

          BitSet shapeBit {};
          axom::utilities::setBitOn(shapeBit, shapesView[index]);
          shapesUsed_reduce |= shapeBit;
        });
      // We redid shapesUsed reduction in case triangles appeared.
      shapesUsed = shapesUsed_reduce.get();
    }
    return shapesUsed;
  }
};
#endif

/*!
 * \brief Base template for handling fields on a strided-structured mesh. The
 *        default is that the mesh is not strided-structured so do nothing.
 *
 * \tparam enabled Whether the mesh is strided-structured.
 * \tparam ExecSpace The execution space.
 * \tparam TopologyView The topology view type.
 *
 * \note This was extracted from ClipField to remove some "if constexpr". Put
 *       it back someday.
 */
template <bool enabled, typename ExecSpace, typename TopologyView>
struct StridedStructuredFields
{
  /*!
   * \brief Slice an element field.
   *
   * \param topologyView The topology view.
   * \param slice Slice data.
   * \param n_field The field being sliced.
   * \param n_newField The node that will contain the new field.
   */
  static bool sliceElementField(
    const TopologyView &AXOM_UNUSED_PARAM(topologyView),
    const axom::mir::utilities::blueprint::SliceData &AXOM_UNUSED_PARAM(slice),
    const conduit::Node &AXOM_UNUSED_PARAM(n_field),
    conduit::Node &AXOM_UNUSED_PARAM(n_newField))
  {
    return false;
  }

  /*!
   * \brief Blend a vertex field.
   *
   * \param topologyView The topology view.
   * \param blend Blend data.
   * \param n_field The field being sliced.
   * \param n_newField The node that will contain the new field.
   */
  static bool blendVertexField(const TopologyView &AXOM_UNUSED_PARAM(topologyView),
                               const axom::mir::utilities::blueprint::BlendData &AXOM_UNUSED_PARAM(blend),
                               const conduit::Node &AXOM_UNUSED_PARAM(n_field),
                               conduit::Node &AXOM_UNUSED_PARAM(n_newField))
  {
    return false;
  }
};

/*!
 * \brief Partial specialization to handle fields on strided-structured mesh.
 *        This is the strided-structured case.
 *
 * \tparam ExecSpace The execution space.
 * \tparam TopologyView The topology view type.
 *
 * \note This was extracted from ClipField to remove some "if constexpr". Put
 *       it back someday.
 */
template <typename ExecSpace, typename TopologyView>
struct StridedStructuredFields<true, ExecSpace, TopologyView>
{
  /*!
   * \brief Slice an element field if the field is strided-structured.
   *
   * \param topologyView The topology view.
   * \param slice Slice data.
   * \param n_field The field being sliced.
   * \param n_newField The node that will contain the new field.
   */
  static bool sliceElementField(const TopologyView &topologyView,
                                const axom::mir::utilities::blueprint::SliceData &slice,
                                const conduit::Node &n_field,
                                conduit::Node &n_newField)
  {
    bool handled = false;
    if(n_field.has_path("offsets") && n_field.has_path("strides"))
    {
      using Indexing = typename TopologyView::IndexingPolicy;
      using IndexingPolicy = axom::mir::utilities::blueprint::SSElementFieldIndexing<Indexing>;
      IndexingPolicy indexing;
      indexing.m_indexing = topologyView.indexing();
      indexing.update(n_field);

      axom::mir::utilities::blueprint::FieldSlicer<ExecSpace, IndexingPolicy> s(indexing);
      s.execute(slice, n_field, n_newField);
      handled = true;
    }
    return handled;
  }

  /*!
   * \brief Blend a vertex field if the field is strided-structured.
   *
   * \param topologyView The topology view.
   * \param blend Blend data.
   * \param n_field The field being sliced.
   * \param n_newField The node that will contain the new field.
   */
  static bool blendVertexField(const TopologyView &topologyView,
                               const axom::mir::utilities::blueprint::BlendData &blend,
                               const conduit::Node &n_field,
                               conduit::Node &n_newField)
  {
    bool handled = false;
    if(n_field.has_path("offsets") && n_field.has_path("strides"))
    {
      // Make node indexing that the field blender can use.
      using Indexing = typename TopologyView::IndexingPolicy;
      using IndexingPolicy = axom::mir::utilities::blueprint::SSVertexFieldIndexing<Indexing>;
      IndexingPolicy indexing;
      indexing.m_topoIndexing = topologyView.indexing().expand();
      indexing.m_fieldIndexing = topologyView.indexing().expand();
      indexing.update(n_field);

      // If the topo and field offsets/strides are different then we need to go through
      // SSVertexFieldIndexing. Otherwise, we can let the normal case further below
      // handle the field.
      if(indexing.m_topoIndexing.m_offsets != indexing.m_fieldIndexing.m_offsets ||
         indexing.m_topoIndexing.m_strides != indexing.m_fieldIndexing.m_strides)
      {
        // Blend the field.
        axom::mir::utilities::blueprint::
          FieldBlender<ExecSpace, axom::mir::utilities::blueprint::SelectSubsetPolicy, IndexingPolicy>
            b(indexing);
        b.execute(blend, n_field, n_newField);
        handled = true;
      }
    }
    return handled;
  }
};

#if defined(AXOM_DEBUG_CLIP_FIELD)
/*!
 * \brief Print device views to std::out after moving data to the host.
 *
 * \param name The name of the view.
 * \param view The view to print.
 */
template <typename ViewType>
void printHost(const std::string &name, const ViewType &deviceView)
{
  using value_type = typename ViewType::value_type;
  int nn = deviceView.size();
  // Move data to host into temp array.
  value_type *host = new value_type[nn];
  axom::copy(host, deviceView.data(), sizeof(value_type) * nn);
  // Print
  std::cout << name << "[" << nn << "] = {";
  for(int ii = 0; ii < nn; ii++)
  {
    if(ii > 0)
    {
      std::cout << ", ";
    }
    std::cout << host[ii];
  }
  std::cout << "}" << std::endl;
  // Cleanup.
  delete[] host;
}
#endif

}  // end namespace internal

//------------------------------------------------------------------------------
/*!
 * \brief This class helps ClipField determine intersection cases and weights
 *        using a field designated by the options.
 */
template <typename ExecSpace, typename ConnectivityType>
class FieldIntersector
{
public:
  using ClipFieldType = float;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /*!
   * \brief This is a view class for FieldIntersector that can be used in device code.
   */
  struct View
  {
    /*!
     * \brief Given a zone index and the node ids that comprise the zone, return
     *        the appropriate clip case, taking into account the clip field and
     *        clip value.
     *
     * \param zoneIndex The zone index.
     * \param nodeIds A view containing node ids for the zone.
     */
    AXOM_HOST_DEVICE
    axom::IndexType determineClipCase(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex),
                                      const ConnectivityView &nodeIds) const
    {
      axom::IndexType clipcase = 0;
      for(IndexType i = 0; i < nodeIds.size(); i++)
      {
        const auto id = nodeIds[i];
        const auto value = m_clipFieldView[id] - m_clipValue;
        clipcase |= (value > 0) ? (1 << i) : 0;
      }
      return clipcase;
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
    ClipFieldType computeWeight(axom::IndexType AXOM_UNUSED_PARAM(zoneIndex),
                                ConnectivityType id0,
                                ConnectivityType id1) const
    {
      const ClipFieldType d0 = m_clipFieldView[id0];
      const ClipFieldType d1 = m_clipFieldView[id1];
      constexpr ClipFieldType tiny = 1.e-09;
      return axom::utilities::clampVal(
        axom::utilities::abs(m_clipValue - d0) / (axom::utilities::abs(d1 - d0) + tiny),
        ClipFieldType(0),
        ClipFieldType(1));
    }

    axom::ArrayView<ClipFieldType> m_clipFieldView {};
    ClipFieldType m_clipValue {};
  };

  /*!
   * \brief Initialize the object from options.
   * \param n_options The node that contains the options.
   * \param n_fields The node that contains fields.
   */
  void initialize(const conduit::Node &n_options, const conduit::Node &n_fields)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Get the clip field and value.
    ClipOptions opts(n_options);
    std::string clipFieldName = opts.clipField();
    m_view.m_clipValue = opts.clipValue();

    // Make sure the clipField is the right data type and store access to it in the view.
    const conduit::Node &n_clip_field = n_fields.fetch_existing(opts.clipField());
    const conduit::Node &n_clip_field_values = n_clip_field["values"];
    SLIC_ASSERT(n_clip_field["association"].as_string() == "vertex");
    if(n_clip_field_values.dtype().id() == bputils::cpp2conduit<ClipFieldType>::id)
    {
      // Make a view.
      m_view.m_clipFieldView = bputils::make_array_view<ClipFieldType>(n_clip_field_values);
    }
    else
    {
      // Convert to ClipFieldType.
      const IndexType n = static_cast<IndexType>(n_clip_field_values.dtype().number_of_elements());
      m_clipFieldData = axom::Array<ClipFieldType>(n, n, allocatorID);
      m_view.m_clipFieldView = m_clipFieldData.view();
      views::Node_to_ArrayView(n_clip_field_values,
                               [&](auto clipFieldViewSrc) { copyValues(clipFieldViewSrc); });
    }
  }

  /*!
   * \brief Determine the name of the topology on which to operate.
   * \param n_input The input mesh node.
   * \param n_options The clipping options.
   * \return The name of the toplogy on which to operate.
   */
  std::string getTopologyName(const conduit::Node &n_input, const conduit::Node &n_options) const
  {
    // Get the clipField's topo name.
    ClipOptions opts(n_options);
    const conduit::Node &n_fields = n_input.fetch_existing("fields");
    const conduit::Node &n_clipField = n_fields.fetch_existing(opts.clipField());
    return n_clipField["topology"].as_string();
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

  /*!
   * \brief Copy values from srcView into m_clipFieldData.
   *
   * \param srcView The source data view.
   */
  template <typename DataView>
  void copyValues(DataView srcView)
  {
    auto clipFieldView = m_clipFieldData.view();
    axom::for_all<ExecSpace>(
      srcView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        clipFieldView[index] = static_cast<ClipFieldType>(srcView[index]);
      });
  }

  axom::Array<ClipFieldType> m_clipFieldData {};
  View m_view {};
};

//------------------------------------------------------------------------------
/*!
 * \accelerated
 * \brief This class clips a topology using a field and puts the new topology into a new Conduit node.
 *
 * \tparam ExecSpace    The execution space where the compute-heavy kernels run.
 * \tparam TopologyView The topology view that can operate on the Blueprint topology.
 * \tparam CoordsetView The coordset view that can operate on the Blueprint coordset.
 * \tparam IntersectPolicy The intersector policy that can helps with cases and weights.
 * \tparam NamingPolicy The policy for making names from arrays of ids.
 */
template <typename ExecSpace,
          typename TopologyView,
          typename CoordsetView,
          typename IntersectPolicy =
            axom::mir::clipping::FieldIntersector<ExecSpace, typename TopologyView::ConnectivityType>,
          typename NamingPolicy = axom::mir::utilities::HashNaming<axom::IndexType>>
class ClipField
{
public:
  using BlendData = axom::mir::utilities::blueprint::BlendData;
  using SliceData = axom::mir::utilities::blueprint::SliceData;
  using ClipTableViews = axom::StackArray<axom::mir::clipping::TableView, 6>;
  using Intersector = IntersectPolicy;

  using BitSet = internal::BitSet;
  using KeyType = typename NamingPolicy::KeyType;
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using BlendGroupBuilderType = BlendGroupBuilder<ExecSpace, typename NamingPolicy::View>;
  using SelectedZones = typename axom::mir::utilities::blueprint::SelectedZones<ExecSpace>;
  using ClipFieldType = float;
  using ZoneType = typename TopologyView::ShapeType;

  /*!
   * \brief Constructor
   *
   * \param topoView A topology view suitable for the supplied topology.
   * \param coordsetView A coordset view suitable for the supplied coordset.
   *
   */
  ClipField(const TopologyView &topoView,
            const CoordsetView &coordsetView,
            const Intersector &intersector = Intersector())
    : m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_intersector(intersector)
    , m_clipTables()
    , m_naming()
  { }

  /*!
   * \brief Allow the user to pass in a NamingPolicy to use when making blend group names.
   *
   * \param naming A new naming policy object. 
   */
  void setNamingPolicy(NamingPolicy &naming) { m_naming = naming; }

  /*!
   * \brief Execute the clipping operation using the data stored in the specified \a clipField.
   *
   * \param[in] n_input The Conduit node that contains the topology, coordsets, and fields.
   * \param[in] n_options A Conduit node that contains clipping options.
   * \param[out] n_output A Conduit node that will hold the clipped output mesh. This should be a different node from \a n_input.
   *
   * \note The clipField field must currently be vertex-associated.
   */
  void execute(const conduit::Node &n_input, const conduit::Node &n_options, conduit::Node &n_output)
  {
    // Get the topo/coordset names in the input.
    ClipOptions opts(n_options);
    const std::string topoName = m_intersector.getTopologyName(n_input, n_options);
    const conduit::Node &n_topo = n_input.fetch_existing("topologies/" + topoName);
    const std::string coordsetName = n_topo["coordset"].as_string();
    const conduit::Node &n_coordset = n_input.fetch_existing("coordsets/" + coordsetName);
    const conduit::Node &n_fields = n_input.fetch_existing("fields");

    conduit::Node &n_newTopo = n_output["topologies/" + opts.topologyName(topoName)];
    conduit::Node &n_newCoordset = n_output["coordsets/" + opts.coordsetName(coordsetName)];
    conduit::Node &n_newFields = n_output["fields"];

    execute(n_topo, n_coordset, n_fields, n_options, n_newTopo, n_newCoordset, n_newFields);
  }

  /*!
   * \brief Execute the clipping operation using the data stored in the specified \a clipField.
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
    namespace bputils = axom::mir::utilities::blueprint;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    AXOM_ANNOTATE_SCOPE("ClipField");

    const std::string newTopologyName = n_newTopo.name();
    // Reset the output nodes just in case they've been reused.
    n_newTopo = conduit::Node();
    n_newCoordset = conduit::Node();
    n_newFields = conduit::Node();

    // Make the selected zones and get the size.
    ClipOptions opts(n_options);
    SelectedZones selectedZones(m_topologyView.numberOfZones(), n_options);
    const auto nzones = selectedZones.view().size();

    // Give the intersector a chance to further initialize.
    {
      AXOM_ANNOTATE_SCOPE("Initialize intersector");
      m_intersector.initialize(n_options, n_fields);
    }

    // Load clip table data and make views.
    m_clipTables.load(m_topologyView.dimension());
    ClipTableViews clipTableViews;
    createClipTableViews(clipTableViews, m_topologyView.dimension());

    // Allocate some memory and store views in ZoneData, FragmentData.
    AXOM_ANNOTATE_BEGIN("allocation");
    axom::Array<int> clipCases(nzones, nzones,
                               allocatorID);  // The clip case for a zone.
    axom::Array<BitSet> pointsUsed(
      nzones,
      nzones,
      allocatorID);  // Which points are used over all selected fragments in a zone
    ZoneData zoneData;
    zoneData.m_clipCasesView = clipCases.view();
    zoneData.m_pointsUsedView = pointsUsed.view();

    // Allocate some memory and store views in NodeData.
    NodeData nodeData;
#if defined(AXOM_REDUCE_BLEND_GROUPS)
    const auto nnodes = m_coordsetView.numberOfNodes();
    axom::Array<int> nodeUsed(nnodes, nnodes, allocatorID);
    nodeData.m_nodeUsedView = nodeUsed.view();
#endif

    // Allocate some memory and store views in FragmentData.
    axom::Array<IndexType> fragments(
      nzones,
      nzones,
      allocatorID);  // The number of fragments (child zones) produced for a zone.
    axom::Array<IndexType> fragmentsSize(
      nzones,
      nzones,
      allocatorID);  // The connectivity size for all selected fragments in a zone.
    axom::Array<IndexType> fragmentOffsets(nzones, nzones, allocatorID);
    axom::Array<IndexType> fragmentSizeOffsets(nzones, nzones, allocatorID);

    FragmentData fragmentData;
    fragmentData.m_fragmentsView = fragments.view();
    fragmentData.m_fragmentsSizeView = fragmentsSize.view();
    fragmentData.m_fragmentOffsetsView = fragmentOffsets.view();
    fragmentData.m_fragmentSizeOffsetsView = fragmentSizeOffsets.view();

    axom::Array<IndexType> blendGroups(nzones,
                                       nzones,
                                       allocatorID);  // Number of blend groups in a zone.
    axom::Array<IndexType> blendGroupsLen(nzones,
                                          nzones,
                                          allocatorID);  // Length of the blend groups in a zone.
    axom::Array<IndexType> blendOffset(nzones,
                                       nzones,
                                       allocatorID);  // Start of zone's blend group indices
    axom::Array<IndexType> blendGroupOffsets(
      nzones,
      nzones,
      allocatorID);  // Start of zone's blend group offsets in definitions.
    AXOM_ANNOTATE_END("allocation");

    // Make sure the naming policy knows the number of nodes.
    m_naming.setMaxId(m_coordsetView.numberOfNodes());

    // Make an object to help manage building the blend groups.
    BlendGroupBuilderType builder;
    builder.setNamingPolicy(m_naming.view());
    builder.setBlendGroupSizes(blendGroups.view(), blendGroupsLen.view());

    // Compute sizes and offsets
    computeSizes(clipTableViews, builder, zoneData, nodeData, fragmentData, opts, selectedZones);
    computeFragmentSizes(fragmentData, selectedZones);
    computeFragmentOffsets(fragmentData);

    // Compute original node count that we're preserving, make node maps.
#if defined(AXOM_REDUCE_BLEND_GROUPS)
    const int compactSize = countOriginalNodes(nodeData);
    axom::Array<IndexType> compactNodes(compactSize, compactSize, allocatorID);
    axom::Array<IndexType> oldNodeToNewNode(nnodes, nnodes, allocatorID);
    nodeData.m_originalIdsView = compactNodes.view();
    nodeData.m_oldNodeToNewNodeView = oldNodeToNewNode.view();
    createNodeMaps(nodeData);

    nodeUsed.clear();
    nodeData.m_nodeUsedView = axom::ArrayView<int>();
#endif

    // Further initialize the blend group builder.
    IndexType blendGroupsSize = 0, blendGroupLenSize = 0;
    builder.computeBlendGroupSizes(blendGroupsSize, blendGroupLenSize);
    builder.setBlendGroupOffsets(blendOffset.view(), blendGroupOffsets.view());
    builder.computeBlendGroupOffsets();

    // Allocate memory for blend groups.
    AXOM_ANNOTATE_BEGIN("allocation2");
    axom::Array<KeyType> blendNames(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<IndexType> blendGroupSizes(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<IndexType> blendGroupStart(blendGroupsSize, blendGroupsSize, allocatorID);
    axom::Array<IndexType> blendIds(blendGroupLenSize, blendGroupLenSize, allocatorID);
    axom::Array<float> blendCoeff(blendGroupLenSize, blendGroupLenSize, allocatorID);

    // Make the blend groups.
    builder.setBlendViews(blendNames.view(),
                          blendGroupSizes.view(),
                          blendGroupStart.view(),
                          blendIds.view(),
                          blendCoeff.view());
    AXOM_ANNOTATE_END("allocation2");
    makeBlendGroups(clipTableViews, builder, zoneData, opts, selectedZones);

    // Make the blend groups unique
    axom::Array<KeyType> uNames;
    axom::Array<axom::IndexType> uIndices;
#if defined(AXOM_REDUCE_BLEND_GROUPS)
    axom::Array<KeyType> newUniqueNames;
    axom::Array<axom::IndexType> newUniqueIndices;
#endif
    {
      AXOM_ANNOTATE_SCOPE("unique");
      axom::mir::utilities::Unique<ExecSpace, KeyType>::execute(builder.blendNames(), uNames, uIndices);
      builder.setUniqueNames(uNames.view(), uIndices.view());

#if defined(AXOM_REDUCE_BLEND_GROUPS)
      // Filter the unique names/indices to remove single node blend groups.
      builder.filterUnique(newUniqueNames, newUniqueIndices);
      uNames.clear();
      uIndices.clear();
#endif
    }

    // Make BlendData.
    bputils::BlendData blend = builder.makeBlendData();
    blend.m_originalIdsView = nodeData.m_originalIdsView;

    // Make the clipped mesh
    makeTopology(clipTableViews,
                 builder,
                 zoneData,
                 nodeData,
                 fragmentData,
                 opts,
                 selectedZones,
                 newTopologyName,
                 n_newTopo,
                 n_newCoordset,
                 n_newFields);

    // Make the coordset
    makeCoordset(blend, n_coordset, n_newCoordset);

    // Get the fields that we want to process.
    std::map<std::string, std::string> fieldsToProcess;
    int numElementFields = 0;
    if(opts.fields(fieldsToProcess))
    {
      // Fields were present in the options. Count the element fields.
      for(auto it = fieldsToProcess.begin(); it != fieldsToProcess.end(); it++)
      {
        const conduit::Node &n_field = n_fields.fetch_existing(it->first);
        if(n_field.fetch_existing("topology").as_string() == n_topo.name())
        {
          numElementFields +=
            (n_field.fetch_existing("association").as_string() == "element") ? 1 : 0;
        }
      }
    }
    else
    {
      // Fields were not present in the options. Select all fields that have the same topology as n_topo.
      for(conduit::index_t i = 0; i < n_fields.number_of_children(); i++)
      {
        const conduit::Node &n_field = n_fields[i];
        if(n_field.fetch_existing("topology").as_string() == n_topo.name())
        {
          numElementFields +=
            (n_field.fetch_existing("association").as_string() == "element") ? 1 : 0;

          fieldsToProcess[n_field.name()] = n_field.name();
        }
      }
    }
    const std::string newNodes = opts.newNodesField();
    if(!newNodes.empty() && n_fields.has_child(newNodes))
    {
      fieldsToProcess[newNodes] = newNodes;
    }

    // Make slice indices if we have element fields.
    bputils::SliceData slice;
    axom::Array<IndexType> sliceIndices;
    if(numElementFields > 0)
    {
      AXOM_ANNOTATE_SCOPE("sliceIndices");
      sliceIndices = axom::Array<IndexType>(fragmentData.m_finalNumZones,
                                            fragmentData.m_finalNumZones,
                                            allocatorID);
      auto sliceIndicesView = sliceIndices.view();

      // Fill in sliceIndicesView.
      const auto selectedZonesView = selectedZones.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          const auto zoneIndex = selectedZonesView[index];
          const auto start = fragmentData.m_fragmentOffsetsView[index];
          for(int i = 0; i < fragmentData.m_fragmentsView[index]; i++)
            sliceIndicesView[start + i] = zoneIndex;
        });
      slice.m_indicesView = sliceIndicesView;
    }

    makeFields(blend, slice, newTopologyName, fieldsToProcess, n_fields, n_newFields);

    makeOriginalElements(fragmentData, opts, selectedZones, n_fields, n_newTopo, n_newFields);

    markNewNodes(blend, newNodes, newTopologyName, n_newFields);
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif
  using FragmentData = internal::FragmentData;

  /*!
   * \brief Contains some per-zone data that we want to hold onto between methods.
   */
  struct ZoneData
  {
    axom::ArrayView<int> m_clipCasesView {};
    axom::ArrayView<BitSet> m_pointsUsedView {};
  };

  /*!
   * \brief Contains some per-node data that we want to hold onto between methods.
   */
  struct NodeData
  {
    axom::ArrayView<int> m_nodeUsedView {};
    axom::ArrayView<IndexType> m_oldNodeToNewNodeView {};
    axom::ArrayView<IndexType> m_originalIdsView {};
  };

  /*!
   * \brief Make a bitset that indicates the parts of the selection that are selected.
   */
  int getSelection(const ClipOptions &opts) const
  {
    int selection = 0;
    if(opts.inside()) axom::utilities::setBitOn(selection, 0);
    if(opts.outside()) axom::utilities::setBitOn(selection, 1);
    SLIC_ASSERT(selection > 0);
    return selection;
  }

  /*!
   * \brief Create views for the clip tables of various shapes.
   *
   * \param[out] views The views array that will contain the table views.
   * \param dimension The dimension the topology (so we can load a subset of tables)
   */
  void createClipTableViews(ClipTableViews &views, int dimension)
  {
    AXOM_ANNOTATE_SCOPE("createClipTableViews");
    if(dimension == -1 || dimension == 2)
    {
      views[internal::getClipTableIndex(views::Tri_ShapeID)] = m_clipTables[ST_TRI].view();
      views[internal::getClipTableIndex(views::Quad_ShapeID)] = m_clipTables[ST_QUA].view();
    }
    if(dimension == -1 || dimension == 3)
    {
      views[internal::getClipTableIndex(views::Tet_ShapeID)] = m_clipTables[ST_TET].view();
      views[internal::getClipTableIndex(views::Pyramid_ShapeID)] = m_clipTables[ST_PYR].view();
      views[internal::getClipTableIndex(views::Wedge_ShapeID)] = m_clipTables[ST_WDG].view();
      views[internal::getClipTableIndex(views::Hex_ShapeID)] = m_clipTables[ST_HEX].view();
    }
  }

  /*!
   * \brief Iterate over zones and their respective fragments to determine sizes
   *        for fragments and blend groups.
   *
   * \param[in] clipTableViews An object that holds views of the clipping table data.
   * \param[in] builder This object holds views to blend group data and helps with building/access.
   * \param[in] zoneData This object holds views to per-zone data.
   * \param[in] nodeData This object holds views to per-node data.
   * \param[in] fragmentData This object holds views to per-fragment data.
   * \param[inout] opts Clipping options.
   *
   * \note Objects that we need to capture into kernels are passed by value (they only contain views anyway). Data can be modified through the views.
   */
  void computeSizes(ClipTableViews clipTableViews,
                    BlendGroupBuilderType builder,
                    ZoneData zoneData,
                    NodeData nodeData,
                    FragmentData fragmentData,
                    const ClipOptions &opts,
                    const SelectedZones &selectedZones) const
  {
    AXOM_ANNOTATE_SCOPE("computeSizes");
    const auto selection = getSelection(opts);

    auto blendGroupsView = builder.state().m_blendGroupsView;
    auto blendGroupsLenView = builder.state().m_blendGroupsLenView;

    // Initialize nodeUsed data for nodes.
    axom::for_all<ExecSpace>(
      nodeData.m_nodeUsedView.size(),
      AXOM_LAMBDA(axom::IndexType index) { nodeData.m_nodeUsedView[index] = 0; });

    const auto deviceIntersector = m_intersector.view();

    const TopologyView deviceTopologyView(m_topologyView);
    const auto selectedZonesView = selectedZones.view();
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const auto zone = deviceTopologyView.zone(zoneIndex);

        // Get the clip case for the current zone.
        const auto clipcase = deviceIntersector.determineClipCase(zoneIndex, zone.getIds());
        zoneData.m_clipCasesView[szIndex] = clipcase;

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = internal::getClipTableIndex(zone.id());
        const auto &ctView = clipTableViews[clipTableIndex];

        int thisBlendGroups = 0;      // The number of blend groups produced in this case.
        int thisBlendGroupLen = 0;    // The total length of the blend groups.
        int thisFragments = 0;        // The number of zone fragments produced in this case.
        int thisFragmentsNumIds = 0;  // The number of points used to make all the fragment zones.
        BitSet ptused = 0;            // A bitset indicating which ST_XX nodes are used.

        auto it = ctView.begin(clipcase);
        const auto end = ctView.end(clipcase);
        for(; it != end; it++)
        {
          // Get the current shape in the clip case.
          const auto fragment = *it;

          if(fragment[0] == ST_PNT)
          {
            if(internal::generatedPointIsSelected(fragment[2], selection))
            {
              const int nIds = static_cast<int>(fragment[3]);

              for(int ni = 0; ni < nIds; ni++)
              {
                const auto pid = fragment[4 + ni];

                // Increase the blend size to include this center point.
                if(pid <= P7)
                {
                  // corner point
                  thisBlendGroupLen++;
                }
                else if(pid >= EA && pid <= EL)
                {
                  // edge point
                  thisBlendGroupLen += 2;
                }
              }

              // This center or face point counts as a blend group.
              thisBlendGroups++;

              // Mark the point used.
              axom::utilities::setBitOn(ptused, N0 + fragment[1]);
            }
          }
          else
          {
            if(internal::shapeIsSelected(fragment[1], selection))
            {
              thisFragments++;
              const int nIdsThisFragment = fragment.size() - 2;
              thisFragmentsNumIds += nIdsThisFragment;

              // Mark the points this fragment used.
              for(int i = 2; i < fragment.size(); i++)
              {
                axom::utilities::setBitOn(ptused, fragment[i]);
              }
            }
          }
        }

        // Save the flags for the points that were used in this zone
        zoneData.m_pointsUsedView[szIndex] = ptused;

#if defined(AXOM_REDUCE_BLEND_GROUPS)
        // NOTE: We are not going to emit blend groups for P0..P7 points.

        // If the zone uses a node, set that node in nodeUsedView.
        for(IndexType pid = P0; pid <= P7; pid++)
        {
          if(axom::utilities::bitIsSet(ptused, pid))
          {
            const auto nodeId = zone.getId(pid);
            // NOTE: Multiple threads may write to this node but they all write the same value.
            nodeData.m_nodeUsedView[nodeId] = 1;
          }
        }
#else
        // Count which points in the original cell are used.
        for(IndexType pid = P0; pid <= P7; pid++)
        {
          const int incr = axom::utilities::bitIsSet(ptused, pid) ? 1 : 0;

          thisBlendGroupLen += incr;  // {p0}
          thisBlendGroups += incr;
        }
#endif

        // Count edges that are used.
        for(IndexType pid = EA; pid <= EL; pid++)
        {
          const int incr = axom::utilities::bitIsSet(ptused, pid) ? 1 : 0;

          thisBlendGroupLen += 2 * incr;  // {p0 p1}
          thisBlendGroups += incr;
        }

        // Save the results.
        fragmentData.m_fragmentsView[szIndex] = thisFragments;
        fragmentData.m_fragmentsSizeView[szIndex] = thisFragmentsNumIds;

        // Set blend group sizes for this zone.
        blendGroupsView[szIndex] = thisBlendGroups;
        blendGroupsLenView[szIndex] = thisBlendGroupLen;
      });  // for_selected_zones

#if defined(AXOM_DEBUG_CLIP_FIELD)
    std::cout << "------------------------ computeSizes ------------------------" << std::endl;
    internal::printHost("fragmentData.m_fragmentsView", fragmentData.m_fragmentsView);
    internal::printHost("fragmentData.m_fragmentsSizeView", fragmentData.m_fragmentsSizeView);
    internal::printHost("blendGroupsView", blendGroupsView);
    internal::printHost("blendGroupsLenView", blendGroupsLenView);
    internal::printHost("zoneData.m_pointsUsedView", zoneData.m_pointsUsedView);
    internal::printHost("zoneData.m_clipCasesView", zoneData.m_clipCasesView);
    std::cout << "--------------------------------------------------------------" << std::endl;
#endif
  }

  /*!
   * \brief Compute the total number of fragments and their size.
   *
   * \param[inout] fragmentData The object that contains data about the zone fragments.
   */
  void computeFragmentSizes(FragmentData &fragmentData, const SelectedZones &selectedZones) const
  {
    AXOM_ANNOTATE_SCOPE("computeFragmentSizes");
    const auto nzones = selectedZones.view().size();

    // Sum the number of fragments.
    RAJA::ReduceSum<reduce_policy, IndexType> fragment_sum(0);
    const auto fragmentsView = fragmentData.m_fragmentsView;
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType szIndex) { fragment_sum += fragmentsView[szIndex]; });
    fragmentData.m_finalNumZones = fragment_sum.get();

    // Sum the fragment connectivity sizes.
    RAJA::ReduceSum<reduce_policy, IndexType> fragment_nids_sum(0);
    const auto fragmentsSizeView = fragmentData.m_fragmentsSizeView;
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType szIndex) { fragment_nids_sum += fragmentsSizeView[szIndex]; });
    fragmentData.m_finalConnSize = fragment_nids_sum.get();
  }

  /*!
   * \brief Compute fragment offsets.
   *
   * \param[inout] fragmentData The object that contains data about the zone fragments.
   */
  void computeFragmentOffsets(FragmentData &fragmentData) const
  {
    AXOM_ANNOTATE_SCOPE("computeFragmentOffsets");
    axom::exclusive_scan<ExecSpace>(fragmentData.m_fragmentsView, fragmentData.m_fragmentOffsetsView);
    axom::exclusive_scan<ExecSpace>(fragmentData.m_fragmentsSizeView,
                                    fragmentData.m_fragmentSizeOffsetsView);

#if defined(AXOM_DEBUG_CLIP_FIELD)
    std::cout << "------------------------ computeFragmentOffsets "
                 "------------------------"
              << std::endl;
    internal::printHost("fragmentData.m_fragmentOffsetsView", fragmentData.m_fragmentOffsetsView);
    internal::printHost("fragmentData.m_fragmentSizeOffsetsView",
                        fragmentData.m_fragmentSizeOffsetsView);
    std::cout << "-------------------------------------------------------------"
                 "-----------"
              << std::endl;
#endif
  }

#if defined(AXOM_REDUCE_BLEND_GROUPS)
  /*!
   * \brief Counts the number of original nodes used by the selected fragments.
   *
   * \param nodeData The node data (passed by value on purpose)
   * \return The number of original nodes used by selected fragments.
   */
  int countOriginalNodes(NodeData nodeData) const
  {
    AXOM_ANNOTATE_SCOPE("countOriginalNodes");
    // Count the number of original nodes we'll use directly.
    RAJA::ReduceSum<reduce_policy, int> nUsed_reducer(0);
    const auto nodeUsedView = nodeData.m_nodeUsedView;
    axom::for_all<ExecSpace>(
      nodeUsedView.size(),
      AXOM_LAMBDA(axom::IndexType index) { nUsed_reducer += nodeUsedView[index]; });
    return nUsed_reducer.get();
  }

  /*!
   * \brief Creates the node lists/maps.
   *
   * \param nodeData The node data that contains views where the node data is stored.
   */
  void createNodeMaps(NodeData nodeData) const
  {
    AXOM_ANNOTATE_SCOPE("createNodeMaps");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Make offsets into a compact array.
    const auto nnodes = nodeData.m_nodeUsedView.size();
    axom::Array<int> nodeOffsets(nnodes, nnodes, allocatorID);
    auto nodeOffsetsView = nodeOffsets.view();
    axom::exclusive_scan<ExecSpace>(nodeData.m_nodeUsedView, nodeOffsetsView);

    // Make the compact node list and oldToNew map.
    axom::for_all<ExecSpace>(
      nnodes,
      AXOM_LAMBDA(axom::IndexType index) {
        IndexType newId = 0;
        if(nodeData.m_nodeUsedView[index] > 0)
        {
          nodeData.m_originalIdsView[nodeOffsetsView[index]] = index;
          newId = nodeOffsetsView[index];
        }
        nodeData.m_oldNodeToNewNodeView[index] = newId;
      });

  #if defined(AXOM_DEBUG_CLIP_FIELD)
    std::cout << "---------------------------- createNodeMaps "
                 "----------------------------"
              << std::endl;
    internal::printHost("nodeData.m_nodeUsedView", nodeData.m_nodeUsedView);
    internal::printHost("nodeData.m_originalIdsView", nodeData.m_originalIdsView);
    internal::printHost("nodeData.m_oldNodeToNewNodeView", nodeData.m_oldNodeToNewNodeView);
    std::cout << "-------------------------------------------------------------"
                 "-----------"
              << std::endl;
  #endif
  }
#endif

  /*!
   * \brief Fill in the data for the blend group views.
   *
   * \param[in] clipTableViews An object that holds views of the clipping table data.
   * \param[in] builder This object holds views to blend group data and helps with building/access.
   * \param[in] zoneData This object holds views to per-zone data.
   * \param[inout] opts Clipping options.
   *
   * \note Objects that we need to capture into kernels are passed by value (they only contain views anyway). Data can be modified through the views.
   */
  void makeBlendGroups(ClipTableViews clipTableViews,
                       BlendGroupBuilderType builder,
                       ZoneData zoneData,
                       const ClipOptions &opts,
                       const SelectedZones &selectedZones) const
  {
    AXOM_ANNOTATE_SCOPE("makeBlendGroups");
    const auto selection = getSelection(opts);

    const auto deviceIntersector = m_intersector.view();
    const TopologyView deviceTopologyView(m_topologyView);
    const auto selectedZonesView = selectedZones.view();
    axom::for_all<ExecSpace>(
      selectedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = selectedZonesView[szIndex];
        const auto zone = deviceTopologyView.zone(zoneIndex);

        // Get the clip case for the current zone.
        const auto clipcase = zoneData.m_clipCasesView[szIndex];

        // Iterate over the shapes in this clip case to determine the number of blend groups.
        const auto clipTableIndex = internal::getClipTableIndex(zone.id());
        const auto &ctView = clipTableViews[clipTableIndex];

        // These are the points used in this zone's fragments.
        const BitSet ptused = zoneData.m_pointsUsedView[szIndex];

        // Get the blend groups for this zone.
        auto groups = builder.blendGroupsForZone(szIndex);

        auto it = ctView.begin(clipcase);
        const auto end = ctView.end(clipcase);
        for(; it != end; it++)
        {
          // Get the current shape in the clip case.
          const auto fragment = *it;

          if(fragment[0] == ST_PNT)
          {
            if(internal::generatedPointIsSelected(fragment[2], selection))
            {
              const int nIds = static_cast<int>(fragment[3]);
              const auto one_over_n = 1.f / static_cast<float>(nIds);

              groups.beginGroup();
              for(int ni = 0; ni < nIds; ni++)
              {
                const auto ptid = fragment[4 + ni];

                // Add the point to the blend group.
                if(ptid <= P7)
                {
                  // corner point.
                  groups.add(zone.getId(ptid), one_over_n);
                }
                else if(ptid >= EA && ptid <= EL)
                {
                  // edge point.
                  const auto edgeIndex = ptid - EA;
                  const auto edge = zone.getEdge(edgeIndex);
                  const auto id0 = zone.getId(edge[0]);
                  const auto id1 = zone.getId(edge[1]);

                  // Figure out the blend for edge.
                  const auto t = deviceIntersector.computeWeight(zoneIndex, id0, id1);

                  groups.add(id0, one_over_n * (1.f - t));
                  groups.add(id1, one_over_n * t);
                }
              }
              groups.endGroup();
            }
          }
        }
#if !defined(AXOM_REDUCE_BLEND_GROUPS)
        // Add blend group for each original point that was used.
        // NOTE - this can add a lot of blend groups with 1 node.
        for(IndexType pid = P0; pid <= P7; pid++)
        {
          if(axom::utilities::bitIsSet(ptused, pid))
          {
            groups.beginGroup();
            groups.add(zone.getId(pid), 1.f);
            groups.endGroup();
          }
        }
#endif
        // Add blend group for each edge point that was used.
        for(IndexType pid = EA; pid <= EL; pid++)
        {
          if(axom::utilities::bitIsSet(ptused, pid))
          {
            const auto edgeIndex = pid - EA;
            const auto edge = zone.getEdge(edgeIndex);
            const auto id0 = zone.getId(edge[0]);
            const auto id1 = zone.getId(edge[1]);

            // Figure out the blend for edge.
            const auto t = deviceIntersector.computeWeight(zoneIndex, id0, id1);

            // Close to the endpoints, just count the edge blend group
            // as an endpoint to ensure better blend group matching later.
            constexpr decltype(t) LOWER = 1.e-4;
            constexpr decltype(t) UPPER = 1. - LOWER;
            groups.beginGroup();
            if(t < LOWER)
            {
              groups.add(id0, 1.f);
            }
            else if(t > UPPER)
            {
              groups.add(id1, 1.f);
            }
            else
            {
              groups.add(id0, 1.f - t);
              groups.add(id1, t);
            }
            groups.endGroup();
          }
        }
      });
  }

  /*!
   * \brief Make the clipped mesh topology.
   *
   * \param[in] clipTableViews An object that holds views of the clipping table data.
   * \param[in] builder This object holds views to blend group data and helps with building/access.
   * \param[in] zoneData This object holds views to per-zone data.
   * \param[in] nodeData This object holds views to per-node data.
   * \param[in] fragmentData This object holds views to per-fragment data.
   * \param[in] opts Clipping options.
   * \param[in] selectedZones The selected zones.
   * \param[in] newTopologyName The name of the new topology.
   * \param[out] n_newTopo The node that will contain the new topology.
   * \param[out] n_newCoordset The node that will contain the new coordset.
   * \param[out] n_newFields The node that will contain the new fields.
   *
   * \note Objects that we need to capture into kernels are passed by value (they only contain views anyway). Data can be modified through the views.
   */
  void makeTopology(ClipTableViews clipTableViews,
                    BlendGroupBuilderType builder,
                    ZoneData zoneData,
                    NodeData nodeData,
                    FragmentData fragmentData,
                    const ClipOptions &opts,
                    const SelectedZones &selectedZones,
                    const std::string &newTopologyName,
                    conduit::Node &n_newTopo,
                    conduit::Node &n_newCoordset,
                    conduit::Node &n_newFields) const
  {
    AXOM_ANNOTATE_SCOPE("makeTopology");
    namespace bputils = axom::mir::utilities::blueprint;
    constexpr auto connTypeID = bputils::cpp2conduit<ConnectivityType>::id;
    const auto selection = getSelection(opts);

    AXOM_ANNOTATE_BEGIN("allocation");
    n_newTopo["type"] = "unstructured";
    n_newTopo["coordset"] = n_newCoordset.name();

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    // _mir_utilities_c2a_begin
    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    const int conduitAllocatorID = c2a.getConduitAllocatorID();

    // Allocate connectivity.
    conduit::Node &n_conn = n_newTopo["elements/connectivity"];
    n_conn.set_allocator(conduitAllocatorID);
    n_conn.set(conduit::DataType(connTypeID, fragmentData.m_finalConnSize));
    auto connView = bputils::make_array_view<ConnectivityType>(n_conn);
    // _mir_utilities_c2a_end

    // Allocate shapes.
    conduit::Node &n_shapes = n_newTopo["elements/shapes"];
    n_shapes.set_allocator(conduitAllocatorID);
    n_shapes.set(conduit::DataType(connTypeID, fragmentData.m_finalNumZones));
    auto shapesView = bputils::make_array_view<ConnectivityType>(n_shapes);

    // Allocate sizes.
    conduit::Node &n_sizes = n_newTopo["elements/sizes"];
    n_sizes.set_allocator(conduitAllocatorID);
    n_sizes.set(conduit::DataType(connTypeID, fragmentData.m_finalNumZones));
    auto sizesView = bputils::make_array_view<ConnectivityType>(n_sizes);

    // Allocate offsets.
    conduit::Node &n_offsets = n_newTopo["elements/offsets"];
    n_offsets.set_allocator(conduitAllocatorID);
    n_offsets.set(conduit::DataType(connTypeID, fragmentData.m_finalNumZones));
    auto offsetsView = bputils::make_array_view<ConnectivityType>(n_offsets);

    // Allocate a color variable to keep track of the "color" of the fragments.
    conduit::Node &n_color = n_newFields[opts.colorField()];
    n_color["topology"] = newTopologyName;
    n_color["association"] = "element";
    conduit::Node &n_color_values = n_color["values"];
    n_color_values.set_allocator(conduitAllocatorID);
    n_color_values.set(conduit::DataType::int32(fragmentData.m_finalNumZones));
    auto colorView = bputils::make_array_view<int>(n_color_values);

    // Fill in connectivity values in case we leave empty slots later.
    axom::for_all<ExecSpace>(
      connView.size(),
      AXOM_LAMBDA(axom::IndexType index) { connView[index] = 0; });

#if defined(AXOM_DEBUG_CLIP_FIELD)
    // Initialize the values beforehand. For debugging.
    axom::for_all<ExecSpace>(
      shapesView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        shapesView[index] = -2;
        sizesView[index] = -3;
        offsetsView[index] = -4;
        colorView[index] = -5;
      });
#endif
    AXOM_ANNOTATE_END("allocation");

    // Here we fill in the new connectivity, sizes, shapes.
    // We get the node ids from the unique blend names, de-duplicating points when making the new connectivity.
    //
    // NOTE: During development, I ran into problems with this kernel not executing
    //       due to point_2_new being too large. The solution was to reduce the values
    //       for EA-EL, N0-N3 to shrink the array to the point where it can fit in
    //       memory available to the thread.
    //
#if defined(AXOM_CLIP_FILTER_DEGENERATES)
    RAJA::ReduceBitOr<reduce_policy, BitSet> degenerates_reduce(0);
#endif
    {
      AXOM_ANNOTATE_SCOPE("build");
      const auto origSize = nodeData.m_originalIdsView.size();

      const TopologyView deviceTopologyView(m_topologyView);
      const auto selectedZonesView = selectedZones.view();
      axom::for_all<ExecSpace>(
        selectedZonesView.size(),
        AXOM_LAMBDA(axom::IndexType szIndex) {
          const auto zoneIndex = selectedZonesView[szIndex];
          const auto zone = deviceTopologyView.zone(zoneIndex);

          // If there are no fragments, return from lambda.
          if(fragmentData.m_fragmentsView[szIndex] == 0) return;

          // Seek to the start of the blend groups for this zone.
          auto groups = builder.blendGroupsForZone(szIndex);

          // Go through the points in the order they would have been added as blend
          // groups, get their blendName, and then overall index of that blendName
          // in uNames, the unique list of new dof names. That will be their index
          // in the final points.
          const BitSet ptused = zoneData.m_pointsUsedView[szIndex];
          ConnectivityType point_2_new[N3 + 1];
          for(BitSet pid = N0; pid <= N3; pid++)
          {
            if(axom::utilities::bitIsSet(ptused, pid))
            {
              point_2_new[pid] = origSize + groups.uniqueBlendGroupIndex();
              groups++;
            }
          }
#if defined(AXOM_REDUCE_BLEND_GROUPS)
          // For single nodes, we did not make a blend group. We look up the new
          // node id from nodeData.m_oldNodeToNewNodeView.
          for(BitSet pid = P0; pid <= P7; pid++)
          {
            if(axom::utilities::bitIsSet(ptused, pid))
            {
              const auto nodeId = zone.getId(pid);
              point_2_new[pid] = nodeData.m_oldNodeToNewNodeView[nodeId];
            }
          }
#else
          for(BitSet pid = P0; pid <= P7; pid++)
          {
            if(axom::utilities::bitIsSet(ptused, pid))
            {
              point_2_new[pid] = origSize + groups.uniqueBlendGroupIndex();
              groups++;
            }
          }
#endif
          for(BitSet pid = EA; pid <= EL; pid++)
          {
            if(axom::utilities::bitIsSet(ptused, pid))
            {
#if defined(AXOM_REDUCE_BLEND_GROUPS)
              // There is a chance that the edge blend group was emitted with a
              // single node if the edge was really close to a corner node.
              if(groups.size() == 1)
              {
                const auto nodeId = groups.id(0);
                point_2_new[pid] = nodeData.m_oldNodeToNewNodeView[nodeId];
              }
              else
              {
                point_2_new[pid] = origSize + groups.uniqueBlendGroupIndex();
              }
#else
              point_2_new[pid] = origSize + groups.uniqueBlendGroupIndex();
#endif
              groups++;
            }
          }

          // This is where the output fragment connectivity start for this zone
          int outputIndex = fragmentData.m_fragmentSizeOffsetsView[szIndex];
          // This is where the output fragment sizes/shapes start for this zone.
          int sizeIndex = fragmentData.m_fragmentOffsetsView[szIndex];
#if defined(AXOM_CLIP_FILTER_DEGENERATES)
          bool degenerates = false;
#endif
          // Iterate over the selected fragments and emit connectivity for them.
          const auto clipcase = zoneData.m_clipCasesView[szIndex];
          const auto clipTableIndex = internal::getClipTableIndex(zone.id());
          const auto ctView = clipTableViews[clipTableIndex];
          auto it = ctView.begin(clipcase);
          const auto end = ctView.end(clipcase);
          for(; it != end; it++)
          {
            // Get the current shape in the clip case.
            const auto fragment = *it;
            const auto fragmentShape = fragment[0];

            if(fragmentShape != ST_PNT)
            {
              if(internal::shapeIsSelected(fragment[1], selection))
              {
                // Output the nodes used in this zone.
                const int fragmentSize = fragment.size();
                offsetsView[sizeIndex] = outputIndex;
                for(int i = 2; i < fragmentSize; i++)
                  connView[outputIndex++] = point_2_new[fragment[i]];

                const auto nIdsThisFragment = fragmentSize - 2;
#if defined(AXOM_CLIP_FILTER_DEGENERATES)
                // Set the output zone size, checking to see whether it is degenerate.
                degenerates |=
                  internal::DegenerateHandler<TopologyView::dimension(), ExecSpace, ConnectivityType>::setSize(
                    fragmentData.m_fragmentsView,
                    connView,
                    sizesView,
                    szIndex,
                    nIdsThisFragment,
                    sizeIndex,
                    outputIndex);
#else
                sizesView[sizeIndex] = nIdsThisFragment;
#endif
                shapesView[sizeIndex] = internal::ST_Index_to_ShapeID(fragmentShape);
                colorView[sizeIndex] = fragment[1] - COLOR0;
                sizeIndex++;
              }
            }
          }

#if defined(AXOM_CLIP_FILTER_DEGENERATES)
          // Reduce overall whether there are degenerates.
          degenerates_reduce |= degenerates;
#endif
        });  // for_selected_zones

#if defined(AXOM_DEBUG_CLIP_FIELD)
      std::cout << "------------------------ makeTopology ------------------------" << std::endl;
      std::cout << "degenerates_reduce=" << degenerates_reduce.get() << std::endl;
      //      internal::printHost("selectedZones", selectedZones.view());
      internal::printHost("m_fragmentsView", fragmentData.m_fragmentsView);
      //      internal::printHost("zoneData.m_clipCasesView", zoneData.m_clipCasesView);
      //      internal::printHost("zoneData.m_pointsUsedView", zoneData.m_pointsUsedView);
      internal::printHost("conn", connView);
      internal::printHost("sizes", sizesView);
      internal::printHost("offsets", offsetsView);
      internal::printHost("shapes", shapesView);
      internal::printHost("color", colorView);
      std::cout << "--------------------------------------------------------------" << std::endl;
#endif
    }

#if defined(AXOM_CLIP_FILTER_DEGENERATES)
    // Filter out shapes that were marked as zero-size, adjusting connectivity and other arrays.
    if(degenerates_reduce.get())
    {
      internal::DegenerateHandler<TopologyView::dimension(), ExecSpace, ConnectivityType>::filterZeroSizes(
        fragmentData,
        n_sizes,
        n_offsets,
        n_shapes,
        n_color_values,
        sizesView,
        offsetsView,
        shapesView,
        colorView);
    }
#endif

    // Figure out which shapes were used.
    BitSet shapesUsed = findUsedShapes(shapesView);

#if defined(AXOM_DEBUG_CLIP_FIELD)
    std::cout << "------------------------ makeTopology ------------------------" << std::endl;
    internal::printHost("selectedZones", selectedZones.view());
    internal::printHost("m_fragmentsView", fragmentData.m_fragmentsView);
    internal::printHost("zoneData.m_clipCasesView", zoneData.m_clipCasesView);
    internal::printHost("zoneData.m_pointsUsedView", zoneData.m_pointsUsedView);
    internal::printHost("conn", connView);
    internal::printHost("sizes", sizesView);
    internal::printHost("offsets", offsetsView);
    internal::printHost("shapes", shapesView);
    internal::printHost("color", colorView);
    std::cout << "--------------------------------------------------------------" << std::endl;
#endif

    // If inside and outside are not selected, remove the color field since we should not need it.
    if(!(opts.inside() && opts.outside()))
    {
      n_newFields.remove(opts.colorField());
    }

#if defined(AXOM_CLIP_FILTER_DEGENERATES)
    // Handle some quad->tri degeneracies, depending on dimension.
    shapesUsed =
      internal::DegenerateHandler<TopologyView::dimension(), ExecSpace, ConnectivityType>::quadtri(
        shapesUsed,
        connView,
        sizesView,
        offsetsView,
        shapesView);
#endif

    // Add shape information to the connectivity.
    SLIC_ASSERT_MSG(shapesUsed != 0, "No shapes were produced!");
    const auto shapeMap = shapeMap_FromFlags(shapesUsed);
    SLIC_ASSERT_MSG(shapeMap.empty() == false, "The shape map is empty!");
    if(axom::utilities::popcount(static_cast<std::uint64_t>(shapesUsed)) > 1)
    {
      n_newTopo["elements/shape"] = "mixed";
      conduit::Node &n_shape_map = n_newTopo["elements/shape_map"];
      for(auto it = shapeMap.cbegin(); it != shapeMap.cend(); it++)
        n_shape_map[it->first] = it->second;
    }
    else
    {
      n_newTopo["elements"].remove("shapes");
      n_newTopo["elements/shape"] = shapeMap.begin()->first;
    }
  }

  /*!
   * \brief Find the shapes that were used.
   *
   * \param shapesView The view that contains the shapes.
   *
   * \return A BitSet where bits are marked for each shape used.
   */
  BitSet findUsedShapes(axom::ArrayView<ConnectivityType> shapesView) const
  {
    AXOM_ANNOTATE_SCOPE("findUsedShapes");

    RAJA::ReduceBitOr<reduce_policy, BitSet> shapesUsed_reduce(0);
    const axom::IndexType nShapes = shapesView.size();
    axom::for_all<ExecSpace>(
      nShapes,
      AXOM_LAMBDA(axom::IndexType index) {
        BitSet shapeBit = 1 << shapesView[index];
        shapesUsed_reduce |= shapeBit;
      });
    BitSet shapesUsed = shapesUsed_reduce.get();
    return shapesUsed;
  }

  /*!
   * \brief Make the new coordset using the blend data and the input coordset/coordsetview.
   *
   * \param blend The BlendData that we need to construct the new coordset.
   * \param n_coordset The input coordset, which is passed for metadata.
   * \param[out] n_newCoordset The new coordset.
   */
  void makeCoordset(const BlendData &blend,
                    const conduit::Node &n_coordset,
                    conduit::Node &n_newCoordset) const
  {
    AXOM_ANNOTATE_SCOPE("makeCoordset");
    // _mir_utilities_coordsetblender_begin
    axom::mir::utilities::blueprint::
      CoordsetBlender<ExecSpace, CoordsetView, axom::mir::utilities::blueprint::SelectSubsetPolicy>
        cb;
    cb.execute(blend, m_coordsetView, n_coordset, n_newCoordset);
    // _mir_utilities_coordsetblender_end
  }

  /*!
   * \brief Make new fields for the output topology.
   *
   * \param blend The BlendData that we need to construct the new vertex fields.
   * \param slice The SliceData we need to construct new element fields.
   * \param topologyName The name of the new field's topology.
   * \param fieldMap A map containing the names of the fields that we'll operate on.
   * \param n_fields The source fields.
   * \param[out] n_out_fields The node that will contain the new fields.
   */
  void makeFields(const BlendData &blend,
                  const SliceData &slice,
                  const std::string &topologyName,
                  const std::map<std::string, std::string> &fieldMap,
                  const conduit::Node &n_fields,
                  conduit::Node &n_out_fields) const
  {
    AXOM_ANNOTATE_SCOPE("makeFields");
    constexpr bool ss = axom::mir::views::view_traits<TopologyView>::supports_strided_structured();

    for(auto it = fieldMap.begin(); it != fieldMap.end(); it++)
    {
      const conduit::Node &n_field = n_fields.fetch_existing(it->first);
      const std::string association = n_field["association"].as_string();
      if(association == "element")
      {
        // Conditionally support strided-structured.
        bool handled =
          internal::StridedStructuredFields<ss, ExecSpace, TopologyView>::sliceElementField(
            m_topologyView,
            slice,
            n_field,
            n_out_fields[it->second]);

        if(!handled)
        {
          axom::mir::utilities::blueprint::FieldSlicer<ExecSpace> s;
          s.execute(slice, n_field, n_out_fields[it->second]);
        }

        n_out_fields[it->second]["topology"] = topologyName;
      }
      else if(association == "vertex")
      {
        // Conditionally support strided-structured.
        bool handled =
          internal::StridedStructuredFields<ss, ExecSpace, TopologyView>::blendVertexField(
            m_topologyView,
            blend,
            n_field,
            n_out_fields[it->second]);

        if(!handled)
        {
          // Blend the field normally.
          axom::mir::utilities::blueprint::FieldBlender<ExecSpace, axom::mir::utilities::blueprint::SelectSubsetPolicy>
            b;
          b.execute(blend, n_field, n_out_fields[it->second]);
        }

        n_out_fields[it->second]["topology"] = topologyName;
      }
    }
  }

  /*!
   * \brief Make an originalElements field so we can know each output zone's original zone number in the input mesh.
   *
   * \param[in] fragmentData This object holds views to per-fragment data.
   * \param[in] opts Clipping options.
   * \param[in] selectedZones The selected zones.
   * \param[in] n_fields The node that contains the input mesh's fields.
   * \param[out] n_newTopo The node that will contain the new topology.
   * \param[out] n_newFields The node that will contain the new fields.
   *
   * \note Objects that we need to capture into kernels are passed by value (they only contain views anyway). Data can be modified through the views.
   */
  void makeOriginalElements(FragmentData fragmentData,
                            const ClipOptions &opts,
                            const SelectedZones &selectedZones,
                            const conduit::Node &n_fields,
                            conduit::Node &n_newTopo,
                            conduit::Node &n_newFields) const
  {
    AXOM_ANNOTATE_SCOPE("makeOriginalElements");
    namespace bputils = axom::mir::utilities::blueprint;
    constexpr auto connTypeID = bputils::cpp2conduit<ConnectivityType>::id;

    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    const int conduitAllocatorID = c2a.getConduitAllocatorID();

    const auto selectedZonesView = selectedZones.view();
    const auto nzones = selectedZonesView.size();
    const std::string originalElements(opts.originalElementsField());

    if(n_fields.has_child(originalElements))
    {
      // originalElements already exists. We need to map it forward.
      const conduit::Node &n_orig = n_fields[originalElements];
      const conduit::Node &n_orig_values = n_orig["values"];
      views::IndexNode_to_ArrayView(n_orig_values, [&](auto origValuesView) {
        using value_type = typename decltype(origValuesView)::value_type;
        conduit::Node &n_origElem = n_newFields[originalElements];
        n_origElem["association"] = "element";
        n_origElem["topology"] = opts.topologyName(n_newTopo.name());
        conduit::Node &n_values = n_origElem["values"];
        n_values.set_allocator(conduitAllocatorID);
        n_values.set(conduit::DataType(n_orig_values.dtype().id(), fragmentData.m_finalNumZones));
        auto valuesView = bputils::make_array_view<value_type>(n_values);
        makeOriginalElements_copy(fragmentData, selectedZones, valuesView, origValuesView);
      });
    }
    else
    {
      // Make a new node and populate originalElement.
      conduit::Node &n_orig = n_newFields[originalElements];
      n_orig["association"] = "element";
      n_orig["topology"] = opts.topologyName(n_newTopo.name());
      conduit::Node &n_values = n_orig["values"];
      n_values.set_allocator(conduitAllocatorID);
      n_values.set(conduit::DataType(connTypeID, fragmentData.m_finalNumZones));
      auto valuesView = bputils::make_array_view<ConnectivityType>(n_values);
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          const int sizeIndex = fragmentData.m_fragmentOffsetsView[index];
          const int nFragments = fragmentData.m_fragmentsView[index];
          const auto zoneIndex = selectedZonesView[index];
          for(int i = 0; i < nFragments; i++) valuesView[sizeIndex + i] = zoneIndex;
        });
    }
  }

  /*!
   * \brief Assist setting original elements that already exist, based on selected zones.
   *
   * \param[in] fragmentData This object holds views to per-fragment data.
   * \param[in] selectedZones The selected zones.
   * \param[out] valuesView The destination values view.
   * \param[in] origValuesView The source values view.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename DataView>
  void makeOriginalElements_copy(FragmentData fragmentData,
                                 const SelectedZones &selectedZones,
                                 DataView valuesView,
                                 DataView origValuesView) const
  {
    const auto selectedZonesView = selectedZones.view();
    const auto nzones = selectedZonesView.size();
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType index) {
        const int sizeIndex = fragmentData.m_fragmentOffsetsView[index];
        const int nFragments = fragmentData.m_fragmentsView[index];
        const auto zoneIndex = selectedZonesView[index];
        for(int i = 0; i < nFragments; i++) valuesView[sizeIndex + i] = origValuesView[zoneIndex];
      });
  }

  /*!
   * \brief Given a flag that includes bitwise-or'd shape ids, make a map that indicates which Conduit shapes are used.
   *
   * \param shapes This is a bitwise-or of various (1 << ShapeID) values.
   *
   * \return A map of Conduit shape name to ShapeID value.
   */
  std::map<std::string, int> shapeMap_FromFlags(std::uint64_t shapes) const
  {
    std::map<std::string, int> sm;

    if(axom::utilities::bitIsSet(shapes, views::Line_ShapeID)) sm["line"] = views::Line_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Tri_ShapeID)) sm["tri"] = views::Tri_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Quad_ShapeID)) sm["quad"] = views::Quad_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Polygon_ShapeID))
      sm["polygon"] = views::Polygon_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Tet_ShapeID)) sm["tet"] = views::Tet_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Pyramid_ShapeID))
      sm["pyramid"] = views::Pyramid_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Wedge_ShapeID)) sm["wedge"] = views::Wedge_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Hex_ShapeID)) sm["hex"] = views::Hex_ShapeID;

    if(axom::utilities::bitIsSet(shapes, views::Polyhedron_ShapeID))
      sm["polyhedron"] = views::Polyhedron_ShapeID;

    return sm;
  }

  /*!
   * \brief If we're making a field that marks the new nodes that were created as
   *        a result of clipping, update those nodes now.
   *
   * \param blend The blend data used to create nodal fields.
   * \param newNodes The name of the new nodes field.
   * \param topoName The name of the output topology.
   * \param[inout] n_newFields The fields node for the output mesh.
   */
  void markNewNodes(const BlendData &blend,
                    const std::string &newNodes,
                    const std::string &topoName,
                    conduit::Node &n_newFields) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("markNewNodes");
    if(!newNodes.empty())
    {
      const auto origSize = blend.m_originalIdsView.size();
      const auto blendSize = blend.m_selectedIndicesView.size();
      const auto outputSize = origSize + blendSize;
      using Precision = int;
      constexpr Precision one = 1;
      constexpr Precision zero = 0;

      if(n_newFields.has_child(newNodes))
      {
        // Update the field. The field would have gone through field blending.
        // We can mark the new nodes with fresh values. This comes up in
        // applications that call Clip multiple times.

        conduit::Node &n_new_nodes = n_newFields.fetch_existing(newNodes);
        conduit::Node &n_new_nodes_values = n_new_nodes["values"];
        auto valuesView = bputils::make_array_view<Precision>(n_new_nodes_values);

        // Update values for the blend groups only.
        axom::for_all<ExecSpace>(
          blendSize,
          AXOM_LAMBDA(axom::IndexType bgid) { valuesView[origSize + bgid] = one; });
      }
      else
      {
        // Make the field for the first time.
        // Allocate Conduit data through Axom.
        bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;
        conduit::Node &n_new_nodes = n_newFields[newNodes];
        n_new_nodes["topology"] = topoName;
        n_new_nodes["association"] = "vertex";
        conduit::Node &n_new_nodes_values = n_new_nodes["values"];
        n_new_nodes_values.set_allocator(c2a.getConduitAllocatorID());
        n_new_nodes_values.set(conduit::DataType(bputils::cpp2conduit<Precision>::id, outputSize));
        auto valuesView = bputils::make_array_view<Precision>(n_new_nodes_values);

        // Fill in values. Everything below origSize is an original node.
        // Everything above is a blended node.
        axom::for_all<ExecSpace>(
          outputSize,
          AXOM_LAMBDA(axom::IndexType index) { valuesView[index] = (index < origSize) ? zero : one; });
      }
    }
  }

// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif
  TopologyView m_topologyView {};
  CoordsetView m_coordsetView {};
  Intersector m_intersector {};
  axom::mir::clipping::ClipTableManager<ExecSpace> m_clipTables {};
  NamingPolicy m_naming {};
};

}  // end namespace clipping
}  // end namespace mir
}  // end namespace axom

#endif
