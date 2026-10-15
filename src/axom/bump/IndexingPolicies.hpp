// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/view_traits.hpp"
#include "axom/fmt.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{

//------------------------------------------------------------------------------
/*!
 * \brief Returns the input index (no changes).
 */
struct DirectIndexing
{
  /*!
   * \brief Return the input index (no changes).
   * \param index The input index.
   * \return The input index.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType index) const { return index; }
};

namespace detail
{
/*!
 * \brief Apply Blueprint offsets and strides to a strided indexing.
 *
 * This preserves the indexing's logical dimensions. A strided layout requires \c strides,
 * while missing \c offsets default to zero.
 *
 * \tparam Indexing A StridedStructuredIndexing of some dimension.
 *
 * \param node A field or topology-dimensions node that contains \c strides and may contain \c offsets.
 *        It can be in device memory.
 * \param[in,out] indexing The indexing to update.
 *
 * \note Executes on the host. Raises SLIC_ERROR if a layout has \c offsets without \c strides,
 *       if metadata does not have one entry per dimension, if an offset is negative,
 *       or if a stride is not positive.
 */
template <typename Indexing>
void setLayoutFromNode(const conduit::Node& node, Indexing& indexing)
{
  using LogicalIndex = typename Indexing::LogicalIndex;
  constexpr int NDIMS = Indexing::dimension();

  const Indexing compact(indexing.logicalDimensions());
  LogicalIndex offsets = compact.m_offsets;
  LogicalIndex strides = compact.m_strides;

  SLIC_ERROR_IF(node.has_path("offsets") && !node.has_path("strides"),
                axom::fmt::format("Blueprint layout '{}' has offsets but no strides. "
                                  "Strided-structured layouts require strides.",
                                  node.path()));

  for(const char* key : {"offsets", "strides"})
  {
    if(node.has_path(key))
    {
      const auto count = node.fetch_existing(key).dtype().number_of_elements();
      SLIC_ERROR_IF(count != NDIMS,
                    axom::fmt::format("'{}' has {} {} entries, but a {}D topology requires {}.",
                                      node.path(),
                                      count,
                                      key,
                                      NDIMS,
                                      NDIMS));
    }
  }
  axom::bump::utilities::fillFromNode(node, "offsets", offsets, true);
  axom::bump::utilities::fillFromNode(node, "strides", strides, true);
  for(int d = 0; d < NDIMS; d++)
  {
    SLIC_ERROR_IF(offsets[d] < 0 || strides[d] < 1,
                  axom::fmt::format("'{}' has offset {} and stride {} in dimension {}. "
                                    "Offsets must be nonnegative and strides must be positive.",
                                    node.path(),
                                    offsets[d],
                                    strides[d],
                                    d));
  }
  indexing = Indexing(compact.logicalDimensions(), offsets, strides);
}
}  // end namespace detail

//------------------------------------------------------------------------------
/*!
 * \brief Map compact zone indices to a strided-structured element field.
 *
 * Bump assigns compact indices from 0 to numberOfZones() - 1. A field on a strided-structured
 * topology may have padding. Blueprint stores that field's \c offsets and \c strides separately
 * from the topology's \c elements/dims metadata.
 *
 * The object is trivially copyable and can be captured by value in device kernels.
 *
 * \tparam Indexing A StridedStructuredIndexing of some dimension.
 */
template <typename Indexing>
struct SSElementFieldIndexing
{
  SSElementFieldIndexing() = default;

  /*!
   * \brief Create an index map from the topology's zone indexing and an element field.
   *
   * \param zoneIndexing The topology's zone indexing. It supplies the zone dimensions.
   * \param field The Conduit node for the element field.
   *
   * \note Executes on the host.
   */
  SSElementFieldIndexing(const Indexing& zoneIndexing, const conduit::Node& field)
    : m_indexing(zoneIndexing)
  {
    update(field);
  }

  /*!
   * \brief Update the field layout.
   *
   * \param field The Conduit node for an element field.
   *
   * \note Preserves the zone dimensions in \c m_indexing. Executes on the host.
   */
  void update(const conduit::Node& field) { detail::setLayoutFromNode(field, m_indexing); }

  /*!
   * \brief Transforms the index from local to global through an indexing object.
   * \param index The local index
   * \return The global index for the field.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType index) const
  {
    return m_indexing.localToGlobal(index);
  }

  Indexing m_indexing {};
};

//------------------------------------------------------------------------------
/*!
 * \brief Map topology node indices to a strided-structured vertex field.
 *
 * Topology node indices address the coordset layout in \c elements/dims.
 * A vertex field can have its own \c offsets and \c strides.
 *
 * The object is trivially copyable and can be captured by value in device kernels.
 *
 * \tparam Indexing A StridedStructuredIndexing of some dimension.
 */
template <typename Indexing>
struct SSVertexFieldIndexing
{
  SSVertexFieldIndexing() = default;

  /*!
   * \brief Create an index map from the topology's node indexing and a vertex field.
   *
   * \param nodeIndexing The topology's node indexing, such as zone indexing's expand().
   * \param field The Conduit node for the vertex field.
   *
   * \note Executes on the host.
   */
  SSVertexFieldIndexing(const Indexing& nodeIndexing, const conduit::Node& field)
    : m_topoIndexing(nodeIndexing)
    , m_fieldIndexing(nodeIndexing)
  {
    update(field);
  }

  /*!
   * \brief Update the field layout.
   *
   * \param field The Conduit node for a vertex field.
   *
   * \note Preserves the node dimensions in \c m_fieldIndexing. Executes on the host.
   */
  void update(const conduit::Node& field) { detail::setLayoutFromNode(field, m_fieldIndexing); }

  /*!
   * \brief Transforms the index from local to global through an indexing object.
   * \param index The global index
   * \return The global index for the field.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType index) const
  {
    // Make the global index into a global logical in the topo.
    const auto topoGlobalLogical = m_topoIndexing.globalToGlobal(index);
    // Make the global logical into a local logical in the topo.
    const auto topoLocalLogical = m_topoIndexing.globalToLocal(topoGlobalLogical);
    // Make the global logical index in the field.
    const auto fieldGlobalLogical = m_fieldIndexing.localToGlobal(topoLocalLogical);
    // Make the global index in the field.
    const auto fieldGlobalIndex = m_fieldIndexing.globalToGlobal(fieldGlobalLogical);
    return fieldGlobalIndex;
  }

  Indexing m_topoIndexing {};
  Indexing m_fieldIndexing {};
};

//------------------------------------------------------------------------------
/*!
 * \brief Make an index map from compact zone indices to an element field.
 *
 * Element fields normally store one value per zone in compact zone order.
 * On a strided-structured topology, a field's Blueprint \c offsets and \c strides can place those
 * values in a larger array. This function maps a compact zone index to that array position.
 *
 * \param topologyView The topology view.
 * \param n_field The Conduit node for the element-associated field.
 *
 * \return An SSElementFieldIndexing that reads the field layout for strided-structured views,
 *         or DirectIndexing for fields stored in compact zone order. Both are trivially copyable
 *         and can be captured by value in device kernels.
 *
 * \note Executes on the host. Raises SLIC_ERROR if the field is not element-associated, if it
 *       has \c offsets or \c strides on a view that is not strided-structured, if its layout
 *       metadata is invalid, or if any zone maps outside the field's values.
 */
template <typename TopologyView>
auto makeElementFieldIndexing(const TopologyView& topologyView, const conduit::Node& n_field)
{
  SLIC_ERROR_IF(!n_field.has_path("association") ||
                  n_field.fetch_existing("association").as_string() != "element",
                axom::fmt::format("Field '{}' is not element-associated.", n_field.name()));

  const auto indexing = [&]() {
    if constexpr(views::view_traits<TopologyView>::supports_strided_structured())
    {
      return SSElementFieldIndexing<typename TopologyView::IndexingPolicy>(topologyView.indexing(),
                                                                           n_field);
    }
    else
    {
      SLIC_ERROR_IF(n_field.has_path("offsets") || n_field.has_path("strides"),
                    axom::fmt::format("Field '{}' has Blueprint offsets or strides, but its "
                                      "topology view does not support strided-structured fields.",
                                      n_field.name()));
      return DirectIndexing {};
    }
  }();

  // Strides are positive and offsets nonnegative,
  // so the first and last zones bound the range of field indices.
  const axom::IndexType nzones = topologyView.numberOfZones();
  if(nzones > 0)
  {
    const conduit::Node& n_values = n_field.fetch_existing("values");
    const conduit::Node& n_component =
      (n_values.number_of_children() > 0) ? n_values.child(0) : n_values;
    const axom::IndexType nvalues = n_component.dtype().number_of_elements();
    const axom::IndexType first = indexing[0];
    const axom::IndexType last = indexing[nzones - 1];
    SLIC_ERROR_IF(
      first < 0 || last >= nvalues,
      axom::fmt::format("Field '{}' has {} values, but its layout maps zones to indices [{}, {}].",
                        n_field.name(),
                        nvalues,
                        first,
                        last));
  }
  return indexing;
}

}  // end namespace bump
}  // end namespace axom
