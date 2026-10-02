// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/StridedStructuredIndexing.hpp"
#include "axom/bump/views/StructuredTopologyView.hpp"
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

//------------------------------------------------------------------------------
/*!
 * \brief Help turn slice data zone indices into strided structured element field indices.
 * \tparam Indexing A StridedStructuredIndexing of some dimension.
 */
template <typename Indexing>
struct SSElementFieldIndexing
{
  /*!
   * \brief Update the indexing offsets/strides from a Conduit node.
   * \param field The Conduit node for a field.
   *
   * \note Executes on the host.
   */
  void update(const conduit::Node& field)
  {
    axom::bump::utilities::fillFromNode(field, "offsets", m_indexing.m_offsets, true);
    axom::bump::utilities::fillFromNode(field, "strides", m_indexing.m_strides, true);
  }

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
 * \brief Help turn blend group node indices (global) into vertex field indices.
 * \tparam Indexing A StridedStructuredIndexing of some dimension.
 */
template <typename Indexing>
struct SSVertexFieldIndexing
{
  /*!
   * \brief Update the indexing offsets/strides from a Conduit node.
   * \param field The Conduit node for a field.
   *
   * \note Executes on the host.
   */
  void update(const conduit::Node& field)
  {
    axom::bump::utilities::fillFromNode(field, "offsets", m_fieldIndexing.m_offsets, true);
    axom::bump::utilities::fillFromNode(field, "strides", m_fieldIndexing.m_strides, true);
  }

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
 * \brief Map compact zone indices to indices in an element-associated field.
 *
 * Bump numbers zones compactly, from 0 to numberOfZones() - 1.
 * An element field on a strided-structured topology can include padding.
 * Since Blueprint stores its \c offsets and \c strides on the field
 * separately from the topology's \c elements/dims metadata,
 * the mapping depends on the field.
 *
 * This template covers views whose element fields are stored in zone order.
 * Its update() rejects fields with \c offsets or \c strides,
 * since that would require strided zone indexing.
 *
 * Use makeElementFieldIndexing() to construct this type. It is trivially
 * copyable and can be captured by value in device kernels.
 *
 * \tparam TopologyView The topology view type.
 */
template <typename TopologyView>
struct ElementFieldIndexing
{
  /*!
   * \brief Initialize the indexing from a field.
   *
   * \param topologyView (unused) The topology view.
   * \param n_field The Conduit node for the element-associated field.
   *
   * \note Executes on the host.
   */
  void update(const TopologyView& AXOM_UNUSED_PARAM(topologyView), const conduit::Node& n_field)
  {
    SLIC_ERROR_IF(n_field.has_path("offsets") || n_field.has_path("strides"),
                  axom::fmt::format("Field '{}' has Blueprint offsets or strides, but its topology "
                                    "view does not support strided-structured element fields.",
                                    n_field.name()));
  }

  /*!
   * \brief Return the field index for a zone.
   * \param zoneIndex The compact zone index.
   * \return \a zoneIndex because the field is stored in zone order.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType zoneIndex) const { return zoneIndex; }
};

/*!
 * \brief Map compact zone indices to indices in an element field on a strided-structured topology.
 *
 * The field's \c offsets and \c strides determine the layout.
 * Missing \c offsets default to zero, and missing \c strides default to the compact zone layout.
 * A field with neither is indexed in zone order, and a field with both matches SSElementFieldIndexing.
 *
 * \tparam IndexT The index type of the strided indexing.
 * \tparam NDIMS  The number of topological dimensions.
 */
template <typename IndexT, int NDIMS>
struct ElementFieldIndexing<views::StructuredTopologyView<views::StridedStructuredIndexing<IndexT, NDIMS>>>
{
  using Indexing = views::StridedStructuredIndexing<IndexT, NDIMS>;
  using TopologyView = views::StructuredTopologyView<Indexing>;
  using LogicalIndex = typename Indexing::LogicalIndex;

  /*!
   * \brief Initialize the indexing from the topology's zone dimensions and the field's layout metadata.
   *
   * \param topologyView The topology view.
   * \param n_field The Conduit node for the element-associated field.
   *
   * \note Executes on the host. The metadata can be in device memory.
   */
  void update(const TopologyView& topologyView, const conduit::Node& n_field)
  {
    const LogicalIndex zoneDims = topologyView.indexing().logicalDimensions();
    LogicalIndex offsets, strides;
    IndexT stride = 1;
    for(int d = 0; d < NDIMS; d++)
    {
      offsets[d] = 0;
      strides[d] = stride;
      stride *= zoneDims[d];
    }

    auto validate_layout_size = [&](const char* key) {
      if(n_field.has_path(key))
      {
        const conduit::Node& n_layout = n_field.fetch_existing(key);
        SLIC_ERROR_IF(
          n_layout.dtype().number_of_elements() != NDIMS,
          axom::fmt::format("Field '{}' has {} {} entries, but a {}D topology requires {}.",
                            n_field.name(),
                            n_layout.dtype().number_of_elements(),
                            key,
                            NDIMS,
                            NDIMS));
      }
    };
    validate_layout_size("offsets");
    validate_layout_size("strides");

    axom::bump::utilities::fillFromNode(n_field, "offsets", offsets, true);
    axom::bump::utilities::fillFromNode(n_field, "strides", strides, true);
    for(int d = 0; d < NDIMS; d++)
    {
      SLIC_ERROR_IF(offsets[d] < 0 || strides[d] < 1,
                    axom::fmt::format("Field '{}' has offset {} and stride {} in dimension {}. "
                                      "Offsets must be nonnegative and strides must be positive.",
                                      n_field.name(),
                                      offsets[d],
                                      strides[d],
                                      d));
    }
    m_indexing = Indexing(zoneDims, offsets, strides);
  }

  /*!
   * \brief Return the field index for a zone.
   * \param zoneIndex The compact zone index.
   * \return The index of \a zoneIndex in the field's values.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType zoneIndex) const
  {
    return m_indexing.localToGlobal(zoneIndex);
  }

  Indexing m_indexing {};
};

/*!
 * \brief Make an ElementFieldIndexing for a topology view and an element field.
 *
 * \param topologyView The topology view.
 * \param n_field The Conduit node for the element-associated field.
 *
 * \return An indexing object that maps compact zone indices to field indices.
 *
 * \note Executes on the host. Raises SLIC_ERROR if the field is not element-associated,
 *       if its layout cannot be honored for this view, or if any zone maps outside the field's values.
 */
template <typename TopologyView>
ElementFieldIndexing<TopologyView> makeElementFieldIndexing(const TopologyView& topologyView,
                                                            const conduit::Node& n_field)
{
  SLIC_ERROR_IF(!n_field.has_path("association") ||
                  n_field.fetch_existing("association").as_string() != "element",
                axom::fmt::format("Field '{}' is not element-associated.", n_field.name()));

  ElementFieldIndexing<TopologyView> indexing;
  indexing.update(topologyView, n_field);

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
