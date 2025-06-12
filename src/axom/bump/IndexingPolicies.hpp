// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_INDEXING_POLICIES_HPP_
#define AXOM_BUMP_INDEXING_POLICIES_HPP_

#include "axom/bump/utilities/conduit_memory.hpp"

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
  void update(const conduit::Node &field)
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
    return m_indexing.LocalToGlobal(index);
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
  void update(const conduit::Node &field)
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
    const auto topoGlobalLogical = m_topoIndexing.GlobalToGlobal(index);
    // Make the global logical into a local logical in the topo.
    const auto topoLocalLogical = m_topoIndexing.GlobalToLocal(topoGlobalLogical);
    // Make the global logical index in the field.
    const auto fieldGlobalLogical = m_fieldIndexing.LocalToGlobal(topoLocalLogical);
    // Make the global index in the field.
    const auto fieldGlobalIndex = m_fieldIndexing.GlobalToGlobal(fieldGlobalLogical);
    return fieldGlobalIndex;
  }

  Indexing m_topoIndexing {};
  Indexing m_fieldIndexing {};
};

}  // end namespace bump
}  // end namespace axom

#endif
