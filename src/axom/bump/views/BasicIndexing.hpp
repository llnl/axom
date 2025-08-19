// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_BASIC_INDEXING_HPP_
#define AXOM_BUMP_BASIC_INDEXING_HPP_

namespace axom
{
namespace bump
{
namespace views
{
/*!
 * \brief A basic indexing class that provides methods for indexing transforms
 *        that have no effect. This is provided for topology view compatibility
 *        in algorithms that might use structured/unstructured topologies where
 *        index transforms are important.
 *        
 * \tparam IndexT The type used for indices.
 */
class BasicIndexing
{
public:
  /// Constructors
  ///@{
  AXOM_HOST_DEVICE
  BasicIndexing() : m_size(0) { }
  AXOM_HOST_DEVICE
  BasicIndexing(axom::IndexType size) : m_size(size) { }
  ///@}

  /*!
   * \brief Return the size
   *
   * \return Returns the size.
   */
  AXOM_HOST_DEVICE axom::IndexType size() const { return m_size; }

  /// No-effect indexing transforms

  AXOM_HOST_DEVICE inline axom::IndexType LocalToLocal(axom::IndexType index) const
  {
    return index;
  }

  AXOM_HOST_DEVICE inline axom::IndexType LocalToGlobal(axom::IndexType index) const
  {
    return index;
  }

  AXOM_HOST_DEVICE inline axom::IndexType GlobalToLocal(axom::IndexType index) const
  {
    return index;
  }

  AXOM_HOST_DEVICE inline axom::IndexType GlobalToGlobal(axom::IndexType index) const
  {
    return index;
  }

  /*!
   * \brief Determines whether the indexing contains the supplied index.
   *
   * \param index The index being tested.
   *
   * \return True if the index is within the index, false otherwise.
   */
  AXOM_HOST_DEVICE inline bool contains(axom::IndexType index) const
  {
    return index >= 0 && index < m_size;
  }

  AXOM_HOST_DEVICE inline axom::IndexType clamp(axom::IndexType index) const
  {
    return axom::utilities::clampVal(index, axom::IndexType(0), axom::IndexType(m_size - 1));
  }

  axom::IndexType m_size {0};
};

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
