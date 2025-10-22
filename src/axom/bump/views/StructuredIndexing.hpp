// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_STRUCTURED_INDEXING_HPP_
#define AXOM_BUMP_STRUCTURED_INDEXING_HPP_

#include "axom/core/StackArray.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/utilities/Utilities.hpp"

namespace axom
{
namespace bump
{
namespace views
{
/*!
 * \brief This class encapsulates a structured mesh size and contains methods to
 *        help with indexing into it.
 *
 * \tparam IndexT The type used for indices.
 * \tparam NDIMS The number of dimensions.
 */
template <typename IndexT, int NDIMS = 3>
class StructuredIndexing
{
public:
  using IndexType = IndexT;
  using LogicalIndex = axom::StackArray<IndexType, NDIMS>;

  AXOM_HOST_DEVICE constexpr static int dimension() { return NDIMS; }

  /*!
   * \brief Return whether the view supports strided structured indexing.
   * \return false
   */
  AXOM_HOST_DEVICE static constexpr bool supports_strided_structured_indexing() { return false; }

  /*!
   * \brief constructor
   *
   * \param dims The dimensions we're indexing.
   */
  AXOM_HOST_DEVICE
  StructuredIndexing() : m_dimensions()
  {
    for(int d = 0; d < NDIMS; d++)
    {
      m_dimensions[d] = 1;
    }
  }

  AXOM_HOST_DEVICE
  StructuredIndexing(const LogicalIndex &dims) : m_dimensions(dims)
  {
#if !defined(AXOM_DEVICE_CODE)
    for(int d = 0; d < NDIMS; d++)
    {
      SLIC_ERROR_IF(m_dimensions[d] < 1, "Dimensions must greater than or equal to 1.");
    }
#endif
  }

  /*!
   * \brief Return the number of points in the index space.
   *
   * \return The number of points in the index space.
   */
  AXOM_HOST_DEVICE
  IndexType size() const
  {
    IndexType sz = 1;
    for(int i = 0; i < NDIMS; i++)
    {
      sz *= m_dimensions[i];
    }
    return sz;
  }

  /*!
   * \brief Return the logical dimensions.
   *
   * \return The logical dimensions.
   */
  AXOM_HOST_DEVICE
  const LogicalIndex &logicalDimensions() const { return m_dimensions; }

  /*!
   * \brief Return the j stride.
   *
   * \return The j stride to move up a row.
   */
  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims >= 2, IndexType>::type jStride() const
  {
    return m_dimensions[0];
  }

  /*!
   * \brief Return the k stride.
   *
   * \return The k stride to move forward a "page".
   */
  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 3, IndexType>::type kStride() const
  {
    return m_dimensions[0] * m_dimensions[1];
  }

  /*!
   * \brief Turn a global logical index into an index.
   * \param global The global logical index to convert.
   * \return The global index.
   */
  AXOM_HOST_DEVICE
  inline IndexType GlobalToGlobal(const LogicalIndex &global) const
  {
    return LogicalIndexToIndex(global);
  }

  /*!
   * \brief Turn a global index into a global logical index.
   * \param global The global index to convert.
   * \return The global logical index.
   */
  AXOM_HOST_DEVICE
  inline LogicalIndex GlobalToGlobal(IndexType global) const { return IndexToLogicalIndex(global); }

  /*!
   * \brief Turn global logical index to local logical index. no-op.
   * \param index The index to convert.
   * \return Same as the input in this case.
   */
  AXOM_HOST_DEVICE
  inline LogicalIndex GlobalToLocal(const LogicalIndex &index) const { return index; }

  /*!
   * \brief Turn global index to local index. no-op.
   * \param index The index to convert.
   * \return Same as the input in this case.
   */
  AXOM_HOST_DEVICE
  inline IndexType GlobalToLocal(IndexType index) const { return index; }

  /*!
   * \brief Turn local logical index to global logical index. no-op.
   * \param index The index to convert.
   * \return Same as the input in this case.
   */
  AXOM_HOST_DEVICE
  inline LogicalIndex LocalToGlobal(const LogicalIndex &index) const { return index; }

  /*!
   * \brief Turn local index to global index. no-op.
   * \param index The index to convert.
   * \return Same as the input in this case.
   */
  AXOM_HOST_DEVICE
  inline IndexType LocalToGlobal(IndexType index) const { return index; }

  /*!
   * \brief Turn an index into a logical index.
   *
   * \param index The index to convert.
   *
   * \return The logical index that corresponds to the \a index.
   */
  /// @{

  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 1, LogicalIndex>::type IndexToLogicalIndex(
    IndexType index) const
  {
    LogicalIndex logical;
    logical[0] = index;
    return logical;
  }

  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 2, LogicalIndex>::type IndexToLogicalIndex(
    IndexType index) const
  {
    LogicalIndex logical;
    const auto nx = m_dimensions[0];
    logical[0] = index % nx;
    logical[1] = index / nx;
    return logical;
  }

  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 3, LogicalIndex>::type IndexToLogicalIndex(
    IndexType index) const
  {
    LogicalIndex logical;
    const auto nx = m_dimensions[0];
    const auto nxy = nx * m_dimensions[1];
    logical[0] = index % nx;
    logical[1] = (index % nxy) / nx;
    logical[2] = index / nxy;
    return logical;
  }

  /// @}

  /*!
   * \brief Turn a logical index into a flat index.
   *
   * \param logical The logical indexto convert to a flat index.
   *
   * \return The index that corresponds to the \a logical index.
   */
  /// @{
  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 1, IndexType>::type LogicalIndexToIndex(
    const LogicalIndex &logical) const
  {
    return logical[0];
  }

  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 2, IndexType>::type LogicalIndexToIndex(
    const LogicalIndex &logical) const
  {
    return logical[1] * m_dimensions[0] + logical[0];
  }

  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 3, IndexType>::type LogicalIndexToIndex(
    const LogicalIndex &logical) const
  {
    return (logical[2] * m_dimensions[1] * m_dimensions[0]) + (logical[1] * m_dimensions[0]) +
      logical[0];
  }

  /// @}

  /*!
   * \brief Determines whether the indexing contains the supplied logical index.
   *
   * \param logical The logical index being tested.
   *
   * \return True if the logical index is within the index, false otherwise.
   */
  AXOM_HOST_DEVICE
  bool contains(const LogicalIndex &logical) const
  {
    bool retval = true;
    for(int i = 0; i < dimension(); i++)
    {
      retval &= (logical[i] >= 0 && logical[i] < m_dimensions[i]);
    }
    return retval;
  }

  /*!
   * \brief Determines whether the indexing contains the supplied index.
   *
   * \param index The index being tested.
   *
   * \return True if the index is within the index, false otherwise.
   */
  AXOM_HOST_DEVICE
  bool contains(IndexType index) const { return contains(IndexToLogicalIndex(index)); }

  /*!
   * \brief Expand the current StructuredIndexing by one in each dimension.
   *
   * \return An expanded StructuredIndexing.
   */
  AXOM_HOST_DEVICE
  StructuredIndexing expand() const
  {
    StructuredIndexing retval(*this);
    for(int i = 0; i < dimension(); i++)
    {
      retval.m_dimensions[i]++;
    }
    return retval;
  }

  /*!
   * \brief Given a logical index that may or may not be in the index space,
   *        return a clamped logical index that is in the index space. The
   *        logical indices are clamped to [0, dimensions-1] in each dimension.
   *
   * \param logical The input logical index being clamped.
   *
   * \retval A new LogicalIndex that is in the index space.
   */
  /// @{
  AXOM_HOST_DEVICE
  LogicalIndex clamp(const LogicalIndex &logical) const
  {
    LogicalIndex retval;
    const IndexType lower(0);
    for(int i = 0; i < dimension(); i++)
    {
      const IndexType upper(m_dimensions[i] - 1);
      retval[i] = axom::utilities::clampVal(logical[i], lower, upper);
    }
    return retval;
  }

  IndexType clamp(IndexType index) const
  {
    return LogicalIndexToIndex(clamp(IndexToLogicalIndex(index)));
  }
  /// @}

private:
  LogicalIndex m_dimensions {};
};

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
