// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_STRIDED_STRUCTURED_INDEXING_HPP_
#define AXOM_BUMP_STRIDED_STRUCTURED_INDEXING_HPP_

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
 * \brief This class encapsulates data for strided structured indexing and provides
 *        methods for creating/manipulating indices.
 *
 * \tparam IndexT The type used for indices.
 * \tparam NDIMS The number of dimensions.
 *
 *
 * Strided structured indexing lets us index part of a larger overall indexing space.
 * We can index it using indices local to the selected window and the class provides a
 * mechanism to return indices in the overall indexing space.
 *
 * Example:
 *
 *   x---x---x---x---x---x---x
 *   |   |   |   |   |   |   |
 *   x---x---*---*---*---*---x   *=real node, x=ignored node, O=origin node 16
 *   |   |   |   |   |   |   |
 *   x---x---*---*---*---*---x   dims={4,3}
 *   |   |   |   |   |   |   |   origin={2,2}
 *   x---x---O---*---*---*---x   stride={1,7}
 *   |   |   |   |   |   |   |
 *   x---x---x---x---x---x---x
 *   |   |   |   |   |   |   |
 *   x---x---x---x---x---x---x
 *
 */
template <typename IndexT, int NDIMS = 3>
struct StridedStructuredIndexing
{
  using IndexType = IndexT;
  using LogicalIndex = axom::StackArray<IndexType, NDIMS>;

  AXOM_HOST_DEVICE constexpr static int dimension() { return NDIMS; }

  /*!
   * \brief Return whether the view supports strided structured indexing.
   * \return true
   */
  AXOM_HOST_DEVICE static constexpr bool supports_strided_structured_indexing() { return true; }

  /*!
   * \brief constructor
   */
  AXOM_HOST_DEVICE
  StridedStructuredIndexing() : m_dimensions(), m_offsets(), m_strides()
  {
    for(int i = 0; i < NDIMS; i++)
    {
      m_dimensions[i] = 1;
      m_offsets[i] = 0;
      m_strides[i] = 1;
    }
  }

  /*!
   * \brief Constructor
   *
   * \param dims The number of zones in each logical dimension.
   * \param offsets The offset of the first zone from the mesh origin in each logical dimension.
   * \param strides The amount to stride when moving to the next element for each logical dimension.
   */
  AXOM_HOST_DEVICE
  StridedStructuredIndexing(const LogicalIndex &dims,
                            const LogicalIndex &offsets,
                            const LogicalIndex &strides)
    : m_dimensions(dims)
    , m_offsets(offsets)
    , m_strides(strides)
  {
#if !defined(AXOM_DEVICE_CODE)
    for(int d = 0; d < NDIMS; d++)
    {
      SLIC_ERROR_IF(m_dimensions[d] < 1, "Dimensions must greater than or equal to 1.");
    }
#endif
  }

  /*!
   * \brief Return the number of values in the index space.
   *
   * \return The number of values in the index space.
   */
  AXOM_HOST_DEVICE
  IndexType size() const
  {
    IndexType sz = 1;
    for(int i = 0; i < NDIMS; i++) sz *= m_dimensions[i];
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
  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims >= 2, IndexType>::type jStride() const
  {
    return m_strides[1];
  }

  /*!
   * \brief Return the k stride.
   *
   * \return The k stride to move forward a "page".
   */
  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 3, IndexType>::type kStride() const
  {
    return m_strides[2];
  }

  /*!
   * \brief Turn a global logical index into an index.
   * \param global The global logical index to convert.
   * \return The global index.
   */
  AXOM_HOST_DEVICE
  IndexType GlobalToGlobal(const LogicalIndex &global) const
  {
    IndexType gl {};
    for(int i = 0; i < NDIMS; i++)
    {
      gl += global[i] * m_strides[i];
    }
    return gl;
  }

  /*!
   * \brief Turn a global index into a global logical index.
   *
   * \param global The index to convert.
   *
   * \return The local index that corresponds to the \a local.
   */
  /// @{
  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 1, LogicalIndex>::type GlobalToGlobal(
    IndexType global) const
  {
    LogicalIndex gl;
    gl[0] = global;
    return global;
  }

  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 2, LogicalIndex>::type GlobalToGlobal(
    IndexType global) const
  {
    LogicalIndex gl;
    gl[0] = global % m_strides[1];
    gl[1] = global / m_strides[1];
    return gl;
  }

  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 3, LogicalIndex>::type GlobalToGlobal(
    IndexType global) const
  {
    LogicalIndex gl;
    gl[0] = global % m_strides[1];
    gl[1] = (global % m_strides[2]) / m_strides[1];
    gl[2] = global / m_strides[2];
    return gl;
  }
  /// @}

  /*!
   * \brief Convert global logical index to a local one.
   * \param local The local logical index.
   * \return local logical index.
   */
  AXOM_HOST_DEVICE
  LogicalIndex GlobalToLocal(const LogicalIndex &global) const
  {
    LogicalIndex local(global);
    for(int i = 0; i < NDIMS; i++)
    {
      local[i] -= m_offsets[i];
    }
    return local;
  }

  /*!
   * \brief Turn a global index into a local index.
   *
   * \param global The index to convert.
   *
   * \return The local index that corresponds to the \a local.
   */
  AXOM_HOST_DEVICE
  IndexType GlobalToLocal(IndexType global) const
  {
    return LogicalIndexToIndex(GlobalToLocal(GlobalToGlobal(global)));
  }

  /*!
   * \brief Convert local logical index to a global one.
   * \param local The local logical index.
   * \return global logical index.
   */
  AXOM_HOST_DEVICE
  LogicalIndex LocalToGlobal(const LogicalIndex &local) const
  {
    LogicalIndex global(local);
    for(int i = 0; i < NDIMS; i++)
    {
      global[i] += m_offsets[i];
    }
    return global;
  }

  /*!
   * \brief Convert local logical index to a global one.
   * \param local The local logical index.
   * \return local logical index.
   */
  AXOM_HOST_DEVICE
  IndexType LocalToGlobal(IndexType local) const
  {
    return GlobalToGlobal(LocalToGlobal(IndexToLogicalIndex(local)));
  }

  /*!
   * \brief Turn a local index into a local logical index.
   *
   * \param index The index to convert.
   *
   * \return The local logical index that corresponds to the \a index.
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
   * \brief Turn a local logical index into a local flat index.
   *
   * \param logical The logical indexto convert to a flat index.
   *
   * \return The index that corresponds to the \a logical index.
   */
  AXOM_HOST_DEVICE
  IndexType LogicalIndexToIndex(const LogicalIndex &logical) const
  {
    IndexType index {};
    IndexType stride {1};
    for(int i = 0; i < NDIMS; i++)
    {
      index += logical[i] * stride;
      stride *= m_dimensions[i];
    }
    return index;
  }

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
   * \brief Expand the current StridedStructuredIndexing by one in each dimension.
   *
   * \return An expanded StridedStructuredIndexing.
   */
  /// @{
  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 1, StridedStructuredIndexing>::type expand() const
  {
    StridedStructuredIndexing retval(*this);
    retval.m_dimensions[0]++;
    return retval;
  }

  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 2, StridedStructuredIndexing>::type expand() const
  {
    StridedStructuredIndexing retval(*this);
    retval.m_dimensions[0]++;
    retval.m_dimensions[1]++;
    retval.m_strides[1]++;
    return retval;
  }

  template <int _ndims = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 3, StridedStructuredIndexing>::type expand() const
  {
    StridedStructuredIndexing retval(*this);
    retval.m_dimensions[0]++;
    retval.m_dimensions[1]++;
    retval.m_dimensions[2]++;
    const auto nx = retval.m_strides[1];
    const auto ny = retval.m_strides[2] / nx;
    retval.m_strides[1] = nx + 1;
    retval.m_strides[2] = (ny + 1) * (nx + 1);
    return retval;
  }

  /// @}

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

  LogicalIndex m_dimensions {};
  LogicalIndex m_offsets {};
  LogicalIndex m_strides {};
};

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
