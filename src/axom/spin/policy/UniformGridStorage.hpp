// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_POLICY_UGRID_STORAGE_HPP
#define AXOM_SPIN_POLICY_UGRID_STORAGE_HPP

#include "axom/core/Array.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/reductions.hpp"
#include "axom/core/execution/scans.hpp"

namespace axom
{
namespace spin
{
namespace policy
{
template <typename T>
struct DynamicGridStorage;

template <typename T>
struct DynamicGridView;

template <typename T>
struct FlatGridStorage;

template <typename T>
struct FlatGridView;

/*!
 * \brief Policy for storing uniform grid as an array of arrays. Supports
 *  efficient single-element insertions and bin clears.
 */
template <typename T>
struct DynamicGridStorage
{
  using ViewType = DynamicGridView<T>;
  using ConstViewType = DynamicGridView<const T>;

  using BinType = axom::Array<T>;
  using ConstBinType = const axom::Array<T>;

  using BinRef = axom::Array<T>&;
  using ConstBinRef = const axom::Array<T>&;

  DynamicGridStorage(int allocID) : m_bins(0, 0, allocID), m_allocatorID(allocID)
  {
#ifdef AXOM_USE_UMPIRE
    SLIC_ASSERT_MSG(axom::detail::getAllocatorSpace(allocID) != MemorySpace::Device,
                    "Umpire device-only memory is not a supported allocator type for dynamic "
                    "UniformGrid storage.");
#endif
  }

  int getAllocatorID() const { return m_allocatorID; }

  void setNumBins(IndexType nbins) { m_bins.resize(nbins); }
  IndexType getNumBins() const { return m_bins.size(); }

  bool isValidIndex(IndexType index) const { return index >= 0 && index < getNumBins(); }

  template <typename ExecSpace>
  void initialize(const axom::ArrayView<const IndexType> binSizes)
  {
    m_bins.clear();
    for(int i = 0; i < binSizes.size(); i++)
    {
      m_bins.emplace_back(binSizes[i], binSizes[i], m_allocatorID);
    }
  }

  void insert(IndexType gridIdx, const T& elem) { m_bins[gridIdx].push_back(elem); }

  void clear(IndexType gridIdx)
  {
    if(isValidIndex(gridIdx))
    {
      m_bins[gridIdx].clear();
    }
  }

  // getters
  BinRef getBinContents(IndexType gridIdx)
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx];
  }
  ConstBinRef getBinContents(IndexType gridIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx];
  }

  T& get(IndexType gridIdx, IndexType binIdx)
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx][binIdx];
  }

  const T& get(IndexType gridIdx, IndexType binIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx][binIdx];
  }

  axom::Array<axom::Array<T>> m_bins;
  int m_allocatorID;
};

template <typename T>
struct DynamicGridView
{
  using BaseT = typename std::remove_const<T>::type;

  using BinStoreType =
    typename std::conditional<std::is_const<T>::value, const axom::Array<BaseT>, axom::Array<BaseT>>::type;

  DynamicGridView(DynamicGridStorage<BaseT>& in) : m_bins(in.m_bins.view()) { }

  DynamicGridView(const DynamicGridStorage<BaseT>& in) : m_bins(in.m_bins.view()) { }

  AXOM_HOST_DEVICE IndexType getNumBins() const { return m_bins.size(); }

  AXOM_HOST_DEVICE bool isValidIndex(IndexType index) const
  {
    return index >= 0 && index < getNumBins();
  }
  // getters
  AXOM_HOST_DEVICE BinStoreType& getBinContents(IndexType gridIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx];
  }

  AXOM_HOST_DEVICE T& get(IndexType gridIdx, IndexType binIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx][binIdx];
  }

  axom::ArrayView<BinStoreType> m_bins;
};

/*!
 * \brief Policy for storing uniform grid as a flat array of type T. Bins are
 *  represented as slices of this flat array.
 */
template <typename T>
struct FlatGridStorage
{
  using BinType = axom::ArrayView<T>;
  using ConstBinType = axom::ArrayView<const T>;

  using ViewType = FlatGridView<T>;
  using ConstViewType = FlatGridView<const T>;

  FlatGridStorage(int allocID)
    : m_binData(0, 0, allocID)
    , m_binOffsets(0, 0, allocID)
    , m_allocatorID(allocID)
  { }

  int getAllocatorID() const { return m_allocatorID; }

  void setNumBins(IndexType nbins) { m_binOffsets.resize(nbins); }
  IndexType getNumBins() const { return m_binOffsets.size(); }

  bool isValidIndex(IndexType index) const { return index >= 0 && index < getNumBins(); }

  template <typename ExecSpace>
  void initialize(const axom::ArrayView<const IndexType> binSizes)
  {
    axom::exclusive_scan<ExecSpace>(binSizes, m_binOffsets);
    axom::ReduceSum<ExecSpace, IndexType> total_elems(0);
    for_all<ExecSpace>(binSizes.size(), AXOM_LAMBDA(IndexType idx) { total_elems += binSizes[idx]; });
    m_binData.resize(total_elems.get());
  }

  void insert(IndexType gridIdx, T elem)
  {
    if(gridIdx + 1 < m_binOffsets.size())
    {
      IndexType flatOffset = m_binOffsets[gridIdx + 1];
      m_binData.insert(flatOffset, elem);
      // Increment offsets of following bins to account for insertion
      for(int i = gridIdx + 1; i < m_binOffsets.size(); i++)
      {
        m_binOffsets[i]++;
      }
    }
    else
    {
      // bin is at the end, just do a push_back
      m_binData.push_back(elem);
    }
  }

  void clear(IndexType gridIdx)
  {
    IndexType offset = m_binOffsets[gridIdx];
    IndexType end =
      (gridIdx + 1 < m_binOffsets.size()) ? m_binOffsets[gridIdx + 1] : m_binData.size();
    m_binData.erase(m_binData.begin() + offset, m_binData.begin() + end);
  }

  // getters
  BinType getBinContents(IndexType gridIdx)
  {
    IndexType offset = m_binOffsets[gridIdx];
    IndexType end =
      (gridIdx + 1 < m_binOffsets.size()) ? m_binOffsets[gridIdx + 1] : m_binData.size();
    return m_binData.view().subspan(offset, end - offset);
  }
  ConstBinType getBinContents(IndexType gridIdx) const
  {
    IndexType offset = m_binOffsets[gridIdx];
    IndexType end =
      (gridIdx + 1 < m_binOffsets.size()) ? m_binOffsets[gridIdx + 1] : m_binData.size();
    return m_binData.view().subspan(offset, end - offset);
  }

  T& get(IndexType gridIdx, IndexType binIdx) { return m_binData[m_binOffsets[gridIdx] + binIdx]; }

  const T& get(IndexType gridIdx, IndexType binIdx) const
  {
    return m_binData[m_binOffsets[gridIdx] + binIdx];
  }

  axom::Array<T> m_binData;
  axom::Array<IndexType> m_binOffsets;
  int m_allocatorID;
};

template <typename T>
struct FlatGridView
{
  using BaseT = typename std::remove_const<T>::type;
  using BinType = axom::ArrayView<T>;

  FlatGridView(FlatGridStorage<BaseT>& in)
    : m_binData(in.m_binData.view())
    , m_binOffsets(in.m_binOffsets.view())
  { }

  FlatGridView(const FlatGridStorage<BaseT>& in)
    : m_binData(in.m_binData.view())
    , m_binOffsets(in.m_binOffsets.view())
  { }

  AXOM_HOST_DEVICE IndexType getNumBins() const { return m_binOffsets.size(); }

  AXOM_HOST_DEVICE bool isValidIndex(IndexType index) const
  {
    return index >= 0 && index < getNumBins();
  }
  // getters
  AXOM_HOST_DEVICE BinType getBinContents(IndexType gridIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    IndexType offset = m_binOffsets[gridIdx];
    IndexType end =
      (gridIdx + 1 < m_binOffsets.size()) ? m_binOffsets[gridIdx + 1] : m_binData.size();
    return m_binData.subspan(offset, end - offset);
  }

  AXOM_HOST_DEVICE T& get(IndexType gridIdx, IndexType binIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_binData[m_binOffsets[gridIdx] + binIdx];
  }

  axom::ArrayView<T> m_binData;
  axom::ArrayView<const IndexType> m_binOffsets;
};

}  // namespace policy
}  // namespace spin
}  // namespace axom
#endif  // AXOM_SPIN_POLICY_UGRID_STORAGE_HPP
