// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Mint_ArrayWrapper_HPP
#define Axom_Mint_ArrayWrapper_HPP

#include "axom/core/Array.hpp"  // to inherit
#include "axom/core/Types.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
  #include "axom/mint/deprecated/SidreMCArray.hpp"
#endif
#include "axom/mint/utils/ExternalArray.hpp"

#include <variant>

namespace axom
{
namespace mint
{
namespace detail
{
/*!
 * \class ArrayWrapper
 *
 * \brief Wrapper around common array operations which holds:
 *   - axom::Array if we own the memory
 *   - sidre::Array if the memory is stored in Sidre
 *   - mint::utilities::ExternalArray if the memory is externally-owned
 */
template <typename T, int DIM = 1>
class ArrayWrapper
{
private:
  using ArrayVariant = std::variant<axom::Array<T, DIM>,
#ifdef AXOM_MINT_USE_SIDRE
                                    axom::sidre::Array<T, DIM>,
#endif
                                    axom::mint::ExternalArray<T, DIM>>;

public:
  ArrayWrapper() = default;

  template <typename ArrayT>
  ArrayWrapper& operator=(ArrayT&& arr)
  {
    m_array = std::forward<ArrayT>(arr);
    return *this;
  }

  void insert(IndexType index, ArrayView<const T, DIM> span)
  {
    std::visit([&](auto& arr) { arr.insert(index, span); }, m_array);
  }
  void append(ArrayView<const T, DIM> span)
  {
    std::visit([&](auto& arr) { arr.append(span); }, m_array);
  }
  void set(const T* elements, IndexType n, IndexType pos)
  {
    std::visit([&](auto& arr) { arr.set(elements, n, pos); }, m_array);
  }

  template <typename... Args>
  void resize(Args... args)
  {
    std::visit([&](auto& arr) { arr.resize(args...); }, m_array);
  }

  void reserve(IndexType newCapacity)
  {
    std::visit([&](auto& arr) { arr.reserve(newCapacity); }, m_array);
  }

  void shrink()
  {
    std::visit([&](auto& arr) { arr.shrink(); }, m_array);
  }

  void setResizeRatio(double ratio)
  {
    std::visit([&](auto& arr) { arr.setResizeRatio(ratio); }, m_array);
  }

  double getResizeRatio() const
  {
    return std::visit([](const auto& arr) -> double { return arr.getResizeRatio(); }, m_array);
  }

  IndexType size() const
  {
    return std::visit([](const auto& arr) -> IndexType { return arr.size(); }, m_array);
  }

  IndexType capacity() const
  {
    return std::visit([](const auto& arr) -> IndexType { return arr.capacity(); }, m_array);
  }

  bool empty() const
  {
    return std::visit([](const auto& arr) -> bool { return arr.empty(); }, m_array);
  }

  T* data()
  {
    return std::visit([](auto& arr) -> T* { return arr.data(); }, m_array);
  }

  const T* data() const
  {
    return std::visit([](const auto& arr) -> const T* { return arr.data(); }, m_array);
  }

  axom::StackArray<IndexType, DIM> shape() const
  {
    return std::visit([](const auto& arr) -> axom::StackArray<IndexType, DIM> { return arr.shape(); },
                      m_array);
  }

  const ArrayVariant& getVariant() const { return m_array; }

private:
  ArrayVariant m_array;
};

}  // namespace detail
}  // namespace mint
}  // namespace axom

#endif
