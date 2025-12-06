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
  /// \brief Default-constructs an empty ArrayWrapper.
  ArrayWrapper() = default;

  /// \brief Assigns the underlying variant with a specific array.
  template <typename ArrayT>
  ArrayWrapper& operator=(ArrayT&& arr)
  {
    m_array = std::forward<ArrayT>(arr);
    return *this;
  }

  /*!
   * \brief Pushes a value to the back of the wrapped array.
   *
   * \param [in] value The value to append.
   *
   * \see axom::Array::push_back
   */
  void push_back(const T& value)
  {
    std::visit([&](auto& arr) { arr.push_back(value); }, m_array);
  }

  /*!
   * \brief Inserts elements into the wrapped array at a position.
   *
   * \param [in] pos    The position at which to begin insertion.
   * \param [in] n      The number of elements to insert.
   * \param [in] values Pointer to the elements to insert; must reference at
   *                    least \p n values.
   *
   * \see axom::Array::insert
   */
  void insert(IndexType pos, IndexType n, const T* values)
  {
    std::visit([&](auto& arr) { arr.insert(pos, n, values); }, m_array);
  }

  /*!
   * \brief Inserts elements from a view into the wrapped array.
   *
   * Inserts the contents of \p span at the given index of the wrapped array.
   *
   * \param [in] index The position at which to begin insertion.
   * \param [in] span  View containing the elements to insert.
   *
   * \see axom::Array::insert
   */
  void insert(IndexType index, ArrayView<const T, DIM> span)
  {
    std::visit([&](auto& arr) { arr.insert(index, span); }, m_array);
  }

  /*!
   * \brief Appends elements from a view to the wrapped array.
   *
   * Elements in \p span are appended at the end of the wrapped array.
   *
   * \param [in] span View containing the elements to append.
   *
   * \see axom::Array::append
   */
  void append(ArrayView<const T, DIM> span)
  {
    std::visit([&](auto& arr) { arr.append(span); }, m_array);
  }

  /*!
   * \brief Sets a range of elements in the wrapped array.
   *
   * Copies \p n elements from \p elements into the wrapped array starting
   * at position \p pos.
   *
   * \param [in] elements Pointer to the source elements.
   * \param [in] n        Number of elements to copy.
   * \param [in] pos      Starting position in the wrapped array.
   *
   * \see axom::Array::set
   */
  void set(const T* elements, IndexType n, IndexType pos)
  {
    std::visit([&](auto& arr) { arr.set(elements, n, pos); }, m_array);
  }

  /*!
   * \brief Resizes the wrapped array.
   *
   * Forwards the given size arguments to the underlying array's \c resize
   * method.
   *
   * \param [in] args New extents for each dimension.
   *
   * \see axom::Array::resize
   */
  template <typename... Args>
  void resize(Args... args)
  {
    std::visit([&](auto& arr) { arr.resize(args...); }, m_array);
  }

  /*!
   * \brief Reserves capacity in the wrapped array.
   *
   * Ensures that the wrapped array can hold at least \p newCapacity
   * elements without reallocation.
   *
   * \param [in] newCapacity The new capacity in number of elements.
   *
   * \see axom::Array::reserve
   */
  void reserve(IndexType newCapacity)
  {
    std::visit([&](auto& arr) { arr.reserve(newCapacity); }, m_array);
  }

  /*!
   * \brief Shrinks the wrapped array's capacity to its size.
   *
   * \see axom::Array::shrink
   */
  void shrink()
  {
    std::visit([&](auto& arr) { arr.shrink(); }, m_array);
  }

  /*!
   * \brief Sets the resize ratio of the wrapped array.
   *
   * Controls the factor by which capacity grows when the wrapped array
   * is reallocated.
   *
   * \param [in] ratio New resize ratio.
   *
   * \see axom::Array::setResizeRatio
   */
  void setResizeRatio(double ratio)
  {
    std::visit([&](auto& arr) { arr.setResizeRatio(ratio); }, m_array);
  }

  /*!
   * \brief Returns the resize ratio used by the wrapped array.
   *
   * \see axom::Array::getResizeRatio
   */
  double getResizeRatio() const
  {
    return std::visit([](const auto& arr) -> double { return arr.getResizeRatio(); }, m_array);
  }

  /*!
   * \brief Returns the number of elements in the wrapped array.
   *
   * \see axom::Array::size
   */
  IndexType size() const
  {
    return std::visit([](const auto& arr) -> IndexType { return arr.size(); }, m_array);
  }

  /*!
   * \brief Returns the capacity of the wrapped array.
   *
   * \see axom::Array::capacity
   */
  IndexType capacity() const
  {
    return std::visit([](const auto& arr) -> IndexType { return arr.capacity(); }, m_array);
  }

  /*!
   * \brief Returns true if the wrapped array contains no elements.
   *
   * \see axom::Array::empty
   */
  bool empty() const
  {
    return std::visit([](const auto& arr) -> bool { return arr.empty(); }, m_array);
  }

  /*!
   * \brief Returns a pointer to the wrapped array's data buffer.
   *
   * \see axom::Array::data
   */
  T* data()
  {
    return std::visit([](auto& arr) -> T* { return arr.data(); }, m_array);
  }

  /*!
   * \brief Returns a const pointer to the wrapped array's data buffer.
   *
   * \see axom::Array::data
   */
  const T* data() const
  {
    return std::visit([](const auto& arr) -> const T* { return arr.data(); }, m_array);
  }

  /*!
   * \brief Returns a pointer to the wrapped data buffer.
   *
   * Equivalent to calling \c data().
   *
   * \see axom::Array::data
   */
  T* operator*() { return data(); }

  /// \overload
  const T* operator*() const { return data(); }

  /*!
   * \brief Returns the logical shape of the wrapped array.
   *
   * \see axom::Array::shape
   */
  axom::StackArray<IndexType, DIM> shape() const
  {
    return std::visit([](const auto& arr) -> axom::StackArray<IndexType, DIM> { return arr.shape(); },
                      m_array);
  }

  /*!
   * \brief Returns the underlying variant object storing the wrapped arrays.
   */
  const ArrayVariant& getVariant() const { return m_array; }

private:
  ArrayVariant m_array;
};

}  // namespace detail
}  // namespace mint
}  // namespace axom

#endif
