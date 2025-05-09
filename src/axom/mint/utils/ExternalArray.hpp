// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_EXTERNALARRAY_HPP_
#define MINT_EXTERNALARRAY_HPP_

#include "axom/core/Array.hpp"  // to inherit
#include "axom/core/Types.hpp"

#include "axom/fmt.hpp"

#include "axom/slic/interface/slic.hpp"  // for slic logging macros

namespace axom
{
namespace mint
{
/*!
 * \class ExternalArray
 *
 * \brief Provides a generic multi-component array, constructed from external
 *  storage.
 *
 *  This ExternalArray class extends axom::Array by storing data in an
 *  externally-owned buffer. This class provides a generic multi-component
 *  array container with dynamic resizing and insertion. Each element in the
 *  array is a tuple consisting of 1 or more components, which are stored
 *  contiguously.
 *
 *  All array operations can be performed as with the base axom::Array class,
 *  with the exception of operations that require reallocation of the
 *  underlying buffer.
 *
 * \note When the Array object is deleted, it does not delete the associated
 *  data.
 *
 * \tparam T the type of the values to hold.
 */
template <typename T, int DIM = 1>
class ExternalArray : public axom::Array<T, DIM>
{
public:
  using BaseClass = axom::Array<T, DIM>;
  static_assert(DIM <= 2, "Only 1- and 2-dimensional external arrays are permitted");
  /*!
   * \brief Default constructor. Disabled.
   */
  ExternalArray() = delete;

  /*!
   * \brief Move constructor.
   * \param [in] other The array to move from
   */
  ExternalArray(ExternalArray&& other) : axom::Array<T, DIM>(std::move(other)) { }

  /// \name ExternalArray constructors
  /// @{

  /*!
   * \brief Generic constructor for an ExternalArray of arbitrary dimension
   *
   * \param [in] data the external data this ExternalArray will wrap.
   * \param [in] shape An array with the "shape" of the ExternalArray
   *
   * \post size() == num_elements
   */
  template <int UDIM = DIM, typename Enable = std::enable_if_t<UDIM != 1>>
  ExternalArray(T* data, const StackArray<IndexType, DIM>& shape, IndexType capacity)
    : axom::Array<T, DIM>()
  {
    SLIC_ASSERT(data != nullptr);

    this->m_shape = shape;
    this->updateStrides();

    SLIC_ERROR_IF(!detail::allNonNegative(shape.m_data),
                  "Dimensions passed as shape must all be non-negative.");

    this->m_num_elements = detail::packProduct(shape.m_data);
    this->m_capacity = capacity;

    if(this->m_num_elements > capacity)
    {
      SLIC_WARNING(
        fmt::format("Attempting to set number of elements greater than the available "
                    "capacity. (elements = {}, capacity = {})",
                    this->m_num_elements,
                    capacity));
      this->m_capacity = this->m_num_elements;
    }

    this->m_data = data;
  }

  /// \overload
  template <int UDIM = DIM, typename Enable = std::enable_if_t<UDIM == 1>>
  ExternalArray(T* data, IndexType size, IndexType capacity) : axom::Array<T, DIM>()
  {
    SLIC_ASSERT(data != nullptr);

    this->m_num_elements = size;
    this->m_capacity = capacity;

    if(this->m_num_elements > capacity)
    {
      SLIC_WARNING(
        fmt::format("Attempting to set number of elements greater than the available "
                    "capacity. (elements = {}, capacity = {})",
                    this->m_num_elements,
                    capacity));
      this->m_capacity = this->m_num_elements;
    }

    this->m_data = data;
  }

  /// @}

  /*!
   * Destructor.
   */
  virtual ~ExternalArray()
  {
    // The below avoids the base class constructor from erasing/destroying the
    // externally-managed buffer.
    this->updateNumElements(0);
    this->m_data = nullptr;
  }

  /*!
   * \brief Move assignment.
   * \param [in] other The ExternalArray to move from
   */
  ExternalArray& operator=(ExternalArray&& other)
  {
    axom::Array<T, DIM>::operator=(std::move(other));
    return *this;
  }

protected:
  /*!
   * \brief Set the number of elements allocated for the data array.
   *
   * \param [in] capacity the new number of elements to allocate.
   */
  virtual void setCapacity(axom::IndexType new_capacity)
  {
    if(this->m_capacity != new_capacity)
    {
      SLIC_ERROR("Cannot modify the capacity of an ExternalArray.");
    }
  }

  /*!
   * \brief Reallocates the data array when the size exceeds the capacity.
   *
   * \param [in] new_num_elements the number of elements which exceeds the current
   *  capacity.
   */
  virtual void dynamicRealloc(axom::IndexType new_num_elements)
  {
    AXOM_UNUSED_VAR(new_num_elements);
    SLIC_ERROR("Cannot increase capacity of an ExternalArray.");
  }
};

}  // namespace mint
}  // namespace axom

#endif
