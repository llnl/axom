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
namespace detail
{
/*!
 * \class ExternalStoragePolicy
 *
 * \brief Externally-managed storage policy.
 *  Reallocation is not supported.
 */
template <typename T>
struct ExternalStoragePolicy
{
  /*!
   * \brief Callback to report changes in shape/size of valid data in Array.
   *
   * \param [in] shape the current dimensions of the array
   * \param [in] size the current number of elements stored in the array
   */
  template <int Dims>
  void onShapeUpdate(StackArray<IndexType, Dims> AXOM_UNUSED_PARAM(shape))
  { }

  /*!
   * \brief Reallocates a buffer. No-op for ExternalArray.
   */
  template <typename Func>
  T* reallocate(T* AXOM_UNUSED_PARAM(old_data),
                int old_capacity,
                int AXOM_UNUSED_PARAM(allocator_id),
                int new_capacity,
                Func&& AXOM_UNUSED_PARAM(nontrivial_move))
  {
    if(old_capacity != new_capacity)
    {
      SLIC_ERROR("Cannot increase capacity of an ExternalArray.");
    }
    return nullptr;
  }

  /*!
   * \brief Frees a buffer. No-op for ExternalArray.
   */
  void deallocate(void* AXOM_UNUSED_PARAM(data)) { }
};

}  // namespace detail

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
class ExternalArray
  : public axom::Array<T, DIM, MemorySpace::Dynamic, detail::ExternalStoragePolicy<T>>
{
public:
  using BaseClass = axom::Array<T, DIM, MemorySpace::Dynamic, detail::ExternalStoragePolicy<T>>;
  static_assert(DIM <= 2, "Only 1- and 2-dimensional external arrays are permitted");
  /*!
   * \brief Default constructor. Disabled.
   */
  ExternalArray() = delete;

  /*!
   * \brief Move constructor.
   * \param [in] other The array to move from
   */
  ExternalArray(ExternalArray&& other) = default;

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
  ExternalArray(T* data, const StackArray<IndexType, DIM>& shape, IndexType capacity) : BaseClass()
  {
    SLIC_ASSERT(data != nullptr);

    this->m_shape = shape;
    this->updateStrides();

    SLIC_ERROR_IF(!axom::detail::allNonNegative(shape.m_data),
                  "Dimensions passed as shape must all be non-negative.");

    this->m_num_elements = axom::detail::packProduct(shape.m_data);
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
  ExternalArray(T* data, IndexType size, IndexType capacity) : BaseClass()
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
  virtual ~ExternalArray() = default;

  /*!
   * \brief Move assignment.
   * \param [in] other The ExternalArray to move from
   */
  ExternalArray& operator=(ExternalArray&& other) = default;
};

}  // namespace mint
}  // namespace axom

#endif
