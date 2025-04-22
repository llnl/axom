// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IndirectionPolicies.hpp
 *
 * \brief Defines several indirection policies for slam
 *
 * Indirection policies encompass the underlying storage for indirection buffers
 * for a SLAM set, relation or map. A valid indirection policy must support the
 * following interface:
 *   * [required]
 *   * type alias IndirectionResult -- the type of the result of an indirection
 *      (const/nonconst and ref/nonref)
 *   * indirection() : IntType  -- returns the value of the element after
 *     indirection
 *   * hasIndirection(): bool -- returns whether there is an indirection
 *     buffer
 *   * isValid() : bool -- indicates whether the Indirection policy of the set
 *      is valid
 *   * [optional]
 *     * operator(): IntType -- alternate accessor for indirection
 *     * data() : ElementType* -- allows direct access to the underlying buffer
 *       (when this exists)
 *
 * \note Slam's Sets, Relations and Maps are not responsible for
 *  allocating/deallocating their own memory
 */

#ifndef SLAM_POLICIES_INDIRECTION_H_
#define SLAM_POLICIES_INDIRECTION_H_

#include "axom/core/Macros.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
namespace detail
{
/*!
 * \class IndexedIndirection
 *
 * \brief Provides a mixin class for a generic indexable indirection policy.
 *
 * \tparam BasePolicy the type of the derived indirection policy class
 *
 * \pre BasePolicy must implement the following methods:
 *   - data(): returns the underlying indirection buffer object
 *   - size(): returns the size of the backing buffer
 *   - checkIndirection(pos): checks whether the indirection and index is valid
 */
template <typename BasePolicy>
struct IndexedIndirection : public BasePolicy
{
  using PositionType = typename BasePolicy::PosType;

  using typename BasePolicy::ConstIndirectionResult;
  using typename BasePolicy::IndirectionResult;
  using ConstResultPtr = std::remove_reference_t<ConstIndirectionResult>*;
  using ResultPtr = std::remove_reference_t<IndirectionResult>*;

  using typename BasePolicy::IndirectionBufferType;
  using typename BasePolicy::IndirectionConstRefType;
  using typename BasePolicy::IndirectionPtrType;
  using typename BasePolicy::IndirectionRefType;

  using BasePolicy::BasePolicy;

  inline bool isValid(PositionType size,
                      PositionType offset,
                      PositionType stride,
                      bool verboseOutput = false) const;

  template <bool DeviceEnable = BasePolicy::DeviceAccessible>
  AXOM_HOST_DEVICE static inline std::enable_if_t<DeviceEnable, ResultPtr> getIndirection(
    IndirectionRefType buf,
    PositionType pos = 0)
  {
    return &buf[pos];
  }

  template <bool DeviceEnable = BasePolicy::DeviceAccessible>
  AXOM_HOST_DEVICE static inline std::enable_if_t<DeviceEnable, ConstResultPtr> getConstIndirection(
    IndirectionConstRefType buf,
    PositionType pos = 0)
  {
    return &buf[pos];
  }

  template <bool DeviceEnable = BasePolicy::DeviceAccessible>
  AXOM_HOST_DEVICE static inline std::enable_if_t<!DeviceEnable, ResultPtr> getIndirection(
    IndirectionRefType buf,
    PositionType pos = 0)
  {
#ifdef AXOM_DEVICE_CODE
    AXOM_UNUSED_VAR(buf);
    AXOM_UNUSED_VAR(pos);
    SLIC_ASSERT_MSG(
      false,
      BasePolicy::Name << " -- Attempting to indirect on an unsupported indirection policy.");

  // Disable no-return warnings from device code
  #if defined(__CUDA_ARCH__)
    __trap();
  #elif defined(__HIP_DEVICE_COMPILE__)
    abort();
  #endif
    return nullptr;
#else
    // Always return a value.
    return &buf[pos];
#endif
  }

  template <bool DeviceEnable = BasePolicy::DeviceAccessible>
  AXOM_HOST_DEVICE static inline std::enable_if_t<!DeviceEnable, ConstResultPtr>
  getConstIndirection(IndirectionConstRefType buf, PositionType pos = 0)
  {
#ifdef AXOM_DEVICE_CODE
    AXOM_UNUSED_VAR(buf);
    AXOM_UNUSED_VAR(pos);
    SLIC_ASSERT_MSG(
      false,
      BasePolicy::Name << " -- Attempting to indirect on an unsupported indirection policy.");

  // Disable no-return warnings from device code
  #if defined(__CUDA_ARCH__)
    __trap();
  #elif defined(__HIP_DEVICE_COMPILE__)
    abort();
  #endif
    return nullptr;
#else
    // Always return a value.
    return &buf[pos];
#endif
  }

  AXOM_HOST_DEVICE inline ConstIndirectionResult indirection(PositionType pos) const
  {
#ifndef AXOM_DEVICE_CODE
    checkIndirection(pos);
#endif
    return *IndexedIndirection::getConstIndirection(BasePolicy::data(), pos);
  }

  AXOM_HOST_DEVICE inline IndirectionResult indirection(PositionType pos)
  {
#ifndef AXOM_DEVICE_CODE
    checkIndirection(pos);
#endif
    return *IndexedIndirection::getIndirection(BasePolicy::data(), pos);
  }

  AXOM_HOST_DEVICE inline ConstIndirectionResult operator()(PositionType pos) const
  {
    return indirection(pos);
  }

  AXOM_HOST_DEVICE inline IndirectionResult operator()(PositionType pos)
  {
    return indirection(pos);
  }

private:
  void checkIndirection(PositionType pos) const
  {
    AXOM_UNUSED_VAR(pos);
    SLIC_ASSERT_MSG(this->hasIndirection(),
                    BasePolicy::Name << " -- Tried to dereference "
                                     << " a null array in an array based indirection set.");
  }
};

template <typename BasePolicy>
bool IndexedIndirection<BasePolicy>::isValid(PositionType size,
                                             PositionType offset,
                                             PositionType stride,
                                             bool verboseOutput) const
{
  AXOM_UNUSED_VAR(verboseOutput);

  // always valid if set has zero size, even if indirection buffer is null
  if(size == 0)
  {
    return true;
  }

  bool bValid = true;

  // Otherwise, check whether the set has elements, but the array ptr is null
  if(!this->hasIndirection())
  {
    SLIC_DEBUG_IF(verboseOutput,
                  "Array-based indirection set with non-zero size "
                    << "(size=" << size << ") requires a valid data buffer,"
                    << "but buffer pointer was null.");

    bValid = false;
  }
  else
  {
    // Verify underlying array has sufficient storage for all set elements
    // Note: it is valid for the data buffer to have extra space
    PositionType firstEltInd = offset;
    PositionType lastEltInd = (size - 1) * stride + offset;
    PositionType vecSize = static_cast<PositionType>(BasePolicy::size());

    bool isRangeValid =
      (0 <= firstEltInd) && (firstEltInd < vecSize) && (0 <= lastEltInd) && (lastEltInd < vecSize);

    if(!isRangeValid)
    {
      SLIC_DEBUG_IF(verboseOutput,
                    "Invalid array-based IndirectionSet -- Data buffer "
                      << "must be large enough to hold all elements of the set. "
                      << "Underlying buffer size is " << vecSize << "."
                      << " Offset of " << offset << " leads to a first index of " << firstEltInd << "."
                      << " Stride of " << stride << " and size of " << size
                      << " leads to a last index of " << lastEltInd << ".");

      bValid = false;
    }
  }

  return bValid;
}

}  // namespace detail

/**
 * \name OrderedSet_Indirection_Policies
 * \brief A few default policies for the indirection of an OrderedSet
 */

/// \{

/**
 * \brief A policy class for sets with no indirection
 */
template <typename PositionType, typename ElementType>
struct NoIndirection
{
  using IndirectionResult = ElementType;
  using ConstIndirectionResult = const ElementType;
  using IndirectionBufferType = struct
  { };
  using IndirectionPtrType = IndirectionBufferType*;

  static constexpr bool DeviceAccessible = true;

  AXOM_HOST_DEVICE NoIndirection() { }

  // This empty .ctor exists to satisfy IndirectionPolicy interface
  NoIndirection(IndirectionBufferType*) { }

  AXOM_HOST_DEVICE inline IndirectionResult indirection(PositionType pos) const
  {
    return static_cast<ElementType>(pos);
  }

  AXOM_HOST_DEVICE inline IndirectionResult operator()(PositionType pos) const
  {
    return indirection(pos);
  }

  AXOM_HOST_DEVICE IndirectionBufferType* ptr() { return nullptr; }

  bool hasIndirection() const { return false; }
  inline bool isValid(PositionType, PositionType, PositionType, bool) const { return true; }
};

template <typename PositionType, typename ElementType>
struct CArrayIndirectionBase
{
  using PosType = PositionType;
  using ElemType = ElementType;

  using IndirectionResult = ElementType&;
  using ConstIndirectionResult = const ElementType&;

  using IndirectionBufferType = ElementType*;
  using IndirectionPtrType = IndirectionBufferType;
  using IndirectionRefType = IndirectionBufferType;
  using IndirectionConstRefType = const IndirectionBufferType&;

  static constexpr bool DeviceAccessible = true;
  static constexpr bool IsMutableBuffer = false;
  static constexpr const char* Name = "SLAM::CArrayIndirection";

  AXOM_HOST_DEVICE CArrayIndirectionBase(IndirectionPtrType buf = nullptr) : m_arrBuf(buf) { }

  AXOM_HOST_DEVICE IndirectionBufferType data() const { return m_arrBuf; }
  AXOM_HOST_DEVICE IndirectionBufferType& ptr() { return m_arrBuf; }

  bool hasIndirection() const { return m_arrBuf != nullptr; }

  constexpr PositionType size() const { return axom::numeric_limits<PositionType>::max(); }

private:
  IndirectionBufferType m_arrBuf;
};

/**
 * \brief A policy class for sets with C-style array-based indirection
 */
template <typename PositionType, typename ElementType>
using CArrayIndirection =
  detail::IndexedIndirection<CArrayIndirectionBase<PositionType, ElementType>>;

template <typename PositionType, typename ElementType>
struct STLVectorIndirectionBase

{
  using PosType = PositionType;
  using ElemType = ElementType;

  using IndirectionResult = ElementType&;
  using ConstIndirectionResult = const ElementType&;

  using IndirectionBufferType = std::vector<ElementType>;
  using IndirectionPtrType = IndirectionBufferType*;
  using IndirectionRefType = IndirectionBufferType&;
  using IndirectionConstRefType = const IndirectionBufferType&;

  static constexpr bool DeviceAccessible = false;
  static constexpr bool IsMutableBuffer = true;
  static constexpr const char* Name = "SLAM::STLVectorIndirection";

  AXOM_HOST_DEVICE STLVectorIndirectionBase(IndirectionBufferType* buf = nullptr)
    : m_vecBuf(buf) { }

  AXOM_HOST_DEVICE IndirectionBufferType& data() { return *m_vecBuf; }
  AXOM_HOST_DEVICE IndirectionBufferType const& data() const { return *m_vecBuf; }

  AXOM_HOST_DEVICE IndirectionPtrType& ptr() { return m_vecBuf; }
  AXOM_HOST_DEVICE IndirectionPtrType const& ptr() const { return m_vecBuf; }

  bool hasIndirection() const { return m_vecBuf != nullptr; }

  PositionType size() const { return m_vecBuf->size(); }

  static IndirectionBufferType create(PositionType size, ConstIndirectionResult value, int allocatorId)
  {
    AXOM_UNUSED_VAR(allocatorId);
    return IndirectionBufferType(size, value);
  }

private:
  IndirectionBufferType* m_vecBuf;
};

/**
 * \brief A policy class for sets with stl vector-based indirection
 */
template <typename PositionType, typename ElementType>
using STLVectorIndirection =
  detail::IndexedIndirection<STLVectorIndirectionBase<PositionType, ElementType>>;

template <typename PositionType, typename ElementType>
struct ArrayIndirectionBase
{
  using PosType = PositionType;
  using ElemType = ElementType;

  using IndirectionResult = ElementType&;
  using ConstIndirectionResult = const ElementType&;

  using IndirectionBufferType = axom::Array<ElementType>;
  using IndirectionPtrType = IndirectionBufferType*;
  using IndirectionRefType = IndirectionBufferType&;
  using IndirectionConstRefType = const IndirectionBufferType&;

  static constexpr bool DeviceAccessible = true;
  static constexpr bool IsMutableBuffer = true;
  static constexpr const char* Name = "SLAM::ArrayIndirection";

  AXOM_HOST_DEVICE ArrayIndirectionBase(IndirectionBufferType* buf = nullptr) : m_vecBuf(buf) { }

  AXOM_HOST_DEVICE IndirectionBufferType& data() { return *m_vecBuf; }
  AXOM_HOST_DEVICE IndirectionBufferType const& data() const { return *m_vecBuf; }

  AXOM_HOST_DEVICE IndirectionPtrType& ptr() { return m_vecBuf; }
  AXOM_HOST_DEVICE IndirectionPtrType const& ptr() const { return m_vecBuf; }

  bool hasIndirection() const { return m_vecBuf != nullptr; }

  PositionType size() const { return m_vecBuf->size(); }

  static IndirectionBufferType create(PositionType size, ConstIndirectionResult value, int allocatorID)
  {
    IndirectionBufferType buf(size, size, allocatorID);
    buf.fill(value);
    return buf;
  }

private:
  IndirectionBufferType* m_vecBuf;
};

/**
 * \brief A policy class for sets with axom::Array-based indirection
 */
template <typename PositionType, typename ElementType>
using ArrayIndirection = detail::IndexedIndirection<ArrayIndirectionBase<PositionType, ElementType>>;

template <typename PositionType, typename ElementType>
struct ArrayViewIndirectionBase
{
  using PosType = PositionType;
  using ElemType = ElementType;

  using IndirectionResult = ElementType&;
  using ConstIndirectionResult = ElementType&;

  using IndirectionBufferType = axom::ArrayView<ElementType>;
  using IndirectionPtrType = IndirectionBufferType;
  using IndirectionRefType = IndirectionBufferType;
  using IndirectionConstRefType = const IndirectionBufferType;

  static constexpr bool DeviceAccessible = true;
  static constexpr bool IsMutableBuffer = false;
  static constexpr const char* Name = "SLAM::ArrayViewIndirection";

  AXOM_HOST_DEVICE ArrayViewIndirectionBase(IndirectionBufferType buf = {}) : m_vecBuf(buf) { }

  ArrayViewIndirectionBase(ArrayIndirection<PositionType, ElementType> ind)
    : m_vecBuf(ind.data().data(), ind.data().size())
  { }

  ArrayViewIndirectionBase(STLVectorIndirection<PositionType, ElementType> ind)
    : m_vecBuf(ind.data().data(), ind.data().size())
  { }

  AXOM_HOST_DEVICE IndirectionBufferType data() { return m_vecBuf; }
  AXOM_HOST_DEVICE const IndirectionBufferType data() const { return m_vecBuf; }

  AXOM_HOST_DEVICE IndirectionBufferType& ptr() { return m_vecBuf; }
  AXOM_HOST_DEVICE const IndirectionBufferType& ptr() const { return m_vecBuf; }

  bool hasIndirection() const { return m_vecBuf.data() != nullptr; }

  PositionType size() const { return m_vecBuf.size(); }

private:
  IndirectionBufferType m_vecBuf;
};

/**
 * \brief A policy class for sets with axom::ArrayView-based indirection
 */
template <typename PositionType, typename ElementType>
using ArrayViewIndirection =
  detail::IndexedIndirection<ArrayViewIndirectionBase<PositionType, ElementType>>;

/// \}

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICIES_INDIRECTION_H_
