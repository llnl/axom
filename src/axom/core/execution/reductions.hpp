// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_REDUCTIONS_HPP_
#define AXOM_CORE_EXECUTION_REDUCTIONS_HPP_

#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

// NOTE: Reduction operations used in Axom have been wrapped here. Additional
//       RAJA reductions may need to be wrapped over time.

#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"

//------------------------------------------------------------------------------
namespace axom
{
// Axom includes RAJA so use RAJA reductions.
template <typename ExecSpace, typename T>
using ReduceSum = RAJA::ReduceSum<typename axom::execution_space<ExecSpace>::reduce_policy, T>;

template <typename ExecSpace, typename T>
using ReduceMin = RAJA::ReduceMin<typename axom::execution_space<ExecSpace>::reduce_policy, T>;

template <typename ExecSpace, typename T>
using ReduceMinLoc = RAJA::ReduceMinLoc<typename axom::execution_space<ExecSpace>::reduce_policy, T>;

template <typename ExecSpace, typename T>
using ReduceMax = RAJA::ReduceMax<typename axom::execution_space<ExecSpace>::reduce_policy, T>;

template <typename ExecSpace, typename T>
using ReduceMaxLoc = RAJA::ReduceMaxLoc<typename axom::execution_space<ExecSpace>::reduce_policy, T>;

template <typename ExecSpace, typename T>
using ReduceBitAnd = RAJA::ReduceBitAnd<typename axom::execution_space<ExecSpace>::reduce_policy, T>;

template <typename ExecSpace, typename T>
using ReduceBitOr = RAJA::ReduceBitOr<typename axom::execution_space<ExecSpace>::reduce_policy, T>;

}  // namespace axom
#else
//------------------------------------------------------------------------------
namespace axom
{
namespace serial
{
namespace reductions
{
// Serial reductions adapted from Ascent.
// https://github.com/Alpine-DAV/ascent/blob/develop/src/libs/ascent/runtimes/expressions/ascent_execution_policies.hpp

/*!
 * \brief A serial implementation of a ReduceSum operation.
 */
template <typename ExecSpace, typename T>
class ReduceSum
{
  static constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

public:
  ReduceSum() : m_value(0), m_value_ptr(&m_value) { }

  ReduceSum(T v_start) : m_value(v_start), m_value_ptr(&m_value) { }

  ReduceSum(const ReduceSum &v)
    : m_value(v.m_value)
    ,                           // will be unused in copies
    m_value_ptr(v.m_value_ptr)  // this is where the magic happens
  { }

  void operator+=(const T value) const { m_value_ptr[0] += value; }

  void sum(const T value) const { m_value_ptr[0] += value; }

  T get() const { return m_value; }

private:
  T m_value;
  T *m_value_ptr;
};

/*!
 * \brief A serial implementation of a ReduceMin operation.
 */
template <typename ExecSpace, typename T>
class ReduceMin
{
  static constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

public:
  ReduceMin() : m_value(std::numeric_limits<T>::max()), m_value_ptr(&m_value) { }

  ReduceMin(T v_start) : m_value(v_start), m_value_ptr(&m_value) { }

  ReduceMin(const ReduceMin &v)
    : m_value(v.m_value)
    ,                           // will be unused in copies
    m_value_ptr(v.m_value_ptr)  // this is where the magic happens
  { }

  void min(const T value) const
  {
    if(value < m_value_ptr[0])
    {
      m_value_ptr[0] = value;
    }
  }

  T get() const { return m_value_ptr[0]; }

private:
  T m_value;
  T *m_value_ptr;
};

/*!
 * \brief A serial implementation of a ReduceMinLoc operation.
 */
template <typename ExecSpace, typename T>
class ReduceMinLoc
{
  static constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

public:
  ReduceMinLoc()
    : m_value(std::numeric_limits<T>::max())
    , m_value_ptr(&m_value)
    , m_index(-1)
    , m_index_ptr(&m_index)
  { }

  ReduceMinLoc(T v_start, axom::IndexType i_start)
    : m_value(v_start)
    , m_value_ptr(&m_value)
    , m_index(i_start)
    , m_index_ptr(&m_index)
  { }

  ReduceMinLoc(const ReduceMinLoc &v)
    : m_value(v.m_value)
    ,  // will be unused in copies
    m_value_ptr(v.m_value_ptr)
    ,  // this is where the magic happens
    m_index(v.m_index)
    ,                           // will be unused in copies
    m_index_ptr(v.m_index_ptr)  // this is where the magic happens
  { }

  inline void minloc(const T v, axom::IndexType i) const
  {
    if(v < m_value_ptr[0])
    {
      m_value_ptr[0] = v;
      m_index_ptr[0] = i;
    }
  };

  inline T get() const { return m_value_ptr[0]; }

  inline axom::IndexType getLoc() const { return m_index_ptr[0]; }

private:
  T m_value;
  T *m_value_ptr;
  axom::IndexType m_index;
  axom::IndexType *m_index_ptr;
};

/*!
 * \brief A serial implementation of a ReduceMax operation.
 */
template <typename ExecSpace, typename T>
class ReduceMax
{
  static constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

public:
  ReduceMax() : m_value(std::numeric_limits<T>::lowest()), m_value_ptr(&m_value) { }

  ReduceMax(T v_start) : m_value(v_start), m_value_ptr(&m_value) { }

  ReduceMax(const ReduceMax &v)
    : m_value(v.m_value)
    ,                           // will be unused in copies
    m_value_ptr(v.m_value_ptr)  // this is where the magic happens
  { }

  // The const crimes we commit here are in the name of [=] capture
  void max(const T value) const
  {
    if(value > m_value_ptr[0])
    {
      m_value_ptr[0] = value;
    }
  }

  T get() const { return m_value_ptr[0]; }

private:
  T m_value;
  T *m_value_ptr;
};

/*!
 * \brief A serial implementation of a ReduceMaxLoc operation.
 */
template <typename ExecSpace, typename T>
class ReduceMaxLoc
{
  static constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

public:
  ReduceMaxLoc()
    : m_value(std::numeric_limits<T>::lowest())
    , m_value_ptr(&m_value)
    , m_index(-1)
    , m_index_ptr(&m_index)
  { }

  ReduceMaxLoc(T v_start, axom::IndexType i_start)
    : m_value(v_start)
    , m_value_ptr(&m_value)
    , m_index(i_start)
    , m_index_ptr(&m_index)
  { }

  ReduceMaxLoc(const ReduceMaxLoc &v)
    : m_value(v.m_value)
    ,  // will be unused in copies
    m_value_ptr(v.m_value_ptr)
    ,  // this is where the magic happens
    m_index(v.m_index)
    ,                           // will be unused in copies
    m_index_ptr(v.m_index_ptr)  // this is where the magic happens
  { }

  // the const crimes we commit here are in the name of [=] capture
  inline void maxloc(const T v, axom::IndexType i) const
  {
    if(v > m_value_ptr[0])
    {
      m_value_ptr[0] = v;
      m_index_ptr[0] = i;
    }
  };

  inline T get() const { return m_value_ptr[0]; }

  inline axom::IndexType getLoc() const { return m_index_ptr[0]; }

private:
  T m_value;
  T *m_value_ptr;
  axom::IndexType m_index;
  axom::IndexType *m_index_ptr;
};

/*!
 * \brief A serial implementation of a ReduceBitAnd operation.
 */
template <typename ExecSpace, typename T>
class ReduceBitAnd
{
  static constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

public:
  ReduceBitAnd() : m_value(0), m_value_ptr(&m_value) { }

  ReduceBitAnd(T v_start) : m_value(v_start), m_value_ptr(&m_value) { }

  ReduceBitAnd(const ReduceBitAnd &v)
    : m_value(v.m_value)
    ,                           // will be unused in copies
    m_value_ptr(v.m_value_ptr)  // this is where the magic happens
  { }

  void operator&=(const T value) const { m_value_ptr[0] &= value; }

  T get() const { return m_value; }

private:
  T m_value;
  T *m_value_ptr;
};

/*!
 * \brief A serial implementation of a ReduceBitOr operation.
 */
template <typename ExecSpace, typename T>
class ReduceBitOr
{
  static constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

public:
  ReduceBitOr() : m_value(0), m_value_ptr(&m_value) { }

  ReduceBitOr(T v_start) : m_value(v_start), m_value_ptr(&m_value) { }

  ReduceBitOr(const ReduceBitOr &v)
    : m_value(v.m_value)
    ,                           // will be unused in copies
    m_value_ptr(v.m_value_ptr)  // this is where the magic happens
  { }

  void operator|=(const T value) const { m_value_ptr[0] |= value; }

  T get() const { return m_value; }

private:
  T m_value;
  T *m_value_ptr;
};

}  // namespace reductions
}  // namespace serial

// Use the serial implementations when we do not have RAJA.
template <typename ExecSpace, typename T>
using ReduceSum = axom::serial::reductions::ReduceSum<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceMin = axom::serial::reductions::ReduceMin<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceMinLoc = axom::serial::reductions::ReduceMinLoc<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceMax = axom::serial::reductions::ReduceMax<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceMaxLoc = axom::serial::reductions::ReduceMaxLoc<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceBitAnd = axom::serial::reductions::ReduceBitAnd<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceBitOr = axom::serial::reductions::ReduceBitOr<ExecSpace, T>;

}  // namespace axom
#endif  // AXOM_HAVE_RAJA

#endif
