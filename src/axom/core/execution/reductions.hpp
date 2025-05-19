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

// NOTE: The usage style for these reductions is slightly different from how one
//       would do it using RAJA directly. This is so we can create an object
//       like this: axom::ReduceSum<axom::SEQ_EXEC, double> instead of needing to
//       do RAJA::ReduceSum<axom::execution_space<axom::SEQ_EXEC>::reduce_policy, double>.
//       The newer usage is more like how Axom does axom::for_all.

#ifdef AXOM_USE_RAJA
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
}  // namespace axom
#else
//------------------------------------------------------------------------------
namespace axom
{
namespace serial
{
// Serial reductions adapted from Ascent.
// https://github.com/Alpine-DAV/ascent/blob/develop/src/libs/ascent/runtimes/expressions/ascent_execution_policies.hpp

/*!
 * \brief A serial implementation of a ReduceSum operation.
 */
template <typename ExecSpace, typename T>
class ReduceSum
{
  static_assert<std::is_same<ExecSpace>, axom::SEQ_EXEC>::value);

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
  static_assert<std::is_same<ExecSpace>, axom::SEQ_EXEC>::value);

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
  static_assert<std::is_same<ExecSpace>, axom::SEQ_EXEC>::value);

public:
  ReduceMinLoc()
    : m_value(std::numeric_limits<T>::max())
    , m_value_ptr(&m_value)
    , m_index(-1)
    , m_index_ptr(&m_index)
  { }

  ReduceMinLoc(T v_start, index_t i_start)
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

  inline void minloc(const T v, index_t i) const
  {
    if(v < m_value_ptr[0])
    {
      m_value_ptr[0] = v;
      m_index_ptr[0] = i;
    }
  };

  inline T get() const { return m_value_ptr[0]; }

  inline index_t getLoc() const { return m_index_ptr[0]; }

private:
  T m_value;
  T *m_value_ptr;
  index_t m_index;
  index_t *m_index_ptr;
};

/*!
 * \brief A serial implementation of a ReduceMax operation.
 */
template <typename ExecSpace, typename T>
class ReduceMax
{
  static_assert<std::is_same<ExecSpace>, axom::SEQ_EXEC>::value);

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
  static_assert<std::is_same<ExecSpace>, axom::SEQ_EXEC>::value);

public:
  ReduceMaxLoc()
    : m_value(std::numeric_limits<T>::lowest())
    , m_value_ptr(&m_value)
    , m_index(-1)
    , m_index_ptr(&m_index)
  { }

  ReduceMaxLoc(T v_start, index_t i_start)
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
  inline void maxloc(const T v, index_t i) const
  {
    if(v > m_value_ptr[0])
    {
      m_value_ptr[0] = v;
      m_index_ptr[0] = i;
    }
  };

  inline T get() const { return m_value_ptr[0]; }

  inline index_t getLoc() const { return m_index_ptr[0]; }

private:
  T m_value;
  T *m_value_ptr;
  index_t m_index;
  index_t *m_index_ptr;
};

}  // namespace serial

// Use the serial implementations when we do not have RAJA.
template <typename ExecSpace, typename T>
using ReduceSum = axom::serial::ReduceSum<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceMin = axom::serial::ReduceMin<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceMinLoc = axom::serial::ReduceMinLoc<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceMax = axom::serial::ReduceMax<ExecSpace, T>;

template <typename ExecSpace, typename T>
using ReduceMaxLoc = axom::serial::ReduceMaxLoc<ExecSpace, T>;

}  // namespace axom
#endif  // AXOM_HAVE_RAJA

#endif
