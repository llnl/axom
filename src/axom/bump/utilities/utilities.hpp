// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_UTILITIES_HPP_
#define AXOM_BUMP_UTILITIES_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <cstdint>
#include <string_view>

namespace axom
{
namespace bump
{
namespace utilities
{
//------------------------------------------------------------------------------
/*!
 * \brief Determines best mask type to use for an ExecSpace.
 *
 * \tparam ExecSpace The execution space.
 * \tparam DefaultMaskType The default mask type to use if we cannot use char.
 */
template <typename ExecSpace, typename DefaultMaskType = axom::IndexType>
struct mask_traits
{
  static constexpr bool is_hip()
  {
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    constexpr int BLOCK_SIZE = 256; // The size here does not matter.
    return std::string_view(axom::execution_space<ExecSpace>::name()) == std::string_view(axom::execution_space<axom::HIP_EXEC<BLOCK_SIZE>>::name());
#else
    return false;
#endif
  }

  using type = typename std::conditional<is_hip(), DefaultMaskType, char>::type;
};

//------------------------------------------------------------------------------
/*!
 * \brief This class and its specializations provide a type trait that lets us
 *        determine the type that should be used to accumulate values when we
 *        do floating point math.
 *
 * \note this belongs in algorithm utilities, maybe core.
 */
template <typename T>
struct accumulation_traits
{
  using type = float;
};

template <>
struct accumulation_traits<double>
{
  using type = double;
};

template <>
struct accumulation_traits<long>
{
  using type = double;
};

template <>
struct accumulation_traits<unsigned long>
{
  using type = double;
};

//------------------------------------------------------------------------------
/*!
 * \brief Base template for computing a shape's area or volume.
 */
template <int NDIMS>
struct ComputeShapeAmount
{ };

/*!
 * \brief 2D specialization for shapes to compute area.
 */
template <>
struct ComputeShapeAmount<2>
{
  template <typename ShapeType>
  static inline AXOM_HOST_DEVICE double execute(const ShapeType &shape)
  {
    return shape.area();
  }
};

/*!
 * \brief 3D specialization for shapes to compute volume.
 */
template <>
struct ComputeShapeAmount<3>
{
  template <typename ShapeType>
  static inline AXOM_HOST_DEVICE double execute(const ShapeType &shape)
  {
    return shape.volume();
  }
};

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom

#endif
