// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_UTILITIES_HPP_
#define AXOM_MIR_DISPATCH_UTILITIES_HPP_

namespace axom
{
namespace mir
{
namespace views
{
#if __cplusplus >= 201703L
// C++17 and later.
template <typename... Dimensions>
constexpr int encode_dimensions(Dimensions... dims)
{
  return (... | dims);
}
#else
template <typename T>
constexpr int encode_dimensions_impl(T arg)
{
  return arg;
}

template <typename T, typename... Dimensions>
constexpr int encode_dimensions_impl(T arg, Dimensions... dims)
{
  return (arg | encode_dimensions_impl(dims...));
}

template <typename... Dimensions>
constexpr int encode_dimensions(Dimensions... dims)
{
  return encode_dimensions_impl(dims...);
}
#endif

template <typename... Dimensions>
constexpr int select_dimensions(Dimensions... dims)
{
  return encode_dimensions((1 << dims)...);
}

constexpr bool dimension_selected(int encoded_dims, int dim) { return encoded_dims & (1 << dim); }

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
