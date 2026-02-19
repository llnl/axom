// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_DIMENSIONS_HPP_
#define AXOM_KLEE_DIMENSIONS_HPP_

#include "axom/fmt.hpp"

namespace axom
{
namespace klee
{
/// The dimensions that are supported for specifying operations in Klee.
enum class Dimensions : int
{
  Unspecified = 0,
  Two = 2,
  Three = 3
};

}  // namespace klee
}  // namespace axom

template <>
struct axom::fmt::formatter<axom::klee::Dimensions> : ostream_formatter
{
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  template <typename FormatContext>
  auto format(const axom::klee::Dimensions& dim, FormatContext& ctx) const
  {
    switch(dim)
    {
    case axom::klee::Dimensions::Two:
      return axom::fmt::format_to(ctx.out(), "Two");
    case axom::klee::Dimensions::Three:
      return axom::fmt::format_to(ctx.out(), "Three");
    case axom::klee::Dimensions::Unspecified:
      return axom::fmt::format_to(ctx.out(), "Unspecified");
    default:
      return axom::fmt::format_to(ctx.out(), "Unknown");
    }
  }
};

#endif  // AXOM_KLEE_DIMENSIONS_HPP_
