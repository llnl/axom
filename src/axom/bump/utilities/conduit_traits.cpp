// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/bump/utilities/conduit_traits.hpp"

namespace axom
{
namespace bump
{
namespace utilities
{

// Static data. These originally appeared as constexpr members in the header file
// but there were linker errors despite constexpr.
const char *cpp2conduit<conduit::int8>::name = "int8";
const char *cpp2conduit<conduit::int16>::name = "int16";
const char *cpp2conduit<conduit::int32>::name = "int32";
const char *cpp2conduit<conduit::int64>::name = "int64";
const char *cpp2conduit<conduit::uint8>::name = "uint8";
const char *cpp2conduit<conduit::uint16>::name = "uint16";
const char *cpp2conduit<conduit::uint32>::name = "uint32";
const char *cpp2conduit<conduit::uint64>::name = "uint64";
const char *cpp2conduit<conduit::float32>::name = "float32";
const char *cpp2conduit<conduit::float64>::name = "float64";

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom
