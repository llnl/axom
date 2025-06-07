// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MAKE_CLIPPER_STRATEGY_HPP
#define AXOM_MAKE_CLIPPER_STRATEGY_HPP

#include "axom/config.hpp"

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/quest/ShapeeMesh.hpp"

namespace axom
{
namespace quest
{
namespace util
{

/*!
  @brief Return a new GeometryClipperStrategy implementation
  to clip the specified klee Geometry.

  This method chooses the correct implementations for known
  klee geometry formats.  It throws an error for unrecognized
  formats.
*/
std::shared_ptr<GeometryClipperStrategy> make_clipper_strategy(
  const axom::klee::Geometry& kleeGeometry);

}  // namespace util
}  // namespace quest
}  // namespace axom

#endif  // AXOM_MAKE_CLIPPER_STRATEGY_HPP
