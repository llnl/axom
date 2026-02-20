// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MAKE_CLIPPER_STRATEGY_HPP
#define AXOM_MAKE_CLIPPER_STRATEGY_HPP

#include "axom/config.hpp"

// MeshClipper depends on sidre
#ifdef AXOM_USE_SIDRE

  #include "axom/klee/Geometry.hpp"
  #include "axom/quest/MeshClipperStrategy.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{
namespace util
{

/*!
 * @brief Return a new MeshClipperStrategy implementation
 * to clip the specified klee Geometry.
 * @param [in] kleeGeometry
 * @param [in] name Name of strategy instance.
 *   If unspecified, the selected implementation will provide a default name.
 *
 * This method chooses the correct implementations for known
 * klee geometry formats.  It issues a warning for unrecognized
 * formats.
 *
 * @return Pointer to new MeshClipperStrategy, or null if the
 * specified geometry is not an axom-provided one.
 */
std::shared_ptr<MeshClipperStrategy> make_clipper_strategy(const axom::klee::Geometry& kleeGeometry,
                                                           const std::string& name = "");

}  // namespace util
}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_USE_SIDRE
#endif  // AXOM_MAKE_CLIPPER_STRATEGY_HPP
