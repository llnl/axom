// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MAKE_CLIPPER_STRATEGY_HPP
#define AXOM_MAKE_CLIPPER_STRATEGY_HPP

#include "axom/config.hpp"

// MeshClipper depends on sidre
#ifdef AXOM_ENABLE_SIDRE

#include "axom/klee/Geometry.hpp"
#include "axom/quest/MeshClipperStrategy.hpp"

namespace axom
{
namespace quest
{
namespace util
{

/*!
 * @brief Return a new MeshClipperStrategy implementation
 * to clip the specified klee Geometry.
 * @param [in] kleeGeometry
 * @param [in] name Name of strategy instance
 *
 * This method chooses the correct implementations for known
 * klee geometry formats.  It throws an error for unrecognized
 * formats.
*/
std::shared_ptr<MeshClipperStrategy> make_clipper_strategy(const axom::klee::Geometry& kleeGeometry,
                                                           const std::string& name = "");

}  // namespace util
}  // namespace quest
}  // namespace axom

#endif  // AXOM_ENABLE_SIDRE
#endif  // AXOM_MAKE_CLIPPER_STRATEGY_HPP
