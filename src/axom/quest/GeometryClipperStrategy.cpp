// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/GeometryClipperStrategy.hpp"

namespace axom
{
namespace quest
{

GeometryClipperStrategy::GeometryClipperStrategy(const klee::Geometry& kGeom)
  : m_info(kGeom.asHierarchy())
{ }

const std::string& GeometryClipperStrategy::name() const
{
  static const std::string n = "UNNAMED";
  return n;
}

}  // end namespace quest
}  // end namespace axom
