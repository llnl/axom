// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/util/make_clipper_strategy.hpp"
#include "axom/quest/Plane3DClipper.hpp"
#include "axom/quest/TetClipper.hpp"
#include "axom/quest/TetMeshClipper.hpp"
#include "axom/quest/HexClipper.hpp"
#include "axom/quest/SphereClipper.hpp"
#include "axom/quest/SorClipper.hpp"
#include "axom/slic/interface/slic_macros.hpp"

namespace axom
{
namespace quest
{
namespace util
{

std::shared_ptr<GeometryClipperStrategy> make_clipper_strategy(
  const axom::klee::Geometry& kleeGeometry)
{
  std::shared_ptr<GeometryClipperStrategy> strategy;

  const std::string& format = kleeGeometry.getFormat();

  if(format == "plane3D")
  {
    strategy.reset(new Plane3DClipper(kleeGeometry, format));
  }
  else if(format == "tet3D")
  {
    strategy.reset(new TetClipper(kleeGeometry, format));
  }
  else if(format == "blueprint-tets")
  {
    strategy.reset(new TetMeshClipper(kleeGeometry, format));
  }
  else if(format == "hex3D")
  {
    strategy.reset(new HexClipper(kleeGeometry, format));
  }
  else if(format == "sphere3D")
  {
    strategy.reset(new SphereClipper(kleeGeometry, format));
  }
  else if(format == "sor3D")
  {
    strategy.reset(new SorClipper(kleeGeometry, format));
  }
  else
  {
    SLIC_ERROR(axom::fmt::format("Unrecognized Klee Geometry format {}.", format));
  }

  return strategy;
}

}  // namespace util
}  // namespace quest
}  // namespace axom
