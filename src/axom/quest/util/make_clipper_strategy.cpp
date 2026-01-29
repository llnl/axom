// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifdef AXOM_USE_SIDRE
  #include "axom/quest/util/make_clipper_strategy.hpp"
  #include "axom/quest/detail/clipping/Plane3DClipper.hpp"
  #include "axom/quest/detail/clipping/TetClipper.hpp"
  #include "axom/quest/detail/clipping/TetMeshClipper.hpp"
  #include "axom/quest/detail/clipping/HexClipper.hpp"
  #include "axom/quest/detail/clipping/SphereClipper.hpp"
  #include "axom/quest/detail/clipping/MonotonicZSORClipper.hpp"
  #include "axom/quest/detail/clipping/SORClipper.hpp"
  #include "axom/slic/interface/slic_macros.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{
namespace util
{

std::shared_ptr<MeshClipperStrategy> make_clipper_strategy(const axom::klee::Geometry& kleeGeometry,
                                                           const std::string& name)
{
  std::shared_ptr<MeshClipperStrategy> strategy;

  const std::string& format = kleeGeometry.getFormat();
  const std::string& instanceName = !name.empty() ? name : kleeGeometry.getFormat();

  if(format == "plane3D")
  {
    strategy.reset(new Plane3DClipper(kleeGeometry, instanceName));
  }
  else if(format == "tet3D")
  {
    strategy.reset(new TetClipper(kleeGeometry, instanceName));
  }
  else if(format == "blueprint-tets")
  {
    strategy.reset(new TetMeshClipper(kleeGeometry, instanceName));
  }
  else if(format == "hex3D")
  {
    strategy.reset(new HexClipper(kleeGeometry, instanceName));
  }
  else if(format == "sphere3D")
  {
    strategy.reset(new SphereClipper(kleeGeometry, instanceName));
  }
  else if(format == "sor3D")
  {
    strategy.reset(new SORClipper(kleeGeometry, instanceName));
  }
  else if(format == "cyl3D")
  {
    strategy.reset(new MonotonicZSORClipper(kleeGeometry, instanceName));
  }
  else if(format == "cone3D")
  {
    strategy.reset(new MonotonicZSORClipper(kleeGeometry, instanceName));
  }
  else
  {
    SLIC_WARNING(axom::fmt::format("Unrecognized Klee Geometry format {}.", format));
  }

  return strategy;
}

}  // namespace util
}  // namespace experimental
}  // namespace quest
}  // namespace axom
#endif  // AXOM_USE_SIDRE
