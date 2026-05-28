// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"

#include "axom/quest/MeshClipperStrategy.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{
namespace detail
{

/*!
 * \brief Compute the positive-side volume of an affine scalar field over a tetrahedron.
 *
 * The field is specified by its values at the tetrahedron vertices. Since the
 * zero isosurface is planar, its clipped-volume fraction can be computed
 * directly from the four values without constructing intersection geometry.
 */
AXOM_HOST_DEVICE inline double clipTetByVertexValues(const double values[4], double tetVolume)
{
  int inside[4];
  int outside[4];
  int insideCount = 0;
  int outsideCount = 0;
  for(int i = 0; i < 4; ++i)
  {
    if(values[i] >= 0.0)
    {
      inside[insideCount++] = i;
    }
    else
    {
      outside[outsideCount++] = i;
    }
  }

  if(insideCount == 0)
  {
    return 0.0;
  }
  if(insideCount == 4)
  {
    return tetVolume;
  }

  if(insideCount == 1 || insideCount == 3)
  {
    const bool invert = insideCount == 3;
    const int isolated = invert ? outside[0] : inside[0];
    const int* others = invert ? inside : outside;
    double fraction = 1.0;
    for(int i = 0; i < 3; ++i)
    {
      fraction *= values[isolated] / (values[isolated] - values[others[i]]);
    }
    return invert ? tetVolume * (1.0 - fraction) : tetVolume * fraction;
  }

  const int i0 = inside[0];
  const int i1 = inside[1];
  const int o0 = outside[0];
  const int o1 = outside[1];
  const double w00 = values[i0] / (values[i0] - values[o0]);
  const double w01 = values[i0] / (values[i0] - values[o1]);
  const double w10 = values[i1] / (values[i1] - values[o0]);
  const double w11 = values[i1] / (values[i1] - values[o1]);
  const double fraction = w00 * w01 + w00 * w11 * (1.0 - w01) + w10 * w11 * (1.0 - w00);
  return tetVolume * fraction;
}

}  // namespace detail
}  // namespace experimental
}  // namespace quest
}  // namespace axom
