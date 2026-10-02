// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file quest_marching_cubes_testing_helpers.hpp
 *
 * @brief Shared setup for the quest::MarchingCubes test suites.
 *
 * The suites share analytic fields, Blueprint memory copies, and structured mesh builders.
 * Each suite defines its own checks.
 */

#pragma once

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"

#include "axom/bump/utilities/conduit_memory.hpp"

#include <conduit/conduit.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace axom::quest::testing::marching_cubes
{

using RuntimePolicy = axom::runtime_policy::Policy;

//! @brief Return each runtime policy enabled in this build.
inline std::vector<RuntimePolicy> enabledPolicies()
{
  std::vector<RuntimePolicy> policies {RuntimePolicy::seq};
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  policies.push_back(RuntimePolicy::omp);
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  policies.push_back(RuntimePolicy::cuda);
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  policies.push_back(RuntimePolicy::hip);
#endif
  return policies;
}

//! @brief Return the printable name of a runtime policy.
inline std::string policyName(RuntimePolicy policy)
{
  return axom::runtime_policy::policyToName(policy);
}

//---------------------------------------------------------------------------
// Memory copies
//---------------------------------------------------------------------------

inline int hostAllocatorID() { return axom::execution_space<axom::SEQ_EXEC>::allocatorID(); }

/*!
 * @brief Copy a Blueprint tree into memory accessible to a runtime policy.
 *
 * Tests build meshes on the host, then copy them before calling MarchingCubes.
 */
inline void copyBlueprintToPolicy(conduit::Node& dst,
                                  const conduit::Node& src,
                                  RuntimePolicy policy,
                                  int allocatorID)
{
  namespace bputils = axom::bump::utilities;

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  if(policy == RuntimePolicy::cuda)
  {
    bputils::copy<axom::CUDA_EXEC<256>>(dst, src, allocatorID);
    return;
  }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  if(policy == RuntimePolicy::hip)
  {
    bputils::copy<axom::HIP_EXEC<256>>(dst, src, allocatorID);
    return;
  }
#endif

  AXOM_UNUSED_VAR(policy);
  AXOM_UNUSED_VAR(allocatorID);
  dst.set(src);
}

//! @brief Copy a Blueprint tree back to host memory for inspection.
inline void copyBlueprintToHost(conduit::Node& dst, const conduit::Node& src)
{
  axom::bump::utilities::copy<axom::SEQ_EXEC>(dst, src, hostAllocatorID());
}

//! @brief Shift the explicit x coordinates of @a mesh by @a offset.
inline void translateExplicitCoordsX(conduit::Node& mesh, double offset)
{
  conduit::Node& x_node = mesh.fetch_existing("coordsets/coords/values/x");
  auto* x = x_node.as_float64_ptr();
  for(conduit::index_t i = 0; i < x_node.dtype().number_of_elements(); ++i)
  {
    x[i] += offset;
  }
}

//---------------------------------------------------------------------------
// Analytic fields
//
// Fields take (x, y, z). A 2D caller passes z == 0, so a RoundField used as a circle
// must also have center z == 0. PlanarField and RoundField use Axom's primal types.
//---------------------------------------------------------------------------

/*!
 * @brief Signed distance to a plane.
 *
 * An axis-aligned plane has no ambiguous cells, so its total measure can be compared directly.
 */
struct PlanarField
{
  using PointType = axom::primal::Point<double, 3>;
  using VectorType = axom::primal::Vector<double, 3>;
  using PlaneType = axom::primal::Plane<double, 3>;

  //! @brief The plane through @a origin with the given @a normal.
  PlanarField(const PointType& origin, const VectorType& normal) : plane(normal, origin) { }

  //! @brief The plane {p : normal . p == offset}.  @a normal need not be unit.
  PlanarField(const VectorType& normal, double offset) : plane(normal, offset) { }

  double operator()(double x, double y, double z) const
  {
    return plane.signedDistance(PointType {x, y, z});
  }

  PlaneType plane;
};

//! @brief Signed distance to a sphere (3D) or circle (2D, with center z == 0).
struct RoundField
{
  using PointType = axom::primal::Point<double, 3>;
  using SphereType = axom::primal::Sphere<double, 3>;

  RoundField(const PointType& center, double radius) : sphere(center, radius) { }

  double operator()(double x, double y, double z) const
  {
    return sphere.computeSignedDistance(PointType {x, y, z});
  }

  SphereType sphere;
};

/*!
 * @brief A gyroid.
 *
 * Its curvature produces non-planar cut polygons for triangulation tests.
 */
struct GyroidField
{
  using PointType = axom::primal::Point<double, 3>;

  explicit GyroidField(double scale_) : scale {scale_, scale_, scale_} { }
  explicit GyroidField(const PointType& scale_) : scale(scale_) { }

  double operator()(double x, double y, double z) const
  {
    const double sx = scale[0] * x, sy = scale[1] * y, sz = scale[2] * z;
    return std::sin(sx) * std::cos(sy) + std::sin(sy) * std::cos(sz) + std::sin(sz) * std::cos(sx);
  }

  PointType scale;
};

//---------------------------------------------------------------------------
// Coordinate warps
//---------------------------------------------------------------------------

//! @brief Identity warp (no displacement); the default coordinate map.
struct NoWarp
{
  void operator()(double& /*x*/, double& /*y*/, double& /*z*/) const { }
};

/*!
 * @brief A smooth sinusoidal coordinate warp.
 *
 * The displacement vanishes on the outer edges and corners.
 * Tests use a small amplitude so the hexahedra remain valid.
 */
struct SinusoidalWarp
{
  double amp;

  void operator()(double& x, double& y, double& z) const
  {
    const double sx = std::sin(M_PI * x), sy = std::sin(M_PI * y), sz = std::sin(M_PI * z);
    x += amp * sy * sz;
    y += amp * sx * sz;
    z += amp * sx * sy;
  }
};

//---------------------------------------------------------------------------
// Mesh builders
//---------------------------------------------------------------------------

/*!
 * @brief Sample @a f at the nodes of an explicit coordset and store it as a vertex field.
 *
 * The field copies the coordset's layout. When the topology has \c elements/dims/offsets
 * and \c strides, the field receives the same \c offsets and \c strides.
 */
template <int DIM, typename Field>
void addVertexField(conduit::Node& mesh, const Field& f, const std::string& fieldName)
{
  static_assert(DIM == 2 || DIM == 3, "DIM must be 2 or 3");

  const conduit::Node& values = mesh.fetch_existing("coordsets/coords/values");
  const auto x = values.fetch_existing("x").as_double_accessor();
  const auto y = values.fetch_existing("y").as_double_accessor();
  const auto z =
    DIM == 3 ? values.fetch_existing("z").as_double_accessor() : conduit::double_accessor {};

  conduit::Node& field = mesh["fields/" + fieldName];
  field["topology"] = "mesh";
  field["association"] = "vertex";
  field["values"].set(conduit::DataType::float64(x.number_of_elements()));
  auto* field_values = field["values"].as_float64_ptr();

  for(conduit::index_t idx = 0; idx < x.number_of_elements(); ++idx)
  {
    const double pz = DIM == 3 ? static_cast<double>(z[idx]) : 0.0;
    field_values[idx] = f(x[idx], y[idx], pz);
  }

  const conduit::Node& dims = mesh.fetch_existing("topologies/mesh/elements/dims");
  for(const char* key : {"offsets", "strides"})
  {
    if(dims.has_child(key))
    {
      field[key].set(dims.fetch_existing(key));
    }
  }
}

/*!
 * @brief Add an element-associated \c int32 field whose value is \a g(i, j, k).
 *
 * @a g receives the logical indices of each real (non-ghost) cell.
 * The field is stored compactly, without padding, in i-fastest order.
 */
template <int DIM, typename CellFunction>
void addCellField(conduit::Node& mesh, const CellFunction& g, const std::string& fieldName)
{
  static_assert(DIM == 2 || DIM == 3, "DIM must be 2 or 3");

  const conduit::Node& dims = mesh.fetch_existing("topologies/mesh/elements/dims");
  const int ni = dims.fetch_existing("i").to_int();
  const int nj = dims.fetch_existing("j").to_int();
  const int nk = (DIM == 3) ? dims.fetch_existing("k").to_int() : 1;

  conduit::Node& field = mesh["fields/" + fieldName];
  field["topology"] = "mesh";
  field["association"] = "element";
  field["values"].set(conduit::DataType::int32(conduit::index_t(ni) * nj * nk));
  auto* field_values = field["values"].as_int32_ptr();

  conduit::index_t idx = 0;
  for(int k = 0; k < nk; ++k)
  {
    for(int j = 0; j < nj; ++j)
    {
      for(int i = 0; i < ni; ++i, ++idx)
      {
        field_values[idx] = g(i, j, k);
      }
    }
  }
}

/*!
 * @brief Build a single-domain structured mesh with an explicit coordset on
 *        [0,1]^DIM with @a n cells per side.
 *
 * Nodes use i-fastest order to match Bump's StructuredIndexing.
 * If @a warp is supplied, the function is evaluated at the warped coordinates.
 */
template <int DIM, typename Field, typename Warp = NoWarp>
void buildStructured(conduit::Node& mesh,
                     int n,
                     const Field& f,
                     const std::string& fieldName,
                     const Warp& warp = Warp {})
{
  static_assert(DIM == 2 || DIM == 3, "DIM must be 2 or 3");
  const int nn = n + 1;
  conduit::index_t N = static_cast<conduit::index_t>(nn) * nn;
  if(DIM == 3)
  {
    N *= nn;
  }

  mesh.reset();

  conduit::Node& cs = mesh["coordsets/coords"];
  cs["type"] = "explicit";
  cs["values/x"].set(conduit::DataType::float64(N));
  cs["values/y"].set(conduit::DataType::float64(N));
  auto* x = cs["values/x"].as_float64_ptr();
  auto* y = cs["values/y"].as_float64_ptr();
  double* z = nullptr;
  if(DIM == 3)
  {
    cs["values/z"].set(conduit::DataType::float64(N));
    z = cs["values/z"].as_float64_ptr();
  }

  conduit::Node& topo = mesh["topologies/mesh"];
  topo["type"] = "structured";
  topo["coordset"] = "coords";
  topo["elements/dims/i"] = n;
  topo["elements/dims/j"] = n;
  if(DIM == 3)
  {
    topo["elements/dims/k"] = n;
  }

  const int nk = (DIM == 3) ? nn : 1;
  conduit::index_t idx = 0;
  for(int k = 0; k < nk; ++k)
  {
    for(int j = 0; j < nn; ++j)
    {
      for(int i = 0; i < nn; ++i, ++idx)
      {
        double px = double(i) / n, py = double(j) / n, pz = (DIM == 3) ? double(k) / n : 0.0;
        warp(px, py, pz);
        x[idx] = px;
        y[idx] = py;
        if(z != nullptr)
        {
          z[idx] = pz;
        }
      }
    }
  }

  addVertexField<DIM>(mesh, f, fieldName);
}

/*!
 * @brief Build a strided structured mesh on [0,1]^DIM with @a n real cells per side
 *        and @a pad ghost nodes on each side of every axis.
 *
 * The coordset and vertex field store (n + 1 + 2*pad)^DIM nodes in i-fastest order.
 * The topology's \c elements/dims/offsets and \c strides select the real nodes.
 * Ghost nodes continue the uniform spacing outside [0,1],
 * so an isocontour that ignores the offsets would include ghost cells.
 */
template <int DIM, typename Field>
void buildStridedStructured(conduit::Node& mesh,
                            int n,
                            int pad,
                            const Field& f,
                            const std::string& fieldName)
{
  static_assert(DIM == 2 || DIM == 3, "DIM must be 2 or 3");
  const int np = n + 1 + 2 * pad;
  conduit::index_t N = static_cast<conduit::index_t>(np) * np;
  if(DIM == 3)
  {
    N *= np;
  }

  mesh.reset();

  conduit::Node& cs = mesh["coordsets/coords"];
  cs["type"] = "explicit";
  cs["values/x"].set(conduit::DataType::float64(N));
  cs["values/y"].set(conduit::DataType::float64(N));
  auto* x = cs["values/x"].as_float64_ptr();
  auto* y = cs["values/y"].as_float64_ptr();
  double* z = nullptr;
  if(DIM == 3)
  {
    cs["values/z"].set(conduit::DataType::float64(N));
    z = cs["values/z"].as_float64_ptr();
  }

  conduit::Node& topo = mesh["topologies/mesh"];
  topo["type"] = "structured";
  topo["coordset"] = "coords";
  topo["elements/dims/i"] = n;
  topo["elements/dims/j"] = n;
  if(DIM == 3)
  {
    topo["elements/dims/k"] = n;
  }

  topo["elements/dims/offsets"].set(conduit::DataType::int32(DIM));
  topo["elements/dims/strides"].set(conduit::DataType::int32(DIM));
  auto* offsets = topo["elements/dims/offsets"].as_int32_ptr();
  auto* strides = topo["elements/dims/strides"].as_int32_ptr();
  int stride = 1;
  for(int d = 0; d < DIM; ++d)
  {
    offsets[d] = pad;
    strides[d] = stride;
    stride *= np;
  }

  const int nk = (DIM == 3) ? np : 1;
  conduit::index_t idx = 0;
  for(int k = 0; k < nk; ++k)
  {
    for(int j = 0; j < np; ++j)
    {
      for(int i = 0; i < np; ++i, ++idx)
      {
        x[idx] = double(i - pad) / n;
        y[idx] = double(j - pad) / n;
        if(z != nullptr)
        {
          z[idx] = double(k - pad) / n;
        }
      }
    }
  }

  addVertexField<DIM>(mesh, f, fieldName);
}

}  // namespace axom::quest::testing::marching_cubes
