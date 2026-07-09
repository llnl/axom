// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_marching_cubes_ray_queries.cpp
 *
 * \brief Builds several 3D marching-cubes test cases directly in memory and
 * checks whether selected rays intersect the generated contour surfaces.
 *
 * This file contains two kinds of coverage:
 *  - a banded structured-mesh case with a nodalized volume-fraction field that
 *    captures the known z-directed query failure
 *  - smaller analytic cases that exercise baseline MarchingCubes behavior on
 *    uniform, warped, planar, and multidomain meshes
 */

#include "axom/config.hpp"

#ifdef AXOM_USE_CONDUIT

  #include "axom/core.hpp"
  #include "axom/fmt.hpp"
  #include "axom/mint.hpp"
  #include "axom/primal.hpp"
  #include "axom/quest/MarchingCubes.hpp"

  #include "conduit_blueprint.hpp"

  #include "gtest/gtest.h"

  #include <array>
  #include <cmath>
  #include <string>
  #include <vector>

namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;

using Point3D = primal::Point<double, 3>;
using Triangle3D = primal::Triangle<double, 3>;
using Vector3D = primal::Vector<double, 3>;
using Ray3D = primal::Ray<double, 3>;
using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

namespace
{

//------------------------------------------------------------------------------
// Test geometry constants
//------------------------------------------------------------------------------
//
// These dimensions and per-k x-band values define a structured mesh with:
//   - 43 x 199 x 3 hexahedral zones
//   - 44 x 200 x 4 vertices
//
// The mesh is intentionally awkward in x: a short left band, a one-cell target
// band, and a long right band. The failing z-rays pass through the sphere in
// this geometry and expose the surface hole in the generated contour.
//
constexpr axom::IndexType NUM_ELEMS_I = 43;
constexpr axom::IndexType NUM_ELEMS_J = 199;
constexpr axom::IndexType NUM_ELEMS_K = 3;

constexpr axom::IndexType NUM_VERTS_I = NUM_ELEMS_I + 1;
constexpr axom::IndexType NUM_VERTS_J = NUM_ELEMS_J + 1;
constexpr axom::IndexType NUM_VERTS_K = NUM_ELEMS_K + 1;

constexpr int VF_SAMPLES_PER_DIM = 8;

// Sphere used by the test case.
const Point3D SPHERE_CENTER {1.5, 0.5, 0.5};
constexpr double SPHERE_RADIUS = 0.45;
constexpr double SPHERE_RADIUS_SQUARED = SPHERE_RADIUS * SPHERE_RADIUS;

constexpr std::array<double, NUM_VERTS_K> Z_VALUES {{0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0}};
constexpr std::array<double, NUM_VERTS_K> X0_VALUES {
  {-0.17, -0.16333333333333333, -0.15666666666666668, -0.15}};
constexpr std::array<double, NUM_VERTS_K> X3_VALUES {
  {-0.01, -0.0033333333333333335, 0.003333333333333333, 0.01}};
constexpr std::array<double, NUM_VERTS_K> X4_VALUES {
  {0.0374, 0.044066666666666664, 0.050733333333333332, 0.0574}};
constexpr std::array<double, NUM_VERTS_K> X5_VALUES {
  {5.0374, 5.044066666666667, 5.050733333333333, 5.0574}};

struct RayCase
{
  const char* name;
  Point3D start;
  Point3D end;
};

struct BlueprintMeshData
{
  // The Conduit mesh stores external pointers into the vectors below, so the
  // owning storage stays in the same object for the lifetime of the contouring
  // call.
  conduit::Node mesh;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> vf;
};

struct StructuredDomainSpec
{
  axom::IndexType ni;
  axom::IndexType nj;
  axom::IndexType nk;
  axom::IndexType domainId;
};

struct MultiDomainBlueprintData
{
  struct DomainArrays
  {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> field;
  };

  conduit::Node mesh;
  std::vector<DomainArrays> domains;
};

using HexVertexIds = std::array<axom::IndexType, 8>;

//! Flatten structured (i,j,k) vertex indices into the explicit coord arrays.
axom::IndexType vertexIndex(axom::IndexType i, axom::IndexType j, axom::IndexType k)
{
  return (k * NUM_VERTS_J + j) * NUM_VERTS_I + i;
}

//! Return the standard hexahedron vertex ordering for a structured cell.
HexVertexIds hexVertexIds(axom::IndexType i, axom::IndexType j, axom::IndexType k)
{
  return {{vertexIndex(i, j, k),
           vertexIndex(i + 1, j, k),
           vertexIndex(i + 1, j + 1, k),
           vertexIndex(i, j + 1, k),
           vertexIndex(i, j, k + 1),
           vertexIndex(i + 1, j, k + 1),
           vertexIndex(i + 1, j + 1, k + 1),
           vertexIndex(i, j + 1, k + 1)}};
}

double interpolate(double a, double b, double t) { return (1.0 - t) * a + t * b; }

axom::IndexType vertexCount(axom::IndexType ni, axom::IndexType nj, axom::IndexType nk)
{
  return (ni + 1) * (nj + 1) * (nk + 1);
}

axom::IndexType vertexIndex(axom::IndexType i,
                            axom::IndexType j,
                            axom::IndexType k,
                            axom::IndexType ni,
                            axom::IndexType nj)
{
  return (k * (nj + 1) + j) * (ni + 1) + i;
}

Point3D affinePoint(const Point3D& origin,
                    const Vector3D& uVector,
                    const Vector3D& vVector,
                    const Vector3D& wVector,
                    double u,
                    double v,
                    double w)
{
  Point3D pt;
  for(int dim = 0; dim < 3; ++dim)
  {
    pt[dim] = origin[dim] + u * uVector[dim] + v * vVector[dim] + w * wVector[dim];
  }
  return pt;
}

double sphereLevelSet(const Point3D& pt, const Point3D& center, double radius)
{
  const double dx = pt[0] - center[0];
  const double dy = pt[1] - center[1];
  const double dz = pt[2] - center[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz) - radius;
}

double planeField(const Point3D& pt, const Vector3D& normal, double offset)
{
  return pt[0] * normal[0] + pt[1] * normal[1] + pt[2] * normal[2] - offset;
}

/*!
 * \brief Return the x-coordinate of vertex \a i in z-plane \a k.
 *
 * The mesh defines three x-bands in each k-plane:
 *   - i in [0,1]   : left band from x0(k) to x3(k)
 *   - i in [1,2]   : thin target band from x3(k) to x4(k)
 *   - i in [2,43]  : long right band from x4(k) to x5(k)
 *
 * We intentionally ignore the tiny y-dependent spline perturbation on x4 here.
 * The current setup already exposes the miss without it, so keeping the
 * geometry piecewise-linear makes the test easier to read and maintain.
 */
double bandedX(axom::IndexType i, axom::IndexType k)
{
  if(i <= 1)
  {
    return interpolate(X0_VALUES[k], X3_VALUES[k], static_cast<double>(i));
  }

  if(i <= 2)
  {
    return interpolate(X3_VALUES[k], X4_VALUES[k], static_cast<double>(i - 1));
  }

  const double t = static_cast<double>(i - 2) / static_cast<double>(NUM_ELEMS_I - 2);
  return interpolate(X4_VALUES[k], X5_VALUES[k], t);
}

/*!
 * \brief Evaluate a point inside a hexahedron using trilinear interpolation.
 *
 * We use this instead of assuming an axis-aligned box because the explicit
 * structured mesh is warped in x across k-planes. That detail matters for
 * reproducing the same sampled volume fractions handed to MarchingCubes.
 */
Point3D trilinearSample(const std::array<Point3D, 8>& corners, double u, double v, double w)
{
  Point3D sample;

  for(int dim = 0; dim < 3; ++dim)
  {
    const double c00 = interpolate(corners[0][dim], corners[1][dim], u);
    const double c10 = interpolate(corners[3][dim], corners[2][dim], u);
    const double c01 = interpolate(corners[4][dim], corners[5][dim], u);
    const double c11 = interpolate(corners[7][dim], corners[6][dim], u);

    const double c0 = interpolate(c00, c10, v);
    const double c1 = interpolate(c01, c11, v);
    sample[dim] = interpolate(c0, c1, w);
  }

  return sample;
}

/*!
 * \brief Estimate a zone-centered sphere volume fraction by deterministic
 * subcell sampling.
 *
 * This mirrors a common "estimate a per-zone volume fraction, then average it
 * onto the zone's eight vertices" handoff used before contouring. The test
 * does not need the exact application-side sampler implementation; it needs a
 * stable, deterministic approximation that produces the same problematic
 * contour.
 */
double estimateZoneVolumeFraction(const BlueprintMeshData& data,
                                  axom::IndexType i,
                                  axom::IndexType j,
                                  axom::IndexType k)
{
  const auto nodeIds = hexVertexIds(i, j, k);

  std::array<Point3D, 8> corners;
  for(int idx = 0; idx < 8; ++idx)
  {
    const auto vertexId = nodeIds[idx];
    corners[idx] = Point3D {data.x[vertexId], data.y[vertexId], data.z[vertexId]};
  }

  int insideCount = 0;
  for(int kk = 0; kk < VF_SAMPLES_PER_DIM; ++kk)
  {
    const double w = (static_cast<double>(kk) + 0.5) / VF_SAMPLES_PER_DIM;
    for(int jj = 0; jj < VF_SAMPLES_PER_DIM; ++jj)
    {
      const double v = (static_cast<double>(jj) + 0.5) / VF_SAMPLES_PER_DIM;
      for(int ii = 0; ii < VF_SAMPLES_PER_DIM; ++ii)
      {
        const double u = (static_cast<double>(ii) + 0.5) / VF_SAMPLES_PER_DIM;
        const Point3D sample = trilinearSample(corners, u, v, w);
        const double dx = sample[0] - SPHERE_CENTER[0];
        const double dy = sample[1] - SPHERE_CENTER[1];
        const double dz = sample[2] - SPHERE_CENTER[2];
        const double distSquared = dx * dx + dy * dy + dz * dz;
        if(distSquared <= SPHERE_RADIUS_SQUARED)
        {
          ++insideCount;
        }
      }
    }
  }

  constexpr double sampleCount = VF_SAMPLES_PER_DIM * VF_SAMPLES_PER_DIM * VF_SAMPLES_PER_DIM;
  return insideCount / sampleCount;
}

//------------------------------------------------------------------------------
// Blueprint mesh construction
//------------------------------------------------------------------------------

void populateCoordinates(BlueprintMeshData& data)
{
  const axom::IndexType vertexCount = NUM_VERTS_I * NUM_VERTS_J * NUM_VERTS_K;
  data.x.resize(vertexCount);
  data.y.resize(vertexCount);
  data.z.resize(vertexCount);

  for(axom::IndexType k = 0; k < NUM_VERTS_K; ++k)
  {
    for(axom::IndexType j = 0; j < NUM_VERTS_J; ++j)
    {
      const double y = static_cast<double>(j) / static_cast<double>(NUM_ELEMS_J);
      for(axom::IndexType i = 0; i < NUM_VERTS_I; ++i)
      {
        const axom::IndexType idx = vertexIndex(i, j, k);
        data.x[idx] = bandedX(i, k);
        data.y[idx] = y;
        data.z[idx] = Z_VALUES[k];
      }
    }
  }
}

void populateNodalVolumeFractions(BlueprintMeshData& data)
{
  const axom::IndexType vertexCount = NUM_VERTS_I * NUM_VERTS_J * NUM_VERTS_K;
  data.vf.resize(vertexCount, 0.0);

  std::vector<double> vfAccum(vertexCount, 0.0);
  std::vector<int> contributionCount(vertexCount, 0);

  for(axom::IndexType k = 0; k < NUM_ELEMS_K; ++k)
  {
    for(axom::IndexType j = 0; j < NUM_ELEMS_J; ++j)
    {
      for(axom::IndexType i = 0; i < NUM_ELEMS_I; ++i)
      {
        const double zoneVf = estimateZoneVolumeFraction(data, i, j, k);
        const auto nodeIds = hexVertexIds(i, j, k);

        for(const auto nodeId : nodeIds)
        {
          vfAccum[nodeId] += zoneVf;
          ++contributionCount[nodeId];
        }
      }
    }
  }

  for(axom::IndexType idx = 0; idx < vertexCount; ++idx)
  {
    EXPECT_GT(contributionCount[idx], 0);
    data.vf[idx] = vfAccum[idx] / static_cast<double>(contributionCount[idx]);
  }
}

void populateBlueprintNode(BlueprintMeshData& data)
{
  const axom::IndexType vertexCount = NUM_VERTS_I * NUM_VERTS_J * NUM_VERTS_K;

  conduit::Node& domain = data.mesh.append();
  domain["coordsets/coords/type"] = "explicit";
  domain["coordsets/coords/values/x"].set_external(conduit::DataType::float64(vertexCount),
                                                   data.x.data());
  domain["coordsets/coords/values/y"].set_external(conduit::DataType::float64(vertexCount),
                                                   data.y.data());
  domain["coordsets/coords/values/z"].set_external(conduit::DataType::float64(vertexCount),
                                                   data.z.data());

  domain["topologies/mesh/type"] = "structured";
  domain["topologies/mesh/coordset"] = "coords";
  domain["topologies/mesh/elements/dims/i"] = NUM_ELEMS_I;
  domain["topologies/mesh/elements/dims/j"] = NUM_ELEMS_J;
  domain["topologies/mesh/elements/dims/k"] = NUM_ELEMS_K;
  domain["state/domain_id"] = 0;

  domain["fields/vf/association"] = "vertex";
  domain["fields/vf/topology"] = "mesh";
  domain["fields/vf/volume_dependent"] = "false";
  domain["fields/vf/values"].set_external(conduit::DataType::float64(vertexCount), data.vf.data());
}

template <typename CoordFunctor, typename FieldFunctor>
void appendStructuredDomain(MultiDomainBlueprintData& data,
                            const StructuredDomainSpec& spec,
                            const CoordFunctor& coordFunctor,
                            const FieldFunctor& fieldFunctor,
                            const std::string& fieldName)
{
  data.domains.emplace_back();
  auto& arrays = data.domains.back();

  const axom::IndexType nVerts = vertexCount(spec.ni, spec.nj, spec.nk);
  arrays.x.resize(nVerts);
  arrays.y.resize(nVerts);
  arrays.z.resize(nVerts);
  arrays.field.resize(nVerts);

  for(axom::IndexType k = 0; k <= spec.nk; ++k)
  {
    const double w = spec.nk > 0 ? static_cast<double>(k) / spec.nk : 0.0;
    for(axom::IndexType j = 0; j <= spec.nj; ++j)
    {
      const double v = spec.nj > 0 ? static_cast<double>(j) / spec.nj : 0.0;
      for(axom::IndexType i = 0; i <= spec.ni; ++i)
      {
        const double u = spec.ni > 0 ? static_cast<double>(i) / spec.ni : 0.0;
        const axom::IndexType idx = vertexIndex(i, j, k, spec.ni, spec.nj);
        const Point3D pt = coordFunctor(u, v, w);

        arrays.x[idx] = pt[0];
        arrays.y[idx] = pt[1];
        arrays.z[idx] = pt[2];
        arrays.field[idx] = fieldFunctor(pt);
      }
    }
  }

  conduit::Node& domain = data.mesh.append();
  domain["coordsets/coords/type"] = "explicit";
  domain["coordsets/coords/values/x"].set_external(conduit::DataType::float64(nVerts),
                                                   arrays.x.data());
  domain["coordsets/coords/values/y"].set_external(conduit::DataType::float64(nVerts),
                                                   arrays.y.data());
  domain["coordsets/coords/values/z"].set_external(conduit::DataType::float64(nVerts),
                                                   arrays.z.data());

  domain["topologies/mesh/type"] = "structured";
  domain["topologies/mesh/coordset"] = "coords";
  domain["topologies/mesh/elements/dims/i"] = spec.ni;
  domain["topologies/mesh/elements/dims/j"] = spec.nj;
  domain["topologies/mesh/elements/dims/k"] = spec.nk;
  domain["state/domain_id"] = spec.domainId;

  domain["fields/" + fieldName + "/association"] = "vertex";
  domain["fields/" + fieldName + "/topology"] = "mesh";
  domain["fields/" + fieldName + "/volume_dependent"] = "false";
  domain["fields/" + fieldName + "/values"].set_external(conduit::DataType::float64(nVerts),
                                                         arrays.field.data());
}

/*!
 * \brief Build the single-domain Blueprint mesh and its nodal "vf" field.
 *
 * The field construction follows the nodalization pattern:
 *   1. compute one zone-centered sphere volume fraction per hex
 *   2. add that value to each of the zone's eight vertices
 *   3. divide each vertex sum by the number of contributing zones
 *
 * That nodal averaging step is important. MarchingCubes contours a vertex
 * field, and this test is specifically about the surface produced from
 * the nodalized volume fractions rather than a direct analytic signed
 * distance or exact sphere indicator.
 */
void buildBlueprintMesh(BlueprintMeshData& data)
{
  // Step 1: fill the explicit structured coordinates.
  populateCoordinates(data);

  // Step 2: sample each zone and finish the MaterialContour-style nodal
  // averaging of those zone-centered values.
  populateNodalVolumeFractions(data);

  // Step 3: build a single-domain explicit-coordinate structured Blueprint
  // mesh from the populated coordinate and field buffers.
  // We use `set_external()` so MarchingCubes reads the same in-memory buffers
  // we just populated.
  populateBlueprintNode(data);
}

/*!
 * \brief Generate the contour mesh from the constructed Blueprint input.
 *
 * The construction is intentionally simple and matches the requested call
 * sequence exactly so the test stays focused on one MarchingCubes path rather
 * than testing alternative setup variations.
 */
void buildContourMesh(const conduit::Node& blueprintMesh,
                      const std::string& fieldName,
                      double contourValue,
                      UMesh& contourMesh)
{
  conduit::Node verifyInfo;
  EXPECT_TRUE(conduit::blueprint::mesh::verify(blueprintMesh, verifyInfo)) << verifyInfo.to_yaml();

  quest::MarchingCubes mc(axom::runtime_policy::Policy::seq,
                          axom::getDefaultAllocatorID(),
                          quest::MarchingCubesDataParallelism::byPolicy);
  mc.setMesh(blueprintMesh, "mesh");
  mc.setFunctionField(fieldName);
  mc.computeIsocontour(contourValue);
  mc.populateContourMesh(contourMesh, "parent_cell", "domain_id");
}

void buildBandedVolumeFractionContourMesh(UMesh& contourMesh)
{
  BlueprintMeshData blueprintMesh;
  buildBlueprintMesh(blueprintMesh);

  buildContourMesh(blueprintMesh.mesh, "vf", 0.5, contourMesh);
}

/*!
 * \brief Check whether a finite ray segment intersects any generated triangle.
 *
 * We deliberately use ray semantics, not segment/triangle intersection,
 * because these tests are about ray queries against the generated surface. The
 * additional `t` bound trims the infinite ray back to the user-supplied start
 * and end points. The `t > eps` check avoids counting an immediate self-hit at
 * the ray origin, which would make the challenging z-directed cases pass for
 * the wrong reason.
 */
bool rayHitsContour(const UMesh& contourMesh, const RayCase& ray)
{
  const Vector3D direction(ray.start, ray.end);
  const double maxT = direction.norm();
  const Ray3D queryRay(ray.start, direction);
  constexpr double eps = 1.0e-12;
  double coords[3] {};

  for(axom::IndexType cellId = 0; cellId < contourMesh.getNumberOfCells(); ++cellId)
  {
    const axom::IndexType* nodeIds = contourMesh.getCellNodeIDs(cellId);
    Triangle3D tri;

    for(int localId = 0; localId < 3; ++localId)
    {
      contourMesh.getNode(nodeIds[localId], coords);
      tri[localId] = Point3D {coords[0], coords[1], coords[2]};
    }

    double t = 0.0;
    if(primal::intersect(tri, queryRay, t) && t > eps && t <= maxT + eps)
    {
      return true;
    }
  }

  return false;
}

std::string rayMessage(const RayCase& ray)
{
  return axom::fmt::format("Expected ray '{}' {} -> {} to intersect the generated contour surface.",
                           ray.name,
                           ray.start,
                           ray.end);
}

template <typename RayContainer>
void expectRaysHit(const UMesh& contourMesh, const RayContainer& rays)
{
  for(const auto& ray : rays)
  {
    EXPECT_TRUE(rayHitsContour(contourMesh, ray)) << rayMessage(ray);
  }
}

const std::array<RayCase, 3>& controlRays()
{
  static const std::array<RayCase, 3> rays {{
    {"control_from_x", Point3D {5.0474, 0.5, 0.5}, Point3D {0.0, 0.5, 0.5}},
    {"control_from_y", Point3D {1.5, 1.0, 0.5}, Point3D {1.5, 0.0, 0.5}},
    {"control_diag", Point3D {2.25, 1.0, 1.0}, Point3D {0.75, 0.0, 0.0}},
  }};
  return rays;
}

const std::array<RayCase, 7>& challengingZRays()
{
  static const std::array<RayCase, 7> rays {{
    {"from_z", Point3D {1.5, 0.5, 1.0}, Point3D {1.5, 0.5, 0.0}},
    {"from_z_inset", Point3D {1.5, 0.5, 0.99}, Point3D {1.5, 0.5, 0.0}},
    {"from_z_xoff", Point3D {1.55, 0.5, 1.0}, Point3D {1.55, 0.5, 0.0}},
    {"from_z_yoff", Point3D {1.5, 0.55, 1.0}, Point3D {1.5, 0.55, 0.0}},
    {"from_z_rev", Point3D {1.5, 0.5, 0.0}, Point3D {1.5, 0.5, 1.0}},
    {"from_z_tiltx", Point3D {1.45, 0.5, 1.0}, Point3D {1.55, 0.5, 0.0}},
    {"from_z_tilty", Point3D {1.5, 0.45, 1.0}, Point3D {1.5, 0.55, 0.0}},
  }};
  return rays;
}

const std::array<RayCase, 4>& uniformSphereRays()
{
  static const std::array<RayCase, 4> rays {{
    {"x_center", Point3D {3.0, 0.5, 0.5}, Point3D {0.0, 0.5, 0.5}},
    {"y_center", Point3D {1.5, 1.0, 0.5}, Point3D {1.5, 0.0, 0.5}},
    {"z_center", Point3D {1.5, 0.5, 1.0}, Point3D {1.5, 0.5, 0.0}},
    {"diag", Point3D {2.1, 1.1, 1.1}, Point3D {0.9, -0.1, -0.1}},
  }};
  return rays;
}

const std::array<RayCase, 3>& warpedSphereRays()
{
  static const std::array<RayCase, 3> rays {{
    {"x_center", Point3D {2.8, 0.55, 0.5}, Point3D {0.2, 0.55, 0.5}},
    {"y_center", Point3D {1.55, 1.0, 0.5}, Point3D {1.55, 0.0, 0.5}},
    {"z_center", Point3D {1.55, 0.55, 1.0}, Point3D {1.55, 0.55, 0.0}},
  }};
  return rays;
}

const std::array<RayCase, 3>& planeCutRays()
{
  static const std::array<RayCase, 3> rays {{
    {"from_x", Point3D {1.0, 0.5, 0.5}, Point3D {0.0, 0.5, 0.5}},
    {"from_y", Point3D {0.4, 1.0, 0.5}, Point3D {0.4, 0.0, 0.5}},
    {"diag", Point3D {1.0, 1.0, 0.5}, Point3D {0.0, 0.0, 0.5}},
  }};
  return rays;
}

const std::array<RayCase, 3>& seamSphereRays()
{
  static const std::array<RayCase, 3> rays {{
    {"across_x", Point3D {3.0, 0.5, 0.5}, Point3D {0.0, 0.5, 0.5}},
    {"through_seam_y", Point3D {1.5, 1.0, 0.5}, Point3D {1.5, 0.0, 0.5}},
    {"through_seam_z", Point3D {1.5, 0.5, 1.0}, Point3D {1.5, 0.5, 0.0}},
  }};
  return rays;
}

MultiDomainBlueprintData makeUniformSphereBlueprint()
{
  MultiDomainBlueprintData data;
  const StructuredDomainSpec spec {32, 24, 24, 0};

  const auto coordFunctor = [](double u, double v, double w) { return Point3D {3.0 * u, v, w}; };
  const auto fieldFunctor = [](const Point3D& pt) {
    return sphereLevelSet(pt, SPHERE_CENTER, SPHERE_RADIUS);
  };

  appendStructuredDomain(data, spec, coordFunctor, fieldFunctor, "field");
  return data;
}

MultiDomainBlueprintData makeWarpedSphereBlueprint()
{
  MultiDomainBlueprintData data;
  const StructuredDomainSpec spec {48, 12, 6, 0};
  const Point3D origin {0.0, 0.0, 0.0};
  const Vector3D uVector {2.9, 0.0, 0.0};
  const Vector3D vVector {0.15, 1.0, 0.0};
  const Vector3D wVector {0.05, 0.05, 1.0};
  const Point3D sphereCenter {1.55, 0.55, 0.5};
  constexpr double warpedRadius = 0.35;

  const auto coordFunctor = [=](double u, double v, double w) {
    return affinePoint(origin, uVector, vVector, wVector, u, v, w);
  };
  const auto fieldFunctor = [=](const Point3D& pt) {
    return sphereLevelSet(pt, sphereCenter, warpedRadius);
  };

  appendStructuredDomain(data, spec, coordFunctor, fieldFunctor, "field");
  return data;
}

MultiDomainBlueprintData makePlaneCutBlueprint()
{
  MultiDomainBlueprintData data;
  const StructuredDomainSpec spec {20, 20, 8, 0};
  const Vector3D normal {1.0, 0.25, 0.0};

  const auto coordFunctor = [](double u, double v, double w) { return Point3D {u, v, w}; };
  const auto fieldFunctor = [=](const Point3D& pt) { return planeField(pt, normal, 0.5); };

  appendStructuredDomain(data, spec, coordFunctor, fieldFunctor, "field");
  return data;
}

MultiDomainBlueprintData makeSeamSphereBlueprint()
{
  MultiDomainBlueprintData data;
  const StructuredDomainSpec leftSpec {16, 24, 24, 0};
  const StructuredDomainSpec rightSpec {16, 24, 24, 1};

  const auto leftCoords = [](double u, double v, double w) { return Point3D {1.5 * u, v, w}; };
  const auto rightCoords = [](double u, double v, double w) { return Point3D {1.5 + 1.5 * u, v, w}; };
  const auto fieldFunctor = [](const Point3D& pt) {
    return sphereLevelSet(pt, SPHERE_CENTER, SPHERE_RADIUS);
  };

  appendStructuredDomain(data, leftSpec, leftCoords, fieldFunctor, "field");
  appendStructuredDomain(data, rightSpec, rightCoords, fieldFunctor, "field");
  return data;
}

}  // namespace

//------------------------------------------------------------------------------
// Tests
//------------------------------------------------------------------------------

TEST(quest_marching_cubes_ray_queries, banded_nodal_volume_fraction_surface)
{
  // First prove the generated contour exists and is queryable at all.
  // If these controls fail, the challenging z-directed rays are not
  // informative.
  UMesh contourMesh(3, mint::TRIANGLE);
  buildBandedVolumeFractionContourMesh(contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  ASSERT_GT(contourMesh.getNumberOfNodes(), 0);

  // These control rays prove the contoured sphere exists and that the basic
  // ray-query setup is working before we check the more demanding z-directed
  // cases.
  expectRaysHit(contourMesh, controlRays());

  // These are the cases for the current z-directed miss. They
  // should hit the sphere surface, but currently miss because the contoured
  // mesh has a hole near the z-directed path through the sphere.
  // Leave these as normal expectations so the test registers as a failing
  // test until the MarchingCubes hole is fixed.
  expectRaysHit(contourMesh, challengingZRays());
}

TEST(quest_marching_cubes_ray_queries, uniform_sphere_surface)
{
  // Baseline closed-surface case on a uniform structured mesh. If this fails,
  // MarchingCubes is broken in a much more general way than the banded mesh
  // case above.
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeUniformSphereBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  expectRaysHit(contourMesh, uniformSphereRays());
}

TEST(quest_marching_cubes_ray_queries, warped_sphere_surface)
{
  // Similar sphere test on an explicit but skewed mesh. This keeps the field
  // analytic while checking that mild geometric warping does not prevent the
  // generated surface from answering simple ray queries.
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeWarpedSphereBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  expectRaysHit(contourMesh, warpedSphereRays());
}

TEST(quest_marching_cubes_ray_queries, planar_surface)
{
  // An open planar cut is a simpler topology than a closed sphere. It is a
  // useful control because the expected surface geometry is easy to reason
  // about and should be robust on a regular grid.
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makePlaneCutBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  expectRaysHit(contourMesh, planeCutRays());
}

TEST(quest_marching_cubes_ray_queries, multidomain_seam_sphere_surface)
{
  // The sphere is split across two domains that meet on a shared x-normal
  // interface. This checks that MarchingCubes can still generate a hittable
  // surface when the sampled field crosses a domain seam.
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeSeamSphereBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  expectRaysHit(contourMesh, seamSphereRays());
}

#endif
