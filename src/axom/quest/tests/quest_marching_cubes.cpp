// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_marching_cubes.cpp
 *
 * \brief Consolidated MarchingCubes coverage for the known banded-mesh
 * reproducer and structured reference cases with exact contour assertions.
 */

#include "axom/config.hpp"

#ifdef AXOM_USE_CONDUIT

  #include "axom/core.hpp"
  #include "axom/fmt.hpp"
  #include "axom/mint.hpp"
  #include "axom/primal.hpp"
  #include "axom/quest/MarchingCubes.hpp"
  #define _MC_LOOKUP_CASES3D
  #define _MC_LOOKUP_NUM_TRIANGLES
  #include "axom/quest/detail/marching_cubes_lookup.hpp"
  #undef _MC_LOOKUP_NUM_TRIANGLES
  #undef _MC_LOOKUP_CASES3D

  #include "conduit_blueprint.hpp"

  #include "gtest/gtest.h"

  #include <algorithm>
  #include <array>
  #include <cmath>
  #include <map>
  #include <set>
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
using DomainIdType = quest::MarchingCubes::DomainIdType;

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
  DomainIdType domainId;
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

struct MaskedSingleCubeBlueprintData
{
  conduit::Node mesh;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> field;
  std::vector<int> mask;
};

using HexVertexIds = std::array<axom::IndexType, 8>;
using EdgeVertexIds = std::array<int, 2>;

constexpr StructuredDomainSpec SINGLE_CUBE_SPEC {2, 2, 2, 0};
constexpr double POINT_TOL = 1.0e-12;

axom::IndexType vertexIndex(axom::IndexType i, axom::IndexType j, axom::IndexType k)
{
  return (k * NUM_VERTS_J + j) * NUM_VERTS_I + i;
}

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

axom::IndexType cellCount(axom::IndexType ni, axom::IndexType nj, axom::IndexType nk)
{
  return ni * nj * nk;
}

axom::IndexType vertexIndex(axom::IndexType i,
                            axom::IndexType j,
                            axom::IndexType k,
                            axom::IndexType ni,
                            axom::IndexType nj)
{
  return (k * (nj + 1) + j) * (ni + 1) + i;
}

axom::IndexType cellIndex(axom::IndexType i,
                          axom::IndexType j,
                          axom::IndexType k,
                          axom::IndexType ni,
                          axom::IndexType nj)
{
  return (k * nj + j) * ni + i;
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
  const axom::IndexType nVerts = NUM_VERTS_I * NUM_VERTS_J * NUM_VERTS_K;
  data.x.resize(nVerts);
  data.y.resize(nVerts);
  data.z.resize(nVerts);

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
  const axom::IndexType nVerts = NUM_VERTS_I * NUM_VERTS_J * NUM_VERTS_K;
  data.vf.resize(nVerts, 0.0);

  std::vector<double> vfAccum(nVerts, 0.0);
  std::vector<int> contributionCount(nVerts, 0);

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

  for(axom::IndexType idx = 0; idx < nVerts; ++idx)
  {
    EXPECT_GT(contributionCount[idx], 0);
    data.vf[idx] = vfAccum[idx] / static_cast<double>(contributionCount[idx]);
  }
}

void populateBlueprintNode(BlueprintMeshData& data)
{
  const axom::IndexType nVerts = NUM_VERTS_I * NUM_VERTS_J * NUM_VERTS_K;

  conduit::Node& domain = data.mesh.append();
  domain["coordsets/coords/type"] = "explicit";
  domain["coordsets/coords/values/x"].set_external(conduit::DataType::float64(nVerts), data.x.data());
  domain["coordsets/coords/values/y"].set_external(conduit::DataType::float64(nVerts), data.y.data());
  domain["coordsets/coords/values/z"].set_external(conduit::DataType::float64(nVerts), data.z.data());

  domain["topologies/mesh/type"] = "structured";
  domain["topologies/mesh/coordset"] = "coords";
  domain["topologies/mesh/elements/dims/i"] = NUM_ELEMS_I;
  domain["topologies/mesh/elements/dims/j"] = NUM_ELEMS_J;
  domain["topologies/mesh/elements/dims/k"] = NUM_ELEMS_K;
  domain["state/domain_id"] = 0;

  domain["fields/vf/association"] = "vertex";
  domain["fields/vf/topology"] = "mesh";
  domain["fields/vf/volume_dependent"] = "false";
  domain["fields/vf/values"].set_external(conduit::DataType::float64(nVerts), data.vf.data());
}

void appendStructuredDomainData(MultiDomainBlueprintData& data,
                                const StructuredDomainSpec& spec,
                                std::vector<double> x,
                                std::vector<double> y,
                                std::vector<double> z,
                                std::vector<double> field,
                                const std::string& fieldName)
{
  const axom::IndexType nVerts = vertexCount(spec.ni, spec.nj, spec.nk);
  EXPECT_EQ(static_cast<axom::IndexType>(x.size()), nVerts);
  EXPECT_EQ(static_cast<axom::IndexType>(y.size()), nVerts);
  EXPECT_EQ(static_cast<axom::IndexType>(z.size()), nVerts);
  EXPECT_EQ(static_cast<axom::IndexType>(field.size()), nVerts);

  data.domains.push_back({std::move(x), std::move(y), std::move(z), std::move(field)});
  auto& arrays = data.domains.back();

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

template <typename CoordFunctor, typename FieldFunctor>
void appendStructuredDomain(MultiDomainBlueprintData& data,
                            const StructuredDomainSpec& spec,
                            const CoordFunctor& coordFunctor,
                            const FieldFunctor& fieldFunctor,
                            const std::string& fieldName)
{
  const axom::IndexType nVerts = vertexCount(spec.ni, spec.nj, spec.nk);
  std::vector<double> x(nVerts);
  std::vector<double> y(nVerts);
  std::vector<double> z(nVerts);
  std::vector<double> field(nVerts);

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

        x[idx] = pt[0];
        y[idx] = pt[1];
        z[idx] = pt[2];
        field[idx] = fieldFunctor(pt);
      }
    }
  }

  appendStructuredDomainData(data,
                             spec,
                             std::move(x),
                             std::move(y),
                             std::move(z),
                             std::move(field),
                             fieldName);
}

void buildBlueprintMesh(BlueprintMeshData& data)
{
  populateCoordinates(data);
  populateNodalVolumeFractions(data);
  populateBlueprintNode(data);
}

void buildContourMesh(const conduit::Node& blueprintMesh,
                      const std::string& fieldName,
                      double contourValue,
                      UMesh& contourMesh,
                      const std::string& maskField = "")
{
  conduit::Node verifyInfo;
  EXPECT_TRUE(conduit::blueprint::mesh::verify(blueprintMesh, verifyInfo)) << verifyInfo.to_yaml();

  quest::MarchingCubes mc(axom::runtime_policy::Policy::seq,
                          axom::getDefaultAllocatorID(),
                          quest::MarchingCubesDataParallelism::byPolicy);
  mc.setMesh(blueprintMesh, "mesh", maskField);
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

std::vector<Point3D> getContourNodes(const UMesh& contourMesh)
{
  std::vector<Point3D> nodes(contourMesh.getNumberOfNodes());
  double coords[3] {};

  for(axom::IndexType nodeId = 0; nodeId < contourMesh.getNumberOfNodes(); ++nodeId)
  {
    contourMesh.getNode(nodeId, coords);
    nodes[nodeId] = Point3D {coords[0], coords[1], coords[2]};
  }

  return nodes;
}

void sortPoints(std::vector<Point3D>& points)
{
  std::sort(points.begin(), points.end(), [](const Point3D& lhs, const Point3D& rhs) {
    if(lhs[0] != rhs[0])
    {
      return lhs[0] < rhs[0];
    }
    if(lhs[1] != rhs[1])
    {
      return lhs[1] < rhs[1];
    }
    return lhs[2] < rhs[2];
  });
}

const axom::IndexType* getParentCellIds(const UMesh& contourMesh)
{
  axom::IndexType numComponents = -1;
  const axom::IndexType* ptr =
    contourMesh.getFieldPtr<axom::IndexType>("parent_cell", mint::CELL_CENTERED, numComponents);
  EXPECT_EQ(numComponents, 1);
  return ptr;
}

const DomainIdType* getDomainIds(const UMesh& contourMesh)
{
  return contourMesh.getFieldPtr<DomainIdType>("domain_id", mint::CELL_CENTERED);
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

MultiDomainBlueprintData makeUniformReferenceBlueprint()
{
  MultiDomainBlueprintData data;
  const StructuredDomainSpec spec {2, 2, 2, 0};
  const axom::IndexType nVerts = vertexCount(spec.ni, spec.nj, spec.nk);

  std::vector<double> x(nVerts);
  std::vector<double> y(nVerts);
  std::vector<double> z(nVerts);
  std::vector<double> field(nVerts, -1.0);

  for(axom::IndexType k = 0; k <= spec.nk; ++k)
  {
    for(axom::IndexType j = 0; j <= spec.nj; ++j)
    {
      for(axom::IndexType i = 0; i <= spec.ni; ++i)
      {
        const axom::IndexType idx = vertexIndex(i, j, k, spec.ni, spec.nj);
        x[idx] = static_cast<double>(i);
        y[idx] = static_cast<double>(j);
        z[idx] = static_cast<double>(k);
      }
    }
  }

  field[vertexIndex(0, 0, 0, spec.ni, spec.nj)] = 1.0;
  field[vertexIndex(2, 2, 2, spec.ni, spec.nj)] = 1.0;

  appendStructuredDomainData(data,
                             spec,
                             std::move(x),
                             std::move(y),
                             std::move(z),
                             std::move(field),
                             "field");
  return data;
}

MultiDomainBlueprintData makeNonuniformSingleTriangleBlueprint()
{
  MultiDomainBlueprintData data;
  const StructuredDomainSpec spec {2, 2, 2, 0};
  const std::array<double, 3> xCoords {{0.0, 2.0, 5.0}};
  const std::array<double, 3> yCoords {{0.0, 1.0, 4.0}};
  const std::array<double, 3> zCoords {{0.0, 3.0, 9.0}};
  const axom::IndexType nVerts = vertexCount(spec.ni, spec.nj, spec.nk);

  std::vector<double> x(nVerts);
  std::vector<double> y(nVerts);
  std::vector<double> z(nVerts);
  std::vector<double> field(nVerts, -2.0);

  for(axom::IndexType k = 0; k <= spec.nk; ++k)
  {
    for(axom::IndexType j = 0; j <= spec.nj; ++j)
    {
      for(axom::IndexType i = 0; i <= spec.ni; ++i)
      {
        const axom::IndexType idx = vertexIndex(i, j, k, spec.ni, spec.nj);
        x[idx] = xCoords[i];
        y[idx] = yCoords[j];
        z[idx] = zCoords[k];
      }
    }
  }

  field[vertexIndex(0, 0, 0, spec.ni, spec.nj)] = 2.0;
  field[vertexIndex(0, 0, 1, spec.ni, spec.nj)] = -4.0;

  appendStructuredDomainData(data,
                             spec,
                             std::move(x),
                             std::move(y),
                             std::move(z),
                             std::move(field),
                             "field");
  return data;
}

MultiDomainBlueprintData makeMultidomainSeamBlueprint()
{
  MultiDomainBlueprintData data;
  const StructuredDomainSpec leftSpec {2, 2, 2, 7};
  const StructuredDomainSpec rightSpec {2, 2, 2, 11};

  const auto leftCoords = [](double u, double v, double w) { return Point3D {u, v, w}; };
  const auto rightCoords = [](double u, double v, double w) { return Point3D {1.0 + u, v, w}; };
  const auto fieldFunctor = [](const Point3D& pt) { return pt[1] + pt[2] - 0.75; };

  appendStructuredDomain(data, leftSpec, leftCoords, fieldFunctor, "field");
  appendStructuredDomain(data, rightSpec, rightCoords, fieldFunctor, "field");
  return data;
}

const std::array<Point3D, 8>& unitCubeVertices()
{
  static const std::array<Point3D, 8> vertices {Point3D {1.0, 0.0, 0.0},
                                                Point3D {1.0, 1.0, 0.0},
                                                Point3D {0.0, 1.0, 0.0},
                                                Point3D {0.0, 0.0, 0.0},
                                                Point3D {1.0, 0.0, 1.0},
                                                Point3D {1.0, 1.0, 1.0},
                                                Point3D {0.0, 1.0, 1.0},
                                                Point3D {0.0, 0.0, 1.0}};
  return vertices;
}

const std::array<EdgeVertexIds, 12>& unitCubeEdgeVertices()
{
  static const std::array<EdgeVertexIds, 12> edges {{
    {{0, 1}},
    {{1, 2}},
    {{3, 2}},
    {{0, 3}},
    {{4, 5}},
    {{5, 6}},
    {{7, 6}},
    {{4, 7}},
    {{0, 4}},
    {{1, 5}},
    {{2, 6}},
    {{3, 7}},
  }};
  return edges;
}

HexVertexIds singleCubeVertexIds()
{
  return {{vertexIndex(1, 0, 0, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj),
           vertexIndex(1, 1, 0, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj),
           vertexIndex(0, 1, 0, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj),
           vertexIndex(0, 0, 0, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj),
           vertexIndex(1, 0, 1, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj),
           vertexIndex(1, 1, 1, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj),
           vertexIndex(0, 1, 1, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj),
           vertexIndex(0, 0, 1, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj)}};
}

Point3D edgeMidpoint(int edgeId)
{
  const auto& vertices = unitCubeVertices();
  const auto& edge = unitCubeEdgeVertices()[edgeId];
  Point3D midpoint;

  for(int dim = 0; dim < 3; ++dim)
  {
    midpoint[dim] = 0.5 * (vertices[edge[0]][dim] + vertices[edge[1]][dim]);
  }

  return midpoint;
}

bool pointsEqualWithinTolerance(const Point3D& lhs, const Point3D& rhs, double tol = POINT_TOL)
{
  return std::abs(lhs[0] - rhs[0]) <= tol && std::abs(lhs[1] - rhs[1]) <= tol &&
    std::abs(lhs[2] - rhs[2]) <= tol;
}

std::vector<Point3D> uniqueSortedPoints(std::vector<Point3D> points, double tol = POINT_TOL)
{
  sortPoints(points);
  points.erase(std::unique(points.begin(),
                           points.end(),
                           [=](const Point3D& lhs, const Point3D& rhs) {
                             return pointsEqualWithinTolerance(lhs, rhs, tol);
                           }),
               points.end());
  return points;
}

std::vector<Point3D> expectedCaseContourNodes(int caseId)
{
  std::set<int> edges;
  for(int idx = 0; idx < 16 && cases3D[caseId][idx] >= 0; ++idx)
  {
    edges.insert(cases3D[caseId][idx]);
  }

  std::vector<Point3D> expectedNodes;
  expectedNodes.reserve(edges.size());
  for(const int edgeId : edges)
  {
    expectedNodes.push_back(edgeMidpoint(edgeId));
  }

  return expectedNodes;
}

MaskedSingleCubeBlueprintData makeSingleCubeCaseBlueprint(int caseId)
{
  MaskedSingleCubeBlueprintData data;
  const axom::IndexType nVerts =
    vertexCount(SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj, SINGLE_CUBE_SPEC.nk);
  const axom::IndexType nCells =
    cellCount(SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj, SINGLE_CUBE_SPEC.nk);
  std::vector<double> x(nVerts);
  std::vector<double> y(nVerts);
  std::vector<double> z(nVerts);
  std::vector<double> field(nVerts, -1.0);
  std::vector<int> mask(nCells, 0);

  for(axom::IndexType k = 0; k <= SINGLE_CUBE_SPEC.nk; ++k)
  {
    for(axom::IndexType j = 0; j <= SINGLE_CUBE_SPEC.nj; ++j)
    {
      for(axom::IndexType i = 0; i <= SINGLE_CUBE_SPEC.ni; ++i)
      {
        const axom::IndexType idx = vertexIndex(i, j, k, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj);
        x[idx] = static_cast<double>(i);
        y[idx] = static_cast<double>(j);
        z[idx] = static_cast<double>(k);
      }
    }
  }

  const auto nodeIds = singleCubeVertexIds();
  for(int vertex = 0; vertex < 8; ++vertex)
  {
    field[nodeIds[vertex]] = (caseId & (1 << vertex)) ? 1.0 : -1.0;
  }

  mask[cellIndex(0, 0, 0, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj)] = 1;

  data.x = std::move(x);
  data.y = std::move(y);
  data.z = std::move(z);
  data.field = std::move(field);
  data.mask = std::move(mask);

  conduit::Node& domain = data.mesh.append();
  domain["coordsets/coords/type"] = "explicit";
  domain["coordsets/coords/values/x"].set_external(conduit::DataType::float64(nVerts), data.x.data());
  domain["coordsets/coords/values/y"].set_external(conduit::DataType::float64(nVerts), data.y.data());
  domain["coordsets/coords/values/z"].set_external(conduit::DataType::float64(nVerts), data.z.data());

  domain["topologies/mesh/type"] = "structured";
  domain["topologies/mesh/coordset"] = "coords";
  domain["topologies/mesh/elements/dims/i"] = SINGLE_CUBE_SPEC.ni;
  domain["topologies/mesh/elements/dims/j"] = SINGLE_CUBE_SPEC.nj;
  domain["topologies/mesh/elements/dims/k"] = SINGLE_CUBE_SPEC.nk;
  domain["state/domain_id"] = 0;

  domain["fields/field/association"] = "vertex";
  domain["fields/field/topology"] = "mesh";
  domain["fields/field/volume_dependent"] = "false";
  domain["fields/field/values"].set_external(conduit::DataType::float64(nVerts), data.field.data());

  domain["fields/mask/association"] = "element";
  domain["fields/mask/topology"] = "mesh";
  domain["fields/mask/values"].set_external(conduit::DataType::int32(nCells), data.mask.data());
  return data;
}

std::vector<Point3D> getUniqueContourNodes(const UMesh& contourMesh)
{
  return uniqueSortedPoints(getContourNodes(contourMesh));
}

void expectPointVectorsEqual(const std::vector<Point3D>& actual,
                             const std::vector<Point3D>& expected,
                             double tol = POINT_TOL)
{
  ASSERT_EQ(actual.size(), expected.size());
  for(std::size_t idx = 0; idx < actual.size(); ++idx)
  {
    EXPECT_NEAR(actual[idx][0], expected[idx][0], tol);
    EXPECT_NEAR(actual[idx][1], expected[idx][1], tol);
    EXPECT_NEAR(actual[idx][2], expected[idx][2], tol);
  }
}

std::vector<Point3D> runSingleCubeCaseAndGetGeometry(int caseId, axom::IndexType& triangleCount)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeSingleCubeCaseBlueprint(caseId);
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh, "mask");
  triangleCount = contourMesh.getNumberOfCells();
  return getUniqueContourNodes(contourMesh);
}

void verifySingleCubeCanonicalCase(int caseId)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeSingleCubeCaseBlueprint(caseId);
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh, "mask");

  ASSERT_EQ(contourMesh.getNumberOfCells(), num_triangles[caseId]);

  std::vector<Point3D> expectedUniqueNodes = expectedCaseContourNodes(caseId);
  sortPoints(expectedUniqueNodes);

  for(const Point3D& node : getContourNodes(contourMesh))
  {
    const bool isExpectedNode = std::any_of(
      expectedUniqueNodes.begin(),
      expectedUniqueNodes.end(),
      [&](const Point3D& expectedNode) { return pointsEqualWithinTolerance(node, expectedNode); });
    EXPECT_TRUE(isExpectedNode) << "Unexpected contour node " << node << " for case " << caseId;
  }

  expectPointVectorsEqual(getUniqueContourNodes(contourMesh), expectedUniqueNodes);

  for(axom::IndexType cellId = 0; cellId < contourMesh.getNumberOfCells(); ++cellId)
  {
    EXPECT_EQ(contourMesh.getCellType(cellId), mint::TRIANGLE);
    EXPECT_EQ(contourMesh.getNumberOfCellNodes(cellId), 3);
  }

  const axom::IndexType* parentCellIds = getParentCellIds(contourMesh);
  const DomainIdType* domainIds = getDomainIds(contourMesh);

  if(contourMesh.getNumberOfCells() > 0)
  {
    ASSERT_NE(parentCellIds, nullptr);
    ASSERT_NE(domainIds, nullptr);
  }

  for(axom::IndexType cellId = 0; cellId < contourMesh.getNumberOfCells(); ++cellId)
  {
    EXPECT_EQ(parentCellIds[cellId], 0);
    EXPECT_EQ(domainIds[cellId], 0);
  }
}

}  // namespace

//------------------------------------------------------------------------------
// Tests
//------------------------------------------------------------------------------

TEST(quest_marching_cubes, banded_nodal_volume_fraction_surface)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  buildBandedVolumeFractionContourMesh(contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  ASSERT_GT(contourMesh.getNumberOfNodes(), 0);

  expectRaysHit(contourMesh, controlRays());

  // These should hit the contoured sphere, but currently miss because the
  // generated surface has a hole near the z-directed path. Keep these as
  // ordinary expectations so the known bug remains tracked.
  expectRaysHit(contourMesh, challengingZRays());
}

TEST(quest_marching_cubes, uniform_structured_reference_surface)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeUniformReferenceBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_EQ(contourMesh.getNumberOfCells(), 2);
  ASSERT_EQ(contourMesh.getNumberOfNodes(), 6);

  for(axom::IndexType cellId = 0; cellId < contourMesh.getNumberOfCells(); ++cellId)
  {
    EXPECT_EQ(contourMesh.getCellType(cellId), mint::TRIANGLE);
    EXPECT_EQ(contourMesh.getNumberOfCellNodes(cellId), 3);
  }

  const axom::IndexType* parentCellIds = getParentCellIds(contourMesh);
  std::set<axom::IndexType> parentCellSet(parentCellIds,
                                          parentCellIds + contourMesh.getNumberOfCells());
  const std::set<axom::IndexType> expectedParents {cellIndex(0, 0, 0, 2, 2),
                                                   cellIndex(1, 1, 1, 2, 2)};
  EXPECT_EQ(parentCellSet, expectedParents);
}

TEST(quest_marching_cubes, uniform_sphere_surface)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeUniformSphereBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  expectRaysHit(contourMesh, uniformSphereRays());
}

TEST(quest_marching_cubes, warped_sphere_surface)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeWarpedSphereBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  expectRaysHit(contourMesh, warpedSphereRays());
}

TEST(quest_marching_cubes, planar_surface)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makePlaneCutBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  expectRaysHit(contourMesh, planeCutRays());
}

TEST(quest_marching_cubes, nonuniform_single_triangle_reference)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeNonuniformSingleTriangleBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_EQ(contourMesh.getNumberOfCells(), 1);
  ASSERT_EQ(contourMesh.getNumberOfNodes(), 3);

  const axom::IndexType* parentCellIds = getParentCellIds(contourMesh);
  const DomainIdType* domainIds = getDomainIds(contourMesh);

  ASSERT_NE(parentCellIds, nullptr);
  ASSERT_NE(domainIds, nullptr);
  EXPECT_EQ(parentCellIds[0], 0);
  EXPECT_EQ(domainIds[0], 0);

  std::vector<Point3D> actualNodes = getContourNodes(contourMesh);
  std::vector<Point3D> expectedNodes {Point3D {0.0, 0.0, 1.0},
                                      Point3D {0.0, 0.5, 0.0},
                                      Point3D {1.0, 0.0, 0.0}};
  sortPoints(actualNodes);
  sortPoints(expectedNodes);

  constexpr double tol = 1.0e-12;
  for(int idx = 0; idx < 3; ++idx)
  {
    EXPECT_NEAR(actualNodes[idx][0], expectedNodes[idx][0], tol);
    EXPECT_NEAR(actualNodes[idx][1], expectedNodes[idx][1], tol);
    EXPECT_NEAR(actualNodes[idx][2], expectedNodes[idx][2], tol);
  }
}

TEST(quest_marching_cubes, multidomain_seam_sphere_surface)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeSeamSphereBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  expectRaysHit(contourMesh, seamSphereRays());
}

TEST(quest_marching_cubes, multidomain_seam_surface)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  auto meshData = makeMultidomainSeamBlueprint();
  buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);

  const axom::IndexType* parentCellIds = getParentCellIds(contourMesh);
  const DomainIdType* domainIds = getDomainIds(contourMesh);

  ASSERT_NE(parentCellIds, nullptr);
  ASSERT_NE(domainIds, nullptr);

  std::set<DomainIdType> observedDomainIds(domainIds, domainIds + contourMesh.getNumberOfCells());
  const std::set<DomainIdType> expectedDomainIds {7, 11};
  EXPECT_EQ(observedDomainIds, expectedDomainIds);

  const std::map<DomainIdType, axom::IndexType> validCellCounts {{7, cellCount(2, 2, 2)},
                                                                 {11, cellCount(2, 2, 2)}};

  for(axom::IndexType cellId = 0; cellId < contourMesh.getNumberOfCells(); ++cellId)
  {
    const auto iter = validCellCounts.find(domainIds[cellId]);
    ASSERT_NE(iter, validCellCounts.end());
    EXPECT_GE(parentCellIds[cellId], 0);
    EXPECT_LT(parentCellIds[cellId], iter->second);
  }
}

TEST(quest_marching_cubes, single_cube_all_256_cases)
{
  for(int caseId = 0; caseId < 256; ++caseId)
  {
    SCOPED_TRACE(axom::fmt::format("caseId={}", caseId));
    verifySingleCubeCanonicalCase(caseId);
  }
}

TEST(quest_marching_cubes, single_cube_complement_symmetry)
{
  for(int caseId = 0; caseId < 256; ++caseId)
  {
    const int complementCaseId = 255 - caseId;
    SCOPED_TRACE(axom::fmt::format("caseId={} complement={}", caseId, complementCaseId));

    axom::IndexType triangleCount = 0;
    axom::IndexType complementTriangleCount = 0;
    const std::vector<Point3D> geometry = runSingleCubeCaseAndGetGeometry(caseId, triangleCount);
    const std::vector<Point3D> complementGeometry =
      runSingleCubeCaseAndGetGeometry(complementCaseId, complementTriangleCount);

    EXPECT_EQ(triangleCount, num_triangles[caseId]);
    EXPECT_EQ(complementTriangleCount, num_triangles[complementCaseId]);

    // Complementary ambiguous cases can triangulate the same contour polygon
    // with different triangle counts, so only compare the geometry set here.
    expectPointVectorsEqual(geometry, complementGeometry);
  }
}

#endif
