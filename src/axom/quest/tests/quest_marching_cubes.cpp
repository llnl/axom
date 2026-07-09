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
  #include <utility>
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

struct StructuredCaseBlueprintData
{
  conduit::Node mesh;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> field;
  std::vector<int> mask;
  double isoValue {0.0};
};

struct UniquePointData
{
  std::vector<Point3D> points;
  std::vector<int> nodeToUnique;
};

struct ScalarRange
{
  double min {0.0};
  double max {0.0};
};

struct MeshInspection
{
  axom::IndexType rawNodeCount {0};
  axom::IndexType triangleCount {0};
  int uniquePointCount {0};
  int collocatedNodeCount {0};
  int duplicateUniquePointCount {0};
  int isolatedPointCount {0};
  int degenerateTriangleCount {0};
  int duplicateTriangleCount {0};
  int connectedComponentCount {0};
  int boundaryEdgeCount {0};

  bool hasOpenBoundary() const { return boundaryEdgeCount > 0; }
};

using HexVertexIds = std::array<axom::IndexType, 8>;
using EdgeVertexIds = std::array<int, 2>;
using CornerValues = std::array<double, 8>;
using TriangleNodeIds = std::array<int, 3>;
using UndirectedEdge = std::array<int, 2>;

constexpr StructuredDomainSpec SINGLE_CUBE_SPEC {2, 2, 2, 0};
constexpr StructuredDomainSpec EMBEDDED_CUBE_SPEC {3, 3, 3, 0};
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

MultiDomainBlueprintData makeUniformSphereBlueprint(
  const StructuredDomainSpec& spec = StructuredDomainSpec {32, 24, 24, 0})
{
  MultiDomainBlueprintData data;

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

const std::array<Point3D, 8>& cgalUnitCubeVertices()
{
  static const std::array<Point3D, 8> vertices {Point3D {0.0, 0.0, 0.0},
                                                Point3D {0.0, 0.0, 1.0},
                                                Point3D {0.0, 1.0, 0.0},
                                                Point3D {0.0, 1.0, 1.0},
                                                Point3D {1.0, 0.0, 0.0},
                                                Point3D {1.0, 0.0, 1.0},
                                                Point3D {1.0, 1.0, 0.0},
                                                Point3D {1.0, 1.0, 1.0}};
  return vertices;
}

HexVertexIds structuredCellVertexIds(axom::IndexType i,
                                     axom::IndexType j,
                                     axom::IndexType k,
                                     const StructuredDomainSpec& spec)
{
  return {{vertexIndex(i + 1, j, k, spec.ni, spec.nj),
           vertexIndex(i + 1, j + 1, k, spec.ni, spec.nj),
           vertexIndex(i, j + 1, k, spec.ni, spec.nj),
           vertexIndex(i, j, k, spec.ni, spec.nj),
           vertexIndex(i + 1, j, k + 1, spec.ni, spec.nj),
           vertexIndex(i + 1, j + 1, k + 1, spec.ni, spec.nj),
           vertexIndex(i, j + 1, k + 1, spec.ni, spec.nj),
           vertexIndex(i, j, k + 1, spec.ni, spec.nj)}};
}

CornerValues cgalToAxomCornerValues(const CornerValues& cgalValues)
{
  return {{cgalValues[4],
           cgalValues[6],
           cgalValues[2],
           cgalValues[0],
           cgalValues[5],
           cgalValues[7],
           cgalValues[3],
           cgalValues[1]}};
}

HexVertexIds singleCubeVertexIds() { return structuredCellVertexIds(0, 0, 0, SINGLE_CUBE_SPEC); }

HexVertexIds embeddedCubeVertexIds()
{
  return structuredCellVertexIds(1, 1, 1, EMBEDDED_CUBE_SPEC);
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

ScalarRange computeScalarRange(const std::vector<double>& field)
{
  EXPECT_FALSE(field.empty());
  if(field.empty())
  {
    return {};
  }

  const auto minmax = std::minmax_element(field.begin(), field.end());
  return {*minmax.first, *minmax.second};
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

StructuredCaseBlueprintData makeSingleCellStructuredCaseBlueprint(const CornerValues& cgalCornerValues,
                                                                  double isoValue)
{
  StructuredCaseBlueprintData data;
  data.isoValue = isoValue;

  const axom::IndexType nVerts =
    vertexCount(SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj, SINGLE_CUBE_SPEC.nk);
  const axom::IndexType nCells =
    cellCount(SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj, SINGLE_CUBE_SPEC.nk);
  data.x.resize(nVerts);
  data.y.resize(nVerts);
  data.z.resize(nVerts);
  data.field.resize(nVerts, isoValue + 1.0);
  data.mask.resize(nCells, 0);

  for(axom::IndexType k = 0; k <= SINGLE_CUBE_SPEC.nk; ++k)
  {
    for(axom::IndexType j = 0; j <= SINGLE_CUBE_SPEC.nj; ++j)
    {
      for(axom::IndexType i = 0; i <= SINGLE_CUBE_SPEC.ni; ++i)
      {
        const axom::IndexType idx = vertexIndex(i, j, k, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj);
        data.x[idx] = static_cast<double>(i);
        data.y[idx] = static_cast<double>(j);
        data.z[idx] = static_cast<double>(k);
      }
    }
  }

  const CornerValues axomCornerValues = cgalToAxomCornerValues(cgalCornerValues);
  const auto nodeIds = singleCubeVertexIds();
  for(int vertex = 0; vertex < 8; ++vertex)
  {
    data.field[nodeIds[vertex]] = axomCornerValues[vertex];
  }

  data.mask[cellIndex(0, 0, 0, SINGLE_CUBE_SPEC.ni, SINGLE_CUBE_SPEC.nj)] = 1;

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

StructuredCaseBlueprintData makeEmbeddedStructuredCaseBlueprint(const CornerValues& cgalCornerValues,
                                                                double isoValue)
{
  StructuredCaseBlueprintData data;
  data.isoValue = isoValue;

  const axom::IndexType nVerts =
    vertexCount(EMBEDDED_CUBE_SPEC.ni, EMBEDDED_CUBE_SPEC.nj, EMBEDDED_CUBE_SPEC.nk);
  data.x.resize(nVerts);
  data.y.resize(nVerts);
  data.z.resize(nVerts);
  data.field.resize(nVerts);

  double maxCornerValue = isoValue;
  double minCornerValue = isoValue;
  for(double value : cgalCornerValues)
  {
    maxCornerValue = std::max(maxCornerValue, value);
    minCornerValue = std::min(minCornerValue, value);
  }
  const double backgroundValue =
    maxCornerValue + std::max(1.0, maxCornerValue - minCornerValue + 1.0e-3);
  std::fill(data.field.begin(), data.field.end(), backgroundValue);

  for(axom::IndexType k = 0; k <= EMBEDDED_CUBE_SPEC.nk; ++k)
  {
    for(axom::IndexType j = 0; j <= EMBEDDED_CUBE_SPEC.nj; ++j)
    {
      for(axom::IndexType i = 0; i <= EMBEDDED_CUBE_SPEC.ni; ++i)
      {
        const axom::IndexType idx =
          vertexIndex(i, j, k, EMBEDDED_CUBE_SPEC.ni, EMBEDDED_CUBE_SPEC.nj);
        data.x[idx] = static_cast<double>(i);
        data.y[idx] = static_cast<double>(j);
        data.z[idx] = static_cast<double>(k);
      }
    }
  }

  const CornerValues axomCornerValues = cgalToAxomCornerValues(cgalCornerValues);
  const auto nodeIds = embeddedCubeVertexIds();
  for(int vertex = 0; vertex < 8; ++vertex)
  {
    data.field[nodeIds[vertex]] = axomCornerValues[vertex];
  }

  conduit::Node& domain = data.mesh.append();
  domain["coordsets/coords/type"] = "explicit";
  domain["coordsets/coords/values/x"].set_external(conduit::DataType::float64(nVerts), data.x.data());
  domain["coordsets/coords/values/y"].set_external(conduit::DataType::float64(nVerts), data.y.data());
  domain["coordsets/coords/values/z"].set_external(conduit::DataType::float64(nVerts), data.z.data());

  domain["topologies/mesh/type"] = "structured";
  domain["topologies/mesh/coordset"] = "coords";
  domain["topologies/mesh/elements/dims/i"] = EMBEDDED_CUBE_SPEC.ni;
  domain["topologies/mesh/elements/dims/j"] = EMBEDDED_CUBE_SPEC.nj;
  domain["topologies/mesh/elements/dims/k"] = EMBEDDED_CUBE_SPEC.nk;
  domain["state/domain_id"] = 0;

  domain["fields/field/association"] = "vertex";
  domain["fields/field/topology"] = "mesh";
  domain["fields/field/volume_dependent"] = "false";
  domain["fields/field/values"].set_external(conduit::DataType::float64(nVerts), data.field.data());
  return data;
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

UniquePointData extractUniquePoints(const UMesh& contourMesh, double tol = POINT_TOL)
{
  struct PointEntry
  {
    Point3D point;
    axom::IndexType nodeId;
  };

  std::vector<PointEntry> entries;
  entries.reserve(contourMesh.getNumberOfNodes());

  double coords[3] {};
  for(axom::IndexType nodeId = 0; nodeId < contourMesh.getNumberOfNodes(); ++nodeId)
  {
    contourMesh.getNode(nodeId, coords);
    entries.push_back({Point3D {coords[0], coords[1], coords[2]}, nodeId});
  }

  std::sort(entries.begin(), entries.end(), [](const PointEntry& lhs, const PointEntry& rhs) {
    if(lhs.point[0] != rhs.point[0])
    {
      return lhs.point[0] < rhs.point[0];
    }
    if(lhs.point[1] != rhs.point[1])
    {
      return lhs.point[1] < rhs.point[1];
    }
    return lhs.point[2] < rhs.point[2];
  });

  UniquePointData uniqueData;
  uniqueData.nodeToUnique.resize(contourMesh.getNumberOfNodes(), -1);

  for(const auto& entry : entries)
  {
    if(uniqueData.points.empty() ||
       !pointsEqualWithinTolerance(entry.point, uniqueData.points.back(), tol))
    {
      uniqueData.points.push_back(entry.point);
    }

    uniqueData.nodeToUnique[entry.nodeId] = static_cast<int>(uniqueData.points.size()) - 1;
  }

  return uniqueData;
}

double twiceTriangleArea(const Point3D& a, const Point3D& b, const Point3D& c)
{
  const double abx = b[0] - a[0];
  const double aby = b[1] - a[1];
  const double abz = b[2] - a[2];

  const double acx = c[0] - a[0];
  const double acy = c[1] - a[1];
  const double acz = c[2] - a[2];

  const double crossX = aby * acz - abz * acy;
  const double crossY = abz * acx - abx * acz;
  const double crossZ = abx * acy - aby * acx;
  return std::sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
}

MeshInspection inspectContourMesh(const UMesh& contourMesh, double tol = POINT_TOL)
{
  MeshInspection inspection;
  inspection.rawNodeCount = contourMesh.getNumberOfNodes();
  inspection.triangleCount = contourMesh.getNumberOfCells();

  const UniquePointData uniqueData = extractUniquePoints(contourMesh, tol);
  inspection.uniquePointCount = static_cast<int>(uniqueData.points.size());
  inspection.collocatedNodeCount =
    static_cast<int>(inspection.rawNodeCount) - inspection.uniquePointCount;
  inspection.duplicateUniquePointCount = inspection.collocatedNodeCount;

  std::vector<int> pointUseCount(uniqueData.points.size(), 0);
  std::map<TriangleNodeIds, int> triangleCounts;
  std::map<UndirectedEdge, std::vector<int>> edgeToTriangles;
  std::vector<std::vector<int>> adjacency;

  double coords[3] {};
  std::vector<Point3D> rawNodes(contourMesh.getNumberOfNodes());
  for(axom::IndexType nodeId = 0; nodeId < contourMesh.getNumberOfNodes(); ++nodeId)
  {
    contourMesh.getNode(nodeId, coords);
    rawNodes[nodeId] = Point3D {coords[0], coords[1], coords[2]};
  }

  int validTriangleCount = 0;
  for(axom::IndexType cellId = 0; cellId < contourMesh.getNumberOfCells(); ++cellId)
  {
    const axom::IndexType* nodeIds = contourMesh.getCellNodeIDs(cellId);
    TriangleNodeIds triangleIds {{uniqueData.nodeToUnique[nodeIds[0]],
                                  uniqueData.nodeToUnique[nodeIds[1]],
                                  uniqueData.nodeToUnique[nodeIds[2]]}};

    for(int nodeId : triangleIds)
    {
      ++pointUseCount[nodeId];
    }

    const Point3D& a = rawNodes[nodeIds[0]];
    const Point3D& b = rawNodes[nodeIds[1]];
    const Point3D& c = rawNodes[nodeIds[2]];
    const bool repeatedVertex = triangleIds[0] == triangleIds[1] ||
      triangleIds[1] == triangleIds[2] || triangleIds[0] == triangleIds[2];
    const bool degenerateByArea = twiceTriangleArea(a, b, c) <= tol;

    TriangleNodeIds canonicalTriangle = triangleIds;
    std::sort(canonicalTriangle.begin(), canonicalTriangle.end());
    if(++triangleCounts[canonicalTriangle] > 1)
    {
      ++inspection.duplicateTriangleCount;
    }

    if(repeatedVertex || degenerateByArea)
    {
      ++inspection.degenerateTriangleCount;
      continue;
    }

    adjacency.emplace_back();
    const int triangleIndex = validTriangleCount++;

    for(int edgeIdx = 0; edgeIdx < 3; ++edgeIdx)
    {
      UndirectedEdge edge {{triangleIds[edgeIdx], triangleIds[(edgeIdx + 1) % 3]}};
      if(edge[0] > edge[1])
      {
        std::swap(edge[0], edge[1]);
      }

      auto& incidentTriangles = edgeToTriangles[edge];
      for(int adjacentTriangle : incidentTriangles)
      {
        adjacency[triangleIndex].push_back(adjacentTriangle);
        adjacency[adjacentTriangle].push_back(triangleIndex);
      }
      incidentTriangles.push_back(triangleIndex);
    }
  }

  inspection.boundaryEdgeCount = 0;
  for(const auto& edgeEntry : edgeToTriangles)
  {
    if(edgeEntry.second.size() == 1)
    {
      ++inspection.boundaryEdgeCount;
    }
  }

  inspection.connectedComponentCount = 0;
  std::vector<bool> visited(adjacency.size(), false);
  for(std::size_t triangleIdx = 0; triangleIdx < adjacency.size(); ++triangleIdx)
  {
    if(visited[triangleIdx])
    {
      continue;
    }

    ++inspection.connectedComponentCount;
    std::vector<int> stack(1, static_cast<int>(triangleIdx));
    visited[triangleIdx] = true;

    while(!stack.empty())
    {
      const int current = stack.back();
      stack.pop_back();

      for(int neighbor : adjacency[current])
      {
        if(!visited[neighbor])
        {
          visited[neighbor] = true;
          stack.push_back(neighbor);
        }
      }
    }
  }

  inspection.isolatedPointCount = 0;
  for(int useCount : pointUseCount)
  {
    if(useCount == 0)
    {
      ++inspection.isolatedPointCount;
    }
  }

  return inspection;
}

bool containsPoint(const std::vector<Point3D>& points, const Point3D& expected, double tol = POINT_TOL)
{
  return std::any_of(points.begin(), points.end(), [&](const Point3D& point) {
    return pointsEqualWithinTolerance(point, expected, tol);
  });
}

void expectContainsPoints(const std::vector<Point3D>& actualPoints,
                          const std::vector<Point3D>& expectedPoints,
                          double tol = POINT_TOL)
{
  for(const Point3D& expectedPoint : expectedPoints)
  {
    EXPECT_TRUE(containsPoint(actualPoints, expectedPoint, tol))
      << "Expected contour to contain point " << expectedPoint;
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

enum class AmbiguousCase
{
  CONTOUR_6_VTS,
  CONTOUR_7_VTS_CONFIG_A,
  MC_4_TUNNEL,
  MC_13_SIMPLE
};

enum class SingularCase
{
  SIMPLE,
  DOUBLE,
  TRIPLE,
  TUNNEL_OPEN,
  TUNNEL,
  TUNNEL_DOUBLE
};

enum class IsoEdgeCase
{
  SIMPLE,
  SEPARATED,
  SINGULAR,
  DOUBLE_SINGULAR_SIMPLE,
  DOUBLE_SINGULAR,
  DOUBLE_SINGULAR_SPLIT,
  WEDGE,
  WEDGE_SINGULAR,
  EDGE_TUNNEL
};

enum class PlaneCase
{
  SIMPLE,
  INTERSECTING,
  PLANE_TUBES,
  CROSS
};

double generatePredefinedInnerAmbiguity(AmbiguousCase topologyCase, CornerValues& values)
{
  switch(topologyCase)
  {
  case AmbiguousCase::CONTOUR_6_VTS:
    values = {{-0.960492426392903,
               0.793207329559554,
               0.916735525189067,
               -0.422761282626275,
               -0.934993247757551,
               -0.850129305868777,
               -0.0367116785741896,
               -0.656740699156587}};
    return 0.0387;
  case AmbiguousCase::CONTOUR_7_VTS_CONFIG_A:
    values = {{10.2967816247556,
               9.45145192686147,
               9.54753711271687,
               10.6482067822841,
               9.81494966341055,
               9.31168538578250,
               9.80950580411527,
               10.7451536262220}};
    return 9.8588;
  case AmbiguousCase::MC_4_TUNNEL:
    values = {{-7.70146936482581,
               -3.21868369245987,
               -5.44023748418735,
               15.6051950593180,
               12.7611835388515,
               -4.46952393442309,
               -11.7240576326183,
               -9.23038948829007}};
    return -1.7660;
  case AmbiguousCase::MC_13_SIMPLE:
    values = {{0.520482995461163,
               -0.839814387388296,
               -0.467491517013617,
               0.937814095887345,
               -0.825777099007084,
               0.506695544835103,
               0.345318915961394,
               -0.861107217966913}};
    return 0.0293;
  }

  values.fill(0.0);
  return 0.0;
}

double generatePredefinedSingular(SingularCase topologyCase, CornerValues& values)
{
  switch(topologyCase)
  {
  case SingularCase::SIMPLE:
    values = {{-1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0}};
    break;
  case SingularCase::DOUBLE:
    values = {{-1.0, -1.0, -0.5, 0.5, -2.0, 1.0, 0.5, -0.5}};
    break;
  case SingularCase::TRIPLE:
    values = {{-1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0}};
    break;
  case SingularCase::TUNNEL_OPEN:
    values = {{-0.5, 0.5, -0.4, 0.46, 0.5, -0.5, 0.71, -0.92}};
    break;
  case SingularCase::TUNNEL:
    values = {{-0.5, 0.5, -0.33, 0.46, 0.5, -0.5, 0.71, -0.92}};
    break;
  case SingularCase::TUNNEL_DOUBLE:
    values = {{-1.2, -0.6, 1.2, 0.9, 9.7, 0.6, -9.7, -0.9}};
    break;
  }

  return 0.0;
}

double generatePredefinedIsoEdge(IsoEdgeCase topologyCase, CornerValues& values)
{
  switch(topologyCase)
  {
  case IsoEdgeCase::SIMPLE:
    values = {{0.0, 0.0, -0.4, -0.3, 0.7, 0.5, -0.71, -0.92}};
    break;
  case IsoEdgeCase::SEPARATED:
    values = {{0.0, 0.0, 0.4, 0.3, 0.7, 0.5, -0.71, -0.92}};
    break;
  case IsoEdgeCase::SINGULAR:
    values = {{0.0, 0.0, 0.4, 0.3, -0.7, 0.5, -0.71, -0.92}};
    break;
  case IsoEdgeCase::DOUBLE_SINGULAR_SIMPLE:
    values = {{0.0, 0.0, -0.8, 0.3, -0.7, 0.5, -0.2, 0.1}};
    break;
  case IsoEdgeCase::DOUBLE_SINGULAR:
    values = {{0.0, 0.0, -0.2, 0.3, -0.7, 0.5, -0.2, -0.92}};
    break;
  case IsoEdgeCase::DOUBLE_SINGULAR_SPLIT:
    values = {{0.0, 0.0, 0.2, -0.3, -0.7, 0.5, -0.2, -0.92}};
    break;
  case IsoEdgeCase::WEDGE:
    values = {{0.0, 0.0, -0.4, -0.3, 0.0, 0.5, -0.2, -0.92}};
    break;
  case IsoEdgeCase::WEDGE_SINGULAR:
    values = {{0.0, 0.0, 0.2, 0.3, 0.0, 0.5, -0.2, -0.92}};
    break;
  case IsoEdgeCase::EDGE_TUNNEL:
    values = {{0.3, 0.4, 0.7, -0.4, -0.7, 0.0, 0.4, 0.0}};
    break;
  }

  return 0.0;
}

double generatePredefinedPlane(PlaneCase topologyCase, CornerValues& values)
{
  switch(topologyCase)
  {
  case PlaneCase::SIMPLE:
    values = {{0.0, 0.0, 0.0, 0.0, 0.7, 0.5, 0.71, 0.92}};
    break;
  case PlaneCase::INTERSECTING:
    values = {{-12.0, -91.0, 12.0, 91.0, 97.0, 9.0, -97.0, -9.0}};
    break;
  case PlaneCase::PLANE_TUBES:
    values = {{-0.5, 0.0, 0.4, 0.0, 0.5, 0.0, -0.5, 0.0}};
    break;
  case PlaneCase::CROSS:
    values = {{1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0}};
    break;
  }

  return 0.0;
}

void runStructuredCase(const StructuredCaseBlueprintData& meshData,
                       UMesh& contourMesh,
                       const std::string& maskField = {})
{
  buildContourMesh(meshData.mesh, "field", meshData.isoValue, contourMesh, maskField);
}

std::vector<Point3D> cgalIsoVertices(const CornerValues& values, double isoValue)
{
  std::vector<Point3D> isoVertices;
  const auto& vertices = cgalUnitCubeVertices();
  for(int idx = 0; idx < 8; ++idx)
  {
    if(axom::utilities::isNearlyEqual(values[idx], isoValue))
    {
      isoVertices.push_back(vertices[idx]);
    }
  }
  return isoVertices;
}

MultiDomainBlueprintData makeCgalGridSphereStyleBlueprint(axom::IndexType n)
{
  MultiDomainBlueprintData data;
  const StructuredDomainSpec spec {n - 1, 2 * (n - 1), 3 * (n - 1), 0};

  const auto coordFunctor = [=](double u, double v, double w) {
    return Point3D {-1.0 + 2.0 * u, -1.0 + 2.0 * v, -1.0 + 2.0 * w};
  };
  const auto fieldFunctor = [=](const Point3D& pt) {
    return std::sqrt(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]);
  };

  appendStructuredDomain(data, spec, coordFunctor, fieldFunctor, "field");
  return data;
}

}  // namespace

//------------------------------------------------------------------------------
// Tests
//------------------------------------------------------------------------------

TEST(quest_marching_cubes, banded_nodal_volume_fraction_surface)
{
  UMesh contourMesh(3, mint::TRIANGLE);
  buildBandedVolumeFractionContourMesh(contourMesh);
  const MeshInspection inspection = inspectContourMesh(contourMesh);

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);
  ASSERT_GT(contourMesh.getNumberOfNodes(), 0);
  EXPECT_EQ(inspection.degenerateTriangleCount, 0);
  EXPECT_EQ(inspection.isolatedPointCount, 0);
  EXPECT_EQ(inspection.boundaryEdgeCount, 0);
  EXPECT_EQ(inspection.connectedComponentCount, 1);

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

TEST(quest_marching_cubes, structured_sphere_outside_scalar_range_yields_empty_mesh)
{
  const auto expectEmptyOutsideRange = [](const MultiDomainBlueprintData& meshData, const char* name) {
    ASSERT_EQ(meshData.domains.size(), 1u);
    const ScalarRange range = computeScalarRange(meshData.domains[0].field);

    for(const double isoValue : {range.min - 1.0, range.max + 1.0})
    {
      SCOPED_TRACE(axom::fmt::format("{} iso={}", name, isoValue));

      UMesh contourMesh(3, mint::TRIANGLE);
      buildContourMesh(meshData.mesh, "field", isoValue, contourMesh);
      EXPECT_EQ(contourMesh.getNumberOfCells(), 0);
      EXPECT_EQ(contourMesh.getNumberOfNodes(), 0);
    }
  };

  const auto uniformSphere = makeUniformSphereBlueprint();
  expectEmptyOutsideRange(uniformSphere, "uniform_sphere");

  const auto warpedSphere = makeWarpedSphereBlueprint();
  expectEmptyOutsideRange(warpedSphere, "warped_sphere");
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

TEST(quest_marching_cubes, cgal_grid_sphere_style_regression)
{
  const std::array<axom::IndexType, 3> resolutions {{3, 6, 10}};
  const std::array<RayCase, 3> rays {{
    {"x_axis", Point3D {1.2, 0.0, 0.0}, Point3D {-1.2, 0.0, 0.0}},
    {"y_axis", Point3D {0.0, 1.2, 0.0}, Point3D {0.0, -1.2, 0.0}},
    {"z_axis", Point3D {0.0, 0.0, 1.2}, Point3D {0.0, 0.0, -1.2}},
  }};

  for(axom::IndexType n : resolutions)
  {
    SCOPED_TRACE(axom::fmt::format("grid_n={}", n));

    UMesh contourMesh(3, mint::TRIANGLE);
    auto meshData = makeCgalGridSphereStyleBlueprint(n);
    buildContourMesh(meshData.mesh, "field", 0.777, contourMesh);

    const MeshInspection inspection = inspectContourMesh(contourMesh);
    EXPECT_GT(contourMesh.getNumberOfCells(), 0);
    EXPECT_EQ(inspection.degenerateTriangleCount, 0);
    EXPECT_EQ(inspection.isolatedPointCount, 0);
    expectRaysHit(contourMesh, rays);
  }
}

TEST(quest_marching_cubes, uniform_sphere_resolution_sweep_regression)
{
  const std::array<StructuredDomainSpec, 3> resolutions {{
    {12, 8, 8, 0},
    {24, 16, 16, 0},
    {36, 24, 24, 0},
  }};

  axom::IndexType prevTriangleCount = -1;
  int prevUniquePointCount = -1;

  for(const StructuredDomainSpec& spec : resolutions)
  {
    SCOPED_TRACE(axom::fmt::format("dims=({}, {}, {})", spec.ni, spec.nj, spec.nk));

    UMesh contourMesh(3, mint::TRIANGLE);
    auto meshData = makeUniformSphereBlueprint(spec);
    buildContourMesh(meshData.mesh, "field", 0.0, contourMesh);

    const MeshInspection inspection = inspectContourMesh(contourMesh);
    EXPECT_GT(contourMesh.getNumberOfCells(), 0);
    EXPECT_GT(inspection.uniquePointCount, 0);
    EXPECT_EQ(inspection.degenerateTriangleCount, 0);
    EXPECT_EQ(inspection.isolatedPointCount, 0);
    EXPECT_EQ(inspection.duplicateTriangleCount, 0);
    EXPECT_FALSE(inspection.hasOpenBoundary());
    EXPECT_EQ(inspection.connectedComponentCount, 1);

    if(prevTriangleCount >= 0)
    {
      EXPECT_GT(inspection.triangleCount, prevTriangleCount);
      EXPECT_GT(inspection.uniquePointCount, prevUniquePointCount);
    }

    prevTriangleCount = inspection.triangleCount;
    prevUniquePointCount = inspection.uniquePointCount;
  }
}

TEST(quest_marching_cubes, single_cube_singular_cases)
{
  const std::array<std::pair<const char*, SingularCase>, 6> cases {{
    {"SIMPLE", SingularCase::SIMPLE},
    {"DOUBLE", SingularCase::DOUBLE},
    {"TRIPLE", SingularCase::TRIPLE},
    {"TUNNEL_OPEN", SingularCase::TUNNEL_OPEN},
    {"TUNNEL", SingularCase::TUNNEL},
    {"TUNNEL_DOUBLE", SingularCase::TUNNEL_DOUBLE},
  }};

  for(const auto& caseEntry : cases)
  {
    SCOPED_TRACE(caseEntry.first);

    CornerValues values {};
    const double isoValue = generatePredefinedSingular(caseEntry.second, values);
    const auto meshData = makeSingleCellStructuredCaseBlueprint(values, isoValue);
    UMesh contourMesh(3, mint::TRIANGLE);
    runStructuredCase(meshData, contourMesh, "mask");
    const MeshInspection inspection = inspectContourMesh(contourMesh);

    EXPECT_GT(contourMesh.getNumberOfCells(), 0);
    EXPECT_EQ(inspection.isolatedPointCount, 0);
    EXPECT_EQ(inspection.degenerateTriangleCount, 0);
    if(caseEntry.second == SingularCase::TUNNEL_OPEN)
    {
      EXPECT_TRUE(inspection.hasOpenBoundary());
    }
  }
}

TEST(quest_marching_cubes, single_cube_iso_edge_cases)
{
  const std::array<std::pair<const char*, IsoEdgeCase>, 4> cases {{
    {"SIMPLE", IsoEdgeCase::SIMPLE},
    {"SEPARATED", IsoEdgeCase::SEPARATED},
    {"WEDGE", IsoEdgeCase::WEDGE},
    {"EDGE_TUNNEL", IsoEdgeCase::EDGE_TUNNEL},
  }};

  for(const auto& caseEntry : cases)
  {
    SCOPED_TRACE(caseEntry.first);

    CornerValues values {};
    const double isoValue = generatePredefinedIsoEdge(caseEntry.second, values);
    const auto meshData = makeSingleCellStructuredCaseBlueprint(values, isoValue);
    UMesh contourMesh(3, mint::TRIANGLE);
    runStructuredCase(meshData, contourMesh, "mask");
    const MeshInspection inspection = inspectContourMesh(contourMesh);
    const std::vector<Point3D> uniquePoints = extractUniquePoints(contourMesh).points;

    EXPECT_GT(contourMesh.getNumberOfCells(), 0);
    EXPECT_EQ(inspection.isolatedPointCount, 0);
    EXPECT_EQ(inspection.degenerateTriangleCount, 0);

    if(caseEntry.second != IsoEdgeCase::SEPARATED)
    {
      expectContainsPoints(uniquePoints, cgalIsoVertices(values, isoValue));
    }
  }
}

TEST(quest_marching_cubes, single_cube_plane_cases)
{
  const std::array<std::pair<const char*, PlaneCase>, 4> cases {{
    {"SIMPLE", PlaneCase::SIMPLE},
    {"INTERSECTING", PlaneCase::INTERSECTING},
    {"PLANE_TUBES", PlaneCase::PLANE_TUBES},
    {"CROSS", PlaneCase::CROSS},
  }};

  for(const auto& caseEntry : cases)
  {
    SCOPED_TRACE(caseEntry.first);

    CornerValues values {};
    const double isoValue = generatePredefinedPlane(caseEntry.second, values);
    const auto meshData = makeSingleCellStructuredCaseBlueprint(values, isoValue);
    UMesh contourMesh(3, mint::TRIANGLE);
    runStructuredCase(meshData, contourMesh, "mask");
    const MeshInspection inspection = inspectContourMesh(contourMesh);
    const std::vector<Point3D> uniquePoints = extractUniquePoints(contourMesh).points;

    if(caseEntry.second == PlaneCase::SIMPLE)
    {
      // Axom's classic MC currently drops the all-iso face case entirely.
      EXPECT_EQ(contourMesh.getNumberOfCells(), 0);
      EXPECT_EQ(inspection.connectedComponentCount, 0);
      EXPECT_FALSE(inspection.hasOpenBoundary());
      continue;
    }

    EXPECT_GT(contourMesh.getNumberOfCells(), 0);
    EXPECT_EQ(inspection.isolatedPointCount, 0);
    EXPECT_EQ(inspection.degenerateTriangleCount, 0);

    if(caseEntry.second == PlaneCase::INTERSECTING)
    {
      expectContainsPoints(uniquePoints, cgalIsoVertices(values, isoValue));
    }
  }
}

TEST(quest_marching_cubes, embedded_case_topology_smoke)
{
  const auto runEmbeddedCase = [](const char* name, const CornerValues& values, double isoValue) {
    SCOPED_TRACE(name);

    const auto meshData = makeEmbeddedStructuredCaseBlueprint(values, isoValue);
    UMesh contourMesh(3, mint::TRIANGLE);
    runStructuredCase(meshData, contourMesh);
    const MeshInspection inspection = inspectContourMesh(contourMesh);

    EXPECT_GT(contourMesh.getNumberOfCells(), 0);
    EXPECT_EQ(inspection.isolatedPointCount, 0);
    EXPECT_EQ(inspection.degenerateTriangleCount, 0);
  };

  CornerValues ambiguousValues {};
  double isoValue = generatePredefinedInnerAmbiguity(AmbiguousCase::CONTOUR_6_VTS, ambiguousValues);
  runEmbeddedCase("AMBIGUOUS_CONTOUR_6_VTS", ambiguousValues, isoValue);

  CornerValues tunnelValues {};
  isoValue = generatePredefinedInnerAmbiguity(AmbiguousCase::MC_4_TUNNEL, tunnelValues);
  runEmbeddedCase("AMBIGUOUS_MC_4_TUNNEL", tunnelValues, isoValue);

  CornerValues singularValues {};
  isoValue = generatePredefinedSingular(SingularCase::TUNNEL, singularValues);
  runEmbeddedCase("SINGULAR_TUNNEL", singularValues, isoValue);

  CornerValues planeValues {};
  isoValue = generatePredefinedPlane(PlaneCase::INTERSECTING, planeValues);
  runEmbeddedCase("PLANE_INTERSECTING", planeValues, isoValue);
}

TEST(quest_marching_cubes, tracked_regression_iso_edge_separated_endpoint_snapping)
{
  CornerValues values {};
  const double isoValue = generatePredefinedIsoEdge(IsoEdgeCase::SEPARATED, values);
  const auto meshData = makeSingleCellStructuredCaseBlueprint(values, isoValue);
  UMesh contourMesh(3, mint::TRIANGLE);
  runStructuredCase(meshData, contourMesh, "mask");
  const std::vector<Point3D> uniquePoints = extractUniquePoints(contourMesh).points;

  ASSERT_GT(contourMesh.getNumberOfCells(), 0);

  // Tracked robustness regression: classic MarchingCubes should not emit
  // an iso-edge contour that misses the snapped endpoint vertices for this
  // CGAL-inspired separated iso-edge case.
  expectContainsPoints(uniquePoints, cgalIsoVertices(values, isoValue));
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
