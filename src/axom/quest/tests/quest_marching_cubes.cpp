// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "quest_marching_cubes.cpp requires conduit"
#endif

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/slic.hpp"

#include "conduit_blueprint.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;

using Point3D = primal::Point<double, 3>;
using Ray3D = primal::Ray<double, 3>;
using Triangle3D = primal::Triangle<double, 3>;
using UMesh = axom::mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

namespace
{

constexpr double X_MIN = 0.0;
constexpr double X_MAX = 5.0474;
constexpr double X_MID = 2.5237;
constexpr double Y_MIN = 0.0;
constexpr double Y_MAX = 1.0;
constexpr double Z_MIN = 0.0;
constexpr double Z_MAX = 1.0;

constexpr axom::IndexType DOM0_NI = 44;
constexpr axom::IndexType DOM0_NJ = 200;
constexpr axom::IndexType DOM0_NK = 3;

constexpr axom::IndexType DOM1_NI = 3;
constexpr axom::IndexType DOM1_NJ = 44;
constexpr axom::IndexType DOM1_NK = 200;

constexpr double SPHERE_CX = 2.5387;
constexpr double SPHERE_CY = 0.505;
constexpr double SPHERE_CZ = 0.34;
constexpr double SPHERE_R = 0.23;
constexpr double CONTOUR_VALUE = 0.5;
// This is intentionally higher than the first cut of the reproducer.
// The original Axom-side mirror under-sampled the per-zone volume fractions,
// which kept the seam-bottom nodal values near 0.75 instead of the
// ~0.774 exported by Ares. That small nodal-field mismatch moved the seam
// contour several 1e-3 in z and prevented a close Ares/Axom comparison.
// Raising the sub-zone sampling density makes the synthetic field converge
// much more closely to the Ares deck's accumulated nodal field.
constexpr int ZONE_SAMPLES_PER_AXIS = 24;
constexpr axom::IndexType PHONY_OFFSET = 2;

constexpr double SEAM_X = 2.5237;
constexpr double LEFT_X = 2.4937;
constexpr double RIGHT_X = 2.5537;
constexpr double LOW_Y = 0.475;
constexpr double HIGH_Y = 0.535;

struct DomainSpec
{
  axom::IndexType ni;
  axom::IndexType nj;
  axom::IndexType nk;
  int domainId;
  bool rotated;
};

struct ProbeSpec
{
  const char* name;
  Point3D point;
};

struct RayHit
{
  bool hit {false};
  Point3D point {Point3D::make_point(0.0, 0.0, 0.0)};
  double rayParam {std::numeric_limits<double>::max()};
  axom::IndexType triangleId {-1};
  axom::IndexType parentCellId {-1};
  axom::IndexType domainId {-1};
};

struct RayResults
{
  RayHit xCenter;
  RayHit yCenter;
  RayHit zTopSeam;
  RayHit zBottomSeam;
  RayHit zBelowSeam;
  RayHit zTopLeft;
  RayHit zBottomLeft;
  RayHit zBelowLeft;
  RayHit zTopRight;
  RayHit zBottomRight;
  RayHit zBelowRight;
  RayHit zTopHighY;
  RayHit zBottomHighY;
  RayHit zTopLowY;
  RayHit zBottomLowY;
};

struct RayCandidate
{
  axom::IndexType triangleId {-1};
  axom::IndexType parentCellId {-1};
  axom::IndexType domainId {-1};
  bool intersects {false};
  double rayParam {-1.0};
  bool aresIntersects {false};
  double aresRayParam {-1.0};
  Point3D pts[3];
};

axom::IndexType flattenNodeIndex(axom::IndexType i,
                                 axom::IndexType j,
                                 axom::IndexType k,
                                 axom::IndexType ni,
                                 axom::IndexType nj)
{
  return i + ni * (j + nj * k);
}

axom::IndexType extrapolatedLogicalIndex(axom::IndexType rawIdx, axom::IndexType activeCount)
{
  const axom::IndexType logicalIdx = rawIdx - PHONY_OFFSET;
  if(logicalIdx < 0)
  {
    return 0;
  }
  if(logicalIdx >= activeCount)
  {
    return activeCount - 1 + (logicalIdx - (activeCount - 1));
  }
  return logicalIdx;
}

std::string describePoint(const Point3D& pt)
{
  std::ostringstream oss;
  oss << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ")";
  return oss.str();
}

template <typename VecType>
std::string describeCoords(const VecType& v)
{
  std::ostringstream oss;
  oss << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
  return oss.str();
}

std::string describeHit(const RayHit& hit)
{
  std::ostringstream oss;
  oss << "hit=" << hit.hit;
  if(hit.hit)
  {
    oss << " point=" << describePoint(hit.point) << " t=" << hit.rayParam << " tri=" << hit.triangleId
        << " parentCell=" << hit.parentCellId << " domain=" << hit.domainId;
  }
  return oss.str();
}

std::string describeTriangle(const Point3D pts[3])
{
  std::ostringstream oss;
  oss << "v0=" << describePoint(pts[0]) << " v1=" << describePoint(pts[1])
      << " v2=" << describePoint(pts[2]);
  return oss.str();
}

bool aresStyleRayIntersects(const Ray3D& ray, const Triangle3D& tri, double& t)
{
  // Keep this helper aligned with the current Ares-side predicate.
  // Earlier debugging showed that using a positive epsilon in the edge sign
  // checks caused an orientation-dependent defect: exact edge hits could be
  // accepted for top-down rays and rejected for bottom-up rays on the same
  // triangle. The Ares code was patched to use exact zero for that
  // classification, and this mirror helper follows that behavior so the test
  // can compare like with like.
  constexpr double zero = 0.0;
  int kx = 0;
  int ky = 0;
  int kz = 0;

  double r[3];
  for(int i = 0; i < 3; ++i)
  {
    r[i] = std::abs(ray.direction()[i]);
  }

  if((r[2] >= r[0]) && (r[2] >= r[1]))
  {
    kz = 2;
  }
  else if((r[1] >= r[0]) && (r[1] >= r[2]))
  {
    kz = 1;
  }
  else
  {
    kz = 0;
  }

  kx = (kz + 1) % 3;
  ky = (kz + 2) % 3;
  if(ray.direction()[kz] < zero)
  {
    axom::utilities::swap(kx, ky);
  }

  const double invDir = 1.0 / ray.direction()[kz];
  const double shearX = invDir * ray.direction()[kx];
  const double shearY = invDir * ray.direction()[ky];

  const auto makeOffset = [&](int vertexIdx, int dim) {
    return tri[vertexIdx][dim] - ray.origin()[dim];
  };

  const double Ax = makeOffset(0, kx) - shearX * makeOffset(0, kz);
  const double Ay = makeOffset(0, ky) - shearY * makeOffset(0, kz);
  const double Bx = makeOffset(1, kx) - shearX * makeOffset(1, kz);
  const double By = makeOffset(1, ky) - shearY * makeOffset(1, kz);
  const double Cx = makeOffset(2, kx) - shearX * makeOffset(2, kz);
  const double Cy = makeOffset(2, ky) - shearY * makeOffset(2, kz);

  const double U = Cx * By - Cy * Bx;
  const double V = Ax * Cy - Ay * Cx;
  const double W = Bx * Ay - By * Ax;

  if((U < zero || V < zero || W < zero) && (U > zero || V > zero || W > zero))
  {
    return false;
  }

  const double det = U + V + W;
  if(axom::utilities::isNearlyEqual(det, zero))
  {
    return false;
  }

  const double Az = invDir * makeOffset(0, kz);
  const double Bz = invDir * makeOffset(1, kz);
  const double Cz = invDir * makeOffset(2, kz);
  t = U * Az + V * Bz + W * Cz;

  if(((t < zero) && !(det < zero)) || ((det < zero) && !(t < zero)))
  {
    return false;
  }

  t /= det;
  return true;
}

std::string sharedNodeKey(double x, double y, double z)
{
  std::ostringstream oss;
  // The strided/phony mesh mirrors Ares Blueprint export, which means the two
  // structured domains each carry their own copies of seam nodes. Ares then
  // sum-reduces those nodal accumulations across block boundaries before
  // running MarchingCubes. If we leave the Axom test with per-domain nodal
  // averages only, the seam values stay too low and the reproduced contour
  // sits measurably below the Ares contour.
  //
  // Quantize coordinates before keying them so physically identical seam nodes
  // from the two domains collapse reliably even if the decimal rendering of
  // the doubles differs in the last few bits.
  constexpr double scale = 1.0e12;
  const long long qx = static_cast<long long>(std::llround(x * scale));
  const long long qy = static_cast<long long>(std::llround(y * scale));
  const long long qz = static_cast<long long>(std::llround(z * scale));
  oss << qx << "|" << qy << "|" << qz;
  return oss.str();
}

void getStructuredFieldLayout(const conduit::Node& dom,
                              axom::IndexType& ni,
                              axom::IndexType& nj,
                              axom::IndexType& nk,
                              axom::IndexType offsets[3],
                              axom::IndexType strides[3]);

axom::IndexType flattenStructuredIndex(axom::IndexType i,
                                       axom::IndexType j,
                                       axom::IndexType k,
                                       const axom::IndexType offsets[3],
                                       const axom::IndexType strides[3]);

bool pointInsideSphere(double x, double y, double z)
{
  const double dx = x - SPHERE_CX;
  const double dy = y - SPHERE_CY;
  const double dz = z - SPHERE_CZ;
  return dx * dx + dy * dy + dz * dz <= SPHERE_R * SPHERE_R;
}

double sampleCellVolumeFraction(double x0, double x1, double y0, double y1, double z0, double z1)
{
  // The Ares deck builds a nodal field by accumulating per-zone volume
  // fractions into the surrounding vertices and then averaging. This helper
  // reproduces the per-zone contribution with a simple Cartesian sub-sampling
  // scheme so the Axom test can generate a comparable field without depending
  // on Ares internals.
  int insideCount = 0;
  const int totalCount = ZONE_SAMPLES_PER_AXIS * ZONE_SAMPLES_PER_AXIS * ZONE_SAMPLES_PER_AXIS;

  for(int kk = 0; kk < ZONE_SAMPLES_PER_AXIS; ++kk)
  {
    const double z = z0 + (static_cast<double>(kk) + 0.5) * (z1 - z0) / ZONE_SAMPLES_PER_AXIS;
    for(int jj = 0; jj < ZONE_SAMPLES_PER_AXIS; ++jj)
    {
      const double y = y0 + (static_cast<double>(jj) + 0.5) * (y1 - y0) / ZONE_SAMPLES_PER_AXIS;
      for(int ii = 0; ii < ZONE_SAMPLES_PER_AXIS; ++ii)
      {
        const double x = x0 + (static_cast<double>(ii) + 0.5) * (x1 - x0) / ZONE_SAMPLES_PER_AXIS;
        insideCount += pointInsideSphere(x, y, z) ? 1 : 0;
      }
    }
  }

  return static_cast<double>(insideCount) / static_cast<double>(totalCount);
}

double* createVertexField(conduit::Node& dom,
                          const std::string& fieldName,
                          axom::IndexType nnodes,
                          bool withPhonies,
                          axom::IndexType rawNi,
                          axom::IndexType rawNj,
                          axom::IndexType rawNk)
{
  AXOM_UNUSED_VAR(rawNk);
  conduit::Node& field = dom[std::string("fields/") + fieldName];
  field["association"] = "vertex";
  field["topology"] = "mesh";
  field["values"].set(conduit::DataType::float64(nnodes));
  if(withPhonies)
  {
    field["offsets"].set(std::vector<std::int32_t> {PHONY_OFFSET, PHONY_OFFSET, PHONY_OFFSET});
    field["strides"].set(std::vector<std::int32_t> {1,
                                                    static_cast<std::int32_t>(rawNi),
                                                    static_cast<std::int32_t>(rawNi * rawNj)});
  }
  return field["values"].as_double_ptr();
}

void setNodeCoordinates(const DomainSpec& spec,
                        axom::IndexType i,
                        axom::IndexType j,
                        axom::IndexType k,
                        double& x,
                        double& y,
                        double& z)
{
  if(!spec.rotated)
  {
    x = X_MIN + (X_MID - X_MIN) * static_cast<double>(i) / (spec.ni - 1);
    y = Y_MIN + (Y_MAX - Y_MIN) * static_cast<double>(j) / (spec.nj - 1);
    z = Z_MIN + (Z_MAX - Z_MIN) * static_cast<double>(k) / (spec.nk - 1);
    return;
  }

  z = Z_MIN + (Z_MAX - Z_MIN) * static_cast<double>(i) / (spec.ni - 1);
  x = X_MAX + (X_MID - X_MAX) * static_cast<double>(j) / (spec.nj - 1);
  y = Y_MAX + (Y_MIN - Y_MAX) * static_cast<double>(k) / (spec.nk - 1);
}

void populateDomain(conduit::Node& dom, const DomainSpec& spec, bool withPhonies)
{
  dom["state/domain_id"] = spec.domainId;

  dom["coordsets/coords/type"] = "explicit";
  dom["topologies/mesh/type"] = "structured";
  dom["topologies/mesh/coordset"] = "coords";
  dom["topologies/mesh/elements/dims/i"] = spec.ni - 1;
  dom["topologies/mesh/elements/dims/j"] = spec.nj - 1;
  dom["topologies/mesh/elements/dims/k"] = spec.nk - 1;

  const axom::IndexType rawNi = withPhonies ? spec.ni + PHONY_OFFSET + 1 : spec.ni;
  const axom::IndexType rawNj = withPhonies ? spec.nj + PHONY_OFFSET + 1 : spec.nj;
  const axom::IndexType rawNk = withPhonies ? spec.nk + PHONY_OFFSET + 1 : spec.nk;
  const axom::IndexType nnodes = rawNi * rawNj * rawNk;

  if(withPhonies)
  {
    dom["topologies/mesh/elements/dims/offsets"].set(
      std::vector<std::int32_t> {PHONY_OFFSET, PHONY_OFFSET, PHONY_OFFSET});
    dom["topologies/mesh/elements/dims/strides"].set(
      std::vector<std::int32_t> {1,
                                 static_cast<std::int32_t>(rawNi),
                                 static_cast<std::int32_t>(rawNi * rawNj)});
  }

  conduit::Node& xvals = dom["coordsets/coords/values/x"];
  conduit::Node& yvals = dom["coordsets/coords/values/y"];
  conduit::Node& zvals = dom["coordsets/coords/values/z"];
  xvals.set(conduit::DataType::float64(nnodes));
  yvals.set(conduit::DataType::float64(nnodes));
  zvals.set(conduit::DataType::float64(nnodes));

  double* x = xvals.as_double_ptr();
  double* y = yvals.as_double_ptr();
  double* z = zvals.as_double_ptr();

  for(axom::IndexType k = 0; k < rawNk; ++k)
  {
    for(axom::IndexType j = 0; j < rawNj; ++j)
    {
      for(axom::IndexType i = 0; i < rawNi; ++i)
      {
        const axom::IndexType idx = flattenNodeIndex(i, j, k, rawNi, rawNj);
        const axom::IndexType ii = withPhonies ? extrapolatedLogicalIndex(i, spec.ni) : i;
        const axom::IndexType jj = withPhonies ? extrapolatedLogicalIndex(j, spec.nj) : j;
        const axom::IndexType kk = withPhonies ? extrapolatedLogicalIndex(k, spec.nk) : k;
        setNodeCoordinates(spec, ii, jj, kk, x[idx], y[idx], z[idx]);
      }
    }
  }

  double* binaryField =
    createVertexField(dom, "indicator_binary", nnodes, withPhonies, rawNi, rawNj, rawNk);
  double* sampledField =
    createVertexField(dom, "indicator_sampled", nnodes, withPhonies, rawNi, rawNj, rawNk);
  double* sampledRawField =
    createVertexField(dom, "indicator_sampled_raw", nnodes, withPhonies, rawNi, rawNj, rawNk);
  double* sampledCountField =
    createVertexField(dom, "indicator_sampled_count", nnodes, withPhonies, rawNi, rawNj, rawNk);

  for(axom::IndexType idx = 0; idx < nnodes; ++idx)
  {
    binaryField[idx] = pointInsideSphere(x[idx], y[idx], z[idx]) ? 1.0 : 0.0;
    sampledField[idx] = 0.0;
    sampledRawField[idx] = 0.0;
    sampledCountField[idx] = 0.0;
  }

  for(axom::IndexType k = 0; k < spec.nk - 1; ++k)
  {
    for(axom::IndexType j = 0; j < spec.nj - 1; ++j)
    {
      for(axom::IndexType i = 0; i < spec.ni - 1; ++i)
      {
        double cellXMin = std::numeric_limits<double>::max();
        double cellYMin = std::numeric_limits<double>::max();
        double cellZMin = std::numeric_limits<double>::max();
        double cellXMax = -std::numeric_limits<double>::max();
        double cellYMax = -std::numeric_limits<double>::max();
        double cellZMax = -std::numeric_limits<double>::max();

        for(axom::IndexType dk = 0; dk <= 1; ++dk)
        {
          for(axom::IndexType dj = 0; dj <= 1; ++dj)
          {
            for(axom::IndexType di = 0; di <= 1; ++di)
            {
              const axom::IndexType nodeIdx =
                flattenNodeIndex(i + di + (withPhonies ? PHONY_OFFSET : 0),
                                 j + dj + (withPhonies ? PHONY_OFFSET : 0),
                                 k + dk + (withPhonies ? PHONY_OFFSET : 0),
                                 rawNi,
                                 rawNj);
              cellXMin = std::min(cellXMin, x[nodeIdx]);
              cellYMin = std::min(cellYMin, y[nodeIdx]);
              cellZMin = std::min(cellZMin, z[nodeIdx]);
              cellXMax = std::max(cellXMax, x[nodeIdx]);
              cellYMax = std::max(cellYMax, y[nodeIdx]);
              cellZMax = std::max(cellZMax, z[nodeIdx]);
            }
          }
        }

        // Each zone contributes the same sampled volume fraction to all eight
        // of its corner nodes, matching the Ares "accumulate then average"
        // pattern used to build the exported contour field.
        const double zoneFrac =
          sampleCellVolumeFraction(cellXMin, cellXMax, cellYMin, cellYMax, cellZMin, cellZMax);

        for(axom::IndexType dk = 0; dk <= 1; ++dk)
        {
          for(axom::IndexType dj = 0; dj <= 1; ++dj)
          {
            for(axom::IndexType di = 0; di <= 1; ++di)
            {
              const axom::IndexType nodeIdx =
                flattenNodeIndex(i + di + (withPhonies ? PHONY_OFFSET : 0),
                                 j + dj + (withPhonies ? PHONY_OFFSET : 0),
                                 k + dk + (withPhonies ? PHONY_OFFSET : 0),
                                 rawNi,
                                 rawNj);
              sampledRawField[nodeIdx] += zoneFrac;
              sampledCountField[nodeIdx] += 1.0;
            }
          }
        }
      }
    }
  }

  for(axom::IndexType idx = 0; idx < nnodes; ++idx)
  {
    sampledField[idx] =
      sampledCountField[idx] > 0.0 ? sampledRawField[idx] / sampledCountField[idx] : 0.0;
  }
}

void reduceSharedVertexFields(conduit::Node& mesh)
{
  struct NodeRef
  {
    int domainIdx;
    axom::IndexType flatIdx;
  };

  std::map<std::string, std::vector<NodeRef>> sharedNodes;
  const int numDomains = conduit::blueprint::mesh::number_of_domains(mesh);
  for(int domainIdx = 0; domainIdx < numDomains; ++domainIdx)
  {
    conduit::Node& dom = mesh.child(domainIdx);
    const double* x = dom.fetch_existing("coordsets/coords/values/x").as_double_ptr();
    const double* y = dom.fetch_existing("coordsets/coords/values/y").as_double_ptr();
    const double* z = dom.fetch_existing("coordsets/coords/values/z").as_double_ptr();

    axom::IndexType ni = 0;
    axom::IndexType nj = 0;
    axom::IndexType nk = 0;
    axom::IndexType offsets[3];
    axom::IndexType strides[3];
    getStructuredFieldLayout(dom, ni, nj, nk, offsets, strides);

    for(axom::IndexType k = 0; k < nk; ++k)
    {
      for(axom::IndexType j = 0; j < nj; ++j)
      {
        for(axom::IndexType i = 0; i < ni; ++i)
        {
          const axom::IndexType idx = flattenStructuredIndex(i, j, k, offsets, strides);
          sharedNodes[sharedNodeKey(x[idx], y[idx], z[idx])].push_back({domainIdx, idx});
        }
      }
    }
  }

  for(const auto& entry : sharedNodes)
  {
    if(entry.second.size() < 2)
    {
      continue;
    }

    // This is the critical Ares-mirroring step.
    //
    // Ares does a communication/sum-reduction of the nodal accumulation and
    // nodal contribution counts across shared block boundaries before
    // contouring. Without reproducing that here, the seam nodes are averaged
    // using only one domain's local contributions, which produced the earlier
    // Axom/Ares mismatch: Axom saw seam-bottom values near 0.75 while Ares
    // exported about 0.7738686, and the contour/ray hits shifted accordingly.
    double rawSum = 0.0;
    double countSum = 0.0;
    for(const auto& ref : entry.second)
    {
      const conduit::Node& dom = mesh.child(ref.domainIdx);
      const double* raw = dom.fetch_existing("fields/indicator_sampled_raw/values").as_double_ptr();
      const double* count =
        dom.fetch_existing("fields/indicator_sampled_count/values").as_double_ptr();
      rawSum += raw[ref.flatIdx];
      countSum += count[ref.flatIdx];
    }

    for(const auto& ref : entry.second)
    {
      conduit::Node& dom = mesh.child(ref.domainIdx);
      double* raw = dom.fetch_existing("fields/indicator_sampled_raw/values").as_double_ptr();
      double* count = dom.fetch_existing("fields/indicator_sampled_count/values").as_double_ptr();
      raw[ref.flatIdx] = rawSum;
      count[ref.flatIdx] = countSum;
    }
  }

  for(int domainIdx = 0; domainIdx < numDomains; ++domainIdx)
  {
    conduit::Node& dom = mesh.child(domainIdx);
    double* sampled = dom.fetch_existing("fields/indicator_sampled/values").as_double_ptr();
    const double* raw = dom.fetch_existing("fields/indicator_sampled_raw/values").as_double_ptr();
    const double* count = dom.fetch_existing("fields/indicator_sampled_count/values").as_double_ptr();
    const axom::IndexType nnodes =
      dom.fetch_existing("fields/indicator_sampled/values").dtype().number_of_elements();

    for(axom::IndexType idx = 0; idx < nnodes; ++idx)
    {
      // Recompute the final nodal field after the shared-node reduction so the
      // Blueprint passed into MarchingCubes matches the post-communication
      // field that Ares actually contours.
      sampled[idx] = count[idx] > 0.0 ? raw[idx] / count[idx] : 0.0;
    }
  }
}

conduit::Node buildRotatedSeamSphereMesh()
{
  conduit::Node mesh;
  populateDomain(mesh.append(), {DOM0_NI, DOM0_NJ, DOM0_NK, 0, false}, false);
  populateDomain(mesh.append(), {DOM1_NI, DOM1_NJ, DOM1_NK, 1, true}, false);
  reduceSharedVertexFields(mesh);
  return mesh;
}

conduit::Node buildRotatedSeamSphereMeshWithPhonies()
{
  conduit::Node mesh;
  populateDomain(mesh.append(), {DOM0_NI, DOM0_NJ, DOM0_NK, 0, false}, true);
  populateDomain(mesh.append(), {DOM1_NI, DOM1_NJ, DOM1_NK, 1, true}, true);
  reduceSharedVertexFields(mesh);
  return mesh;
}

void getStructuredFieldLayout(const conduit::Node& dom,
                              axom::IndexType& ni,
                              axom::IndexType& nj,
                              axom::IndexType& nk,
                              axom::IndexType offsets[3],
                              axom::IndexType strides[3])
{
  const conduit::Node& dims = dom.fetch_existing("topologies/mesh/elements/dims");
  ni = dims.fetch_existing("i").to_int64() + 1;
  nj = dims.fetch_existing("j").to_int64() + 1;
  nk = dims.fetch_existing("k").to_int64() + 1;

  offsets[0] = offsets[1] = offsets[2] = 0;
  strides[0] = 1;
  strides[1] = ni;
  strides[2] = ni * nj;

  if(dims.has_child("offsets"))
  {
    const std::int32_t* vals = dims.fetch_existing("offsets").as_int32_ptr();
    offsets[0] = vals[0];
    offsets[1] = vals[1];
    offsets[2] = vals[2];
  }

  if(dims.has_child("strides"))
  {
    const std::int32_t* vals = dims.fetch_existing("strides").as_int32_ptr();
    strides[0] = vals[0];
    strides[1] = vals[1];
    strides[2] = vals[2];
  }
}

axom::IndexType flattenStructuredIndex(axom::IndexType i,
                                       axom::IndexType j,
                                       axom::IndexType k,
                                       const axom::IndexType offsets[3],
                                       const axom::IndexType strides[3])
{
  return (offsets[0] + i) * strides[0] + (offsets[1] + j) * strides[1] +
    (offsets[2] + k) * strides[2];
}

void logStructuredFieldSummary(const conduit::Node& mesh, const std::string& fieldName)
{
  const int numDomains = conduit::blueprint::mesh::number_of_domains(mesh);
  for(int domainIdx = 0; domainIdx < numDomains; ++domainIdx)
  {
    const conduit::Node& dom = mesh.child(domainIdx);
    const conduit::Node& xvals = dom.fetch_existing("coordsets/coords/values/x");
    const conduit::Node& yvals = dom.fetch_existing("coordsets/coords/values/y");
    const conduit::Node& zvals = dom.fetch_existing("coordsets/coords/values/z");
    const conduit::Node& fieldVals =
      dom.fetch_existing(std::string("fields/") + fieldName + "/values");
    axom::IndexType ni = 0;
    axom::IndexType nj = 0;
    axom::IndexType nk = 0;
    axom::IndexType offsets[3];
    axom::IndexType strides[3];
    getStructuredFieldLayout(dom, ni, nj, nk, offsets, strides);

    const double* x = xvals.as_double_ptr();
    const double* y = yvals.as_double_ptr();
    const double* z = zvals.as_double_ptr();
    const double* values = fieldVals.as_double_ptr();

    double minVal = std::numeric_limits<double>::max();
    double maxVal = -std::numeric_limits<double>::max();
    axom::IndexType sampleI = 0;
    axom::IndexType sampleJ = 0;
    axom::IndexType sampleK = 0;
    double bestMetric = std::numeric_limits<double>::max();

    for(axom::IndexType k = 0; k < nk; ++k)
    {
      for(axom::IndexType j = 0; j < nj; ++j)
      {
        for(axom::IndexType i = 0; i < ni; ++i)
        {
          const axom::IndexType idx = flattenStructuredIndex(i, j, k, offsets, strides);
          minVal = std::min(minVal, values[idx]);
          maxVal = std::max(maxVal, values[idx]);

          const double metric =
            std::abs(x[idx] - SEAM_X) + std::abs(y[idx] - SPHERE_CY) + std::abs(z[idx] - Z_MIN);
          if(metric < bestMetric)
          {
            bestMetric = metric;
            sampleI = i;
            sampleJ = j;
            sampleK = k;
          }
        }
      }
    }

    std::ostringstream oss;
    oss << "[mc-debug] testField field=" << fieldName
        << " domain=" << dom.fetch_existing("state/domain_id").to_int() << " dims=(" << ni << ", "
        << nj << ", " << nk << ") offsets=(" << offsets[0] << ", " << offsets[1] << ", "
        << offsets[2] << ") strides=(" << strides[0] << ", " << strides[1] << ", " << strides[2]
        << ") min=" << minVal << " max=" << maxVal;
    SLIC_INFO(oss.str());

    for(axom::IndexType lineK = 0; lineK < nk; ++lineK)
    {
      const axom::IndexType idx = flattenStructuredIndex(sampleI, sampleJ, lineK, offsets, strides);
      std::ostringstream line;
      line << "[mc-debug]   testField sampleLine field=" << fieldName
           << " domain=" << dom.fetch_existing("state/domain_id").to_int() << " logical=("
           << sampleI << ", " << sampleJ << ", " << lineK << ") coord=(" << x[idx] << ", " << y[idx]
           << ", " << z[idx] << ") value=" << values[idx];
      SLIC_INFO(line.str());
    }

    AXOM_UNUSED_VAR(sampleK);
  }
}

std::vector<ProbeSpec> makeDebugProbes()
{
  return {
    {"seam_bottom", Point3D::make_point(SEAM_X, SPHERE_CY, Z_MIN)},
    {"seam_mid", Point3D::make_point(SEAM_X, SPHERE_CY, 0.5)},
    {"seam_top", Point3D::make_point(SEAM_X, SPHERE_CY, Z_MAX)},
    {"left_bottom", Point3D::make_point(LEFT_X, SPHERE_CY, Z_MIN)},
    {"right_bottom", Point3D::make_point(RIGHT_X, SPHERE_CY, Z_MIN)},
    {"highy_bottom", Point3D::make_point(SEAM_X, HIGH_Y, Z_MIN)},
    {"lowy_bottom", Point3D::make_point(SEAM_X, LOW_Y, Z_MIN)},
  };
}

const conduit::Node& findDomainById(const conduit::Node& mesh, int domainId)
{
  const int numDomains = conduit::blueprint::mesh::number_of_domains(mesh);
  for(int domainIdx = 0; domainIdx < numDomains; ++domainIdx)
  {
    const conduit::Node& dom = mesh.child(domainIdx);
    if(dom.fetch_existing("state/domain_id").to_int() == domainId)
    {
      return dom;
    }
  }
  SLIC_ERROR(axom::fmt::format("Could not find domain {}", domainId));
  return mesh.child(0);
}

void logParentCellDetails(const conduit::Node& mesh,
                          const std::string& fieldName,
                          const std::string& label,
                          int domainId,
                          axom::IndexType parentCellId)
{
  const conduit::Node& dom = findDomainById(mesh, domainId);
  const double* x = dom.fetch_existing("coordsets/coords/values/x").as_double_ptr();
  const double* y = dom.fetch_existing("coordsets/coords/values/y").as_double_ptr();
  const double* z = dom.fetch_existing("coordsets/coords/values/z").as_double_ptr();
  const double* values =
    dom.fetch_existing(std::string("fields/") + fieldName + "/values").as_double_ptr();

  axom::IndexType ni = 0;
  axom::IndexType nj = 0;
  axom::IndexType nk = 0;
  axom::IndexType offsets[3];
  axom::IndexType strides[3];
  getStructuredFieldLayout(dom, ni, nj, nk, offsets, strides);

  const axom::IndexType nzI = ni - 1;
  const axom::IndexType nzJ = nj - 1;
  const axom::IndexType zonePlane = nzI * nzJ;
  const axom::IndexType k = parentCellId / zonePlane;
  const axom::IndexType rem = parentCellId % zonePlane;
  const axom::IndexType j = rem / nzI;
  const axom::IndexType i = rem % nzI;

  std::ostringstream hdr;
  hdr << "[mc-debug] cellDetail label=" << label << " field=" << fieldName << " domain=" << domainId
      << " parentCell=" << parentCellId << " logicalZone=(" << i << ", " << j << ", " << k << ")";
  SLIC_INFO(hdr.str());

  int aboveIsoCount = 0;
  for(axom::IndexType dk = 0; dk <= 1; ++dk)
  {
    for(axom::IndexType dj = 0; dj <= 1; ++dj)
    {
      for(axom::IndexType di = 0; di <= 1; ++di)
      {
        const axom::IndexType nodeI = i + di;
        const axom::IndexType nodeJ = j + dj;
        const axom::IndexType nodeK = k + dk;
        const axom::IndexType idx = flattenStructuredIndex(nodeI, nodeJ, nodeK, offsets, strides);
        const bool aboveIso = values[idx] >= CONTOUR_VALUE;
        aboveIsoCount += aboveIso ? 1 : 0;

        std::ostringstream oss;
        oss << "[mc-debug]   cellCorner domain=" << domainId << " parentCell=" << parentCellId
            << " logical=(" << nodeI << ", " << nodeJ << ", " << nodeK << ") coord=(" << x[idx]
            << ", " << y[idx] << ", " << z[idx] << ") value=" << values[idx]
            << " aboveIso=" << aboveIso;
        SLIC_INFO(oss.str());
      }
    }
  }

  SLIC_INFO(
    axom::fmt::format("[mc-debug]   cellCaseSummary domain={} parentCell={} aboveIsoCorners={}",
                      domainId,
                      parentCellId,
                      aboveIsoCount));
}

void logStructuredFieldProbes(const conduit::Node& mesh,
                              const std::string& fieldName,
                              const std::string& label)
{
  const auto probes = makeDebugProbes();
  const int numDomains = conduit::blueprint::mesh::number_of_domains(mesh);
  for(int domainIdx = 0; domainIdx < numDomains; ++domainIdx)
  {
    const conduit::Node& dom = mesh.child(domainIdx);
    const double* x = dom.fetch_existing("coordsets/coords/values/x").as_double_ptr();
    const double* y = dom.fetch_existing("coordsets/coords/values/y").as_double_ptr();
    const double* z = dom.fetch_existing("coordsets/coords/values/z").as_double_ptr();
    const double* values =
      dom.fetch_existing(std::string("fields/") + fieldName + "/values").as_double_ptr();

    axom::IndexType ni = 0;
    axom::IndexType nj = 0;
    axom::IndexType nk = 0;
    axom::IndexType offsets[3];
    axom::IndexType strides[3];
    getStructuredFieldLayout(dom, ni, nj, nk, offsets, strides);

    for(const auto& probe : probes)
    {
      axom::IndexType bestI = 0;
      axom::IndexType bestJ = 0;
      axom::IndexType bestK = 0;
      double bestMetric = std::numeric_limits<double>::max();
      for(axom::IndexType k = 0; k < nk; ++k)
      {
        for(axom::IndexType j = 0; j < nj; ++j)
        {
          for(axom::IndexType i = 0; i < ni; ++i)
          {
            const axom::IndexType idx = flattenStructuredIndex(i, j, k, offsets, strides);
            const double metric = std::abs(x[idx] - probe.point[0]) +
              std::abs(y[idx] - probe.point[1]) + std::abs(z[idx] - probe.point[2]);
            if(metric < bestMetric)
            {
              bestMetric = metric;
              bestI = i;
              bestJ = j;
              bestK = k;
            }
          }
        }
      }

      const axom::IndexType idx = flattenStructuredIndex(bestI, bestJ, bestK, offsets, strides);
      std::ostringstream oss;
      oss << "[mc-debug] probe label=" << label << " field=" << fieldName
          << " domain=" << dom.fetch_existing("state/domain_id").to_int() << " probe=" << probe.name
          << " target=" << describePoint(probe.point) << " logical=(" << bestI << ", " << bestJ
          << ", " << bestK << ") coord=(" << x[idx] << ", " << y[idx] << ", " << z[idx]
          << ") value=" << values[idx] << " manhattanDist=" << bestMetric;
      SLIC_INFO(oss.str());
    }
  }
}

void logStructuredFieldLine(const conduit::Node& mesh,
                            const std::string& fieldName,
                            const std::string& label,
                            const std::string& lineName,
                            const Point3D& target,
                            int varyingDim)
{
  const int numDomains = conduit::blueprint::mesh::number_of_domains(mesh);
  for(int domainIdx = 0; domainIdx < numDomains; ++domainIdx)
  {
    const conduit::Node& dom = mesh.child(domainIdx);
    const double* x = dom.fetch_existing("coordsets/coords/values/x").as_double_ptr();
    const double* y = dom.fetch_existing("coordsets/coords/values/y").as_double_ptr();
    const double* z = dom.fetch_existing("coordsets/coords/values/z").as_double_ptr();
    const double* values =
      dom.fetch_existing(std::string("fields/") + fieldName + "/values").as_double_ptr();

    axom::IndexType ni = 0;
    axom::IndexType nj = 0;
    axom::IndexType nk = 0;
    axom::IndexType offsets[3];
    axom::IndexType strides[3];
    getStructuredFieldLayout(dom, ni, nj, nk, offsets, strides);

    axom::IndexType logical[3] = {0, 0, 0};
    double bestMetric = std::numeric_limits<double>::max();

    for(axom::IndexType k = 0; k < nk; ++k)
    {
      for(axom::IndexType j = 0; j < nj; ++j)
      {
        for(axom::IndexType i = 0; i < ni; ++i)
        {
          const axom::IndexType idx = flattenStructuredIndex(i, j, k, offsets, strides);
          const double coords[3] = {x[idx], y[idx], z[idx]};
          double metric = 0.0;
          for(int dim = 0; dim < 3; ++dim)
          {
            if(dim != varyingDim)
            {
              metric += std::abs(coords[dim] - target[dim]);
            }
          }
          if(metric < bestMetric)
          {
            bestMetric = metric;
            logical[0] = i;
            logical[1] = j;
            logical[2] = k;
          }
        }
      }
    }

    const axom::IndexType counts[3] = {ni, nj, nk};
    std::ostringstream hdr;
    hdr << "[mc-debug] line label=" << label << " field=" << fieldName
        << " domain=" << dom.fetch_existing("state/domain_id").to_int() << " line=" << lineName
        << " varyingDim=" << varyingDim << " target=" << describePoint(target) << " anchor=("
        << logical[0] << ", " << logical[1] << ", " << logical[2] << ") metric=" << bestMetric;
    SLIC_INFO(hdr.str());

    for(axom::IndexType v = 0; v < counts[varyingDim]; ++v)
    {
      axom::IndexType ijk[3] = {logical[0], logical[1], logical[2]};
      ijk[varyingDim] = v;
      const axom::IndexType idx = flattenStructuredIndex(ijk[0], ijk[1], ijk[2], offsets, strides);
      std::ostringstream oss;
      oss << "[mc-debug]   lineSample label=" << label << " field=" << fieldName
          << " domain=" << dom.fetch_existing("state/domain_id").to_int() << " line=" << lineName
          << " logical=(" << ijk[0] << ", " << ijk[1] << ", " << ijk[2] << ") coord=(" << x[idx]
          << ", " << y[idx] << ", " << z[idx] << ") value=" << values[idx];
      SLIC_INFO(oss.str());
    }
  }
}

RayHit findFirstRayHit(const UMesh& mesh, const Ray3D& ray)
{
  RayHit result;
  const axom::IndexType ncells = mesh.getNumberOfCells();
  const axom::IndexType* cellIds = mesh.getFieldPtr<axom::IndexType>("cellId", mint::CELL_CENTERED);
  const axom::IndexType* domainIds =
    mesh.getFieldPtr<axom::IndexType>("domainId", mint::CELL_CENTERED);

  for(axom::IndexType icell = 0; icell < ncells; ++icell)
  {
    const axom::IndexType* nodeIds = mesh.getCellNodeIDs(icell);

    Point3D pts[3];
    for(int n = 0; n < 3; ++n)
    {
      mesh.getNode(nodeIds[n], pts[n].data());
    }

    const Triangle3D tri {pts[0], pts[1], pts[2]};
    double t = 0.0;
    if(primal::intersect(tri, ray, t) && t >= 0.0 && t < result.rayParam)
    {
      result.hit = true;
      result.rayParam = t;
      result.point = ray.at(t);
      result.triangleId = icell;
      result.parentCellId = cellIds[icell];
      result.domainId = domainIds[icell];
    }
  }

  return result;
}

std::vector<RayCandidate> findRayCandidates(const UMesh& mesh, const Ray3D& ray)
{
  std::vector<RayCandidate> candidates;
  const axom::IndexType ncells = mesh.getNumberOfCells();
  const axom::IndexType* cellIds = mesh.getFieldPtr<axom::IndexType>("cellId", mint::CELL_CENTERED);
  const axom::IndexType* domainIds =
    mesh.getFieldPtr<axom::IndexType>("domainId", mint::CELL_CENTERED);

  constexpr double EPS = 1e-12;
  for(axom::IndexType icell = 0; icell < ncells; ++icell)
  {
    const axom::IndexType* nodeIds = mesh.getCellNodeIDs(icell);

    RayCandidate cand;
    cand.triangleId = icell;
    cand.parentCellId = cellIds[icell];
    cand.domainId = domainIds[icell];
    for(int n = 0; n < 3; ++n)
    {
      mesh.getNode(nodeIds[n], cand.pts[n].data());
    }

    double minX = cand.pts[0][0];
    double maxX = cand.pts[0][0];
    double minY = cand.pts[0][1];
    double maxY = cand.pts[0][1];
    double minZ = cand.pts[0][2];
    double maxZ = cand.pts[0][2];
    for(int n = 1; n < 3; ++n)
    {
      minX = std::min(minX, cand.pts[n][0]);
      maxX = std::max(maxX, cand.pts[n][0]);
      minY = std::min(minY, cand.pts[n][1]);
      maxY = std::max(maxY, cand.pts[n][1]);
      minZ = std::min(minZ, cand.pts[n][2]);
      maxZ = std::max(maxZ, cand.pts[n][2]);
    }

    const auto origin = ray.origin();
    const auto dir = ray.direction();
    bool plausible = true;
    if(std::abs(dir[0]) < EPS)
    {
      plausible = plausible && (origin[0] >= minX - EPS && origin[0] <= maxX + EPS);
    }
    if(std::abs(dir[1]) < EPS)
    {
      plausible = plausible && (origin[1] >= minY - EPS && origin[1] <= maxY + EPS);
    }
    if(std::abs(dir[2]) < EPS)
    {
      plausible = plausible && (origin[2] >= minZ - EPS && origin[2] <= maxZ + EPS);
    }

    if(!plausible)
    {
      continue;
    }

    const Triangle3D tri {cand.pts[0], cand.pts[1], cand.pts[2]};
    double t = 0.0;
    cand.intersects = primal::intersect(tri, ray, t) && t >= 0.0;
    cand.rayParam = cand.intersects ? t : -1.0;
    double aresT = 0.0;
    cand.aresIntersects = aresStyleRayIntersects(ray, tri, aresT) && aresT >= 0.0;
    cand.aresRayParam = cand.aresIntersects ? aresT : -1.0;
    candidates.push_back(cand);
  }

  std::sort(candidates.begin(), candidates.end(), [](const RayCandidate& a, const RayCandidate& b) {
    if(a.parentCellId != b.parentCellId)
    {
      return a.parentCellId < b.parentCellId;
    }
    return a.triangleId < b.triangleId;
  });
  return candidates;
}

void logRayCandidates(const UMesh& mesh,
                      const Ray3D& ray,
                      const std::string& label,
                      const std::string& rayName)
{
  const auto candidates = findRayCandidates(mesh, ray);
  std::ostringstream hdr;
  hdr << "[mc-debug] " << label << " candidateRay name=" << rayName
      << " origin=" << describePoint(ray.origin()) << " dir=" << describeCoords(ray.direction())
      << " candidates=" << candidates.size();
  SLIC_INFO(hdr.str());

  for(std::size_t i = 0; i < candidates.size(); ++i)
  {
    const auto& cand = candidates[i];
    std::ostringstream oss;
    oss << "[mc-debug]   candidate[" << i << "] tri=" << cand.triangleId
        << " parentCell=" << cand.parentCellId << " domain=" << cand.domainId
        << " intersects=" << cand.intersects << " t=" << cand.rayParam
        << " aresIntersects=" << cand.aresIntersects << " aresT=" << cand.aresRayParam << " "
        << describeTriangle(cand.pts);
    SLIC_INFO(oss.str());
  }
}

void logContourSummary(const UMesh& mesh, const std::string& label)
{
  primal::BoundingBox<double, 3> bbox;
  Point3D pt;
  for(axom::IndexType inode = 0; inode < mesh.getNumberOfNodes(); ++inode)
  {
    mesh.getNode(inode, pt.data());
    bbox.addPoint(pt);
  }

  std::ostringstream oss;
  oss << "[mc-debug] contour label=" << label << " cells=" << mesh.getNumberOfCells()
      << " nodes=" << mesh.getNumberOfNodes() << " bboxMin=" << describePoint(bbox.getMin())
      << " bboxMax=" << describePoint(bbox.getMax());
  SLIC_INFO(oss.str());
}

RayResults runMarchingCubesCase(const conduit::Node& mesh,
                                const std::string& fieldName,
                                const std::string& label)
{
  quest::MarchingCubes marchingCubes(axom::runtime_policy::Policy::seq,
                                     axom::execution_space<axom::SEQ_EXEC>::allocatorID(),
                                     quest::MarchingCubesDataParallelism::byPolicy);

  marchingCubes.setMesh(mesh, "mesh");
  marchingCubes.setFunctionField(fieldName);
  marchingCubes.computeIsocontour(CONTOUR_VALUE);

  EXPECT_GT(marchingCubes.getContourCellCount(), 0) << label;
  EXPECT_GT(marchingCubes.getContourNodeCount(), 0) << label;

  UMesh contourMesh(3, mint::TRIANGLE);
  marchingCubes.populateContourMesh(contourMesh, "cellId", "domainId");
  logContourSummary(contourMesh, label);

  const Ray3D xCenter {{X_MAX, SPHERE_CY, SPHERE_CZ}, {-1.0, 0.0, 0.0}};
  const Ray3D yCenter {{SPHERE_CX, Y_MAX, SPHERE_CZ}, {0.0, -1.0, 0.0}};
  const Ray3D zTopSeam {{SEAM_X, SPHERE_CY, Z_MAX}, {0.0, 0.0, -1.0}};
  const Ray3D zBottomSeam {{SEAM_X, SPHERE_CY, Z_MIN}, {0.0, 0.0, 1.0}};
  const Ray3D zBelowSeam {{SEAM_X, SPHERE_CY, -0.02}, {0.0, 0.0, 1.0}};
  const Ray3D zTopLeft {{LEFT_X, SPHERE_CY, Z_MAX}, {0.0, 0.0, -1.0}};
  const Ray3D zBottomLeft {{LEFT_X, SPHERE_CY, Z_MIN}, {0.0, 0.0, 1.0}};
  const Ray3D zBelowLeft {{LEFT_X, SPHERE_CY, -0.02}, {0.0, 0.0, 1.0}};
  const Ray3D zTopRight {{RIGHT_X, SPHERE_CY, Z_MAX}, {0.0, 0.0, -1.0}};
  const Ray3D zBottomRight {{RIGHT_X, SPHERE_CY, Z_MIN}, {0.0, 0.0, 1.0}};
  const Ray3D zBelowRight {{RIGHT_X, SPHERE_CY, -0.02}, {0.0, 0.0, 1.0}};
  const Ray3D zTopHighY {{SEAM_X, HIGH_Y, Z_MAX}, {0.0, 0.0, -1.0}};
  const Ray3D zBottomHighY {{SEAM_X, HIGH_Y, Z_MIN}, {0.0, 0.0, 1.0}};
  const Ray3D zTopLowY {{SEAM_X, LOW_Y, Z_MAX}, {0.0, 0.0, -1.0}};
  const Ray3D zBottomLowY {{SEAM_X, LOW_Y, Z_MIN}, {0.0, 0.0, 1.0}};

  RayResults results;
  results.xCenter = findFirstRayHit(contourMesh, xCenter);
  results.yCenter = findFirstRayHit(contourMesh, yCenter);
  results.zTopSeam = findFirstRayHit(contourMesh, zTopSeam);
  results.zBottomSeam = findFirstRayHit(contourMesh, zBottomSeam);
  results.zBelowSeam = findFirstRayHit(contourMesh, zBelowSeam);
  results.zTopLeft = findFirstRayHit(contourMesh, zTopLeft);
  results.zBottomLeft = findFirstRayHit(contourMesh, zBottomLeft);
  results.zBelowLeft = findFirstRayHit(contourMesh, zBelowLeft);
  results.zTopRight = findFirstRayHit(contourMesh, zTopRight);
  results.zBottomRight = findFirstRayHit(contourMesh, zBottomRight);
  results.zBelowRight = findFirstRayHit(contourMesh, zBelowRight);
  results.zTopHighY = findFirstRayHit(contourMesh, zTopHighY);
  results.zBottomHighY = findFirstRayHit(contourMesh, zBottomHighY);
  results.zTopLowY = findFirstRayHit(contourMesh, zTopLowY);
  results.zBottomLowY = findFirstRayHit(contourMesh, zBottomLowY);

  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_x_center " + describeHit(results.xCenter));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_y_center " + describeHit(results.yCenter));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_top_seam " + describeHit(results.zTopSeam));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_bottom_seam " +
            describeHit(results.zBottomSeam));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_below_seam " +
            describeHit(results.zBelowSeam));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_top_left " + describeHit(results.zTopLeft));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_bottom_left " +
            describeHit(results.zBottomLeft));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_below_left " +
            describeHit(results.zBelowLeft));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_top_right " + describeHit(results.zTopRight));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_bottom_right " +
            describeHit(results.zBottomRight));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_below_right " +
            describeHit(results.zBelowRight));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_top_highy " + describeHit(results.zTopHighY));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_bottom_highy " +
            describeHit(results.zBottomHighY));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_top_lowy " + describeHit(results.zTopLowY));
  SLIC_INFO(std::string("[mc-debug] ") + label + " ray_z_bottom_lowy " +
            describeHit(results.zBottomLowY));

  logRayCandidates(contourMesh, zTopSeam, label, "ray_z_top_seam");
  logRayCandidates(contourMesh, zBottomSeam, label, "ray_z_bottom_seam");
  logRayCandidates(contourMesh, zBelowSeam, label, "ray_z_below_seam");
  logRayCandidates(contourMesh, zTopLeft, label, "ray_z_top_left");
  logRayCandidates(contourMesh, zBottomLeft, label, "ray_z_bottom_left");
  logRayCandidates(contourMesh, zTopRight, label, "ray_z_top_right");
  logRayCandidates(contourMesh, zBottomRight, label, "ray_z_bottom_right");
  logRayCandidates(contourMesh, zTopHighY, label, "ray_z_top_highy");
  logRayCandidates(contourMesh, zBottomHighY, label, "ray_z_bottom_highy");
  logRayCandidates(contourMesh, zTopLowY, label, "ray_z_top_lowy");
  logRayCandidates(contourMesh, zBottomLowY, label, "ray_z_bottom_lowy");

  const RayHit detailedHits[] = {results.zTopSeam,
                                 results.zBottomSeam,
                                 results.zTopLeft,
                                 results.zBottomLeft,
                                 results.zTopRight,
                                 results.zBottomRight,
                                 results.zTopHighY,
                                 results.zBottomHighY,
                                 results.zTopLowY,
                                 results.zBottomLowY};
  std::vector<std::pair<int, axom::IndexType>> loggedCells;
  for(const auto& hit : detailedHits)
  {
    if(!hit.hit)
    {
      continue;
    }

    const auto key = std::make_pair(static_cast<int>(hit.domainId), hit.parentCellId);
    if(std::find(loggedCells.begin(), loggedCells.end(), key) == loggedCells.end())
    {
      loggedCells.push_back(key);
      logParentCellDetails(mesh, fieldName, label, key.first, key.second);
    }
  }

  return results;
}

void expectLateralHits(const RayResults& results)
{
  EXPECT_TRUE(results.xCenter.hit) << describeHit(results.xCenter);
  EXPECT_TRUE(results.yCenter.hit) << describeHit(results.yCenter);

  if(results.xCenter.hit)
  {
    EXPECT_GE(results.xCenter.point[0], 2.60);
    EXPECT_LE(results.xCenter.point[0], 2.70);
    EXPECT_NEAR(results.xCenter.point[1], SPHERE_CY, 1e-6);
    EXPECT_NEAR(results.xCenter.point[2], SPHERE_CZ, 1e-6);
  }

  if(results.yCenter.hit)
  {
    EXPECT_NEAR(results.yCenter.point[0], SPHERE_CX, 1e-6);
    EXPECT_GE(results.yCenter.point[1], 0.60);
    EXPECT_LE(results.yCenter.point[1], 0.65);
    EXPECT_NEAR(results.yCenter.point[2], SPHERE_CZ, 1e-6);
  }
}

}  // namespace

TEST(quest_marching_cubes, rotated_seam_sampled_field_matches_ares_reproducer)
{
  conduit::Node compactMesh = buildRotatedSeamSphereMesh();
  conduit::Node compactInfo;
  EXPECT_TRUE(conduit::blueprint::mesh::verify(compactMesh, compactInfo)) << compactInfo.to_yaml();

  logStructuredFieldSummary(compactMesh, "indicator_sampled");
  logStructuredFieldProbes(compactMesh, "indicator_sampled", "compact");
  logStructuredFieldLine(compactMesh,
                         "indicator_sampled",
                         "compact",
                         "seam_zline",
                         Point3D::make_point(SEAM_X, SPHERE_CY, Z_MIN),
                         2);

  const RayResults compactResults =
    runMarchingCubesCase(compactMesh, "indicator_sampled", "rotated_seam_sampled_compact");
  SLIC_INFO(std::string("[mc-debug] compact bottom seam ") + describeHit(compactResults.zBottomSeam));

  conduit::Node mesh = buildRotatedSeamSphereMeshWithPhonies();
  conduit::Node info;
  EXPECT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();

  logStructuredFieldSummary(mesh, "indicator_binary");
  logStructuredFieldSummary(mesh, "indicator_sampled_raw");
  logStructuredFieldSummary(mesh, "indicator_sampled_count");
  logStructuredFieldSummary(mesh, "indicator_sampled");
  logStructuredFieldProbes(mesh, "indicator_sampled_raw", "strided");
  logStructuredFieldProbes(mesh, "indicator_sampled_count", "strided");
  logStructuredFieldProbes(mesh, "indicator_sampled", "strided");
  logStructuredFieldLine(mesh,
                         "indicator_sampled",
                         "strided",
                         "domain0_seam_zline",
                         Point3D::make_point(SEAM_X, SPHERE_CY, Z_MIN),
                         2);
  logStructuredFieldLine(mesh,
                         "indicator_sampled",
                         "strided",
                         "domain1_seam_yline_z0",
                         Point3D::make_point(SEAM_X, SPHERE_CY, Z_MIN),
                         2);

  const RayResults sampledResults =
    runMarchingCubesCase(mesh, "indicator_sampled", "rotated_seam_sampled_strided");

  expectLateralHits(sampledResults);

  // Final mirrored expectation:
  // after fixing the Ares probe edge-hit asymmetry and matching Ares's nodal
  // field construction more closely on the Axom side, both codes show the same
  // qualitative behavior for this setup. The seam-directed rays from above,
  // below, and slightly below all intersect the upper seam surface rather than
  // a lower-cap crossing.
  EXPECT_TRUE(sampledResults.zTopSeam.hit) << describeHit(sampledResults.zTopSeam);
  EXPECT_TRUE(sampledResults.zBottomSeam.hit) << describeHit(sampledResults.zBottomSeam);
  EXPECT_TRUE(sampledResults.zBelowSeam.hit) << describeHit(sampledResults.zBelowSeam);

  EXPECT_TRUE(sampledResults.zTopLeft.hit) << describeHit(sampledResults.zTopLeft);
  EXPECT_TRUE(sampledResults.zBottomLeft.hit) << describeHit(sampledResults.zBottomLeft);
  EXPECT_TRUE(sampledResults.zBelowLeft.hit) << describeHit(sampledResults.zBelowLeft);

  EXPECT_TRUE(sampledResults.zTopRight.hit) << describeHit(sampledResults.zTopRight);
  EXPECT_TRUE(sampledResults.zBottomRight.hit) << describeHit(sampledResults.zBottomRight);
  EXPECT_TRUE(sampledResults.zBelowRight.hit) << describeHit(sampledResults.zBelowRight);

  EXPECT_TRUE(sampledResults.zTopHighY.hit) << describeHit(sampledResults.zTopHighY);
  EXPECT_TRUE(sampledResults.zBottomHighY.hit) << describeHit(sampledResults.zBottomHighY);

  EXPECT_TRUE(sampledResults.zTopLowY.hit) << describeHit(sampledResults.zTopLowY);
  EXPECT_TRUE(sampledResults.zBottomLowY.hit) << describeHit(sampledResults.zBottomLowY);

  if(sampledResults.zTopSeam.hit)
  {
    EXPECT_NEAR(sampledResults.zTopSeam.point[0], SEAM_X, 1e-6);
    EXPECT_NEAR(sampledResults.zTopSeam.point[1], SPHERE_CY, 1e-6);
    EXPECT_GE(sampledResults.zTopSeam.point[2], 0.42);
    EXPECT_LE(sampledResults.zTopSeam.point[2], 0.44);
  }

  if(sampledResults.zBottomSeam.hit)
  {
    EXPECT_NEAR(sampledResults.zBottomSeam.point[0], SEAM_X, 1e-6);
    EXPECT_NEAR(sampledResults.zBottomSeam.point[1], SPHERE_CY, 1e-6);
    EXPECT_NEAR(sampledResults.zBottomSeam.point[2], sampledResults.zTopSeam.point[2], 1e-9);
  }

  if(sampledResults.zBelowSeam.hit)
  {
    EXPECT_NEAR(sampledResults.zBelowSeam.point[0], SEAM_X, 1e-6);
    EXPECT_NEAR(sampledResults.zBelowSeam.point[1], SPHERE_CY, 1e-6);
    EXPECT_NEAR(sampledResults.zBelowSeam.point[2], sampledResults.zTopSeam.point[2], 1e-9);
  }

  if(sampledResults.zTopLeft.hit)
  {
    EXPECT_NEAR(sampledResults.zTopLeft.point[0], LEFT_X, 1e-6);
    EXPECT_NEAR(sampledResults.zTopLeft.point[1], SPHERE_CY, 1e-6);
    EXPECT_GE(sampledResults.zTopLeft.point[2], 0.39);
    EXPECT_LE(sampledResults.zTopLeft.point[2], 0.42);
  }

  if(sampledResults.zBottomLeft.hit)
  {
    EXPECT_GT(sampledResults.zBottomLeft.point[2], 0.35);
    EXPECT_NEAR(sampledResults.zBottomLeft.point[2], sampledResults.zTopLeft.point[2], 1e-9);
  }

  if(sampledResults.zBelowLeft.hit)
  {
    EXPECT_GT(sampledResults.zBelowLeft.point[2], 0.35);
    EXPECT_NEAR(sampledResults.zBelowLeft.point[2], sampledResults.zTopLeft.point[2], 1e-9);
  }

  if(sampledResults.zTopRight.hit)
  {
    EXPECT_NEAR(sampledResults.zTopRight.point[0], RIGHT_X, 1e-6);
    EXPECT_NEAR(sampledResults.zTopRight.point[1], SPHERE_CY, 1e-6);
    EXPECT_GE(sampledResults.zTopRight.point[2], 0.41);
    EXPECT_LE(sampledResults.zTopRight.point[2], 0.43);
  }

  if(sampledResults.zBottomRight.hit)
  {
    EXPECT_GT(sampledResults.zBottomRight.point[2], 0.35);
    EXPECT_NEAR(sampledResults.zBottomRight.point[2], sampledResults.zTopRight.point[2], 1e-9);
  }

  if(sampledResults.zBelowRight.hit)
  {
    EXPECT_GT(sampledResults.zBelowRight.point[2], 0.35);
    EXPECT_NEAR(sampledResults.zBelowRight.point[2], sampledResults.zTopRight.point[2], 1e-9);
  }

  if(sampledResults.zTopHighY.hit)
  {
    EXPECT_NEAR(sampledResults.zTopHighY.point[0], SEAM_X, 1e-6);
    EXPECT_NEAR(sampledResults.zTopHighY.point[1], HIGH_Y, 1e-6);
    EXPECT_GE(sampledResults.zTopHighY.point[2], 0.41);
    EXPECT_LE(sampledResults.zTopHighY.point[2], 0.43);
  }

  if(sampledResults.zBottomHighY.hit)
  {
    EXPECT_NEAR(sampledResults.zBottomHighY.point[0], SEAM_X, 1e-6);
    EXPECT_NEAR(sampledResults.zBottomHighY.point[1], HIGH_Y, 1e-6);
    EXPECT_NEAR(sampledResults.zBottomHighY.point[2], sampledResults.zTopHighY.point[2], 1e-9);
  }

  if(sampledResults.zTopLowY.hit)
  {
    EXPECT_NEAR(sampledResults.zTopLowY.point[0], SEAM_X, 1e-6);
    EXPECT_NEAR(sampledResults.zTopLowY.point[1], LOW_Y, 1e-6);
    EXPECT_GE(sampledResults.zTopLowY.point[2], 0.41);
    EXPECT_LE(sampledResults.zTopLowY.point[2], 0.43);
  }

  if(sampledResults.zBottomLowY.hit)
  {
    EXPECT_NEAR(sampledResults.zBottomLowY.point[0], SEAM_X, 1e-6);
    EXPECT_NEAR(sampledResults.zBottomLowY.point[1], LOW_Y, 1e-6);
    EXPECT_NEAR(sampledResults.zBottomLowY.point[2], sampledResults.zTopLowY.point[2], 1e-9);
  }
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
