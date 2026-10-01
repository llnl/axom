// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file quest_marching_cubes.cpp
 *
 * @brief Tests for quest::MarchingCubes on structured Blueprint meshes.
 *
 * Most checks use linear fields. Marching Cubes interpolates linearly along
 * straight cell edges, so every contour node of a linear field lies exactly on
 * its zero plane, even on warped meshes. The expected facet counts for
 * axis-aligned planes follow from the number of cells the plane crosses.
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/quest/MarchingCubes.hpp"

#include "quest_marching_cubes_testing_helpers.hpp"

#include "gtest/gtest.h"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <map>
#include <string>
#include <vector>

namespace mctest = axom::quest::testing::marching_cubes;

namespace
{
using RuntimePolicy = axom::runtime_policy::Policy;
using DataParallelism = axom::quest::MarchingCubesDataParallelism;
using Vec3 = mctest::PlanarField::VectorType;

constexpr double POSITION_TOL = 1e-12;

//---------------------------------------------------------------------------
// Runtime policies and memory
//---------------------------------------------------------------------------

std::vector<RuntimePolicy> enabledPolicies()
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

int allocatorForPolicy(RuntimePolicy policy)
{
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  if(policy == RuntimePolicy::cuda)
  {
    return axom::execution_space<axom::CUDA_EXEC<256>>::allocatorID();
  }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  if(policy == RuntimePolicy::hip)
  {
    return axom::execution_space<axom::HIP_EXEC<256>>::allocatorID();
  }
#endif
  AXOM_UNUSED_VAR(policy);
  return mctest::hostAllocatorID();
}

std::string policyName(RuntimePolicy policy) { return axom::runtime_policy::policyToName(policy); }

//---------------------------------------------------------------------------
// Meshes
//---------------------------------------------------------------------------

//! @brief Shift the explicit coordinates of @a dom by @a dx along x.
void translateX(conduit::Node& dom, double dx)
{
  auto* x = dom["coordsets/coords/values/x"].as_float64_ptr();
  const auto n = dom["coordsets/coords/values/x"].dtype().number_of_elements();
  for(conduit::index_t i = 0; i < n; ++i)
  {
    x[i] += dx;
  }
}

/*!
 * @brief Build a multi-domain mesh from @a ndom unit structured domains placed
 *        side by side along x, sampling the planar field @a f on each.
 *
 * Domain @a d covers [d, d+1] x [0,1]^(DIM-1). When @a setDomainIds is true,
 * domain @a d gets \c state/domain_id = 10 + d, except the last domain, which
 * keeps the default so both id sources are exercised.
 */
template <int DIM>
void buildMultiDomain(conduit::Node& mdMesh,
                      int ndom,
                      int n,
                      const mctest::PlanarField& f,
                      bool setDomainIds)
{
  mdMesh.reset();
  for(int d = 0; d < ndom; ++d)
  {
    conduit::Node& dom = mdMesh.append();
    // The field is sampled after translation, so it sees global coordinates.
    mctest::buildStructured<DIM>(dom, n, [](double, double, double) { return 0.0; }, "fcn");
    translateX(dom, double(d));
    mctest::addVertexField<DIM>(dom, f, "fcn");
    if(setDomainIds && d + 1 < ndom)
    {
      dom["state/domain_id"] = 10 + d;
    }
  }
}

//! @brief Wrap a single domain as a one-domain multi-domain mesh.
void wrapAsMultiDomain(conduit::Node& mdMesh, const conduit::Node& dom)
{
  mdMesh.reset();
  mdMesh.append().set_external(dom);
}

/*!
 * @brief Host-side access to a domain's real (non-ghost) nodes by logical index.
 *
 * Honors \c elements/dims/offsets and \c strides when present.
 */
template <int DIM>
struct DomainGeometry
{
  explicit DomainGeometry(const conduit::Node& dom)
  {
    const conduit::Node& dims = dom.fetch_existing("topologies/mesh/elements/dims");
    cells[0] = dims.fetch_existing("i").to_int();
    cells[1] = dims.fetch_existing("j").to_int();
    cells[2] = DIM == 3 ? dims.fetch_existing("k").to_int() : 1;

    int stride = 1;
    for(int d = 0; d < DIM; ++d)
    {
      offsets[d] = 0;
      strides[d] = stride;
      stride *= cells[d] + 1;
    }
    if(dims.has_child("offsets"))
    {
      const auto acc = dims.fetch_existing("offsets").as_int_accessor();
      for(int d = 0; d < DIM; ++d)
      {
        offsets[d] = acc[d];
      }
    }
    if(dims.has_child("strides"))
    {
      const auto acc = dims.fetch_existing("strides").as_int_accessor();
      for(int d = 0; d < DIM; ++d)
      {
        strides[d] = acc[d];
      }
    }

    const conduit::Node& values = dom.fetch_existing("coordsets/coords/values");
    coords[0] = values.fetch_existing("x").as_float64_ptr();
    coords[1] = values.fetch_existing("y").as_float64_ptr();
    coords[2] = DIM == 3 ? values.fetch_existing("z").as_float64_ptr() : nullptr;
  }

  axom::IndexType cellCount() const
  {
    return axom::IndexType(cells[0]) * cells[1] * (DIM == 3 ? cells[2] : 1);
  }

  //! @brief Logical indices of the cell with i-fastest flat index @a cellId.
  std::array<int, 3> cellIndices(axom::IndexType cellId) const
  {
    std::array<int, 3> ijk {0, 0, 0};
    ijk[0] = int(cellId % cells[0]);
    ijk[1] = int((cellId / cells[0]) % cells[1]);
    ijk[2] = DIM == 3 ? int(cellId / (axom::IndexType(cells[0]) * cells[1])) : 0;
    return ijk;
  }

  double nodeCoord(int axis, int i, int j, int k) const
  {
    const axom::IndexType idx = axom::IndexType(i + offsets[0]) * strides[0] +
      axom::IndexType(j + offsets[1]) * strides[1] +
      (DIM == 3 ? axom::IndexType(k + offsets[2]) * strides[2] : 0);
    return coords[axis][idx];
  }

  //! @brief Bounding box of cell @a cellId, from its corner nodes.
  axom::primal::BoundingBox<double, DIM> cellBox(axom::IndexType cellId) const
  {
    const auto ijk = cellIndices(cellId);
    axom::primal::BoundingBox<double, DIM> box;
    for(int dk = 0; dk < (DIM == 3 ? 2 : 1); ++dk)
    {
      for(int dj = 0; dj < 2; ++dj)
      {
        for(int di = 0; di < 2; ++di)
        {
          axom::primal::Point<double, DIM> p;
          for(int a = 0; a < DIM; ++a)
          {
            p[a] = nodeCoord(a, ijk[0] + di, ijk[1] + dj, ijk[2] + dk);
          }
          box.addPoint(p);
        }
      }
    }
    return box;
  }

  int cells[3];
  int offsets[3];
  int strides[3];
  const double* coords[3];
};

//---------------------------------------------------------------------------
// Contour output on the host
//---------------------------------------------------------------------------

struct HostContour
{
  axom::IndexType facetCount {0};
  axom::IndexType nodeCount {0};
  axom::Array<axom::IndexType, 2> corners;
  axom::Array<double, 2> coords;
  axom::Array<axom::IndexType> parents;
  axom::Array<axom::IndexType> domains;
};

HostContour hostContour(const axom::quest::MarchingCubes& mc)
{
  const int host = mctest::hostAllocatorID();
  HostContour out;
  out.facetCount = mc.getContourCellCount();
  out.nodeCount = mc.getContourNodeCount();
  out.corners = axom::Array<axom::IndexType, 2>(mc.getContourFacetCorners(), host);
  out.coords = axom::Array<double, 2>(mc.getContourNodeCoords(), host);
  out.parents = axom::Array<axom::IndexType>(mc.getContourFacetParents(), host);
  out.domains = axom::Array<axom::IndexType>(mc.getContourFacetDomainIds(), host);
  return out;
}

/*!
 * @brief Check the array shapes and the legacy backend's connectivity contract.
 *
 * The legacy backend does not weld vertices: each facet owns DIM nodes,
 * and each node belongs to exactly one facet.
 */
template <int DIM>
void checkUnweldedConnectivity(const HostContour& c)
{
  EXPECT_EQ(c.nodeCount, DIM * c.facetCount);
  ASSERT_EQ(c.corners.shape()[0], c.facetCount);
  ASSERT_EQ(c.coords.shape()[0], c.nodeCount);
  ASSERT_EQ(c.parents.size(), c.facetCount);
  ASSERT_EQ(c.domains.size(), c.facetCount);

  std::vector<int> uses(c.nodeCount, 0);
  for(axom::IndexType f = 0; f < c.facetCount; ++f)
  {
    for(int d = 0; d < DIM; ++d)
    {
      const auto node = c.corners(f, d);
      ASSERT_GE(node, 0);
      ASSERT_LT(node, c.nodeCount);
      ++uses[node];
    }
  }
  EXPECT_TRUE(std::all_of(uses.begin(), uses.end(), [](int u) { return u == 1; }));
}

//! @brief Check that every contour node lies on the plane @a f.
template <int DIM>
void checkNodesOnPlane(const HostContour& c, const mctest::PlanarField& f)
{
  for(axom::IndexType n = 0; n < c.nodeCount; ++n)
  {
    const double z = DIM == 3 ? c.coords(n, 2) : 0.0;
    ASSERT_NEAR(f(c.coords(n, 0), c.coords(n, 1), z), 0.0, POSITION_TOL) << "node " << n;
  }
}

/*!
 * @brief Check that every facet's corners lie in its parent cell.
 *
 * @a domains maps each facet domain id to the host copy of that domain.
 * This ties together the parent-id, domain-id, connectivity, and coordinate
 * arrays, so an offset error in any of them fails the check.
 */
template <int DIM>
void checkFacetsInParentCells(const HostContour& c,
                              const std::map<axom::IndexType, const conduit::Node*>& domains)
{
  std::map<axom::IndexType, DomainGeometry<DIM>> geometry;
  for(const auto& [id, dom] : domains)
  {
    geometry.emplace(id, DomainGeometry<DIM>(*dom));
  }

  for(axom::IndexType f = 0; f < c.facetCount; ++f)
  {
    const auto it = geometry.find(c.domains[f]);
    ASSERT_NE(it, geometry.end()) << "facet " << f << " has unknown domain id " << c.domains[f];
    const auto& geom = it->second;
    ASSERT_GE(c.parents[f], 0);
    ASSERT_LT(c.parents[f], geom.cellCount());

    auto box = geom.cellBox(c.parents[f]);
    box.expand(POSITION_TOL);
    for(int d = 0; d < DIM; ++d)
    {
      axom::primal::Point<double, DIM> p;
      for(int a = 0; a < DIM; ++a)
      {
        p[a] = c.coords(c.corners(f, d), a);
      }
      ASSERT_TRUE(box.contains(p)) << "facet " << f << " corner " << d << " at " << p
                                   << " is outside parent cell " << c.parents[f] << " " << box;
    }
  }
}

//! @brief Map domain ids to the domains of a multi-domain mesh, as MarchingCubes assigns them.
std::map<axom::IndexType, const conduit::Node*> domainsById(const conduit::Node& mdMesh)
{
  std::map<axom::IndexType, const conduit::Node*> rval;
  for(conduit::index_t d = 0; d < mdMesh.number_of_children(); ++d)
  {
    const conduit::Node& dom = mdMesh.child(d);
    const axom::IndexType id = dom.has_path("state/domain_id")
      ? dom.fetch_existing("state/domain_id").to_int64()
      : axom::IndexType(d);
    rval[id] = &dom;
  }
  return rval;
}

/*!
 * @brief Run MarchingCubes on host mesh @a hostMesh with each policy and scan strategy,
 *        calling @a check with the host contour.
 */
template <typename CheckFunction>
void forEachConfiguration(const conduit::Node& hostMesh,
                          const std::string& maskField,
                          double contourValue,
                          CheckFunction&& check)
{
  for(auto policy : enabledPolicies())
  {
    const int allocatorID = allocatorForPolicy(policy);
    conduit::Node mesh;
    mctest::copyBlueprintToPolicy(mesh, hostMesh, policy, allocatorID);

    for(auto parallelism : {DataParallelism::hybridParallel, DataParallelism::fullParallel})
    {
      SCOPED_TRACE(axom::fmt::format("policy {}, data parallelism {}",
                                     policyName(policy),
                                     static_cast<int>(parallelism)));
      axom::quest::MarchingCubes mc(policy, allocatorID, parallelism);
      mc.setMesh(mesh, "mesh", maskField);
      mc.setFunctionField("fcn");
      mc.computeIsocontour(contourValue);
      check(mc, hostContour(mc));
    }
  }
}

}  // end anonymous namespace

//---------------------------------------------------------------------------
// Tests
//---------------------------------------------------------------------------

template <int DIM>
void testAxisAlignedPlane()
{
  constexpr int n = 8;
  constexpr double xPlane = 0.3;  // Between nodes 2/8 and 3/8.
  const mctest::PlanarField f(Vec3 {1., 0., 0.}, xPlane);

  conduit::Node dom, mdMesh;
  mctest::buildStructured<DIM>(dom, n, f, "fcn");
  wrapAsMultiDomain(mdMesh, dom);

  // The plane crosses one column of cells. In 2D each crossed cell yields one segment;
  // in 3D each crossed hex yields a quadrilateral split into two triangles.
  const axom::IndexType expectedFacets = DIM == 2 ? n : 2 * n * n;
  const int crossedColumn = int(xPlane * n);

  forEachConfiguration(mdMesh, "", 0.0, [&](const auto&, const HostContour& c) {
    EXPECT_EQ(c.facetCount, expectedFacets);
    checkUnweldedConnectivity<DIM>(c);
    checkNodesOnPlane<DIM>(c, f);
    checkFacetsInParentCells<DIM>(c, domainsById(mdMesh));

    const DomainGeometry<DIM> geom(dom);
    for(axom::IndexType i = 0; i < c.facetCount; ++i)
    {
      EXPECT_EQ(geom.cellIndices(c.parents[i])[0], crossedColumn);
      EXPECT_EQ(c.domains[i], 0);
    }
  });
}

TEST(quest_marching_cubes, axis_aligned_plane_2d) { testAxisAlignedPlane<2>(); }
TEST(quest_marching_cubes, axis_aligned_plane_3d) { testAxisAlignedPlane<3>(); }

/*!
 * @brief A smooth in-plane warp for 2D meshes.
 *
 * SinusoidalWarp displaces x and y only through sin(pi z),
 * so it is the identity when z == 0. This warp vanishes on the boundary of [0,1]^2.
 */
struct PlanarSinusoidalWarp
{
  double amp;

  void operator()(double& x, double& y, double& /*z*/) const
  {
    const double dx = amp * std::sin(M_PI * x) * std::sin(2 * M_PI * y);
    const double dy = amp * std::sin(2 * M_PI * x) * std::sin(M_PI * y);
    x += dx;
    y += dy;
  }
};

template <int DIM>
void testObliquePlaneOnWarpedMesh()
{
  constexpr int n = 10;
  const mctest::PlanarField f(mctest::PlanarField::PointType {0.4, 0.55, 0.45},
                              Vec3 {1., 2., DIM == 3 ? 3. : 0.});

  conduit::Node dom, mdMesh;
  if constexpr(DIM == 2)
  {
    mctest::buildStructured<DIM>(dom, n, f, "fcn", PlanarSinusoidalWarp {0.03});
  }
  else
  {
    mctest::buildStructured<DIM>(dom, n, f, "fcn", mctest::SinusoidalWarp {0.03});
  }
  wrapAsMultiDomain(mdMesh, dom);

  forEachConfiguration(mdMesh, "", 0.0, [&](const auto&, const HostContour& c) {
    EXPECT_GT(c.facetCount, 0);
    checkUnweldedConnectivity<DIM>(c);
    checkNodesOnPlane<DIM>(c, f);
    checkFacetsInParentCells<DIM>(c, domainsById(mdMesh));
  });
}

TEST(quest_marching_cubes, oblique_plane_warped_mesh_2d) { testObliquePlaneOnWarpedMesh<2>(); }
TEST(quest_marching_cubes, oblique_plane_warped_mesh_3d) { testObliquePlaneOnWarpedMesh<3>(); }

template <int DIM>
void testMultiDomain()
{
  constexpr int n = 6;
  constexpr int ndom = 3;
  // An oblique plane that crosses all three domains.
  const mctest::PlanarField f(mctest::PlanarField::PointType {1.5, 0.5, 0.5},
                              Vec3 {1., -2.5, DIM == 3 ? 0.5 : 0.});

  conduit::Node mdMesh;
  buildMultiDomain<DIM>(mdMesh, ndom, n, f, true);

  // Reference: contour each domain by itself.
  std::vector<axom::IndexType> perDomainFacets;
  for(int d = 0; d < ndom; ++d)
  {
    conduit::Node one;
    wrapAsMultiDomain(one, mdMesh.child(d));
    axom::quest::MarchingCubes mc(RuntimePolicy::seq,
                                  mctest::hostAllocatorID(),
                                  DataParallelism::byPolicy);
    mc.setMesh(one, "mesh");
    mc.setFunctionField("fcn");
    mc.computeIsocontour(0.0);
    perDomainFacets.push_back(mc.getContourCellCount());
    EXPECT_GT(perDomainFacets.back(), 0) << "domain " << d << " is not crossed";
  }

  const auto byId = domainsById(mdMesh);
  forEachConfiguration(mdMesh, "", 0.0, [&](const auto&, const HostContour& c) {
    checkUnweldedConnectivity<DIM>(c);
    checkNodesOnPlane<DIM>(c, f);
    checkFacetsInParentCells<DIM>(c, byId);

    // Facets are grouped by domain in domain order.
    const std::vector<axom::IndexType> ids {10, 11, 2};
    axom::IndexType f0 = 0;
    for(int d = 0; d < ndom; ++d)
    {
      for(axom::IndexType i = f0; i < f0 + perDomainFacets[d]; ++i)
      {
        ASSERT_EQ(c.domains[i], ids[d]) << "facet " << i;
      }
      f0 += perDomainFacets[d];
    }
    EXPECT_EQ(c.facetCount, f0);
  });
}

TEST(quest_marching_cubes, multi_domain_2d) { testMultiDomain<2>(); }
TEST(quest_marching_cubes, multi_domain_3d) { testMultiDomain<3>(); }

//! @brief Convert an \c int32 index array at @a path in @a node to \c int64.
void convertToInt64(conduit::Node& node, const std::string& path)
{
  conduit::Node converted;
  node.fetch_existing(path).to_int64_array(converted);
  node[path].set(converted);
}

/*!
 * @brief Contour a strided mesh with ghost layers.
 *
 * With @a int64Metadata, the layout offsets and strides are \c int64 instead of \c int32,
 * covering both index types that MeshViewUtil accepts. On device policies, copyBlueprintToPolicy()
 * moves these arrays to device memory.
 */
template <int DIM>
void testStridedMesh(bool int64Metadata)
{
  constexpr int n = 8;
  constexpr int pad = 2;
  constexpr double xPlane = 0.3;
  const mctest::PlanarField f(Vec3 {1., 0., 0.}, xPlane);

  conduit::Node dom, mdMesh;
  mctest::buildStridedStructured<DIM>(dom, n, pad, f, "fcn");
  if(int64Metadata)
  {
    for(const char* path : {"topologies/mesh/elements/dims/offsets",
                            "topologies/mesh/elements/dims/strides",
                            "fields/fcn/offsets",
                            "fields/fcn/strides"})
    {
      convertToInt64(dom, path);
      ASSERT_TRUE(dom.fetch_existing(path).dtype().is_int64());
    }
  }
  wrapAsMultiDomain(mdMesh, dom);

  // The ghost layers continue the grid, so any ghost cell that entered the
  // contour would add facets beyond the real column.
  const axom::IndexType expectedFacets = DIM == 2 ? n : 2 * n * n;
  forEachConfiguration(mdMesh, "", 0.0, [&](const auto&, const HostContour& c) {
    EXPECT_EQ(c.facetCount, expectedFacets);
    checkUnweldedConnectivity<DIM>(c);
    checkNodesOnPlane<DIM>(c, f);
    checkFacetsInParentCells<DIM>(c, domainsById(mdMesh));
  });
}

TEST(quest_marching_cubes, strided_mesh_2d) { testStridedMesh<2>(false); }
TEST(quest_marching_cubes, strided_mesh_3d) { testStridedMesh<3>(false); }
TEST(quest_marching_cubes, strided_mesh_int64_metadata_2d) { testStridedMesh<2>(true); }
TEST(quest_marching_cubes, strided_mesh_int64_metadata_3d) { testStridedMesh<3>(true); }

template <int DIM>
void testMask()
{
  constexpr int n = 8;
  constexpr double xPlane = 0.3;
  const mctest::PlanarField f(Vec3 {1., 0., 0.}, xPlane);

  conduit::Node dom, mdMesh;
  mctest::buildStructured<DIM>(dom, n, f, "fcn");
  // Mask value 1 in the lower half of the domain (j < n/2), 0 elsewhere.
  mctest::addCellField<DIM>(dom, [](int, int j, int) { return j < n / 2 ? 1 : 0; }, "mask");
  wrapAsMultiDomain(mdMesh, dom);
  const DomainGeometry<DIM> geom(dom);

  const axom::IndexType allFacets = DIM == 2 ? n : 2 * n * n;
  for(int maskValue : {1, 0})
  {
    SCOPED_TRACE(axom::fmt::format("mask value {}", maskValue));
    for(auto policy : enabledPolicies())
    {
      const int allocatorID = allocatorForPolicy(policy);
      conduit::Node mesh;
      mctest::copyBlueprintToPolicy(mesh, mdMesh, policy, allocatorID);

      axom::quest::MarchingCubes mc(policy, allocatorID, DataParallelism::byPolicy);
      mc.setMesh(mesh, "mesh", "mask");
      mc.setFunctionField("fcn");
      mc.setMaskValue(maskValue);
      mc.computeIsocontour(0.0);
      const HostContour c = hostContour(mc);

      EXPECT_EQ(c.facetCount, allFacets / 2) << policyName(policy);
      checkUnweldedConnectivity<DIM>(c);
      checkNodesOnPlane<DIM>(c, f);
      for(axom::IndexType i = 0; i < c.facetCount; ++i)
      {
        const bool lowerHalf = geom.cellIndices(c.parents[i])[1] < n / 2;
        EXPECT_EQ(lowerHalf, maskValue == 1) << "facet " << i;
      }
    }
  }
}

TEST(quest_marching_cubes, mask_2d) { testMask<2>(); }
TEST(quest_marching_cubes, mask_3d) { testMask<3>(); }

template <int DIM>
void testAppendClearAndRelinquish()
{
  constexpr int n = 8;
  const mctest::PlanarField f(Vec3 {1., 0., 0.}, 0.0);
  const double c1 = 0.3, c2 = 0.55;  // Isovalues x = 0.3 and x = 0.55.
  const axom::IndexType perContour = DIM == 2 ? n : 2 * n * n;

  conduit::Node dom, mdMesh;
  mctest::buildStructured<DIM>(dom, n, f, "fcn");
  wrapAsMultiDomain(mdMesh, dom);

  for(auto policy : enabledPolicies())
  {
    SCOPED_TRACE(policyName(policy));
    const int allocatorID = allocatorForPolicy(policy);
    conduit::Node mesh;
    mctest::copyBlueprintToPolicy(mesh, mdMesh, policy, allocatorID);

    axom::quest::MarchingCubes mc(policy, allocatorID, DataParallelism::byPolicy);
    mc.setMesh(mesh, "mesh");
    mc.setFunctionField("fcn");

    // Each computeIsocontour() call appends to the output.
    mc.computeIsocontour(c1);
    mc.computeIsocontour(c2);
    {
      const HostContour c = hostContour(mc);
      ASSERT_EQ(c.facetCount, 2 * perContour);
      checkUnweldedConnectivity<DIM>(c);
      checkFacetsInParentCells<DIM>(c, domainsById(mdMesh));
      for(axom::IndexType i = 0; i < c.facetCount; ++i)
      {
        const double expected = i < perContour ? c1 : c2;
        for(int d = 0; d < DIM; ++d)
        {
          ASSERT_NEAR(c.coords(c.corners(i, d), 0), expected, POSITION_TOL) << "facet " << i;
        }
      }
    }

    // clearOutput() discards the output; the next call starts over.
    mc.clearOutput();
    EXPECT_EQ(mc.getContourCellCount(), 0);
    EXPECT_EQ(mc.getContourNodeCount(), 0);
    mc.computeIsocontour(c2);
    {
      const HostContour c = hostContour(mc);
      ASSERT_EQ(c.facetCount, perContour);
      checkUnweldedConnectivity<DIM>(c);
      checkNodesOnPlane<DIM>(c, mctest::PlanarField(Vec3 {1., 0., 0.}, c2));
    }

    // relinquishContourData() hands over the arrays and empties the object.
    axom::Array<axom::IndexType, 2> corners;
    axom::Array<double, 2> coords;
    axom::Array<axom::IndexType> parents, domains;
    mc.relinquishContourData(corners, coords, parents, domains);
    EXPECT_EQ(corners.shape()[0], perContour);
    EXPECT_EQ(coords.shape()[0], DIM * perContour);
    EXPECT_EQ(parents.size(), perContour);
    EXPECT_EQ(domains.size(), perContour);
    EXPECT_EQ(mc.getContourCellCount(), 0);
    EXPECT_EQ(mc.getContourNodeCount(), 0);
    EXPECT_EQ(mc.getContourFacetCorners().size(), 0);
    EXPECT_EQ(mc.getContourNodeCoords().size(), 0);

    // The object remains usable after relinquishing its output.
    mc.computeIsocontour(c1);
    {
      const HostContour c = hostContour(mc);
      ASSERT_EQ(c.facetCount, perContour);
      checkUnweldedConnectivity<DIM>(c);
      checkNodesOnPlane<DIM>(c, mctest::PlanarField(Vec3 {1., 0., 0.}, c1));
    }
  }
}

TEST(quest_marching_cubes, append_clear_relinquish_2d) { testAppendClearAndRelinquish<2>(); }
TEST(quest_marching_cubes, append_clear_relinquish_3d) { testAppendClearAndRelinquish<3>(); }

TEST(quest_marching_cubes, scan_strategies_agree_on_gyroid)
{
  constexpr int DIM = 3;
  constexpr int n = 12;
  const mctest::GyroidField f(2. * M_PI);

  conduit::Node dom, mdMesh;
  mctest::buildStructured<DIM>(dom, n, f, "fcn");
  wrapAsMultiDomain(mdMesh, dom);

  // Collect sorted parent ids from every configuration; all must agree.
  std::vector<std::vector<axom::IndexType>> parentLists;
  forEachConfiguration(mdMesh, "", 0.1, [&](const auto&, const HostContour& c) {
    EXPECT_GT(c.facetCount, 0);
    checkUnweldedConnectivity<DIM>(c);
    checkFacetsInParentCells<DIM>(c, domainsById(mdMesh));
    std::vector<axom::IndexType> parents(c.parents.begin(), c.parents.end());
    std::sort(parents.begin(), parents.end());
    parentLists.push_back(parents);
  });

  for(std::size_t i = 1; i < parentLists.size(); ++i)
  {
    EXPECT_EQ(parentLists[i], parentLists[0]) << "configuration " << i;
  }
}

template <int DIM>
void testMintOutput()
{
  constexpr int n = 6;
  const mctest::PlanarField f(mctest::PlanarField::PointType {1.0, 0.5, 0.5},
                              Vec3 {1., -2.5, DIM == 3 ? 0.5 : 0.});

  conduit::Node mdMesh;
  buildMultiDomain<DIM>(mdMesh, 2, n, f, true);

  for(auto policy : enabledPolicies())
  {
    SCOPED_TRACE(policyName(policy));
    const int allocatorID = allocatorForPolicy(policy);
    conduit::Node mesh;
    mctest::copyBlueprintToPolicy(mesh, mdMesh, policy, allocatorID);

    axom::quest::MarchingCubes mc(policy, allocatorID, DataParallelism::byPolicy);
    mc.setMesh(mesh, "mesh");
    mc.setFunctionField("fcn");
    mc.computeIsocontour(0.0);
    const HostContour c = hostContour(mc);
    ASSERT_GT(c.facetCount, 0);

    const auto cellType = DIM == 2 ? axom::mint::SEGMENT : axom::mint::TRIANGLE;
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> contour(DIM, cellType);
    mc.populateContourMesh(contour, "parentCell", "parentDomain");

    ASSERT_EQ(contour.getNumberOfCells(), c.facetCount);
    ASSERT_EQ(contour.getNumberOfNodes(), c.nodeCount);
    const auto* parentCell =
      contour.getFieldPtr<axom::IndexType>("parentCell", axom::mint::CELL_CENTERED);
    const auto* parentDomain =
      contour.getFieldPtr<axom::quest::MarchingCubes::DomainIdType>("parentDomain",
                                                                    axom::mint::CELL_CENTERED);
    for(axom::IndexType i = 0; i < c.facetCount; ++i)
    {
      EXPECT_EQ(parentCell[i], c.parents[i]);
      EXPECT_EQ(parentDomain[i], c.domains[i]);

      axom::IndexType cellNodes[3];
      contour.getCellNodeIDs(i, cellNodes);
      for(int d = 0; d < DIM; ++d)
      {
        EXPECT_EQ(cellNodes[d], c.corners(i, d));
      }
    }
    for(axom::IndexType node = 0; node < c.nodeCount; ++node)
    {
      double xyz[3];
      contour.getNode(node, xyz);
      for(int a = 0; a < DIM; ++a)
      {
        EXPECT_EQ(xyz[a], c.coords(node, a));
      }
    }
  }
}

TEST(quest_marching_cubes, mint_output_2d) { testMintOutput<2>(); }
TEST(quest_marching_cubes, mint_output_3d) { testMintOutput<3>(); }

//---------------------------------------------------------------------------
// Use before setMesh() and input validation
//---------------------------------------------------------------------------

TEST(quest_marching_cubes, calls_before_set_mesh)
{
  using axom::quest::MarchingCubes;

  // Construct in storage that holds a nonzero byte pattern, so a member that
  // the constructor leaves uninitialized reads as garbage instead of zero.
  alignas(MarchingCubes) unsigned char storage[sizeof(MarchingCubes)];
  std::fill(std::begin(storage), std::end(storage), static_cast<unsigned char>(0x5a));
  auto* mc = new(storage)
    MarchingCubes(RuntimePolicy::seq, mctest::hostAllocatorID(), DataParallelism::byPolicy);

  EXPECT_EQ(mc->getContourCellCount(), 0);
  EXPECT_EQ(mc->getContourNodeCount(), 0);
  mc->setFunctionField("fcn");
  mc->computeIsocontour(0.5);
  EXPECT_EQ(mc->getContourCellCount(), 0);
  EXPECT_EQ(mc->getContourNodeCount(), 0);
  mc->clearOutput();
  EXPECT_EQ(mc->getContourCellCount(), 0);

  mc->~MarchingCubes();
}

namespace
{
//! @brief A small 2D structured test domain with a vertex field "fcn" and a cell field "mask".
void buildValidationDomain(conduit::Node& dom)
{
  mctest::buildStructured<2>(dom, 4, mctest::PlanarField(Vec3 {1., 0., 0.}, 0.3), "fcn");
  mctest::addCellField<2>(dom, [](int, int, int) { return 1; }, "mask");
}

/*!
 * @brief Expect setMesh() (and setFunctionField(), if @a fcnField is not empty)
 *        to report a SLIC error for @a dom.
 *
 * The error must come from SLIC in every build type, not from a Conduit
 * exception or a debug-only assertion.
 */
void expectInputError(const conduit::Node& dom,
                      const std::string& maskField = {},
                      const std::string& fcnField = {})
{
  conduit::Node mdMesh;
  wrapAsMultiDomain(mdMesh, dom);
  axom::quest::MarchingCubes mc(RuntimePolicy::seq,
                                mctest::hostAllocatorID(),
                                DataParallelism::byPolicy);

  axom::slic::ScopedAbortToThrow abortGuard;
  if(fcnField.empty())
  {
    EXPECT_THROW(mc.setMesh(mdMesh, "mesh", maskField), axom::slic::SlicAbortException);
  }
  else
  {
    mc.setMesh(mdMesh, "mesh", maskField);
    EXPECT_THROW(mc.setFunctionField(fcnField), axom::slic::SlicAbortException);
  }
}
}  // end anonymous namespace

TEST(quest_marching_cubes, rejects_unstructured_topology)
{
  conduit::Node dom;
  buildValidationDomain(dom);
  dom["topologies/mesh/type"] = "unstructured";
  expectInputError(dom);
}

TEST(quest_marching_cubes, rejects_missing_coordset)
{
  conduit::Node dom;
  buildValidationDomain(dom);
  dom["topologies/mesh/coordset"] = "no_such_coordset";
  expectInputError(dom);
}

TEST(quest_marching_cubes, rejects_one_dimensional_mesh)
{
  conduit::Node dom;
  buildValidationDomain(dom);
  // Conduit derives a structured topology's dimension from its coordset.
  dom["topologies/mesh/elements/dims"].remove("j");
  dom["coordsets/coords/values"].remove("y");
  expectInputError(dom);
}

TEST(quest_marching_cubes, rejects_interleaved_coordinates)
{
  conduit::Node dom;
  buildValidationDomain(dom);

  // Rewrite the coordinates as one interleaved xyxy... buffer.
  const conduit::Node& values = dom["coordsets/coords/values"];
  const auto n = values["x"].dtype().number_of_elements();
  std::vector<double> xy(2 * n);
  for(conduit::index_t i = 0; i < n; ++i)
  {
    xy[2 * i] = values["x"].as_float64_ptr()[i];
    xy[2 * i + 1] = values["y"].as_float64_ptr()[i];
  }
  constexpr conduit::index_t stride = 2 * sizeof(double);
  dom["coordsets/coords/values/x"].set_external(conduit::DataType::float64(n, 0, stride), xy.data());
  dom["coordsets/coords/values/y"].set_external(conduit::DataType::float64(n, sizeof(double), stride),
                                                xy.data());
  ASSERT_TRUE(conduit::blueprint::mcarray::is_interleaved(dom["coordsets/coords/values"]));

  expectInputError(dom);
}

TEST(quest_marching_cubes, rejects_missing_mask_field)
{
  conduit::Node dom;
  buildValidationDomain(dom);
  expectInputError(dom, "no_such_mask");
}

TEST(quest_marching_cubes, rejects_missing_function_field)
{
  conduit::Node dom;
  buildValidationDomain(dom);
  expectInputError(dom, {}, "no_such_field");
}

TEST(quest_marching_cubes, rejects_element_function_field)
{
  conduit::Node dom;
  buildValidationDomain(dom);
  // "mask" exists but is element-associated; the contour field must be nodal.
  expectInputError(dom, {}, "mask");
}

//---------------------------------------------------------------------------
// Single-domain input
//---------------------------------------------------------------------------

namespace
{
//! @brief Expect two host contours to be identical, entry by entry.
void expectSameContour(const HostContour& a, const HostContour& b)
{
  ASSERT_EQ(a.facetCount, b.facetCount);
  ASSERT_EQ(a.nodeCount, b.nodeCount);
  EXPECT_TRUE(std::equal(a.corners.begin(), a.corners.end(), b.corners.begin()));
  EXPECT_TRUE(std::equal(a.coords.begin(), a.coords.end(), b.coords.begin()));
  EXPECT_TRUE(std::equal(a.parents.begin(), a.parents.end(), b.parents.begin()));
  EXPECT_TRUE(std::equal(a.domains.begin(), a.domains.end(), b.domains.begin()));
}

//! @brief Contour @a mesh (already in @a policy's memory) and return the host result.
HostContour contourOf(const conduit::Node& mesh, RuntimePolicy policy, double contourValue)
{
  axom::quest::MarchingCubes mc(policy, allocatorForPolicy(policy), DataParallelism::byPolicy);
  mc.setMesh(mesh, "mesh");
  mc.setFunctionField("fcn");
  mc.computeIsocontour(contourValue);
  return hostContour(mc);
}
}  // end anonymous namespace

/*!
 * @brief A single domain passed directly gives the same contour as the same
 *        domain wrapped in a one-domain multi-domain mesh.
 */
template <int DIM>
void testSingleDomainInput(bool setDomainId)
{
  constexpr int n = 8;
  const mctest::PlanarField f(mctest::PlanarField::PointType {0.4, 0.55, 0.45},
                              Vec3 {1., 2., DIM == 3 ? 3. : 0.});

  conduit::Node hostDom, hostMd;
  mctest::buildStructured<DIM>(hostDom, n, f, "fcn", mctest::SinusoidalWarp {0.03});
  if(setDomainId)
  {
    hostDom["state/domain_id"] = 7;
  }
  wrapAsMultiDomain(hostMd, hostDom);

  for(auto policy : enabledPolicies())
  {
    SCOPED_TRACE(policyName(policy));
    const int allocatorID = allocatorForPolicy(policy);
    conduit::Node dom, md;
    mctest::copyBlueprintToPolicy(dom, hostDom, policy, allocatorID);
    mctest::copyBlueprintToPolicy(md, hostMd, policy, allocatorID);

    const HostContour single = contourOf(dom, policy, 0.0);
    const HostContour wrapped = contourOf(md, policy, 0.0);
    ASSERT_GT(single.facetCount, 0);
    expectSameContour(single, wrapped);
    checkUnweldedConnectivity<DIM>(single);
    checkNodesOnPlane<DIM>(single, f);
    checkFacetsInParentCells<DIM>(single, domainsById(hostMd));
    for(axom::IndexType i = 0; i < single.facetCount; ++i)
    {
      ASSERT_EQ(single.domains[i], setDomainId ? 7 : 0);
    }
  }
}

TEST(quest_marching_cubes, single_domain_input_2d) { testSingleDomainInput<2>(false); }
TEST(quest_marching_cubes, single_domain_input_3d) { testSingleDomainInput<3>(false); }
TEST(quest_marching_cubes, single_domain_input_with_domain_id_3d)
{
  testSingleDomainInput<3>(true);
}

//! @brief One MarchingCubes object can switch between single- and multi-domain inputs.
TEST(quest_marching_cubes, switch_between_single_and_multi_domain)
{
  constexpr int DIM = 3;
  constexpr int n = 6;
  const mctest::PlanarField f(mctest::PlanarField::PointType {1.5, 0.5, 0.5}, Vec3 {1., -2.5, 0.5});

  conduit::Node multi;
  buildMultiDomain<DIM>(multi, 3, n, f, true);
  const conduit::Node& single = multi.child(1);  // state/domain_id == 11

  const HostContour multiRef = contourOf(multi, RuntimePolicy::seq, 0.0);
  const HostContour singleRef = contourOf(single, RuntimePolicy::seq, 0.0);
  ASSERT_GT(singleRef.facetCount, 0);
  ASSERT_GT(multiRef.facetCount, singleRef.facetCount);

  axom::quest::MarchingCubes mc(RuntimePolicy::seq,
                                mctest::hostAllocatorID(),
                                DataParallelism::byPolicy);
  auto run = [&](const conduit::Node& mesh) {
    mc.clearOutput();
    mc.setMesh(mesh, "mesh");
    mc.setFunctionField("fcn");
    mc.computeIsocontour(0.0);
    return hostContour(mc);
  };

  expectSameContour(run(single), singleRef);
  expectSameContour(run(multi), multiRef);
  expectSameContour(run(single), singleRef);
  const HostContour again = run(single);
  for(axom::IndexType i = 0; i < again.facetCount; ++i)
  {
    ASSERT_EQ(again.domains[i], 11);
  }
}

//! @brief A multi-domain mesh with no local domains is valid and yields an empty contour.
TEST(quest_marching_cubes, empty_multi_domain_mesh)
{
  conduit::Node empty;
  axom::quest::MarchingCubes mc(RuntimePolicy::seq,
                                mctest::hostAllocatorID(),
                                DataParallelism::byPolicy);
  mc.setMesh(empty, "mesh");
  mc.setFunctionField("fcn");
  mc.computeIsocontour(0.0);
  EXPECT_EQ(mc.getContourCellCount(), 0);
  EXPECT_EQ(mc.getContourNodeCount(), 0);
}

//! @brief A single domain that lacks the named topology is an error, not a multi-domain mesh.
TEST(quest_marching_cubes, rejects_single_domain_without_topology)
{
  conduit::Node dom;
  buildValidationDomain(dom);
  axom::quest::MarchingCubes mc(RuntimePolicy::seq,
                                mctest::hostAllocatorID(),
                                DataParallelism::byPolicy);
  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(mc.setMesh(dom, "no_such_topology"), axom::slic::SlicAbortException);
}

//---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
