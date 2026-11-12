// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file quest_marching_cubes_equivalence.cpp
 *
 * @brief Compares the legacy and Bump MarchingCubes backends on structured meshes.
 *
 * Legacy output duplicates vertices by facet.
 * Bump output welds vertices and can contain polygons that the adaptor triangulates.
 * The checks compare:
 *   - E1: parent cell sets
 *   - E2: vertex sets
 *   - E3: area for planar 3D fixtures and length for unambiguous 2D fixtures
 */

#include "axom/config.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "quest_marching_cubes_equivalence.cpp requires conduit"
#endif
#ifndef AXOM_USE_BUMP
  #error "quest_marching_cubes_equivalence.cpp requires bump"
#endif

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/spin/UniformGrid.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/util/mesh_helpers.hpp"

#include "conduit_blueprint.hpp"

#include "axom/quest/tests/quest_marching_cubes_testing_helpers.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace
{
namespace mctest = axom::quest::testing::marching_cubes;

using mctest::copyBlueprintToPolicy;
using mctest::forEachEnabledPolicy;
using mctest::forEachHostPolicy;
using mctest::GyroidField;
using mctest::hostAllocatorID;
using mctest::PlanarField;
using mctest::policyName;
using mctest::RoundField;

using RuntimePolicy = axom::runtime_policy::Policy;

/*!
 * @brief Build the same box as mctest::buildStructured<3> with a uniform coordset and topology.
 *
 * @note @a n must be a power of two so the uniform and explicit coordinates agree bit for bit.
 */
template <typename Field>
void buildUniform3D(conduit::Node& mesh, int n, const Field& f, const std::string& fieldName)
{
  const int nn = n + 1;
  const conduit::index_t N = static_cast<conduit::index_t>(nn) * nn * nn;
  mesh.reset();

  conduit::Node& cs = mesh["coordsets/coords"];
  cs["type"] = "uniform";
  cs["dims/i"] = nn;
  cs["dims/j"] = nn;
  cs["dims/k"] = nn;
  cs["origin/x"] = 0.0;
  cs["origin/y"] = 0.0;
  cs["origin/z"] = 0.0;
  cs["spacing/dx"] = 1.0 / n;
  cs["spacing/dy"] = 1.0 / n;
  cs["spacing/dz"] = 1.0 / n;

  conduit::Node& topo = mesh["topologies/mesh"];
  topo["type"] = "uniform";
  topo["coordset"] = "coords";

  conduit::Node& fld = mesh["fields/" + fieldName];
  fld["topology"] = "mesh";
  fld["association"] = "vertex";
  fld["values"].set(conduit::DataType::float64(N));
  auto* fv = fld["values"].as_float64_ptr();
  conduit::index_t idx = 0;
  for(int k = 0; k < nn; ++k)
  {
    for(int j = 0; j < nn; ++j)
    {
      for(int i = 0; i < nn; ++i, ++idx)
      {
        fv[idx] = f(double(i) / n, double(j) / n, double(k) / n);
      }
    }
  }
}

//---------------------------------------------------------------------------
// Comparable backend results
//---------------------------------------------------------------------------

struct BackendResult
{
  std::set<axom::IndexType> crossingCells;               //!< E1
  std::vector<axom::primal::Point<double, 3>> vertices;  //!< E2 (z==0 in 2D)
  double measure {0.0};                                  //!< E3: area in 3D, length in 2D
  axom::IndexType facetCount {0};
  axom::IndexType nodeCount {0};
};

template <int DIM>
BackendResult extractBackendResult(const axom::quest::MarchingCubes& mc,
                                   axom::IndexType facet_begin = 0,
                                   axom::IndexType node_begin = 0)
{
  BackendResult r;
  const axom::IndexType facet_end = mc.getContourCellCount();
  const axom::IndexType node_end = mc.getContourNodeCount();
  r.facetCount = facet_end - facet_begin;
  r.nodeCount = node_end - node_begin;

  const axom::Array<double, 2> coordsHost(mc.getContourNodeCoords(), hostAllocatorID());
  const axom::Array<axom::IndexType, 2> cornersHost(mc.getContourFacetCorners(), hostAllocatorID());
  const axom::Array<axom::IndexType> parentsHost(mc.getContourFacetParents(), hostAllocatorID());
  const auto coords = coordsHost.view();
  const auto corners = cornersHost.view();
  const auto parents = parentsHost.view();

  for(axom::IndexType f = facet_begin; f < facet_end; ++f)
  {
    r.crossingCells.insert(parents[f]);
  }

  for(axom::IndexType v = node_begin; v < node_end; ++v)
  {
    axom::primal::Point<double, 3> p {};
    p[0] = coords(v, 0);
    p[1] = coords(v, 1);
    p[2] = (DIM == 3) ? coords(v, 2) : 0.0;
    r.vertices.push_back(p);
  }

  for(axom::IndexType f = facet_begin; f < facet_end; ++f)
  {
    if(DIM == 3)
    {
      axom::primal::Point<double, 3> a {}, b {}, c {};
      for(int d = 0; d < 3; ++d)
      {
        a[d] = coords(corners(f, 0), d);
        b[d] = coords(corners(f, 1), d);
        c[d] = coords(corners(f, 2), d);
      }
      const auto u = axom::primal::Vector<double, 3>(a, b);
      const auto w = axom::primal::Vector<double, 3>(a, c);
      r.measure += 0.5 * axom::primal::Vector<double, 3>::cross_product(u, w).norm();
    }
    else
    {
      axom::primal::Point<double, 2> a {}, b {};
      for(int d = 0; d < 2; ++d)
      {
        a[d] = coords(corners(f, 0), d);
        b[d] = coords(corners(f, 1), d);
      }
      r.measure += axom::primal::Vector<double, 2>(a, b).norm();
    }
  }

  return r;
}

template <int DIM>
BackendResult runBackend(const conduit::Node& mesh,
                         const std::string& field_name,
                         double contour_value,
                         RuntimePolicy policy,
                         bool use_bump)
{
  namespace quest = axom::quest;

  const int allocator_id = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocator_id, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setBackend(use_bump ? quest::MarchingCubesBackend::bump : quest::MarchingCubesBackend::legacy);

  conduit::Node exec_mesh;
  copyBlueprintToPolicy(exec_mesh, mesh, policy, allocator_id);
  mc.setMesh(exec_mesh, "mesh");
  mc.setFunctionField(field_name);
  mc.computeIsocontour(contour_value);

  return extractBackendResult<DIM>(mc);
}

template <int DIM>
std::array<BackendResult, 3> runAccumulatedBackend(const conduit::Node& mesh,
                                                   const std::array<std::string, 3>& field_names,
                                                   RuntimePolicy policy,
                                                   bool use_bump)
{
  namespace quest = axom::quest;

  const int allocator_id = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocator_id, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setBackend(use_bump ? quest::MarchingCubesBackend::bump : quest::MarchingCubesBackend::legacy);

  conduit::Node exec_mesh;
  copyBlueprintToPolicy(exec_mesh, mesh, policy, allocator_id);
  mc.setMesh(exec_mesh, "mesh");

  std::array<BackendResult, 3> results;
  axom::IndexType facet_begin = 0;
  axom::IndexType node_begin = 0;
  for(int field = 0; field < 3; ++field)
  {
    mc.setFunctionField(field_names[field]);
    mc.computeIsocontour(0.0);
    results[field] = extractBackendResult<DIM>(mc, facet_begin, node_begin);
    facet_begin = mc.getContourFacetCount();
    node_begin = mc.getContourNodeCount();
  }
  return results;
}

//---------------------------------------------------------------------------
// Check each point set against the other within a fixed tolerance.
//---------------------------------------------------------------------------

//! @brief Return the largest nearest-candidate distance and its source point index.
double oneSidedHausdorff(const std::vector<axom::primal::Point<double, 3>>& from,
                         const std::vector<axom::primal::Point<double, 3>>& to,
                         double tolerance,
                         std::size_t& argMax)
{
  argMax = 0;
  if(from.empty() || to.empty())
  {
    return from.empty() ? 0.0 : std::numeric_limits<double>::infinity();
  }

  using Grid = axom::spin::UniformGrid<std::size_t, 3>;
  Grid::BoxType bounds(to.data(), static_cast<int>(to.size()));
  bounds.expand(tolerance);
  const int resolution[3] {16, 16, 16};
  Grid grid(bounds, resolution, hostAllocatorID());
  for(std::size_t vertex = 0; vertex < to.size(); ++vertex)
  {
    grid.insert(Grid::BoxType(to[vertex]), vertex);
  }

  const auto grid_query = grid.getQueryObject();
  double worst = 0.0;
  for(std::size_t vertex = 0; vertex < from.size(); ++vertex)
  {
    Grid::BoxType query_bounds(from[vertex]);
    query_bounds.expand(tolerance);
    double nearest = std::numeric_limits<double>::infinity();
    grid_query.visitCandidates(query_bounds, [&](std::size_t candidate) {
      const double distance = axom::primal::Vector<double, 3>(from[vertex], to[candidate]).norm();
      nearest = std::min(nearest, distance);
    });
    if(nearest > worst)
    {
      worst = nearest;
      argMax = vertex;
    }
  }
  return worst;
}

void compareVertexSets(const BackendResult& reference, const BackendResult& result, double tolerance)
{
  std::size_t worst_vertex = 0;
  EXPECT_LT(oneSidedHausdorff(result.vertices, reference.vertices, tolerance, worst_vertex), tolerance)
    << "result vertex " << worst_vertex << " has no reference counterpart within tolerance";
  EXPECT_LT(oneSidedHausdorff(reference.vertices, result.vertices, tolerance, worst_vertex), tolerance)
    << "reference vertex " << worst_vertex << " has no result counterpart within tolerance";
}

//---------------------------------------------------------------------------
// Backend comparison
//---------------------------------------------------------------------------

void compareBackendResults(const BackendResult& legacy,
                           const BackendResult& bump,
                           const std::string& label,
                           bool compare_measure,
                           double vertex_tolerance = 1.0e-5,
                           double measure_relative_tolerance = 1.0e-5)
{
  SCOPED_TRACE(label);
  SLIC_INFO(axom::fmt::format(
    "[{}] legacy: {} facets / {} nodes / {} cells; bump: {} facets / {} nodes / {} cells",
    label,
    legacy.facetCount,
    legacy.nodeCount,
    legacy.crossingCells.size(),
    bump.facetCount,
    bump.nodeCount,
    bump.crossingCells.size()));

  ASSERT_GT(legacy.facetCount, 0) << "[" << label << "] legacy produced an empty contour";
  ASSERT_GT(bump.facetCount, 0) << "[" << label << "] bump produced an empty contour";
  EXPECT_GT(legacy.measure, 0.0);
  EXPECT_GT(bump.measure, 0.0);

  EXPECT_EQ(legacy.crossingCells, bump.crossingCells)
    << "E1 [" << label << "]: backends report different parent cells.";

  compareVertexSets(legacy, bump, vertex_tolerance);

  if(compare_measure)
  {
    const double relative_difference =
      std::abs(bump.measure - legacy.measure) / std::max(legacy.measure, 1.0e-300);
    EXPECT_LT(relative_difference, measure_relative_tolerance)
      << "E3 [" << label << "]: comparable contours have different measures";
  }
}

template <int DIM, typename Field>
void compareBackends(int n,
                     const Field& f,
                     double contourVal,
                     RuntimePolicy policy,
                     const std::string& label,
                     bool compare_measure = true,
                     double vertexTol = 1.0e-5,
                     // Bump interpolates edge crossings in float.
                     // Passing cases have measured relative errors from 5e-9 to 3e-6.
                     double measureRelTol = 1.0e-5)
{
  const std::string fieldName = "fcn";
  conduit::Node mesh;
  mctest::buildStructured<DIM>(mesh, n, f, fieldName);

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();

  const auto legacy = runBackend<DIM>(mesh, fieldName, contourVal, policy, /*useBump=*/false);
  const auto bump = runBackend<DIM>(mesh, fieldName, contourVal, policy, /*useBump=*/true);

  compareBackendResults(legacy, bump, label, compare_measure, vertexTol, measureRelTol);
}

/*!
 * @brief An isovalue outside the data range is valid input, not an error.
 *
 * This covers Blueprint population and relinquishment after an extraction with no crossing cells.
 */
void test_empty_contour(RuntimePolicy policy)
{
  namespace quest = axom::quest;
  conduit::Node mesh;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  mctest::buildStructured<3>(mesh, 6, f, "fcn");

  const int allocatorID = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocatorID, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setBackend(quest::MarchingCubesBackend::bump);
  mc.setMesh(mesh, "mesh");
  mc.setFunctionField("fcn");
  mc.computeIsocontour(0.0);
  ASSERT_GT(mc.getContourCellCount(), 0);
  mc.clearOutput();
  mc.computeIsocontour(1000.0);  // far outside the range of the signed distance

  EXPECT_EQ(mc.getContourCellCount(), 0);
  EXPECT_EQ(mc.getContourNodeCount(), 0);

  // Neither of these may throw or abort.
  conduit::Node bp;
  mc.populateContourMeshBlueprint(bp);
  ASSERT_EQ(bp.number_of_children(), 1);
  EXPECT_FALSE(bp.child(0).has_path("topologies"));
  conduit::Node bpTri;
  mc.populateContourMeshBlueprint(bpTri, /*triangulate=*/true);

  conduit::Node relinquished;
  mc.relinquishContourDataBlueprint(relinquished);
  SUCCEED();
}

/*!
 * @brief Compare uniform and rectilinear input with equivalent explicit input.
 *
 * The legacy backend cannot read uniform or rectilinear meshes, so the explicit
 * Bump result is the reference for all three representations of the same box.
 */
void test_uniform_and_rectilinear(RuntimePolicy policy)
{
  const int n = 8;  // power of two: see buildUniform3D's note on coordinate agreement
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  const std::string fieldName = "fcn";

  conduit::Node structured, uniform, rectilinear;
  mctest::buildStructured<3>(structured, n, f, fieldName);
  buildUniform3D(uniform, n, f, fieldName);
  rectilinear.set(uniform);
  conduit::blueprint::mesh::coordset::uniform::to_rectilinear(uniform["coordsets/coords"],
                                                              rectilinear["coordsets/coords"]);
  rectilinear["topologies/mesh/type"] = "rectilinear";

  for(const auto& m : {&uniform, &rectilinear})
  {
    conduit::Node info;
    ASSERT_TRUE(conduit::blueprint::mesh::verify(*m, info)) << info.to_yaml();
  }

  const auto refRun = runBackend<3>(structured, fieldName, 0.0, policy, /*useBump=*/true);
  ASSERT_GT(refRun.facetCount, 0);

  const auto uniformRun = runBackend<3>(uniform, fieldName, 0.0, policy, /*useBump=*/true);
  const auto rectRun = runBackend<3>(rectilinear, fieldName, 0.0, policy, /*useBump=*/true);

  SLIC_INFO(
    axom::fmt::format("[uniform/rectilinear] structured-explicit: {} facets; uniform: {} facets; "
                      "rectilinear: {} facets",
                      refRun.facetCount,
                      uniformRun.facetCount,
                      rectRun.facetCount));

  for(const auto& kv :
      {std::make_pair("uniform", &uniformRun), std::make_pair("rectilinear", &rectRun)})
  {
    const std::string what = kv.first;
    SCOPED_TRACE(what);
    const BackendResult& r = *kv.second;
    ASSERT_GT(r.facetCount, 0);
    EXPECT_EQ(r.crossingCells, refRun.crossingCells)
      << what << " cut a different cell set than the equivalent structured-explicit mesh";
    compareVertexSets(refRun, r, 1.0e-8);
    EXPECT_NEAR(r.measure, refRun.measure, 1.0e-12 * std::max(refRun.measure, 1.0))
      << what << " produced a different contour measure";
  }
}

/*!
 * @brief Check that MarchingCubes rejects a float32 function field.
 *
 * Both backends require float64 function fields.
 */
void test_float32_field_rejected(RuntimePolicy policy)
{
  namespace quest = axom::quest;
  const int n = 6;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  conduit::Node mesh;
  mctest::buildStructured<3>(mesh, n, f, "fcn");

  // Rewrite the function field as float32 while preserving its values.
  {
    const conduit::Node& n_old = mesh.fetch_existing("fields/fcn/values");
    const auto acc = n_old.as_double_accessor();
    const conduit::index_t N = n_old.dtype().number_of_elements();
    std::vector<float> tmp(static_cast<std::size_t>(N));
    for(conduit::index_t i = 0; i < N; ++i)
    {
      tmp[static_cast<std::size_t>(i)] = static_cast<float>(acc[i]);
    }
    mesh["fields/fcn/values"].set(tmp.data(), static_cast<conduit::index_t>(tmp.size()));
  }
  ASSERT_TRUE(mesh.fetch_existing("fields/fcn/values").dtype().is_float32());

  const int allocatorID = axom::policyToDefaultAllocatorID(policy);
  quest::MarchingCubes mc(policy, allocatorID, quest::MarchingCubesDataParallelism::byPolicy);
  mc.setBackend(quest::MarchingCubesBackend::bump);

  axom::slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(
    {
      mc.setMesh(mesh, "mesh");
      mc.setFunctionField("fcn");
      mc.computeIsocontour(0.0);
    },
    axom::slic::SlicAbortException);
}

//! @brief Run the Bump backend and require a field-layout error.
void expectBumpFieldLayoutRejected(const conduit::Node& mesh, RuntimePolicy policy)
{
  axom::slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(
    {
      const int allocatorID = axom::policyToDefaultAllocatorID(policy);
      axom::quest::MarchingCubes mc(policy,
                                    allocatorID,
                                    axom::quest::MarchingCubesDataParallelism::byPolicy);
      mc.setBackend(axom::quest::MarchingCubesBackend::bump);
      mc.setMesh(mesh, "mesh");
      mc.setFunctionField("fcn");
      mc.computeIsocontour(0.0);
    },
    axom::slic::SlicAbortException);
}

/*!
 * @brief Check that MarchingCubes rejects unsupported strided field layouts.
 *
 * Bump's flat field view cannot represent a permuted field or one whose offsets
 * and strides differ from the topology.
 */
void test_invalid_field_layouts_rejected(RuntimePolicy policy)
{
  constexpr int n = 6;
  constexpr int pad = 2;
  constexpr int nnPad = n + 1 + 2 * pad;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};

  conduit::Node permuted;
  mctest::buildStridedStructured<3>(permuted, n, pad, f, "fcn");
  permuted["fields/fcn/strides"].set(std::vector<conduit::int32> {nnPad * nnPad, nnPad, 1});
  expectBumpFieldLayoutRejected(permuted, policy);

  conduit::Node mismatched;
  mctest::buildStridedStructured<3>(mismatched, n, pad, f, "fcn");
  mismatched["fields/fcn/offsets"].set(std::vector<conduit::int32> {pad + 1, pad, pad});
  expectBumpFieldLayoutRejected(mismatched, policy);
}

/*!
 * @brief Preserve parent ids when the input has an "originalElements" field.
 *
 * TableBasedExtractor otherwise propagates the input field instead of writing
 * source zone indices. MarchingCubes uses a private field name to avoid this collision.
 */
void test_original_elements_collision(RuntimePolicy policy)
{
  const int n = 8;
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  const std::string fieldName = "fcn";

  conduit::Node clean, poisoned;
  mctest::buildStructured<3>(clean, n, f, fieldName);
  poisoned.set(clean);

  // Use negative decoy values that cannot be valid zone indices.
  for(const char* field_name :
      {"originalElements", "__axom_mc_originalElements", "__axom_mc_originalElements_"})
  {
    const conduit::index_t nCells = static_cast<conduit::index_t>(n) * n * n;
    conduit::Node& fld = poisoned[std::string("fields/") + field_name];
    fld["topology"] = "mesh";
    fld["association"] = "element";
    fld["values"].set(conduit::DataType::int64(nCells));
    auto* v = fld["values"].as_int64_ptr();
    for(conduit::index_t i = 0; i < nCells; ++i)
    {
      v[i] = -7;  // if these leak through as parent ids, the check below fails
    }
  }
  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(poisoned, info)) << info.to_yaml();

  const auto cleanRun = runBackend<3>(clean, fieldName, 0.0, policy, /*useBump=*/true);
  const auto poisonedRun = runBackend<3>(poisoned, fieldName, 0.0, policy, /*useBump=*/true);

  ASSERT_GT(cleanRun.facetCount, 0);
  EXPECT_EQ(poisonedRun.crossingCells, cleanRun.crossingCells)
    << "a pre-existing 'originalElements' field on the input changed the reported parent ids";
  for(const auto id : poisonedRun.crossingCells)
  {
    EXPECT_GE(id, 0) << "decoy 'originalElements' values leaked through as parent cell ids";
  }
}

/*!
 * @brief Compare strided and compact representations of the same mesh.
 *
 * Ghost values continue the field outside the real zone range,
 * so ignoring the topology offsets produces extra facets.
 */
void test_strided_structured(RuntimePolicy policy)
{
  const int n = 8;
  const int pad = 2;  // ghost layers
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  const std::string fieldName = "fcn";

  conduit::Node compact, strided;
  mctest::buildStructured<3>(compact, n, f, fieldName);
  mctest::buildStridedStructured<3>(strided, n, pad, f, fieldName);

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(strided, info)) << info.to_yaml();

  // Legacy results establish that the compact and strided fixtures agree.
  // Parent-id bounds catch ghost cells even when facet counts happen to match.
  const auto legacyCompact = runBackend<3>(compact, fieldName, 0.0, policy, /*useBump=*/false);
  const auto legacyStrided = runBackend<3>(strided, fieldName, 0.0, policy, /*useBump=*/false);
  const auto bumpStrided = runBackend<3>(strided, fieldName, 0.0, policy, /*useBump=*/true);

  SLIC_INFO(axom::fmt::format(
    "[strided] legacy/compact={} facets, legacy/strided={} facets, bump/strided={} facets",
    legacyCompact.facetCount,
    legacyStrided.facetCount,
    bumpStrided.facetCount));

  ASSERT_GT(legacyCompact.facetCount, 0);
  ASSERT_EQ(legacyStrided.crossingCells, legacyCompact.crossingCells)
    << "legacy results differ between strided and compact input";

  EXPECT_EQ(bumpStrided.crossingCells, legacyCompact.crossingCells)
    << "Bump strided input cut a different cell set than compact input";

  const axom::IndexType nCells = static_cast<axom::IndexType>(n) * n * n;
  for(const auto id : bumpStrided.crossingCells)
  {
    EXPECT_GE(id, 0) << "parent id below the real zone range (ghost leak)";
    EXPECT_LT(id, nCells) << "parent id above the real zone range (ghost leak)";
  }
}

//---------------------------------------------------------------------------
// Test bodies
//---------------------------------------------------------------------------

void test_planar_3d(RuntimePolicy policy)
{
  // An axis-aligned plane has no ambiguous cells, so compare its area.
  PlanarField f {{0.0, 0.0, 1.0}, 0.5};
  compareBackends<3>(8, f, 0.0, policy, "planar3d");
}

void test_sidre_native_layout(RuntimePolicy policy)
{
  constexpr int cells_per_axis = 4;
  constexpr int cell_count = cells_per_axis * cells_per_axis * cells_per_axis;
  constexpr double plane_height = 0.375;
  constexpr conduit::int32 mask_value = 3;
  const PlanarField field {{0.0, 0.0, 1.0}, plane_height};
  conduit::Node mesh;
  mctest::buildStructured<3>(mesh, cells_per_axis, field, "fcn");

  conduit::Node& mask = mesh["fields/mask"];
  mask["association"] = "element";
  mask["topology"] = "mesh";
  mask["values"].set(conduit::DataType::int32(cell_count));
  auto* mask_values = mask["values"].as_int32_ptr();
  std::set<axom::IndexType> expected_parents;
  for(int cell = 0; cell < cell_count; ++cell)
  {
    const bool selected = cell % cells_per_axis < cells_per_axis / 2;
    mask_values[cell] = selected ? mask_value : 0;
    if(selected && cell / (cells_per_axis * cells_per_axis) == 1)
    {
      expected_parents.insert(cell);
    }
  }

  axom::sidre::DataStore datastore;
  auto* mesh_group = datastore.getRoot()->createGroup("mesh");
  mesh_group->setDefaultAllocator(hostAllocatorID());
  ASSERT_TRUE(mesh_group->importConduitTree(mesh));
  conduit::Node native_layout;
  ASSERT_TRUE(mesh_group->createNativeLayout(native_layout));

  for(const auto* path : {"coordsets/coords/values/x",
                          "coordsets/coords/values/y",
                          "coordsets/coords/values/z",
                          "fields/fcn/values",
                          "fields/mask/values"})
  {
    SCOPED_TRACE(path);
    const auto* view = mesh_group->getView(path);
    ASSERT_NE(view, nullptr);
    EXPECT_EQ(native_layout[path].data_ptr(), view->getVoidPtr());
    EXPECT_NE(native_layout[path].data_ptr(), mesh[path].data_ptr());
  }

  for(const bool use_bump : {false, true})
  {
    SCOPED_TRACE(use_bump ? "bump" : "legacy");
    std::array<BackendResult, 2> results;
    int result_index = 0;
    for(const auto* input : {&mesh, &native_layout})
    {
      SCOPED_TRACE(result_index == 0 ? "conduit" : "sidre");
      axom::quest::MarchingCubes mc(policy,
                                    hostAllocatorID(),
                                    axom::quest::MarchingCubesDataParallelism::byPolicy);
      mc.setBackend(use_bump ? axom::quest::MarchingCubesBackend::bump
                             : axom::quest::MarchingCubesBackend::legacy);
      mc.setMesh(*input, "mesh", "mask");
      mc.setMaskValue(mask_value);
      mc.setFunctionField("fcn");
      mc.computeIsocontour(0.0);
      auto& result = results[result_index++];
      result = extractBackendResult<3>(mc);
      ASSERT_GT(result.facetCount, 0);
      EXPECT_EQ(result.crossingCells, expected_parents);
      EXPECT_NEAR(result.measure, 0.5, 1.0e-12);
      for(const auto& vertex : result.vertices)
      {
        EXPECT_NEAR(vertex[2], plane_height, 1.0e-12);
      }
    }

    const auto& reference = results[0];
    const auto& sidre = results[1];
    EXPECT_EQ(sidre.facetCount, reference.facetCount);
    EXPECT_EQ(sidre.nodeCount, reference.nodeCount);
    compareVertexSets(reference, sidre, 1.0e-8);
  }
}

void test_oblique_planar_3d(RuntimePolicy policy)
{
  // An oblique plane exercises more case-table entries without ambiguity.
  const double s = 1.0 / std::sqrt(1.0 + 0.16 + 1.44);
  PlanarField f {{1.0 * s, 0.4 * s, 1.2 * s}, 1.3 * s};
  compareBackends<3>(8, f, 0.0, policy, "oblique_planar3d");
}

void test_round_3d(RuntimePolicy policy)
{
  RoundField f {{0.5, 0.5, 0.5}, 0.25};
  compareBackends<3>(12, f, 0.0, policy, "round3d", false);
}

void test_gyroid_3d(RuntimePolicy policy)
{
  GyroidField f {3.0 * M_PI};
  compareBackends<3>(10, f, 0.0, policy, "gyroid3d", false);
}

void test_planar_2d(RuntimePolicy policy)
{
  PlanarField f {{0.0, 1.0, 0.0}, 0.5};
  compareBackends<2>(8, f, 0.0, policy, "planar2d");
}

void test_round_2d(RuntimePolicy policy)
{
  RoundField f {{0.5, 0.5, 0.0}, 0.25};
  compareBackends<2>(12, f, 0.0, policy, "round2d");
}

template <int DIM>
void test_accumulated_fields(RuntimePolicy policy)
{
  constexpr int n = 12;
  PlanarField plane {{0.47, 0.43, 0.39}, {1.0, 0.4, 1.2}};
  RoundField round {{0.5, 0.5, DIM == 3 ? 0.5 : 0.0}, 0.27};
  GyroidField gyroid {{3.0, 3.0, DIM == 3 ? 1.5 : 0.0}};
  const std::array<std::string, 3> field_names {{"plane", "round", "gyroid"}};

  conduit::Node mesh;
  mctest::buildStructured<DIM>(mesh, n, plane, field_names[0]);
  mctest::addVertexField<DIM>(mesh, round, field_names[1]);
  mctest::addVertexField<DIM>(mesh, gyroid, field_names[2]);

  const auto legacy = runAccumulatedBackend<DIM>(mesh, field_names, policy, false);
  const auto bump = runAccumulatedBackend<DIM>(mesh, field_names, policy, true);
  const std::array<bool, 3> compare_measure {{true, DIM == 2, false}};

  for(int field = 0; field < 3; ++field)
  {
    compareBackendResults(legacy[field],
                          bump[field],
                          "accumulated_" + field_names[field],
                          compare_measure[field]);
  }
}

/*!
 * @brief Compare float and double classification near the isovalue.
 *
 * Values just above the isovalue in double can round to the isovalue in float.
 * Structured and unstructured forms of the same mesh must still select the
 * same crossing cells.
 */
void test_float_ulp_band_falsification(RuntimePolicy policy)
{
  const int n = 4;
  const std::string fieldName = "fcn";
  const double contourVal = 1.0;

  conduit::Node mesh;
  PlanarField f {{0.0, 0.0, 1.0}, 0.5};
  mctest::buildStructured<3>(mesh, n, f, fieldName);

  // Put one cell's corners just above the isovalue in double but equal to it in float.
  auto* fv = mesh["fields/" + fieldName + "/values"].as_float64_ptr();
  const conduit::index_t N = mesh["fields/" + fieldName + "/values"].dtype().number_of_elements();
  const int nn = n + 1;
  auto nodeAt = [&](int i, int j, int k) { return i + j * nn + k * nn * nn; };
  for(conduit::index_t i = 0; i < N; ++i)
  {
    fv[i] = contourVal + 1.0;
  }
  // Add a separate contour so equal empty results cannot pass the test.
  for(int k = 3; k <= n; ++k)
  {
    for(int j = 0; j < nn; ++j)
    {
      for(int i = 0; i < nn; ++i)
      {
        fv[nodeAt(i, j, k)] = contourVal - 1.0;
      }
    }
  }
  // Four corners of cell (0,0,0) lie in the double-to-float rounding gap.
  const double tiny = std::nextafter(contourVal, 2.0) - contourVal;  // one double ULP
  fv[nodeAt(0, 0, 0)] = contourVal + tiny;
  fv[nodeAt(1, 0, 0)] = contourVal + tiny;
  fv[nodeAt(0, 1, 0)] = contourVal + tiny;
  fv[nodeAt(1, 1, 0)] = contourVal + tiny;

  // Confirm the rounding behavior required by the test.
  ASSERT_GT(fv[nodeAt(0, 0, 0)], contourVal);
  ASSERT_FALSE(static_cast<float>(fv[nodeAt(0, 0, 0)]) > static_cast<float>(contourVal))
    << "perturbed value did not round to float(contourVal)";

  // Reuse the coordinates and field with an unstructured hex topology.
  conduit::Node unstructuredMesh;
  unstructuredMesh.set(mesh);
  axom::quest::util::convert_blueprint_structured_explicit_to_unstructured(unstructuredMesh,
                                                                           "mesh",
                                                                           RuntimePolicy::seq);
  conduit::Node uinfo;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(unstructuredMesh, uinfo)) << uinfo.to_yaml();

  const auto structuredRun = runBackend<3>(mesh, fieldName, contourVal, policy, /*useBump=*/true);
  const auto unstructuredRun =
    runBackend<3>(unstructuredMesh, fieldName, contourVal, policy, /*useBump=*/true);

  SLIC_INFO(axom::fmt::format(
    "[ulp_band] bump/structured: {} facets; bump/unstructured (same geometry): {} facets",
    structuredRun.facetCount,
    unstructuredRun.facetCount));

  ASSERT_GT(structuredRun.facetCount, 0) << "test mesh must contain a contour";
  ASSERT_GT(unstructuredRun.facetCount, 0);
  std::set<axom::IndexType> expected_parents;
  for(int cell = 2 * n * n; cell < 3 * n * n; ++cell)
  {
    expected_parents.insert(cell);
  }
  EXPECT_EQ(structuredRun.crossingCells, expected_parents);
  EXPECT_EQ(structuredRun.crossingCells, unstructuredRun.crossingCells)
    << "the two paths cut different cell sets on identical geometry";
  compareVertexSets(structuredRun, unstructuredRun, 1.0e-5);
}

}  // namespace

TEST(quest_marching_cubes_equivalence, planar_3d) { forEachHostPolicy(test_planar_3d); }
TEST(quest_marching_cubes_equivalence, sidre_native_layout)
{
  forEachHostPolicy(test_sidre_native_layout);
}
TEST(quest_marching_cubes_equivalence, oblique_planar_3d_seq)
{
  test_oblique_planar_3d(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, round_3d) { forEachEnabledPolicy(test_round_3d); }
TEST(quest_marching_cubes_equivalence, gyroid_3d) { forEachEnabledPolicy(test_gyroid_3d); }
TEST(quest_marching_cubes_equivalence, planar_2d_seq) { test_planar_2d(RuntimePolicy::seq); }
TEST(quest_marching_cubes_equivalence, round_2d) { forEachHostPolicy(test_round_2d); }
TEST(quest_marching_cubes_equivalence, accumulated_fields_2d)
{
  forEachEnabledPolicy([](RuntimePolicy policy) { test_accumulated_fields<2>(policy); });
}
TEST(quest_marching_cubes_equivalence, accumulated_fields_3d)
{
  forEachEnabledPolicy([](RuntimePolicy policy) { test_accumulated_fields<3>(policy); });
}
TEST(quest_marching_cubes_equivalence, uniform_and_rectilinear_seq)
{
  test_uniform_and_rectilinear(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, strided_structured)
{
  forEachEnabledPolicy(test_strided_structured);
}
TEST(quest_marching_cubes_equivalence, float32_field_rejected_seq)
{
  test_float32_field_rejected(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, invalid_field_layouts_rejected_seq)
{
  test_invalid_field_layouts_rejected(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, original_elements_collision_seq)
{
  test_original_elements_collision(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, empty_contour_seq)
{
  test_empty_contour(RuntimePolicy::seq);
}
TEST(quest_marching_cubes_equivalence, float_ulp_band_falsification_seq)
{
  test_float_ulp_band_falsification(RuntimePolicy::seq);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
