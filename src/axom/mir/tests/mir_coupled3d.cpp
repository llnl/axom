// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/fmt.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>
#include <cstdlib>

namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_coupled3d");
}

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

/// Spacing on left and right of real zones.
constexpr int NLEFT = 2;
constexpr int NRIGHT = 1;

//------------------------------------------------------------------------------
/*!
 * \brief Make an explicit coordset in 3D.
 *
 * \param n_mesh The Conduit node that will contain the output mesh.
 * \param coordsetName The name of the coordset to make.
 * \param extents The mesh coordinate extents {xmin, xmax, ymin, ymax, zmin, zmax}.
 * \param dims The total number of NODES in each dimension.
 */
void make_explicit_coordset(conduit::Node &n_mesh,
                            const std::string &coordsetName,
                            double extents[6],
                            int dims[3])
{
  SLIC_ASSERT(dims[0] > 0 && dims[1] > 0 && dims[2] > 0);
  const int nnodes = dims[0] * dims[1] * dims[2];
  conduit::Node &n_coordset = n_mesh["coordsets/" + coordsetName];
  n_coordset["type"] = "explicit";
  n_coordset["values/x"].set(conduit::DataType::float64(nnodes));
  n_coordset["values/y"].set(conduit::DataType::float64(nnodes));
  n_coordset["values/z"].set(conduit::DataType::float64(nnodes));
  double *x = n_coordset["values/x"].as_double_ptr();
  double *y = n_coordset["values/y"].as_double_ptr();
  double *z = n_coordset["values/z"].as_double_ptr();
  int index = 0;
  for(int k = 0; k < dims[2]; k++)
  {
    const double tk = static_cast<double>(k) / static_cast<double>(dims[2] - 1);
    const double zc = extents[4] + tk * (extents[5] - extents[4]);
    for(int j = 0; j < dims[1]; j++)
    {
      double ty = static_cast<double>(j) / static_cast<double>(dims[1] - 1);
      const double yc = extents[2] + ty * (extents[3] - extents[2]);
      for(int i = 0; i < dims[0]; i++, index++)
      {
        double tx = static_cast<double>(i) / static_cast<double>(dims[0] - 1);
        const double xc = extents[0] + tx * (extents[1] - extents[0]);
        x[index] = xc;
        y[index] = yc;
        z[index] = zc;
      }
    }
  }
}

/*!
 * \brief Make a strided structured mesh in 3D.
 *
 * \param n_mesh The Conduit node that will contain the output mesh.
 * \param topoName The name of the topology to make.
 * \param coordsetName The name of the coordset to make.
 * \param extents The mesh coordinate extents {xmin, xmax, ymin, ymax, zmin, zmax}.
 * \param dims The total number of NODES (real+phony) in each dimension.
 *
 * \note extents and dims include phonies.
 */
void make_mesh(conduit::Node &n_mesh,
               const std::string &topoName,
               const std::string &coordsetName,
               const int dims[3])
{
  SLIC_ASSERT(dims[0] > 0 && dims[1] > 0 && dims[2] > 0);
  const int real_zone_dims[] = {dims[0] - 1 - NLEFT - NRIGHT,
                                dims[1] - 1 - NLEFT - NRIGHT,
                                dims[2] - 1 - NLEFT - NRIGHT};
  // Offsets and strides appear nodal
  const int offsets[] = {NLEFT, NLEFT, NLEFT};
  int strides[3];
  strides[0] = 1;
  strides[1] = dims[0];
  strides[2] = dims[0] * dims[1];

  // Make the strided-structured topology,
  conduit::Node &n_topo = n_mesh["topologies/" + topoName];
  n_topo["type"] = "structured";
  n_topo["coordset"] = coordsetName;
  n_topo["elements/dims/i"] = real_zone_dims[0];
  n_topo["elements/dims/j"] = real_zone_dims[1];
  n_topo["elements/dims/k"] = real_zone_dims[2];
  n_topo["elements/dims/offsets"].set(offsets, 3);
  n_topo["elements/dims/strides"].set(strides, 3);
}

/*!
 * \brief Make a matset that shapes in a CU wall and ball with other regions filled with AIR.
 *        The volume_fractions and material_ids are sized to 2x the whole mesh size, including
 *        the phony zones. The actual matset is limited to just the real zones in the
 *        strided-structured mesh.
 *
 * \param n_mesh The Conduit node that will contain the mesh.
 * \param topologyName The name of the topology associated with the matset.
 * \param matsetName The name of the matset to create.
 * \param extents The total mesh extents, including phonies.
 * \param dims The total mesh node dimensions, including phonies.
 * \param ballCenter The center of the ball.
 * \param ballRadius The radius of the ball.
 * \param wallX The X value above which zones are filled with CU.
 */
void make_matset(conduit::Node &n_mesh,
                 const std::string &topologyName,
                 const std::string &matsetName,
                 const double extents[6],
                 const int dims[3],
                 const double ballCenter[3],
                 double ballRadius,
                 const double wallX)
{
  /// Sample a zone and determine vfCU and vfAIR, returning number of materials 1 or 2.
  auto ballVF = [&](const double zExt[6], double &vfCU, double &vfAIR) -> int {
    const int nSamples = 10;
    const double br2 = ballRadius * ballRadius;
    const int nTotalSamples = nSamples * nSamples * nSamples;
    int sums[2] = {0, 0};
    for(int k = 0; k < nSamples; k++)
    {
      const double tk = static_cast<double>(k) / static_cast<double>(nSamples - 1);
      const double zc = zExt[4] + tk * (zExt[5] - zExt[4]);
      const double dz = ballCenter[2] - zc;

      for(int j = 0; j < nSamples; j++)
      {
        const double ty = static_cast<double>(j) / static_cast<double>(nSamples - 1);
        const double yc = zExt[2] + ty * (zExt[3] - zExt[2]);
        const double dy = ballCenter[1] - yc;

        for(int i = 0; i < nSamples; i++)
        {
          const double tx = static_cast<double>(i) / static_cast<double>(nSamples - 1);
          const double xc = zExt[0] + tx * (zExt[1] - zExt[0]);
          const double dx = ballCenter[0] - xc;

          const double dist2 = dx * dx + dy * dy + dz * dz;
          const int mi = dist2 <= br2 ? 1 : 0;
          sums[mi]++;
        }
      }
    }
    vfAIR = static_cast<double>(sums[0]) / static_cast<double>(nTotalSamples);
    vfCU = static_cast<double>(sums[1]) / static_cast<double>(nTotalSamples);
    const int nmats = (sums[0] == 0 || sums[1] == 0) ? 1 : 2;
    return nmats;
  };

  const int AIR = 0;
  const int CU = 1;

  // zonal dims (real + phony)
  const int zdims[] = {dims[0] - 1, dims[1] - 1, dims[2] - 1};
  const int totalZones = zdims[0] * zdims[1] * zdims[2];

  double dx = (extents[1] - extents[0]) / zdims[0];
  double dy = (extents[3] - extents[2]) / zdims[1];
  double dz = (extents[5] - extents[4]) / zdims[2];

  // Allocate these 2 arrays over all zones and we'll index into them.
  std::vector<double> volume_fractions(totalZones * 2, 0.);
  const int INVALID_MATERIAL = -1;
  std::vector<int> material_ids(totalZones * 2, INVALID_MATERIAL);
  // Defined for only the real zones.
  std::vector<int> indices;
  std::vector<int> sizes;
  std::vector<int> offsets;

  for(int k = 0; k < zdims[2]; k++)
  {
    for(int j = 0; j < zdims[1]; j++)
    {
      for(int i = 0; i < zdims[0]; i++)
      {
        const int globalZoneIndex = (k * zdims[0] * zdims[1]) + (j * zdims[0]) + i;

        // Define the material only on the real zones.
        if((i >= NLEFT && i < (zdims[0] - NRIGHT)) && (j >= NLEFT && j < (zdims[1] - NRIGHT)) &&
           (k >= NLEFT && k < (zdims[2] - NRIGHT)))
        {
          double zoneExtents[] = {extents[0] + i * dx,
                                  extents[0] + (i + 1) * dx,
                                  extents[2] + j * dy,
                                  extents[2] + (j + 1) * dy,
                                  extents[4] + k * dz,
                                  extents[4] + (k + 1) * dz};

          const double midX = (zoneExtents[0] + zoneExtents[1]) / 2.;
          double vfCU = 0., vfAIR = 1.;
          int nmats = 1;
          if(midX > wallX)
          {
            // Wall
            vfCU = 1.;
            vfAIR = 0.;
          }
          else
          {
            // Ball
            nmats = ballVF(zoneExtents, vfCU, vfAIR);
          }

          const int index = globalZoneIndex * 2;
          if(nmats == 1)
          {
            volume_fractions[index] = (vfCU > 0.) ? vfCU : vfAIR;
            material_ids[index] = (vfCU > 0.) ? CU : AIR;

            offsets.push_back(indices.size());
            indices.push_back(index);
            sizes.push_back(1);
          }
          else
          {
            volume_fractions[index] = vfAIR;
            volume_fractions[index + 1] = vfCU;

            material_ids[index] = AIR;
            material_ids[index + 1] = CU;

            offsets.push_back(indices.size());
            indices.push_back(index);
            indices.push_back(index + 1);
            sizes.push_back(2);
          }
        }
      }
    }
  }

  conduit::Node &n_matset = n_mesh["matsets/" + matsetName];
  n_matset["topology"] = topologyName;
  n_matset["material_map/AIR"] = AIR;
  n_matset["material_map/CU"] = CU;
  n_matset["volume_fractions"].set(volume_fractions);
  n_matset["material_ids"].set(material_ids);
  n_matset["indices"].set(indices);
  n_matset["sizes"].set(sizes);
  n_matset["offsets"].set(offsets);
}

/*!
 * \brief Take the real_dims and real_extents and add phonies.
 *
 * \param[in] real_dims The real node dimensions.
 * \param[in] real_extents The box that defines the extents for the real nodes.
 * \param[out] total_dims The total number of node dimensions, including phonies.
 * \param[out] total_extents The box that defines the extents for all zones, including phonies.
 */
void adjust_sizes(int real_dims[3], double real_extents[6], int total_dims[3], double total_extents[6])
{
  // Size of single zone.
  const double dx = (real_extents[1] - real_extents[0]) / (real_dims[0] - 1);
  const double dy = (real_extents[3] - real_extents[2]) / (real_dims[1] - 1);
  const double dz = (real_extents[5] - real_extents[4]) / (real_dims[2] - 1);

  total_dims[0] = real_dims[0] + NLEFT + NRIGHT;
  total_dims[1] = real_dims[1] + NLEFT + NRIGHT;
  total_dims[2] = real_dims[2] + NLEFT + NRIGHT;

  total_extents[0] = real_extents[0] - NLEFT * dx;
  total_extents[1] = real_extents[1] + NRIGHT * dx;
  total_extents[2] = real_extents[2] - NLEFT * dy;
  total_extents[3] = real_extents[3] + NRIGHT * dy;
  total_extents[4] = real_extents[4] - NLEFT * dz;
  total_extents[5] = real_extents[5] + NRIGHT * dz;
}

/*!
 * \brief Make the coarse mesh.
 *
 * \param n_mesh The Conduit node that will contain the mesh.
 * \param real_extents The box that defines the extents for the real nodes.
 * \param real_dims The number of real nodes in each dimension.
 */
void make_coarse(conduit::Node &n_mesh, double real_extents[6], int real_dims[3])
{
  SLIC_ASSERT(real_dims[0] > 0 && real_dims[1] > 0 && real_dims[2] > 0);

  int total_dims[3];
  double total_extents[6];
  adjust_sizes(real_dims, real_extents, total_dims, total_extents);

  make_explicit_coordset(n_mesh, "coarse_coords", total_extents, total_dims);
  make_mesh(n_mesh, "coarse", "coarse_coords", total_dims);

  double ballRadius = 0.3 * (real_extents[1] - real_extents[0]);
  double ballCenter[3];
  ballCenter[0] = real_extents[0] + 0.4 * (real_extents[1] - real_extents[0]);
  ballCenter[1] = real_extents[2] + 0.5 * (real_extents[3] - real_extents[2]);
  ballCenter[2] = real_extents[4] + 0.5 * (real_extents[5] - real_extents[4]);
  double wallX = real_extents[1] - 0.07 * (real_extents[1] - real_extents[0]);

  make_matset(n_mesh, "coarse", "coarse_matset", total_extents, total_dims, ballCenter, ballRadius, wallX);
}

/*!
 * \brief Make the fine mesh.
 *
 * \param n_mesh The Conduit node that will contain the mesh.
 * \param real_extents The box that defines the extents for the real nodes.
 * \param real_dims The number of real nodes in each dimension.
 * \param refinement The refinement ration in each dimension.
 */
void make_fine(conduit::Node &n_mesh, double real_extents[6], int real_dims[3], int refinement[3])
{
  SLIC_ASSERT(real_dims[0] > 0 && real_dims[1] > 0 && real_dims[2] > 0);
  SLIC_ASSERT(refinement[0] > 0 && refinement[1] > 0 && refinement[2] > 0);

  int nx = (real_dims[0] - 1) * refinement[0] + 1;
  int ny = (real_dims[1] - 1) * refinement[1] + 1;
  int nz = (real_dims[2] - 1) * refinement[2] + 1;
  int rdims[] = {nx, ny, nz};

  int total_dims[3];
  double total_extents[6];
  adjust_sizes(rdims, real_extents, total_dims, total_extents);

  make_explicit_coordset(n_mesh, "fine_coords", total_extents, total_dims);
  make_mesh(n_mesh, "fine", "fine_coords", total_dims);
}

//------------------------------------------------------------------------------
/*!
 * \brief Test coupling Elvira MIR to TopologyMapper to reconstruct material zones
 *        on a coarse mesh and then make a new material that indicates the overlap
 *        on a finer mesh.
 */
template <typename ExecSpace>
class test_coupling
{
public:
  static void test(const std::string &name)
  {
    // Make the input mesh.
    conduit::Node n_mesh;
    initialize(n_mesh);

    // host->device
    conduit::Node n_dev;
    utils::copy<ExecSpace>(n_dev, n_mesh);

    // Do MIR on the coarse mesh. The new objects will be added to n_mesh.
    mir("coarse", n_dev, "postmir", n_dev);

    // Map MIR output in n_mesh onto the fine mesh as a new matset.
    mapping(n_dev, n_dev);

    // As a check, run the generated fine matset through elvira again to make clean zones.
    mir("fine", n_dev, "check", n_dev);

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, n_dev);

    TestApp.saveVisualization(name, hostResult);

    // Remove the check and postmir meshes from the baseline to make it smaller.
    hostResult.remove("coordsets/check_coords");
    hostResult.remove("topologies/check");
    hostResult.remove("matsets/check_matset");
    hostResult.remove("coordsets/postmir_coords");
    hostResult.remove("topologies/postmir");
    hostResult.remove("matsets/postmir_matset");

    // Handle baseline comparison.
    EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostResult));
  }

private:
  /// Make the meshes
  static void initialize(conduit::Node &n_mesh)
  {
    // Unit cube with different numbers of zones (and refinements) in each dimension.
    double extents[] = {0., 1., 0., 1., 0., 1.};
    int dims[] = {11, 16, 13};
    int refinement[] = {4, 3, 2};

    make_coarse(n_mesh, extents, dims);
    make_fine(n_mesh, extents, dims, refinement);
  }

  /// Perform MIR on input mesh and make new output mesh.
  static void mir(const std::string &input_prefix,
                  conduit::Node &n_input,
                  const std::string &output_prefix,
                  conduit::Node &n_output)
  {
    SLIC_INFO(axom::fmt::format("mir {} to {}", input_prefix, output_prefix));

    // Wrap the input mesh in views.
    const conduit::Node &n_coordset = n_input[axom::fmt::format("coordsets/{}_coords", input_prefix)];
    const conduit::Node &n_topology = n_input[axom::fmt::format("topologies/{}", input_prefix)];
    const conduit::Node &n_matset = n_input[axom::fmt::format("matsets/{}_matset", input_prefix)];

    auto coordsetView = views::make_explicit_coordset<double, 3>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    // Get the indexing policy from the TopologyView type.
    auto topologyView = views::make_strided_structured_topology<3>::view(n_topology);
    using TopologyView = decltype(topologyView);
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    const int MAXMATS = 20;
    auto matsetView = views::make_unibuffer_matset<int, double, MAXMATS>::view(n_matset);
    using MatsetView = decltype(matsetView);

    // Do MIR on the mesh.
    using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView, 8>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node options;
    // Select that matset we'll operate on.
    options["matset"] = axom::fmt::format("{}_matset", input_prefix);

    // Change the names of the topology, coordset, and matset in the output.
    options["topologyName"] = output_prefix;
    options["coordsetName"] = axom::fmt::format("{}_coords", output_prefix);
    options["matsetName"] = axom::fmt::format("{}_matset", output_prefix);
    m.execute(n_input, options, n_output);
  }

  /// Map material from postmir mesh onto fine mesh to make fine matset.
  static void mapping(conduit::Node &n_src, conduit::Node &n_target)
  {
    SLIC_INFO("mapping postmir to fine");

    // Wrap the source mesh from (coarse MIR output).
    const conduit::Node &n_src_coordset = n_src["coordsets/postmir_coords"];
    const conduit::Node &n_src_topology = n_src["topologies/postmir"];
    const conduit::Node &n_src_matset = n_src["matsets/postmir_matset"];

    auto srcCoordsetView = views::make_explicit_coordset<double, 3>::view(n_src_coordset);
    using SrcCoordsetView = decltype(srcCoordsetView);

    // 3D Elvira makes polyhedral meshes
    auto srcTopologyView =
      views::make_unstructured_polyhedral_topology<axom::IndexType>::view(n_src_topology);
    using SrcTopologyView = decltype(srcTopologyView);

    constexpr int MAXMATS = 20;
    auto srcMatsetView = views::make_unibuffer_matset<int, double, MAXMATS>::view(n_src_matset);
    using SrcMatsetView = decltype(srcMatsetView);

    // Wrap the target mesh (fine)
    const conduit::Node &n_target_coordset = n_target["coordsets/fine_coords"];
    const conduit::Node &n_target_topology = n_target["topologies/fine"];

    auto targetCoordsetView = views::make_explicit_coordset<double, 3>::view(n_target_coordset);
    using TargetCoordsetView = decltype(targetCoordsetView);

    auto targetTopologyView = views::make_strided_structured_topology<3>::view(n_target_topology);
    using TargetTopologyView = decltype(targetTopologyView);

    // Make new a new matset on the target topology to record material overlaps.
    using Mapper = axom::bump::
      TopologyMapper<ExecSpace, SrcTopologyView, SrcCoordsetView, SrcMatsetView, TargetTopologyView, TargetCoordsetView>;
    Mapper mapper(srcTopologyView,
                  srcCoordsetView,
                  srcMatsetView,
                  targetTopologyView,
                  targetCoordsetView);
    conduit::Node n_opts;
    // Select the matset on the post-MIR mesh.
    n_opts["source/matsetName"] = "postmir_matset";
    // Set the name of the topology to use for the target mesh.
    n_opts["target/topologyName"] = "fine";
    // Set the name of the matset to create on the target mesh.
    n_opts["target/matsetName"] = "fine_matset";
    mapper.execute(n_src, n_opts, n_target);
  }
};

//------------------------------------------------------------------------------
TEST(mir_coupling, coupling_3d_seq)
{
  AXOM_ANNOTATE_SCOPE("coupling_3d_seq");
  test_coupling<seq_exec>::test("coupling_3d");
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_coupling, coupling_3d_omp)
{
  AXOM_ANNOTATE_SCOPE("coupling_3d_omp");
  test_coupling<omp_exec>::test("coupling_3d");
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_coupling, coupling_3d_cuda)
{
  AXOM_ANNOTATE_SCOPE("coupling_3d_cuda");
  test_coupling<cuda_exec>::test("coupling_3d");
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_coupling, coupling_3d_hip)
{
  AXOM_ANNOTATE_SCOPE("coupling_3d_hip");
  test_coupling<hip_exec>::test("coupling_3d");
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
