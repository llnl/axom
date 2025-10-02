// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EXAMPLES_TUTORIAL_SIMPLE_RUNMIR_HPP
#define AXOM_MIR_EXAMPLES_TUTORIAL_SIMPLE_RUNMIR_HPP
#include "axom/config.hpp"
#include "axom/core.hpp"  // for axom macros
#include "axom/slic.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"  // for Mir classes & functions

#include <conduit/conduit_node.hpp>

//--------------------------------------------------------------------------------
/*!
 * \brief Run MIR on the tri input mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 *
 * \param hostMesh A conduit node that contains the test mesh.
 * \param options A conduit node that contains the test mesh.
 * \param hostResult A conduit node that will contain the MIR results.
 */
template <typename ExecSpace>
int runMIR_tri(const conduit::Node &hostMesh, const conduit::Node &options, conduit::Node &hostResult)
{
  AXOM_ANNOTATE_SCOPE("runMIR_tri");
  namespace utils = axom::bump::utilities;
  namespace views = axom::bump::views;
  std::string shape = hostMesh["topologies/mesh/elements/shape"].as_string();
  SLIC_INFO(axom::fmt::format("Using policy {}", axom::execution_space<ExecSpace>::name()));

  // host->device
  conduit::Node deviceMesh;
  {
    AXOM_ANNOTATE_SCOPE("host->device");
    utils::copy<ExecSpace>(deviceMesh, hostMesh);
  }

  conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
  conduit::Node &n_topo = deviceMesh["topologies/mesh"];
  conduit::Node &n_matset = deviceMesh["matsets/mat"];
  auto connView = utils::make_array_view<int>(n_topo["elements/connectivity"]);

  // Make matset view. (There's often 1 more material so add 1)
  constexpr int MAXMATERIALS = 12;
  using MatsetView = views::UnibufferMaterialView<int, float, MAXMATERIALS + 1>;
  MatsetView matsetView;
  matsetView.set(utils::make_array_view<int>(n_matset["material_ids"]),
                 utils::make_array_view<float>(n_matset["volume_fractions"]),
                 utils::make_array_view<int>(n_matset["sizes"]),
                 utils::make_array_view<int>(n_matset["offsets"]),
                 utils::make_array_view<int>(n_matset["indices"]));

  // Make Coord/Topo views.
  conduit::Node deviceResult;
  auto coordsetView = views::make_explicit_coordset<float, 2>::view(n_coordset);
  using CoordsetView = decltype(coordsetView);
  using TopologyView = views::UnstructuredTopologySingleShapeView<views::TriShape<int>>;
  TopologyView topologyView(connView);

  using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
  MIR m(topologyView, coordsetView, matsetView);
  m.execute(deviceMesh, options, deviceResult);

  // device->host
  {
    AXOM_ANNOTATE_SCOPE("device->host");
    utils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);
  }

  return 0;
}

//--------------------------------------------------------------------------------
/*!
 * \brief Run MIR on the quad input mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 *
 * \param hostMesh A conduit node that contains the test mesh.
 * \param options A conduit node that contains the test mesh.
 * \param hostResult A conduit node that will contain the MIR results.
 */
template <typename ExecSpace>
int runMIR_quad(const conduit::Node &hostMesh, const conduit::Node &options, conduit::Node &hostResult)
{
  AXOM_ANNOTATE_SCOPE("runMIR_quad");
  namespace utils = axom::bump::utilities;
  namespace views = axom::bump::views;
  SLIC_INFO(axom::fmt::format("Using policy {}", axom::execution_space<ExecSpace>::name()));

  // host->device
  conduit::Node deviceMesh;
  {
    AXOM_ANNOTATE_SCOPE("host->device");
    utils::copy<ExecSpace>(deviceMesh, hostMesh);
  }

  conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
  conduit::Node &n_topo = deviceMesh["topologies/mesh"];
  conduit::Node &n_matset = deviceMesh["matsets/mat"];
  auto connView = utils::make_array_view<int>(n_topo["elements/connectivity"]);

  // Make matset view. (There's often 1 more material so add 1)
  constexpr int MAXMATERIALS = 12;
  using MatsetView = views::UnibufferMaterialView<int, float, MAXMATERIALS + 1>;
  MatsetView matsetView;
  matsetView.set(utils::make_array_view<int>(n_matset["material_ids"]),
                 utils::make_array_view<float>(n_matset["volume_fractions"]),
                 utils::make_array_view<int>(n_matset["sizes"]),
                 utils::make_array_view<int>(n_matset["offsets"]),
                 utils::make_array_view<int>(n_matset["indices"]));

  // Make Coord/Topo views.
  conduit::Node deviceResult;
  auto coordsetView = views::make_explicit_coordset<float, 2>::view(n_coordset);
  using CoordsetView = decltype(coordsetView);
  using TopologyView = views::UnstructuredTopologySingleShapeView<views::QuadShape<int>>;
  TopologyView topologyView(connView);

  using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
  MIR m(topologyView, coordsetView, matsetView);
  m.execute(deviceMesh, options, deviceResult);

  // device->host
  {
    AXOM_ANNOTATE_SCOPE("device->host");
    utils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);
  }
  return 0;
}

//--------------------------------------------------------------------------------
/*!
 * \brief Run MIR on the hex input mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 *
 * \param hostMesh A conduit node that contains the test mesh.
 * \param options A conduit node that contains the test mesh.
 * \param hostResult A conduit node that will contain the MIR results.
 */
template <typename ExecSpace>
int runMIR_hex(const conduit::Node &hostMesh, const conduit::Node &options, conduit::Node &hostResult)
{
  AXOM_ANNOTATE_SCOPE("runMIR_hex");
  namespace utils = axom::bump::utilities;
  namespace views = axom::bump::views;
  SLIC_INFO(axom::fmt::format("Using policy {}", axom::execution_space<ExecSpace>::name()));

  // host->device
  conduit::Node deviceMesh;
  {
    AXOM_ANNOTATE_SCOPE("host->device");
    utils::copy<ExecSpace>(deviceMesh, hostMesh);
  }

  conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
  conduit::Node &n_topo = deviceMesh["topologies/mesh"];
  conduit::Node &n_matset = deviceMesh["matsets/mat"];
  auto connView = utils::make_array_view<int>(n_topo["elements/connectivity"]);

  // Make matset view. (There's often 1 more material so add 1)
  constexpr int MAXMATERIALS = 12;
  using MatsetView = views::UnibufferMaterialView<int, float, MAXMATERIALS + 1>;
  MatsetView matsetView;
  matsetView.set(utils::make_array_view<int>(n_matset["material_ids"]),
                 utils::make_array_view<float>(n_matset["volume_fractions"]),
                 utils::make_array_view<int>(n_matset["sizes"]),
                 utils::make_array_view<int>(n_matset["offsets"]),
                 utils::make_array_view<int>(n_matset["indices"]));

  // Make Coord/Topo views.
  conduit::Node deviceResult;
  auto coordsetView = views::make_explicit_coordset<float, 3>::view(n_coordset);
  using CoordsetView = decltype(coordsetView);
  using TopologyView = views::UnstructuredTopologySingleShapeView<views::HexShape<int>>;
  TopologyView topologyView(connView);

  using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
  MIR m(topologyView, coordsetView, matsetView);
  m.execute(deviceMesh, options, deviceResult);

  // device->host
  {
    AXOM_ANNOTATE_SCOPE("device->host");
    utils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);
  }

  return 0;
}

// Prototypes.
int runMIR_seq(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_omp(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_cuda(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_hip(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);

#endif
