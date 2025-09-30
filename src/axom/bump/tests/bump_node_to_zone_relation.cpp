// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump.hpp"

#include "axom/bump/tests/blueprint_testing_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"

#include <iostream>
#include <algorithm>

namespace bump = axom::bump;
namespace utils = axom::bump::utilities;

//------------------------------------------------------------------------------
template <typename ExecSpace, typename IndexT = int>
struct test_node_to_zone_relation_builder
{
  static void test(const conduit::Node &hostMesh)
  {
    // host -> device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);
    // _bump_utilities_n2zrel_begin
    const conduit::Node &deviceTopo = deviceMesh["topologies/mesh"];
    const conduit::Node &deviceCoordset = deviceMesh["coordsets/coords"];

    // Run the algorithm on the device
    conduit::Node deviceRelation;
    bump::NodeToZoneRelationBuilder<ExecSpace> n2z;
    n2z.execute(deviceTopo, deviceCoordset, deviceRelation);
    // _bump_utilities_n2zrel_end

    // device -> host
    conduit::Node hostRelation;
    utils::copy<seq_exec>(hostRelation, deviceRelation);

    // Expected answers
    // clang-format off
    const int zones[] = {
     0,
     0, 1,
     1, 2,
     2,
     0, 3,        // NOTE: these are sorted here
     0, 1, 3, 4,
     1, 2, 4, 5,
     2, 5,
     3,
     3, 4,
     4, 5,
     5
    };
    const int sizes[] = {1, 2, 2, 1, 2, 4, 4, 2, 1, 2, 2, 1};
    const int offsets[] = {0, 1, 3, 5, 6, 8, 12, 16, 18, 19, 21, 23};
    // clang-format on

    // Compare answers.
    compareRelation(hostRelation,
                    axom::ArrayView<const int>(zones, sizeof(zones) / sizeof(int)),
                    axom::ArrayView<const int>(sizes, sizeof(sizes) / sizeof(int)),
                    axom::ArrayView<const int>(offsets, sizeof(offsets) / sizeof(int)));
  }

  static void compareRelation(const conduit::Node &hostRelation,
                              const axom::ArrayView<const int> &zones,
                              const axom::ArrayView<const int> &sizes,
                              const axom::ArrayView<const int> &offsets)
  {
    const auto zonesView = utils::make_array_view<IndexT>(hostRelation["zones"]);
    const auto sizesView = utils::make_array_view<IndexT>(hostRelation["sizes"]);
    const auto offsetsView = utils::make_array_view<IndexT>(hostRelation["offsets"]);
    EXPECT_EQ(sizesView.size(), sizes.size());
    EXPECT_EQ(offsetsView.size(), offsets.size());
    for(axom::IndexType i = 0; i < sizesView.size(); i++)
    {
      EXPECT_EQ(sizes[i], sizesView[i]);
      EXPECT_EQ(offsets[i], offsetsView[i]);
    }
    for(axom::IndexType i = 0; i < sizesView.size(); i++)
    {
      // Sort the result so we can compare to the expected answer.
      IndexT *begin = zonesView.data() + offsetsView[i];
      IndexT *end = zonesView.data() + offsetsView[i] + sizesView[i];
      std::sort(begin, end);

      for(int j = 0; j < sizesView[i]; j++)
      {
        EXPECT_EQ(zones[offsets[i] + j], zonesView[offsetsView[i] + j]);
      }
    }
  }
};

TEST(bump_node_to_zone_relation, n2zrel_unstructured_seq)
{
  /*
    8---9--10--11
    |   |   |   |
    4---5---6---7
    |   |   |   |
    0---1---2---3
    */
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::blueprint::testing::data::braid("quads", dims, mesh);
  test_node_to_zone_relation_builder<seq_exec>::test(mesh);
}

#if defined(AXOM_USE_OPENMP)
TEST(bump_node_to_zone_relation, n2zrel_unstructured_omp)
{
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::blueprint::testing::data::braid("quads", dims, mesh);
  test_node_to_zone_relation_builder<omp_exec>::test(mesh);
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(bump_node_to_zone_relation, n2zrel_unstructured_cuda)
{
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::blueprint::testing::data::braid("quads", dims, mesh);
  test_node_to_zone_relation_builder<cuda_exec>::test(mesh);
}
#endif

#if defined(AXOM_USE_HIP)
TEST(bump_node_to_zone_relation, n2zrel_unstructured_hip)
{
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::blueprint::testing::data::braid("quads", dims, mesh);
  test_node_to_zone_relation_builder<hip_exec>::test(mesh);
}
#endif

TEST(bump_node_to_zone_relation, n2zrel_rectilinear_seq)
{
  /*
    8---9--10--11
    |   |   |   |
    4---5---6---7
    |   |   |   |
    0---1---2---3
    */
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::blueprint::testing::data::braid("rectilinear", dims, mesh);
  test_node_to_zone_relation_builder<seq_exec, conduit::index_t>::test(mesh);
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_node_to_zone_relation, n2zrel_rectilinear_omp)
{
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::blueprint::testing::data::braid("rectilinear", dims, mesh);
  test_node_to_zone_relation_builder<omp_exec, conduit::index_t>::test(mesh);
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_node_to_zone_relation, n2zrel_rectilinear_cuda)
{
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::blueprint::testing::data::braid("rectilinear", dims, mesh);
  test_node_to_zone_relation_builder<cuda_exec, conduit::index_t>::test(mesh);
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_node_to_zone_relation, n2zrel_rectilinear_hip)
{
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::blueprint::testing::data::braid("rectilinear", dims, mesh);
  test_node_to_zone_relation_builder<hip_exec, conduit::index_t>::test(mesh);
}
#endif

template <typename ExecSpace, typename IndexT = int>
struct test_node_to_zone_relation_builder_polyhedral
  : public test_node_to_zone_relation_builder<ExecSpace, IndexT>
{
  using SuperClass = test_node_to_zone_relation_builder<ExecSpace, IndexT>;

  static void test(const conduit::Node &hostMesh)
  {
    // host -> device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);
    const conduit::Node &deviceTopo = deviceMesh["topologies/mesh"];
    const conduit::Node &deviceCoordset = deviceMesh["coordsets/coords"];

    // Run the algorithm on the device
    conduit::Node deviceRelation;
    bump::NodeToZoneRelationBuilder<ExecSpace> n2z;
    n2z.execute(deviceTopo, deviceCoordset, deviceRelation);

    // device -> host
    conduit::Node hostRelation;
    utils::copy<seq_exec>(hostRelation, deviceRelation);

    // Expected answers
    // clang-format off
    const int zones[] = {
     /*node 0*/ 0,
     /*node 1*/ 0, 1,
     /*node 2*/ 1,
     /*node 3*/ 0, 2,
     /*node 4*/ 0, 1, 2, 3,
     /*node 5*/ 1, 3,
     /*node 6*/ 2,
     /*node 7*/ 2, 3,
     /*node 8*/ 3,

     /*node 9*/ 0, 4,
     /*node 10*/ 0, 1, 4, 5,
     /*node 11*/ 1, 5,
     /*node 12*/ 0, 2, 4, 6,
     /*node 13*/ 0, 1, 2, 3, 4, 5, 6, 7,
     /*node 14*/ 1, 3, 5, 7,
     /*node 15*/ 2, 6,
     /*node 16*/ 2, 3, 6, 7,
     /*node 17*/ 3, 7,

     /*node 18*/ 4,
     /*node 19*/ 4, 5,
     /*node 20*/ 5,
     /*node 21*/ 4, 6,
     /*node 22*/ 4, 5, 6, 7,
     /*node 23*/ 5, 7,
     /*node 24*/ 6,
     /*node 25*/ 6, 7,
     /*node 26*/ 7
    };
    const int sizes[] = {1, 2, 1, 2, 4, 2, 1, 2, 1,
                         2, 4, 2, 4, 8, 4, 2, 4, 2,
                         1, 2, 1, 2, 4, 2, 1, 2, 1
                        };
    const int offsets[] = {0, 1, 3, 4, 6, 10, 12, 13, 15,
                           16, 18, 22, 24, 28, 36, 40, 42, 46,
                           48, 49, 51, 52, 54, 58, 60, 61, 63};
    // clang-format on

    // Compare answers.
    SuperClass::compareRelation(hostRelation,
                                axom::ArrayView<const int>(zones, sizeof(zones) / sizeof(int)),
                                axom::ArrayView<const int>(sizes, sizeof(sizes) / sizeof(int)),
                                axom::ArrayView<const int>(offsets, sizeof(offsets) / sizeof(int)));
  }

  static void create(conduit::Node &mesh)
  {
    conduit::blueprint::mesh::examples::basic("polyhedra", 3, 3, 3, mesh);
    // Make sure all the types are the same.
    conduit::blueprint::mesh::utils::convert(
      mesh,
      conduit::DataType::int32(),
      std::vector<std::string> {{"topologies/mesh/elements/connectivity",
                                 "topologies/mesh/elements/sizes",
                                 "topologies/mesh/elements/offsets",
                                 "topologies/mesh/subelements/connectivity",
                                 "topologies/mesh/subelements/sizes",
                                 "topologies/mesh/subelements/offsets"}});
  }
};

TEST(bump_node_to_zone_relation, n2zrel_polyhedral_seq)
{
  conduit::Node mesh;
  test_node_to_zone_relation_builder_polyhedral<seq_exec, conduit::int32>::create(mesh);
  test_node_to_zone_relation_builder_polyhedral<seq_exec, conduit::int32>::test(mesh);
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_node_to_zone_relation, n2zrel_polyhedral_omp)
{
  conduit::Node mesh;
  test_node_to_zone_relation_builder_polyhedral<omp_exec, conduit::int32>::create(mesh);
  test_node_to_zone_relation_builder_polyhedral<omp_exec, conduit::int32>::test(mesh);
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_node_to_zone_relation, n2zrel_polyhedral_cuda)
{
  conduit::Node mesh;
  test_node_to_zone_relation_builder_polyhedral<cuda_exec, conduit::int32>::create(mesh);
  test_node_to_zone_relation_builder_polyhedral<cuda_exec, conduit::int32>::test(mesh);
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_node_to_zone_relation, n2zrel_polyhedral_hip)
{
  conduit::Node mesh;
  test_node_to_zone_relation_builder_polyhedral<hip_exec, conduit::int32>::create(mesh);
  test_node_to_zone_relation_builder_polyhedral<hip_exec, conduit::int32>::test(mesh);
}
#endif

//------------------------------------------------------------------------------
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1);
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::CLI::App app;
#if defined(AXOM_USE_CALIPER)
  std::string annotationMode("none");
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif
  bool handlerEnabled = false;
  app.add_flag("--handler", handlerEnabled, "Enable Conduit handler.");

  // Parse command line options.
  try
  {
    app.parse(argc, argv);

#if defined(AXOM_USE_CALIPER)
    axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(annotationMode);
#endif

    axom::slic::SimpleLogger logger;  // create & initialize test logger,
    if(handlerEnabled)
    {
      conduit::utils::set_error_handler(conduit_debug_err_handler);
    }

    result = RUN_ALL_TESTS();
  }
  catch(axom::CLI::CallForHelp &e)
  {
    std::cout << app.help() << std::endl;
    result = 0;
  }
  catch(axom::CLI::ParseError &e)
  {
    // Handle other parsing errors
    std::cerr << e.what() << std::endl;
    result = app.exit(e);
  }

  return result;
}
