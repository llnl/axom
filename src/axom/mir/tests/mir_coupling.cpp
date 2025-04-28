// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/fmt.hpp"
#include "axom/mir.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>
#include <cstdlib>

namespace bputils = axom::mir::utilities::blueprint;

//------------------------------------------------------------------------------

// Uncomment to generate baselines
//#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
//#define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_coupled");
}

//------------------------------------------------------------------------------
/*!

The area defined with '#' characters is meant to be a strided structured mesh.

coarse

%42-------%43-------%44-------%45-------%46-------%47-------%48
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|   (3)   |   (3)   |     (3) |    (3)  |    (3)  |   (3)   |
%35-------%36-------%37#######%38#######%39#######%40-------%41
|         |         #\        #         #         #         |
|         |         #  \  (3) #         #         #         |
|         |         #    \    #         #         #         |
|         |         #      \  #         #         #         |
|   (2)   |   (2)   #  (2)   \#    (3)  #    (3)  #   (3)   |
%28-------%29------%30#######%31#######%32#######%33-------%34
|         |         #\ (2)    |\   (3)  |    (3)  #   (3)   |
|         |         #  \      |  \      |         #         |
|         |         #    \    |    - - -|- - - - -#- - - - -|
|         |         #      \  |         |         #         |
|   (1)   |   (1)   #  (1)   \|    (2)  |    (2)  #   (2)   |
%21-------%22-------%23#######%24#######%25#######%26-------%27
|         |   (0)   #\ (1)    |\   (2)  |    (2)  #   (2)   |
|         |         #  \      |(1)\     |         #         |
|         |         #    \- - |- - - - -|- - - - -#- - - - -|
|         |         #         |   (0)   |    (0)  #   (0)   |
|   (0)   |         #  (0)    |         |         #         |
%14-------%15-------%16#######%17#######%18#######%19-------%20
|   (0)   |   (0)   |  (0)    |   (0)   |    (0)  |   (0)   |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
%7--------%8--------%9--------%10-------%11-------%12-------%13
|   (0)   |   (0)   |  (0)    |   (0)   |    (0)  |   (0)   |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
%---------%---------%---------%---------%---------%---------%
0         1         2         3         4         5         6


fine - refines coarse with equal sized quads.

%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
0    1    2    3    4    5    6    7    8    9    10   11   12

*/
// NOTE: uniform coordsets used here for simplicity. Other coordset types can be used too.
const char *yaml = R"(
coordsets:
  coarse_coords:
    type: uniform
    dims:
      i: 7
      j: 7
    origin:
      x: 0.
      y: 0.
    spacing:
      dx: 1.
      dy: 1.
  fine_coords:
    type: uniform
    dims:
      i: 13
      j: 13
    origin:
      x: 0.
      y: 0.
    spacing:
      dx: 0.5
      dy: 0.5
topologies:
  coarse_strided:
    type: structured
    coordset: coarse_coords
    elements:
      dims:
       i: 3
       j: 3
       offsets: [2, 2]
       strides: [1, 7]
  coarse:
    type: structured
    coordset: coarse_coords
    elements:
      dims:
       i: 6
       j: 6
  fine:
    type: structured
    coordset: fine_coords
    elements:
      dims:
       i: 12
       j: 12
matsets:
  coarse_matset:
   topology: coarse
   material_map:
     mat0: 0
     mat1: 1
     mat2: 2
     mat3: 3
   material_ids: [0, 0, 0, 0, 0, 0,               0, 0, 0, 0, 0, 0,           0, 0, 0,1, 0,1,2, 0,2, 0,2,                                 1, 1, 1,2, 2,3, 2,3, 2,3,                           2, 2, 2,3, 3, 3, 3,             3, 3, 3, 3, 3, 3]
   volume_fractions: [1., 1., 1., 1., 1., 1.,     1., 1., 1., 1., 1., 1.,     1., 1., 0.625,0.375, 0.5,0.125,0.375, 0.5,0.5, 0.5,0.5,     1., 1., 0.5,0.5, 0.625,0.375, 0.5,0.5, 0.5,0.5,     1., 1., 0.5,0.5, 1., 1., 1.,    1., 1., 1., 1., 1., 1.]
   indices: [0, 1, 2, 3, 4, 5,                    6, 7, 8, 9, 10, 11,         12, 13, 14,15, 16,17,18, 19,20, 21,22,                      23, 24, 25,26, 27,28, 29,30, 31,32,                 33, 34, 35,36, 37, 38, 39,      40, 41, 42, 43, 44, 45]
   sizes: [1, 1, 1, 1, 1, 1,     1, 1, 1, 1, 1, 1,     1, 1, 2, 3, 2, 2,    1, 1, 2, 2, 2, 2,    1, 1, 2, 1, 1, 1,    1, 1, 1, 1, 1, 1]
   offsets: [0, 1, 2, 3, 4, 5,     6, 7, 8, 9, 10, 11,     12, 13, 14, 16, 19, 21,     23, 24, 25, 27, 29, 31,     33, 34, 35, 37, 38, 39,     40, 41, 42, 43, 44, 45]
)";

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
  static void test2D(const std::string &name, bool selectedZones = false, bool stridedStructured = false)
  {
    // Make the 2D input mesh.
    conduit::Node n_mesh;
    initialize(stridedStructured, n_mesh);

    // host->device
    conduit::Node n_dev;
    axom::mir::utilities::blueprint::copy<ExecSpace>(n_dev, n_mesh);

    // Do 2D MIR on the coarse mesh. The new objects will be added to n_mesh.
    std::string coarseName(stridedStructured ? "coarse_strided" : "coarse");
    mir2D(coarseName, n_dev, "postmir", n_dev, selectedZones, stridedStructured);

    // Map MIR output in n_mesh onto the fine mesh as a new matset.
    mapping2D(n_dev, n_dev, selectedZones, stridedStructured);

    if(!stridedStructured)
    {
      // As a check, run the generated fine matset through elvira again to make clean zones.
      mir2D("fine", n_dev, "check", n_dev, false, false);
    }

    // device->host
    conduit::Node hostResult;
    bputils::copy<seq_exec>(hostResult, n_dev);

#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  #if defined(AXOM_USE_HDF5)
    conduit::relay::io::blueprint::save_mesh(hostResult, name, "hdf5");
  #endif
    conduit::relay::io::save(hostResult, name + ".yaml", "yaml");
#endif

    // Handle baseline comparison.
    const auto paths = baselinePaths<ExecSpace>();
    std::string baselineName(yamlRoot(name));
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostResult);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostResult));
#endif
  }

private:
  static void initialize(bool stridedStructured, conduit::Node &n_mesh)
  {
    // Make the 2D input mesh.
    n_mesh.parse(yaml);

    if(stridedStructured)
    {
      // Adjust the mesh so we have the right names.
      n_mesh["coordsets/coarse_strided_coords"].set(n_mesh["coordsets/coarse_coords"]);
      n_mesh["topologies/coarse_strided/coordset"] = "coarse_strided_coords";
      n_mesh["matsets/coarse_strided_matset"].set(n_mesh["matsets/coarse_matset"]);
      n_mesh["matsets/coarse_strided_matset/topology"] = "coarse_strided";
      n_mesh.remove("coordsets/coarse_coords");
      n_mesh.remove("topologies/coarse");
      n_mesh.remove("matsets/coarse_matset");
    }
    else
    {
      n_mesh.remove("topologies/coarse_strided");
    }
  }

  static void mir2D(const std::string &input_prefix,
                    conduit::Node &n_input,
                    const std::string &output_prefix,
                    conduit::Node &n_output,
                    bool selectedZones,
                    bool stridedStructured)
  {
    namespace bputils = axom::mir::utilities::blueprint;

    // Wrap the coarse mesh in views.
    const conduit::Node &n_topology = n_input[axom::fmt::format("topologies/{}", input_prefix)];

    if(stridedStructured)
    {
      auto topologyView = axom::mir::views::make_strided_structured<2>::view(n_topology);
      mir2D(topologyView, input_prefix, n_input, output_prefix, n_output, selectedZones, stridedStructured);
    }
    else
    {
      auto topologyView = axom::mir::views::make_structured<2>::view(n_topology);
      mir2D(topologyView, input_prefix, n_input, output_prefix, n_output, selectedZones, stridedStructured);
    }
  }

  template <typename TopologyView>
  static void mir2D(TopologyView topologyView,
                    const std::string &input_prefix,
                    conduit::Node &n_input,
                    const std::string &output_prefix,
                    conduit::Node &n_output,
                    bool selectedZones,
                    bool stridedStructured)
  {
    namespace bputils = axom::mir::utilities::blueprint;

    // Wrap the coarse mesh in views.
    const conduit::Node &n_coordset = n_input[axom::fmt::format("coordsets/{}_coords", input_prefix)];
    const conduit::Node &n_matset = n_input[axom::fmt::format("matsets/{}_matset", input_prefix)];

    auto coordsetView = axom::mir::views::make_uniform_coordset<2>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    // Get the indexing policy from the TopologyView type.
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    auto matsetView =
      axom::mir::views::make_unibuffer_matset<std::int64_t, double, 4>::view(n_matset);
    using MatsetView = decltype(matsetView);

    // Do MIR on the mesh.
    using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node options;
    // Select that matset we'll operate on.
    options["matset"] = axom::fmt::format("{}_matset", input_prefix);
    // Change the names of the topology, coordset, and matset in the output.
    options["topologyName"] = output_prefix;
    options["coordsetName"] = axom::fmt::format("{}_coords", output_prefix);
    options["matsetName"] = axom::fmt::format("{}_matset", output_prefix);
    if(selectedZones)
    {
      selectCoarseZones2D(stridedStructured, options);
    }
    m.execute(n_input, options, n_output);
  }

  /// Add a list of selected zones to the options going into 2D mir.
  static void selectCoarseZones2D(bool stridedStructured, conduit::Node &n_options)
  {
    if(stridedStructured)
    {
      n_options["selectedZones"].set(std::vector<axom::IndexType> {0, 1, 3, 4, 6, 7});
    }
    else
    {
      n_options["selectedZones"].set(
        std::vector<axom::IndexType> {13, 14, 15, 19, 20, 21, 25, 26, 27});
    }
  }

  static void selectFineZones2D(bool stridedStructured, conduit::Node &n_options)
  {
    // These selected zones are on the fine mesh and are used for TopologyMapper.
    conduit::Node &n_selectedZones = n_options["target/selectedZones"];

    if(stridedStructured)
    {
      n_selectedZones.set(std::vector<axom::IndexType> {52,  53,  54,  55,  64,  65,  66,  67,
                                                        76,  77,  78,  79,  88,  89,  90,  91,
                                                        100, 101, 102, 103, 112, 113, 114, 115});
    }
    else
    {
      n_selectedZones.set(std::vector<axom::IndexType> {
        50, 51, 52, 53, 54, 55, 62, 63, 64,  65,  66,  67,  74,  75,  76,  77,  78,  79,
        86, 87, 88, 89, 90, 91, 98, 99, 100, 101, 102, 103, 110, 111, 112, 113, 114, 115});
    }
  }

  static void mapping2D(conduit::Node &n_src,
                        conduit::Node &n_target,
                        bool selectedZones,
                        bool stridedStructured)
  {
    namespace bputils = axom::mir::utilities::blueprint;

    // Wrap the source mesh from (coarse MIR output).
    const conduit::Node &n_src_coordset = n_src["coordsets/postmir_coords"];
    const conduit::Node &n_src_topology = n_src["topologies/postmir"];
    const conduit::Node &n_src_matset = n_src["matsets/postmir_matset"];

    auto srcCoordsetView = axom::mir::views::make_explicit_coordset<double, 2>::view(n_src_coordset);
    using SrcCoordsetView = decltype(srcCoordsetView);

    // 2D Elvira makes polygonal meshes
    using SrcShapeType = axom::mir::views::PolygonShape<axom::IndexType>;
    auto srcTopologyView =
      axom::mir::views::make_unstructured_single_shape<SrcShapeType>::view(n_src_topology);
    using SrcTopologyView = decltype(srcTopologyView);

    auto srcMatsetView =
      axom::mir::views::make_unibuffer_matset<std::int64_t, double, 4>::view(n_src_matset);
    using SrcMatsetView = decltype(srcMatsetView);

    // Wrap the target mesh (fine)
    const conduit::Node &n_target_coordset = n_target["coordsets/fine_coords"];
    const conduit::Node &n_target_topology = n_target["topologies/fine"];

    auto targetCoordsetView = axom::mir::views::make_uniform_coordset<2>::view(n_target_coordset);
    using TargetCoordsetView = decltype(targetCoordsetView);

    auto targetTopologyView = axom::mir::views::make_structured<2>::view(n_target_topology);
    using TargetTopologyView = decltype(targetTopologyView);

    // Make new a new matset on the target topology to record material overlaps.
    using Mapper =
      bputils::TopologyMapper<ExecSpace, SrcTopologyView, SrcCoordsetView, SrcMatsetView, TargetTopologyView, TargetCoordsetView>;
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
    if(selectedZones)
    {
      selectFineZones2D(stridedStructured, n_opts);
    }
    mapper.execute(n_src, n_opts, n_target);
  }
};

//------------------------------------------------------------------------------
TEST(mir_coupling, coupling_2D_sz0_ss0_seq)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss0_seq");
  test_coupling<seq_exec>::test2D("coupling_2D_sz0_ss0", false, false);
}
TEST(mir_coupling, coupling_2D_sz0_ss1_seq)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss1_seq");
  test_coupling<seq_exec>::test2D("coupling_2D_sz0_ss1", false, true);
}
TEST(mir_coupling, coupling_2D_sz1_ss0_seq)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss0_seq");
  test_coupling<seq_exec>::test2D("coupling_2D_sz1_ss0", true, false);
}
TEST(mir_coupling, coupling_2D_sz1_ss1_seq)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss1_seq");
  test_coupling<seq_exec>::test2D("coupling_2D_sz1_ss1", true, true);
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_coupling, coupling_2D_sz0_ss0_omp)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss0_omp");
  test_coupling<omp_exec>::test2D("coupling_2D_sz0_ss0", false, false);
}
TEST(mir_coupling, coupling_2D_sz0_ss1_omp)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss1_omp");
  test_coupling<omp_exec>::test2D("coupling_2D_sz0_ss1", false, true);
}
TEST(mir_coupling, coupling_2D_sz1_ss0_omp)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss0_omp");
  test_coupling<omp_exec>::test2D("coupling_2D_sz1_ss0", true, false);
}
TEST(mir_coupling, coupling_2D_sz1_ss1_omp)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss1_omp");
  test_coupling<omp_exec>::test2D("coupling_2D_sz1_ss1", true, true);
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_coupling, coupling_2D_sz0_ss0_cuda)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss0_cuda");
  test_coupling<cuda_exec>::test2D("coupling_2D_sz0_ss0", false, false);
}
TEST(mir_coupling, coupling_2D_sz0_ss1_cuda)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss1_cuda");
  test_coupling<cuda_exec>::test2D("coupling_2D_sz0_ss1", false, true);
}
TEST(mir_coupling, coupling_2D_sz1_ss0_cuda)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss0_cuda");
  test_coupling<cuda_exec>::test2D("coupling_2D_sz1_ss0", true, false);
}
TEST(mir_coupling, coupling_2D_sz1_ss1_cuda)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss1_cuda");
  test_coupling<cuda_exec>::test2D("coupling_2D_sz1_ss1", true, true);
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_coupling, coupling_2D_sz0_ss0_hip)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss0_hip");
  test_coupling<hip_exec>::test2D("coupling_2D_sz0_ss0", false, false);
}
TEST(mir_coupling, coupling_2D_sz0_ss1_hip)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss1_hip");
  test_coupling<hip_exec>::test2D("coupling_2D_sz0_ss1", false, true);
}
TEST(mir_coupling, coupling_2D_sz1_ss0_hip)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss0_hip");
  test_coupling<hip_exec>::test2D("coupling_2D_sz1_ss0", true, false);
}
TEST(mir_coupling, coupling_2D_sz1_ss1_hip)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss1_hip");
  test_coupling<hip_exec>::test2D("coupling_2D_sz1_ss1", true, true);
}
#endif

//------------------------------------------------------------------------------
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1)
    ;
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
