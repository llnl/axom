// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/primal.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

namespace bump = axom::bump;
namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_coordset_extents
{
  static constexpr double eps = 1.e-10;

  static void test_uniform_2d()
  {
    const char *yaml = R"(
coords:
  type: uniform
  dims:
    i: 7
    j: 7
  origin:
    x: 1.
    y: 2.
  spacing:
    dx: 0.5
    dy: 0.5
)";

    conduit::Node n_coordset;
    initialize(yaml, n_coordset);

    auto coordsetView = views::make_uniform_coordset<2>::view(n_coordset["coords"]);
    using CoordsetView = decltype(coordsetView);
    bump::CoordsetExtents<ExecSpace, CoordsetView> exts(coordsetView);
    double extents[4];
    exts.execute(extents);

    const double expectedExtents[] = {1., 4., 2., 5.};
    EXPECT_NEAR(extents[0], expectedExtents[0], eps);
    EXPECT_NEAR(extents[1], expectedExtents[1], eps);
    EXPECT_NEAR(extents[2], expectedExtents[2], eps);
    EXPECT_NEAR(extents[3], expectedExtents[3], eps);
  }

  static void test_uniform_3d()
  {
    const char *yaml = R"(
coords:
  type: uniform
  dims:
    i: 7
    j: 7
    k: 7
  origin:
    x: -2.
    y: -1.
    z: 0.
  spacing:
    dx: 0.5
    dy: 1.
    dz: 1.5
)";

    conduit::Node n_coordset;
    initialize(yaml, n_coordset);

    auto coordsetView = views::make_uniform_coordset<3>::view(n_coordset["coords"]);
    using CoordsetView = decltype(coordsetView);
    bump::CoordsetExtents<ExecSpace, CoordsetView> exts(coordsetView);
    double extents[6];
    exts.execute(extents);

    const double expectedExtents[] = {-2., 1., -1., 5., 0., 9.};
    EXPECT_NEAR(extents[0], expectedExtents[0], eps);
    EXPECT_NEAR(extents[1], expectedExtents[1], eps);
    EXPECT_NEAR(extents[2], expectedExtents[2], eps);
    EXPECT_NEAR(extents[3], expectedExtents[3], eps);
    EXPECT_NEAR(extents[4], expectedExtents[4], eps);
    EXPECT_NEAR(extents[5], expectedExtents[5], eps);
  }

  static void test_rectilinear_2d()
  {
    const char *yaml = R"(
coords:
  type: rectilinear
  values:
    x: [-1., 0., 1., 2., 3., 5.]
    y: [2., 4., 6., 8., 10.]
)";

    conduit::Node n_coordset;
    initialize(yaml, n_coordset);

    auto coordsetView = views::make_rectilinear_coordset<double, 2>::view(n_coordset["coords"]);
    using CoordsetView = decltype(coordsetView);
    bump::CoordsetExtents<ExecSpace, CoordsetView> exts(coordsetView);
    double extents[4];
    exts.execute(extents);

    const double expectedExtents[] = {-1., 5., 2., 10.};
    EXPECT_NEAR(extents[0], expectedExtents[0], eps);
    EXPECT_NEAR(extents[1], expectedExtents[1], eps);
    EXPECT_NEAR(extents[2], expectedExtents[2], eps);
    EXPECT_NEAR(extents[3], expectedExtents[3], eps);
  }

  static void test_rectilinear_3d()
  {
    const char *yaml = R"(
coords:
  type: rectilinear
  values:
    x: [-1., 0., 1., 2., 3., 5.]
    y: [2., 4., 6., 8., 10.]
    z: [-3., 0., 3., 6., 9., 12.]
)";

    conduit::Node n_coordset;
    initialize(yaml, n_coordset);

    auto coordsetView = views::make_rectilinear_coordset<double, 3>::view(n_coordset["coords"]);
    using CoordsetView = decltype(coordsetView);
    bump::CoordsetExtents<ExecSpace, CoordsetView> exts(coordsetView);
    double extents[6];
    exts.execute(extents);

    const double expectedExtents[] = {-1., 5., 2., 10., -3., 12.};
    EXPECT_NEAR(extents[0], expectedExtents[0], eps);
    EXPECT_NEAR(extents[1], expectedExtents[1], eps);
    EXPECT_NEAR(extents[2], expectedExtents[2], eps);
    EXPECT_NEAR(extents[3], expectedExtents[3], eps);
    EXPECT_NEAR(extents[4], expectedExtents[4], eps);
    EXPECT_NEAR(extents[5], expectedExtents[5], eps);
  }

  static void test_explicit_2d()
  {
    const char *yaml = R"(
coords:
  type: explicit
  values:
    x: [-1., 0., 1., 2., -1., 0., 1.1, 2.2, -1., 0., 1., 2., -1., 0., 1., 2.2]
    y: [0., 0., 0., -0.01, 1., 1., 1., 1., 2., 2., 2., 2., 3., 3.1, 3., 3.]
)";

    conduit::Node n_coordset;
    initialize(yaml, n_coordset);

    auto coordsetView = views::make_explicit_coordset<double, 2>::view(n_coordset["coords"]);
    using CoordsetView = decltype(coordsetView);
    bump::CoordsetExtents<ExecSpace, CoordsetView> exts(coordsetView);
    double extents[4];
    exts.execute(extents);

    const double expectedExtents[] = {-1., 2.2, -0.01, 3.1};
    EXPECT_NEAR(extents[0], expectedExtents[0], eps);
    EXPECT_NEAR(extents[1], expectedExtents[1], eps);
    EXPECT_NEAR(extents[2], expectedExtents[2], eps);
    EXPECT_NEAR(extents[3], expectedExtents[3], eps);
  }

  static void test_explicit_3d()
  {
    const char *yaml = R"(
coords:
  type: explicit
  values:
    x: [-1., 0., 1., 2., -1., 0., 1.1, 2.2, -1., 0., 1., 2., -1., 0., 1., 2.2]
    y: [0., 0., 0., -0.01, 1., 1., 1., 1., 2., 2., 2., 2., 3., 3.1, 3., 3.]
    z: [0., 0., -0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.1]
)";

    conduit::Node n_coordset;
    initialize(yaml, n_coordset);

    auto coordsetView = views::make_explicit_coordset<double, 3>::view(n_coordset["coords"]);
    using CoordsetView = decltype(coordsetView);
    bump::CoordsetExtents<ExecSpace, CoordsetView> exts(coordsetView);
    double extents[6];
    exts.execute(extents);

    const double expectedExtents[] = {-1., 2.2, -0.01, 3.1, -0.1, 0.1};
    EXPECT_NEAR(extents[0], expectedExtents[0], eps);
    EXPECT_NEAR(extents[1], expectedExtents[1], eps);
    EXPECT_NEAR(extents[2], expectedExtents[2], eps);
    EXPECT_NEAR(extents[3], expectedExtents[3], eps);
    EXPECT_NEAR(extents[4], expectedExtents[4], eps);
    EXPECT_NEAR(extents[5], expectedExtents[5], eps);
  }

  static void initialize(const char *yaml, conduit::Node &n_device)
  {
    conduit::Node n_coordset;
    n_coordset.parse(yaml);
    utils::copy<ExecSpace>(n_device, n_coordset);
  }
};

//------------------------------------------------------------------------------
TEST(bump_coordset_extents, uniform2d_seq)
{
  AXOM_ANNOTATE_SCOPE("uniform2d_seq");
  test_coordset_extents<seq_exec>::test_uniform_2d();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_coordset_extents, uniform2d_omp)
{
  AXOM_ANNOTATE_SCOPE("uniform2d_omp");
  test_coordset_extents<omp_exec>::test_uniform_2d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_coordset_extents, uniform2d_cuda)
{
  AXOM_ANNOTATE_SCOPE("uniform2d_cuda");
  test_coordset_extents<cuda_exec>::test_uniform_2d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_coordset_extents, uniform2d_hip)
{
  AXOM_ANNOTATE_SCOPE("uniform2d_hip");
  test_coordset_extents<hip_exec>::test_uniform_2d();
}
#endif

//------------------------------------------------------------------------------
TEST(bump_coordset_extents, uniform3d_seq)
{
  AXOM_ANNOTATE_SCOPE("uniform3d_seq");
  test_coordset_extents<seq_exec>::test_uniform_3d();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_coordset_extents, uniform3d_omp)
{
  AXOM_ANNOTATE_SCOPE("uniform3d_omp");
  test_coordset_extents<omp_exec>::test_uniform_3d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_coordset_extents, uniform3d_cuda)
{
  AXOM_ANNOTATE_SCOPE("uniform3d_cuda");
  test_coordset_extents<cuda_exec>::test_uniform_3d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_coordset_extents, uniform3d_hip)
{
  AXOM_ANNOTATE_SCOPE("uniform3d_hip");
  test_coordset_extents<hip_exec>::test_uniform_3d();
}
#endif

//------------------------------------------------------------------------------
TEST(bump_coordset_extents, rectilinear2d_seq)
{
  AXOM_ANNOTATE_SCOPE("rectilinear2d_seq");
  test_coordset_extents<seq_exec>::test_rectilinear_2d();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_coordset_extents, rectilinear2d_omp)
{
  AXOM_ANNOTATE_SCOPE("rectilinear2d_omp");
  test_coordset_extents<omp_exec>::test_rectilinear_2d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_coordset_extents, rectilinear2d_cuda)
{
  AXOM_ANNOTATE_SCOPE("rectilinear2d_cuda");
  test_coordset_extents<cuda_exec>::test_rectilinear_2d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_coordset_extents, rectilinear2d_hip)
{
  AXOM_ANNOTATE_SCOPE("rectilinear2d_hip");
  test_coordset_extents<hip_exec>::test_rectilinear_2d();
}
#endif

//------------------------------------------------------------------------------
TEST(bump_coordset_extents, rectilinear3d_seq)
{
  AXOM_ANNOTATE_SCOPE("rectilinear3d_seq");
  test_coordset_extents<seq_exec>::test_rectilinear_3d();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_coordset_extents, rectilinear3d_omp)
{
  AXOM_ANNOTATE_SCOPE("rectilinear3d_omp");
  test_coordset_extents<omp_exec>::test_rectilinear_3d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_coordset_extents, rectilinear3d_cuda)
{
  AXOM_ANNOTATE_SCOPE("rectilinear3d_cuda");
  test_coordset_extents<cuda_exec>::test_rectilinear_3d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_coordset_extents, rectilinear3d_hip)
{
  AXOM_ANNOTATE_SCOPE("rectilinear3d_hip");
  test_coordset_extents<hip_exec>::test_rectilinear_3d();
}
#endif

//------------------------------------------------------------------------------
TEST(bump_coordset_extents, explicit2d_seq)
{
  AXOM_ANNOTATE_SCOPE("explicit2d_seq");
  test_coordset_extents<seq_exec>::test_explicit_2d();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_coordset_extents, explicit2d_omp)
{
  AXOM_ANNOTATE_SCOPE("explicit2d_omp");
  test_coordset_extents<omp_exec>::test_explicit_2d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_coordset_extents, explicit2d_cuda)
{
  AXOM_ANNOTATE_SCOPE("explicit2d_cuda");
  test_coordset_extents<cuda_exec>::test_explicit_2d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_coordset_extents, explicit2d_hip)
{
  AXOM_ANNOTATE_SCOPE("explicit2d_hip");
  test_coordset_extents<hip_exec>::test_explicit_2d();
}
#endif

//------------------------------------------------------------------------------
TEST(bump_coordset_extents, explicit3d_seq)
{
  AXOM_ANNOTATE_SCOPE("explicit3d_seq");
  test_coordset_extents<seq_exec>::test_explicit_3d();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_coordset_extents, explicit3d_omp)
{
  AXOM_ANNOTATE_SCOPE("explicit3d_omp");
  test_coordset_extents<omp_exec>::test_explicit_3d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_coordset_extents, explicit3d_cuda)
{
  AXOM_ANNOTATE_SCOPE("explicit3d_cuda");
  test_coordset_extents<cuda_exec>::test_explicit_3d();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_coordset_extents, explicit3d_hip)
{
  AXOM_ANNOTATE_SCOPE("explicit3d_hip");
  test_coordset_extents<hip_exec>::test_explicit_3d();
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

  // Define command line options.
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
