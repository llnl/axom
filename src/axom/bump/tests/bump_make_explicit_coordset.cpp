// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/bump/MakeExplicitCoordset.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/dispatch_coordset.hpp"

#include "conduit.hpp"

namespace
{

using seq_exec = axom::SEQ_EXEC;

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
using omp_exec = axom::OMP_EXEC;
#else
using omp_exec = seq_exec;
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
constexpr int CUDA_BLOCK_SIZE = 256;
using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
#else
using cuda_exec = seq_exec;
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
constexpr int HIP_BLOCK_SIZE = 64;
using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
#else
using hip_exec = seq_exec;
#endif

namespace bump = axom::bump;
namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

struct CoordsetData
{
  int dimension {0};
  axom::IndexType size {0};
  axom::StackArray<axom::Array<double>, 3> components;
};

CoordsetData collectCoordsetData(const conduit::Node& coordset)
{
  CoordsetData data;
  views::dispatch_coordset(coordset, [&](auto coordsetView) {
    data.dimension = coordsetView.dimension();
    data.size = coordsetView.size();

    for(int dim = 0; dim < data.dimension; ++dim)
    {
      data.components[dim].resize(data.size);
    }

    for(axom::IndexType i = 0; i < data.size; ++i)
    {
      const auto pt = coordsetView[i];
      for(int dim = 0; dim < data.dimension; ++dim)
      {
        data.components[dim][i] = pt[dim];
      }
    }
  });

  return data;
}

void expectExplicitCoordsetMatches(const conduit::Node& coordset, const CoordsetData& expected)
{
  ASSERT_TRUE(coordset.has_path("type"));
  EXPECT_EQ(coordset["type"].as_string(), "explicit");

  views::dispatch_explicit_coordset(coordset, [&](auto coordsetView) {
    EXPECT_EQ(coordsetView.dimension(), expected.dimension);
    EXPECT_EQ(coordsetView.size(), expected.size);

    for(axom::IndexType i = 0; i < expected.size; ++i)
    {
      const auto pt = coordsetView[i];
      for(int dim = 0; dim < expected.dimension; ++dim)
      {
        EXPECT_NEAR(pt[dim], expected.components[dim][i], 1e-12);
      }
    }
  });
}

template <typename ExecSpace>
void checkMakeExplicitCoordset(const char* yaml)
{
  conduit::Node hostCoordset;
  hostCoordset.parse(yaml);
  const CoordsetData expected = collectCoordsetData(hostCoordset);

  conduit::Node deviceCoordset;
  utils::copy<ExecSpace>(deviceCoordset, hostCoordset);

  bump::MakeExplicitCoordset<ExecSpace>::execute(deviceCoordset);

  conduit::Node convertedCoordset;
  utils::copy<seq_exec>(convertedCoordset, deviceCoordset);
  expectExplicitCoordsetMatches(convertedCoordset, expected);
}

template <typename ExecSpace>
struct test_make_explicit_coordset
{
  static void uniform_2d()
  {
    const char* yaml = R"xx(
type: uniform
dims:
  i: 4
  j: 3
origin:
  x: 1.5
  y: -2.
spacing:
  dx: 0.25
  dy: 1.5
)xx";

    checkMakeExplicitCoordset<ExecSpace>(yaml);
  }

  static void rectilinear_3d()
  {
    const char* yaml = R"xx(
type: rectilinear
values:
  x: [-1., 0.5, 2.]
  y: [2., 5.]
  z: [0., 1.5, 3.5]
)xx";

    checkMakeExplicitCoordset<ExecSpace>(yaml);
  }

  static void explicit_noop()
  {
    const char* yaml = R"xx(
type: explicit
values:
  x: [0., 1., 3.]
  y: [2., 4., 8.]
  z: [-1., -2., -4.]
)xx";

    checkMakeExplicitCoordset<ExecSpace>(yaml);
  }
};

TEST(bump_make_explicit_coordset, uniform_2d_seq)
{
  test_make_explicit_coordset<seq_exec>::uniform_2d();
}

TEST(bump_make_explicit_coordset, rectilinear_3d_seq)
{
  test_make_explicit_coordset<seq_exec>::rectilinear_3d();
}

TEST(bump_make_explicit_coordset, explicit_noop_seq)
{
  test_make_explicit_coordset<seq_exec>::explicit_noop();
}

#if defined(AXOM_USE_OPENMP)
TEST(bump_make_explicit_coordset, uniform_2d_omp)
{
  test_make_explicit_coordset<omp_exec>::uniform_2d();
}

TEST(bump_make_explicit_coordset, rectilinear_3d_omp)
{
  test_make_explicit_coordset<omp_exec>::rectilinear_3d();
}

TEST(bump_make_explicit_coordset, explicit_noop_omp)
{
  test_make_explicit_coordset<omp_exec>::explicit_noop();
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(bump_make_explicit_coordset, uniform_2d_cuda)
{
  test_make_explicit_coordset<cuda_exec>::uniform_2d();
}

TEST(bump_make_explicit_coordset, rectilinear_3d_cuda)
{
  test_make_explicit_coordset<cuda_exec>::rectilinear_3d();
}

TEST(bump_make_explicit_coordset, explicit_noop_cuda)
{
  test_make_explicit_coordset<cuda_exec>::explicit_noop();
}
#endif

#if defined(AXOM_USE_HIP)
TEST(bump_make_explicit_coordset, uniform_2d_hip)
{
  test_make_explicit_coordset<hip_exec>::uniform_2d();
}

TEST(bump_make_explicit_coordset, rectilinear_3d_hip)
{
  test_make_explicit_coordset<hip_exec>::rectilinear_3d();
}

TEST(bump_make_explicit_coordset, explicit_noop_hip)
{
  test_make_explicit_coordset<hip_exec>::explicit_noop();
}
#endif

}  // namespace
