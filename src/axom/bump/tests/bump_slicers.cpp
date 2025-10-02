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
template <typename ExecSpace>
struct test_matset_slice
{
  static void test()
  {
    conduit::Node hostMatset;
    create(hostMatset);

    // host->device
    conduit::Node deviceMatset;
    utils::copy<ExecSpace>(deviceMatset, hostMatset);

    axom::Array<axom::IndexType> ids {{1, 3, 5}};
    axom::Array<axom::IndexType> selectedZones(3, 3, axom::execution_space<ExecSpace>::allocatorID());
    axom::copy(selectedZones.data(), ids.data(), 3 * sizeof(axom::IndexType));

    using MatsetView = axom::bump::views::UnibufferMaterialView<conduit::int64, conduit::float64, 3>;
    MatsetView matsetView;
    matsetView.set(utils::make_array_view<conduit::int64>(deviceMatset["material_ids"]),
                   utils::make_array_view<conduit::float64>(deviceMatset["volume_fractions"]),
                   utils::make_array_view<conduit::int64>(deviceMatset["sizes"]),
                   utils::make_array_view<conduit::int64>(deviceMatset["offsets"]),
                   utils::make_array_view<conduit::int64>(deviceMatset["indices"]));

    // Slice it.
    bump::MatsetSlicer<ExecSpace, MatsetView> slicer(matsetView);
    conduit::Node newDeviceMatset;
    bump::SliceData slice;
    slice.m_indicesView = selectedZones.view();
    slicer.execute(slice, deviceMatset, newDeviceMatset);

    // device->host
    conduit::Node newHostMatset;
    utils::copy<axom::SEQ_EXEC>(newHostMatset, newDeviceMatset);

    // Expected answers.
    const axom::Array<conduit::int64> sizes {{2, 1, 2}};
    const axom::Array<conduit::int64> offsets {{0, 2, 3}};
    const axom::Array<conduit::int64> indices {{0, 1, 2, 3, 4}};
    const axom::Array<conduit::int64> material_ids {{1, 2, 2, 2, 3}};
    const axom::Array<conduit::float64> volume_fractions {{0.5, 0.5, 1.0, 0.8, 0.2}};

    EXPECT_EQ(conduit::DataType::INT64_ID, newHostMatset["material_ids"].dtype().id());
    EXPECT_EQ(conduit::DataType::FLOAT64_ID, newHostMatset["volume_fractions"].dtype().id());

    EXPECT_TRUE(
      compare_views(sizes.view(), utils::make_array_view<conduit::int64>(newHostMatset["sizes"])));
    EXPECT_TRUE(compare_views(offsets.view(),
                              utils::make_array_view<conduit::int64>(newHostMatset["offsets"])));
    EXPECT_TRUE(compare_views(indices.view(),
                              utils::make_array_view<conduit::int64>(newHostMatset["indices"])));
    EXPECT_TRUE(
      compare_views(material_ids.view(),
                    utils::make_array_view<conduit::int64>(newHostMatset["material_ids"])));
    EXPECT_TRUE(
      compare_views(volume_fractions.view(),
                    utils::make_array_view<conduit::float64>(newHostMatset["volume_fractions"])));
  }

  static void create(conduit::Node &matset)
  {
    /*
      8-------9------10------11
      |  2/1  | 1/0.1 | 2/0.8 |
      |       | 2/0.5 | 3/0.2 |
      |       | 3/0.4 |       |
      4-------5-------6-------7
      |       | 1/0.5 | 1/0.2 |
      |  1/1  | 2/0.5 | 2/0.8 |
      |       |       |       |
      0-------1-------2-------3
      */
    const char *yaml = R"xx(
topology: mesh
material_map:
  a: 1
  b: 2
  c: 3
material_ids: [1, 1,2, 1,2, 2, 1,2,3, 2,3]
volume_fractions: [1., 0.5,0.5, 0.2,0.8, 1., 0.1,0.5,0.4, 0.8,0.2]
indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
sizes: [1, 2, 2, 1, 3, 2]
offsets: [0, 1, 3, 5, 6, 9]
)xx";

    matset.parse(yaml);
  }
};

TEST(bump_slicers, matsetslice_seq) { test_matset_slice<seq_exec>::test(); }

#if defined(AXOM_USE_OPENMP)
TEST(bump_slicers, matsetslice_omp) { test_matset_slice<omp_exec>::test(); }
#endif

#if defined(AXOM_USE_CUDA)
TEST(bump_slicers, matsetslice_cuda) { test_matset_slice<cuda_exec>::test(); }
#endif

#if defined(AXOM_USE_HIP)
TEST(bump_slicers, matsetslice_hip) { test_matset_slice<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace, typename Func>
void test_coordsetslicer(const conduit::Node &hostCoordset, Func &&makeView)
{
  axom::Array<axom::IndexType> ids {{0, 1, 2, 4, 5, 6}};

  const auto nnodes = ids.size();
  axom::Array<axom::IndexType> selectedNodes(nnodes,
                                             nnodes,
                                             axom::execution_space<ExecSpace>::allocatorID());
  axom::copy(selectedNodes.data(), ids.data(), nnodes * sizeof(axom::IndexType));

  bump::SliceData slice;
  slice.m_indicesView = selectedNodes.view();

  // host->device
  conduit::Node deviceCoordset;
  utils::copy<ExecSpace>(deviceCoordset, hostCoordset);

  // Make a view.
  auto coordsetView = makeView(deviceCoordset);
  using CoordsetView = decltype(coordsetView);

  // Pull out selected nodes
  bump::CoordsetSlicer<ExecSpace, CoordsetView> slicer(coordsetView);
  conduit::Node newDeviceCoordset;
  slicer.execute(slice, deviceCoordset, newDeviceCoordset);

  // device->host
  conduit::Node newHostCoordset;
  utils::copy<axom::SEQ_EXEC>(newHostCoordset, newDeviceCoordset);

  // We get an explicit coordset out of the slicer.
  const axom::Array<conduit::float64> x {{0., 1., 2., 0., 1., 2.}};
  const axom::Array<conduit::float64> y {{0., 0., 0., 1., 1., 1.}};
  EXPECT_TRUE(
    compare_views(x.view(), utils::make_array_view<conduit::float64>(newHostCoordset["values/x"])));
  EXPECT_TRUE(
    compare_views(y.view(), utils::make_array_view<conduit::float64>(newHostCoordset["values/y"])));
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct coordsetslicer_explicit
{
  static void test()
  {
    /*
      8---9--10--11
      |   |   |   |
      4---5---6---7
      |   |   |   |
      0---1---2---3
      */
    const char *yaml = R"xx(
type: explicit
values:
  x: [0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.]
  y: [0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.]
)xx";

    conduit::Node coordset;
    coordset.parse(yaml);

    auto makeView = [](const conduit::Node &deviceCoordset) {
      return axom::bump::views::make_explicit_coordset<conduit::float64, 2>::view(deviceCoordset);
    };

    test_coordsetslicer<seq_exec>(coordset, makeView);
  }
};

TEST(bump_slicers, coordsetslicer_explicit_seq) { coordsetslicer_explicit<seq_exec>::test(); }

#if defined(AXOM_USE_OPENMP)
TEST(bump_slicers, coordsetslicer_explicit_omp) { coordsetslicer_explicit<omp_exec>::test(); }
#endif

#if defined(AXOM_USE_CUDA)
TEST(bump_slicers, coordsetslicer_explicit_cuda) { coordsetslicer_explicit<cuda_exec>::test(); }
#endif

#if defined(AXOM_USE_HIP)
TEST(bump_slicers, coordsetslicer_explicit_hip) { coordsetslicer_explicit<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct coordsetslicer_rectilinear
{
  static void test()
  {
    /*
      8---9--10--11
      |   |   |   |
      4---5---6---7
      |   |   |   |
      0---1---2---3
      */
    const char *yaml = R"xx(
type: rectilinear
values:
  x: [0., 1., 2., 3.]
  y: [0., 1., 2.]
)xx";

    conduit::Node coordset;
    coordset.parse(yaml);

    auto makeView = [](const conduit::Node &deviceCoordset) {
      return axom::bump::views::make_rectilinear_coordset<conduit::float64, 2>::view(deviceCoordset);
    };
    test_coordsetslicer<ExecSpace>(coordset, makeView);
  }
};

TEST(bump_slicers, coordsetslicer_rectilinear_seq) { coordsetslicer_rectilinear<seq_exec>::test(); }

#if defined(AXOM_USE_OPENMP)
TEST(bump_slicers, coordsetslicer_rectilinear_omp) { coordsetslicer_rectilinear<omp_exec>::test(); }
#endif

#if defined(AXOM_USE_CUDA)
TEST(bump_slicers, coordsetslicer_rectilinear_cuda)
{
  coordsetslicer_rectilinear<cuda_exec>::test();
}
#endif

#if defined(AXOM_USE_HIP)
TEST(bump_slicers, coordsetslicer_rectilinear_hip) { coordsetslicer_rectilinear<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct coordsetslicer_uniform
{
  static void test()
  {
    /*
      8---9--10--11
      |   |   |   |
      4---5---6---7
      |   |   |   |
      0---1---2---3
      */
    const char *yaml = R"xx(
type: uniform
dims:
  i: 4
  j: 3
)xx";

    conduit::Node coordset;
    coordset.parse(yaml);

    auto makeView = [](const conduit::Node &deviceCoordset) {
      return axom::bump::views::make_uniform_coordset<2>::view(deviceCoordset);
    };
    test_coordsetslicer<ExecSpace>(coordset, makeView);
  }
};

TEST(bump_slicers, coordsetslicer_uniform_seq) { coordsetslicer_uniform<seq_exec>::test(); }

#if defined(AXOM_USE_OPENMP)
TEST(bump_slicers, coordsetslicer_uniform_omp) { coordsetslicer_uniform<omp_exec>::test(); }
#endif

#if defined(AXOM_USE_CUDA)
TEST(bump_slicers, coordsetslicer_uniform_cuda) { coordsetslicer_uniform<cuda_exec>::test(); }
#endif

#if defined(AXOM_USE_HIP)
TEST(bump_slicers, coordsetslicer_uniform_hip) { coordsetslicer_uniform<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_fieldslicer
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    // _bump_utilities_fieldslicer_begin
    std::vector<axom::IndexType> indices {0, 1, 2, 7, 8, 9};
    axom::Array<axom::IndexType> sliceIndices(indices.size(),
                                              indices.size(),
                                              axom::execution_space<ExecSpace>::allocatorID());
    axom::copy(sliceIndices.data(), indices.data(), sizeof(axom::IndexType) * indices.size());

    bump::SliceData slice;
    slice.m_indicesView = sliceIndices.view();

    conduit::Node slicedMesh;
    bump::FieldSlicer<ExecSpace> fs;
    fs.execute(slice, deviceMesh["fields/scalar"], slicedMesh["fields/scalar"]);
    fs.execute(slice, deviceMesh["fields/vector"], slicedMesh["fields/vector"]);
    // _bump_utilities_fieldslicer_end

    // device->host
    conduit::Node hostSlicedMesh;
    utils::copy<axom::SEQ_EXEC>(hostSlicedMesh, slicedMesh);

    std::vector<double> resultX {0., 1., 2., 7., 8., 9.};
    std::vector<double> resultY {0., 10., 20., 70., 80., 90.};

    EXPECT_EQ(hostSlicedMesh["fields/scalar/topology"].as_string(), "mesh");
    EXPECT_EQ(hostSlicedMesh["fields/scalar/association"].as_string(), "element");
    EXPECT_EQ(hostSlicedMesh["fields/scalar/values"].dtype().number_of_elements(), indices.size());
    for(size_t i = 0; i < indices.size(); i++)
    {
      const auto acc = hostSlicedMesh["fields/scalar/values"].as_double_accessor();
      EXPECT_EQ(acc[i], resultX[i]);
    }

    EXPECT_EQ(hostSlicedMesh["fields/vector/topology"].as_string(), "mesh");
    EXPECT_EQ(hostSlicedMesh["fields/vector/association"].as_string(), "element");
    EXPECT_EQ(hostSlicedMesh["fields/vector/values/x"].dtype().number_of_elements(), indices.size());
    EXPECT_EQ(hostSlicedMesh["fields/vector/values/y"].dtype().number_of_elements(), indices.size());
    for(size_t i = 0; i < indices.size(); i++)
    {
      const auto x = hostSlicedMesh["fields/vector/values/x"].as_double_accessor();
      const auto y = hostSlicedMesh["fields/vector/values/y"].as_double_accessor();
      EXPECT_EQ(x[i], resultX[i]);
      EXPECT_EQ(y[i], resultY[i]);
    }
  }

  static void create(conduit::Node &fields)
  {
    const char *yaml = R"xx(
fields:
  scalar: 
    topology: mesh
    association: element
    values: [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]
  vector: 
    topology: mesh
    association: element
    values:
      x: [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]
      y: [0., 10., 20., 30., 40., 50., 60., 70., 80., 90.]
)xx";
    fields.parse(yaml);
  }
};

TEST(bump_slicers, fieldslicer_seq) { test_fieldslicer<seq_exec>::test(); }
#if defined(AXOM_USE_OPENMP)
TEST(bump_slicers, fieldslicer_omp) { test_fieldslicer<omp_exec>::test(); }
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_slicers, fieldslicer_cuda) { test_fieldslicer<cuda_exec>::test(); }
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_slicers, fieldslicer_hip) { test_fieldslicer<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
