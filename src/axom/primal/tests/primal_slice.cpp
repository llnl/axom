// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"

#include "axom/core/Array.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"

#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/operators/slice.hpp"

namespace primal = axom::primal;

namespace
{

constexpr double EPS = 1e-10;

using Point3D = primal::Point<double, 3>;

template <int N, typename PolygonType>
void expect_vertex_set(const PolygonType& poly, const axom::StackArray<Point3D, N>& expected)
{
  ASSERT_EQ(poly.numVertices(), N);

  bool matched[N] = {false};

  for(int i = 0; i < N; ++i)
  {
    bool found = false;
    for(int j = 0; j < N; ++j)
    {
      if(!matched[j] && poly[i].isNearlyEqual(expected[j], EPS))
      {
        matched[j] = true;
        found = true;
        break;
      }
    }

    EXPECT_TRUE(found) << "Unexpected polygon vertex: " << poly[i];
  }
}

template <typename PolygonType, typename PlaneType>
void expect_normal_aligned(const PolygonType& poly, const PlaneType& plane)
{
  ASSERT_GE(poly.numVertices(), 3);
  EXPECT_GT(poly.normal().dot(plane.getNormal()), 0.);
}

template <typename ExecSpace>
void check_slice_policy()
{
  using TetType = primal::Tetrahedron<double, 3>;
  using PlaneType = primal::Plane<double, 3>;
  using PolygonType = primal::Polygon<double, 3, primal::PolygonArray::Static, 4>;

  const TetType tet {Point3D {-1., -1., -1.},
                     Point3D {1., 1., -1.},
                     Point3D {-1., -1., 1.},
                     Point3D {-1., 1., 1.}};
  const PlaneType plane({0., 0., 1.}, 0.);

  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int kernel_allocator = axom::execution_space<ExecSpace>::allocatorID();

  axom::Array<PolygonType> polys_device(1, 1, kernel_allocator);
  auto polys_view = polys_device.view();

  axom::Array<double> areas_device(1, 1, kernel_allocator);
  auto areas_view = areas_device.view();

  axom::for_all<ExecSpace>(
    1,
    AXOM_LAMBDA(int i) {
      polys_view[i] = primal::slice<double, primal::PolygonArray::Static, 4>(tet, plane);
      areas_view[i] = polys_view[i].area();
    });

  axom::Array<PolygonType> polys_host(polys_device, host_allocator);
  axom::Array<double> areas_host(areas_device, host_allocator);

  EXPECT_EQ(polys_host[0].numVertices(), 4);
  EXPECT_NEAR(areas_host[0], 1., EPS);
  expect_vertex_set(
    polys_host[0],
    axom::StackArray<Point3D, 4> {
      {Point3D {-1., -1., 0.}, Point3D {-1., 0., 0.}, Point3D {0., 0., 0.}, Point3D {0., 1., 0.}}});
  expect_normal_aligned(polys_host[0], plane);
}

template <typename ExecSpace>
void check_slice_degenerate_policy()
{
  using TetType = primal::Tetrahedron<double, 3>;
  using PlaneType = primal::Plane<double, 3>;
  using PolygonType = primal::Polygon<double, 3, primal::PolygonArray::Static, 4>;

  const TetType tet {Point3D {0., 0., 0.},
                     Point3D {1., 0., 0.},
                     Point3D {0., 1., 0.},
                     Point3D {0., 0., 1.}};

  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int kernel_allocator = axom::execution_space<ExecSpace>::allocatorID();

  axom::Array<PolygonType> polys_device(3, 3, kernel_allocator);
  auto polys_view = polys_device.view();

  axom::Array<double> areas_device(3, 3, kernel_allocator);
  auto areas_view = areas_device.view();

  axom::for_all<ExecSpace>(
    3,
    AXOM_LAMBDA(int i) {
      PlaneType plane;

      // Edge case: the plane intersects the tet on a vertex.
      if(i == 0)
      {
        plane = PlaneType({1., 1., 1.}, 0.);
      }

      // Edge case: the plane intersects the tet on an edge, here the edge from
      // (0,0,0) to (1,0,0).
      if(i == 1)
      {
        plane = PlaneType({0., 1., 1.}, 0.);
      }

      // Edge case: the plane intersects the tet on a face, here the face z = 0.
      if(i == 2)
      {
        plane = PlaneType({0., 0., 1.}, 0.);
      }

      polys_view[i] = primal::slice<double, primal::PolygonArray::Static, 4>(tet, plane);
      areas_view[i] = polys_view[i].area();
    });

  axom::Array<PolygonType> polys_host(polys_device, host_allocator);
  axom::Array<double> areas_host(areas_device, host_allocator);

  EXPECT_EQ(polys_host[0].numVertices(), 1);
  EXPECT_NEAR(areas_host[0], 0., EPS);
  expect_vertex_set(polys_host[0], axom::StackArray<Point3D, 1> {{Point3D {0., 0., 0.}}});

  EXPECT_EQ(polys_host[1].numVertices(), 2);
  EXPECT_NEAR(areas_host[1], 0., EPS);
  expect_vertex_set(polys_host[1],
                    axom::StackArray<Point3D, 2> {{Point3D {0., 0., 0.}, Point3D {1., 0., 0.}}});

  EXPECT_EQ(polys_host[2].numVertices(), 3);
  EXPECT_NEAR(areas_host[2], 0.5, EPS);
  expect_vertex_set(polys_host[2],
                    axom::StackArray<Point3D, 3> {
                      {Point3D {0., 0., 0.}, Point3D {1., 0., 0.}, Point3D {0., 1., 0.}}});
  expect_normal_aligned(polys_host[2], PlaneType({0., 0., 1.}, 0.));
}

}  // namespace

TEST(primal_slice, tet_plane_slice_dynamic)
{
  using TetType = primal::Tetrahedron<double, 3>;
  using PlaneType = primal::Plane<double, 3>;

  const TetType tet {Point3D {0., 0., 0.},
                     Point3D {1., 0., 0.},
                     Point3D {0., 1., 0.},
                     Point3D {0., 0., 1.}};

  const auto poly = primal::slice(tet, PlaneType({0., 0., 1.}, 0.25));
  EXPECT_EQ(poly.numVertices(), 3);
  EXPECT_NEAR(poly.area(), 0.28125, EPS);
  expect_vertex_set(poly,
                    axom::StackArray<Point3D, 3> {
                      {Point3D {0., 0., 0.25}, Point3D {0.75, 0., 0.25}, Point3D {0., 0.75, 0.25}}});
  expect_normal_aligned(poly, PlaneType({0., 0., 1.}, 0.25));

  const auto empty_poly = primal::slice(tet, PlaneType({0., 0., 1.}, 2.));
  EXPECT_EQ(empty_poly.numVertices(), 0);
}

TEST(primal_slice, tet_plane_slice_degenerate_dynamic)
{
  using TetType = primal::Tetrahedron<double, 3>;
  using PlaneType = primal::Plane<double, 3>;

  const TetType tet {Point3D {0., 0., 0.},
                     Point3D {1., 0., 0.},
                     Point3D {0., 1., 0.},
                     Point3D {0., 0., 1.}};

  // Edge case: the plane intersects the tet on a vertex.
  const auto vertex_poly = primal::slice(tet, PlaneType({1., 1., 1.}, 0.));
  EXPECT_EQ(vertex_poly.numVertices(), 1);
  EXPECT_NEAR(vertex_poly.area(), 0., EPS);
  expect_vertex_set(vertex_poly, axom::StackArray<Point3D, 1> {{Point3D {0., 0., 0.}}});

  // Edge case: the plane intersects the tet on an edge.
  const auto edge_poly = primal::slice(tet, PlaneType({0., 1., 1.}, 0.));
  EXPECT_EQ(edge_poly.numVertices(), 2);
  EXPECT_NEAR(edge_poly.area(), 0., EPS);
  expect_vertex_set(edge_poly,
                    axom::StackArray<Point3D, 2> {{Point3D {0., 0., 0.}, Point3D {1., 0., 0.}}});

  // Edge case: the plane intersects the tet on a face.
  const auto face_poly = primal::slice(tet, PlaneType({0., 0., 1.}, 0.));
  EXPECT_EQ(face_poly.numVertices(), 3);
  EXPECT_NEAR(face_poly.area(), 0.5, EPS);
  expect_vertex_set(face_poly,
                    axom::StackArray<Point3D, 3> {
                      {Point3D {0., 0., 0.}, Point3D {1., 0., 0.}, Point3D {0., 1., 0.}}});
  expect_normal_aligned(face_poly, PlaneType({0., 0., 1.}, 0.));
}

TEST(primal_slice, tet_plane_slice_orientation_follows_plane)
{
  using TetType = primal::Tetrahedron<double, 3>;
  using PlaneType = primal::Plane<double, 3>;

  const TetType tet {Point3D {-1., -1., -1.},
                     Point3D {1., 1., -1.},
                     Point3D {-1., -1., 1.},
                     Point3D {-1., 1., 1.}};

  // Flipping the slicing plane should preserve the vertex set while reversing
  // the polygon orientation so the Polygon orientation matches that of the
  // Plane.
  const PlaneType plane_pos({0., 0., 1.}, 0.);
  const PlaneType plane_neg({0., 0., -1.}, 0.);

  const auto poly_pos = primal::slice(tet, plane_pos);
  const auto poly_neg = primal::slice(tet, plane_neg);

  expect_vertex_set(
    poly_pos,
    axom::StackArray<Point3D, 4> {
      {Point3D {-1., -1., 0.}, Point3D {-1., 0., 0.}, Point3D {0., 0., 0.}, Point3D {0., 1., 0.}}});
  expect_vertex_set(
    poly_neg,
    axom::StackArray<Point3D, 4> {
      {Point3D {-1., -1., 0.}, Point3D {-1., 0., 0.}, Point3D {0., 0., 0.}, Point3D {0., 1., 0.}}});
  expect_normal_aligned(poly_pos, plane_pos);
  expect_normal_aligned(poly_neg, plane_neg);
  EXPECT_LT(poly_pos.normal().dot(poly_neg.normal()), 0.);
}

TEST(primal_slice, tet_plane_slice_near_boundary_tolerance)
{
  using TetType = primal::Tetrahedron<double, 3>;
  using PlaneType = primal::Plane<double, 3>;

  const TetType tet {Point3D {0., 0., 0.},
                     Point3D {1., 0., 0.},
                     Point3D {0., 1., 0.},
                     Point3D {0., 0., 1.}};

  // A plane within the current fuzzy zero tolerance of the base face behaves
  // like the coincident-face case and returns that face.
  const auto fuzzy_face = primal::slice(tet, PlaneType({0., 0., 1.}, 5e-11));
  EXPECT_EQ(fuzzy_face.numVertices(), 3);
  EXPECT_NEAR(fuzzy_face.area(), 0.5, EPS);
  expect_vertex_set(fuzzy_face,
                    axom::StackArray<Point3D, 3> {
                      {Point3D {0., 0., 0.}, Point3D {1., 0., 0.}, Point3D {0., 1., 0.}}});

  // Moving the plane clearly away from that fuzzy interval produces the
  // expected small triangle strictly above the base face instead.
  const auto nearby_slice = primal::slice(tet, PlaneType({0., 0., 1.}, 1e-6));
  EXPECT_EQ(nearby_slice.numVertices(), 3);
  EXPECT_NEAR(nearby_slice.area(), 0.4999990000005, EPS);
  expect_vertex_set(
    nearby_slice,
    axom::StackArray<Point3D, 3> {
      {Point3D {0., 0., 1e-6}, Point3D {0.999999, 0., 1e-6}, Point3D {0., 0.999999, 1e-6}}});
  expect_normal_aligned(nearby_slice, PlaneType({0., 0., 1.}, 1e-6));
}

TEST(primal_slice, tet_plane_slice_degenerate_tet_repeated_vertex)
{
  using TetType = primal::Tetrahedron<double, 3>;
  using PlaneType = primal::Plane<double, 3>;

  const TetType degenerate_tet {Point3D {0., 0., 0.},
                                Point3D {1., 0., 0.},
                                Point3D {0., 1., 0.},
                                Point3D {0., 0., 0.}};

  // A tetrahedron with a repeated vertex degenerates to a triangle; slicing by
  // that triangle's supporting plane should still return the unique face
  // vertices rather than duplicates.
  const auto poly = primal::slice(degenerate_tet, PlaneType({0., 0., 1.}, 0.));
  EXPECT_EQ(poly.numVertices(), 3);
  EXPECT_NEAR(poly.area(), 0.5, EPS);
  expect_vertex_set(poly,
                    axom::StackArray<Point3D, 3> {
                      {Point3D {0., 0., 0.}, Point3D {1., 0., 0.}, Point3D {0., 1., 0.}}});
  expect_normal_aligned(poly, PlaneType({0., 0., 1.}, 0.));
}

TEST(primal_slice, tet_plane_slice_invalid_plane)
{
  using TetType = primal::Tetrahedron<double, 3>;
  using PlaneType = primal::Plane<double, 3>;

  const TetType tet {Point3D {0., 0., 0.},
                     Point3D {1., 0., 0.},
                     Point3D {0., 1., 0.},
                     Point3D {0., 0., 1.}};

  // An invalid plane with zero normal does not define a meaningful slice, so
  // the result should be empty.
  const auto poly = primal::slice(tet, PlaneType());
  EXPECT_EQ(poly.numVertices(), 0);
}

TEST(primal_slice, tet_plane_slice_seq) { check_slice_policy<axom::SEQ_EXEC>(); }

TEST(primal_slice, tet_plane_slice_degenerate_seq)
{
  check_slice_degenerate_policy<axom::SEQ_EXEC>();
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)

  #ifdef AXOM_USE_OPENMP
TEST(primal_slice, tet_plane_slice_omp) { check_slice_policy<axom::OMP_EXEC>(); }

TEST(primal_slice, tet_plane_slice_degenerate_omp)
{
  check_slice_degenerate_policy<axom::OMP_EXEC>();
}
  #endif

  #ifdef AXOM_USE_CUDA
AXOM_CUDA_TEST(primal_slice, tet_plane_slice_cuda) { check_slice_policy<axom::CUDA_EXEC<256>>(); }

AXOM_CUDA_TEST(primal_slice, tet_plane_slice_degenerate_cuda)
{
  check_slice_degenerate_policy<axom::CUDA_EXEC<256>>();
}
  #endif

  #ifdef AXOM_USE_HIP
TEST(primal_slice, tet_plane_slice_hip) { check_slice_policy<axom::HIP_EXEC<256>>(); }

TEST(primal_slice, tet_plane_slice_degenerate_hip)
{
  check_slice_degenerate_policy<axom::HIP_EXEC<256>>();
}
  #endif

#endif

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
