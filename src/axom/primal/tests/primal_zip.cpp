// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/utils/ZipIndexable.hpp"
#include "axom/primal/utils/ZipPoint.hpp"
#include "axom/primal/utils/ZipVector.hpp"
#include "axom/primal/utils/ZipBoundingBox.hpp"
#include "axom/primal/utils/ZipRay.hpp"

#include "gtest/gtest.h"

using namespace axom;

//------------------------------------------------------------------------------

template <typename PrimitiveType>
void check_default_ctor()
{
  using ZipType = primal::ZipIndexable<PrimitiveType>;

  ZipType zip;
}

template <typename ExecSpace, typename PrimitiveType>
void check_zip_points_3d()
{
  using ZipType = primal::ZipIndexable<PrimitiveType>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // create arrays of data
  constexpr int N = 8;
  double* x = axom::allocate<double>(N);
  double* y = axom::allocate<double>(N);
  double* z = axom::allocate<double>(N);

  bool* valid = axom::allocate<bool>(N);
  bool valid_host[N];

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      x[idx] = (idx % 2) + 1;
      y[idx] = ((idx / 2) % 2) + 1;
      z[idx] = (idx / 4) + 1;
    });

  ZipType it {{x, y, z}};

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) { valid[idx] = (it[idx] == PrimitiveType {x[idx], y[idx], z[idx]}); });

  axom::copy(&valid_host, valid, N * sizeof(bool));

  for(int i = 0; i < N; i++)
  {
    EXPECT_TRUE(valid_host[i]);
  }

  axom::deallocate(x);
  axom::deallocate(y);
  axom::deallocate(z);
  axom::deallocate(valid);
  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecSpace>
void check_zip_points_2d_from_3d()
{
  using PointType = primal::Point<double, 2>;
  using ZipType = primal::ZipIndexable<PointType>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // create arrays of data
  constexpr int N = 4;
  double* x = axom::allocate<double>(N);
  double* y = axom::allocate<double>(N);
  double* z = nullptr;

  bool* valid = axom::allocate<bool>(N);
  bool valid_host[N];

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      x[idx] = (idx % 2) + 1;
      y[idx] = (idx / 2) + 1;
    });

  ZipType it = {{x, y, z}};

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) { valid[idx] = (it[idx] == PointType {x[idx], y[idx]}); });

  axom::copy(&valid_host, valid, N * sizeof(bool));

  for(int i = 0; i < N; i++)
  {
    EXPECT_TRUE(valid_host[i]);
  }

  axom::deallocate(x);
  axom::deallocate(y);
  axom::deallocate(valid);
  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecSpace>
void check_zip_vectors_2d_from_3d()
{
  using PointType = primal::Vector<double, 2>;
  using ZipType = primal::ZipIndexable<PointType>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // create arrays of data
  constexpr int N = 4;
  double* x = axom::allocate<double>(N);
  double* y = axom::allocate<double>(N);
  double* z = nullptr;

  bool* valid = axom::allocate<bool>(N);
  bool valid_host[N];

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      x[idx] = (idx % 2) + 1;
      y[idx] = (idx / 2) + 1;
    });

  ZipType it {{x, y, z}};

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) { valid[idx] = (it[idx] == PointType {x[idx], y[idx]}); });

  axom::copy(&valid_host, valid, N * sizeof(bool));

  for(int i = 0; i < N; i++)
  {
    EXPECT_TRUE(valid_host[i]);
  }

  axom::deallocate(x);
  axom::deallocate(y);
  axom::deallocate(valid);
  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecSpace>
void check_zip_bbs_3d()
{
  using BoxType = primal::BoundingBox<double, 3>;
  using PointType = typename BoxType::PointType;
  using ZipType = primal::ZipIndexable<BoxType>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // create arrays of data
  constexpr int N = 8;
  double* xmin = axom::allocate<double>(N);
  double* ymin = axom::allocate<double>(N);
  double* zmin = axom::allocate<double>(N);

  double* xmax = axom::allocate<double>(N);
  double* ymax = axom::allocate<double>(N);
  double* zmax = axom::allocate<double>(N);

  bool* valid = axom::allocate<bool>(N);
  bool valid_host[N];

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      xmin[idx] = (idx % 2) + 1;
      ymin[idx] = ((idx / 2) % 2) + 1;
      zmin[idx] = (idx / 4) + 1;

      xmax[idx] = xmin[idx] + 1.0;
      ymax[idx] = ymin[idx] + 1.0;
      zmax[idx] = zmin[idx] + 1.0;
    });

  ZipType it {{xmin, ymin, zmin}, {xmax, ymax, zmax}};

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      BoxType actual;
      actual.addPoint(PointType {xmin[idx], ymin[idx], zmin[idx]});
      actual.addPoint(PointType {xmax[idx], ymax[idx], zmax[idx]});
      valid[idx] = (it[idx] == actual);
    });

  axom::copy(&valid_host, valid, N * sizeof(bool));

  for(int i = 0; i < N; i++)
  {
    EXPECT_TRUE(valid_host[i]);
  }

  axom::deallocate(xmin);
  axom::deallocate(ymin);
  axom::deallocate(zmin);

  axom::deallocate(xmax);
  axom::deallocate(ymax);
  axom::deallocate(zmax);

  axom::deallocate(valid);
  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecSpace>
void check_zip_bbs_2d_from_3d()
{
  using BoxType = primal::BoundingBox<double, 2>;
  using PointType = typename BoxType::PointType;
  using ZipType = primal::ZipIndexable<BoxType>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // create arrays of data
  constexpr int N = 4;
  double* xmin = axom::allocate<double>(N);
  double* ymin = axom::allocate<double>(N);
  double* zmin = nullptr;

  double* xmax = axom::allocate<double>(N);
  double* ymax = axom::allocate<double>(N);
  double* zmax = nullptr;

  bool* valid = axom::allocate<bool>(N);
  bool valid_host[N];

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      xmin[idx] = (idx % 2) + 1;
      ymin[idx] = (idx / 2) + 1;

      xmax[idx] = xmin[idx] + 1.0;
      ymax[idx] = ymin[idx] + 1.0;
    });

  ZipType it {{xmin, ymin, zmin}, {xmax, ymax, zmax}};

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      BoxType actual;
      actual.addPoint(PointType {xmin[idx], ymin[idx]});
      actual.addPoint(PointType {xmax[idx], ymax[idx]});
      valid[idx] = (it[idx] == actual);
    });

  axom::copy(&valid_host, valid, N * sizeof(bool));

  for(int i = 0; i < N; i++)
  {
    EXPECT_TRUE(valid_host[i]);
  }

  axom::deallocate(xmin);
  axom::deallocate(ymin);

  axom::deallocate(xmax);
  axom::deallocate(ymax);

  axom::deallocate(valid);
  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecSpace>
void check_zip_rays_3d()
{
  using RayType = primal::Ray<double, 3>;
  using PointType = typename RayType::PointType;
  using VectorType = typename RayType::VectorType;
  using ZipType = primal::ZipIndexable<RayType>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // create arrays of data
  constexpr int N = 4;
  double* xo = axom::allocate<double>(N);
  double* yo = axom::allocate<double>(N);
  double* zo = axom::allocate<double>(N);

  double* xd = axom::allocate<double>(N);
  double* yd = axom::allocate<double>(N);
  double* zd = axom::allocate<double>(N);

  bool* valid = axom::allocate<bool>(N);
  bool valid_host[N];

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      if(idx < 2)
      {
        xo[idx] = yo[idx] = zo[idx] = 1.0;
      }
      else
      {
        xo[idx] = yo[idx] = zo[idx] = -1.0;
      }

      if(idx % 2 == 0)
      {
        xd[idx] = yd[idx] = zd[idx] = 1.0;
      }
      else
      {
        xd[idx] = yd[idx] = zd[idx] = -1.0;
      }
    });

  ZipType it {{xo, yo, zo}, {xd, yd, zd}};

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      PointType orig {xo[idx], yo[idx], zo[idx]};
      VectorType dir {xd[idx], yd[idx], zd[idx]};
      RayType actual(orig, dir);
      valid[idx] =
        (it[idx].origin() == actual.origin()) && (it[idx].direction() == actual.direction());
    });

  axom::copy(&valid_host, valid, N * sizeof(bool));

  for(int i = 0; i < N; i++)
  {
    EXPECT_TRUE(valid_host[i]);
  }

  axom::deallocate(xo);
  axom::deallocate(yo);
  axom::deallocate(zo);

  axom::deallocate(xd);
  axom::deallocate(yd);
  axom::deallocate(zd);

  axom::deallocate(valid);
  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecSpace>
void check_zip_rays_2d_from_3d()
{
  using RayType = primal::Ray<double, 2>;
  using PointType = typename RayType::PointType;
  using VectorType = typename RayType::VectorType;
  using ZipType = primal::ZipIndexable<RayType>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // create arrays of data
  constexpr int N = 4;
  double* xo = axom::allocate<double>(N);
  double* yo = axom::allocate<double>(N);
  double* zo = nullptr;

  double* xd = axom::allocate<double>(N);
  double* yd = axom::allocate<double>(N);
  double* zd = nullptr;

  bool* valid = axom::allocate<bool>(N);
  bool valid_host[N];

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      if(idx < 2)
      {
        xo[idx] = yo[idx] = 1.0;
      }
      else
      {
        xo[idx] = yo[idx] = -1.0;
      }

      if(idx % 2 == 0)
      {
        xd[idx] = yd[idx] = 1.0;
      }
      else
      {
        xd[idx] = yd[idx] = -1.0;
      }
    });

  ZipType it {{xo, yo, zo}, {xd, yd, zd}};

  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(int idx) {
      PointType orig {xo[idx], yo[idx]};
      VectorType dir {xd[idx], yd[idx]};
      RayType actual(orig, dir);
      valid[idx] =
        (it[idx].origin() == actual.origin()) && (it[idx].direction() == actual.direction());
    });

  axom::copy(&valid_host, valid, N * sizeof(bool));

  for(int i = 0; i < N; i++)
  {
    EXPECT_TRUE(valid_host[i]);
  }

  axom::deallocate(xo);
  axom::deallocate(yo);

  axom::deallocate(xd);
  axom::deallocate(yd);

  axom::deallocate(valid);
  axom::setDefaultAllocator(current_allocator);
}

TEST(primal_zip, default_ctor)
{
  check_default_ctor<primal::Point<double, 1>>();
  check_default_ctor<primal::Point<double, 2>>();
  check_default_ctor<primal::Point<double, 3>>();

  check_default_ctor<primal::Vector<double, 1>>();
  check_default_ctor<primal::Vector<double, 2>>();
  check_default_ctor<primal::Vector<double, 3>>();

  check_default_ctor<primal::Ray<double, 1>>();
  check_default_ctor<primal::Ray<double, 2>>();
  check_default_ctor<primal::Ray<double, 3>>();

  check_default_ctor<primal::BoundingBox<double, 1>>();
  check_default_ctor<primal::BoundingBox<double, 2>>();
  check_default_ctor<primal::BoundingBox<double, 3>>();
}

TEST(primal_zip, zip_points_3d)
{
  using PointType = primal::Point<double, 3>;
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_points_3d<ExecSpace, PointType>();
}

TEST(primal_zip, zip_points_2d_from_3d)
{
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_points_2d_from_3d<ExecSpace>();
}

TEST(primal_zip, zip_vectors_3d)
{
  using PointType = primal::Vector<double, 3>;
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_points_3d<ExecSpace, PointType>();
}

TEST(primal_zip, zip_vectors_2d_from_3d)
{
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_vectors_2d_from_3d<ExecSpace>();
}

TEST(primal_zip, zip_bbs_3d)
{
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_bbs_3d<ExecSpace>();
}

TEST(primal_zip, zip_bbs_2d_from_3d)
{
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_bbs_2d_from_3d<ExecSpace>();
}

TEST(primal_zip, zip_rays_3d)
{
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_rays_3d<ExecSpace>();
}

TEST(primal_zip, zip_rays_2d_from_3d)
{
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_rays_2d_from_3d<ExecSpace>();
}

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
AXOM_CUDA_TEST(primal_zip, zip_points_3d_cuda)
{
  using PointType = primal::Point<double, 3>;
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_points_3d<ExecSpace, PointType>();
}

TEST(primal_zip, zip_points_2d_from_3d_cuda)
{
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_points_2d_from_3d<ExecSpace>();
}

AXOM_CUDA_TEST(primal_zip, zip_vectors_3d_cuda)
{
  using PointType = primal::Vector<double, 3>;
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_points_3d<ExecSpace, PointType>();
}

TEST(primal_zip, zip_vectors_2d_from_3d_cuda)
{
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_vectors_2d_from_3d<ExecSpace>();
}

AXOM_CUDA_TEST(primal_zip, zip_bbs_3d_cuda)
{
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_bbs_3d<ExecSpace>();
}

TEST(primal_zip, zip_bbs_2d_from_3d_cuda)
{
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_bbs_2d_from_3d<ExecSpace>();
}

TEST(primal_zip, zip_rays_3d_cuda)
{
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_rays_3d<ExecSpace>();
}

TEST(primal_zip, zip_rays_2d_from_3d_cuda)
{
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_rays_2d_from_3d<ExecSpace>();
}
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(primal_zip, zip_points_3d_hip)
{
  using PointType = primal::Point<double, 3>;
  using ExecSpace = axom::HIP_EXEC<256>;

  check_zip_points_3d<ExecSpace, PointType>();
}

TEST(primal_zip, zip_points_2d_from_3d_hip)
{
  using ExecSpace = axom::HIP_EXEC<256>;

  check_zip_points_2d_from_3d<ExecSpace>();
}

TEST(primal_zip, zip_vectors_3d_hip)
{
  using PointType = primal::Vector<double, 3>;
  using ExecSpace = axom::HIP_EXEC<256>;

  check_zip_points_3d<ExecSpace, PointType>();
}

TEST(primal_zip, zip_vectors_2d_from_3d_hip)
{
  using ExecSpace = axom::HIP_EXEC<256>;

  check_zip_vectors_2d_from_3d<ExecSpace>();
}

TEST(primal_zip, zip_bbs_3d_hip)
{
  using ExecSpace = axom::HIP_EXEC<256>;

  check_zip_bbs_3d<ExecSpace>();
}

TEST(primal_zip, zip_bbs_2d_from_3d_hip)
{
  using ExecSpace = axom::HIP_EXEC<256>;

  check_zip_bbs_2d_from_3d<ExecSpace>();
}

TEST(primal_zip, zip_rays_3d_hip)
{
  using ExecSpace = axom::HIP_EXEC<256>;

  check_zip_rays_3d<ExecSpace>();
}

TEST(primal_zip, zip_rays_2d_from_3d_hip)
{
  using ExecSpace = axom::HIP_EXEC<256>;

  check_zip_rays_2d_from_3d<ExecSpace>();
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
