// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"                          // compile-time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits
#include "axom/core/memory_management.hpp"          // for alloc() /free()

// Mint includes
#include "axom/mint/config.hpp"               // mint compile-time definitions
#include "axom/mint/execution/interface.hpp"  // for_all()

// Slic includes
#include "axom/slic.hpp"  // for SLIC macros

#include "mint_test_utilities.hpp"

// gtest includes
#include "gtest/gtest.h"  // for gtest

namespace axom
{
namespace mint
{
//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_nodes_idx(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  constexpr int MAGIC_VAL = 42;

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? 20 : -1;
  const IndexType Nk = (dimension == 3) ? 20 : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numNodes = test_mesh->getNumberOfNodes();
  axom::Array<int> field_d(numNodes, numNodes, device_allocator);

  auto field_v = field_d.view();

  for_all_nodes<ExecPolicy>(
    test_mesh,
    AXOM_LAMBDA(IndexType nodeIdx) { field_v[nodeIdx] = MAGIC_VAL; });

  // Copy data back to host
  axom::Array<int> field_h = axom::Array<int>(field_d, host_allocator);

  // Create mesh field from buffer
  int* n1_field = test_mesh->template createField<int>("n1", NODE_CENTERED, field_h.data());

  for(IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_EQ(n1_field[inode], MAGIC_VAL);
  }

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType>
void check_for_all_nodes_ij()
{
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name()
                      << ", mesh_type=" << internal::mesh_type<MeshType>::name());

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  constexpr IndexType N = 20;
  const double lo[] = {-10, -10};
  const double hi[] = {10, 10};
  UniformMesh uniform_mesh(lo, hi, N, N);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numNodes = test_mesh->getNumberOfNodes();

  axom::Array<IndexType> icoords_d(numNodes, numNodes, device_allocator);
  axom::Array<IndexType> jcoords_d(numNodes, numNodes, device_allocator);

  auto icoords_v = icoords_d.view();
  auto jcoords_v = jcoords_d.view();

  for_all_nodes<ExecPolicy, xargs::ij>(
    test_mesh,
    AXOM_LAMBDA(IndexType nodeIdx, IndexType i, IndexType j) {
      icoords_v[nodeIdx] = i;
      jcoords_v[nodeIdx] = j;
    });

  // Copy data back to host
  axom::Array<IndexType> icoords_h = axom::Array<IndexType>(icoords_d, host_allocator);
  axom::Array<IndexType> jcoords_h = axom::Array<IndexType>(jcoords_d, host_allocator);

  IndexType inode = 0;
  for(IndexType j = 0; j < N; ++j)
  {
    for(IndexType i = 0; i < N; ++i)
    {
      EXPECT_EQ(icoords_h[inode], i);
      EXPECT_EQ(jcoords_h[inode], j);
      ++inode;

    }  // END for all i
  }  // END for all j

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType>
void check_for_all_nodes_ijk()
{
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name()
                      << ", mesh_type=" << internal::mesh_type<MeshType>::name());

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  constexpr IndexType N = 20;
  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, N, N, N);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numNodes = test_mesh->getNumberOfNodes();

  axom::Array<IndexType> icoords_d(numNodes, numNodes, device_allocator);
  axom::Array<IndexType> jcoords_d(numNodes, numNodes, device_allocator);
  axom::Array<IndexType> kcoords_d(numNodes, numNodes, device_allocator);

  auto icoords_v = icoords_d.view();
  auto jcoords_v = jcoords_d.view();
  auto kcoords_v = kcoords_d.view();

  for_all_nodes<ExecPolicy, xargs::ijk>(
    test_mesh,
    AXOM_LAMBDA(IndexType nodeIdx, IndexType i, IndexType j, IndexType k) {
      icoords_v[nodeIdx] = i;
      jcoords_v[nodeIdx] = j;
      kcoords_v[nodeIdx] = k;
    });

  // Copy data back to host
  axom::Array<IndexType> icoords_h = axom::Array<IndexType>(icoords_d, host_allocator);
  axom::Array<IndexType> jcoords_h = axom::Array<IndexType>(jcoords_d, host_allocator);
  axom::Array<IndexType> kcoords_h = axom::Array<IndexType>(kcoords_d, host_allocator);

  IndexType inode = 0;
  for(IndexType k = 0; k < N; ++k)
  {
    for(IndexType j = 0; j < N; ++j)
    {
      for(IndexType i = 0; i < N; ++i)
      {
        EXPECT_EQ(icoords_h[inode], i);
        EXPECT_EQ(jcoords_h[inode], j);
        EXPECT_EQ(kcoords_h[inode], k);
        ++inode;

      }  // END for all i
    }  // END for all j
  }  // END for all k

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_nodes_xyz()
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name() << ", mesh_type=" << mesh_name);

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  constexpr IndexType N = 20;
  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, N, N, N);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  // STEP 1: generate test coordinate arrays
  const IndexType numNodes = test_mesh->getNumberOfNodes();

  axom::Array<double> x_d(numNodes, numNodes, device_allocator);
  axom::Array<double> y_d(numNodes, numNodes, device_allocator);
  axom::Array<double> z_d(numNodes, numNodes, device_allocator);

  auto x_v = x_d.view();
  auto y_v = y_d.view();
  auto z_v = z_d.view();

  for_all_nodes<ExecPolicy, xargs::xyz>(
    test_mesh,
    AXOM_LAMBDA(IndexType idx, double xx, double yy, double zz) {
      x_v[idx] = xx;
      y_v[idx] = yy;
      z_v[idx] = zz;
    });

  // Copy data back to host
  axom::Array<double> x_h = axom::Array<double>(x_d, host_allocator);
  axom::Array<double> y_h = axom::Array<double>(y_d, host_allocator);
  axom::Array<double> z_h = axom::Array<double>(z_d, host_allocator);

  // STEP 2:check coordinate arrays
  for(int inode = 0; inode < numNodes; ++inode)
  {
    double node[3];
    test_mesh->getNode(inode, node);
    EXPECT_DOUBLE_EQ(x_h[inode], node[X_COORDINATE]);
    EXPECT_DOUBLE_EQ(y_h[inode], node[Y_COORDINATE]);
    EXPECT_DOUBLE_EQ(z_h[inode], node[Z_COORDINATE]);
  }  // END for all nodes

  // STEP 3: clean up
  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_nodes_xy()
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name() << ", mesh_type=" << mesh_name);

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  constexpr IndexType N = 20;
  const double lo[] = {-10, -10};
  const double hi[] = {10, 10};
  UniformMesh uniform_mesh(lo, hi, N, N);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  // STEP 1: generate test coordinate arrays
  const IndexType numNodes = test_mesh->getNumberOfNodes();

  axom::Array<double> x_d(numNodes, numNodes, device_allocator);
  axom::Array<double> y_d(numNodes, numNodes, device_allocator);

  auto x_v = x_d.view();
  auto y_v = y_d.view();

  for_all_nodes<ExecPolicy, xargs::xy>(
    test_mesh,
    AXOM_LAMBDA(IndexType idx, double xx, double yy) {
      x_v[idx] = xx;
      y_v[idx] = yy;
    });

  // Copy data back to host
  axom::Array<double> x_h = axom::Array<double>(x_d, host_allocator);
  axom::Array<double> y_h = axom::Array<double>(y_d, host_allocator);

  // STEP 2:check coordinate arrays
  for(int inode = 0; inode < numNodes; ++inode)
  {
    double node[2];
    test_mesh->getNode(inode, node);
    EXPECT_DOUBLE_EQ(x_h[inode], node[X_COORDINATE]);
    EXPECT_DOUBLE_EQ(y_h[inode], node[Y_COORDINATE]);
  }  // END for all nodes

  // STEP 3: clean up
  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_nodes_x()
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name() << ", mesh_type=" << mesh_name);

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  constexpr IndexType Ni = 20;
  const double lo[] = {-10};
  const double hi[] = {10};
  UniformMesh uniform_mesh(lo, hi, Ni);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  if(MeshType != PARTICLE_MESH)
  {
    EXPECT_EQ(test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells());
  }

  // STEP 1: generate test coordinate arrays
  const IndexType numNodes = uniform_mesh.getNumberOfNodes();

  axom::Array<double> x_d(numNodes, numNodes, device_allocator);
  auto x_v = x_d.view();

  for_all_nodes<ExecPolicy, xargs::x>(
    test_mesh,
    AXOM_LAMBDA(IndexType idx, double xx) { x_v[idx] = xx; });

  // Copy data back to host
  axom::Array<double> x_h = axom::Array<double>(x_d, host_allocator);

  // STEP 2:check coordinate arrays
  for(int inode = 0; inode < numNodes; ++inode)
  {
    double node[1];
    uniform_mesh.getNode(inode, node);
    EXPECT_DOUBLE_EQ(x_h[inode], node[X_COORDINATE]);
  }  // END for all nodes

  // STEP 3: clean up
  delete test_mesh;
  test_mesh = nullptr;
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_xyz)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_xyz<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xyz<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xyz<seq_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xyz<seq_exec, PARTICLE_MESH>();
  check_for_all_nodes_xyz<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xyz<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_xyz<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xyz<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xyz<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xyz<openmp_exec, PARTICLE_MESH>();
  check_for_all_nodes_xyz<openmp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xyz<openmp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_for_all_nodes_xyz<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xyz<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xyz<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xyz<cuda_exec, PARTICLE_MESH>();
  check_for_all_nodes_xyz<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xyz<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_for_all_nodes_xyz<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xyz<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xyz<hip_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xyz<hip_exec, PARTICLE_MESH>();
  check_for_all_nodes_xyz<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xyz<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_xy)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_xy<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xy<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xy<seq_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xy<seq_exec, PARTICLE_MESH>();
  check_for_all_nodes_xy<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xy<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_xy<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xy<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xy<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xy<openmp_exec, PARTICLE_MESH>();
  check_for_all_nodes_xy<openmp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xy<openmp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_for_all_nodes_xy<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xy<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xy<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xy<cuda_exec, PARTICLE_MESH>();
  check_for_all_nodes_xy<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xy<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_for_all_nodes_xy<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xy<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xy<hip_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xy<hip_exec, PARTICLE_MESH>();
  check_for_all_nodes_xy<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xy<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_x)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_x<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_x<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_x<seq_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_x<seq_exec, PARTICLE_MESH>();
  check_for_all_nodes_x<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_x<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_x<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_x<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_x<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_x<openmp_exec, PARTICLE_MESH>();
  check_for_all_nodes_x<openmp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_x<openmp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_for_all_nodes_x<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_x<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_x<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_x<cuda_exec, PARTICLE_MESH>();
  check_for_all_nodes_x<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_x<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_for_all_nodes_x<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_x<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_x<hip_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_x<hip_exec, PARTICLE_MESH>();
  check_for_all_nodes_x<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_x<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_ijk)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_ijk<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ijk<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ijk<seq_exec, STRUCTURED_RECTILINEAR_MESH>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_ijk<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ijk<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ijk<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_for_all_nodes_ijk<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ijk<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ijk<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_for_all_nodes_ijk<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ijk<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ijk<hip_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_ij)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_ij<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ij<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ij<seq_exec, STRUCTURED_RECTILINEAR_MESH>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_ij<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ij<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ij<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_for_all_nodes_ij<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ij<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ij<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_for_all_nodes_ij<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ij<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ij<hip_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_index)
{
  constexpr int NDIMS = 3;
  for(int i = 1; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_nodes_idx<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_nodes_idx<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_nodes_idx<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_nodes_idx<seq_exec, PARTICLE_MESH>(i);
    check_for_all_nodes_idx<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_nodes_idx<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_nodes_idx<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_nodes_idx<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_nodes_idx<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_nodes_idx<omp_exec, PARTICLE_MESH>(i);
    check_for_all_nodes_idx<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_nodes_idx<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    check_for_all_nodes_idx<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_nodes_idx<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_nodes_idx<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_nodes_idx<cuda_exec, PARTICLE_MESH>(i);
    check_for_all_nodes_idx<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_nodes_idx<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    check_for_all_nodes_idx<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_nodes_idx<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_nodes_idx<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_nodes_idx<hip_exec, PARTICLE_MESH>(i);
    check_for_all_nodes_idx<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_nodes_idx<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

  }  // END for all dimensions
}

} /* namespace mint */
} /* namespace axom */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
