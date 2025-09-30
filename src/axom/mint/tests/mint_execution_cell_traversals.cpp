// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/numerics/Matrix.hpp"

// Mint includes
#include "axom/mint/config.hpp"
#include "axom/mint/execution/interface.hpp"

// Slic includes
#include "axom/slic.hpp"

#include "mint_test_utilities.hpp"

// gtest includes
#include "gtest/gtest.h"

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
void check_for_all_cells_idx(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numCells = test_mesh->getNumberOfCells();

  axom::Array<IndexType> field_d(numCells, numCells, device_allocator);

  auto field_v = field_d.view();

  for_all_cells<ExecPolicy>(test_mesh, AXOM_LAMBDA(IndexType cellID) { field_v[cellID] = cellID; });

  // Copy field back to host
  axom::Array<IndexType> field_h = axom::Array<IndexType>(field_d, host_allocator);

  // Create mesh field from buffer
  IndexType* c1_field =
    test_mesh->template createField<IndexType>("c1", CELL_CENTERED, field_h.data());

  for(IndexType cellID = 0; cellID < numCells; ++cellID)
  {
    EXPECT_EQ(c1_field[cellID], cellID);
  }

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType>
void check_for_all_cells_ij()
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

  axom::Array<IndexType> icoords_d(N * N, N * N, device_allocator);
  axom::Array<IndexType> jcoords_d(N * N, N * N, device_allocator);

  auto icoords_v = icoords_d.view();
  auto jcoords_v = jcoords_d.view();

  for_all_cells<ExecPolicy, xargs::ij>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellIdx, IndexType i, IndexType j) {
      icoords_v[cellIdx] = i;
      jcoords_v[cellIdx] = j;
    });

  // Copy data back to host
  axom::Array<IndexType> icoords_h = axom::Array<IndexType>(icoords_d, host_allocator);
  axom::Array<IndexType> jcoords_h = axom::Array<IndexType>(jcoords_d, host_allocator);

  // Create mesh fields from buffers
  IndexType* i_field =
    test_mesh->template createField<IndexType>("i", CELL_CENTERED, icoords_h.data());
  IndexType* j_field =
    test_mesh->template createField<IndexType>("j", CELL_CENTERED, jcoords_h.data());

  IndexType icell = 0;
  for(IndexType j = 0; j < (N - 1); ++j)
  {
    for(IndexType i = 0; i < (N - 1); ++i)
    {
      EXPECT_EQ(i_field[icell], i);
      EXPECT_EQ(j_field[icell], j);
      ++icell;
    }  // END for all i
  }  // END for all j

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType>
void check_for_all_cells_ijk()
{
  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name()
                      << ", mesh_type=" << internal::mesh_type<MeshType>::name());

  constexpr IndexType N = 20;
  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, N, N, N);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const int size = N * N * N;
  axom::Array<IndexType> icoords_d(size, size, device_allocator);
  axom::Array<IndexType> jcoords_d(size, size, device_allocator);
  axom::Array<IndexType> kcoords_d(size, size, device_allocator);

  auto icoords_v = icoords_d.view();
  auto jcoords_v = jcoords_d.view();
  auto kcoords_v = kcoords_d.view();

  for_all_cells<ExecPolicy, xargs::ijk>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellIdx, IndexType i, IndexType j, IndexType k) {
      icoords_v[cellIdx] = i;
      jcoords_v[cellIdx] = j;
      kcoords_v[cellIdx] = k;
    });

  // Copy data back to host
  axom::Array<IndexType> icoords_h = axom::Array<IndexType>(icoords_d, host_allocator);
  axom::Array<IndexType> jcoords_h = axom::Array<IndexType>(jcoords_d, host_allocator);
  axom::Array<IndexType> kcoords_h = axom::Array<IndexType>(kcoords_d, host_allocator);

  // Create mesh fields from buffers
  IndexType* i_field =
    test_mesh->template createField<IndexType>("i", CELL_CENTERED, icoords_h.data());
  IndexType* j_field =
    test_mesh->template createField<IndexType>("j", CELL_CENTERED, jcoords_h.data());

  IndexType* k_field =
    test_mesh->template createField<IndexType>("k", CELL_CENTERED, kcoords_h.data());

  IndexType icell = 0;
  for(IndexType k = 0; k < (N - 1); ++k)
  {
    for(IndexType j = 0; j < (N - 1); ++j)
    {
      for(IndexType i = 0; i < (N - 1); ++i)
      {
        EXPECT_EQ(i_field[icell], i);
        EXPECT_EQ(j_field[icell], j);
        EXPECT_EQ(k_field[icell], k);
        ++icell;
      }  // END for all i
    }  // END for all j
  }  // END for all k

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_cell_nodes(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numCells = test_mesh->getNumberOfCells();

  axom::Array<IndexType> conn_d(numCells * MAX_CELL_NODES, numCells * MAX_CELL_NODES, device_allocator);

  auto conn_v = conn_d.view();

  for_all_cells<ExecPolicy, xargs::nodeids>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellID, const IndexType* nodes, IndexType N) {
      for(int i = 0; i < N; ++i)
      {
        conn_v[cellID * MAX_CELL_NODES + i] = nodes[i];
      }  // END for all cell nodes
    });

  // Copy data back to host
  axom::Array<IndexType> conn_h = axom::Array<IndexType>(conn_d, host_allocator);

  // Create mesh field from buffer
  IndexType* conn_field =
    test_mesh->template createField<IndexType>("conn", CELL_CENTERED, conn_h.data(), MAX_CELL_NODES);

  IndexType cellNodes[MAX_CELL_NODES];
  for(IndexType cellID = 0; cellID < numCells; ++cellID)
  {
    const IndexType N = test_mesh->getCellNodeIDs(cellID, cellNodes);
    for(int i = 0; i < N; ++i)
    {
      EXPECT_EQ(conn_field[cellID * MAX_CELL_NODES + i], cellNodes[i]);
    }
  }  // END for all cells

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_cell_coords(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numCells = test_mesh->getNumberOfCells();

  axom::Array<IndexType> conn_d(numCells * MAX_CELL_NODES, numCells * MAX_CELL_NODES, device_allocator);
  axom::Array<double> coords_d(numCells * dimension * MAX_CELL_NODES,
                               numCells * dimension * MAX_CELL_NODES,
                               device_allocator);

  auto conn_v = conn_d.view();
  auto coords_v = coords_d.view();

  for_all_cells<ExecPolicy, xargs::coords>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellID, const numerics::Matrix<double>& coordsMatrix, const IndexType* nodes) {
      const IndexType numNodes = coordsMatrix.getNumColumns();
      for(int i = 0; i < numNodes; ++i)
      {
        conn_v[cellID * MAX_CELL_NODES + i] = nodes[i];

        for(int dim = 0; dim < dimension; ++dim)
        {
          coords_v[cellID * dimension * MAX_CELL_NODES + i * dimension + dim] = coordsMatrix(dim, i);
        }
      }  // END for all cell nodes
    });

  // Copy data back to host
  axom::Array<IndexType> conn_h = axom::Array<IndexType>(conn_d, host_allocator);
  axom::Array<double> coords_h = axom::Array<double>(coords_d, host_allocator);

  // Create mesh fields from buffers
  IndexType* conn_field =
    test_mesh->template createField<IndexType>("conn", CELL_CENTERED, conn_h.data(), MAX_CELL_NODES);
  double* coords_field = test_mesh->template createField<double>("coords",
                                                                 CELL_CENTERED,
                                                                 coords_h.data(),
                                                                 dimension * MAX_CELL_NODES);

  double nodeCoords[3];
  IndexType cellNodes[MAX_CELL_NODES];
  for(IndexType cellID = 0; cellID < numCells; ++cellID)
  {
    const IndexType numNodes = test_mesh->getCellNodeIDs(cellID, cellNodes);
    for(int i = 0; i < numNodes; ++i)
    {
      EXPECT_EQ(conn_field[cellID * MAX_CELL_NODES + i], cellNodes[i]);

      for(int dim = 0; dim < dimension; ++dim)
      {
        test_mesh->getNode(cellNodes[i], nodeCoords);
        EXPECT_NEAR(coords_field[cellID * dimension * MAX_CELL_NODES + i * dimension + dim],
                    nodeCoords[dim],
                    1e-8);
      }
    }
  }  // END for all cells

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_cell_faces(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numCells = test_mesh->getNumberOfCells();
  axom::Array<IndexType> cell_faces_d(numCells * MAX_CELL_FACES,
                                      numCells * MAX_CELL_FACES,
                                      device_allocator);

  auto cell_faces_v = cell_faces_d.view();

  for_all_cells<ExecPolicy, xargs::faceids>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellID, const IndexType* faces, IndexType N) {
      for(int i = 0; i < N; ++i)
      {
        cell_faces_v[cellID * MAX_CELL_FACES + i] = faces[i];
      }
    });

  // Copy data back to host
  axom::Array<IndexType> cell_faces_h = axom::Array<IndexType>(cell_faces_d, host_allocator);

  // Create mesh fields from buffers
  IndexType* cell_faces_field = test_mesh->template createField<IndexType>("cellFaces",
                                                                           CELL_CENTERED,
                                                                           cell_faces_h.data(),
                                                                           MAX_CELL_NODES);

  IndexType faces[MAX_CELL_FACES];
  for(IndexType cellID = 0; cellID < numCells; ++cellID)
  {
    const IndexType N = test_mesh->getCellFaceIDs(cellID, faces);

    for(IndexType i = 0; i < N; ++i)
    {
      EXPECT_EQ(cell_faces_field[cellID * MAX_CELL_FACES + i], faces[i]);
    }
  }

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_nodeids)
{
  constexpr int NDIMS = 3;
  for(int i = 1; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_cell_nodes<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_nodes<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_nodes<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_nodes<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_nodes<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_cell_nodes<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_nodes<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_nodes<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_nodes<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_nodes<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    check_for_all_cell_nodes<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_nodes<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_nodes<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_nodes<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_nodes<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    check_for_all_cell_nodes<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_nodes<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_nodes<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_nodes<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_nodes<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

  }  // END for all dimensions
}

AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_coords)
{
  constexpr int NDIMS = 3;
  for(int i = 1; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_cell_coords<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_coords<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_coords<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_coords<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_coords<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_cell_coords<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_coords<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_coords<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_coords<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_coords<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    check_for_all_cell_coords<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_coords<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_coords<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_coords<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_coords<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    check_for_all_cell_coords<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_coords<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_coords<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_coords<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_coords<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

  }  // END for all dimensions
}

AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_faceids)
{
  constexpr int NDIMS = 3;
  for(int i = 2; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_cell_faces<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_faces<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_faces<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_cell_faces<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_faces<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_faces<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    check_for_all_cell_faces<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_faces<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_faces<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    check_for_all_cell_faces<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_faces<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_faces<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_faces<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_faces<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

  }  // END for all dimensions
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_ij)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_cells_ij<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ij<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ij<seq_exec, STRUCTURED_RECTILINEAR_MESH>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using omp_exec = axom::OMP_EXEC;
  check_for_all_cells_ij<omp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ij<omp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ij<omp_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_for_all_cells_ij<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ij<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ij<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_for_all_cells_ij<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ij<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ij<hip_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_ijk)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_cells_ijk<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ijk<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ijk<seq_exec, STRUCTURED_RECTILINEAR_MESH>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using omp_exec = axom::OMP_EXEC;
  check_for_all_cells_ijk<omp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ijk<omp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ijk<omp_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_for_all_cells_ijk<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ijk<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ijk<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_for_all_cells_ijk<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ijk<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ijk<hip_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_index)
{
  constexpr int NDIMS = 3;
  for(int i = 1; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_cells_idx<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cells_idx<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cells_idx<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cells_idx<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cells_idx<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_cells_idx<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cells_idx<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cells_idx<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cells_idx<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cells_idx<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    check_for_all_cells_idx<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cells_idx<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cells_idx<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cells_idx<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cells_idx<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    check_for_all_cells_idx<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cells_idx<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cells_idx<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cells_idx<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cells_idx<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);
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
