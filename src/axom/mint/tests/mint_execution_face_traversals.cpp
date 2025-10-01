// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"                          // compile-time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits

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
void check_for_all_faces(int dimension)
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

  const double lo[] = {-10, -9, -8};
  const double hi[] = {10, 9, 8};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numFaces = test_mesh->getNumberOfFaces();

  axom::Array<IndexType> field_d(numFaces, numFaces, device_allocator);

  auto field_v = field_d.view();

  for_all_faces<ExecPolicy>(test_mesh, AXOM_LAMBDA(IndexType faceID) { field_v[faceID] = faceID; });

  // Copy field back to host
  axom::Array<IndexType> field_h = axom::Array<IndexType>(field_d, host_allocator);

  // Create mesh field from buffer
  IndexType* f1_field =
    test_mesh->template createField<IndexType>("f1", FACE_CENTERED, field_h.data());

  for(IndexType faceID = 0; faceID < numFaces; ++faceID)
  {
    EXPECT_EQ(f1_field[faceID], faceID);
  }

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_face_nodes(int dimension)
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

  const double lo[] = {-10, -9, -8};
  const double hi[] = {10, 9, 8};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numFaces = test_mesh->getNumberOfFaces();

  axom::Array<IndexType> conn_d(numFaces * MAX_FACE_NODES, numFaces * MAX_FACE_NODES, device_allocator);

  auto conn_v = conn_d.view();

  for_all_faces<ExecPolicy, xargs::nodeids>(
    test_mesh,
    AXOM_LAMBDA(IndexType faceID, const IndexType* nodes, IndexType N) {
      for(int i = 0; i < N; ++i)
      {
        conn_v[faceID * MAX_FACE_NODES + i] = nodes[i];
      }  // END for all face nodes
    });

  // Copy field back to host
  axom::Array<IndexType> conn_h = axom::Array<IndexType>(conn_d, host_allocator);

  // Create mesh field from buffer
  IndexType* conn_field =
    test_mesh->template createField<IndexType>("f1", FACE_CENTERED, conn_h.data());

  IndexType faceNodes[MAX_FACE_NODES];
  for(IndexType faceID = 0; faceID < numFaces; ++faceID)
  {
    const IndexType N = test_mesh->getFaceNodeIDs(faceID, faceNodes);
    for(int i = 0; i < N; ++i)
    {
      EXPECT_EQ(conn_field[faceID * MAX_FACE_NODES + i], faceNodes[i]);
    }
  }  // END for all cells

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_face_coords(int dimension)
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

  const double lo[] = {-10, -9, -8};
  const double hi[] = {10, 9, 8};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numFaces = test_mesh->getNumberOfFaces();
  axom::Array<IndexType> conn_d(numFaces * MAX_FACE_NODES, numFaces * MAX_FACE_NODES, device_allocator);
  axom::Array<double> coords_d(numFaces * dimension * MAX_FACE_NODES,
                               numFaces * dimension * MAX_FACE_NODES,
                               device_allocator);

  auto conn_v = conn_d.view();
  auto coords_v = coords_d.view();

  for_all_faces<ExecPolicy, xargs::coords>(
    test_mesh,
    AXOM_LAMBDA(IndexType faceID, const numerics::Matrix<double>& coordsMatrix, const IndexType* nodes) {
      const IndexType numNodes = coordsMatrix.getNumColumns();
      for(int i = 0; i < numNodes; ++i)
      {
        conn_v[faceID * MAX_FACE_NODES + i] = nodes[i];

        for(int dim = 0; dim < dimension; ++dim)
        {
          coords_v[faceID * dimension * MAX_FACE_NODES + i * dimension + dim] = coordsMatrix(dim, i);
        }
      }  // END for all face nodes
    });

  // Copy data back to host
  axom::Array<IndexType> conn_h = axom::Array<IndexType>(conn_d, host_allocator);
  axom::Array<double> coords_h = axom::Array<double>(coords_d, host_allocator);

  // Create mesh fields from buffers
  IndexType* conn_field =
    test_mesh->template createField<IndexType>("conn", FACE_CENTERED, conn_h.data(), MAX_FACE_NODES);
  double* coords_field = test_mesh->template createField<double>("coords",
                                                                 FACE_CENTERED,
                                                                 coords_h.data(),
                                                                 dimension * MAX_FACE_NODES);

  double nodeCoords[3];
  IndexType faceNodes[MAX_FACE_NODES];
  for(IndexType faceID = 0; faceID < numFaces; ++faceID)
  {
    const IndexType numNodes = test_mesh->getFaceNodeIDs(faceID, faceNodes);
    for(int i = 0; i < numNodes; ++i)
    {
      EXPECT_EQ(conn_field[faceID * MAX_FACE_NODES + i], faceNodes[i]);

      for(int dim = 0; dim < dimension; ++dim)
      {
        test_mesh->getNode(faceNodes[i], nodeCoords);
        EXPECT_NEAR(coords_field[faceID * dimension * MAX_FACE_NODES + i * dimension + dim],
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
void check_for_all_face_cells(int dimension)
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

  const double lo[] = {-10, -9, -8};
  const double hi[] = {10, 9, 8};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh = dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numFaces = test_mesh->getNumberOfFaces();
  axom::Array<IndexType> face_cells_d(numFaces * 2, numFaces * 2, device_allocator);

  auto face_cells_v = face_cells_d.view();

  for_all_faces<ExecPolicy, xargs::cellids>(
    test_mesh,
    AXOM_LAMBDA(IndexType faceID, IndexType cellIDOne, IndexType cellIDTwo) {
      face_cells_v[2 * faceID + 0] = cellIDOne;
      face_cells_v[2 * faceID + 1] = cellIDTwo;
    });

  // Copy field back to host
  axom::Array<IndexType> face_cells_h = axom::Array<IndexType>(face_cells_d, host_allocator);

  // Create mesh field from buffer
  IndexType* face_cells_field =
    test_mesh->template createField<IndexType>("f1", FACE_CENTERED, face_cells_h.data(), 2);

  for(IndexType faceID = 0; faceID < numFaces; ++faceID)
  {
    IndexType cellIDOne, cellIDTwo;
    test_mesh->getFaceCellIDs(faceID, cellIDOne, cellIDTwo);

    EXPECT_EQ(face_cells_field[2 * faceID + 0], cellIDOne);
    EXPECT_EQ(face_cells_field[2 * faceID + 1], cellIDTwo);
  }

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

AXOM_CUDA_TEST(mint_execution_face_traversals, for_all_face_nodeids)
{
  for(int dim = 2; dim <= 3; ++dim)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_face_nodes<seq_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_nodes<seq_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_nodes<seq_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_nodes<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_nodes<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_face_nodes<omp_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_nodes<omp_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_nodes<omp_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_nodes<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_nodes<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    check_for_all_face_nodes<cuda_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_nodes<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_nodes<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_nodes<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_nodes<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    check_for_all_face_nodes<hip_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_nodes<hip_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_nodes<hip_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_nodes<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_nodes<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

  }  // END for all dimensions
}

AXOM_CUDA_TEST(mint_execution_face_traversals, for_all_face_coords)
{
  for(int dim = 2; dim <= 3; ++dim)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_face_coords<seq_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_coords<seq_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_coords<seq_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_coords<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_coords<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_face_coords<omp_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_coords<omp_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_coords<omp_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_coords<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_coords<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    check_for_all_face_coords<cuda_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_coords<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_coords<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_coords<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_coords<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    check_for_all_face_coords<hip_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_coords<hip_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_coords<hip_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_coords<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_coords<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

  }  // END for all dimensions
}

AXOM_CUDA_TEST(mint_execution_face_traversals, for_all_face_cellids)
{
  for(int dim = 2; dim <= 3; ++dim)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_face_cells<seq_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_cells<seq_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_cells<seq_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_nodes<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_nodes<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_face_cells<omp_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_cells<omp_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_cells<omp_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_nodes<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_nodes<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    check_for_all_face_cells<cuda_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_cells<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_cells<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_nodes<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_nodes<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    check_for_all_face_cells<hip_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_face_cells<hip_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_face_cells<hip_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_face_nodes<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_face_nodes<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

  }  // END for all dimensions
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_face_traversals, for_all_faces_index)
{
  for(int dim = 2; dim <= 3; ++dim)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_faces<seq_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_faces<seq_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_faces<seq_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_faces<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_faces<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_faces<omp_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_faces<omp_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_faces<omp_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_faces<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_faces<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    check_for_all_faces<cuda_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_faces<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_faces<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_faces<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_faces<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    check_for_all_faces<hip_exec, STRUCTURED_UNIFORM_MESH>(dim);
    check_for_all_faces<hip_exec, STRUCTURED_CURVILINEAR_MESH>(dim);
    check_for_all_faces<hip_exec, STRUCTURED_RECTILINEAR_MESH>(dim);
    check_for_all_faces<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(dim);
    check_for_all_faces<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(dim);

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
