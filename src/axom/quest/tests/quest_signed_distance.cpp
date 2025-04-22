// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"

// _quest_distance_cpp_include_start
#include "axom/quest/SignedDistance.hpp"
// _quest_distance_cpp_include_end

#include "quest_test_utilities.hpp"

// Google Test includes
#include "gtest/gtest.h"

// C/C++ includes
#include <cmath>

// Aliases
namespace mint = axom::mint;
namespace quest = axom::quest;
namespace primal = axom::primal;
using UMesh = axom::mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

// Note: the following macro is intended for debugging purposes. Uncomment the
// following line to enable VTK dumps from this test.
// #define QUEST_SIGNED_DISTANCE_TEST_DUMP_VTK

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
/*!
 * \brief Returns the bounding box of the mesh.
 * \param [in] mesh pointer to the mesh instance.
 * \return bb bounding box of the mesh
 */
primal::BoundingBox<double, 3> getBounds(const axom::mint::Mesh* mesh)
{
  SLIC_ASSERT(mesh != nullptr);

  primal::BoundingBox<double, 3> bb;
  primal::Point<double, 3> pt;

  const int nnodes = mesh->getNumberOfNodes();
  for(int inode = 0; inode < nnodes; ++inode)
  {
    mesh->getNode(inode, pt.data());
    bb.addPoint(pt);
  }

  return (bb);
}

/*!
 * \brief Generates a uniform mesh surrounding the given triangle mesh.
 * \param [in] mesh pointer to the input mesh.
 * \param [in] umesh pointer to the uniform mesh;
 */
void getUniformMesh(const UMesh* mesh, mint::UniformMesh*& umesh)
{
  SLIC_ASSERT(mesh != nullptr);
  SLIC_ASSERT(umesh == nullptr);

  constexpr int N = 16;  // number of points along each dimension

  primal::BoundingBox<double, 3> bb = getBounds(mesh);
  bb.expand(2.0);

  const double* lo = bb.getMin().data();
  const double* hi = bb.getMax().data();

  // construct an N x N x N grid
  umesh = new mint::UniformMesh(lo, hi, N, N, N);
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
TEST(quest_signed_distance, sphere_test)
{
  constexpr double l1norm_expected = 6.7051997372579715;
  constexpr double l2norm_expected = 2.5894400431865519;
  constexpr double linf_expected = 0.00532092;
  constexpr double TOL = 1.e-3;

  constexpr double SPHERE_RADIUS = 0.5;
  constexpr int SPHERE_THETA_RES = 25;
  constexpr int SPHERE_PHI_RES = 25;
  const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};

  primal::Sphere<double, 3> analytic_sphere(SPHERE_RADIUS);

  SLIC_INFO("Constructing sphere mesh...");
  UMesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
  quest::utilities::getSphereSurfaceMesh(surface_mesh,
                                         SPHERE_CENTER,
                                         SPHERE_RADIUS,
                                         SPHERE_THETA_RES,
                                         SPHERE_PHI_RES);

  SLIC_INFO("Generating uniform mesh...");
  mint::UniformMesh* umesh = nullptr;
  getUniformMesh(surface_mesh, umesh);

  double* phi_computed = umesh->createField<double>("phi_computed", mint::NODE_CENTERED);
  double* phi_expected = umesh->createField<double>("phi_expected", mint::NODE_CENTERED);
  double* phi_diff = umesh->createField<double>("phi_diff", mint::NODE_CENTERED);
  double* phi_err = umesh->createField<double>("phi_err", mint::NODE_CENTERED);

  const int nnodes = umesh->getNumberOfNodes();

  SLIC_INFO("Generate BVHTree...");
  constexpr bool is_watertight = true;
  constexpr bool compute_signs = true;
  axom::quest::SignedDistance<3> signed_distance(surface_mesh, is_watertight, compute_signs);

  SLIC_INFO("Compute signed distance...");

  double l1norm = 0.0;
  double l2norm = 0.0;
  double linf = axom::numeric_limits<double>::min();

  for(int inode = 0; inode < nnodes; ++inode)
  {
    primal::Point<double, 3> pt;
    umesh->getNode(inode, pt.data());

    phi_computed[inode] = signed_distance.computeDistance(pt);
    phi_expected[inode] = analytic_sphere.computeSignedDistance(pt);
    EXPECT_NEAR(phi_computed[inode], phi_expected[inode], 1.e-2);

    // compute error
    phi_diff[inode] = phi_computed[inode] - phi_expected[inode];
    phi_err[inode] = std::fabs(phi_diff[inode]);

    // update norms
    l1norm += phi_err[inode];
    l2norm += phi_diff[inode];
    linf = (phi_err[inode] > linf) ? phi_err[inode] : linf;

  }  // END for all nodes

  l2norm = std::sqrt(l2norm);

  SLIC_INFO("l1 = " << l1norm);
  SLIC_INFO("l2 = " << l2norm);
  SLIC_INFO("linf = " << linf);

#ifdef QUEST_SIGNED_DISTANCE_TEST_DUMP_VTK
  mint::write_vtk(umesh, "uniform_mesh.vtk");
  mint::write_vtk(surface_mesh, "sphere_mesh.vtk");
#endif

  EXPECT_NEAR(l1norm_expected, l1norm, TOL);
  EXPECT_NEAR(l2norm_expected, l2norm, TOL);
  EXPECT_NEAR(linf_expected, linf, TOL);

  delete surface_mesh;
  delete umesh;

  SLIC_INFO("Done.");
}

//------------------------------------------------------------------------------
TEST(quest_signed_distance, sphere_test_with_normals)
{
  using PointType = primal::Point<double, 3>;
  using VectorType = primal::Vector<double, 3>;
  using SphereType = primal::Sphere<double, 3>;

  constexpr double l1norm_expected = 6.7051997372579715;
  constexpr double l2norm_expected = 2.5894400431865519;
  constexpr double linf_expected = 0.00532092;
  constexpr double TOL = 1.e-3;

  constexpr double SPHERE_RADIUS = 0.5;
  constexpr int SPHERE_THETA_RES = 25;
  constexpr int SPHERE_PHI_RES = 25;
  const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};

  SphereType analytic_sphere(SPHERE_RADIUS);

  SLIC_INFO("Constructing sphere mesh...");
  UMesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
  quest::utilities::getSphereSurfaceMesh(surface_mesh,
                                         SPHERE_CENTER,
                                         SPHERE_RADIUS,
                                         SPHERE_THETA_RES,
                                         SPHERE_PHI_RES);

  SLIC_INFO("Generating uniform mesh...");
  mint::UniformMesh* umesh = nullptr;
  getUniformMesh(surface_mesh, umesh);

  double* phi_computed = umesh->createField<double>("phi_computed", mint::NODE_CENTERED);
  double* phi_expected = umesh->createField<double>("phi_expected", mint::NODE_CENTERED);
  double* phi_diff = umesh->createField<double>("phi_diff", mint::NODE_CENTERED);
  double* phi_err = umesh->createField<double>("phi_err", mint::NODE_CENTERED);

  double* cp_computed_x = umesh->createField<double>("cp_computed_x", mint::NODE_CENTERED);
  double* cp_computed_y = umesh->createField<double>("cp_computed_y", mint::NODE_CENTERED);
  double* cp_computed_z = umesh->createField<double>("cp_computed_z", mint::NODE_CENTERED);

  double* normal_computed_x = umesh->createField<double>("normal_computed_x", mint::NODE_CENTERED);
  double* normal_computed_y = umesh->createField<double>("normal_computed_y", mint::NODE_CENTERED);
  double* normal_computed_z = umesh->createField<double>("normal_computed_z", mint::NODE_CENTERED);

  const int nnodes = umesh->getNumberOfNodes();

  SLIC_INFO("Generate BVHTree...");
  constexpr bool is_watertight = true;
  constexpr bool compute_signs = true;
  axom::quest::SignedDistance<3> signed_distance(surface_mesh, is_watertight, compute_signs);

  SLIC_INFO("Compute signed distance...");

  double l1norm = 0.0;
  double l2norm = 0.0;
  double linf = axom::numeric_limits<double>::min();

  for(int inode = 0; inode < nnodes; ++inode)
  {
    PointType pt;
    umesh->getNode(inode, pt.data());

    PointType cp;
    VectorType normal;

    phi_computed[inode] = signed_distance.computeDistance(pt, cp, normal);

    // Check that the computed closest point on the discretized sphere is close
    // to where it would be on an analytic sphere
    cp_computed_x[inode] = cp[0];
    cp_computed_y[inode] = cp[1];
    cp_computed_z[inode] = cp[2];
    const auto cp_expected = primal::closest_point(pt, analytic_sphere);
    EXPECT_NEAR(0., primal::squared_distance(cp_expected, cp), 1e-2);

    // Check that the computed (pseudo)-normal on the discretized sphere is close
    // to where it would be on an analytic sphere, i.e. the dot product
    // between the two should be close to 1
    normal_computed_x[inode] = normal[0];
    normal_computed_y[inode] = normal[1];
    normal_computed_z[inode] = normal[2];
    const auto normal_expected = VectorType(cp_expected).unitVector();
    EXPECT_NEAR(1., normal.dot(normal_expected), 1e-2);

    // Check that the distance (squared) is the same as the distance
    // between the query point and the returned closest point
    EXPECT_NEAR(phi_computed[inode] * phi_computed[inode], primal::squared_distance(pt, cp), 1e-8);

    phi_expected[inode] = analytic_sphere.computeSignedDistance(pt);
    EXPECT_NEAR(phi_computed[inode], phi_expected[inode], 1.e-2);

    // compute error
    phi_diff[inode] = phi_computed[inode] - phi_expected[inode];
    phi_err[inode] = std::fabs(phi_diff[inode]);

    // update norms
    l1norm += phi_err[inode];
    l2norm += phi_diff[inode];
    linf = (phi_err[inode] > linf) ? phi_err[inode] : linf;

  }  // END for all nodes

  l2norm = std::sqrt(l2norm);

  SLIC_INFO("l1 = " << l1norm);
  SLIC_INFO("l2 = " << l2norm);
  SLIC_INFO("linf = " << linf);

#ifdef QUEST_SIGNED_DISTANCE_TEST_DUMP_VTK
  mint::write_vtk(umesh, "uniform_mesh_with_normals.vtk");
  mint::write_vtk(surface_mesh, "sphere_mesh_with_normals.vtk");
#endif

  EXPECT_NEAR(l1norm_expected, l1norm, TOL);
  EXPECT_NEAR(l2norm_expected, l2norm, TOL);
  EXPECT_NEAR(linf_expected, linf, TOL);

  delete surface_mesh;
  delete umesh;

  SLIC_INFO("Done.");
}
//------------------------------------------------------------------------------
template <typename ExecSpace>
void run_vectorized_sphere_test()
{
  using PointType = primal::Point<double, 3>;

  int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  int kernel_allocator = axom::execution_space<ExecSpace>::allocatorID();

  //Use unified memory on device
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernel_allocator = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  }
#endif

  axom::setDefaultAllocator(kernel_allocator);

  constexpr double l1norm_expected = 6.7051997372579715;
  constexpr double l2norm_expected = 2.5894400431865519;
  constexpr double linf_expected = 0.00532092;
  constexpr double TOL = 1.e-3;

  constexpr double SPHERE_RADIUS = 0.5;
  constexpr int SPHERE_THETA_RES = 25;
  constexpr int SPHERE_PHI_RES = 25;
  const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};

  primal::Sphere<double, 3> analytic_sphere(SPHERE_RADIUS);

  SLIC_INFO("Constructing sphere mesh...");
  UMesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
  quest::utilities::getSphereSurfaceMesh(surface_mesh,
                                         SPHERE_CENTER,
                                         SPHERE_RADIUS,
                                         SPHERE_THETA_RES,
                                         SPHERE_PHI_RES);

  SLIC_INFO("Generating uniform mesh...");
  mint::UniformMesh* umesh = nullptr;
  getUniformMesh(surface_mesh, umesh);

  double* phi_computed = umesh->createField<double>("phi_computed", mint::NODE_CENTERED);
  double* phi_expected = umesh->createField<double>("phi_expected", mint::NODE_CENTERED);
  double* phi_diff = umesh->createField<double>("phi_diff", mint::NODE_CENTERED);
  double* phi_err = umesh->createField<double>("phi_err", mint::NODE_CENTERED);

  const int nnodes = umesh->getNumberOfNodes();

  SLIC_INFO("Generate BVHTree...");
  constexpr bool is_watertight = true;
  constexpr bool compute_signs = true;
  quest::SignedDistance<3, ExecSpace> signed_distance(surface_mesh, is_watertight, compute_signs);

  SLIC_INFO("Compute signed distance...");

  double l1norm = 0.0;
  double l2norm = 0.0;
  double linf = axom::numeric_limits<double>::min();

  axom::Array<PointType> queryPts = axom::Array<PointType>(nnodes, nnodes, host_allocator);
  for(int inode = 0; inode < nnodes; inode++)
  {
    umesh->getNode(inode, queryPts[inode].data());
  }

  // Copy query points to device
  axom::Array<PointType> queryPtsDevice = axom::Array<PointType>(queryPts, kernel_allocator);
  double* phi_computed_device = axom::allocate<double>(nnodes, kernel_allocator);

  signed_distance.computeDistances(nnodes, queryPtsDevice.view(), phi_computed_device);

  // Copy output to host
  axom::copy(phi_computed, phi_computed_device, nnodes * sizeof(double));

  for(int inode = 0; inode < nnodes; ++inode)
  {
    phi_expected[inode] = analytic_sphere.computeSignedDistance(queryPts[inode]);
    EXPECT_NEAR(phi_computed[inode], phi_expected[inode], 1.e-2);

    // compute error
    phi_diff[inode] = phi_computed[inode] - phi_expected[inode];
    phi_err[inode] = std::fabs(phi_diff[inode]);

    // update norms
    l1norm += phi_err[inode];
    l2norm += phi_diff[inode];
    linf = (phi_err[inode] > linf) ? phi_err[inode] : linf;

  }  // END for all nodes

  l2norm = std::sqrt(l2norm);

  SLIC_INFO("l1 = " << l1norm);
  SLIC_INFO("l2 = " << l2norm);
  SLIC_INFO("linf = " << linf);

#ifdef QUEST_SIGNED_DISTANCE_TEST_DUMP_VTK
  mint::write_vtk(umesh, "uniform_mesh.vtk");
  mint::write_vtk(surface_mesh, "sphere_mesh.vtk");
#endif

  EXPECT_NEAR(l1norm_expected, l1norm, TOL);
  EXPECT_NEAR(l2norm_expected, l2norm, TOL);
  EXPECT_NEAR(linf_expected, linf, TOL);

  axom::deallocate(phi_computed_device);

  delete surface_mesh;
  delete umesh;

  axom::setDefaultAllocator(host_allocator);

  SLIC_INFO("Done.");
}

//------------------------------------------------------------------------------
TEST(quest_signed_distance, sphere_vec_test) { run_vectorized_sphere_test<axom::SEQ_EXEC>(); }

//------------------------------------------------------------------------------
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
TEST(quest_signed_distance, sphere_vec_omp_test) { run_vectorized_sphere_test<axom::OMP_EXEC>(); }
#endif  // AXOM_USE_OPENMP

//------------------------------------------------------------------------------
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
TEST(quest_signed_distance, sphere_vec_device_test)
{
  constexpr int BLOCK_SIZE = 256;
  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  run_vectorized_sphere_test<exec>();
}

//------------------------------------------------------------------------------
TEST(quest_signed_distance, sphere_vec_device_custom_alloc)
{
  constexpr int BLOCK_SIZE = 256;

  #if defined(__CUDACC__)
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;
  #elif defined(__HIPCC__)
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;
  #else
  using exec = axom::SEQ_EXEC;
  #endif

  using PointType = primal::Point<double, 3>;

  const int host_allocator = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  constexpr bool on_device = axom::execution_space<exec>::onDevice();
  const int kernel_allocator = on_device
    ? axom::getUmpireResourceAllocatorID(umpire::resource::Unified)
    : axom::execution_space<exec>::allocatorID();
  axom::setDefaultAllocator(kernel_allocator);

  constexpr double l1norm_expected = 6.7051997372579715;
  constexpr double l2norm_expected = 2.5894400431865519;
  constexpr double linf_expected = 0.00532092;
  constexpr double TOL = 1.e-3;

  constexpr double SPHERE_RADIUS = 0.5;
  constexpr int SPHERE_THETA_RES = 25;
  constexpr int SPHERE_PHI_RES = 25;
  const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};

  primal::Sphere<double, 3> analytic_sphere(SPHERE_RADIUS);

  SLIC_INFO("Constructing sphere mesh...");
  UMesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
  quest::utilities::getSphereSurfaceMesh(surface_mesh,
                                         SPHERE_CENTER,
                                         SPHERE_RADIUS,
                                         SPHERE_THETA_RES,
                                         SPHERE_PHI_RES);

  SLIC_INFO("Generating uniform mesh...");
  mint::UniformMesh* umesh = nullptr;
  getUniformMesh(surface_mesh, umesh);

  double* phi_computed = umesh->createField<double>("phi_computed", mint::NODE_CENTERED);
  double* phi_expected = umesh->createField<double>("phi_expected", mint::NODE_CENTERED);
  double* phi_diff = umesh->createField<double>("phi_diff", mint::NODE_CENTERED);
  double* phi_err = umesh->createField<double>("phi_err", mint::NODE_CENTERED);

  const int nnodes = umesh->getNumberOfNodes();

  SLIC_INFO("Generate BVHTree...");
  // _quest_distance_cpp_init_start
  // Set execution space
  constexpr int BlockSize = 256;
  #if defined(__CUDACC__)
  using ExecSpace = axom::CUDA_EXEC<BlockSize>;
  #elif defined(__HIPCC__)
  using ExecSpace = axom::HIP_EXEC<BlockSize>;
  #else
  using ExecSpace = axom::SEQ_EXEC;
  #endif

  // Create a custom allocator
  constexpr size_t PoolSize = 1024 * 1024 * 1024;
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator device_allocator =
    rm.makeAllocator<umpire::strategy::QuickPool>("DEVICE_POOL",
                                                  rm.getAllocator(umpire::resource::Device),
                                                  PoolSize);
  int device_pool_id = device_allocator.getId();

  // Set SignedDistance options
  constexpr bool is_watertight = true;
  constexpr bool compute_signs = true;
  quest::SignedDistance<3, ExecSpace> signed_distance(surface_mesh,
                                                      is_watertight,
                                                      compute_signs,
                                                      device_pool_id);
  // _quest_distance_cpp_init_end

  SLIC_INFO("Compute signed distance...");

  double l1norm = 0.0;
  double l2norm = 0.0;
  double linf = axom::numeric_limits<double>::min();

  axom::Array<PointType> queryPts = axom::Array<PointType>(nnodes, nnodes, host_allocator);
  for(int inode = 0; inode < nnodes; inode++)
  {
    umesh->getNode(inode, queryPts[inode].data());
  }

  // Copy query points to device
  axom::Array<PointType> queryPtsDevice(queryPts, kernel_allocator);
  double* phi_computed_device = axom::allocate<double>(nnodes, kernel_allocator);

  signed_distance.computeDistances(nnodes, queryPtsDevice.view(), phi_computed_device);

  // Copy output to host
  axom::copy(phi_computed, phi_computed_device, nnodes * sizeof(double));

  for(int inode = 0; inode < nnodes; ++inode)
  {
    phi_expected[inode] = analytic_sphere.computeSignedDistance(queryPts[inode]);
    EXPECT_NEAR(phi_computed[inode], phi_expected[inode], 1.e-2);

    // compute error
    phi_diff[inode] = phi_computed[inode] - phi_expected[inode];
    phi_err[inode] = std::fabs(phi_diff[inode]);

    // update norms
    l1norm += phi_err[inode];
    l2norm += phi_diff[inode];
    linf = (phi_err[inode] > linf) ? phi_err[inode] : linf;

  }  // END for all nodes

  l2norm = std::sqrt(l2norm);

  SLIC_INFO("l1 = " << l1norm);
  SLIC_INFO("l2 = " << l2norm);
  SLIC_INFO("linf = " << linf);

  #ifdef QUEST_SIGNED_DISTANCE_TEST_DUMP_VTK
  mint::write_vtk(umesh, "uniform_mesh.vtk");
  mint::write_vtk(surface_mesh, "sphere_mesh.vtk");
  #endif

  EXPECT_NEAR(l1norm_expected, l1norm, TOL);
  EXPECT_NEAR(l2norm_expected, l2norm, TOL);
  EXPECT_NEAR(linf_expected, linf, TOL);

  axom::deallocate(phi_computed_device);

  delete surface_mesh;
  delete umesh;

  axom::setDefaultAllocator(host_allocator);
  SLIC_INFO("Done.");
}
#endif  // defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
