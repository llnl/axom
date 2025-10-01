// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/quest/interface/signed_distance.hpp"

// Axom includes
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/SignedDistance.hpp"

#include "axom/primal.hpp"
#include "axom/mint/mesh/Mesh.hpp"

#include "axom/slic/interface/slic.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

#include <string>

namespace axom
{
namespace quest
{
//------------------------------------------------------------------------------
// INTERNAL DATA STRUCTURES
//------------------------------------------------------------------------------
namespace
{
constexpr int INIT_FAILED = -1;
constexpr int INIT_SUCCESS = 0;

using ExecSeq = axom::SEQ_EXEC;
using SignedDistance3D = SignedDistance<3>;
using SignedDistance2D = SignedDistance<2>;

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
using ExecOMP = axom::OMP_EXEC;
using SignedDistance3DOMP = SignedDistance<3, ExecOMP>;
using SignedDistance2DOMP = SignedDistance<2, ExecOMP>;
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA)
using ExecGPU = axom::CUDA_EXEC<256>;
using SignedDistance3DGPU = SignedDistance<3, ExecGPU>;
using SignedDistance2DGPU = SignedDistance<2, ExecGPU>;
#endif

#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA)
using ExecGPU = axom::HIP_EXEC<256>;
using SignedDistance3DGPU = SignedDistance<3, ExecGPU>;
using SignedDistance2DGPU = SignedDistance<2, ExecGPU>;
#endif

/*!
 * \brief Holds the options for the SignedDistance query.
 */
static struct parameters_t
{
  int dimension; /*!< the dimension, 2 or 3 */

  bool verbose;              /*!< logger verbosity */
  bool is_closed_surface;    /*!< indicates if the input is a closed surface */
  bool use_shared_memory;    /*!< use Umpire shared memory for the surface mesh */
  bool compute_sign;         /*!< indicates if sign should be computed */
  int allocator_id;          /*!< the allocator ID to create BVH with (-1 for default) */
  SignedDistExec exec_space; /*!< indicates the execution space to run in */

  /*!
   * \brief Default Constructor. Sets default values for the parameters.
   */
  parameters_t()
    : dimension(3)
    , verbose(false)
    , is_closed_surface(true)
    , use_shared_memory(false)
    , compute_sign(true)
    , allocator_id(-1)
    , exec_space(SignedDistExec::CPU)
  { }

} Parameters;

// TODO: note the SignedDistance query is currently only supported in 3-D
static SignedDistance3D* s_query = nullptr;
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
static SignedDistance3DOMP* s_query_omp = nullptr;
#endif
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
static SignedDistance3DGPU* s_query_gpu = nullptr;
#endif
static mint::Mesh* s_surface_mesh = nullptr;
static bool s_must_delete_mesh = false;
static bool s_must_finalize_logger = false;
static bool s_logger_is_initialized = false;

#if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
static unsigned char* s_shared_mesh_buffer = nullptr;
#else
static std::string s_shared_memory_requirements(
  "Shared memory requires MPI and an Umpire library built with shared memory support");
#endif

}  // end anonymous namespace

//------------------------------------------------------------------------------
// SIGNED DISTANCE QUERY INTERFACE IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
int signed_distance_init(const std::string& file, MPI_Comm comm)
{
  // STEP 0: initialize logger
  internal::logger_init(s_logger_is_initialized, s_must_finalize_logger, Parameters.verbose, comm);

  SLIC_ASSERT(s_query == nullptr);

  if(Parameters.dimension != 3)
  {
    SLIC_WARNING("the SignedDistance Query is currently only supported in 3D");
    return INIT_FAILED;
  }

  // STEP 0: read the STL mesh
  int rc = INIT_FAILED;

#if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
  if(Parameters.use_shared_memory)
  {
    rc = internal::read_stl_mesh_shared(file, comm, s_shared_mesh_buffer, s_surface_mesh);
  }
  else
  {
    rc = internal::read_stl_mesh(file, s_surface_mesh, comm);
  }
#else
  SLIC_WARNING_IF(Parameters.use_shared_memory, s_shared_memory_requirements);

  rc = internal::read_stl_mesh(file, s_surface_mesh, comm);
#endif

  if(rc != 0)
  {
    SLIC_WARNING("reading mesh from [" << file << "] failed!");
    return INIT_FAILED;
  }

  // STEP 1: initialized the signed distance query
  s_must_delete_mesh = true;
  rc = signed_distance_init(s_surface_mesh, comm);
  return rc;
}

//------------------------------------------------------------------------------
int signed_distance_init(const mint::Mesh* m, MPI_Comm comm)
{
  internal::logger_init(s_logger_is_initialized, s_must_finalize_logger, Parameters.verbose, comm);

  SLIC_ERROR_IF(signed_distance_initialized(),
                "signed distance query has already been initialized!");
  SLIC_ERROR_IF(m->getDimension() != 3, "signed distance query currently only support 3-D meshes");
  SLIC_ERROR_IF(m->getMeshType() != mint::UNSTRUCTURED_MESH,
                "signed distance query currently only supports unstructured meshes");
  SLIC_ERROR_IF(m->hasMixedCellTypes() == true,
                "signed distance query does not support meshes with mixed shape topology");
  SLIC_ERROR_IF(m->getCellType() != mint::TRIANGLE,
                "signed distance currently only support 3D triangular surface meshes");

  if(s_surface_mesh != m)
  {
    SLIC_ASSERT(s_surface_mesh == nullptr);
    s_surface_mesh = const_cast<mint::Mesh*>(m);
    s_must_delete_mesh = false;
  }

  int allocatorID = Parameters.allocator_id;
  switch(Parameters.exec_space)
  {
  case SignedDistExec::CPU:
    if(allocatorID == -1)
    {
      allocatorID = axom::execution_space<ExecSeq>::allocatorID();
    }
    s_query = new SignedDistance3D(s_surface_mesh,
                                   Parameters.is_closed_surface,
                                   Parameters.compute_sign,
                                   allocatorID);
    break;
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  case SignedDistExec::OpenMP:
    if(allocatorID == -1)
    {
      allocatorID = axom::execution_space<ExecOMP>::allocatorID();
    }
    s_query_omp = new SignedDistance3DOMP(s_surface_mesh,
                                          Parameters.is_closed_surface,
                                          Parameters.compute_sign,
                                          allocatorID);
    break;
#endif
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
  case SignedDistExec::GPU:
    if(allocatorID == -1)
    {
      allocatorID = axom::execution_space<ExecGPU>::allocatorID();
    }
    s_query_gpu = new SignedDistance3DGPU(s_surface_mesh,
                                          Parameters.is_closed_surface,
                                          Parameters.compute_sign,
                                          allocatorID);
    break;
#endif
  default:
    SLIC_ERROR("Unsupported execution space");
    return INIT_FAILED;
  }

  return INIT_SUCCESS;
}

//------------------------------------------------------------------------------
bool signed_distance_initialized()
{
  switch(Parameters.exec_space)
  {
  case SignedDistExec::CPU:
    return (s_query != nullptr);
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  case SignedDistExec::OpenMP:
    return (s_query_omp != nullptr);
#endif
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
  case SignedDistExec::GPU:
    return (s_query_gpu != nullptr);
#endif
  default:
    SLIC_ERROR("Unsupported execution space");
    return false;
  }
}

//------------------------------------------------------------------------------
void signed_distance_get_mesh_bounds(double* lo, double* hi)
{
  SLIC_ERROR_IF(
    !signed_distance_initialized(),
    "signed distance query must be initialized prior to" << "calling get_mesh_bounds()");
  SLIC_ERROR_IF(lo == nullptr, "supplied buffer is null");
  SLIC_ERROR_IF(hi == nullptr, "supplied buffer is null");

  internal::compute_mesh_bounds(s_surface_mesh, lo, hi);
}

//------------------------------------------------------------------------------
void signed_distance_set_dimension(int dim)
{
  SLIC_ERROR_IF(dim != 3, "The signed distance query only support 3D");
  SLIC_ERROR_IF(signed_distance_initialized(),
                "signed distance query already initialized; setting option has no effect!");

  Parameters.dimension = dim;
}

//------------------------------------------------------------------------------
void signed_distance_set_closed_surface(bool status)
{
  SLIC_ERROR_IF(signed_distance_initialized(),
                "signed distance query already initialized; setting option has no effect!");

  Parameters.is_closed_surface = status;
}

//------------------------------------------------------------------------------
void signed_distance_set_compute_signs(bool computeSign)
{
  SLIC_ERROR_IF(signed_distance_initialized(),
                "signed distance query already initialized; setting option has no effect!");

  Parameters.compute_sign = computeSign;
}

//------------------------------------------------------------------------------
void signed_distance_set_allocator(int allocatorID)
{
  SLIC_ERROR_IF(signed_distance_initialized(),
                "signed distance query already initialized; setting option has no effect!");

  Parameters.allocator_id = allocatorID;
}

//------------------------------------------------------------------------------
void signed_distance_set_verbose(bool status)
{
  SLIC_ERROR_IF(signed_distance_initialized(),
                "signed distance query already initialized; setting option has no effect!");

  Parameters.verbose = status;
}

//------------------------------------------------------------------------------
void signed_distance_use_shared_memory(bool status)
{
  SLIC_ERROR_IF(signed_distance_initialized(),
                "signed distance query already initialized; setting option has no effect!");

  Parameters.use_shared_memory = status;

#if !defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
  SLIC_WARNING_IF(Parameters.use_shared_memory,
                  s_shared_memory_requirements + " Option is ignored!");
#endif
}

//------------------------------------------------------------------------------
void signed_distance_set_execution_space(SignedDistExec exec_space)
{
  SLIC_ERROR_IF(signed_distance_initialized(),
                "signed distance query already initialized; setting option has no effect!");

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  if(exec_space == SignedDistExec::OpenMP)
  {
    SLIC_ERROR("Signed distance query not compiled with OpenMP support");
  }
#endif

#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
  if(exec_space == SignedDistExec::GPU)
  {
    SLIC_ERROR("Signed distance query not compiled with GPU support");
  }
#endif

  Parameters.exec_space = exec_space;
}

//------------------------------------------------------------------------------
double signed_distance_evaluate(double x, double y, double z)
{
  SLIC_ERROR_IF(!signed_distance_initialized(),
                "signed distance query must be initialized prior to calling evaluate()!");

  double phi = 0.0;
  switch(Parameters.exec_space)
  {
  case SignedDistExec::CPU:
    phi = s_query->computeDistance(x, y, z);
    break;
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  case SignedDistExec::OpenMP:
    phi = s_query_omp->computeDistance(x, y, z);
    break;
#endif
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
  case SignedDistExec::GPU:
    phi = s_query_gpu->computeDistance(x, y, z);
    break;
#endif
  default:
    SLIC_ERROR("Unsupported execution space");
    break;
  }
  return (phi);
}

//------------------------------------------------------------------------------
double signed_distance_evaluate(double x,
                                double y,
                                double z,
                                double& cp_x,
                                double& cp_y,
                                double& cp_z,
                                double& n_x,
                                double& n_y,
                                double& n_z)
{
  SLIC_ERROR_IF(!signed_distance_initialized(),
                "signed distance query must be initialized prior to calling evaluate()!");

  SLIC_ERROR_IF(Parameters.dimension != 3,
                "This overload of signed_distance_evaluate is only available "
                "for 3D queries");

  using PointType = primal::Point<double, 3>;
  using VectorType = primal::Vector<double, 3>;

  const PointType query {x, y, z};
  PointType closest_pt;
  VectorType normal;

  double phi = 0.0;
  switch(Parameters.exec_space)
  {
  case SignedDistExec::CPU:
    phi = s_query->computeDistance(query, closest_pt, normal);
    break;
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  case SignedDistExec::OpenMP:
    phi = s_query_omp->computeDistance(query, closest_pt, normal);
    break;
#endif
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
  case SignedDistExec::GPU:
    phi = s_query_gpu->computeDistance(query, closest_pt, normal);
    break;
#endif
  default:
    SLIC_ERROR("Unsupported execution space");
    break;
  }

  // copy result back to the output parameters
  cp_x = closest_pt[0];
  cp_y = closest_pt[1];
  cp_z = closest_pt[2];
  n_x = normal[0];
  n_y = normal[1];
  n_z = normal[2];

  return (phi);
}

//------------------------------------------------------------------------------
void signed_distance_evaluate(const double* x, const double* y, const double* z, int npoints, double* phi)
{
  SLIC_ERROR_IF(!signed_distance_initialized(),
                "signed distance query must be initialized prior to calling evaluate()!");
  SLIC_ERROR_IF(x == nullptr, "x-coords array is null");
  SLIC_ERROR_IF(y == nullptr, "y-coords array is null");
  SLIC_ERROR_IF(z == nullptr, "z-coords array is null");
  SLIC_ERROR_IF(phi == nullptr, "output phi array is null");

  using PointType = primal::Point<double, 3>;
  using ZipPoint = primal::ZipIndexable<PointType>;

  ZipPoint it {{x, y, z}};

  switch(Parameters.exec_space)
  {
  case SignedDistExec::CPU:
    s_query->computeDistances(npoints, it, phi);
    break;
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  case SignedDistExec::OpenMP:
    s_query_omp->computeDistances(npoints, it, phi);
    break;
#endif
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
  case SignedDistExec::GPU:
    s_query_gpu->computeDistances(npoints, it, phi);
    break;
#endif
  default:
    SLIC_ERROR("Unsupported execution space");
    break;
  }
}

//------------------------------------------------------------------------------
void signed_distance_finalize()
{
  if(s_query != nullptr)
  {
    delete s_query;
    s_query = nullptr;
  }

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  if(s_query_omp != nullptr)
  {
    delete s_query_omp;
    s_query_omp = nullptr;
  }
#endif

#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA)
  if(s_query_gpu != nullptr)
  {
    delete s_query_gpu;
    s_query_gpu = nullptr;
  }
#endif

  if(s_surface_mesh != nullptr && s_must_delete_mesh)
  {
    delete s_surface_mesh;
  }

  s_surface_mesh = nullptr;

  SLIC_ASSERT(!signed_distance_initialized());
  internal::logger_finalize(s_must_finalize_logger);

#if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
  deallocate(s_shared_mesh_buffer);
#endif
}

}  // end namespace quest
}  // end namespace axom
