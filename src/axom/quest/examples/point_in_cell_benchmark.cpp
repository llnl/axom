// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/quest/PointInCell.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#ifdef AXOM_USE_MFEM
  #include "axom/quest/detail/PointInCellMeshWrapper_mfem.hpp"
#else
  #error "Quest's PointInCell on mfem meshes requires mfem library."
#endif

// namespace aliases
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;
namespace utilities = axom::utilities;

namespace
{
#ifdef AXOM_USE_GPU
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
using GPUExec = axom::HIP_EXEC<256>;
  #elif defined(AXOM_RUNTIME_POLICY_USE_CUDA)
using GPUExec = axom::CUDA_EXEC<256>;
  #endif
#endif

enum class ExecPolicy
{
  CPU,
  OpenMP,
  GPU
};

// clang-format off
const std::map<std::string, ExecPolicy> validExecPolicies
{
    {"seq", ExecPolicy::CPU}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  , {"omp", ExecPolicy::OpenMP}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP) || defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  , {"gpu", ExecPolicy::GPU}
#endif
};

const std::map<std::string, mfem::InverseElementTransformation::InitGuessType> validGuessTypes
{
  {"center", mfem::InverseElementTransformation::Center},
  {"phys"  , mfem::InverseElementTransformation::ClosestPhysNode},
  {"ref"   , mfem::InverseElementTransformation::ClosestRefNode}
};

const std::map<std::string, mfem::InverseElementTransformation::SolverType> validSolverProjectionTypes
{
  {"none",     mfem::InverseElementTransformation::Newton},
  {"segment",  mfem::InverseElementTransformation::NewtonSegmentProject},
  {"boundary", mfem::InverseElementTransformation::NewtonElementProject}
};
// clang-format on

struct Arguments
{
  std::string file_name;
  int num_rand_pts {10000};
  int num_bins {25};
  ExecPolicy exec_space {ExecPolicy::CPU};
  bool should_verify_points {false};

  int verbosity {-1};
  int init_guess_type {mfem::InverseElementTransformation::ClosestPhysNode};
  int init_guess_order {-1};
  int solver_projection_type {mfem::InverseElementTransformation::NewtonElementProject};

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    std::string pol_info = "Sets execution space of the PointInCell query.\n";
    pol_info += "Set to 'seq' to use sequential execution policy.";
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    pol_info += "\nSet to 'omp' to use an OpenMP execution policy.";
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP) || defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    pol_info += "\nSet to 'gpu' to use a GPU execution policy.";
#endif

    app.add_option("-f,--file", this->file_name, "specifies the input mesh file")
      ->check(axom::CLI::ExistingFile)
      ->required();

    app.add_option("-n,--num_pts", this->num_rand_pts)
      ->description("the number of points to query")
      ->capture_default_str();

    app.add_option("-b,--num-bins", this->num_bins)
      ->description("the number of bins to construct for each dimension")
      ->capture_default_str();

    app.add_flag("--verify", this->should_verify_points)
      ->description("verify query by reconstructing points after locating them")
      ->capture_default_str();

    app.add_option("-e, --exec_space", this->exec_space, pol_info)
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(validExecPolicies));

    auto* xform_opts = app.add_option_group("inverse transformation options")
                         ->description(
                           "Extra parameters for controlling the per-element "
                           "point-in-cell queries");

    xform_opts->add_option("-v, --verbosity", this->verbosity)
      ->description(
        "controls verbosity of output for point location. "
        "\n-1: minimal"
        "\n 0: print only errors"
        "\n 1: print first and last iteration"
        "\n 2: print every iteration"
        "\n 3: print every iteration including point coordinates")
      ->capture_default_str()
      ->check(axom::CLI::Range(-1, 3));

    xform_opts->add_option("--guess-type", this->init_guess_type)
      ->description(
        "controls the initial guess type"
        "\n 'center': Use center of element in reference space"
        "\n 'phys': Use closest point to a grid in physical space"
        "\n 'ref': Use closest point to a grid in reference space")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(validGuessTypes));

    xform_opts->add_option("--guess-order", this->init_guess_order)
      ->description("controls the number of grid points for the initial guess")
      ->capture_default_str();

    xform_opts->add_option("--projection-type", this->solver_projection_type)
      ->description(
        "controls handling of points that leave reference space element"
        "\n 'none': Allow search to continue outside element"
        "\n 'segment': Project along current segment leaving the element"
        "\n 'boundary': Project to closest reference space point on boundary")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(validSolverProjectionTypes));

    app.get_formatter()->column_width(40);

    // could throw an exception
    app.parse(argc, argv);
  }
};

}  // namespace

template <int NDIMS>
primal::Point<double, NDIMS> get_rand_pt(const primal::BoundingBox<double, NDIMS>& bbox)
{
  primal::Point<double, NDIMS> rnd_pt;
  for(int i = 0; i < NDIMS; i++)
  {
    rnd_pt[i] = axom::utilities::random_real(bbox.getMin()[i], bbox.getMax()[i]);
  }
  return rnd_pt;
}

template <typename ExecSpace, int DIM>
void benchmark_point_in_cell(mfem::Mesh& mesh, const Arguments& args)
{
  using PointType = primal::Point<double, DIM>;
  using BoxType = primal::BoundingBox<double, DIM>;
  using mesh_tag = axom::quest::quest_point_in_cell_mfem_tag;
  using IndexType = typename quest::PointInCellTraits<mesh_tag>::IndexType;
  constexpr auto NO_CELL = quest::PointInCellTraits<mesh_tag>::NO_CELL;

  const int npts = args.num_rand_pts;
  const int nbins = args.num_bins;
  const bool verifyPoints = args.should_verify_points;

  // Get ids of necessary allocators
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecSpace>::allocatorID();

  BoxType meshBb;
  {
    mfem::Vector meshMin, meshMax;
    mesh.GetBoundingBox(meshMin, meshMax);
    meshBb.addPoint(PointType(meshMin.GetData()));
    meshBb.addPoint(PointType(meshMax.GetData()));
  }
  SLIC_DEBUG("Mesh bounding box " << meshBb);

  axom::Array<PointType> pts_h(npts, npts, host_allocator);

  // Generate random points
  utilities::Timer timer(true);
  for(int i = 0; i < npts; i++)
  {
    pts_h[i] = get_rand_pt(meshBb);
  }
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Constructed {:L} random points in {} s.",
                              npts,
                              timer.elapsed()));

  primal::Point<int, DIM> bins(nbins);

  // Initialize the spatial index
  timer.start();
  quest::PointInCell<mesh_tag, ExecSpace> query(&mesh, bins.data());
  query.setPrintLevel(args.verbosity);
  query.setInitialGuessType(args.init_guess_type);
  query.setInitialGridOrder(args.init_guess_order);
  query.setSolverProjectionType(args.solver_projection_type);

  SLIC_INFO(axom::fmt::format("Initialized point-in-cell query in {} s.", timer.elapsed()));

  // Run query
  axom::Array<IndexType> outCellIds_d(npts, npts, device_allocator);
  axom::Array<PointType> outIsoParams_d(npts, npts, device_allocator);

  auto outCellIds_v = outCellIds_d.view();
  auto outIsoParams_v = outIsoParams_d.view();

  axom::Array<PointType> pts_d = axom::Array<PointType>(pts_h, device_allocator);

  timer.start();
  query.locatePoints(pts_d.view(), outCellIds_v.data(), outIsoParams_v.data());
  double time = timer.elapsed();
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Ran query on {:L} points in {} s -- rate: {:L} q/s",
                              npts,
                              time,
                              npts / time));

  // Copy back to host
  axom::Array<IndexType> outCellIds_h =
    on_device ? axom::Array<IndexType>(outCellIds_d, host_allocator) : std::move(outCellIds_d);
  axom::Array<PointType> outIsoParams_h =
    on_device ? axom::Array<PointType>(outIsoParams_d, host_allocator) : std::move(outIsoParams_d);

  // Verify the results by reconstructing physical points from refrerence coordinates
  if(verifyPoints)
  {
    constexpr double EPS = 1e-16;
    int num_not_found = 0;
    int num_wrong = 0;
    int num_found = 0;

    for(int i = 0; i < npts; ++i)
    {
      if(outCellIds_h[i] != NO_CELL)
      {
        PointType reconstructed;
        query.reconstructPoint(outCellIds_h[i], outIsoParams_h[i].data(), reconstructed.data());
        if(primal::squared_distance(pts_h[i], reconstructed) > EPS)
        {
          ++num_wrong;
          SLIC_DEBUG(
            axom::fmt::format("Incorrect reconstruction: Original point {}; "
                              "reconstructed point {} found in cell {} w/ isoparametric "
                              "coordinates {}; distance between these is {}",
                              pts_h[i],
                              reconstructed[i],
                              outCellIds_h[i],
                              outIsoParams_h[i],
                              sqrt(primal::squared_distance(pts_h[i], reconstructed))));
        }
        else
        {
          ++num_found;
        }
      }
      else
      {
        ++num_not_found;
        SLIC_DEBUG(axom::fmt::format("Did not reconstruct point {}", pts_h[i]));
      }
    }

    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Correctly reconstructed {:L} of {:L} points ({:.3f}%).\n"
                                "\t{:L} points were not reconstructed; "
                                "{:L} points were incorrectly reconstructed",
                                num_found,
                                npts,
                                num_found * 100. / npts,
                                num_not_found,
                                num_wrong));
  }
}

template <typename ExecSpace>
void benchmark_point_in_cell(mfem::Mesh& mesh, const Arguments& args)
{
  switch(mesh.SpaceDimension())
  {
  case 2:
    benchmark_point_in_cell<ExecSpace, 2>(mesh, args);
    break;
  case 3:
    benchmark_point_in_cell<ExecSpace, 3>(mesh, args);
    break;
  default:
    SLIC_ERROR("Unsupported dimension (" << mesh.SpaceDimension() << ")");
    break;
  }
}

int main(int argc, char** argv)
{
  slic::SimpleLogger logger(slic::message::Info);

  Arguments args;
  axom::CLI::App app {"Driver for signed distance query"};

  try
  {
    args.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    retval = app.exit(e);
    return retval;
  }

  if(args.verbosity >= 0)
  {
    slic::setLoggingMsgLevel(slic::message::Debug);
  }

  // Open MFEM mesh
  utilities::Timer timer(true);
  mfem::Mesh testMesh(args.file_name.c_str());
  testMesh.EnsureNodes();
  SLIC_INFO(
    axom::fmt::format("Initialized MFEM mesh '{}' in {} s. \n"
                      "\tThe mesh has {} {}d elements and its order is {}.",
                      args.file_name,
                      timer.elapsed(),
                      testMesh.GetNE(),
                      testMesh.Dimension(),
                      testMesh.GetNodalFESpace()->FEColl()->GetOrder()));

  // Run the benchmark
  switch(args.exec_space)
  {
  case ExecPolicy::CPU:
    benchmark_point_in_cell<axom::SEQ_EXEC>(testMesh, args);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case ExecPolicy::OpenMP:
    benchmark_point_in_cell<axom::OMP_EXEC>(testMesh, args);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP) || defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case ExecPolicy::GPU:
    benchmark_point_in_cell<GPUExec>(testMesh, args);
    break;
#endif
  default:
    SLIC_ERROR("Unsupported execution space.");
    return 1;
  }

  return 0;
}
