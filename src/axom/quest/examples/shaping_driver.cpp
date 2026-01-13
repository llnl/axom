// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file shaping_driver.cpp
 * \brief Driver application for shaping material volume fractions onto a simulation mesh
 */

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/klee.hpp"
#include "axom/quest.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

// NOTE: The shaping driver requires Axom to be configured with conduit or mfem and
// the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION CMake option
#if !defined(AXOM_USE_MFEM) && !defined(AXOM_USE_CONDUIT)
  #error Shaping functionality requires Axom to be configured with Conduit or MFEM and the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION option
#endif

#include "mfem.hpp"

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// C/C++ includes
#include <string>
#include <vector>
#include <memory>

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;
namespace sidre = axom::sidre;

using VolFracSampling = quest::shaping::VolFracSampling;
using SamplingMethod = quest::SamplingShaper::SamplingMethod;

//------------------------------------------------------------------------------

/// Struct to help choose our shaping method: sampling or intersection for now
enum class ShapingMethod : int
{
  Sampling,
  Intersection
};

using RuntimePolicy = axom::runtime_policy::Policy;

/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string meshFile;

  // Inline mesh parameters
  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int boxDim {-1};

  std::string shapeFile;
  klee::ShapeSet shapeSet;

  ShapingMethod shapingMethod {ShapingMethod::Sampling};
  SamplingMethod samplingMethod {SamplingMethod::InOut};
  RuntimePolicy policy {RuntimePolicy::seq};
  int quadratureOrder {5};
  int outputOrder {2};
  int samplesPerKnotSpan {25};
  int refinementLevel {7};
  double weldThresh {1e-9};
  double percentError {-1.};
  std::string annotationMode {"none"};

  std::string backgroundMaterial;

  VolFracSampling vfSampling {VolFracSampling::SAMPLE_AT_QPTS};

private:
  bool m_verboseOutput {false};

public:
  bool isVerbose() const { return m_verboseOutput; }

  /// Generate an mfem Cartesian mesh, scaled to the bounding box range
  mfem::Mesh* createBoxMesh()
  {
    mfem::Mesh* mesh = nullptr;

    switch(boxDim)
    {
    case 2:
    {
      using BBox2D = primal::BoundingBox<double, 2>;
      using Pt2D = primal::Point<double, 2>;
      auto res = axom::NumericArray<int, 2>(boxResolution.data());
      auto bbox = BBox2D(Pt2D(boxMins.data()), Pt2D(boxMaxs.data()));

      SLIC_INFO_ROOT(
        axom::fmt::format("Creating inline box mesh of resolution {} and bounding box {}", res, bbox));

      mesh = quest::util::make_cartesian_mfem_mesh_2D(bbox, res, outputOrder);
    }
    break;
    case 3:
    {
      using BBox3D = primal::BoundingBox<double, 3>;
      using Pt3D = primal::Point<double, 3>;
      auto res = axom::NumericArray<int, 3>(boxResolution.data());
      auto bbox = BBox3D(Pt3D(boxMins.data()), Pt3D(boxMaxs.data()));

      SLIC_INFO_ROOT(
        axom::fmt::format("Creating inline box mesh of resolution {} and bounding box {}", res, bbox));

      mesh = quest::util::make_cartesian_mfem_mesh_3D(bbox, res, outputOrder);
    }
    break;
    default:
      SLIC_ERROR_ROOT("Only 2D and 3D meshes are currently supported.");
      break;
    }

    // Handle conversion to parallel mfem mesh
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    {
      int* partitioning = nullptr;
      int part_method = 0;
      mfem::Mesh* pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *mesh, partitioning, part_method);
      delete[] partitioning;
      delete mesh;
      mesh = pmesh;
    }
#endif

    return mesh;
  }

  std::unique_ptr<sidre::MFEMSidreDataCollection> loadComputationalMesh()
  {
    constexpr bool dc_owns_data = true;
    mfem::Mesh* mesh = meshFile.empty() ? createBoxMesh() : nullptr;
    std::string name = meshFile.empty() ? "mesh" : getDCMeshName();

    auto dc = std::unique_ptr<sidre::MFEMSidreDataCollection>(
      new sidre::MFEMSidreDataCollection(name, mesh, dc_owns_data));
    dc->SetComm(MPI_COMM_WORLD);

    if(!meshFile.empty())
    {
      dc->Load(meshFile, "sidre_hdf5");
    }

    return dc;
  }

  std::string getDCMeshName() const
  {
    using axom::utilities::string::removeSuffix;

    // Remove the parent directories and file suffix
    std::string name = axom::Path(meshFile).baseName();
    name = removeSuffix(name, ".root");

    return name;
  }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-i,--shape-file", shapeFile)
      ->description("Path to input shape file")
      ->check(axom::CLI::ExistingFile)
      ->required();

    app.add_flag("-v,--verbose,!--no-verbose", m_verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("-n,--segments-per-knot-span", samplesPerKnotSpan)
      ->description("(2D only) Number of linear segments to generate per NURBS knot span")
      ->capture_default_str()
      ->check(axom::CLI::PositiveNumber);

    app.add_option("-t,--weld-threshold", weldThresh)
      ->description("Threshold for welding")
      ->check(axom::CLI::NonNegativeNumber)
      ->capture_default_str();

    app.add_option("-e,--percent-error", percentError)
      ->description(
        "Percent error used for calculating curve refinement and revolved volume.\n"
        "If this value is provided then dynamic curve refinement will be used\n"
        "instead of segment-based curve refinement.")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();

    std::map<std::string, ShapingMethod> methodMap {{"sampling", ShapingMethod::Sampling},
                                                    {"intersection", ShapingMethod::Intersection}};
    app.add_option("--method", shapingMethod)
      ->description("Determines the shaping method -- either sampling or intersection")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(methodMap, axom::CLI::ignore_case));

    std::map<std::string, SamplingMethod> sMethodMap {
      {"inout", SamplingMethod::InOut},
      {"windingnumber", SamplingMethod::WindingNumber}};
    app.add_option("--sampling", samplingMethod)
      ->description(
        "Determines the sampling method for the sampling shaper -- either inout or windingnumber")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(sMethodMap, axom::CLI::ignore_case));

#ifdef AXOM_USE_CALIPER
    app.add_option("--caliper", annotationMode)
      ->description(
        "caliper annotation mode. Valid options include 'none' and 'report'. "
        "Use 'help' to see full list.")
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    // use either an input mesh file or a simple inline Cartesian mesh
    {
      auto* mesh_file = app.add_option("-m,--mesh-file", meshFile)
                          ->description(
                            "Path to computational mesh. \n"
                            "Alternatively, use the `inline_mesh` subcommand.")
                          ->check(axom::CLI::ExistingFile);

      auto* inline_mesh_subcommand = app.add_subcommand("inline_mesh")
                                       ->description("Options for setting up a simple inline mesh")
                                       ->fallthrough();

      inline_mesh_subcommand->add_option("--min", boxMins)
        ->description("Min bounds for box mesh (x,y[,z])")
        ->expected(2, 3)
        ->required();
      inline_mesh_subcommand->add_option("--max", boxMaxs)
        ->description("Max bounds for box mesh (x,y[,z])")
        ->expected(2, 3)
        ->required();

      inline_mesh_subcommand->add_option("--res", boxResolution)
        ->description("Resolution of the box mesh (i,j[,k])")
        ->expected(2, 3)
        ->required();

      auto* inline_mesh_dim = inline_mesh_subcommand->add_option("-d,--dimension", boxDim)
                                ->description("Dimension of the box mesh")
                                ->check(axom::CLI::PositiveNumber)
                                ->required();

      // we want either the mesh_file or an inline mesh
      mesh_file->excludes(inline_mesh_dim);
      inline_mesh_dim->excludes(mesh_file);
    }

    app.add_option("--background-material", backgroundMaterial)
      ->description("Sets the name of the background material");

    // parameters that only apply to the sampling method
    {
      auto* sampling_options =
        app.add_option_group("sampling", "Options related to sampling-based queries");

      sampling_options->add_option("-o,--order", outputOrder)
        ->description("Order of the output grid function")
        ->capture_default_str()
        ->check(axom::CLI::NonNegativeNumber);

      sampling_options->add_option("-q,--quadrature-order", quadratureOrder)
        ->description(
          "Quadrature order for sampling the inout field. \n"
          "Determines number of samples per element in determining volume fraction field")
        ->capture_default_str()
        ->check(axom::CLI::PositiveNumber);

      std::map<std::string, VolFracSampling> vfsamplingMap {
        {"qpts", VolFracSampling::SAMPLE_AT_QPTS},
        {"dofs", VolFracSampling::SAMPLE_AT_DOFS}};
      sampling_options->add_option("-s,--sampling-type", vfSampling)
        ->description(
          "Sampling strategy. \n"
          "Sampling either at quadrature points or collocated with degrees of freedom")
        ->capture_default_str()
        ->transform(axom::CLI::CheckedTransformer(vfsamplingMap, axom::CLI::ignore_case));
    }

    // parameters that only apply to the intersection method
    {
      auto* intersection_options =
        app.add_option_group("intersection", "Options related to intersection-based queries");

      intersection_options->add_option("-r, --refinements", refinementLevel)
        ->description("(3D only) Number of refinements to perform for revolved contour")
        ->capture_default_str()
        ->check(axom::CLI::NonNegativeNumber);

      std::stringstream pol_sstr;
      pol_sstr << "Set runtime policy for intersection-based sampling method.";
      pol_sstr << "\nSet to 'seq' or 0 to use the sequential policy.";
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
      pol_sstr << "\nSet to 'omp' or 1 to use the RAJA OpenMP policy.";
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
      pol_sstr << "\nSet to 'cuda' or 2 to use the RAJA CUDA policy.";
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_HIP)
      pol_sstr << "\nSet to 'hip' or 3 to use the RAJA HIP policy.";
#endif

      intersection_options->add_option("-p, --policy", policy, pol_sstr.str())
        ->capture_default_str()
        ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));
    }
    app.get_formatter()->column_width(50);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug : slic::message::Info);
  }
};

/**
 * \brief Print some info about the mesh
 *
 * \note In MPI-based configurations, this is a collective call, but only prints on rank 0
 */
void printMeshInfo(mfem::Mesh* mesh, const std::string& prefixMessage = "")
{
  namespace primal = axom::primal;

  int myRank = 0;
#ifdef AXOM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

  int numElements = mesh->GetNE();

  mfem::Vector mins, maxs;
#ifdef MFEM_USE_MPI
  auto* pmesh = dynamic_cast<mfem::ParMesh*>(mesh);
  if(pmesh != nullptr)
  {
    pmesh->GetBoundingBox(mins, maxs);
    numElements = pmesh->ReduceInt(numElements);
    myRank = pmesh->GetMyRank();
  }
  else
#endif
  {
    mesh->GetBoundingBox(mins, maxs);
  }

  if(myRank == 0)
  {
    switch(mesh->Dimension())
    {
    case 2:
      SLIC_INFO(axom::fmt::format(
        axom::utilities::locale(),
        "{} mesh has {:L} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 2>(primal::Point<double, 2>(mins.GetData()),
                                       primal::Point<double, 2>(maxs.GetData()))));
      break;
    case 3:
      SLIC_INFO(axom::fmt::format(
        axom::utilities::locale(),
        "{} mesh has {:L} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 3>(primal::Point<double, 3>(mins.GetData()),
                                       primal::Point<double, 3>(maxs.GetData()))));
      break;
    }
  }

  slic::flushStreams();
}

/// \brief Utility function to initialize the logger
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  slic::LogStream* logStream {nullptr};

#ifdef AXOM_USE_MPI
  int num_ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  slic::setIsRoot(my_rank == 0);

  if(num_ranks > 1)
  {
    std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
    const int RLIMIT = 8;
    logStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
    logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
  }
  else
#endif  // AXOM_USE_MPI
  {
    std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
    logStream = new slic::GenericOutputStream(&std::cout, fmt);
  }

  slic::addStreamToAllMsgLevels(logStream);
}

/// \brief Utility function to finalize the logger
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);
  const int my_rank = mpi_raii_wrapper.my_rank();

  initializeLogger();

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  Input params;
  axom::CLI::App app {"Driver for Klee shaping query"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(my_rank == 0)
    {
      retval = app.exit(e);
    }
    finalizeLogger();

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    exit(retval);
  }

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(params.annotationMode);

  AXOM_ANNOTATE_BEGIN("quest shaping example");
  AXOM_ANNOTATE_BEGIN("init");

  //---------------------------------------------------------------------------
  // Load the klee shape file and extract some information
  //---------------------------------------------------------------------------
  try
  {
    AXOM_ANNOTATE_SCOPE("read Klee shape set");
    params.shapeSet = klee::readShapeSet(params.shapeFile);

    slic::flushStreams();
  }
  catch(klee::KleeError& error)
  {
    std::vector<std::string> errs;
    for(auto verificationError : error.getErrors())
    {
      errs.push_back(axom::fmt::format(" - '{}': {}",
                                       static_cast<std::string>(verificationError.path),
                                       verificationError.message));
    }

    SLIC_WARNING(
      axom::fmt::format("Error during parsing klee input. Found the following errors:\n{}",
                        axom::fmt::join(errs, "\n")));

    finalizeLogger();
    exit(1);
  }

  AXOM_ANNOTATE_BEGIN("load mesh");
  //---------------------------------------------------------------------------
  // Load the computational mesh
  //---------------------------------------------------------------------------
  auto originalMeshDC = params.loadComputationalMesh();

  //---------------------------------------------------------------------------
  // Set up DataCollection for shaping
  //---------------------------------------------------------------------------
  mfem::Mesh* shapingMesh = nullptr;
  constexpr bool dc_owns_data = true;
  sidre::MFEMSidreDataCollection shapingDC("shaping", shapingMesh, dc_owns_data);
  {
    shapingDC.SetMeshNodesName("positions");

    auto* pmesh = dynamic_cast<mfem::ParMesh*>(originalMeshDC->GetMesh());
    shapingMesh =
      (pmesh != nullptr) ? new mfem::ParMesh(*pmesh) : new mfem::Mesh(*originalMeshDC->GetMesh());
    shapingDC.SetMesh(shapingMesh);
  }
  AXOM_ANNOTATE_END("load mesh");
  printMeshInfo(shapingDC.GetMesh(), "After loading");

  //---------------------------------------------------------------------------
  // Initialize the shaping query object
  //---------------------------------------------------------------------------
  AXOM_ANNOTATE_BEGIN("setup shaping problem");
  quest::Shaper* shaper = nullptr;
  switch(params.shapingMethod)
  {
  case ShapingMethod::Sampling:
    shaper = new quest::SamplingShaper(params.policy,
                                       axom::policyToDefaultAllocatorID(params.policy),
                                       params.shapeSet,
                                       &shapingDC);
    break;
  case ShapingMethod::Intersection:
    shaper = new quest::IntersectionShaper(params.policy,
                                           axom::policyToDefaultAllocatorID(params.policy),
                                           params.shapeSet,
                                           &shapingDC);
    break;
  }
  SLIC_ASSERT_MSG(shaper != nullptr, "Invalid shaping method selected!");

  // Set generic parameters for the base Shaper instance
  shaper->setSamplesPerKnotSpan(params.samplesPerKnotSpan);
  shaper->setVertexWeldThreshold(params.weldThresh);
  shaper->setVerbosity(params.isVerbose());
  if(params.percentError > 0.)
  {
    shaper->setPercentError(params.percentError);
    shaper->setRefinementType(quest::DiscreteShape::RefinementDynamic);
  }

  // Associate any fields that begin with "vol_frac" with "material" so when
  // the data collection is written, a matset will be created.
  shaper->getDC()->AssociateMaterialSet("vol_frac", "material");

  // Set specific parameters for a SamplingShaper, if appropriate
  if(auto* samplingShaper = dynamic_cast<quest::SamplingShaper*>(shaper))
  {
    samplingShaper->setSamplingType(params.vfSampling);
    samplingShaper->setQuadratureOrder(params.quadratureOrder);
    samplingShaper->setVolumeFractionOrder(params.outputOrder);
    samplingShaper->setSamplingMethod(params.samplingMethod);

    // register point projectors
    if(shapingDC.GetMesh()->Dimension() == 3)
    {
      samplingShaper->setPointProjector32([](primal::Point<double, 3> pt) {
        const double& x = pt[0];
        const double& y = pt[1];
        const double& z = pt[2];
        return primal::Point<double, 2> {z, sqrt(x * x + y * y)};
      });
    }
    else if(shapingDC.GetMesh()->Dimension() == 2)
    {
      samplingShaper->setPointProjector23([](primal::Point<double, 2> pt) {
        const double& x = pt[0];
        const double& y = pt[1];
        return primal::Point<double, 3> {x, y, 0.};
      });
    }
  }

  // Set specific parameters here for IntersectionShaper
  if(auto* intersectionShaper = dynamic_cast<quest::IntersectionShaper*>(shaper))
  {
    intersectionShaper->setLevel(params.refinementLevel);

    if(!params.backgroundMaterial.empty())
    {
      intersectionShaper->setFreeMaterialName(params.backgroundMaterial);
    }
  }

  //---------------------------------------------------------------------------
  // Project initial volume fractions, if applicable
  //---------------------------------------------------------------------------
  if(auto* samplingShaper = dynamic_cast<quest::SamplingShaper*>(shaper))
  {
    AXOM_ANNOTATE_SCOPE("import initial volume fractions");
    std::map<std::string, mfem::GridFunction*> initial_grid_functions;

    // Generate a background material (w/ volume fractions set to 1) if user provided a name
    if(!params.backgroundMaterial.empty())
    {
      auto material = params.backgroundMaterial;
      auto name = axom::fmt::format("vol_frac_{}", material);

      const int order = params.outputOrder;
      const int dim = shapingMesh->Dimension();
      const auto basis = mfem::BasisType::Positive;

      auto* coll = new mfem::L2_FECollection(order, dim, basis);
      auto* fes = new mfem::FiniteElementSpace(shapingDC.GetMesh(), coll);
      const int sz = fes->GetVSize();

      auto* view = shapingDC.AllocNamedBuffer(name, sz);
      auto* volFrac = new mfem::GridFunction(fes, view->getArray());
      volFrac->MakeOwner(coll);

      (*volFrac) = 1.;

      shapingDC.RegisterField(name, volFrac);

      initial_grid_functions[material] = shapingDC.GetField(name);
    }

    // Project provided volume fraction grid functions as quadrature point data
    samplingShaper->importInitialVolumeFractions(initial_grid_functions);
  }
  AXOM_ANNOTATE_END("setup shaping problem");
  AXOM_ANNOTATE_END("init");

  //---------------------------------------------------------------------------
  // Process each of the shapes
  //---------------------------------------------------------------------------
  SLIC_INFO(axom::fmt::format("{:=^80}", "Sampling InOut fields for shapes"));
  AXOM_ANNOTATE_BEGIN("shaping");
  for(const auto& shape : params.shapeSet.getShapes())
  {
    const std::string shapeFormat = shape.getGeometry().getFormat();
    SLIC_INFO(
      axom::fmt::format("{:-^80}",
                        axom::fmt::format("Processing shape '{}' of material '{}' (format '{}')",
                                          shape.getName(),
                                          shape.getMaterial(),
                                          shapeFormat)));

    const klee::Dimensions shapeDim = shape.getGeometry().getInputDimensions();

    // Apply error checking
#ifndef AXOM_USE_C2C
    SLIC_ERROR_IF(shapeDim == klee::Dimensions::Two && shapeFormat == "c2c",
                  "Shaping with contour files requires an Axom configured with "
                  "the C2C library");
#endif

    // Load the shape from file. This also applies any transformations.
    shaper->loadShape(shape);
    slic::flushStreams();

    // Generate a spatial index over the shape
    shaper->prepareShapeQuery(shapeDim, shape);
    slic::flushStreams();

    // Query the mesh against this shape
    shaper->runShapeQuery(shape);
    slic::flushStreams();

    // Apply the replacement rules for this shape against the existing materials
    shaper->applyReplacementRules(shape);
    slic::flushStreams();

    // Finalize data structures associated with this shape and spatial index
    shaper->finalizeShapeQuery();
    slic::flushStreams();
  }
  AXOM_ANNOTATE_END("shaping");

  //---------------------------------------------------------------------------
  // After shaping in all shapes, generate/adjust the material volume fractions
  //---------------------------------------------------------------------------
  AXOM_ANNOTATE_BEGIN("adjust");
  SLIC_INFO(axom::fmt::format("{:=^80}", "Generating volume fraction fields for materials"));

  shaper->adjustVolumeFractions();

  //---------------------------------------------------------------------------
  // Compute and print volumes of each material's volume fraction
  //---------------------------------------------------------------------------
  using axom::utilities::string::startsWith;
  for(auto& kv : shaper->getDC()->GetFieldMap())
  {
    if(startsWith(kv.first, "vol_frac_"))
    {
      const auto mat_name = kv.first.substr(9);
      auto* gf = kv.second;

      mfem::ConstantCoefficient one(1.0);
      mfem::LinearForm vol_form(gf->FESpace());
      vol_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
      vol_form.Assemble();

      const double volume = shaper->allReduceSum(*gf * vol_form);

      SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                  "Volume of material '{}' is {:.6Lf}",
                                  mat_name,
                                  volume));
    }
  }
  AXOM_ANNOTATE_END("adjust");

  //---------------------------------------------------------------------------
  // Save meshes and fields
  //---------------------------------------------------------------------------
  if(params.isVerbose())
  {
    if(auto* samplingShaper = dynamic_cast<quest::SamplingShaper*>(shaper))
    {
      SLIC_INFO(axom::fmt::format("{:-^80}", ""));
      samplingShaper->printRegisteredFieldNames(" -- after shaping");
    }
  }

#ifdef MFEM_USE_MPI
  {
    AXOM_ANNOTATE_SCOPE("save shaping results");
    shaper->getDC()->Save();
  }
#endif

  delete shaper;

  //---------------------------------------------------------------------------
  // Cleanup and exit
  //---------------------------------------------------------------------------
  SLIC_INFO(axom::fmt::format("{:-^80}", ""));
  slic::flushStreams();

  AXOM_ANNOTATE_END("quest shaping example");

  finalizeLogger();

  return 0;
}
