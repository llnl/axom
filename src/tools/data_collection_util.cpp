// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

// MFEM is required
#ifndef AXOM_USE_MFEM
  #error This example requires an Axom configured with MFEM
#endif
#include "mfem.hpp"

// MPI is optional
#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

namespace slic = axom::slic;
namespace sidre = axom::sidre;
namespace primal = axom::primal;
namespace quest = axom::quest;

//------------------------------------------------------------------------------

/// Struct to help choose if we're loading a box (Cartesian) mesh or a file
enum class MeshForm : int
{
  Box,
  File
};

/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string dcName;

  int uniformRefinements {0};
  int polynomialOrder {2};

  MeshForm meshForm {MeshForm::Box};

  // Box mesh parameters
  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int boxDim {-1};

  // File mesh parameters
  std::string mfemFile;
  std::vector<double> fileScale;
  std::vector<double> fileTranslate;

private:
  bool m_verboseOutput {false};

  /**
   * Fixes up parameters associated with a box mesh based on user-provided input
   * \throws axom::CLI::ValidationError if the options are invalid or insufficient
   */
  void fixBoxParams()
  {
    const bool dimProvided = boxDim > 0;
    const bool rangeProvided = !(boxMins.empty() && boxMaxs.empty());
    const bool resProvided = !boxResolution.empty();

    // First check: we need at least one of: dimension, mins/maxs, resolution
    if(!(dimProvided || rangeProvided || resProvided))
    {
      throw axom::CLI::ValidationError(
        "box",
        "Box mesh must have at least one of: dimension, mins and maxs or "
        "resolution to determine the dimension");
    }

    int szMins = boxMins.size();
    int szMaxs = boxMaxs.size();
    int szRes = boxResolution.size();

    // Error checking on provided dimension
    if(dimProvided)
    {
      if(boxDim <= 1 || boxDim > 3)
      {
        throw axom::CLI::ValidationError(
          "box",
          axom::fmt::format("Invalid dimension: {}. Only 2D and 3D supported", boxDim));
      }

      // Ensure that range and resolution have the right number of entries
      if(rangeProvided && (szMins != boxDim || szMaxs != boxDim))
      {
        throw axom::CLI::ValidationError(
          "box",
          "Bounding box has different dimension than provided dimension");
      }
      if(resProvided && (szRes != boxDim))
      {
        throw axom::CLI::ValidationError(
          "box",
          "Box resolution has different dimension than provided dimension");
      }
    }

    // Error checking on provided range; note: accounts for previous checks on dim
    if(rangeProvided && (szMins != szMaxs))
    {
      throw axom::CLI::ValidationError("box", "Bounding box mins and maxs has different dimensions");
    }

    // Error checking on provided range and resolution; note: accounts for previous checks on dim and range
    if(rangeProvided && resProvided && (szMins != szRes))
    {
      throw axom::CLI::ValidationError(
        "box",
        "Bounding box mins and maxs has different dimensions than resolution");
    }

    // Now that we've checked errors, let's fill in missing data
    const int dim = dimProvided ? boxDim : (resProvided ? szRes : szMins);
    if(!dimProvided)
    {
      boxDim = dim;
    }

    if(!rangeProvided)  // set range to unit box
    {
      boxMins.resize(dim);
      boxMaxs.resize(dim);
      for(int i = 0; i < dim; ++i)
      {
        boxMins[i] = 0.;
        boxMaxs[i] = 1.;
      }
    }

    if(!resProvided)  // set resolution to 1
    {
      boxResolution.resize(dim);
      for(int i = 0; i < dim; ++i)
      {
        boxResolution[i] = 1;
      }
    }
  }

public:
  Input() = default;

  bool isVerbose() const { return m_verboseOutput; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    // Options that are always available
    app.add_option("-o, --output-name", dcName)
      ->description(
        "Name of the output mesh. Defaults to box_{2,3}d for box meshes, "
        "and the mfem mesh name when loading from an mfem file");

    app.add_flag("-v,--verbose", m_verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("-l,--ref-level", uniformRefinements)
      ->description("The number of uniform refinement levels to apply to the mesh")
      ->capture_default_str()
      ->check(axom::CLI::NonNegativeNumber);

    // Parameter to determine if we're using a file or a box mesh
    std::map<std::string, MeshForm> meshFormMap {{"box", MeshForm::Box}, {"file", MeshForm::File}};
    app.add_option("-f,--mesh-form", meshForm)
      ->description("Determines the input type -- either box or file")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(meshFormMap, axom::CLI::ignore_case));

    // Parameters for the box mesh option
    auto* box_options = app.add_option_group("box", "Options for setting up a box mesh");
    auto* minbb = box_options->add_option("--min", boxMins)
                    ->description("Min bounds for box mesh (x,y[,z])")
                    ->expected(2, 3);
    auto* maxbb = box_options->add_option("--max", boxMaxs)
                    ->description("Max bounds for box mesh (x,y[,z])")
                    ->expected(2, 3);
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    box_options->add_option("--res", boxResolution)
      ->description("Resolution of the box mesh (i,j[,k])")
      ->expected(2, 3);

    box_options->add_option("-d,--dimension", boxDim)
      ->description("Dimension of the box mesh")
      ->check(axom::CLI::PositiveNumber);

    box_options->add_option("-p,--polynomial-order", polynomialOrder)
      ->description("polynomial order of the generated mesh")
      ->capture_default_str()
      ->check(axom::CLI::NonNegativeNumber);

    // Parameters for the 'file' option
    auto* file_options = app.add_option_group("file", "Options for loading from an mfem mesh file");

    file_options->add_option("-m, --mfem-file", mfemFile)
      ->description("Path to a mesh file in the mfem format")
      ->check(axom::CLI::ExistingFile);

    file_options->add_option("--scale", fileScale)
      ->description(
        "Scale factors for the file mesh. Can either be one value or "
        "mesh.Dimension() values")
      ->expected(1, 3);
    file_options->add_option("--translate", fileTranslate)
      ->description("Translation for the file mesh (sx,sy[,sz])")
      ->expected(2, 3);

    app.get_formatter()->column_width(50);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug : slic::message::Info);

    // Fix up some parameters
    if(meshForm == MeshForm::Box)
    {
      fixBoxParams();  // Note: throws on error conditions
      if(dcName.empty())
      {
        dcName = axom::fmt::format("box_{}d", boxDim);
      }
    }
    else
    {
      if(mfemFile.empty())
      {
        throw axom::CLI::ValidationError("file", "'mfemFile' required when `--mesh-form == File`");
      }

      if(dcName.empty())
      {
        using axom::utilities::string::endsWith;
        dcName = axom::Path(mfemFile).baseName();

        const std::string suffix = ".mesh";
        if(endsWith(dcName, suffix))
        {
          dcName = dcName.substr(0, dcName.size() - suffix.size());
        }
      }
    }
  }
};

mfem::Mesh* createBoxMesh(const Input& params)
{
  // Create a background mesh
  // Generate an mfem Cartesian mesh, scaled to the bounding box range

  const int dim = params.boxDim;
  auto& lo = params.boxMins;
  auto& hi = params.boxMaxs;

  mfem::Mesh* mesh = nullptr;

  switch(dim)
  {
  case 2:
  {
    using Pt2D = primal::Point<double, 2>;
    auto res = axom::NumericArray<int, 2>(params.boxResolution.data());
    auto bbox = primal::BoundingBox<double, 2>(Pt2D(lo.data()), Pt2D(hi.data()));

    SLIC_INFO(
      axom::fmt::format("Creating a box mesh of resolution {} and bounding box {}", res, bbox));

    mesh = quest::util::make_cartesian_mfem_mesh_2D(bbox, res, params.polynomialOrder);
  }
  break;
  case 3:
  {
    using Pt3D = primal::Point<double, 3>;
    auto res = axom::NumericArray<int, 3>(params.boxResolution.data());
    auto bbox = primal::BoundingBox<double, 3>(Pt3D(lo.data()), Pt3D(hi.data()));

    SLIC_INFO(
      axom::fmt::format("Creating a box mesh of resolution {} and bounding box {}", res, bbox));

    mesh = quest::util::make_cartesian_mfem_mesh_3D(bbox, res, params.polynomialOrder);
  }
  break;
  default:
    SLIC_ERROR("Only 2D and 3D meshes are currently supported.");
    break;
  }

  return mesh;
}

mfem::Mesh* loadFileMesh(const Input& params)
{
  mfem::Mesh* mesh = new mfem::Mesh(params.mfemFile.c_str(), 1, 1);
  mesh->EnsureNodes();

  const bool shouldScale = !params.fileScale.empty();
  const bool shouldTranslate = !params.fileTranslate.empty();
  bool shouldXForm = shouldScale || shouldTranslate;
  if(shouldXForm)
  {
    const int dim = mesh->Dimension();

    // Set up scaling vector sc
    mfem::Vector sc(dim);
    {
      const int numProvidedScales = params.fileScale.size();
      switch(numProvidedScales)
      {
      case 0:
        for(int d = 0; d < dim; ++d)
        {
          sc(d) = 1.;
        }
        break;
      case 1:
        for(int d = 0; d < dim; ++d)
        {
          sc(d) = params.fileScale[0];
        }
        break;
      default:
        SLIC_ERROR_IF(dim != numProvidedScales,
                      axom::fmt::format("Incorrect number of scale values. Expected {} got {}",
                                        dim,
                                        numProvidedScales));
        for(int d = 0; d < dim; ++d)
        {
          sc(d) = params.fileScale[d];
        }
        break;
      }
    }

    // Set up translation vector tr
    mfem::Vector tr(dim);
    {
      const int numProvidedTranslations = params.fileTranslate.size();
      switch(numProvidedTranslations)
      {
      case 0:
        for(int d = 0; d < dim; ++d)
        {
          tr(d) = 0.;
        }
        break;
      default:
        SLIC_ERROR_IF(
          dim != numProvidedTranslations,
          axom::fmt::format("Incorrect number of translations values. Expected {} got {}",
                            dim,
                            numProvidedTranslations));
        for(int d = 0; d < dim; ++d)
        {
          tr(d) = params.fileTranslate[d];
        }
        break;
      }
    }

    // Get approximate center of mesh for scaling
    mfem::Vector ctr(dim);
    {
      mfem::Vector mins, maxs;
      mesh->GetBoundingBox(mins, maxs);
      for(int d = 0; d < dim; ++d)
      {
        ctr(d) = (mins(d) + maxs(d)) / 2.;
      }
    }

    // Perform transformation on nodal dofs: scaling about center followed by translation
    auto* nodes = mesh->GetNodes();
    auto* fes = nodes->FESpace();
    const int NDOFS = fes->GetNDofs();
    for(int i = 0; i < NDOFS; ++i)
    {
      for(int d = 0; d < dim; ++d)
      {
        double val = (*nodes)(fes->DofToVDof(i, d));
        val = sc(d) * (val - ctr(d)) + ctr(d) + tr(d);
        (*nodes)(fes->DofToVDof(i, d)) = val;
      }
    }
  }

  return mesh;
}

/**
 * \brief Print some info about the mesh
 *
 * \note In MPI configurations, this is a collective call, but only prints on rank 0
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
        "{} mesh has {} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 2>(primal::Point<double, 2>(mins.GetData()),
                                       primal::Point<double, 2>(maxs.GetData()))));
      break;
    case 3:
      SLIC_INFO(axom::fmt::format(
        "{} mesh has {} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 3>(primal::Point<double, 3>(mins.GetData()),
                                       primal::Point<double, 3>(maxs.GetData()))));
      break;
    }
  }

  slic::flushStreams();
}

/*!
 * \brief Utility function to initialize the logger
 */
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  slic::LogStream* logStream {nullptr};

#ifdef AXOM_USE_MPI
  int num_ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
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

/*!
 * \brief Utility function to finalize the logger
 */
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
#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#else
  int my_rank = 0;
  int num_ranks = 1;
#endif

  initializeLogger();

  // Set up and parse command line arguments
  Input params;
  axom::CLI::App app {"Utility tool to create a blueprint compliant data store"};

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
    MPI_Finalize();
#endif
    exit(retval);
  }

  // Create the mfem sidre data collection
  sidre::MFEMSidreDataCollection dc(params.dcName, nullptr, true);

  // Create or load the serial mfem mesh
  mfem::Mesh* mesh {nullptr};
  switch(params.meshForm)
  {
  case MeshForm::Box:
    mesh = createBoxMesh(params);
    break;

  case MeshForm::File:
    mesh = loadFileMesh(params);
    break;
  }

  // Apply uniform refinements
  for(int i = 0; i < params.uniformRefinements; ++i)
  {
    mesh->UniformRefinement();
  }

  // Print out some info about the mesh, including the approximate mesh bounding box
  // before parallel decomposition
  printMeshInfo(mesh, "After transformations,");

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

  // Add mfem mesh to data collection
  dc.SetMeshNodesName("positions");
#ifdef MFEM_USE_MPI
  dc.SetMesh(MPI_COMM_WORLD, mesh);
#else
  dc.SetMesh(mesh);
#endif

  // TODO: Optionally convert to low order mesh ?
  SLIC_INFO(axom::fmt::format("Saving mesh file '{}' in directory '{}'",
                              dc.GetCollectionName(),
                              axom::utilities::filesystem::getCWD()));
  dc.Save();

  // Cleanup and exit
  finalizeLogger();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
