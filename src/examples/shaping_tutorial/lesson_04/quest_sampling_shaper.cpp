// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
/*!
 * \file quest_sampling_shaper.cpp
 * \brief An example that shapes a Klee input onto a compuational mesh
 *
 * This example demonstrates how to use the Quest SamplingShaper to shape
 * a Klee input onto a computational mesh. The program:
 * 1. Creates a structured mesh from Inlet mesh metadata
 * 2. Loads Klee shapes from a YAML file
 * 3. Generates volume fraction fields for each material based on the input shapes
 * 4. Outputs the generated mesh to disk
 * 
 * This example supports both serial and MPI execution.
 * 
 * Example run:
 * > [srun -n8] ./quest_sampling_shaper -m mesh_metadata.lua -k shapes.yaml [-v]
 * 
 */
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/klee.hpp"
#include "axom/quest.hpp"
#include "axom/inlet.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include "mfem.hpp"

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// NOTE: The shaping driver requires Axom to be configured with conduit and mfem
#if !defined(AXOM_USE_MFEM) && !defined(AXOM_USE_CONDUIT)
  #error Shaping functionality requires Axom to be configured with Conduit and MFEM
#endif

#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace slic = axom::slic;
namespace klee = axom::klee;
namespace quest = axom::quest;
namespace primal = axom::primal;
namespace inlet = axom::inlet;
namespace sidre = axom::sidre;

// Mesh metadata
struct MeshMetadata
{
public:
  int dim;
  axom::Array<double> bb_min;
  axom::Array<double> bb_max;
  axom::Array<int> resolution;

  std::string background_material;
  int volume_fraction_order {2};
  int mesh_order {1};
  int quadrature_order {5};
  quest::SamplingShaper::SamplingMethod sampling_method {quest::SamplingShaper::SamplingMethod::InOut};

public:
  /// Returns the axis-aligned bounding box as a templated primal BoundingBox primitive
  template <int DIM>
  axom::primal::BoundingBox<double, DIM> getBoundingBox() const
  {
    static_assert(DIM == 2 || DIM == 3, "Invalid dimension");
    SLIC_ASSERT_MSG(DIM == dim, "Template dimension must match MeshMetadata dimension");

    if constexpr(DIM == 2)
    {
      return axom::primal::BoundingBox<double, DIM>({bb_min[0], bb_min[1]}, {bb_max[0], bb_max[1]});
    }
    else
    {
      return axom::primal::BoundingBox<double, DIM>({bb_min[0], bb_min[1], bb_min[2]},
                                                    {bb_max[0], bb_max[1], bb_max[2]});
    }
  }

  /// Defines the input schema for our mesh metadata, w/ some validation checks
  static void defineSchema(inlet::Container& mesh_schema)
  {
    mesh_schema.addInt("dim", "Dimension (2 or 3)").required().range(2, 3);

    auto& bb = mesh_schema.addStruct("bounding_box", "Mesh bounding box").required();

    auto& min = bb.addStruct("min", "Minimum coordinates").required();
    min.addDouble("x", "Minimum x coordinate").required();
    min.addDouble("y", "Minimum y coordinate").required();
    min.addDouble("z", "Minimum z coordinate (only specify when dim is 3)");

    auto& max = bb.addStruct("max", "Maximum coordinates").required();
    max.addDouble("x", "Maximum x coordinate").required();
    max.addDouble("y", "Maximum y coordinate").required();
    max.addDouble("z", "Maximum z coordinate (only specify when dim is 3)");

    auto& res = mesh_schema.addStruct("resolution", "Mesh resolution").required();
    res.addInt("x", "Resolution in x direction").required().range(1, std::numeric_limits<int>::max());
    res.addInt("y", "Resolution in y direction").required().range(1, std::numeric_limits<int>::max());
    res.addInt("z", "Resolution in z direction (only specify when dim is 3)")
      .range(1, std::numeric_limits<int>::max());

    mesh_schema.addString("background_material", "Optional background material");
    mesh_schema.addInt("volume_fraction_order", "Order for volume fraction fields (>= 1)")
      .range(1, std::numeric_limits<int>::max());
    mesh_schema.addInt("mesh_order", "Order for mesh nodes (>= 1)")
      .range(1, std::numeric_limits<int>::max());
    mesh_schema.addInt("quadrature_order", "Order for quadrature (>= 1)")
      .range(1, std::numeric_limits<int>::max());

    mesh_schema.addString("sampling_method", "Sampling method ('inout' or 'winding')")
      .validValues({"inout", "winding"});

    // Validate bounding box min/max
    bb.registerVerifier([](const inlet::Container& input) {
      bool valid = true;
      for(const std::string axis : {"x", "y", "z"})
      {
        const std::string min_str = axom::fmt::format("min/{}", axis);
        const std::string max_str = axom::fmt::format("max/{}", axis);
        if(axis == "z" && (!input.contains(min_str) && !input.contains(max_str)))
        {
          continue;
        }

        if(const double min_val = input[min_str], max_val = input[max_str]; min_val >= max_val)
        {
          SLIC_WARNING(axom::fmt::format("Invalid bounding box range for {}-coordinate: {} >= {}",
                                         axis,
                                         min_val,
                                         max_val));
          valid = false;
        }
      }
      return valid;
    });

    // Validate presence/absence of z fields based on dim
    mesh_schema.registerVerifier([](const inlet::Container& input) {
      const int dim = input["dim"];
      bool valid = true;

      for(const auto& field : {"bounding_box/min/z", "bounding_box/max/z", "resolution/z"})
      {
        if(dim == 3)
        {
          if(!input.contains(field))
          {
            SLIC_WARNING(
              axom::fmt::format("Z-coordinate for '{}' is required when dimension is 3", field));
            valid = false;
          }
        }
        else if(dim == 2)
        {
          if(input.contains(field))
          {
            SLIC_WARNING(
              axom::fmt::format("Z-coordinate for '{}' should not be provided when dimension is 2",
                                field));
            valid = false;
          }
        }
      }

      return valid;
    });
  }
};

/// Generate a MeshMetadata instance from the in-memory inlet representation
template <>
struct FromInlet<MeshMetadata>
{
  MeshMetadata operator()(const inlet::Container& input_data)
  {
    MeshMetadata result;

    result.dim = input_data["dim"];

    result.bb_min.resize(result.dim);
    result.bb_max.resize(result.dim);
    result.resolution.resize(result.dim);

    auto bb = input_data["bounding_box"];
    result.bb_min[0] = bb["min/x"];
    result.bb_min[1] = bb["min/y"];

    result.bb_max[0] = bb["max/x"];
    result.bb_max[1] = bb["max/y"];

    auto res = input_data["resolution"];
    result.resolution[0] = res["x"];
    result.resolution[1] = res["y"];

    if(result.dim == 3)
    {
      result.bb_min[2] = bb["min/z"];
      result.bb_max[2] = bb["max/z"];
      result.resolution[2] = res["z"];
    }

    if(input_data.contains("background_material"))
    {
      result.background_material = static_cast<std::string>(input_data["background_material"]);
    }

    if(input_data.contains("volume_fraction_order"))
    {
      result.volume_fraction_order = static_cast<int>(input_data["volume_fraction_order"]);
    }

    if(input_data.contains("quadrature_order"))
    {
      result.quadrature_order = static_cast<int>(input_data["quadrature_order"]);
    }

    if(input_data.contains("sampling_method"))
    {
      const auto str = static_cast<std::string>(input_data["sampling_method"]);
      if(str == "inout")
      {
        result.sampling_method = quest::SamplingShaper::SamplingMethod::InOut;
      }
      else if(str == "winding")
      {
        result.sampling_method = quest::SamplingShaper::SamplingMethod::WindingNumber;
      }
    }

    return result;
  }
};

/// Helper function to generate a 2D or 3D mfem cartesian mesh (partitioned across ranks in MPI configurations)
mfem::Mesh* createCartesianMesh(const MeshMetadata& meta, int nodal_order)
{
  mfem::Mesh* mesh = nullptr;

  switch(meta.dim)
  {
  case 2:
  {
    const axom::NumericArray<int, 2> res {meta.resolution[0], meta.resolution[1]};
    const auto bbox = meta.getBoundingBox<2>();

    SLIC_INFO_ROOT(axom::fmt::format("Creating 2D Cartesian mesh of res {} and bbox {}", res, bbox));
    mesh = quest::util::make_cartesian_mfem_mesh_2D(bbox, res, nodal_order);
  }
  break;
  case 3:
  {
    const axom::NumericArray<int, 3> res {meta.resolution[0], meta.resolution[1], meta.resolution[2]};
    const auto bbox = meta.getBoundingBox<3>();

    SLIC_INFO_ROOT(axom::fmt::format("Creating 3D Cartesian mesh of res {} and bbox {}", res, bbox));
    mesh = quest::util::make_cartesian_mfem_mesh_3D(bbox, res, nodal_order);
  }
  break;
  default:
    SLIC_ERROR("Only 2D and 3D meshes are supported");
    break;
  }

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

int main(int argc, char** argv)
{
  // use Axom's utility classes to initialize MPI and a logger
  // these will be automatically finalized at the end of the function
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);
  axom::slic::SimpleLogger raii_logger;
  axom::slic::setIsRoot(mpi_raii_wrapper.my_rank() == 0);

  // --------------------------------------------------------------------------
  // CLI for input files
  // --------------------------------------------------------------------------
  axom::CLI::App app {"Shaping pipeline using separate Inlet mesh metadata and Klee shapes"};
  std::string inputFilename;  // Mesh metadata Inlet Lua
  std::string kleeFilename;   // Klee shape set YAML
  bool verbose = false;

  app.add_option("-m,--mesh_file", inputFilename)
    ->description("Mesh metadata Inlet Lua file")
    ->required()
    ->check(axom::CLI::ExistingFile);
  app.add_option("-k,--klee_file", kleeFilename)
    ->description("Klee shape set YAML file")
    ->required()
    ->check(axom::CLI::ExistingFile);
  app.add_flag("-v,--verbose", verbose)->description("Enable verbose (debug) logging");

  try
  {
    app.parse(argc, argv);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = app.exit(e);
#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    return retval;
  }

  slic::setLoggingMsgLevel(verbose ? slic::message::Debug : slic::message::Info);

  // --------------------------------------------------------------------------
  // Parse and validate inlet input into MeshMetadata
  // --------------------------------------------------------------------------
  MeshMetadata meta = [&]() -> MeshMetadata {
    std::unique_ptr<inlet::Reader> reader = std::make_unique<inlet::LuaReader>();
    reader->parseFile(inputFilename);
    inlet::Inlet input(std::move(reader));

    auto& mesh_schema = input.addStruct("mesh", "Mesh metadata").required();
    MeshMetadata::defineSchema(mesh_schema);

    SLIC_ERROR_IF(!input.verify(), "Input validation failed.");

    // Parse mesh metadata and shaping params
    return input["mesh"].get<MeshMetadata>();
  }();

  // --------------------------------------------------------------------------
  // Set up computational mesh from MeshMetadata, store within a sidre-basead
  // DataCollection following the mesh blueprint conventions
  // --------------------------------------------------------------------------
  mfem::Mesh* mesh = createCartesianMesh(meta, meta.mesh_order);
  constexpr bool dc_owns_data = true;  // Note: dc takes ownership of mesh
  sidre::MFEMSidreDataCollection dc("shaping", nullptr, dc_owns_data);
#ifdef AXOM_USE_MPI
  dc.SetMesh(MPI_COMM_WORLD, mesh);
#else
  dc.SetMesh(mesh);
#endif
  dc.SetMeshNodesName("positions");

  // Associate any fields that begin with "vol_frac" with "material" so when
  // the data collection is written, a matset will be created.
  dc.AssociateMaterialSet("vol_frac", "material");

  // --------------------------------------------------------------------------
  // Load and validate Klee shape set
  // --------------------------------------------------------------------------
  klee::ShapeSet shapeSet;
  try
  {
    shapeSet = klee::readShapeSet(kleeFilename);
  }
  catch(klee::KleeError& error)
  {
    std::vector<std::string> errs;
    for(const auto& verificationError : error.getErrors())
    {
      errs.push_back(axom::fmt::format(" - '{}': {}",
                                       static_cast<std::string>(verificationError.path),
                                       verificationError.message));
    }
    SLIC_WARNING(axom::fmt::format("Error parsing klee input:\n{}", axom::fmt::join(errs, "\n")));
    return 1;
  }

  // --------------------------------------------------------------------------
  // Setup sample-based shaper
  // --------------------------------------------------------------------------
  using RuntimePolicy = axom::runtime_policy::Policy;
  RuntimePolicy policy = RuntimePolicy::seq;

  auto shaper = std::make_unique<quest::SamplingShaper>(policy,
                                                        axom::policyToDefaultAllocatorID(policy),
                                                        shapeSet,
                                                        &dc);
  shaper->setVerbosity(verbose);
  shaper->setQuadratureOrder(meta.quadrature_order);
  shaper->setVolumeFractionOrder(meta.volume_fraction_order);

  shaper->setSamplingMethod(meta.sampling_method);
  if(meta.sampling_method == quest::SamplingShaper::SamplingMethod::InOut)
  {
    shaper->setSamplesPerKnotSpan(50);
  }

  // Project initial volume fractions, if applicable
  if(!meta.background_material.empty())
  {
    std::map<std::string, mfem::GridFunction*> initial_grid_functions;

    // Generate a background material (volume fraction set to 1) if provided
    auto material = meta.background_material;
    auto name = axom::fmt::format("vol_frac_{}", material);

    const int order = meta.volume_fraction_order;
    const int dim = meta.dim;
    const auto basis = mfem::BasisType::Positive;

    auto* coll = new mfem::L2_FECollection(order, dim, basis);
    auto* fes = new mfem::FiniteElementSpace(dc.GetMesh(), coll);
    const int sz = fes->GetVSize();

    auto* view = dc.AllocNamedBuffer(name, sz);
    auto* volFrac = new mfem::GridFunction(fes, view->getArray());
    volFrac->MakeOwner(coll);

    (*volFrac) = 1.0;

    dc.RegisterField(name, volFrac);

    initial_grid_functions[material] = dc.GetField(name);
    // Inform the shaper about any initial volume fraction fields
    shaper->importInitialVolumeFractions(initial_grid_functions);
  }

  // --------------------------------------------------------------------------
  // Run shaping pipeline
  // --------------------------------------------------------------------------
  SLIC_INFO_ROOT(axom::fmt::format("{:=^80}", "Shaping"));
  for(const auto& shape : shapeSet.getShapes())
  {
    const std::string shapeFormat = shape.getGeometry().getFormat();
    SLIC_INFO_ROOT(
      axom::fmt::format("{:-^80}",
                        axom::fmt::format("Processing shape '{}' of material '{}' (format '{}')",
                                          shape.getName(),
                                          shape.getMaterial(),
                                          shapeFormat)));

    const klee::Dimensions shapeDim = shape.getGeometry().getInputDimensions();

    shaper->loadShape(shape);
    shaper->prepareShapeQuery(shapeDim, shape);
    shaper->runShapeQuery(shape);
    shaper->applyReplacementRules(shape);
    shaper->finalizeShapeQuery();
    slic::flushStreams();
  }

  //---------------------------------------------------------------------------
  // After processing all shapes, generate/adjust the material volume fractions
  //---------------------------------------------------------------------------
  SLIC_INFO_ROOT(axom::fmt::format("{:=^80}", "Generating volume fraction fields for materials"));

  shaper->adjustVolumeFractions();

  // --------------------------------------------------------------------------
  // Write data collection to disk
  // --------------------------------------------------------------------------
#ifdef MFEM_USE_MPI
  {
    dc.Save();
  }
#endif

  return 0;
}
