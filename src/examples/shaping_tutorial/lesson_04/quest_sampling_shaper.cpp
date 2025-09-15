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

// Mesh metadata (from examples, adapted)
struct MeshMetadata
{
  int dim;
  axom::Array<double> bb_min;
  axom::Array<double> bb_max;
  axom::Array<int> resolution;

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

    // Validate bounding box min/max
    bb.registerVerifier([](const inlet::Container& input) {
      bool valid = true;
      for(const std::string& axis : {"x", "y", "z"})
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

    return result;
  }
};

mfem::Mesh* createCartesianMesh(const MeshMetadata& meta, int outputOrder = 1)
{
  mfem::Mesh* mesh = nullptr;

  switch(meta.dim)
  {
  case 2:
  {
    const axom::NumericArray<int, 2> res {meta.resolution[0], meta.resolution[1]};
    const auto bbox = meta.getBoundingBox<2>();

    SLIC_INFO(axom::fmt::format("Creating 2D Cartesian mesh of res {} and bbox {}", res, bbox));
    mesh = quest::util::make_cartesian_mfem_mesh_2D(bbox, res, outputOrder);
  }
  break;
  case 3:
  {
    const axom::NumericArray<int, 3> res {meta.resolution[0], meta.resolution[1], meta.resolution[2]};
    const auto bbox = meta.getBoundingBox<3>();

    SLIC_INFO(axom::fmt::format("Creating 3D Cartesian mesh of res {} and bbox {}", res, bbox));
    mesh = quest::util::make_cartesian_mfem_mesh_3D(bbox, res, outputOrder);
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
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);
  axom::slic::SimpleLogger raii_logger;

  // --------------------------------------------------------------------------
  // CLI for input files
  // --------------------------------------------------------------------------
  axom::CLI::App app {"Shaping pipeline using separate Inlet mesh metadata and Klee shapes"};
  std::string inputFilename;  // Mesh metadata Inlet Lua
  std::string kleeFilename;   // Klee shape set YAML
  bool verbose = false;
  int outputOrder = 1;

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

  slic::setAbortOnWarning(true);
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
  // Set up computational mesh from MeshMetadata
  // --------------------------------------------------------------------------
  mfem::Mesh* mesh = createCartesianMesh(meta, outputOrder);
  constexpr bool dc_owns_data = true;  // Note: dc takes ownership of mesh
  sidre::MFEMSidreDataCollection dc("shaping", nullptr, dc_owns_data);
#ifdef AXOM_USE_MPI
  dc.SetMesh(MPI_COMM_WORLD, mesh);
#else
  dc.SetMesh(mesh);
#endif
  dc.SetMeshNodesName("positions");

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
  // Setup shaper
  // --------------------------------------------------------------------------
  using RuntimePolicy = axom::runtime_policy::Policy;
  RuntimePolicy policy = RuntimePolicy::seq;

  auto shaper = std::make_unique<quest::SamplingShaper>(policy,
                                                        axom::policyToDefaultAllocatorID(policy),
                                                        shapeSet,
                                                        &dc);
  // TODO: Set some additional parameters....

  // --------------------------------------------------------------------------
  // Run pipeline
  // --------------------------------------------------------------------------
  // Process each shape
  SLIC_INFO(axom::fmt::format("{:=^80}", "Shaping"));
  for(const auto& shape : shapeSet.getShapes())
  {
    const std::string shapeFormat = shape.getGeometry().getFormat();
    SLIC_INFO_ROOT(
      axom::fmt::format("{:-^80}",
                        axom::fmt::format("Processing shape '{}' of material '{}' (format '{}')",
                                          shape.getName(),
                                          shape.getMaterial(),
                                          shapeFormat)));

    shaper->loadShape(shape);
    shaper->prepareShapeQuery(shape.getGeometry().getInputDimensions(), shape);
    shaper->runShapeQuery(shape);
    shaper->applyReplacementRules(shape);
    shaper->finalizeShapeQuery();
    slic::flushStreams();
  }

  return 0;
}
