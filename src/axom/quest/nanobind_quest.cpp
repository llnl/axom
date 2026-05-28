// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "axom/config.hpp"
#include "axom/core/utilities/RAII.hpp"
#include "axom/fmt.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/YAMLReader.hpp"
#ifdef AXOM_USE_LUA
  #include "axom/inlet/LuaReader.hpp"
#endif
#include "axom/klee/KleeError.hpp"
#include "axom/klee/io/IO.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/quest/SamplingShaper.hpp"
#include "axom/quest/util/mesh_helpers.hpp"
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"
#include "axom/slic.hpp"

#include "mfem.hpp"

#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace nb = nanobind;

namespace
{
using BoundingBox2D = axom::primal::BoundingBox<double, 2>;
using BoundingBox3D = axom::primal::BoundingBox<double, 3>;
using Point2D = axom::primal::Point<double, 2>;
using Point3D = axom::primal::Point<double, 3>;

struct MeshMetadata
{
  int dim {2};
  std::vector<double> bb_min {0.0, 0.0};
  std::vector<double> bb_max {1.0, 1.0};
  std::vector<int> resolution {10, 10};

  std::string background_material;
  int volume_fraction_order {2};
  int mesh_order {1};
  int quadrature_order {5};
  std::string sampling_method {"inout"};

  BoundingBox2D getBoundingBox2D() const
  {
    return BoundingBox2D(Point2D {bb_min[0], bb_min[1]}, Point2D {bb_max[0], bb_max[1]});
  }

  BoundingBox3D getBoundingBox3D() const
  {
    return BoundingBox3D(Point3D {bb_min[0], bb_min[1], bb_min[2]},
                         Point3D {bb_max[0], bb_max[1], bb_max[2]});
  }
};

bool hasSuffix(const std::string& path, const std::string& suffix)
{
  return path.size() >= suffix.size() &&
    path.compare(path.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::unique_ptr<axom::inlet::Reader> makeReader(const std::string& path)
{
  if(hasSuffix(path, ".yaml") || hasSuffix(path, ".yml"))
  {
    return std::make_unique<axom::inlet::YAMLReader>();
  }

#ifdef AXOM_USE_LUA
  if(hasSuffix(path, ".lua"))
  {
    return std::make_unique<axom::inlet::LuaReader>();
  }
#endif

  throw std::runtime_error(
    axom::fmt::format("Unsupported mesh metadata extension for '{}'. Expected .yaml, .yml or .lua.",
                      path));
}

void defineMeshSchema(axom::inlet::Container& mesh_schema)
{
  mesh_schema.addInt("dim", "Dimension (2 or 3)").range(2, 3);

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
}

MeshMetadata loadMeshMetadata(const std::string& path)
{
  auto reader = makeReader(path);
  if(!reader->parseFile(path))
  {
    throw std::runtime_error(axom::fmt::format("Failed to parse '{}'.", path));
  }

  axom::inlet::Inlet inlet(std::move(reader));
  auto& mesh_schema = inlet.addStruct("mesh", "Mesh metadata").required();
  defineMeshSchema(mesh_schema);

  std::vector<axom::inlet::VerificationError> errors;
  inlet.verify(&errors);
  if(!errors.empty())
  {
    std::vector<std::string> messages;
    for(const auto& err : errors)
    {
      messages.push_back(axom::fmt::format("{}: {}", static_cast<std::string>(err.path), err.message));
    }
    throw std::runtime_error(
      axom::fmt::format("Mesh metadata validation failed:\n{}", axom::fmt::join(messages, "\n")));
  }

  MeshMetadata result;
  auto mesh = inlet["mesh"];

  if(mesh.contains("dim"))
  {
    result.dim = static_cast<int>(mesh["dim"]);
  }

  result.bb_min.resize(result.dim);
  result.bb_max.resize(result.dim);
  result.resolution.resize(result.dim);

  auto bb = mesh["bounding_box"];
  result.bb_min[0] = bb["min/x"];
  result.bb_min[1] = bb["min/y"];
  result.bb_max[0] = bb["max/x"];
  result.bb_max[1] = bb["max/y"];

  auto res = mesh["resolution"];
  result.resolution[0] = res["x"];
  result.resolution[1] = res["y"];

  if(result.dim == 3)
  {
    result.bb_min[2] = bb["min/z"];
    result.bb_max[2] = bb["max/z"];
    result.resolution[2] = res["z"];
  }

  if(mesh.contains("background_material"))
  {
    result.background_material = static_cast<std::string>(mesh["background_material"]);
  }
  if(mesh.contains("volume_fraction_order"))
  {
    result.volume_fraction_order = static_cast<int>(mesh["volume_fraction_order"]);
  }
  if(mesh.contains("mesh_order"))
  {
    result.mesh_order = static_cast<int>(mesh["mesh_order"]);
  }
  if(mesh.contains("quadrature_order"))
  {
    result.quadrature_order = static_cast<int>(mesh["quadrature_order"]);
  }
  if(mesh.contains("sampling_method"))
  {
    result.sampling_method = static_cast<std::string>(mesh["sampling_method"]);
  }

  return result;
}

mfem::Mesh* createCartesianMesh(const MeshMetadata& meta)
{
  mfem::Mesh* mesh = nullptr;

  if(meta.dim == 2)
  {
    const axom::NumericArray<int, 2> res {meta.resolution[0], meta.resolution[1]};
    mesh =
      axom::quest::util::make_cartesian_mfem_mesh_2D(meta.getBoundingBox2D(), res, meta.mesh_order);
  }
  else if(meta.dim == 3)
  {
    const axom::NumericArray<int, 3> res {meta.resolution[0], meta.resolution[1], meta.resolution[2]};
    mesh =
      axom::quest::util::make_cartesian_mfem_mesh_3D(meta.getBoundingBox3D(), res, meta.mesh_order);
  }
  else
  {
    throw std::runtime_error("Only 2D and 3D meshes are supported.");
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

axom::klee::ShapeSet readShapeSet(const std::string& path)
{
  try
  {
    return axom::klee::readShapeSet(path);
  }
  catch(const axom::klee::KleeError& error)
  {
    std::vector<std::string> errs;
    for(const auto& verificationError : error.getErrors())
    {
      errs.push_back(axom::fmt::format("{}: {}",
                                       static_cast<std::string>(verificationError.path),
                                       verificationError.message));
    }
    throw std::runtime_error(
      axom::fmt::format("Error parsing klee input:\n{}", axom::fmt::join(errs, "\n")));
  }
}

void ensureMPIInitialized()
{
#ifdef AXOM_USE_MPI
  static std::unique_ptr<axom::utilities::raii::MPIWrapper> mpi_wrapper;
  if(!mpi_wrapper)
  {
    int argc = 1;
    char program_name[] = "pyquest";
    char* argv[] = {program_name, nullptr};
    mpi_wrapper = std::make_unique<axom::utilities::raii::MPIWrapper>(argc, argv);
  }
#endif
}

struct ShapingResult
{
  explicit ShapingResult(std::unique_ptr<axom::sidre::MFEMSidreDataCollection> dc)
    : collection(std::move(dc))
  { }

  void save() { collection->Save(); }

  axom::sidre::Group* getBlueprintGroup() { return collection->GetBPGroup(); }

  axom::sidre::DataStore* getDataStore() { return collection->GetBPGroup()->getDataStore(); }

  std::string getCollectionName() const { return collection->GetCollectionName(); }

  std::unique_ptr<axom::sidre::MFEMSidreDataCollection> collection;
};

ShapingResult runShaping(const std::string& meshFile,
                         const std::string& kleeFile,
                         const std::string& outputName,
                         bool verbose)
{
  ensureMPIInitialized();

  axom::slic::SimpleLogger logger;
  axom::slic::setLoggingMsgLevel(verbose ? axom::slic::message::Debug : axom::slic::message::Info);

  const MeshMetadata meta = loadMeshMetadata(meshFile);
  auto shapeSet = readShapeSet(kleeFile);

  auto* mesh = createCartesianMesh(meta);
  auto dc = std::make_unique<axom::sidre::MFEMSidreDataCollection>(outputName, nullptr, true);
#if defined(AXOM_USE_MPI)
  dc->SetMesh(MPI_COMM_WORLD, mesh);
#else
  dc->SetMesh(mesh);
#endif
  dc->SetMeshNodesName("positions");
  dc->AssociateMaterialSet("vol_frac", "material");

  using RuntimePolicy = axom::runtime_policy::Policy;
  const RuntimePolicy policy = RuntimePolicy::seq;

  auto shaper =
    std::make_unique<axom::quest::SamplingShaper>(policy,
                                                  axom::policyToDefaultAllocatorID(policy),
                                                  shapeSet,
                                                  dc.get());

  shaper->setVerbosity(verbose);
  shaper->setQuadratureOrder(meta.quadrature_order);
  shaper->setVolumeFractionOrder(meta.volume_fraction_order);

  if(meta.sampling_method == "winding")
  {
    shaper->setSamplingMethod(axom::quest::SamplingShaper::SamplingMethod::WindingNumber);
  }
  else
  {
    shaper->setSamplingMethod(axom::quest::SamplingShaper::SamplingMethod::InOut);
    shaper->setSamplesPerKnotSpan(50);
  }

  if(!meta.background_material.empty())
  {
    std::map<std::string, mfem::GridFunction*> initial_grid_functions;

    const auto material = meta.background_material;
    const auto name = axom::fmt::format("vol_frac_{}", material);
    const int order = meta.volume_fraction_order;
    const int dim = meta.dim;
    const auto basis = mfem::BasisType::Positive;

    auto* coll = new mfem::L2_FECollection(order, dim, basis);
    auto* fes = new mfem::FiniteElementSpace(dc->GetMesh(), coll);
    const int sz = fes->GetVSize();

    auto* view = dc->AllocNamedBuffer(name, sz);
    auto* volFrac = new mfem::GridFunction(fes, view->getArray());
    volFrac->MakeOwner(coll);
    (*volFrac) = 1.0;

    dc->RegisterField(name, volFrac);
    initial_grid_functions[material] = dc->GetField(name);
    shaper->importInitialVolumeFractions(initial_grid_functions);
  }

  for(const auto& shape : shapeSet.getShapes())
  {
    const auto shapeDim = shape.getGeometry().getInputDimensions();
    shaper->loadShape(shape);
    shaper->prepareShapeQuery(shapeDim, shape);
    shaper->runShapeQuery(shape);
    shaper->applyReplacementRules(shape);
    shaper->finalizeShapeQuery();
  }

  shaper->adjustVolumeFractions();
  return ShapingResult(std::move(dc));
}
}  // namespace

NB_MODULE(pyquest, m)
{
  m.doc() = "Tutorial-facing Python bindings for Axom Quest shaping.";

  nb::module_::import_("pysidre");

  nb::enum_<axom::quest::SamplingShaper::SamplingMethod>(m, "SamplingMethod")
    .value("InOut", axom::quest::SamplingShaper::SamplingMethod::InOut)
    .value("WindingNumber", axom::quest::SamplingShaper::SamplingMethod::WindingNumber)
    .export_values();

  nb::class_<ShapingResult>(m, "ShapingResult")
    .def("save", &ShapingResult::save)
    .def("getBlueprintGroup", &ShapingResult::getBlueprintGroup, nb::rv_policy::reference_internal)
    .def("getDataStore", &ShapingResult::getDataStore, nb::rv_policy::reference_internal)
    .def("getCollectionName", &ShapingResult::getCollectionName);

  m.def("runShaping",
        &runShaping,
        nb::arg("mesh_file"),
        nb::arg("klee_file"),
        nb::arg("output_name") = "shaping",
        nb::arg("verbose") = false,
        "Run the shaping tutorial pipeline and return an owned result object.");
}
