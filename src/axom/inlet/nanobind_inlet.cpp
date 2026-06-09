// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/Container.hpp"
#include "axom/inlet/Field.hpp"
#include "axom/inlet/Proxy.hpp"
#include "axom/inlet/Reader.hpp"
#include "axom/inlet/VerifiableScalar.hpp"
#include "axom/inlet/YAMLReader.hpp"
#ifdef AXOM_USE_LUA
  #include "axom/inlet/LuaReader.hpp"
#endif
#include "axom/fmt.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

#include <limits>
#include <memory>
#include <sstream>
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
    if(dim != 2)
    {
      throw std::runtime_error("MeshMetadata does not describe a 2D mesh.");
    }

    return BoundingBox2D(Point2D {bb_min[0], bb_min[1]}, Point2D {bb_max[0], bb_max[1]});
  }

  BoundingBox3D getBoundingBox3D() const
  {
    if(dim != 3)
    {
      throw std::runtime_error("MeshMetadata does not describe a 3D mesh.");
    }

    return BoundingBox3D(Point3D {bb_min[0], bb_min[1], bb_min[2]},
                         Point3D {bb_max[0], bb_max[1], bb_max[2]});
  }
};

bool hasSuffix(const std::string& path, const std::string& suffix)
{
  return path.size() >= suffix.size() &&
    path.compare(path.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::unique_ptr<axom::inlet::Reader> makeReader(const std::string& path);

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

  bb.registerVerifier([](const axom::inlet::Container& input) {
    bool valid = true;
    for(const char* axis : {"x", "y", "z"})
    {
      const std::string min_str = axom::fmt::format("min/{}", axis);
      const std::string max_str = axom::fmt::format("max/{}", axis);
      if(std::string(axis) == "z" && (!input.contains(min_str) || !input.contains(max_str)))
      {
        continue;
      }

      const double min_val = input[min_str];
      const double max_val = input[max_str];
      if(min_val >= max_val)
      {
        valid = false;
      }
    }
    return valid;
  });

  mesh_schema.registerVerifier([](const axom::inlet::Container& input) {
    const int dim = input.contains("dim") ? static_cast<int>(input["dim"]) : 2;
    bool valid = true;

    for(const auto& field : {"bounding_box/min/z", "bounding_box/max/z", "resolution/z"})
    {
      if(dim == 3)
      {
        valid = valid && input.contains(field);
      }
      else if(dim == 2)
      {
        valid = valid && !input.contains(field);
      }
    }

    return valid;
  });
}

void defineValidatedMeshSchema(axom::inlet::Container& mesh_schema)
{
  auto& bb = mesh_schema.addStruct("bounding_box", "Mesh bounding box").required();

  auto& min = bb.addStruct("min", "Minimum coordinates").required();
  min.addDouble("x", "Minimum x coordinate").required();
  min.addDouble("y", "Minimum y coordinate").required();

  auto& max = bb.addStruct("max", "Maximum coordinates").required();
  max.addDouble("x", "Maximum x coordinate").required();
  max.addDouble("y", "Maximum y coordinate").required();

  auto& res = mesh_schema.addStruct("resolution", "Mesh resolution").required();
  res.addInt("x", "Resolution in x direction").required().range(1, std::numeric_limits<int>::max());
  res.addInt("y", "Resolution in y direction").required().range(1, std::numeric_limits<int>::max());

  bb.registerVerifier([](const axom::inlet::Container& input) {
    const double min_x = input["min/x"];
    const double max_x = input["max/x"];
    const double min_y = input["min/y"];
    const double max_y = input["max/y"];
    return (min_x < max_x) && (min_y < max_y);
  });
}

MeshMetadata meshMetadataFromProxy(const axom::inlet::Proxy& mesh)
{
  MeshMetadata result;

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

std::vector<std::string> validateMeshMetadata(const std::string& path)
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

  std::vector<std::string> messages;
  messages.reserve(errors.size());
  for(const auto& err : errors)
  {
    messages.push_back(axom::fmt::format("{}: {}", static_cast<std::string>(err.path), err.message));
  }

  return messages;
}

MeshMetadata loadMeshMetadata(const std::string& path)
{
  auto errors = validateMeshMetadata(path);
  if(!errors.empty())
  {
    throw std::runtime_error(
      axom::fmt::format("Mesh metadata validation failed:\n{}", axom::fmt::join(errors, "\n")));
  }

  auto reader = makeReader(path);
  if(!reader->parseFile(path))
  {
    throw std::runtime_error(axom::fmt::format("Failed to parse '{}'.", path));
  }

  axom::inlet::Inlet inlet(std::move(reader));
  auto& mesh_schema = inlet.addStruct("mesh", "Mesh metadata").required();
  defineMeshSchema(mesh_schema);
  return meshMetadataFromProxy(inlet["mesh"]);
}
}  // namespace

NB_MODULE(pyinlet, m)
{
  m.doc() = "Tutorial-facing Python bindings for parsing shaping tutorial inputs with Axom Inlet.";

  nb::module_::import_("pyprimal");

  nb::class_<axom::inlet::Reader>(m, "Reader");

  nb::class_<axom::inlet::YAMLReader, axom::inlet::Reader>(m, "YAMLReader")
    .def(nb::new_([]() { return new axom::inlet::YAMLReader(); }))
    .def("parseFile", &axom::inlet::YAMLReader::parseFile, nb::arg("file_path"))
    .def("parseString", &axom::inlet::YAMLReader::parseString, nb::arg("input_string"));

#ifdef AXOM_USE_LUA
  nb::class_<axom::inlet::LuaReader, axom::inlet::Reader>(m, "LuaReader")
    .def(nb::new_([]() { return new axom::inlet::LuaReader(); }))
    .def("parseFile", &axom::inlet::LuaReader::parseFile, nb::arg("file_path"))
    .def("parseString", &axom::inlet::LuaReader::parseString, nb::arg("input_string"));
#endif

  nb::class_<axom::inlet::Proxy>(m, "Proxy")
    .def("contains", &axom::inlet::Proxy::contains, nb::arg("name"))
    .def("__getitem__",
         [](const axom::inlet::Proxy& self, const std::string& name) { return self[name]; })
    .def("getMeshMetadata", &meshMetadataFromProxy)
    .def("__bool__",
         [](const axom::inlet::Proxy& self) { return static_cast<bool>(self.get<bool>()); })
    .def("__int__", [](const axom::inlet::Proxy& self) { return static_cast<int>(self.get<int>()); })
    .def("__float__",
         [](const axom::inlet::Proxy& self) { return static_cast<double>(self.get<double>()); })
    .def("__str__", [](const axom::inlet::Proxy& self) {
      return static_cast<std::string>(self.get<std::string>());
    });

  nb::class_<axom::inlet::VerifiableScalar>(m, "VerifiableScalar")
    .def("required",
         nb::overload_cast<bool>(&axom::inlet::VerifiableScalar::required),
         nb::arg("is_required") = true,
         nb::rv_policy::reference_internal)
    .def("isRequired", &axom::inlet::VerifiableScalar::isRequired)
    .def("defaultValue",
         nb::overload_cast<const std::string&>(&axom::inlet::VerifiableScalar::defaultValue),
         nb::arg("value"),
         nb::rv_policy::reference_internal)
    .def("defaultValue",
         nb::overload_cast<bool>(&axom::inlet::VerifiableScalar::defaultValue),
         nb::arg("value"),
         nb::rv_policy::reference_internal)
    .def("defaultValue",
         nb::overload_cast<int>(&axom::inlet::VerifiableScalar::defaultValue),
         nb::arg("value"),
         nb::rv_policy::reference_internal)
    .def("defaultValue",
         nb::overload_cast<double>(&axom::inlet::VerifiableScalar::defaultValue),
         nb::arg("value"),
         nb::rv_policy::reference_internal)
    .def("range",
         nb::overload_cast<int, int>(&axom::inlet::VerifiableScalar::range),
         nb::arg("start_val"),
         nb::arg("end_val"),
         nb::rv_policy::reference_internal)
    .def("range",
         nb::overload_cast<double, double>(&axom::inlet::VerifiableScalar::range),
         nb::arg("start_val"),
         nb::arg("end_val"),
         nb::rv_policy::reference_internal)
    .def("validValues",
         nb::overload_cast<const std::vector<int>&>(&axom::inlet::VerifiableScalar::validValues),
         nb::arg("values"),
         nb::rv_policy::reference_internal)
    .def("validValues",
         nb::overload_cast<const std::vector<double>&>(&axom::inlet::VerifiableScalar::validValues),
         nb::arg("values"),
         nb::rv_policy::reference_internal)
    .def(
      "validValues",
      nb::overload_cast<const std::vector<std::string>&>(&axom::inlet::VerifiableScalar::validValues),
      nb::arg("values"),
      nb::rv_policy::reference_internal)
    .def(
      "registerVerifier",
      [](axom::inlet::VerifiableScalar& self, nb::callable callback) -> axom::inlet::VerifiableScalar& {
        self.registerVerifier([callback](const axom::inlet::Field& field) {
          nb::gil_scoped_acquire gil;
          return nb::cast<bool>(callback(nb::cast(&field, nb::rv_policy::reference)));
        });
        return self;
      },
      nb::arg("verifier"),
      nb::rv_policy::reference_internal);

  nb::class_<axom::inlet::Field, axom::inlet::VerifiableScalar>(m, "Field")
    .def("required",
         nb::overload_cast<bool>(&axom::inlet::Field::required),
         nb::arg("is_required") = true,
         nb::rv_policy::reference_internal)
    .def("isRequired", &axom::inlet::Field::isRequired)
    .def("defaultValue",
         nb::overload_cast<const std::string&>(&axom::inlet::Field::defaultValue),
         nb::arg("value"),
         nb::rv_policy::reference_internal)
    .def("defaultValue",
         nb::overload_cast<bool>(&axom::inlet::Field::defaultValue),
         nb::arg("value"),
         nb::rv_policy::reference_internal)
    .def("defaultValue",
         nb::overload_cast<int>(&axom::inlet::Field::defaultValue),
         nb::arg("value"),
         nb::rv_policy::reference_internal)
    .def("defaultValue",
         nb::overload_cast<double>(&axom::inlet::Field::defaultValue),
         nb::arg("value"),
         nb::rv_policy::reference_internal)
    .def("range",
         nb::overload_cast<int, int>(&axom::inlet::Field::range),
         nb::arg("start_val"),
         nb::arg("end_val"),
         nb::rv_policy::reference_internal)
    .def("range",
         nb::overload_cast<double, double>(&axom::inlet::Field::range),
         nb::arg("start_val"),
         nb::arg("end_val"),
         nb::rv_policy::reference_internal)
    .def("validValues",
         nb::overload_cast<const std::vector<int>&>(&axom::inlet::Field::validValues),
         nb::arg("values"),
         nb::rv_policy::reference_internal)
    .def("validValues",
         nb::overload_cast<const std::vector<double>&>(&axom::inlet::Field::validValues),
         nb::arg("values"),
         nb::rv_policy::reference_internal)
    .def("validValues",
         nb::overload_cast<const std::vector<std::string>&>(&axom::inlet::Field::validValues),
         nb::arg("values"),
         nb::rv_policy::reference_internal);

  nb::class_<axom::inlet::Container>(m, "Container")
    .def(
      "addStruct",
      nb::overload_cast<const std::string&, const std::string&>(&axom::inlet::Container::addStruct),
      nb::arg("name"),
      nb::arg("description") = "",
      nb::rv_policy::reference_internal)
    .def("addInt",
         nb::overload_cast<const std::string&, const std::string&>(&axom::inlet::Container::addInt),
         nb::arg("name"),
         nb::arg("description") = "",
         nb::rv_policy::reference_internal)
    .def(
      "addDouble",
      nb::overload_cast<const std::string&, const std::string&>(&axom::inlet::Container::addDouble),
      nb::arg("name"),
      nb::arg("description") = "",
      nb::rv_policy::reference_internal)
    .def(
      "addString",
      nb::overload_cast<const std::string&, const std::string&>(&axom::inlet::Container::addString),
      nb::arg("name"),
      nb::arg("description") = "",
      nb::rv_policy::reference_internal)
    .def("required",
         nb::overload_cast<bool>(&axom::inlet::Container::required),
         nb::arg("is_required") = true,
         nb::rv_policy::reference_internal)
    .def("isRequired", &axom::inlet::Container::isRequired)
    .def(
      "registerVerifier",
      [](axom::inlet::Container& self, nb::callable callback) -> axom::inlet::Container& {
        self.registerVerifier([callback](const axom::inlet::Container& container,
                                         std::vector<axom::inlet::VerificationError>*) {
          nb::gil_scoped_acquire gil;
          return nb::cast<bool>(callback(nb::cast(&container, nb::rv_policy::reference)));
        });
        return self;
      },
      nb::arg("verifier"),
      nb::rv_policy::reference_internal)
    .def("verify", [](const axom::inlet::Container& self) { return self.verify(); })
    .def("contains", &axom::inlet::Container::contains, nb::arg("name"))
    .def("__getitem__",
         [](const axom::inlet::Container& self, const std::string& name) { return self[name]; });

  nb::class_<axom::inlet::Inlet>(m, "Inlet")
    .def(nb::init<std::unique_ptr<axom::inlet::Reader>, bool, bool>(),
         nb::arg("reader"),
         nb::arg("doc_enabled") = true,
         nb::arg("reconstruct") = false)
    .def("addStruct",
         nb::overload_cast<const std::string&, const std::string&>(&axom::inlet::Inlet::addStruct),
         nb::arg("name"),
         nb::arg("description") = "",
         nb::rv_policy::reference_internal)
    .def("verify", [](const axom::inlet::Inlet& self) { return self.verify(); })
    .def("contains", &axom::inlet::Inlet::contains, nb::arg("name"))
    .def("__getitem__",
         [](const axom::inlet::Inlet& self, const std::string& name) { return self[name]; });

  m.def(
    "open_inlet",
    [](const std::string& path, bool doc_enabled, bool reconstruct) {
      auto reader = makeReader(path);
      if(!reader->parseFile(path))
      {
        throw std::runtime_error(axom::fmt::format("Failed to parse '{}'.", path));
      }
      return axom::inlet::Inlet(std::move(reader), doc_enabled, reconstruct);
    },
    nb::arg("path"),
    nb::arg("doc_enabled") = true,
    nb::arg("reconstruct") = false);

  m.def("define_validated_mesh_schema",
        &defineValidatedMeshSchema,
        nb::arg("mesh_schema"),
        "Define the 2D validated mesh schema.");
  m.def("define_improved_mesh_schema",
        &defineMeshSchema,
        nb::arg("mesh_schema"),
        "Define the 2D/3D improved mesh schema.");

  nb::class_<MeshMetadata>(m, "MeshMetadata")
    .def(nb::init<>())
    .def_rw("dim", &MeshMetadata::dim)
    .def_rw("bb_min", &MeshMetadata::bb_min)
    .def_rw("bb_max", &MeshMetadata::bb_max)
    .def_rw("resolution", &MeshMetadata::resolution)
    .def_rw("background_material", &MeshMetadata::background_material)
    .def_rw("volume_fraction_order", &MeshMetadata::volume_fraction_order)
    .def_rw("mesh_order", &MeshMetadata::mesh_order)
    .def_rw("quadrature_order", &MeshMetadata::quadrature_order)
    .def_rw("sampling_method", &MeshMetadata::sampling_method)
    .def("getBoundingBox2D", &MeshMetadata::getBoundingBox2D)
    .def("getBoundingBox3D", &MeshMetadata::getBoundingBox3D)
    .def("__repr__", [](const MeshMetadata& self) {
      return axom::fmt::format(
        "MeshMetadata(dim={}, bb_min={}, bb_max={}, resolution={}, "
        "background_material='{}', volume_fraction_order={}, mesh_order={}, "
        "quadrature_order={}, sampling_method='{}')",
        self.dim,
        self.bb_min,
        self.bb_max,
        self.resolution,
        self.background_material,
        self.volume_fraction_order,
        self.mesh_order,
        self.quadrature_order,
        self.sampling_method);
    });

  m.def("validate_mesh_metadata",
        &validateMeshMetadata,
        nb::arg("path"),
        "Validate a shaping tutorial mesh metadata file and return the list of error messages.");
  m.def("load_mesh_metadata",
        &loadMeshMetadata,
        nb::arg("path"),
        "Load and validate a shaping tutorial mesh metadata file.");
}
