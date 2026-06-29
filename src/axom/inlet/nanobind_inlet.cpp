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

#include <limits>
#include <string>
#include <vector>

namespace nb = nanobind;

namespace
{ }  // namespace

NB_MODULE(pyinlet, m)
{
  m.doc() = "Tutorial-facing Python bindings for parsing shaping tutorial inputs with Axom Inlet.";

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
}
