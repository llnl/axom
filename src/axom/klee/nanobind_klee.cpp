// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "axom/klee/Geometry.hpp"
#include "axom/klee/KleeError.hpp"
#include "axom/klee/Shape.hpp"
#include "axom/klee/ShapeSet.hpp"
#include "axom/klee/Units.hpp"
#include "axom/klee/io/IO.hpp"
#include "axom/fmt.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace nb = nanobind;

namespace
{
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
}  // namespace

NB_MODULE(pyklee, m)
{
  m.doc() = "Tutorial-facing Python bindings for Axom Klee shape-set loading and inspection.";

  nb::enum_<axom::klee::Dimensions>(m, "Dimensions")
    .value("Unspecified", axom::klee::Dimensions::Unspecified)
    .value("Two", axom::klee::Dimensions::Two)
    .value("Three", axom::klee::Dimensions::Three)
    .export_values();

  nb::enum_<axom::klee::LengthUnit>(m, "LengthUnit")
    .value("km", axom::klee::LengthUnit::km)
    .value("m", axom::klee::LengthUnit::m)
    .value("dm", axom::klee::LengthUnit::dm)
    .value("cm", axom::klee::LengthUnit::cm)
    .value("mm", axom::klee::LengthUnit::mm)
    .value("um", axom::klee::LengthUnit::um)
    .value("nm", axom::klee::LengthUnit::nm)
    .value("angstrom", axom::klee::LengthUnit::angstrom)
    .value("miles", axom::klee::LengthUnit::miles)
    .value("feet", axom::klee::LengthUnit::feet)
    .value("inches", axom::klee::LengthUnit::inches)
    .value("mils", axom::klee::LengthUnit::mils)
    .value("unspecified", axom::klee::LengthUnit::unspecified)
    .export_values();

  nb::class_<axom::klee::TransformableGeometryProperties>(m, "TransformableGeometryProperties")
    .def(nb::init<>())
    .def_rw("dimensions", &axom::klee::TransformableGeometryProperties::dimensions)
    .def_rw("units", &axom::klee::TransformableGeometryProperties::units);

  nb::class_<axom::klee::Geometry>(m, "Geometry")
    .def("getFormat", &axom::klee::Geometry::getFormat)
    .def("getPath", &axom::klee::Geometry::getPath)
    .def("getInputDimensions", &axom::klee::Geometry::getInputDimensions)
    .def("getOutputDimensions", &axom::klee::Geometry::getOutputDimensions)
    .def("getStartProperties",
         &axom::klee::Geometry::getStartProperties,
         nb::rv_policy::reference_internal)
    .def("getEndProperties", &axom::klee::Geometry::getEndProperties)
    .def("hasGeometry", &axom::klee::Geometry::hasGeometry);

  nb::class_<axom::klee::Shape>(m, "Shape")
    .def("getName", &axom::klee::Shape::getName)
    .def("getMaterial", &axom::klee::Shape::getMaterial)
    .def("replaces", &axom::klee::Shape::replaces, nb::arg("material"))
    .def("getGeometry", &axom::klee::Shape::getGeometry, nb::rv_policy::reference_internal)
    .def("getMaterialsReplaced", &axom::klee::Shape::getMaterialsReplaced)
    .def("getMaterialsNotReplaced", &axom::klee::Shape::getMaterialsNotReplaced);

  nb::class_<axom::klee::ShapeSet>(m, "ShapeSet")
    .def("getShapes", &axom::klee::ShapeSet::getShapes, nb::rv_policy::reference_internal)
    .def("getPath", &axom::klee::ShapeSet::getPath)
    .def("getDimensions", &axom::klee::ShapeSet::getDimensions);

  m.def("readShapeSet", &readShapeSet, nb::arg("path"), "Load and validate a Klee shape-set file.");
}
