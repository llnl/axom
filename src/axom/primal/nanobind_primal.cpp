// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace nb = nanobind;

namespace axom
{
namespace primal
{
namespace
{
template <int NDIMS>
using BoundingBoxD = BoundingBox<double, NDIMS>;

template <int NDIMS>
using PointD = Point<double, NDIMS>;

template <int NDIMS>
PointD<NDIMS> vectorToPoint(const std::vector<double>& values, const char* which)
{
  if(static_cast<int>(values.size()) != NDIMS)
  {
    throw std::runtime_error(std::string(which) + " point must have exactly " +
                             std::to_string(NDIMS) + " coordinates.");
  }

  PointD<NDIMS> pt;
  for(int i = 0; i < NDIMS; ++i)
  {
    pt[i] = values[i];
  }
  return pt;
}

template <int NDIMS>
std::vector<double> pointToVector(const PointD<NDIMS>& pt)
{
  std::vector<double> values(NDIMS);
  for(int i = 0; i < NDIMS; ++i)
  {
    values[i] = pt[i];
  }
  return values;
}

template <int NDIMS>
void bindBoundingBox(nb::module_& m, const char* name)
{
  using Box = BoundingBoxD<NDIMS>;

  nb::class_<Box>(m, name)
    .def(nb::init<>())
    .def_static(
      "fromPoints",
      [](const std::vector<double>& lower, const std::vector<double>& upper) {
        return Box(vectorToPoint<NDIMS>(lower, "Lower"), vectorToPoint<NDIMS>(upper, "Upper"));
      },
      nb::arg("lower"),
      nb::arg("upper"))
    .def("getMin",
         [](const Box& self) { return pointToVector<NDIMS>(self.getMin()); },
         "Return the lower corner of the bounding box as a Python list.")
    .def("getMax",
         [](const Box& self) { return pointToVector<NDIMS>(self.getMax()); },
         "Return the upper corner of the bounding box as a Python list.")
    .def("contains",
         [](const Box& self, const std::vector<double>& point) {
           return self.contains(vectorToPoint<NDIMS>(point, "Point"));
         },
         nb::arg("point"),
         "Return true when the given point lies in the box.")
    .def("containsBox",
         [](const Box& self, const Box& other) { return self.contains(other); },
         nb::arg("other"),
         "Return true when the other bounding box is fully contained in this box.")
    .def("intersectsWith",
         [](const Box& self, const Box& other) { return self.intersectsWith(other); },
         nb::arg("other"),
         "Return true when the given bounding box intersects this box.")
    .def("isValid", &Box::isValid, "Return true when the bounding box is valid.")
    .def("__repr__",
         [](const Box& self) {
           std::ostringstream oss;
           oss << self;
           return oss.str();
         });
}
}  // namespace
}  // namespace primal
}  // namespace axom

NB_MODULE(pyprimal, m)
{
  m.doc() = "Python bindings for a small tutorial-facing slice of Axom Primal.";

  axom::primal::bindBoundingBox<2>(m, "BoundingBox2D");
  axom::primal::bindBoundingBox<3>(m, "BoundingBox3D");
}
