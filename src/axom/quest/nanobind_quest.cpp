// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/utilities/RAII.hpp"
#include "axom/fmt.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/quest/SamplingShaper.hpp"
#include "axom/quest/util/mesh_helpers.hpp"
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"
#include "axom/slic.hpp"

#include "mfem.hpp"

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace
{
using BoundingBox2D = axom::primal::BoundingBox<double, 2>;
using BoundingBox3D = axom::primal::BoundingBox<double, 3>;
using Point2D = axom::primal::Point<double, 2>;
using Point3D = axom::primal::Point<double, 3>;
using RuntimePolicy = axom::runtime_policy::Policy;
using SamplingMethod = axom::quest::SamplingShaper::SamplingMethod;

void ensureAxomRuntime()
{
#ifdef AXOM_USE_MPI
  static axom::utilities::raii::MPIWrapper* mpi_wrapper = nullptr;
#endif
  static axom::slic::SimpleLogger* logger = nullptr;

#ifdef AXOM_USE_MPI
  if(!mpi_wrapper)
  {
    int argc = 1;
    char program_name[] = "pyquest";
    char* argv[] = {program_name, nullptr};
    mpi_wrapper = new axom::utilities::raii::MPIWrapper(argc, argv);
  }
#endif

  if(!logger)
  {
    logger = new axom::slic::SimpleLogger();
  }

#ifdef AXOM_USE_MPI
  axom::slic::setIsRoot(mpi_wrapper->my_rank() == 0);
#else
  axom::slic::setIsRoot(true);
#endif
}

mfem::Mesh* createCartesianMesh(int dim,
                                const std::vector<double>& bb_min,
                                const std::vector<double>& bb_max,
                                const std::vector<int>& resolution,
                                int mesh_order)
{
  mfem::Mesh* mesh = nullptr;
  switch(dim)
  {
  case 2:
  {
    const BoundingBox2D bbox(Point2D {bb_min[0], bb_min[1]}, Point2D {bb_max[0], bb_max[1]});
    const axom::NumericArray<int, 2> res {resolution[0], resolution[1]};
    mesh = axom::quest::util::make_cartesian_mfem_mesh_2D(bbox, res, mesh_order);
  }
  break;
  case 3:
  {
    const BoundingBox3D bbox(Point3D {bb_min[0], bb_min[1], bb_min[2]},
                             Point3D {bb_max[0], bb_max[1], bb_max[2]});
    const axom::NumericArray<int, 3> res {resolution[0], resolution[1], resolution[2]};
    mesh = axom::quest::util::make_cartesian_mfem_mesh_3D(bbox, res, mesh_order);
  }
  break;
  default:
    throw std::invalid_argument(axom::fmt::format("'dim' must be 2 or 3; got {}.", dim));
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

std::unique_ptr<axom::sidre::MFEMSidreDataCollection> createDataCollection(
  int dim,
  const std::vector<double>& bb_min,
  const std::vector<double>& bb_max,
  const std::vector<int>& resolution,
  int mesh_order,
  const std::string& output_name)
{
  ensureAxomRuntime();

  auto* mesh = createCartesianMesh(dim, bb_min, bb_max, resolution, mesh_order);
  auto dc = std::make_unique<axom::sidre::MFEMSidreDataCollection>(output_name, nullptr, true);

#if defined(AXOM_USE_MPI)
  dc->SetMesh(MPI_COMM_WORLD, mesh);
#else
  dc->SetMesh(mesh);
#endif

  dc->SetMeshNodesName("positions");
  dc->AssociateMaterialSet("vol_frac", "material");
  return dc;
}

void importBackgroundMaterial(axom::quest::SamplingShaper& shaper,
                              axom::sidre::MFEMSidreDataCollection& collection,
                              int mesh_dim,
                              const std::string& material,
                              int volume_fraction_order)
{
  std::map<std::string, mfem::GridFunction*> initial_grid_functions;

  const auto name = axom::fmt::format("vol_frac_{}", material);
  const auto basis = mfem::BasisType::Positive;

  auto* coll = new mfem::L2_FECollection(volume_fraction_order, mesh_dim, basis);
  auto* fes = new mfem::FiniteElementSpace(collection.GetMesh(), coll);
  const int size = fes->GetVSize();

  auto* view = collection.AllocNamedBuffer(name, size);
  auto* vol_frac = new mfem::GridFunction(fes, view->getArray());
  vol_frac->MakeOwner(coll);
  (*vol_frac) = 1.0;

  collection.RegisterField(name, vol_frac);
  initial_grid_functions[material] = collection.GetField(name);
  shaper.importInitialVolumeFractions(initial_grid_functions);
}

void configureVerbosity(axom::quest::SamplingShaper& shaper, bool verbose)
{
  ensureAxomRuntime();
  axom::slic::setLoggingMsgLevel(verbose ? axom::slic::message::Debug : axom::slic::message::Info);
  shaper.setVerbosity(verbose);
}

class PySamplingShaper : public axom::quest::SamplingShaper
{
public:
  PySamplingShaper(const axom::klee::ShapeSet& shape_set,
                   int dim,
                   const std::vector<double>& bb_min,
                   const std::vector<double>& bb_max,
                   const std::vector<int>& resolution,
                   int mesh_order = 1,
                   const std::string& output_name = "shaping")
    : PySamplingShaper(std::make_shared<axom::klee::ShapeSet>(shape_set),
                       createDataCollection(dim, bb_min, bb_max, resolution, mesh_order, output_name),
                       dim,
                       RuntimePolicy::seq)
  { }

  void save() { m_collection->Save(); }

  void importBackgroundMaterial(const std::string& material, int volume_fraction_order)
  {
    ::importBackgroundMaterial(*this, *m_collection, m_mesh_dim, material, volume_fraction_order);
  }

private:
  PySamplingShaper(std::shared_ptr<axom::klee::ShapeSet> shape_set,
                   std::unique_ptr<axom::sidre::MFEMSidreDataCollection> collection,
                   int mesh_dim,
                   RuntimePolicy policy)
    : axom::quest::SamplingShaper(policy,
                                  axom::policyToDefaultAllocatorID(policy),
                                  *shape_set,
                                  collection.get())
    , m_shape_set(std::move(shape_set))
    , m_collection(std::move(collection))
    , m_mesh_dim(mesh_dim)
  { }

  std::shared_ptr<axom::klee::ShapeSet> m_shape_set;
  std::unique_ptr<axom::sidre::MFEMSidreDataCollection> m_collection;
  int m_mesh_dim;
};
}  // namespace

NB_MODULE(pyquest, m)
{
  m.doc() = "Tutorial-facing Python bindings for Axom Quest shaping.";

  nb::module_::import_("pyklee");

  nb::enum_<SamplingMethod>(m, "SamplingMethod")
    .value("InOut", SamplingMethod::InOut)
    .value("WindingNumber", SamplingMethod::WindingNumber)
    .export_values();

  nb::class_<PySamplingShaper>(m, "SamplingShaper")
    .def(nb::init<const axom::klee::ShapeSet&,
                  int,
                  const std::vector<double>&,
                  const std::vector<double>&,
                  const std::vector<int>&,
                  int,
                  const std::string&>(),
         nb::arg("shape_set"),
         nb::arg("dim"),
         nb::arg("bb_min"),
         nb::arg("bb_max"),
         nb::arg("resolution"),
         nb::arg("mesh_order") = 1,
         nb::arg("output_name") = "shaping")
    .def("setVerbosity",
         [](PySamplingShaper& self, bool verbose) { configureVerbosity(self, verbose); })
    .def("setSamplesPerKnotSpan", &PySamplingShaper::setSamplesPerKnotSpan, nb::arg("samples"))
    .def("setSamplingMethod", &PySamplingShaper::setSamplingMethod, nb::arg("sampling_method"))
    .def("setSamplingResolution",
         nb::overload_cast<int>(&PySamplingShaper::setSamplingResolution),
         nb::arg("resolution"))
    .def("setVolumeFractionOrder", &PySamplingShaper::setVolumeFractionOrder, nb::arg("order"))
    .def("loadShape", &PySamplingShaper::loadShape, nb::arg("shape"))
    .def("prepareShapeQuery",
         &PySamplingShaper::prepareShapeQuery,
         nb::arg("shape_dimension"),
         nb::arg("shape"))
    .def("runShapeQuery", &PySamplingShaper::runShapeQuery, nb::arg("shape"))
    .def("applyReplacementRules", &PySamplingShaper::applyReplacementRules, nb::arg("shape"))
    .def("finalizeShapeQuery", &PySamplingShaper::finalizeShapeQuery)
    .def("adjustVolumeFractions", &PySamplingShaper::adjustVolumeFractions)
    .def("importBackgroundMaterial",
         &PySamplingShaper::importBackgroundMaterial,
         nb::arg("material"),
         nb::arg("volume_fraction_order"))
    .def("save", &PySamplingShaper::save);
}
