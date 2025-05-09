// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \brief Unit tests for quest's SamplingShaper class and associated replacement rules.
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/klee.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"
#include "axom/sidre.hpp"
#include "axom/slic.hpp"
#include "axom/quest/SamplingShaper.hpp"
#include "axom/quest/util/mesh_helpers.hpp"

#ifndef AXOM_USE_MFEM
  #error "Quest's SamplingShaper tests on mfem meshes require mfem library."
#else
  #include "mfem.hpp"
#endif

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace sidre = axom::sidre;
namespace slic = axom::slic;

namespace
{
const std::string unit_circle_contour =
  "piece = circle(origin=(0cm, 0cm), radius=1cm, start=0deg, end=360deg)";

const std::string unit_semicircle_contour = R"(
  piece = circle(origin=(0cm, 0cm), radius=1cm, start=0deg, end=180deg)
  piece = line(end=(0cm, 1cm)))";

// Set the following to true for verbose output and for saving vis files
constexpr bool very_verbose_output = false;

/// RAII utility class to write a file at construction time and remove it
/// once the instance is out of scope
class ScopedTemporaryFile
{
public:
  ScopedTemporaryFile(const std::string& filename, const std::string& contents)
    : m_filename(filename)
  {
    std::ofstream ofs(m_filename.c_str(), std::ios::out);
    ofs << contents;
  }

  ~ScopedTemporaryFile() { EXPECT_EQ(axom::utilities::filesystem::removeFile(m_filename), 0); }

  const std::string& getFileName() const { return m_filename; }

  std::string getFileContents() const
  {
    std::stringstream buffer;

    std::ifstream ifs(m_filename.c_str(), std::ios::in);
    if(ifs.is_open())
    {
      buffer << ifs.rdbuf();
    }

    return buffer.str();
  }

private:
  const std::string m_filename;
};
}  // namespace

/// Test fixture for SamplingShaper tests on MFEM meshes
class SamplingShaperTest : public ::testing::Test
{
public:
  SamplingShaperTest() : m_dc("test", nullptr, true) { }

  virtual ~SamplingShaperTest() { }

  void SetUp() override { }

  sidre::MFEMSidreDataCollection& getDC() { return m_dc; }
  mfem::Mesh& getMesh() { return *m_dc.GetMesh(); }

  /// parse and validate the Klee shapefile; fail the test if invalid
  void validateShapeFile(const std::string& shapefile)
  {
    axom::klee::ShapeSet shapeSet;

    try
    {
      shapeSet = axom::klee::readShapeSet(shapefile);
    }
    catch(axom::klee::KleeError& error)
    {
      std::vector<std::string> errs;
      for(auto verificationError : error.getErrors())
      {
        errs.push_back(axom::fmt::format(" - '{}': {}",
                                         static_cast<std::string>(verificationError.path),
                                         verificationError.message));
      }

      if(!errs.empty())
      {
        SLIC_WARNING(
          axom::fmt::format("Error during parsing klee input file '{}'. "
                            "Found the following errors:\n{}",
                            shapefile,
                            axom::fmt::join(errs, "\n")));
        FAIL();
      }
    }
  }

  /// Initializes the Shaper instance over a shapefile and optionally sets up initial "preshaped" volume fractions
  void initializeShaping(const std::string& shapefile,
                         const std::map<std::string, mfem::GridFunction*>& init_vf_map = {})
  {
    SLIC_INFO_IF(very_verbose_output, axom::fmt::format("Reading shape set from {}", shapefile));
    m_shapeSet = std::make_unique<klee::ShapeSet>(klee::readShapeSet(shapefile));

    SLIC_INFO_IF(very_verbose_output, axom::fmt::format("Shaping materials..."));
    m_shaper = std::make_unique<quest::SamplingShaper>(*m_shapeSet, &m_dc);
    m_shaper->setVerbosity(very_verbose_output);

    if(!init_vf_map.empty())
    {
      m_shaper->importInitialVolumeFractions(init_vf_map);
    }

    if(very_verbose_output)
    {
      m_shaper->printRegisteredFieldNames("*** After importing volume fractions");
    }
  }

  /// Runs the shaping query over a shapefile; must be called after initializeShaping()
  void runShaping()
  {
    EXPECT_NE(nullptr, m_shaper) << "Shaper needs to be initialized via initializeShaping()";

    // Define lambda to override default dimensions, when necessary
    auto getShapeDim = [defaultDim = m_shapeSet->getDimensions()](const auto& shape) {
      static std::map<std::string, klee::Dimensions> format_dim = {{"c2c", klee::Dimensions::Two},
                                                                   {"stl", klee::Dimensions::Three}};

      const auto& format_str = shape.getGeometry().getFormat();
      return format_dim.find(format_str) != format_dim.end() ? format_dim[format_str] : defaultDim;
    };

    for(const auto& shape : m_shapeSet->getShapes())
    {
      SLIC_INFO_IF(
        very_verbose_output,
        axom::fmt::format("\tshape {} -> material {}", shape.getName(), shape.getMaterial()));

      const auto shapeDim = getShapeDim(shape);

      m_shaper->loadShape(shape);
      m_shaper->prepareShapeQuery(shapeDim, shape);
      m_shaper->runShapeQuery(shape);
      m_shaper->applyReplacementRules(shape);
      m_shaper->finalizeShapeQuery();
    }

    m_shaper->adjustVolumeFractions();

    if(very_verbose_output)
    {
      m_shaper->printRegisteredFieldNames("*** After shaping volume fractions");
    }
  }

  // Computes the total volume of the associated volume fraction grid function
  double gridFunctionVolume(const std::string& name)
  {
    mfem::GridFunction* gf = m_dc.GetField(name);

    mfem::ConstantCoefficient one(1.0);
    mfem::LinearForm vol_form(gf->FESpace());
    vol_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    vol_form.Assemble();

    return *gf * vol_form;
  }

  /// Registers and allocates a volume fraction grid function within the datastore
  mfem::GridFunction* registerVolFracGridFunction(const std::string& name, int vfOrder = 2)
  {
    SLIC_ASSERT(!m_dc.HasField(name));

    auto& mesh = getMesh();
    const int dim = mesh.Dimension();

    // create grid function
    auto* coll = new mfem::L2_FECollection(vfOrder, dim, mfem::BasisType::Positive);
    auto* fes = new mfem::FiniteElementSpace(&mesh, coll);
    auto* vf = new mfem::GridFunction(fes);
    vf->MakeOwner(coll);

    // allocate grid function via sidre
    const int sz = fes->GetVSize();
    mfem::Vector v(m_dc.AllocNamedBuffer(name, sz)->getData(), sz);
    vf->MakeRef(fes, v, 0);

    // register grid function w/ data collection
    m_dc.RegisterField(name, vf);

    return vf;
  }

  /** 
   * \brief Initializes the values of the DOFs of a volume fraction grid function
   * using a provided lambda w/ parameters for the cell index, DOF position and cell attribute
   * The signature of DOFInitializer is [](int idx, Point<double,DIM>>& pt, int attribute) -> double
   */
  template <int DIM, typename DOFInitializer>
  void initializeVolFracGridFunction(mfem::GridFunction* vf, DOFInitializer&& dof_initializer)
  {
    auto& mesh = this->getMesh();
    const int dim = mesh.Dimension();
    const int NE = mesh.GetNE();

    // Assume all elements have the same integration rule
    const auto* fes = vf->FESpace();
    auto* fe = fes->GetFE(0);
    auto& ir = fe->GetNodes();
    const int nq = ir.GetNPoints();

    // Get positions of DOFs
    mfem::DenseTensor pos_coef(dim, nq, NE);
    {
      const auto* geomFactors = mesh.GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);

      // Rearrange positions
      for(int i = 0; i < NE; ++i)
      {
        for(int j = 0; j < dim; ++j)
        {
          for(int k = 0; k < nq; ++k)
          {
            pos_coef(j, k, i) = geomFactors->X((i * nq * dim) + (j * nq) + k);
          }
        }
      }
    }

    // Initialize volume fraction DOFs using passed in lambda based on cell index, DOF position and attribute
    mfem::Vector res(nq);
    mfem::Array<int> dofs;
    for(int idx = 0; idx < NE; ++idx)
    {
      const int attr = mesh.GetAttribute(idx);

      mfem::DenseMatrix& m = pos_coef(idx);
      for(int p = 0; p < nq; ++p)
      {
        const primal::Point<double, DIM> pt(m.GetColumn(p), dim);
        res(p) = dof_initializer(idx, pt, attr);
      }

      fes->GetElementDofs(idx, dofs);
      vf->SetSubVector(dofs, res);
    }
  }

  /// Helper to check integrated volume of a volume fraction grid fucntion
  void checkExpectedVolumeFractions(const std::string& material_name,
                                    double expected_volume,
                                    double EPS = 1e-2)
  {
    auto vf_name = axom::fmt::format("vol_frac_{}", material_name);

    EXPECT_TRUE(m_dc.HasField(vf_name))
      << axom::fmt::format("Did not have expected volume fraction '{:.4}' for material '{}'",
                           material_name,
                           vf_name);

    const double actual_volume = this->gridFunctionVolume(vf_name);
    SLIC_INFO(axom::fmt::format("Shaped volume fraction of '{}' is {:.4}  (expected: {:.4})",
                                material_name,
                                actual_volume,
                                expected_volume));

    EXPECT_NEAR(expected_volume, actual_volume, EPS);
  }

protected:
  sidre::MFEMSidreDataCollection m_dc;
  std::unique_ptr<klee::ShapeSet> m_shapeSet;
  std::unique_ptr<quest::SamplingShaper> m_shaper;
};

/// Test fixture for SamplingShaper tests on 2D MFEM meshes
class SamplingShaperTest2D : public SamplingShaperTest
{
public:
  using Point2D = primal::Point<double, 2>;
  using BBox2D = primal::BoundingBox<double, 2>;

public:
  virtual ~SamplingShaperTest2D() { }

  void SetUp() override
  {
    const int polynomialOrder = 2;
    const BBox2D bbox({-2, -2}, {2, 2});
    const axom::NumericArray<int, 2> celldims {64, 64};

    // memory for mesh will be managed by data collection
    auto* mesh = quest::util::make_cartesian_mfem_mesh_2D(bbox, celldims, polynomialOrder);

    // Set element attributes based on quadrant where centroid is located
    // These will be used later in some cases when setting volume fractions
    mfem::Array<int> v;
    const int NE = mesh->GetNE();
    for(int i = 0; i < NE; ++i)
    {
      mesh->GetElementVertices(i, v);
      BBox2D elem_bbox;
      for(int j = 0; j < v.Size(); ++j)
      {
        elem_bbox.addPoint(Point2D(mesh->GetVertex(v[j]), 2));
      }

      const auto centroid = elem_bbox.getCentroid();
      if(centroid[0] >= 0 && centroid[1] >= 0)
      {
        mesh->SetAttribute(i, 1);
      }
      else if(centroid[0] >= 0 && centroid[1] < 0)
      {
        mesh->SetAttribute(i, 2);
      }
      else if(centroid[0] < 0 && centroid[1] >= 0)
      {
        mesh->SetAttribute(i, 3);
      }
      else
      {
        mesh->SetAttribute(i, 4);
      }
    }

    m_dc.SetOwnData(true);
    m_dc.SetMeshNodesName("positions");
    m_dc.SetMesh(mesh);

#ifdef AXOM_USE_MPI
    m_dc.SetComm(MPI_COMM_WORLD);
#endif
  }

  BBox2D meshBoundingBox()
  {
    mfem::Vector bbmin, bbmax;
    getMesh().GetBoundingBox(bbmin, bbmax);

    return BBox2D(Point2D(bbmin.GetData()), Point2D(bbmax.GetData()));
  }
};

/// Test fixture for SamplingShaper tests on 3D MFEM meshes
class SamplingShaperTest3D : public SamplingShaperTest
{
public:
  using Point3D = primal::Point<double, 3>;
  using BBox3D = primal::BoundingBox<double, 3>;

public:
  virtual ~SamplingShaperTest3D() { }

  void SetUp() override
  {
    const int polynomialOrder = 2;
    const BBox3D bbox({-2, -2, -2}, {2, 2, 2});
    const axom::NumericArray<int, 3> celldims {8, 8, 8};

    // memory for mesh will be managed by data collection
    auto* mesh = quest::util::make_cartesian_mfem_mesh_3D(bbox, celldims, polynomialOrder);

    // Set element attributes based on octant where centroid is located
    // These will be used later in some cases when setting volume fractions
    mfem::Array<int> v;
    const int NE = mesh->GetNE();
    for(int i = 0; i < NE; ++i)
    {
      mesh->GetElementVertices(i, v);
      BBox3D elem_bbox;
      for(int j = 0; j < v.Size(); ++j)
      {
        elem_bbox.addPoint(Point3D(mesh->GetVertex(v[j]), 3));
      }
      const auto centroid = elem_bbox.getCentroid();

      int attr = 0;
      attr |= (centroid[0] < 0) ? 1 << 0 : 0;
      attr |= (centroid[1] < 0) ? 1 << 1 : 0;
      attr |= (centroid[2] < 0) ? 1 << 2 : 0;
      mesh->SetAttribute(i, attr);
    }

    m_dc.SetOwnData(true);
    m_dc.SetMeshNodesName("positions");
    m_dc.SetMesh(mesh);

#ifdef AXOM_USE_MPI
    m_dc.SetComm(MPI_COMM_WORLD);
#endif
  }

  BBox3D meshBoundingBox()
  {
    mfem::Vector bbmin, bbmax;
    getMesh().GetBoundingBox(bbmin, bbmax);

    return BBox3D(Point3D(bbmin.GetData()), Point3D(bbmax.GetData()));
  }
};

//-----------------------------------------------------------------------------

TEST(ScopedTemporaryFile, basic)
{
  using axom::utilities::filesystem::pathExists;

  const std::string filename = "previously_nonexistent.file";
  const std::string contents = "file contents!";

  // File does not exist before entering scope
  EXPECT_FALSE(pathExists(filename));

  // File is created and exists within the scope
  {
    ScopedTemporaryFile test_file(filename, contents);

    EXPECT_EQ(test_file.getFileName(), filename);
    ASSERT_TRUE(pathExists(test_file.getFileName()));

    // check contents
    EXPECT_EQ(contents, test_file.getFileContents());
  }

  // File no longer exists outside the scope
  EXPECT_FALSE(pathExists(filename));
}

//-----------------------------------------------------------------------------

TEST_F(SamplingShaperTest2D, check_mesh)
{
  auto& mesh = this->getMesh();

  const int NE = mesh.GetNE();
  SLIC_INFO(axom::fmt::format("The mesh has {} elements", NE));
  EXPECT_GT(NE, 0);

  const int NV = mesh.GetNV();
  SLIC_INFO(axom::fmt::format("The mesh has {} vertices", NV));
  EXPECT_GT(NV, 0);

  const auto bbox = this->meshBoundingBox();
  SLIC_INFO(axom::fmt::format("The mesh bounding box is: {}", bbox));
  EXPECT_TRUE(bbox.isValid());
}

//-----------------------------------------------------------------------------

TEST_F(SamplingShaperTest2D, basic_circle)
{
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: circle_shape
  material: {}
  geometry:
    format: c2c
    path: {}
)";

  const std::string circle_material = "circleMat";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname), unit_circle_contour);

  ScopedTemporaryFile shape_file(
    axom::fmt::format("{}.yaml", testname),
    axom::fmt::format(shape_template, circle_material, contour_file.getFileName()));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName());
  this->runShaping();

  // check that the result has a volume fraction field associated with the circle material
  constexpr double expected_volume = M_PI;
  this->checkExpectedVolumeFractions(circle_material, expected_volume);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest2D, basic_circle_projector)
{
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: circle_shape
  material: {}
  geometry:
    format: c2c
    path: {}
)";

  const std::string circle_material = "circleMat";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname), unit_circle_contour);

  ScopedTemporaryFile shape_file(
    axom::fmt::format("{}.yaml", testname),
    axom::fmt::format(shape_template, circle_material, contour_file.getFileName()));

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName());

  // check that we can set several projectors in 2D and 3D
  // uses simplest projectors, e.g. identity in 2D and 3D
  this->m_shaper->setPointProjector([](const Point3D& pt) {
    return Point3D {pt[0], pt[1], pt[2]};
  });
  this->m_shaper->setPointProjector([](const Point2D& pt) { return Point2D {pt[0], pt[1]}; });
  this->m_shaper->setPointProjector([](const Point3D& pt) { return Point2D {pt[0], pt[1]}; });
  this->m_shaper->setPointProjector([](const Point2D& pt) { return Point3D {pt[0], pt[1], 0}; });

  this->runShaping();

  // check that the result has a volume fraction field associated with the circle material
  constexpr double expected_volume = M_PI;
  this->checkExpectedVolumeFractions(circle_material, expected_volume);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest2D, circle_projector_anisotropic)
{
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: circle_shape
  material: {}
  geometry:
    format: c2c
    path: {}
)";

  const std::string circle_material = "circleMat";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname), unit_circle_contour);

  ScopedTemporaryFile shape_file(
    axom::fmt::format("{}.yaml", testname),
    axom::fmt::format(shape_template, circle_material, contour_file.getFileName()));

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName());

  // check that we can set several projectors in 2D and 3D
  // creating an ellipse by scaling input x and y by scale_a and scale_b
  constexpr double scale_a = 3. / 2.;
  constexpr double scale_b = 3. / 4.;
  this->m_shaper->setPointProjector([](const Point2D& pt) {
    return Point2D {pt[0] / scale_a, pt[1] / scale_b};
  });
  // check that we can register another projector that's not used
  this->m_shaper->setPointProjector([](const Point3D&) { return Point3D {0., 0.}; });

  this->runShaping();

  // check that the result has a volume fraction field associated with the circle material
  constexpr double expected_volume = M_PI * scale_a * scale_b;
  this->checkExpectedVolumeFractions(circle_material, expected_volume);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest2D, disk_via_replacement)
{
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2
units: cm

shapes:
- name: circle_outer
  material: outer
  geometry:
    format: c2c
    path: {0}
    units: cm
- name: void_inner
  material: inner
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: .5
)";

  const std::string outer_material = "outer";
  const std::string inner_material = "inner";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname), unit_circle_contour);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template, contour_file.getFileName()));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());

  this->initializeShaping(shape_file.getFileName());
  this->runShaping();

  // check that the result has a volume fraction field associated with the circle material
  constexpr double expected_inner_area = .5 * .5 * M_PI;
  constexpr double expected_outer_area = M_PI - expected_inner_area;
  this->checkExpectedVolumeFractions(outer_material, expected_outer_area);
  this->checkExpectedVolumeFractions(inner_material, expected_inner_area);

  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest2D, disk_via_replacement_with_background)
{
  using Point2D = typename SamplingShaperTest2D::Point2D;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2
units: cm

shapes:
- name: background
  material: {1}
  geometry:
    format: none
- name: circle_outer
  material: {2}
  geometry:
    format: c2c
    path: {0}
    units: cm
- name: void_inner
  material: {3}
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: .5
)";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname), unit_circle_contour);

  // Set background material to 'void' (which is not present elsewhere)
  {
    ScopedTemporaryFile shape_file(
      axom::fmt::format("{}.yaml", testname),
      axom::fmt::format(shape_template, contour_file.getFileName(), "void", "disk", "hole"));

    // Create an initial background material set to 1 everywhere
    std::map<std::string, mfem::GridFunction*> initialGridFunctions;
    {
      auto* vf = this->registerVolFracGridFunction("init_vf_bg");
      this->initializeVolFracGridFunction<2>(vf,
                                             [](int, const Point2D&, int) -> double { return 1.; });
      initialGridFunctions["void"] = vf;
    }

    this->validateShapeFile(shape_file.getFileName());
    this->initializeShaping(shape_file.getFileName(), initialGridFunctions);
    this->runShaping();

    // check that the result has a volume fraction field associated with the circle material
    constexpr double expected_hole_area = .5 * .5 * M_PI;
    constexpr double expected_disk_area = M_PI - expected_hole_area;
    const auto range = this->meshBoundingBox().range();
    const double expected_bg_area = range[0] * range[1] - M_PI;

    this->checkExpectedVolumeFractions("disk", expected_disk_area);
    this->checkExpectedVolumeFractions("hole", expected_hole_area);
    this->checkExpectedVolumeFractions("void", expected_bg_area);
  }

  // clean up data collection
  for(const auto& name : {"vol_frac_void", "vol_frac_hole", "vol_frac_disk", "init_vf_bg"})
  {
    this->getDC().DeregisterField(name);
  }

  // Set background and inner hole materials to 'void'
  {
    ScopedTemporaryFile shape_file(
      axom::fmt::format("{}.yaml", testname),
      axom::fmt::format(shape_template, contour_file.getFileName(), "void", "disk", "void"));

    // Create an initial background material set to 1 everywhere
    std::map<std::string, mfem::GridFunction*> initialGridFunctions;
    {
      auto* vf = this->registerVolFracGridFunction("init_vf_bg");
      this->initializeVolFracGridFunction<2>(vf,
                                             [](int, const Point2D&, int) -> double { return 1.; });
      initialGridFunctions["void"] = vf;
    }

    this->validateShapeFile(shape_file.getFileName());
    this->initializeShaping(shape_file.getFileName(), initialGridFunctions);
    this->runShaping();

    // check that the result has a volume fraction field associated with the circle material
    constexpr double expected_disk_area = M_PI - .5 * .5 * M_PI;
    const auto range = this->meshBoundingBox().range();
    const double expected_void_area = range[0] * range[1] - expected_disk_area;

    this->checkExpectedVolumeFractions("disk", expected_disk_area);
    this->checkExpectedVolumeFractions("void", expected_void_area);
  }

  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest2D, preshaped_materials)
{
  using Point2D = typename SamplingShaperTest2D::Point2D;

  const std::string& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2
units: cm

shapes:
- name: background
  material: {0}
  geometry:
    format: none
# Left replaces background void
- name: left_side
  material: {1}
  geometry:
    format: none
# Odd cells replace background void, but not left
- name: odd_cells
  material: {2}
  geometry:
    format: none
  does_not_replace: [{1}]
)";

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template, "void", "left", "odds"));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  // Create an initial background material set to 1 everywhere
  std::map<std::string, mfem::GridFunction*> initialGridFunctions;
  {
    auto* vf = this->registerVolFracGridFunction("init_vf_bg");
    this->initializeVolFracGridFunction<2>(vf, [](int, const Point2D&, int) -> double { return 1.; });
    initialGridFunctions["void"] = vf;

    // Note: element attributes were set earlier based on quadrant of cell's centroid (1, 2, 3 and 4)
    vf = this->registerVolFracGridFunction("init_vf_left");
    this->initializeVolFracGridFunction<2>(vf, [](int, const Point2D&, int attr) -> double {
      return (attr == 3 || attr == 4) ? 1. : 0;
    });
    initialGridFunctions["left"] = vf;

    vf = this->registerVolFracGridFunction("init_vf_odds");
    this->initializeVolFracGridFunction<2>(vf, [](int idx, const Point2D&, int) -> double {
      return idx % 2 == 1 ? 1. : 0;
    });
    initialGridFunctions["odds"] = vf;
  }

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName(), initialGridFunctions);
  this->runShaping();

  // check that the result has a volume fraction field associated with the circle material
  const auto range = this->meshBoundingBox().range();

  // Left covers half the mesh and is not replaced
  const double expected_left_area = range[0] * range[1] / 2.;
  // Odds should cover half of the right side of the mesh
  const double expected_odds_area = range[0] * range[1] / 4.;
  // The rest should be void
  const double expected_void_area = range[0] * range[1] - expected_left_area - expected_odds_area;

  this->checkExpectedVolumeFractions("left", expected_left_area);
  this->checkExpectedVolumeFractions("odds", expected_odds_area);
  this->checkExpectedVolumeFractions("void", expected_void_area);

  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest2D, disk_with_multiple_preshaped_materials)
{
  using Point2D = typename SamplingShaperTest2D::Point2D;

  const std::string& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2
units: cm

shapes:
- name: background
  material: void
  geometry:
    format: none
- name: circle_outer
  material: disk
  geometry:
    format: c2c
    path: {0}
    units: cm
- name: circle_inner
  material: hole
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: .5
  replaces: [void, disk]
- name: left_side
  material: left
  geometry:
    format: none
  does_not_replace: [disk]
- name: odd_cells
  material: odds
  geometry:
    format: none
  replaces: [void, hole]
)";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname), unit_circle_contour);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template, contour_file.getFileName()));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  // Create an initial background material set to 1 everywhere
  std::map<std::string, mfem::GridFunction*> initialGridFunctions;
  {
    // initial background void material is set everywhere
    auto* vf = this->registerVolFracGridFunction("init_vf_bg");
    this->initializeVolFracGridFunction<2>(vf, [](int, const Point2D&, int) -> double { return 1.; });
    initialGridFunctions["void"] = vf;

    // initial left material is set based on mesh attributes
    // Note: element attributes were set earlier based on quadrant of cell's centroid (1, 2, 3 and 4)
    vf = this->registerVolFracGridFunction("init_vf_left");
    this->initializeVolFracGridFunction<2>(vf, [](int, const Point2D&, int attr) -> double {
      return (attr == 3 || attr == 4) ? 1. : 0;
    });
    initialGridFunctions["left"] = vf;

    // initial "odds" material is based on the parity of the element indices
    vf = this->registerVolFracGridFunction("init_vf_odds");
    this->initializeVolFracGridFunction<2>(vf, [](int idx, const Point2D&, int) -> double {
      return idx % 2 == 1 ? 1. : 0;
    });
    initialGridFunctions["odds"] = vf;
  }

  // For this example, we keep the full disk between radii 1 and .5
  // The 'left' material is set for all cells to the left of the y-axis
  // but does not replace the disk material
  // The interior hole is within radius .5, but is replaced by left and odds
  // The odds material is all cells w/ odd index, but not covering 'left' or 'disk'
  // The end result for the void background is everything that's left
  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName(), initialGridFunctions);
  this->runShaping();

  // check that the result has the correct volume fractions
  const auto range = this->meshBoundingBox().range();
  const auto total_area = range[0] * range[1];
  const auto left_orig = total_area / 2;
  constexpr auto hole_orig = .5 * .5 * M_PI;

  constexpr double expected_disk_area = M_PI - hole_orig;
  const double expected_left_area = left_orig - expected_disk_area / 2;
  const double expected_hole_area = 0.19683;  // from program output
  const double expected_odds_area = 3.41180;  // from program output
  const double expected_void_area = 3.21497;  // from program output

  this->checkExpectedVolumeFractions("left", expected_left_area);
  this->checkExpectedVolumeFractions("disk", expected_disk_area);
  this->checkExpectedVolumeFractions("odds", expected_odds_area);
  this->checkExpectedVolumeFractions("hole", expected_hole_area);
  this->checkExpectedVolumeFractions("void", expected_void_area);

  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest2D, check_underscores)
{
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  constexpr double radius = 1.5;

  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: {2}
  material: {3}
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: {1}
- name: {4}
  material: {5}
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: {1}
- name: {6}
  material: {7}
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: {1}
)";

  const std::string shape_name {"shape"};
  const std::string mat_name {"mat"};

  const std::string underscored_shape_name {"underscored_shape"};
  const std::string underscored_mat_name {"underscored_mat"};

  const std::string double_underscored_shape_name {"double_underscored_shape"};
  const std::string double_underscored_mat_name {"double_underscored_mat"};

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname),
                                   unit_semicircle_contour);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template,
                                                   contour_file.getFileName(),
                                                   radius,
                                                   shape_name,
                                                   mat_name,
                                                   underscored_shape_name,
                                                   underscored_mat_name,
                                                   double_underscored_shape_name,
                                                   double_underscored_mat_name));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName());

  this->runShaping();

  // Collect and print registered fields
  std::vector<std::string> regFields;
  for(const auto& pr : this->getDC().GetFieldMap())
  {
    regFields.push_back(pr.first);
  }
  SLIC_INFO(axom::fmt::format("Registered fields: {}", axom::fmt::join(regFields, ", ")));

  // check that output materials are present
  EXPECT_TRUE(this->getDC().HasField(axom::fmt::format("vol_frac_{}", mat_name)));
  EXPECT_TRUE(this->getDC().HasField(axom::fmt::format("vol_frac_{}", underscored_mat_name)));
  EXPECT_TRUE(this->getDC().HasField(axom::fmt::format("vol_frac_{}", double_underscored_mat_name)));

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest2D, contour_and_stl_2D)
{
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  constexpr double radius = 1.5;

  const std::string shape_template = R"(
dimensions: 2

shapes:
# preshape a background material; dimension should be default
- name: background
  material: {3}
  geometry:
    format: none
# shape in a revolved sphere given as a c2c contour
- name: circle_shape
  material: {4}
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: {2}
# shape in a sphere given as an stl surface mesh
- name: sphere_shape
  material: {5}
  geometry:
    format: stl
    path: {1}
    units: cm
)";

  const std::string background_material = "luminiferous_ether";
  const std::string circle_material = "steel";
  const std::string sphere_material = "vaccum";
  const std::string sphere_path = axom::fmt::format("{}/quest/unit_sphere.stl", AXOM_DATA_DIR);

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname), unit_circle_contour);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template,
                                                   contour_file.getFileName(),
                                                   sphere_path,
                                                   radius,
                                                   background_material,
                                                   circle_material,
                                                   sphere_material));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  // Create an initial background material set to 1 everywhere
  std::map<std::string, mfem::GridFunction*> initialGridFunctions;
  {
    auto* vf = this->registerVolFracGridFunction("init_vf_bg");
    this->initializeVolFracGridFunction<2>(vf, [](int, const Point2D&, int) -> double { return 1.; });
    initialGridFunctions[background_material] = vf;
  }

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName(), initialGridFunctions);

  // set projector from 2D mesh points to 3D query points within STL
  this->m_shaper->setPointProjector([](Point2D pt) { return Point3D {pt[0], pt[1], 0.}; });

  this->m_shaper->setQuadratureOrder(8);

  this->runShaping();

  // Check that the result has a volume fraction field associated with circle and sphere materials
  constexpr double exp_volume_contour = M_PI * radius * radius;
  constexpr double exp_volume_sphere = M_PI * 1. * 1.;
  this->checkExpectedVolumeFractions(circle_material, exp_volume_contour - exp_volume_sphere, 3e-2);
  this->checkExpectedVolumeFractions(sphere_material, exp_volume_sphere, 3e-2);

  for(const auto& vf_name : {background_material, circle_material, sphere_material})
  {
    EXPECT_TRUE(this->getDC().HasField(axom::fmt::format("vol_frac_{}", vf_name)));
  }

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

//-----------------------------------------------------------------------------

TEST_F(SamplingShaperTest3D, basic_tet)
{
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 3

shapes:
- name: tet_shape
  material: {}
  geometry:
    format: stl
    path: {}
)";

  const std::string tet_material = "steel";
  const std::string tet_path = axom::fmt::format("{}/quest/tetrahedron.stl", AXOM_DATA_DIR);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template, tet_material, tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Bounding box of 3D input mesh: \n" << this->meshBoundingBox());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName());
  this->runShaping();

  // Check that the result has a volume fraction field associated with the tetrahedron material
  // The tet lives in cube of edge length 2 (and volume 8) and is defined by opposite corners.
  // It occupies 1/3 of the cube's volume
  constexpr double expected_volume = 8. / 3.;
  this->checkExpectedVolumeFractions(tet_material, expected_volume);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest3D, tet_preshaped)
{
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 3

shapes:
- name: octant_0
  material: octant0
  geometry:
    format: none
- name: octant_1
  material: octant1
  geometry:
    format: none
- name: octant_2
  material: octant2
  geometry:
    format: none
- name: octant_3
  material: octant3
  geometry:
    format: none
- name: octant_4
  material: octant4
  geometry:
    format: none
- name: octant_5
  material: octant5
  geometry:
    format: none
- name: octant_6
  material: octant6
  geometry:
    format: none
- name: octant_7
  material: octant7
  geometry:
    format: none
- name: tet_shape
  material: {}
  geometry:
    format: stl
    path: {}
)";

  const std::string tet_material = "steel";
  const std::string tet_path = axom::fmt::format("{}/quest/tetrahedron.stl", AXOM_DATA_DIR);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template, tet_material, tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());

  // Create initial background materials based on octant attributes
  // Octants were offset by one since mfem doesn't allow setting attribute to zero
  std::map<std::string, mfem::GridFunction*> initialGridFunctions;
  {
    for(int attr_i = 0; attr_i < 8; ++attr_i)
    {
      auto* vf = this->registerVolFracGridFunction(axom::fmt::format("init_vf_octant_{}", attr_i));
      this->initializeVolFracGridFunction<3>(vf, [attr_i](int, const Point3D&, int attr) -> double {
        return attr == attr_i ? 1 : 0;
      });
      initialGridFunctions[axom::fmt::format("octant{}", attr_i)] = vf;
    }
  }

  this->initializeShaping(shape_file.getFileName(), initialGridFunctions);
  this->runShaping();

  // Check that the result has a volume fraction field associated with the tetrahedron material
  // The tet lives in cube of edge length 2 (and volume 8) and is defined by opposite corners.
  // It occupies 1/3 of the cube's volume
  constexpr double tet_volume = 8. / 3.;
  this->checkExpectedVolumeFractions(tet_material, tet_volume);

  // The background mesh is a cube of edge length 4 centered around the origin
  // Each octant's volume is 8 and its vf gets overlaid by a piece of the tet
  constexpr double missing_half = 8. - 1. / 2.;
  constexpr double missing_sixth = 8. - 1. / 6.;
  this->checkExpectedVolumeFractions("octant0", missing_half);
  this->checkExpectedVolumeFractions("octant1", missing_sixth);
  this->checkExpectedVolumeFractions("octant2", missing_sixth);
  this->checkExpectedVolumeFractions("octant3", missing_half);
  this->checkExpectedVolumeFractions("octant4", missing_sixth);
  this->checkExpectedVolumeFractions("octant5", missing_half);
  this->checkExpectedVolumeFractions("octant6", missing_half);
  this->checkExpectedVolumeFractions("octant7", missing_sixth);

  constexpr double total_volume = 4 * 4 * 4;
  EXPECT_EQ(total_volume, tet_volume + 4 * missing_sixth + 4 * missing_half);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest3D, tet_preshaped_with_replacements)
{
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  // Use somewhat complex rules: tet's material will replace octants 0-3, but not 4-7
  const std::string shape_template = R"(
dimensions: 3

shapes:
- name: octant_0
  material: octant0
  geometry:
    format: none
- name: octant_1
  material: octant1
  geometry:
    format: none

- name: octant_4
  material: octant4
  geometry:
    format: none
- name: octant_5
  material: octant5
  geometry:
    format: none

- name: tet_shape
  material: steel
  geometry:
    format: stl
    path: {}
  replaces: [octant0,octant1]

- name: octant_6
  material: octant6
  geometry:
    format: none
- name: octant_7
  material: octant7
  geometry:
    format: none

- name: octant_2
  material: octant2
  geometry:
    format: none
  does_not_replace: [steel]
- name: octant_3
  material: octant3
  geometry:
    format: none
  does_not_replace: [steel]
)";

  const std::string tet_path = axom::fmt::format("{}/quest/tetrahedron.stl", AXOM_DATA_DIR);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template, tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());

  // Create initial background materials based on octant attributes
  std::map<std::string, mfem::GridFunction*> initialGridFunctions;
  {
    for(int attr_i = 0; attr_i < 8; ++attr_i)
    {
      auto* vf = this->registerVolFracGridFunction(axom::fmt::format("init_vf_octant_{}", attr_i));
      this->initializeVolFracGridFunction<3>(vf, [attr_i](int, const Point3D&, int attr) -> double {
        return attr == attr_i ? 1 : 0;
      });
      initialGridFunctions[axom::fmt::format("octant{}", attr_i)] = vf;
    }
  }

  this->initializeShaping(shape_file.getFileName(), initialGridFunctions);
  this->runShaping();

  // Check that the result has a volume fraction field associated with the tetrahedron material
  // The tet has volume 8/3, but only half of it is replaced
  constexpr double tet_volume = 8. / 3.;
  this->checkExpectedVolumeFractions("steel", tet_volume / 2.);

  // The background mesh is a cube of edge length 4 centered around the origin
  // octants 0-3 are replaced by the tet, but 4-7 are not
  constexpr double missing_half = 8. - 1. / 2.;
  constexpr double missing_sixth = 8. - 1. / 6.;
  this->checkExpectedVolumeFractions("octant0", missing_half);
  this->checkExpectedVolumeFractions("octant1", missing_sixth);
  this->checkExpectedVolumeFractions("octant2", missing_sixth);
  this->checkExpectedVolumeFractions("octant3", missing_half);
  this->checkExpectedVolumeFractions("octant4", 8.);
  this->checkExpectedVolumeFractions("octant5", 8.);
  this->checkExpectedVolumeFractions("octant6", 8.);
  this->checkExpectedVolumeFractions("octant7", 8.);

  constexpr double total_volume = 4 * 4 * 4;
  EXPECT_EQ(total_volume, tet_volume / 2 + 2 * (missing_sixth + missing_half) + 8 * 4);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest3D, tet_identity_projector)
{
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 3

shapes:
- name: tet_shape
  material: {}
  geometry:
    format: stl
    path: {}
)";

  const std::string tet_material = "steel";
  const std::string tet_path = axom::fmt::format("{}/quest/tetrahedron.stl", AXOM_DATA_DIR);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template, tet_material, tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName());

  // check that we can set several projectors in 2D and 3D
  // uses simplest projectors, e.g. identity in 2D and 3D
  this->m_shaper->setPointProjector([](const Point3D& pt) {
    return Point3D {pt[0], pt[1], pt[2]};
  });
  this->m_shaper->setPointProjector([](const Point2D& pt) { return Point2D {pt[0], pt[1]}; });
  this->m_shaper->setPointProjector([](const Point3D& pt) { return Point2D {pt[0], pt[1]}; });
  this->m_shaper->setPointProjector([](const Point2D& pt) { return Point3D {pt[0], pt[1], 0}; });

  this->runShaping();

  // Check that the result has a volume fraction field associated with the tetrahedron material
  // The tet lives in cube of edge length 2 (and volume 8) and is defined by opposite corners.
  // It occupies 1/3 of the cube's volume
  constexpr double expected_volume = 8. / 3.;
  this->checkExpectedVolumeFractions(tet_material, expected_volume);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest3D, tet_doubling_projector)
{
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 3

shapes:
- name: tet_shape
  material: {}
  geometry:
    format: stl
    path: {}
)";

  const std::string tet_material = "steel";
  const std::string tet_path = axom::fmt::format("{}/quest/tetrahedron.stl", AXOM_DATA_DIR);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template, tet_material, tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName());

  // scale input points by a factor of 1/2 in each dimension
  this->m_shaper->setPointProjector([](const Point3D& pt) {
    return Point3D {pt[0] / 2, pt[1] / 2, pt[2] / 2};
  });

  // for good measure, add a 3D->2D projector that will not be used
  this->m_shaper->setPointProjector([](const Point3D&) { return Point2D {0, 0}; });

  this->runShaping();

  // Check that the result has a volume fraction field associated with the tetrahedron material
  // Scaling by a factor of 1/2 in each dimension should multiply the total volume by a factor of 8
  constexpr double orig_tet_volume = 8. / 3.;
  this->checkExpectedVolumeFractions(tet_material, 8 * orig_tet_volume);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest3D, circle_2D_projector)
{
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  constexpr double radius = 1.5;

  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: circle_shape
  material: {}
  geometry:
    format: c2c
    path: {}
    units: cm
    operators:
      - scale: {}
)";

  const std::string circle_material = "steel";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname),
                                   unit_semicircle_contour);

  ScopedTemporaryFile shape_file(
    axom::fmt::format("{}.yaml", testname),
    axom::fmt::format(shape_template, circle_material, contour_file.getFileName(), radius));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName());

  // set projector from 3D points to axisymmetric plane
  this->m_shaper->setPointProjector([](Point3D pt) {
    const double& x = pt[0];
    const double& y = pt[1];
    const double& z = pt[2];
    return Point2D {z, sqrt(x * x + y * y)};
  });

  // we need a higher quadrature order to resolve this shape at the (low) testing resolution
  this->m_shaper->setQuadratureOrder(8);

  this->runShaping();

  // Check that the result has a volume fraction field associated with the circle material
  constexpr double exp_volume = 4. / 3. * M_PI * radius * radius * radius;
  this->checkExpectedVolumeFractions(circle_material, exp_volume, 3e-2);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest3D, contour_and_stl_3D)
{
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  constexpr double radius = 1.5;

  const std::string shape_template = R"(
dimensions: 2

shapes:
# shape in a revolved sphere given as a c2c contour
- name: circle_shape
  material: {3}
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: {2}
# shape in a sphere given as an stl surface mesh
- name: sphere_shape
  material: {4}
  geometry:
    format: stl
    path: {1}
    units: cm
)";

  const std::string circle_material = "steel";
  const std::string sphere_material = "void";
  const std::string sphere_path = axom::fmt::format("{}/quest/unit_sphere.stl", AXOM_DATA_DIR);

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname),
                                   unit_semicircle_contour);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template,
                                                   contour_file.getFileName(),
                                                   sphere_path,
                                                   radius,
                                                   circle_material,
                                                   sphere_material));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());
  this->initializeShaping(shape_file.getFileName());

  // set projector from 3D points to axisymmetric plane
  this->m_shaper->setPointProjector([](Point3D pt) {
    const double& x = pt[0];
    const double& y = pt[1];
    const double& z = pt[2];
    return Point2D {z, sqrt(x * x + y * y)};
  });

  // we need a higher quadrature order to resolve this shape at the (low) testing resolution
  this->m_shaper->setQuadratureOrder(8);

  this->runShaping();

  // Check that the result has a volume fraction field associated with sphere and circle materials
  constexpr double exp_volume_contour = 4. / 3. * M_PI * radius * radius * radius;
  constexpr double exp_volume_sphere = 4. / 3. * M_PI * 1. * 1. * 1.;
  this->checkExpectedVolumeFractions(circle_material, exp_volume_contour - exp_volume_sphere, 3e-2);
  this->checkExpectedVolumeFractions(sphere_material, exp_volume_sphere, 3e-2);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;
#ifdef AXOM_USE_MPI
  // This is needed because of Axom's c2c reader.
  MPI_Init(&argc, &argv);
  // See if this aborts right away.
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#endif

  ::testing::InitGoogleTest(&argc, argv);
  slic::SimpleLogger logger(slic::message::Info);

  result = RUN_ALL_TESTS();
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif
  return result;
}
