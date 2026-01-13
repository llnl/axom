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
namespace fs = axom::utilities::filesystem;

namespace
{
const std::string unit_circle_contour =
  "piece = circle(origin=(0cm, 0cm), radius=1cm, start=0deg, end=360deg)";

const std::string unit_semicircle_contour = R"(
  piece = circle(origin=(0cm, 0cm), radius=1cm, start=0deg, end=180deg)
  piece = line(end=(0cm, 1cm)))";

const std::string proe_tet_fmt_str = R"(
4 1
1 {} {} {}
2 {} {} {}
3 {} {} {}
4 {} {} {}
1 1 2 3 4
)";

// Set the following to true for verbose output and for saving vis files
constexpr bool very_verbose_output = false;

// Utility function to slice a tetrahedron along a plane
primal::Polygon<double, 3> slice(const primal::Tetrahedron<double, 3>& tet,
                                 const primal::Plane<double, 3>& plane)
{
  primal::Polygon<double, 3> intersectionPolygon;

  // find intersection vertices
  for(int i = 0; i < 4; ++i)
  {
    for(int j = i + 1; j < 4; ++j)
    {
      primal::Segment<double, 3> edge(tet[i], tet[j]);
      double t {};
      if(primal::intersect(plane, edge, t))
      {
        intersectionPolygon.addVertex(edge.at(t));
      }
    }
  }
  SLIC_ASSERT(intersectionPolygon.numVertices() <= 4);

  // fix the polygon if it bowties
  if(intersectionPolygon.numVertices() == 4)
  {
    // note: using BezierCurve since Axom doesn't currently have intersect(segment, segment)
    primal::BezierCurve<double, 2> seg1(1);
    seg1[0] = Point2D(intersectionPolygon[0][0], intersectionPolygon[0][1]);
    seg1[1] = Point2D(intersectionPolygon[1][0], intersectionPolygon[1][1]);
    primal::BezierCurve<double, 2> seg2(1);
    seg2[0] = Point2D(intersectionPolygon[2][0], intersectionPolygon[2][1]);
    seg2[1] = Point2D(intersectionPolygon[3][0], intersectionPolygon[3][1]);
    axom::Array<double> sp, tp;

    if(!primal::intersect(seg1, seg2, sp, tp))
    {
      axom::utilities::swap(intersectionPolygon[2], intersectionPolygon[3]);
    }
  }
  return intersectionPolygon;
}

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

  axom::quest::SamplingShaper* getSamplingShaper() { return m_shaper.get(); }

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
    const auto policy = axom::runtime_policy::Policy::seq;
    const auto alloc = axom::policyToDefaultAllocatorID(policy);
    m_shaper = std::make_unique<quest::SamplingShaper>(policy, alloc, *m_shapeSet, &m_dc);
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

  void resetShaping()
  {
    std::vector<std::string> dereg;
    for(const auto& kv : m_dc.GetFieldMap())
    {
      if(axom::utilities::string::startsWith(kv.first, "vol_"))
      {
        dereg.push_back(kv.first);
      }
    }
    for(const auto& fld : dereg)
    {
      m_dc.DeregisterField(fld);
    }

    m_shaper.reset();
    m_shapeSet.reset();
  }

  /// Runs the shaping query over a shapefile; must be called after initializeShaping()
  void runShaping()
  {
    EXPECT_NE(nullptr, m_shaper) << "Shaper needs to be initialized via initializeShaping()";

    // Define lambda to override default dimensions, when necessary
    auto getShapeDim = [](const auto& shape) {
      static std::map<std::string, klee::Dimensions> format_dim = {{"c2c", klee::Dimensions::Two},
                                                                   {"stl", klee::Dimensions::Three}};

      const auto& shape_dim = shape.getGeometry().getInputDimensions();
      const auto& format_str = shape.getGeometry().getFormat();
      return format_dim.find(format_str) != format_dim.end() ? format_dim[format_str] : shape_dim;
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

    mfem::GridFunction* vf =
      axom::quest::shaping::getOrAllocateL2GridFunction(&m_dc,
                                                        name,
                                                        vfOrder,
                                                        dim,
                                                        mfem::BasisType::Positive);

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

/// Test fixture for SamplingShaper tests on 2D MFEM meshes
class SampleTester2D : public SamplingShaperTest
{
public:
  using Point2D = primal::Point<double, 2>;
  using BBox2D = primal::BoundingBox<double, 2>;

public:
  virtual ~SampleTester2D() { }

  void SetUp() override
  {
    // create a single element mesh in the unit 2D square
    const int polynomialOrder = 1;
    const BBox2D bbox({0, 0}, {1, 1});
    const axom::NumericArray<int, 2> celldims {1, 1};

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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_circle_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(
    axom::fmt::format(axom::fmt::runtime(shape_template), circle_material, contour_file.getPath()));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());
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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_circle_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(
    axom::fmt::format(axom::fmt::runtime(shape_template), circle_material, contour_file.getPath()));

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  // check that we can set several projectors in 2D and 3D
  // uses simplest projectors, e.g. identity in 2D and 3D
  this->m_shaper->setPointProjector33([](const Point3D& pt) { return Point3D {pt[0], pt[1], pt[2]}; });
  this->m_shaper->setPointProjector22([](const Point2D& pt) { return Point2D {pt[0], pt[1]}; });
  this->m_shaper->setPointProjector32([](const Point3D& pt) { return Point2D {pt[0], pt[1]}; });
  this->m_shaper->setPointProjector23([](const Point2D& pt) { return Point3D {pt[0], pt[1], 0}; });

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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_circle_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(
    axom::fmt::format(axom::fmt::runtime(shape_template), circle_material, contour_file.getPath()));

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  // check that we can set several projectors in 2D and 3D
  // creating an ellipse by scaling input x and y by scale_a and scale_b
  constexpr double scale_a = 3. / 2.;
  constexpr double scale_b = 3. / 4.;
  this->m_shaper->setPointProjector22(
    [](const Point2D& pt) { return Point2D {pt[0] / scale_a, pt[1] / scale_b}; });
  // check that we can register another projector that's not used
  this->m_shaper->setPointProjector33([](const Point3D&) { return Point3D {0., 0.}; });

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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_circle_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), contour_file.getPath()));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());

  this->initializeShaping(shape_file.getPath());
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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_circle_contour);

  // Set background material to 'void' (which is not present elsewhere)
  {
    fs::TempFile shape_file(testname, ".yaml");
    shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template),
                                       contour_file.getPath(),
                                       "void",
                                       "disk",
                                       "hole"));

    // Create an initial background material set to 1 everywhere
    std::map<std::string, mfem::GridFunction*> initialGridFunctions;
    {
      auto* vf = this->registerVolFracGridFunction("init_vf_bg");
      this->initializeVolFracGridFunction<2>(vf,
                                             [](int, const Point2D&, int) -> double { return 1.; });
      initialGridFunctions["void"] = vf;
    }

    this->validateShapeFile(shape_file.getPath());
    this->initializeShaping(shape_file.getPath(), initialGridFunctions);
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
    fs::TempFile shape_file(testname, ".yaml");
    shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template),
                                       contour_file.getPath(),
                                       "void",
                                       "disk",
                                       "void"));

    // Create an initial background material set to 1 everywhere
    std::map<std::string, mfem::GridFunction*> initialGridFunctions;
    {
      auto* vf = this->registerVolFracGridFunction("init_vf_bg");
      this->initializeVolFracGridFunction<2>(vf,
                                             [](int, const Point2D&, int) -> double { return 1.; });
      initialGridFunctions["void"] = vf;
    }

    this->validateShapeFile(shape_file.getPath());
    this->initializeShaping(shape_file.getPath(), initialGridFunctions);
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

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), "void", "left", "odds"));

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

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath(), initialGridFunctions);
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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_circle_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), contour_file.getPath()));

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
  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath(), initialGridFunctions);
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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_semicircle_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template),
                                     contour_file.getPath(),
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

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_circle_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template),
                                     contour_file.getPath(),
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

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath(), initialGridFunctions);

  // set projector from 2D mesh points to 3D query points within STL
  this->m_shaper->setPointProjector23([](Point2D pt) { return Point3D {pt[0], pt[1], 0.}; });

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

TEST_F(SamplingShaperTest2D, contour_and_mfem_2D)
{
  using Point2D = primal::Point<double, 2>;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  // Shape file
  const std::string shape_template = R"(
dimensions: 2

shapes:
# Background material
- name: bg
  material: luminiferous_ether
  geometry:
    format: none
# shape in a square using mfem
- name: square1
  material: square
  geometry:
    format: mfem
    path: {0}
    units: cm
)";

  // MFEM file
  const std::string mfem_square_contour = R"(
MFEM NURBS mesh v1.0

# MFEM Geometry Types (see fem/geom.hpp):
#
# SEGMENT = 1 | SQUARE = 3 | CUBE = 5
#
# element: <attr> 1 <v0> <v1>
# edge: <idx++> 0 1  <-- idx increases by one each time
# knotvector: <order> <num_ctrl_pts> [knots]; sizeof(knots) is 1+order+num_ctrl_pts
# weights: array of weights corresponding to the NURBS element
# FES: list of control points; vertex control points at top, then interior control points

dimension
1

elements
4
1 1 0 1
1 1 2 3
1 1 4 5
1 1 6 7

boundary
0

edges
4
0 0 1
1 0 1
2 0 1
3 0 1

vertices
8

knotvectors
4
3 4 0 0 0 0 1 1 1 1
3 4 0 0 0 0 1 1 1 1
3 4 0 0 0 0 1 1 1 1
3 4 0 0 0 0 1 1 1 1

weights
1 1
1 1
1 1
1 1
1 1
1 1
1 1
1 1

FiniteElementSpace
FiniteElementCollection: NURBS
VDim: 2
Ordering: 1
-1.0 -1.0 1.0 -1.0
1.0 -1.0 1.0 1.0
1.0 1.0 -1.0 1.0
-1.0 1.0 -1.0 -1.0
0.6 -1.0 -0.4 -1.0
1.0 0.6 1.0 -0.4
-0.4 1.0 0.6 1.0
-1.0 -0.4 -1.0 0.6
)";

  const std::string background_material = "luminiferous_ether";
  const std::string square_material = "square";

  fs::TempFile contour_file(testname, ".mesh");
  contour_file.write(mfem_square_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), contour_file.getPath()));

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

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath(), initialGridFunctions);
  // Use WindingNumber shaping!
  this->m_shaper->setSamplingMethod(quest::SamplingShaper::SamplingMethod::WindingNumber);

  this->m_shaper->setQuadratureOrder(8);
  this->runShaping();

  // Check that the result has a volume fraction field associated with square materials
  constexpr double exp_volume_square = 4.;
  this->checkExpectedVolumeFractions(square_material, exp_volume_square, 1.e-4);

  for(const auto& vf_name : {background_material, square_material})
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

TEST_F(SamplingShaperTest3D, basic_tet_boundary)
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

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), tet_material, tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Bounding box of 3D input mesh: \n" << this->meshBoundingBox());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());
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

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), tet_material, tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());

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

  this->initializeShaping(shape_file.getPath(), initialGridFunctions);
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

TEST_F(SamplingShaperTest3D, tet_boundary_preshaped_with_replacements)
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

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());

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

  this->initializeShaping(shape_file.getPath(), initialGridFunctions);
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

TEST_F(SamplingShaperTest3D, tet_boundary_identity_projector)
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

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), tet_material, tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  // check that we can set several projectors in 2D and 3D
  // uses simplest projectors, e.g. identity in 2D and 3D
  this->m_shaper->setPointProjector33([](const Point3D& pt) { return Point3D {pt[0], pt[1], pt[2]}; });
  this->m_shaper->setPointProjector22([](const Point2D& pt) { return Point2D {pt[0], pt[1]}; });
  this->m_shaper->setPointProjector32([](const Point3D& pt) { return Point2D {pt[0], pt[1]}; });
  this->m_shaper->setPointProjector23([](const Point2D& pt) { return Point3D {pt[0], pt[1], 0}; });

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

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), tet_material, tet_path));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  // scale input points by a factor of 1/2 in each dimension
  this->m_shaper->setPointProjector33(
    [](const Point3D& pt) { return Point3D {pt[0] / 2, pt[1] / 2, pt[2] / 2}; });

  // for good measure, add a 3D->2D projector that will not be used
  this->m_shaper->setPointProjector32([](const Point3D&) { return Point2D {0, 0}; });

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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_semicircle_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template),
                                     circle_material,
                                     contour_file.getPath(),
                                     radius));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  // set projector from 3D points to axisymmetric plane
  this->m_shaper->setPointProjector32([](Point3D pt) {
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

  fs::TempFile contour_file(testname, ".contour");
  contour_file.write(unit_semicircle_contour);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template),
                                     contour_file.getPath(),
                                     sphere_path,
                                     radius,
                                     circle_material,
                                     sphere_material));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  // set projector from 3D points to axisymmetric plane
  this->m_shaper->setPointProjector32([](Point3D pt) {
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

TEST_F(SamplingShaperTest2D, shape_proe_tet_with_2D_projection)
{
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;

  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 3

shapes:
- name: tet_shape
  material: {1}
  geometry:
    format: proe
    path: {0}
    units: cm
)";

  // regular tet w/ vertices at corners of cube
  const auto tet = primal::Tetrahedron<double, 3> {Point3D {-1, -1, -1},
                                                   Point3D {1, 1, -1},
                                                   Point3D {-1, -1, 1},
                                                   Point3D {-1, 1, 1}};

  // Check that the volume of the tetrahedron is 4/3
  double tetVolume = tet.volume();
  constexpr double expectedTetVolume = 4.0 / 3.0;
  EXPECT_NEAR(tetVolume, expectedTetVolume, 1e-6);
  SLIC_INFO(axom::fmt::format("Computed tetrahedron volume: {}", tetVolume));

  // Write out tet as a proe file
  // clang-format off
  const std::string proe_tet_str 
    = axom::fmt::format(axom::fmt::runtime(proe_tet_fmt_str), tet[0][0], tet[0][1], tet[0][2],
                                          tet[1][0], tet[1][1], tet[1][2],
                                          tet[2][0], tet[2][1], tet[2][2],
                                          tet[3][0], tet[3][1], tet[3][2]);
  // clang-format on

  const std::string tet_material = "steel";
  fs::TempFile tet_file(testname, ".proe");
  tet_file.write(proe_tet_str);

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(
    axom::fmt::format(axom::fmt::runtime(shape_template), tet_file.getPath(), tet_material));

  if(very_verbose_output)
  {
    SLIC_INFO("Tet file: \n" << tet_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getPath());

  // exercise the point projector by running at several XY planes
  for(const double z : {-1.2, -.9, -.75, -.1, 0., .25, 1. / 3., .5})
  {
    this->resetShaping();

    this->initializeShaping(shape_file.getPath());

    primal::Plane<double, 3> plane({0, 0, 1}, z);
    const auto polygon = slice(tet, plane);
    const double intersectionArea = polygon.area();
    SLIC_INFO(axom::fmt::format("Area of intersection polygon: {}", intersectionArea));

    // set projector from 2D points to 3-space, z-coord is lambda captured
    this->m_shaper->setPointProjector23([z](Point2D pt) -> Point3D {
      const double& x = pt[0];
      const double& y = pt[1];
      return Point3D {x, y, z};
    });

    // we need a higher quadrature order to resolve this shape at the (low) testing resolution
    this->m_shaper->setQuadratureOrder(8);

    this->runShaping();

    // Check that the result has a volume fraction field associated with sphere and circle materials
    this->checkExpectedVolumeFractions(tet_material, intersectionArea, 3e-2);

    // Save meshes and fields
    if(very_verbose_output)
    {
      this->getDC().Save(axom::fmt::format("{}_{}", testname, z),
                         axom::sidre::Group::getDefaultIOProtocol());
    }
  }
}

//-----------------------------------------------------------------------------

TEST_F(SampleTester2D, check_bbox_inouts)
{
  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: {}
  material: {}
  geometry:
    format: c2c
    path: {}
    units: cm
)";

  // inputs: {x_min} {y_min} {x_max} {y_max}
  const std::string rectangle_contour = R"(
point = start
piece = line(start=({0}cm, {1}cm), end=({0}cm, {3}cm))
piece = line()
piece = line(start=({2}cm, {3}cm), end=({2}cm, {1}cm))
piece = line(end=start)
)";

  const std::string rect_shape = "rectShape";
  const std::string rect_material = "rectMat";

  using Pt = axom::primal::Point<double, 2>;
  using BBox = axom::primal::BoundingBox<double, 2>;

  // Test a few rectangular shapes
  constexpr double EPS = 1e-4;
  for(auto bb : {BBox {Pt {0., 0.}, Pt {1., 1.}},

                 BBox {Pt {0., 0.}, Pt {.5, .5}},
                 BBox {Pt {0., 0.}, Pt {.5, 1.}},
                 BBox {Pt {0., 0.}, Pt {1., .5}},

                 BBox {Pt {.5, 0.}, Pt {1., 1.}},
                 BBox {Pt {0., .5}, Pt {1., 1.}},
                 BBox {Pt {.5, .5}, Pt {1., 1.}}})
  {
    bb.expand(EPS);  // expand slightly to catch quadrature points on boundary of bbox

    const std::string contour_str = axom::fmt::format(axom::fmt::runtime(rectangle_contour),
                                                      bb.getMin()[0],
                                                      bb.getMin()[1],
                                                      bb.getMax()[0],
                                                      bb.getMax()[1]);

    const std::string testname =
      axom::fmt::format("{}_{}_{}_{}_{}",
                        ::testing::UnitTest::GetInstance()->current_test_info()->name(),
                        bb.getMin()[0],
                        bb.getMin()[1],
                        bb.getMax()[0],
                        bb.getMax()[1]);

    fs::TempFile contour_file(testname, ".contour");
    contour_file.write(contour_str);

    fs::TempFile shape_file(testname, ".yaml");
    shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template),
                                       rect_shape,
                                       rect_material,
                                       contour_file.getPath()));

    if(very_verbose_output)
    {
      SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
      SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
    }

    this->validateShapeFile(shape_file.getPath());

    for(int qorder : {3, 4, 5, 6, 7, 8, 9})
    {
      this->resetShaping();
      this->initializeShaping(shape_file.getPath());

      this->m_shaper->setVolumeFractionOrder(0);
      this->m_shaper->setQuadratureOrder(qorder);

      this->runShaping();

      // grab the inout quadrature data
      auto& mesh = this->getMesh();
      EXPECT_EQ(mesh.GetNE(), 1);

      auto* inout =
        this->getSamplingShaper()->getShapeQFunction(axom::fmt::format("inout_{}", rect_shape));
      EXPECT_NE(inout, nullptr);

      auto* qfs = dynamic_cast<mfem::QuadratureSpace*>(inout->GetSpace());
      EXPECT_NE(qfs, nullptr);
      auto* T = qfs->GetTransformation(0);
      const mfem::IntegrationRule& ir = qfs->GetElementIntRule(0);

      mfem::Vector inout_data;
      inout->GetValues(0, inout_data);
      EXPECT_EQ(inout_data.Size(), ir.GetNPoints());

      // check that inout values match corresponding bbox containment
      for(int i = 0; i < ir.GetNPoints(); i++)
      {
        const mfem::IntegrationPoint& ip = ir.IntPoint(i);
        mfem::Vector ip_phys(2);
        T->Transform(ip, ip_phys);
        const Pt p {ip_phys[1], ip_phys[0]};

        EXPECT_TRUE(inout_data[i] == 1. || inout_data[i] == 0.);
        const bool is_in = inout_data[i] == 1 ? true : false;
        const bool exp_in = bb.contains(p) ? true : false;
        EXPECT_EQ(is_in, exp_in)
          << axom::fmt::format(
               "Qorder: {}, Quadrature point {}:  physical position: {}, inout: {}, expected: {}. ",
               qorder,
               i,
               p,
               is_in,
               exp_in)
          << axom::fmt::format("Bounding box is: {}", bb);
      }

      // Save meshes and fields
      if(very_verbose_output)
      {
        this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
      }
    }
  }
}

//-----------------------------------------------------------------------------

TEST_F(SamplingShaperTest2D, loadShape_missing_c2c_file_aborts)
{
  // Tests Klee shape file referencing non-existant c2c file; should fail
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: missing_c2c
  material: mat
  geometry:
    format: c2c
    path: {}
)";

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), "missing.contour"));

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  EXPECT_TRUE(m_shapeSet);
  EXPECT_TRUE(m_shaper);
  EXPECT_FALSE(m_shapeSet->getShapes().empty());

  const auto& shape = m_shapeSet->getShapes().front();
  slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(m_shaper->loadShape(shape), slic::SlicAbortException);
}

TEST_F(SamplingShaperTest2D, loadShape_missing_mfem_mesh_file_aborts)
{
  // Tests Klee shape file referencing non-existant mfem file; should fail
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: missing_mfem
  material: mat
  geometry:
    format: mfem
    path: {}
)";

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), "missing.mesh"));

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  EXPECT_TRUE(m_shapeSet);
  EXPECT_TRUE(m_shaper);
  EXPECT_FALSE(m_shapeSet->getShapes().empty());

  const auto& shape = m_shapeSet->getShapes().front();
  slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(m_shaper->loadShape(shape), slic::SlicAbortException);
}

TEST_F(SamplingShaperTest2D, loadShape_missing_mfem_mesh_file_windingnumber_aborts)
{
  // Tests Klee shape file referencing non-existant mfem file; should fail
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: missing_mfem_wn
  material: mat
  geometry:
    format: mfem
    path: {}
)";

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), "missing.mesh"));

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  EXPECT_TRUE(m_shapeSet);
  EXPECT_TRUE(m_shaper);
  m_shaper->setSamplingMethod(quest::SamplingShaper::SamplingMethod::WindingNumber);
  EXPECT_FALSE(m_shapeSet->getShapes().empty());

  const auto& shape = m_shapeSet->getShapes().front();
  slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(m_shaper->loadShape(shape), slic::SlicAbortException);
}

//-----------------------------------------------------------------------------

TEST_F(SamplingShaperTest3D, loadShape_missing_stl_file_aborts)
{
  // Tests Klee shape file referencing non-existant stl mesh; should fail
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 3

shapes:
- name: missing_stl
  material: mat
  geometry:
    format: stl
    path: {}
)";

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), "missing.stl"));

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  EXPECT_TRUE(m_shapeSet);
  EXPECT_TRUE(m_shaper);
  EXPECT_FALSE(m_shapeSet->getShapes().empty());

  const auto& shape = m_shapeSet->getShapes().front();
  slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(m_shaper->loadShape(shape), slic::SlicAbortException);
}

TEST_F(SamplingShaperTest3D, loadShape_missing_proe_file_aborts)
{
  // Tests Klee shape file referencing non-existant pro-e file; should fail
  const auto& testname = ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 3

shapes:
- name: missing_proe
  material: mat
  geometry:
    format: proe
    path: {}
)";

  fs::TempFile shape_file(testname, ".yaml");
  shape_file.write(axom::fmt::format(axom::fmt::runtime(shape_template), "missing.proe"));

  this->validateShapeFile(shape_file.getPath());
  this->initializeShaping(shape_file.getPath());

  EXPECT_TRUE(m_shapeSet);
  EXPECT_TRUE(m_shaper);
  EXPECT_FALSE(m_shapeSet->getShapes().empty());

  const auto& shape = m_shapeSet->getShapes().front();
  slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(m_shaper->loadShape(shape), slic::SlicAbortException);
}

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);

  ::testing::InitGoogleTest(&argc, argv);
  slic::SimpleLogger logger(slic::message::Info);

  const int result = RUN_ALL_TESTS();
  return result;
}
