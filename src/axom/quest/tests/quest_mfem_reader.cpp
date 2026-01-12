// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#if !defined(AXOM_USE_MFEM) || !defined(AXOM_USE_SIDRE)
  #error These tests should only be included when Axom is configured with MFEM and SIDRE
#endif

#include "axom/quest/io/MFEMReader.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include "mfem.hpp"

// gtest includes
#include "gtest/gtest.h"

#include <fstream>
#include <map>
#include <string>

// namespace aliases
namespace primal = axom::primal;
namespace quest = axom::quest;

//------------------------------------------------------------------------------
std::string pjoin(const std::string &str) { return str; }

std::string pjoin(const char *str) { return std::string(str); }

template <typename... Args>
std::string pjoin(const std::string &str, Args... args)
{
  return axom::utilities::filesystem::joinPath(str, pjoin(args...));
}

template <typename... Args>
std::string pjoin(const char *str, Args... args)
{
  return axom::utilities::filesystem::joinPath(std::string(str), pjoin(args...));
}
//------------------------------------------------------------------------------

namespace
{
using BezierCurve2D = primal::BezierCurve<double, 2>;
using Point2D = primal::Point<double, 2>;

void write_mesh_from_bezier_curves(const std::string &mesh_path,
                                   const axom::Array<BezierCurve2D> &bezier_curves,
                                   const axom::Array<int> &attributes,
                                   mfem::FiniteElementCollection &fec)
{
  ASSERT_EQ(bezier_curves.size(), attributes.size());
  const int num_curves = bezier_curves.size();
  ASSERT_GT(num_curves, 0);

  constexpr int VDIM = 2;

  mfem::Mesh mesh(/*Dim*/ 1,
                  /*NVert*/ VDIM * num_curves,
                  /*NElem*/ num_curves,
                  /*NBdrElem*/ 0,
                  /*spaceDim*/ VDIM);

  for(int i = 0; i < num_curves; ++i)
  {
    const auto &curve = bezier_curves[i];
    ASSERT_EQ(curve.getOrder(), fec.GetOrder());

    const auto &p0 = curve.getInitPoint();
    const auto &p1 = curve.getEndPoint();

    const double v0[] = {p0[0], p0[1]};
    const double v1[] = {p1[0], p1[1]};

    mesh.AddVertex(v0);
    mesh.AddVertex(v1);
    mesh.AddSegment(VDIM * i, VDIM * i + 1, attributes[i]);
  }

  mesh.FinalizeTopology(/*generate_bdr*/ true);
  mesh.Finalize(/*refine*/ false, /*fix_orientation*/ true);

  mfem::FiniteElementSpace fes(&mesh,
                               &fec,
                               /*vdim*/ mesh.SpaceDimension(),
                               mfem::Ordering::byVDIM);

  mfem::GridFunction nodes(&fes);
  nodes = 0.0;

  mfem::Array<int> dofs, vdofs_c;
  for(int e = 0; e < mesh.GetNE(); ++e)
  {
    // utility lambda to help with indexing; mfem stores endpoint and the interior points
    const int order = fes.GetOrder(e);
    auto mfemLocalToBezier = [order](int j) -> int {
      if(j == 0) return 0;      // start vertex
      if(j == 1) return order;  // end vertex
      return j - 1;             // interior dofs
    };

    fes.GetElementDofs(e, dofs);
    EXPECT_EQ(order, bezier_curves[e].getOrder());
    EXPECT_EQ(dofs.Size(), order + 1);

    for(int i = 0; i <= order; ++i)
    {
      const auto &cp = bezier_curves[e][mfemLocalToBezier(i)];
      nodes(fes.DofToVDof(dofs[i], 0)) = cp[0];
      nodes(fes.DofToVDof(dofs[i], 1)) = cp[1];
    }
  }

  mesh.NewNodes(nodes, /*make_owner*/ false);

  std::ofstream ofs(mesh_path);
  ASSERT_TRUE(ofs.good());
  ofs.precision(17);
  mesh.Print(ofs);
}
}  // namespace

TEST(quest_mfem_reader, read_nurbs_curves)
{
  const std::string fileName =
    pjoin(AXOM_DATA_DIR, "contours", "heroic_roses", "mfem", "blue0.mesh");

  quest::MFEMReader reader;
  reader.setFileName(fileName);

  // Read as 9 NURBS curves
  axom::Array<primal::NURBSCurve<double, 2>> curves;
  axom::Array<int> attributes;
  EXPECT_EQ(reader.read(curves, attributes), quest::MFEMReader::READ_SUCCESS);
  EXPECT_EQ(curves.size(), 9);
  ASSERT_EQ(attributes.size(), 9);
  for(int i = 0; i < attributes.size(); ++i)
  {
    EXPECT_EQ(attributes[i], 1);
  }
}

TEST(quest_mfem_reader, read_curved_polygon)
{
  const std::string fileName =
    pjoin(AXOM_DATA_DIR, "contours", "heroic_roses", "mfem", "blue0.mesh");

  quest::MFEMReader reader;
  reader.setFileName(fileName);

  // Read as 1 CurvedPolygon with 9 edges
  axom::Array<primal::CurvedPolygon<axom::primal::NURBSCurve<double, 2>>> polys;
  EXPECT_EQ(reader.read(polys), 0);
  EXPECT_EQ(polys.size(), 1);
  EXPECT_EQ(polys[0].numEdges(), 9);

  // Read as CurvedPolygon
  polys.clear();
  const std::string fileNameB =
    pjoin(AXOM_DATA_DIR, "contours", "heroic_roses", "mfem_cp", "black.mesh");
  reader.setFileName(fileNameB);
  EXPECT_EQ(reader.read(polys), 0);
  EXPECT_EQ(polys.size(), 73);
  // Pick some curved polygons and check lengths.
  EXPECT_EQ(polys[0].numEdges(), 11);
  EXPECT_EQ(polys[20].numEdges(), 20);
  EXPECT_EQ(polys[40].numEdges(), 4);
  EXPECT_EQ(polys[72].numEdges(), 2);
}

TEST(quest_mfem_reader, preserves_rational_weights)
{
  const std::string fileName =
    pjoin(AXOM_DATA_DIR, "contours", "heroic_roses", "mfem", "brightgreen_over.mesh");

  quest::MFEMReader reader;
  reader.setFileName(fileName);

  // test read on NURBSCurve array and check that at least one curve is rational
  {
    axom::Array<primal::NURBSCurve<double, 2>> curves;
    EXPECT_EQ(reader.read(curves), quest::MFEMReader::READ_SUCCESS);
    ASSERT_GT(curves.size(), 0);

    bool any_rational = false;
    for(const auto &curve : curves)
    {
      if(curve.isRational())
      {
        any_rational = true;
        break;
      }
    }
    EXPECT_TRUE(any_rational);
  }

  // do the same with the extracted Curved polygon array
  {
    axom::Array<primal::CurvedPolygon<axom::primal::NURBSCurve<double, 2>>> polys;
    EXPECT_EQ(reader.read(polys), quest::MFEMReader::READ_SUCCESS);
    ASSERT_GT(polys.size(), 0);

    bool any_rational = false;
    for(const auto &poly : polys)
    {
      for(const auto &cur : poly.getEdges())
      {
        if(cur.isRational())
        {
          any_rational = true;
          break;
        }
      }
    }
    EXPECT_TRUE(any_rational);
  }
}

TEST(quest_mfem_reader, read_bernstein_basis_roundtrip_bezier_order3)
{
  axom::utilities::filesystem::TempFile tmp_mesh("bernstein_basis", ".mesh");

  // write out mfem file containing a few Bezier curves; store curves and attributes in `expected`, keyed by attribute
  std::map<int, BezierCurve2D> expected;
  {
    axom::Array<BezierCurve2D> input_curves = {
      BezierCurve2D({Point2D {0.0, 0.0}, Point2D {0.25, 0.5}, Point2D {0.75, 0.5}, Point2D {1.0, 0.0}},
                    3),
      BezierCurve2D({Point2D {2.0, 0.0}, Point2D {2.2, -0.2}, Point2D {2.8, 0.2}, Point2D {3.0, 0.0}},
                    3),
      BezierCurve2D({Point2D {-1.0, 1.0}, Point2D {-0.5, 1.5}, Point2D {0.5, 0.5}, Point2D {1.0, 1.0}},
                    3)};
    axom::Array<int> input_attributes {30, 10, 20};

    constexpr int order = 3;
    mfem::H1Pos_FECollection fec(order, /*dim*/ 1);
    write_mesh_from_bezier_curves(tmp_mesh.getPath(), input_curves, input_attributes, fec);

    for(int i = 0; i < input_curves.size(); ++i)
    {
      expected.emplace(input_attributes[i], input_curves[i]);
    }
  }

  // we should be able to successfully load the file and compare to originals
  quest::MFEMReader reader;
  reader.setFileName(tmp_mesh.getPath());

  axom::Array<primal::NURBSCurve<double, 2>> curves;
  axom::Array<int> attributes;
  EXPECT_EQ(reader.read(curves, attributes), quest::MFEMReader::READ_SUCCESS);

  ASSERT_EQ(curves.size(), expected.size());
  ASSERT_EQ(attributes.size(), expected.size());

  for(int i = 0; i < curves.size(); ++i)
  {
    const int attr = attributes[i];
    const auto expected_it = expected.find(attr);
    ASSERT_TRUE(expected_it != expected.end());
    const auto &expected_curve = expected_it->second;

    EXPECT_EQ(curves[i].getDegree(), expected_curve.getOrder());
    EXPECT_EQ(curves[i].getNumControlPoints(), expected_curve.getOrder() + 1);
    for(int j = 0; j < expected_curve.getOrder() + 1; ++j)
    {
      EXPECT_NEAR(curves[i][j][0], expected_curve[j][0], 1e-12);
      EXPECT_NEAR(curves[i][j][1], expected_curve[j][1], 1e-12);
    }
  }
}

TEST(quest_mfem_reader, read_non_bernstein_basis_rejected)
{
  axom::utilities::filesystem::TempFile tmp_mesh("non_bernstein_basis", ".mesh");

  // create file containing curves in non-Bernstein basis
  {
    axom::Array<BezierCurve2D> input_curves = {BezierCurve2D(
      {Point2D {0.0, 0.0}, Point2D {0.25, 0.5}, Point2D {0.75, 0.5}, Point2D {1.0, 0.0}},
      3)};
    axom::Array<int> input_attributes = {7};

    constexpr int order = 3;
    mfem::H1_FECollection fec(order, /*dim*/ 1);
    write_mesh_from_bezier_curves(tmp_mesh.getPath(), input_curves, input_attributes, fec);
  }

  // attempt to read it in; should fail
  quest::MFEMReader reader;
  reader.setFileName(tmp_mesh.getPath());

  axom::Array<primal::NURBSCurve<double, 2>> curves;
  axom::Array<int> attributes;
  EXPECT_EQ(reader.read(curves, attributes), quest::MFEMReader::READ_FAILED);
}

TEST(quest_mfem_reader, read_curved_polygon_noncontiguous_attributes)
{
  axom::utilities::filesystem::TempFile tmp_mesh("noncontiguous_attributes", ".mesh");

  // This test uses non-contiguous attributes.
  // For testing, we're setting the y-coordinate to be the same as the attribute
  constexpr int attr10 {10};
  constexpr int attr20 {20};
  constexpr int order = 1;
  constexpr int dim = 1;

  // write out a mesh containing two segments w/ non-contiguous attributes
  {
    axom::Array<BezierCurve2D> input_curves = {
      BezierCurve2D({Point2D {0, attr20}, Point2D {1, attr20}}, order),
      BezierCurve2D({Point2D {0, attr10}, Point2D {1, attr10}}, order)};

    axom::Array<int> input_attributes {attr20, attr10};

    mfem::H1Pos_FECollection fec(order, dim);
    write_mesh_from_bezier_curves(tmp_mesh.getPath(), input_curves, input_attributes, fec);
  }

  quest::MFEMReader reader;
  reader.setFileName(tmp_mesh.getPath());

  // check that we can successfully read in the attributes to CurvedPolygon array
  {
    axom::Array<primal::CurvedPolygon<axom::primal::NURBSCurve<double, 2>>> polys;
    axom::Array<int> attributes;
    EXPECT_EQ(reader.read(polys, attributes), quest::MFEMReader::READ_SUCCESS);

    ASSERT_EQ(polys.size(), 2);
    ASSERT_EQ(attributes.size(), 2);

    EXPECT_EQ(polys[0].numEdges(), 1);
    EXPECT_EQ(polys[1].numEdges(), 1);

    // note: the curves are added to a map, so the order is not guaranteed
    // let's check that the attributes and geometry match expectations
    // the y-coordinates of the edges start and end vertex should equal the attribute
    for(int i : {0, 1})
    {
      const auto &curve = polys[i][0];
      switch(attributes[i])
      {
      case attr10:
        EXPECT_EQ(curve[0][1], attr10);
        EXPECT_EQ(curve[1][1], attr10);
        break;
      case attr20:
        EXPECT_EQ(curve[0][1], attr20);
        EXPECT_EQ(curve[1][1], attr20);
        break;
      default:
        FAIL() << "Got unexpected attribute for polygon " << i << ": " << attributes[i] << "\n";
        break;
      }
    }
  }

  // check that we can successfully read in the attributes to NURBSCurve array
  {
    axom::Array<primal::NURBSCurve<double, 2>> curves;
    axom::Array<int> attributes;
    EXPECT_EQ(reader.read(curves, attributes), quest::MFEMReader::READ_SUCCESS);
    ASSERT_EQ(curves.size(), 2);
    ASSERT_EQ(attributes.size(), 2);

    for(int i : {0, 1})
    {
      const auto &curve = curves[i];
      switch(attributes[i])
      {
      case attr10:
        EXPECT_EQ(curve[0][1], attr10);
        EXPECT_EQ(curve[1][1], attr10);
        break;
      case attr20:
        EXPECT_EQ(curve[0][1], attr20);
        EXPECT_EQ(curve[1][1], attr20);
        break;
      default:
        FAIL() << "Got unexpected attribute for curve " << i << ": " << attributes[i] << "\n";
        break;
      }
    }
  }
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
