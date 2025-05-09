// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_winding_number.cpp
 * \brief Example that computes the winding number of a grid of points
 * against a collection of 2D parametric rational curves.
 * Supports MFEM meshes in the cubic positive Bernstein basis or the (rational)
 * NURBS basis.
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "mfem.hpp"

namespace primal = axom::primal;
using Point2D = primal::Point<double, 2>;
using BezierCurve2D = primal::BezierCurve<double, 2>;
using BoundingBox2D = primal::BoundingBox<double, 2>;

/*!
 * Given an mfem mesh, convert element with id \a elem_id to a (rational) BezierCurve
 * \pre Assumes the elements of the mfem mesh are in the positive (Bernstein)
 * basis, or in the NURBS basis
 */
BezierCurve2D segment_to_curve(const mfem::Mesh* mesh, int elem_id)
{
  const auto* fes = mesh->GetNodes()->FESpace();
  const auto* fec = fes->FEColl();

  const bool isBernstein = dynamic_cast<const mfem::H1Pos_FECollection*>(fec) != nullptr;
  const bool isNURBS = dynamic_cast<const mfem::NURBSFECollection*>(fec) != nullptr;

  SLIC_ERROR_IF(!(isBernstein || isNURBS),
                "MFEM mesh elements must be in either the Bernstein or NURBS basis");

  const int NE = isBernstein ? mesh->GetNE() : fes->GetNURBSext()->GetNP();
  SLIC_ERROR_IF(NE < elem_id, axom::fmt::format("Mesh does not have {} elements", elem_id));

  const int order = isBernstein ? fes->GetOrder(elem_id) : mesh->NURBSext->GetOrders()[elem_id];
  SLIC_ERROR_IF(order != 3,
                axom::fmt::format("This example currently requires the input mfem mesh to "
                                  "contain cubic elements, but the order of element {} is {}",
                                  elem_id,
                                  order));

  mfem::Array<int> dofs;
  mfem::Array<int> vdofs;

  mfem::Vector dvec;
  mfem::Vector v;

  fes->GetElementDofs(elem_id, dofs);
  fes->GetElementVDofs(elem_id, vdofs);
  mesh->GetNodes()->GetSubVector(vdofs, v);

  // Currently hard-coded for 3rd order. This can easily be extended to arbitrary order
  axom::Array<Point2D> points(4, 4);
  if(isBernstein)
  {
    points[0] = Point2D {v[0], v[0 + 4]};
    points[1] = Point2D {v[2], v[2 + 4]};
    points[2] = Point2D {v[3], v[3 + 4]};
    points[3] = Point2D {v[1], v[1 + 4]};

    return BezierCurve2D(points, fec->GetOrder());
  }
  else  // isNURBS
  {
    // temporary assumption is that there are no interior knots
    // i.e. the NURBS curve is essentially a rational Bezier curve

    points[0] = Point2D {v[0], v[0 + 4]};
    points[1] = Point2D {v[1], v[1 + 4]};
    points[2] = Point2D {v[2], v[2 + 4]};
    points[3] = Point2D {v[3], v[3 + 4]};

    fes->GetNURBSext()->GetWeights().GetSubVector(dofs, dvec);
    axom::Array<double> weights {dvec[0], dvec[1], dvec[2], dvec[3]};

    return BezierCurve2D(points, weights, fec->GetOrder());
  }
}

bool check_mesh_valid(const mfem::Mesh* mesh)
{
  const auto* fes = mesh->GetNodes()->FESpace();
  if(fes == nullptr)
  {
    SLIC_WARNING("MFEM mesh finite element space was null");
    return false;
  }

  const auto* fec = fes->FEColl();
  if(fec == nullptr)
  {
    SLIC_WARNING("MFEM mesh finite element collection was null");
    return false;
  }

  const bool isBernstein = dynamic_cast<const mfem::H1Pos_FECollection*>(fec) != nullptr;
  const bool isNURBS = dynamic_cast<const mfem::NURBSFECollection*>(fec) != nullptr;
  const bool isValidFEC = isBernstein || isNURBS;

  // TODO: Convert from Lagrange to Bernstein, if/when necessary
  if(!isValidFEC)
  {
    SLIC_WARNING(
      "Example only currently supports 1D NURBS meshes "
      "or meshes with nodes in the Bernstein basis");
    return false;
  }

  if(fes->GetVDim() != 2)
  {
    SLIC_WARNING("Example only currently supports 2D meshes");
    return false;
  }

  const int NE = isBernstein ? mesh->GetNE() : fes->GetNURBSext()->GetNP();
  int order = -1;
  if(isBernstein)
  {
    order = NE > 0 ? fes->GetOrder(0) : 3;
  }
  else  // isNURBS
  {
    //SLIC_INFO("nurbsext order :" << mesh->NURBSext->GetOrder());
    order = NE > 0 ? mesh->NURBSext->GetOrders()[0] : 3;
  }

  if(order != 3)
  {
    SLIC_WARNING(
      axom::fmt::format("This example currently requires the input mfem mesh to contain cubic "
                        "elements, but the provided mesh has order {}",
                        order));
    return false;
  }

  return true;
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger raii_logger;

  axom::CLI::App app {
    "Load mesh containing collection of curves"
    " and optionally generate a query mesh of winding numbers."};

  std::string inputFile;
  std::string outputPrefix = {"winding"};

  bool verbose {false};

  // Query mesh parameters
  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int queryOrder {1};

  app.add_option("-i,--input", inputFile)
    ->description("MFEM mesh containing contours (1D segments)")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_option("-o,--output-prefix", outputPrefix)
    ->description(
      "Prefix for output 2D query mesh (in MFEM format) mesh containing "
      "winding number calculations")
    ->capture_default_str();

  app.add_flag("-v,--verbose", verbose, "verbose output")->capture_default_str();

  auto* query_mesh_subcommand =
    app.add_subcommand("query_mesh")->description("Options for setting up a query mesh")->fallthrough();
  query_mesh_subcommand->add_option("--min", boxMins)
    ->description("Min bounds for box mesh (x,y)")
    ->expected(2)
    ->required();
  query_mesh_subcommand->add_option("--max", boxMaxs)
    ->description("Max bounds for box mesh (x,y)")
    ->expected(2)
    ->required();
  query_mesh_subcommand->add_option("--res", boxResolution)
    ->description("Resolution of the box mesh (i,j)")
    ->expected(2)
    ->required();
  query_mesh_subcommand->add_option("--order", queryOrder)
    ->description("polynomial order of the query mesh")
    ->check(axom::CLI::PositiveNumber);

  CLI11_PARSE(app, argc, argv);

  mfem::Mesh mesh(inputFile);
  SLIC_INFO(
    axom::fmt::format("Curve mesh has a topological dimension of {}d, "
                      "has {} vertices and {} elements",
                      mesh.Dimension(),
                      mesh.GetNV(),
                      mesh.GetNE()));

  if(!check_mesh_valid(&mesh))
  {
    return 1;
  }

  axom::Array<int> segments;
  axom::Array<BezierCurve2D> curves;

  // Loop through mesh elements, retaining the (curved) 1D segments
  for(int i = 0; i < mesh.GetNE(); ++i)
  {
    auto* el = mesh.GetElement(i);
    if(el->GetGeometryType() == mfem::Geometry::SEGMENT)
    {
      segments.push_back(i);
    }
  }

  // Extract the curves and compute their bounding boxes along the way
  BoundingBox2D bbox;
  for(int i = 0; i < segments.size(); ++i)
  {
    auto curve = segment_to_curve(&mesh, i);
    SLIC_INFO_IF(verbose, axom::fmt::format("Element {}: {}", i, curve));

    bbox.addBox(curve.boundingBox());

    curves.emplace_back(std::move(curve));
  }

  SLIC_INFO(axom::fmt::format("Curve mesh bounding box: {}", bbox));

  // Early return if user didn't set up a query mesh
  if(boxMins.empty())
  {
    return 0;
  }

  // Generate a Cartesian (high order) mesh for the query points
  const auto query_res = axom::NumericArray<int, 2>(boxResolution.data());
  const auto query_box = BoundingBox2D(Point2D(boxMins.data()), Point2D(boxMaxs.data()));

  auto query_mesh = std::unique_ptr<mfem::Mesh>(
    axom::quest::util::make_cartesian_mfem_mesh_2D(query_box, query_res, queryOrder));
  auto fec = mfem::H1_FECollection(queryOrder, 2);
  auto fes = mfem::FiniteElementSpace(query_mesh.get(), &fec, 1);
  auto winding = mfem::GridFunction(&fes);
  auto inout = mfem::GridFunction(&fes);
  auto nodes_fes = query_mesh->GetNodalFESpace();

  // Query the winding numbers at each degree of freedom (DoF) of the query mesh.
  // The loop below independently checks (and adaptively refines) every curve for each query points.
  // A more efficient algorithm can de defined that caches the refined curves to avoid
  // extra refinements. We will add this in a follow-up PR.
  for(int nidx = 0; nidx < nodes_fes->GetNDofs(); ++nidx)
  {
    Point2D q;
    query_mesh->GetNode(nidx, q.data());

    double wn {};
    for(const auto& c : curves)
    {
      wn += axom::primal::winding_number(q, c);
    }

    winding[nidx] = wn;
    inout[nidx] = std::round(wn);
  }

  // Save the query mesh and fields to disk using a format that can be viewed in VisIt
  mfem::VisItDataCollection windingDC(outputPrefix, query_mesh.get());
  windingDC.RegisterField("winding", &winding);
  windingDC.RegisterField("inout", &inout);
  windingDC.Save();

  SLIC_INFO(axom::fmt::format("Outputting generated mesh '{}' to '{}'",
                              windingDC.GetCollectionName(),
                              axom::utilities::filesystem::getCWD()));

  return 0;
}
