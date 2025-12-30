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
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "mfem.hpp"

namespace primal = axom::primal;
using Point2D = primal::Point<double, 2>;
using NURBSCurve2D = primal::NURBSCurve<double, 2>;
using CurveGWNCache = primal::detail::NURBSCurveGWNCache<double>;
using BoundingBox2D = primal::BoundingBox<double, 2>;

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

  axom::Array<NURBSCurve2D> curves;
  axom::quest::MFEMReader mfem_reader;
  mfem_reader.setFileName(inputFile);

  int ret = mfem_reader.read(curves);
  if(ret != axom::quest::MFEMReader::READ_SUCCESS)
  {
    return 1;
  }

  // Extract the curves and compute their bounding boxes along the way
  BoundingBox2D bbox;
  int count {0};
  for(const auto& cur : curves)
  {
    SLIC_INFO_IF(verbose, axom::fmt::format("Element {}: {}", count++, cur));

    bbox.addBox(cur.boundingBox());

    // // Add curves to GWN Cache objects that dynamically store intermediate
    // //  curve subdivisions to be reused across query points
    // curves.emplace_back(CurveGWNCache(std::move(curve)));
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

  // Create grid functions for the winding number and inout fields
  auto fec = mfem::H1_FECollection(queryOrder, 2);
  auto fes = mfem::FiniteElementSpace(query_mesh.get(), &fec, 1);
  auto winding = mfem::GridFunction(&fes);
  auto inout = mfem::GridFunction(&fes);

  // Add utility function to convert an index into a query point
  const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();
  auto query_point = [&query_mesh](int idx) -> Point2D {
    Point2D pt;
    query_mesh->GetNode(idx, pt.data());
    return pt;
  };

  // Query the winding numbers at each degree of freedom (DoF) of the query mesh.
  // The loop below independently checks (and adaptively refines) every curve for each query point.
  for(int nidx = 0; nidx < num_query_points; ++nidx)
  {
    const Point2D q = query_point(nidx);
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
