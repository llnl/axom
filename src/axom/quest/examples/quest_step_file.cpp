// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "axom/quest/io/STEPReader.hpp"

#include "opencascade/BRepAdaptor_Curve.hxx"
#include "opencascade/BRepBuilderAPI_MakeFace.hxx"
#include "opencascade/BRepBuilderAPI_NurbsConvert.hxx"
#include "opencascade/BRepMesh_IncrementalMesh.hxx"
#include "opencascade/BRepTools.hxx"
#include "opencascade/BRep_Tool.hxx"
#include "opencascade/Geom_BSplineSurface.hxx"
#include "opencascade/Geom_RectangularTrimmedSurface.hxx"
#include "opencascade/Geom_Surface.hxx"
#include "opencascade/Poly_Triangulation.hxx"
#include "opencascade/Precision.hxx"
#include "opencascade/TopExp.hxx"
#include "opencascade/TopExp_Explorer.hxx"
#include "opencascade/TopLoc_Location.hxx"
#include "opencascade/TopoDS.hxx"
#include "opencascade/TopoDS_Face.hxx"
#include "opencascade/TopoDS_Shape.hxx"
#include <iostream>

/**
 * /file quest_step_file.cpp
 * /brief Example that loads in a STEP file and converts the surface patches and curves to Axom's NURBS representations
 *
 * This example reads in STEP files representing trimmed NURBS meshes using Open Cascade, 
 * converts the patches and trimming curves to Axom's NURBSPatch and NURBSCurve primitives, 
 * and generates various outputs including SVG and STL files.
 *
 * /note This example requires Axom to be configured with Open Cascade enabled.
 */

using PatchData = axom::quest::internal::PatchData;
using NURBSPatch = axom::quest::STEPReader::NURBSPatch;

/**
 * Class that generates SVG files over the parametric space of trimmed NURBS patches
 *
 * Uses a <rect> for the bounding box in parameter space; 
 * adds a <line> for each knot vector in u- and v-
 * and a <path> for each oriented trimming curve
 */
class PatchParametricSpaceProcessor
{
public:
  PatchParametricSpaceProcessor() { }

  void setOutputDirectory(const std::string& dir) { m_outputDirectory = dir; }
  void setUnits(const std::string& units) { m_units = units; }
  void setVerbosity(bool verbosityFlag) { m_verbose = verbosityFlag; }
  void setNumFillZeros(int num)
  {
    if(num >= 0)
    {
      m_numFillZeros = num;
    }
  }

  void generateSVGForPatch(int patchIndex, const PatchData& patchData, const NURBSPatch& patch)
  {
    const auto& parametricBBox = patchData.parametricBBox;

    SLIC_INFO_IF(m_verbose,
                 axom::fmt::format("Parametric BBox for patch {}: {}", patchIndex, parametricBBox));

    const auto& curves = patch.getTrimmingCurves();
    axom::fmt::memory_buffer svgContent;

    // Create a new bounding box by scaling and translating the parametricBBox
    auto scaledParametricBBox = parametricBBox;
    scaledParametricBBox.scale(1.25);

    SLIC_INFO_IF(m_verbose,
                 axom::fmt::format("Scaled and translated parametric BBox for patch {}: {}",
                                   patchIndex,
                                   scaledParametricBBox));

    // add the SVG header
    axom::fmt::format_to(std::back_inserter(svgContent),
                         "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' \n"
                         "     width='{0}{2}' height='{1}{2}' \n"
                         "     viewBox='{3} {4} {0} {1}' >\n",
                         scaledParametricBBox.range()[0],
                         scaledParametricBBox.range()[1],
                         m_units,
                         scaledParametricBBox.getMin()[0],
                         scaledParametricBBox.getMin()[1]);

    // add some CSS styles
    axom::fmt::format_to(std::back_inserter(svgContent), R"raw(
    <style>
      path {{ fill:none; stroke:black; stroke-width:.03; marker-end:url(#arrow); paint-order:fill stroke markers; stroke-linejoin:round; stroke-linecap:round; }}
      rect {{ fill: white; stroke: gray; stroke-width: 0.05; }}
      .u-line {{ fill: none; stroke: gray; stroke-width: 0.01; }}
      .v-line {{ fill: none; stroke: gray; stroke-width: 0.01; }}
    </style>
    )raw");

    // add a marker for the arrow's head to indicate the orientation
    axom::fmt::format_to(std::back_inserter(svgContent), R"raw(
    <defs>
      <marker id='arrow' style='overflow:visible' orient='auto-start-reverse'
          refX='0' refY='0'
          markerWidth='3.3239999' markerHeight='3.8427744'
          viewBox='0 0 5.3244081 6.1553851'>
        <path
            transform='scale(0.8)'
            style='fill:context-stroke;fill-rule:evenodd;stroke:none'
            d='M 5.77,0 L -2.88,4.5 L -1.44,0 L -2.88,-4.5 Z' />
      </marker>
    </defs>
    )raw");

    // add a rectangle for the parametric bounding box and a comment for its bounding boxes
    axom::fmt::format_to(std::back_inserter(svgContent),
                         "  <!-- Bounding box of ({},{})-degree patch in parametric space: {}; \n"
                         "       BBox in physical space: {} -->\n",
                         patch.getDegree_u(),
                         patch.getDegree_v(),
                         patchData.parametricBBox,
                         patchData.physicalBBox);

    axom::fmt::format_to(std::back_inserter(svgContent),
                         "  <rect x='{}' y='{}' width='{}' height='{}' />\n",
                         parametricBBox.getMin()[0],
                         parametricBBox.getMin()[1],
                         parametricBBox.range()[0],
                         parametricBBox.range()[1]);

    // add lines for the u- and v- knots
    axom::fmt::format_to(std::back_inserter(svgContent), "  <!-- Lines for u- and v- knots -->\n");

    auto unique_knots_and_multiplicities = [](const axom::Array<double>& knots_vector) {
      axom::Array<std::pair<double, int>> uniqueCounts;
      if(knots_vector.size() == 0)
      {
        return uniqueCounts;
      }

      double currentValue = knots_vector[0];
      int count = 1;

      for(int i = 1; i < knots_vector.size(); ++i)
      {
        if(knots_vector[i] == currentValue)
        {
          ++count;
        }
        else
        {
          uniqueCounts.emplace_back(currentValue, count);
          currentValue = knots_vector[i];
          count = 1;
        }
      }
      uniqueCounts.emplace_back(currentValue, count);

      return uniqueCounts;
    };

    for(const auto& u : unique_knots_and_multiplicities(patch.getKnots_u().getArray()))
    {
      axom::fmt::format_to(std::back_inserter(svgContent),
                           "  <line class='u-line mult-{}' x1='{}' y1='{}' x2='{}' y2='{}' />\n",
                           u.second,
                           u.first,
                           parametricBBox.getMin()[1],
                           u.first,
                           parametricBBox.getMax()[1]);
    }

    for(const auto& v : unique_knots_and_multiplicities(patch.getKnots_v().getArray()))
    {
      axom::fmt::format_to(std::back_inserter(svgContent),
                           "  <line class='v-line mult-{}' x1='{}' y1='{}' x2='{}' y2='{}' />\n",
                           v.second,
                           parametricBBox.getMin()[0],
                           v.first,
                           parametricBBox.getMax()[0],
                           v.first);
    }

    // add a path for each trimming curve
    // add lines for the u- and v- knots
    axom::fmt::format_to(std::back_inserter(svgContent),
                         "  <!-- Paths for patch trimming curves -->\n");
    for(const auto& curve : curves)
    {
      std::string pathData = nurbsCurveToSVGPath(curve);
      axom::fmt::format_to(std::back_inserter(svgContent), "{}\n", pathData);
    }

    // close the image and write to disk
    axom::fmt::format_to(std::back_inserter(svgContent), "</svg>");

    std::string svgFilename = axom::utilities::filesystem::joinPath(
      m_outputDirectory,
      axom::fmt::format("patch_{:0{}}.svg", patchIndex, m_numFillZeros));
    std::ofstream svgFile(svgFilename);
    if(svgFile.is_open())
    {
      svgFile << axom::fmt::to_string(svgContent);
      svgFile.close();
      SLIC_INFO_IF(m_verbose, "SVG file generated: " << svgFilename);
    }
    else
    {
      std::cerr << "Error: Unable to open file " << svgFilename << " for writing." << std::endl;
    }
  }

private:
  /**
   * Utility function to represent a NURBSCurve as an SVG path
   *
   * Since SVG only represents polynomial Bezier splines up to order 3,
   * this function discretizes rational curves and linear curves with order above three
   * to a polyline representation
   */
  std::string nurbsCurveToSVGPath(const axom::primal::NURBSCurve<double, 2>& curve)
  {
    using PointType = axom::primal::Point<double, 2>;

    const int degree = curve.getDegree();
    const auto& knotVector = curve.getKnots();
    const bool isRational = curve.isRational();

    axom::fmt::memory_buffer svgPath;
    axom::fmt::format_to(std::back_inserter(svgPath),
                         "  <path class='{} degree-{}' d='",
                         isRational ? "rational" : "non-rational",
                         degree);

    if(curve.isRational() || degree > 3)
    {
      const int numSamples = 100;
      const double tMin = knotVector[0];
      const double tMax = knotVector[knotVector.getNumKnots() - 1];

      for(int i = 0; i <= numSamples; ++i)
      {
        const double t = axom::utilities::lerp(tMin, tMax, static_cast<double>(i) / numSamples);

        PointType pt = curve.evaluate(t);
        if(i == 0)
        {
          axom::fmt::format_to(std::back_inserter(svgPath), "M {} {} ", pt[0], pt[1]);
        }
        else
        {
          axom::fmt::format_to(std::back_inserter(svgPath), "L {} {} ", pt[0], pt[1]);
        }
      }
    }
    else
    {
      auto bezierCurves = curve.extractBezier();
      for(const auto& bezier : bezierCurves)
      {
        const auto& bezierControlPoints = bezier.getControlPoints();
        if(degree == 2)
        {
          for(int i = 0; i < bezierControlPoints.size(); ++i)
          {
            const PointType& pt = bezierControlPoints[i];
            if(i == 0)
            {
              axom::fmt::format_to(std::back_inserter(svgPath), "M {} {} ", pt[0], pt[1]);
            }
            else if(i == 2)
            {
              axom::fmt::format_to(std::back_inserter(svgPath),
                                   "Q {} {} {} {} ",
                                   bezierControlPoints[1][0],
                                   bezierControlPoints[1][1],
                                   pt[0],
                                   pt[1]);
            }
          }
        }
        else if(degree == 3)
        {
          for(int i = 0; i < bezierControlPoints.size(); ++i)
          {
            const PointType& pt = bezierControlPoints[i];
            if(i == 0)
            {
              axom::fmt::format_to(std::back_inserter(svgPath), "M {} {} ", pt[0], pt[1]);
            }
            else if(i == 3)
            {
              axom::fmt::format_to(std::back_inserter(svgPath),
                                   "C {} {} {} {} {} {} ",
                                   bezierControlPoints[1][0],
                                   bezierControlPoints[1][1],
                                   bezierControlPoints[2][0],
                                   bezierControlPoints[2][1],
                                   pt[0],
                                   pt[1]);
            }
          }
        }
        else
        {
          for(int i = 0; i < bezierControlPoints.size(); ++i)
          {
            const PointType& pt = bezierControlPoints[i];
            if(i == 0)
            {
              axom::fmt::format_to(std::back_inserter(svgPath), "M {} {} ", pt[0], pt[1]);
            }
            else
            {
              axom::fmt::format_to(std::back_inserter(svgPath), "L {} {} ", pt[0], pt[1]);
            }
          }
        }
      }
    }

    // add the closing tags for the path
    axom::fmt::format_to(std::back_inserter(svgPath), "' />");

    return axom::fmt::to_string(svgPath);
  }

private:
  std::string m_outputDirectory {"."};
  std::string m_units {"mm"};
  bool m_verbose {false};
  int m_numFillZeros {0};
};

/**
 * Utility class to assist with triangulating STEP files
 * 
 * This class uses Open Cascade's triangulation functionality to generate triangle meshes
 * Supported triangulations:
 *  - triangulateTrimmedPatches: Generate a triangulation of the entire (trimmed) mesh and return as a mint mesh
 *  - triangulateUntrimmedPatches: Generate a triangulation of the entire (trimmed) mesh and return as a mint mesh
 *  The output mesh has a field 'patch_index' that associates each triangle with the index of its original patch in the input mesh
 */
class PatchTriangulator
{
public:
  PatchTriangulator() = delete;

  PatchTriangulator(const TopoDS_Shape& shape, double deflection, double angularDeflection)
    : m_shape(shape)
    , m_deflection(deflection)
    , m_angularDeflection(angularDeflection)
  {
    const bool is_relative = false;
    BRepTools::Clean(shape);
    BRepMesh_IncrementalMesh mesh(m_shape, m_deflection, is_relative, m_angularDeflection);

    if(!mesh.IsDone())
    {
      throw std::runtime_error("Mesh generation failed.");
    }
  }

  /// Utility function to triangulate each patch, ignoring the trimming curves, and write it to disk as as STL mesh
  void triangulateUntrimmedPatches(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& output_mesh)
  {
    std::vector<int> patch_id;

    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More(); faceExp.Next(), ++patchIndex)
    {
      TopoDS_Face face = TopoDS::Face(faceExp.Current());

      // Get the underlying surface of the face
      opencascade::handle<Geom_Surface> surface = BRep_Tool::Surface(face);

      // Optionally, you can create a rectangular trimmed surface if needed
      // Here, we assume the surface is already suitable for creating a new face
      Standard_Real u1, u2, v1, v2;
      surface->Bounds(u1, u2, v1, v2);
      opencascade::handle<Geom_RectangularTrimmedSurface> untrimmedSurface =
        new Geom_RectangularTrimmedSurface(surface, u1, u2, v1, v2);

      // Create a new face from the untrimmed surface
      TopoDS_Face newFace = BRepBuilderAPI_MakeFace(untrimmedSurface, Precision::Confusion());

      // Mesh the new face
      BRepMesh_IncrementalMesh mesh(newFace, m_deflection, Standard_False, m_angularDeflection);

      // Now you can access the triangulation of the new face
      TopLoc_Location loc;
      opencascade::handle<Poly_Triangulation> triangulation = BRep_Tool::Triangulation(newFace, loc);

      if(triangulation.IsNull())
      {
        SLIC_WARNING(
          axom::fmt::format("Error: Triangulation could not be generated for untrimmed patch {}",
                            patchIndex));
        break;
      }

      const int numTriangles = triangulation->NbTriangles();
      auto trsf = loc.Transformation();

      for(int i = 1; i <= numTriangles; ++i)
      {
        Poly_Triangle triangle = triangulation->Triangle(i);
        int n1, n2, n3;
        triangle.Get(n1, n2, n3);

        gp_Pnt p1 = triangulation->Node(n1).Transformed(trsf);
        gp_Pnt p2 = triangulation->Node(n2).Transformed(trsf);
        gp_Pnt p3 = triangulation->Node(n3).Transformed(trsf);

        axom::IndexType v1 = output_mesh.appendNode(p1.X(), p1.Y(), p1.Z());
        axom::IndexType v2 = output_mesh.appendNode(p2.X(), p2.Y(), p2.Z());
        axom::IndexType v3 = output_mesh.appendNode(p3.X(), p3.Y(), p3.Z());

        axom::IndexType cell[3] = {v1, v2, v3};
        output_mesh.appendCell(cell);
        patch_id.push_back(patchIndex);
      }
    }

    // Add a field to store the patch index for each cell
    auto* patchIndexField = output_mesh.createField<int>("patch_index", axom::mint::CELL_CENTERED);

    for(axom::IndexType i = 0; i < output_mesh.getNumberOfCells(); ++i)
    {
      patchIndexField[i] = patch_id[i];
    }
  }

  /// Triangulates the entire mesh as a single VTK file
  /// The association to the original patches is tracked via the patch_index field
  void triangulateTrimmedPatches(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& output_mesh)
  {
    std::vector<int> patch_id;

    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More(); faceExp.Next(), ++patchIndex)
    {
      TopoDS_Face face = TopoDS::Face(faceExp.Current());

      // Create a triangulation of this patch
      TopLoc_Location loc;
      opencascade::handle<Poly_Triangulation> triangulation = BRep_Tool::Triangulation(face, loc);

      if(triangulation.IsNull())
      {
        SLIC_WARNING(axom::fmt::format("Error: Triangulation could not be generated for patch {}",
                                       patchIndex));
        continue;
      }

      const int numTriangles = triangulation->NbTriangles();
      auto trsf = loc.Transformation();

      for(int i = 1; i <= numTriangles; ++i)
      {
        Poly_Triangle triangle = triangulation->Triangle(i);
        int n1, n2, n3;
        triangle.Get(n1, n2, n3);

        gp_Pnt p1 = triangulation->Node(n1).Transformed(trsf);
        gp_Pnt p2 = triangulation->Node(n2).Transformed(trsf);
        gp_Pnt p3 = triangulation->Node(n3).Transformed(trsf);

        axom::IndexType v1 = output_mesh.appendNode(p1.X(), p1.Y(), p1.Z());
        axom::IndexType v2 = output_mesh.appendNode(p2.X(), p2.Y(), p2.Z());
        axom::IndexType v3 = output_mesh.appendNode(p3.X(), p3.Y(), p3.Z());

        axom::IndexType cell[3] = {v1, v2, v3};
        output_mesh.appendCell(cell);
        patch_id.push_back(patchIndex);
      }
    }

    // Add a field to store the patch index for each cell
    auto* patchIndexField = output_mesh.createField<int>("patch_index", axom::mint::CELL_CENTERED);

    for(axom::IndexType i = 0; i < output_mesh.getNumberOfCells(); ++i)
    {
      patchIndexField[i] = patch_id[i];
    }
  }

private:
  TopoDS_Shape m_shape;
  double m_deflection;
  double m_angularDeflection;
};

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  // Set up and parse command line options
  axom::CLI::App app {"Quest Step File Example"};

  std::string filename;
  app.add_option("-f,--file", filename, "Input file")->required();

  bool verbosity {false};
  app.add_flag("-v,--verbose", verbosity)->description("Enable verbose output")->capture_default_str();

  double deflection {.1};
  app.add_option("--deflection", deflection)
    ->description("Max distance between actual geometry and triangulated geometry")
    ->capture_default_str();

  double angular_deflection {0.5};
  app.add_option("--angular_deflection", angular_deflection)
    ->description("Angular deflection between adjacent normals when triangulating surfaces")
    ->capture_default_str();

  std::string output_dir = "step_output";
  app.add_option("-o,--output_dir", output_dir)
    ->description("Output directory for generated meshes")
    ->capture_default_str()
    ->check([](const std::string& dir) -> std::string {
      if(dir.find_first_of("\\:*?\"<>|") != std::string::npos)
      {
        return std::string("Output directory contains invalid characters.");
      }
      return std::string();
    });

  CLI11_PARSE(app, argc, argv);

  // Ensure output directory exists
  if(!axom::utilities::filesystem::pathExists(output_dir))
  {
    axom::utilities::filesystem::makeDirsForPath(output_dir);
  }

  // Load and process file
  SLIC_INFO("Processing file: " << filename);
  SLIC_INFO_IF(verbosity, "Current working directory: " << axom::utilities::filesystem::getCWD());

  using NURBSPatch = axom::primal::NURBSPatch<double, 3>;
  using PatchArray = axom::Array<NURBSPatch>;

  axom::quest::STEPReader stepReader;
  stepReader.setFileName(filename);
  stepReader.setVerbosity(verbosity);

  int res = stepReader.read();

  if(res != 0)
  {
    std::cerr << "Error: The shape is invalid or empty." << std::endl;
    return 1;
  }

  ////---------------------------
  PatchArray& patches = stepReader.getPatchArray();
  const int numPatches = patches.size();
  const int numFillZeros = static_cast<int>(std::log10(numPatches)) + 1;

  // Generate outputs
  PatchParametricSpaceProcessor patchProcessor;
  patchProcessor.setUnits(stepReader.getFileUnits());
  patchProcessor.setVerbosity(verbosity);
  patchProcessor.setOutputDirectory(output_dir);
  patchProcessor.setNumFillZeros(numFillZeros);
  SLIC_INFO(
    axom::fmt::format("Generating SVG meshes for patches and their trimming "
                      "curves in '{}' directory",
                      output_dir));
  for(const auto& [index, patchData] : stepReader.getPatchDataMap())
  {
    patchProcessor.generateSVGForPatch(index, patchData, patches[index]);
  }

  auto& nurbs_shape = stepReader.getShape();
  PatchTriangulator patchTriangulator(nurbs_shape, deflection, angular_deflection);
  SLIC_INFO(
    axom::fmt::format("Generating triangles meshes for trimmed and untrimmed "
                      "patches in '{}' directory",
                      output_dir));
  {
    // Create an unstructured mesh with 3D vertices and triangular cells
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TRIANGLE);

    patchTriangulator.triangulateTrimmedPatches(mesh);

    const std::string filename =
      axom::utilities::filesystem::joinPath(output_dir, "triangulated_mesh.vtk");
    axom::mint::write_vtk(&mesh, filename);
    SLIC_INFO_IF(verbosity, "VTK triangle mesh of entire model generated: " << filename);
  }
  {
    // Create an unstructured mesh with 3D vertices and triangular cells
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TRIANGLE);
    patchTriangulator.triangulateUntrimmedPatches(mesh);

    const std::string filename =
      axom::utilities::filesystem::joinPath(output_dir, "untrimmed_mesh.vtk");
    axom::mint::write_vtk(&mesh, filename);
    SLIC_INFO_IF(verbosity, "VTK triangle mesh of entire (untrimmed) model generated: " << filename);
  }

  return 0;
}
