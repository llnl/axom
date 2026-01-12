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

#include "axom/quest.hpp"

#include <iostream>

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

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

namespace slic = axom::slic;

namespace
{
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

  void generateSVGForPatch(int patchIndex, const NURBSPatch& patch)
  {
    using Point2D = axom::primal::Point<double, 2>;

    axom::primal::BoundingBox<double, 2> parametricBBox;
    parametricBBox.addPoint(Point2D {patch.getMinKnot_u(), patch.getMinKnot_v()});
    parametricBBox.addPoint(Point2D {patch.getMaxKnot_u(), patch.getMaxKnot_v()});

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

    const auto physicalBBox = patch.boundingBox();

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
                         "\n  <!-- Bounding box of ({},{})-degree patch in parametric space: {}; \n"
                         "       BBox in physical space: {} -->\n",
                         patch.getDegree_u(),
                         patch.getDegree_v(),
                         parametricBBox,
                         physicalBBox);

    axom::fmt::format_to(std::back_inserter(svgContent),
                         "  <rect x='{}' y='{}' width='{}' height='{}' />\n",
                         parametricBBox.getMin()[0],
                         parametricBBox.getMin()[1],
                         parametricBBox.range()[0],
                         parametricBBox.range()[1]);

    // add lines for the u- and v- knots

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

    axom::fmt::format_to(std::back_inserter(svgContent), "\n  <!-- Lines for u- knots -->\n");
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

    axom::fmt::format_to(std::back_inserter(svgContent), "\n  <!-- Lines for v- knots -->\n");
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
                         "\n  <!-- Paths for patch trimming curves -->\n");
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

#ifdef AXOM_USE_MPI

// utility function to help with MPI_Allreduce calls
template <typename T>
T allreduce_val(T localVal, MPI_Op op)
{
  T result = 0;
  MPI_Allreduce(&localVal, &result, 1, axom::mpi_traits<T>::type, op, MPI_COMM_WORLD);
  return result;
}

// utility function to help with MPI_Allreduce calls on booleans
bool allreduce_bool(bool localBool, MPI_Op op)
{
  const int localInt = localBool ? 1 : 0;
  int result = 0;
  MPI_Allreduce(&localInt, &result, 1, axom::mpi_traits<int>::type, op, MPI_COMM_WORLD);
  return result ? true : false;
}

// utility function to compare a value across ranks
template <typename T>
bool compare_across_ranks(T localVal, const std::string& check_str)
{
  const T minVal = allreduce_val(localVal, MPI_MIN);
  const T maxVal = allreduce_val(localVal, MPI_MAX);
  SLIC_WARNING_ROOT_IF(
    minVal != maxVal,
    axom::fmt::format("validation failed: {} is not consistent across ranks. min={}, max={}",
                      check_str,
                      minVal,
                      maxVal));

  return minVal == maxVal;
}

// utility function to compare a bool across ranks
bool compare_bool_across_ranks(bool localVal, const std::string& check_str)
{
  const int ival = localVal ? 1 : 0;
  const int all_true = allreduce_val(ival, MPI_LAND);
  const int all_false = !allreduce_val(ival, MPI_LOR);
  const bool consistent = all_true || all_false;

  SLIC_WARNING_ROOT_IF(
    !consistent,
    axom::fmt::format("patch validation failed: {} is not consistent across ranks.", check_str));

  return consistent;
}

/// Validates consistency of patch data across all ranks
bool validate_patches(const axom::Array<NURBSPatch>& patches)
{
  bool is_valid = true;

  // First, check that all ranks have the same number of patches
  axom::IndexType localNumPatches = patches.size();
  if(!compare_across_ranks(localNumPatches, "number of patches"))
  {
    return false;
  }
  // If no patches (and consistent), nothing more to check
  if(localNumPatches == 0)
  {
    return is_valid;
  }

  for(const auto& patch : patches)
  {
    // check that all patches are valid
    if(!allreduce_bool(patch.isValidNURBS(), MPI_LAND))
    {
      SLIC_WARNING("patch validation failed: patch is not valid on this rank");
      is_valid = false;
    }

    // check that orders and num control points and rationality are consistent
    if(!compare_across_ranks(patch.getOrder_u(), "order u"))
    {
      is_valid = false;
    }
    if(!compare_across_ranks(patch.getOrder_v(), "order v"))
    {
      is_valid = false;
    }
    if(!compare_across_ranks(patch.getNumControlPoints_u(), "num control points u"))
    {
      is_valid = false;
    }
    if(!compare_across_ranks(patch.getNumControlPoints_v(), "num control points v"))
    {
      is_valid = false;
    }

    if(!compare_bool_across_ranks(patch.isRational(), "patch rationality"))
    {
      is_valid = false;
    }

    // check trimming curve validity and consistency
    if(!compare_across_ranks(patch.getNumTrimmingCurves(), "num trimming curves"))
    {
      is_valid = false;
    }
    else
    {
      for(const auto& cur : patch.getTrimmingCurves())
      {
        if(!allreduce_bool(cur.isValidNURBS(), MPI_LAND))
        {
          SLIC_WARNING("patch validation failed: trimming curve is not valid on this rank");
          is_valid = false;
        }

        if(!compare_across_ranks(cur.getDegree(), "trimming curve degree"))
        {
          is_valid = false;
        }

        if(!compare_across_ranks(cur.getNumKnots(), "num trimming curve knots"))
        {
          is_valid = false;
        }

        if(!compare_across_ranks(cur.getNumControlPoints(), "num trimming curve control points"))
        {
          is_valid = false;
        }

        if(!compare_bool_across_ranks(cur.isRational(), "trimming curve rationality"))
        {
          is_valid = false;
        }
      }
    }
  }

  return is_valid;
}

/// Simple check that the triangle meshes are valid and consistent on all ranks
bool validate_triangle_mesh(const axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh)
{
  bool is_valid = true;

  // Check that all ranks have the same number of cells
  if(!compare_across_ranks(mesh.getNumberOfCells(), "number of cells"))
  {
    is_valid = false;
  }

  // Check that all ranks have the same number of nodes
  if(!compare_across_ranks(mesh.getNumberOfNodes(), "number of nodes"))
  {
    is_valid = false;
  }

  // Check consistency of presence of cell-centered "patch_index" field on all ranks
  const bool has_patch_index = mesh.hasField("patch_index", axom::mint::CELL_CENTERED);
  if(!compare_bool_across_ranks(has_patch_index, "presence of 'patch_index' cell field"))
  {
    is_valid = false;
  }
  if(is_valid && has_patch_index)
  {
    auto* patch_index_ptr = mesh.getFieldPtr<int>("patch_index", axom::mint::CELL_CENTERED);
    const bool is_patch_index_ptr_null = (patch_index_ptr == nullptr);
    if(!compare_bool_across_ranks(is_patch_index_ptr_null, "'patch_index' field nullptr"))
    {
      is_valid = false;
    }
  }

  return is_valid;
}
#endif

/// Sets up a parallel logger for this example, taking care of initialization and finalization
struct RAIILogger
{
  RAIILogger(int my_rank)
  {
    slic::initialize();
    slic::setIsRoot(my_rank == 0);

    slic::LogStream* logStream;

#ifdef AXOM_USE_MPI
    std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
    const int RLIMIT = 8;
    logStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
    logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
#else
    std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
    logStream = new slic::GenericOutputStream(&std::cout, fmt);
#endif  // AXOM_USE_MPI

    slic::addStreamToAllMsgLevels(logStream);
  }

  ~RAIILogger()
  {
    if(slic::isInitialized())
    {
      slic::flushStreams();
      slic::finalize();
    }
  }

  void setLoggingLevel(slic::message::Level level) { slic::setLoggingMsgLevel(level); }
};

enum class TriangleMeshOutputType
{
  NONE,
  VTK,
  STL
};

const std::map<std::string, TriangleMeshOutputType> validTriangleMeshOutputs {
  {"none", TriangleMeshOutputType::NONE},
  {"vtk", TriangleMeshOutputType::VTK},
  {"stl", TriangleMeshOutputType::STL}};

}  // namespace

int main(int argc, char** argv)
{
  constexpr static int RETURN_VALID = 0;
  constexpr static int RETURN_INVALID = 1;
  int rc = RETURN_VALID;

  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);

  const bool is_root = mpi_raii_wrapper.my_rank() == 0;
  RAIILogger raii_logger(mpi_raii_wrapper.my_rank());
  raii_logger.setLoggingLevel(slic::message::Info);

  //---------------------------------------------------------------------------
  // Set up and parse command line options
  //---------------------------------------------------------------------------
  axom::CLI::App app {"Quest Step File Example"};

  std::string filename;
  app.add_option("-f,--file", filename, "Input STEP file")->required();

  bool verbosity {false};
  app.add_flag("-v,--verbose", verbosity)->description("Enable verbose output")->capture_default_str();

  bool validate_model {true};
  app.add_flag("--validate", validate_model)
    ->description(
      axom::fmt::format("Validate the model while reading it in? (default: {})", validate_model))
    ->capture_default_str();

  std::string output_dir = "step_output";
  app.add_option("-o,--output-dir", output_dir)
    ->description("Output directory for generated meshes")
    ->capture_default_str()
    ->check([](const std::string& dir) -> std::string {
      if(dir.find_first_of("\\:*?\"<>|") != std::string::npos)
      {
        return std::string("Output directory contains invalid characters.");
      }
      return std::string();
    });

  TriangleMeshOutputType output_trimmed {TriangleMeshOutputType::VTK};
  app.add_option("--output-trimmed", output_trimmed)
    ->description("Output format for trimmed model triangulation: 'none', 'vtk', 'stl'")
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(validTriangleMeshOutputs));

  TriangleMeshOutputType output_untrimmed {TriangleMeshOutputType::NONE};
  app.add_option("--output-untrimmed", output_untrimmed)
    ->description("Output format for untrimmed model triangulation: 'none', 'vtk', 'stl'")
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(validTriangleMeshOutputs));

  double deflection {.1};
  app.add_option("--deflection", deflection)
    ->description("Max distance between actual geometry and triangulated geometry")
    ->capture_default_str();

  bool relative_deflection {false};
  app.add_flag("--relative", relative_deflection)
    ->description("Use relative deflection instead of absolute?")
    ->capture_default_str();

  double angular_deflection {0.5};
  app.add_option("--angular-deflection", angular_deflection)
    ->description("Angular deflection between adjacent normals when triangulating surfaces")
    ->capture_default_str();

  bool output_svg {false};
  app.add_flag("--output-svg", output_svg)
    ->description("Generate SVG files for each NURBS patch?")
    ->capture_default_str();

  app.get_formatter()->column_width(50);

  try
  {
    app.parse(argc, argv);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(is_root)
    {
      retval = app.exit(e);
    }

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    exit(retval);
#else
    return retval;
#endif
  }

  // Ensure output directory exists
  if(is_root && !axom::utilities::filesystem::pathExists(output_dir))
  {
    axom::utilities::filesystem::makeDirsForPath(output_dir);
  }

  //---------------------------------------------------------------------------
  // Load and process file
  //---------------------------------------------------------------------------
  SLIC_INFO_ROOT("Processing file: " << filename);
  SLIC_INFO_ROOT_IF(verbosity,
                    "Current working directory: " << axom::utilities::filesystem::getCWD());

  using NURBSPatch = axom::primal::NURBSPatch<double, 3>;
  using PatchArray = axom::Array<NURBSPatch>;

#ifdef AXOM_USE_MPI
  axom::quest::PSTEPReader stepReader(MPI_COMM_WORLD);
#else
  axom::quest::STEPReader stepReader;
#endif
  stepReader.setFileName(filename);
  stepReader.setVerbosity(verbosity);

  int res = stepReader.read(validate_model);
  if(res != 0)
  {
    SLIC_WARNING_ROOT("Error: The shape is invalid or empty.");
    return RETURN_INVALID;
  }

  if(is_root)
  {
    SLIC_INFO(axom::fmt::format("STEP file units: '{}'", stepReader.getFileUnits()));
    SLIC_INFO(stepReader.getBRepStats());
  }

  PatchArray& patches = stepReader.getPatchArray();

#ifdef AXOM_USE_MPI
  if(validate_model && !validate_patches(patches))
  {
    rc = RETURN_INVALID;
  }
#endif

  //---------------------------------------------------------------------------
  // Triangulate model
  //---------------------------------------------------------------------------
  using axom::utilities::filesystem::joinPath;

  // Create an unstructured triangle mesh of the model (using trimmed patches)
  if(output_trimmed != TriangleMeshOutputType::NONE)
  {
    SLIC_INFO_ROOT(
      axom::fmt::format("Generating triangulation of model in '{}' directory", output_dir));
    constexpr bool extract_trimmed_surface = true;

    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TRIANGLE);
    stepReader.getTriangleMesh(&mesh,
                               deflection,
                               angular_deflection,
                               relative_deflection,
                               extract_trimmed_surface);

#ifdef AXOM_USE_MPI
    if(validate_model && !validate_triangle_mesh(mesh))
    {
      rc = RETURN_INVALID;
    }
#endif

    if(is_root)
    {
      const std::string format = output_trimmed == TriangleMeshOutputType::VTK ? "vtk" : "stl";
      const std::string fname = axom::fmt::format("triangulated_mesh.{}", format);
      const std::string output_file = joinPath(output_dir, fname);

      if(output_trimmed == TriangleMeshOutputType::VTK)
      {
        axom::mint::write_vtk(&mesh, output_file);
      }
      else
      {
        axom::quest::STLWriter writer(output_file, true);
        writer.write(&mesh);
      }

      SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                  "Generated {} triangle mesh of trimmed model with {} deflection "
                                  "{} and {} angular deflection containing {:L} triangles: '{}'",
                                  format,
                                  false ? "relative" : "absolute",
                                  deflection,
                                  angular_deflection,
                                  mesh.getNumberOfCells(),
                                  output_file));
    }
  }

  // Create an unstructured triangle mesh of the model's untrimmed patches (mostly to understand the model better)
  if(output_untrimmed != TriangleMeshOutputType::NONE)
  {
    SLIC_INFO_ROOT(axom::fmt::format(
      "Generating triangulation of model (ignoring trimming curves) in '{}' directory",
      output_dir));

    constexpr bool extract_trimmed_surface = false;

    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TRIANGLE);
    stepReader.getTriangleMesh(&mesh,
                               deflection,
                               angular_deflection,
                               relative_deflection,
                               extract_trimmed_surface);

#ifdef AXOM_USE_MPI
    if(validate_model && !validate_triangle_mesh(mesh))
    {
      rc = RETURN_INVALID;
    }
#endif

    if(is_root)
    {
      const std::string format = output_untrimmed == TriangleMeshOutputType::VTK ? "vtk" : "stl";
      const std::string fname = axom::fmt::format("untrimmed_mesh.{}", format);
      const std::string output_file = joinPath(output_dir, fname);

      if(output_untrimmed == TriangleMeshOutputType::VTK)
      {
        axom::mint::write_vtk(&mesh, output_file);
      }
      else
      {
        axom::quest::STLWriter writer(output_file, true);
        writer.write(&mesh);
      }

      SLIC_INFO(
        axom::fmt::format(axom::utilities::locale(),
                          "Generated {} triangle mesh of untrimmed model with {} deflection {} "
                          "and {} angular deflection containing {:L} triangles: '{}'",
                          format,
                          false ? "relative" : "absolute",
                          deflection,
                          angular_deflection,
                          mesh.getNumberOfCells(),
                          output_file));
    }
  }

  //---------------------------------------------------------------------------
  // Optionally output an SVG for each patch, only on root rank
  //---------------------------------------------------------------------------
  if(output_svg && is_root)
  {
    SLIC_INFO(axom::fmt::format(
      "Generating SVG meshes for patches and their trimming curves in '{}' directory",
      output_dir));

    const int numPatches = patches.size();
    const int numFillZeros = static_cast<int>(std::log10(numPatches)) + 1;

    PatchParametricSpaceProcessor patchProcessor;
    patchProcessor.setUnits(stepReader.getFileUnits());
    patchProcessor.setVerbosity(verbosity);
    patchProcessor.setOutputDirectory(output_dir);
    patchProcessor.setNumFillZeros(numFillZeros);

    for(int index = 0; index < numPatches; ++index)
    {
      patchProcessor.generateSVGForPatch(index, patches[index]);
    }
  }

  return rc;
}
