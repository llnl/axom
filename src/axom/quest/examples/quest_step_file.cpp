// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
#include <fstream>
#include <map>
#include <set>

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

/**
 * Class that writes a patch-wise MFEM NURBS mesh containing a patch's trimming curves
 *
 * Each trimming curve is output as a separate 1D NURBS "patch" embedded in 2D space (u,v).
 *
 * Note: Support for reading these meshes was added to mfem after mfem@4.9.0
 * and is not in the current release of VisIt (visit@3.4.2)
 */
class PatchMFEMTrimmingCurveWriter
{
public:
  PatchMFEMTrimmingCurveWriter() { }

  void setOutputDirectory(const std::string& dir) { m_outputDirectory = dir; }
  void setVerbosity(bool verbosityFlag) { m_verbose = verbosityFlag; }
  void setNumFillZeros(int num)
  {
    if(num >= 0)
    {
      m_numFillZeros = num;
    }
  }

  void writeMFEMForPatch(int patchId,
                         const NURBSPatch& patch,
                         axom::ArrayView<const int> trimmingCurveWireIds) const
  {
    const auto& curves = patch.getTrimmingCurves();
    const int numCurves = curves.size();
    if(numCurves == 0)
    {
      return;
    }

    SLIC_WARNING_IF(trimmingCurveWireIds.size() != numCurves,
                    axom::fmt::format(
                      "Trimming curve wire id list size mismatch for patch {}: ids={}, curves={}. "
                      "Falling back to a single wire id.",
                      patchId,
                      trimmingCurveWireIds.size(),
                      numCurves));

    using axom::utilities::filesystem::joinPath;

    const std::string outDir = joinPath(m_outputDirectory, "mfem_trim_curves");
    if(!axom::utilities::filesystem::pathExists(outDir))
    {
      axom::utilities::filesystem::makeDirsForPath(outDir);
    }

    const std::string meshFilename =
      joinPath(outDir, axom::fmt::format("trim_curves_patch_{:0{}}.mesh", patchId, m_numFillZeros));

    std::ofstream meshFile(meshFilename);
    if(!meshFile.is_open())
    {
      SLIC_WARNING(axom::fmt::format("Unable to open '{}' for writing.", meshFilename));
      return;
    }

    axom::fmt::memory_buffer content;

    axom::fmt::format_to(std::back_inserter(content), "MFEM NURBS mesh v1.0\n\n");
    axom::fmt::format_to(std::back_inserter(content),
                         "# Trim curves for STEP NURBSPatch {}\n",
                         patchId);
    axom::fmt::format_to(std::back_inserter(content),
                         "# Parametric bbox: [{:.17g}, {:.17g}] x [{:.17g}, {:.17g}]\n",
                         patch.getMinKnot_u(),
                         patch.getMaxKnot_u(),
                         patch.getMinKnot_v(),
                         patch.getMaxKnot_v());
    axom::fmt::format_to(std::back_inserter(content), "# Number of trimming curves: {}\n\n", numCurves);

    axom::fmt::format_to(std::back_inserter(content), "dimension\n1\n\n");

    // One element per trimming curve; each element uses its own vertex pair.
    axom::fmt::format_to(std::back_inserter(content), "elements\n{}\n", numCurves);
    for(int i = 0; i < numCurves; ++i)
    {
      // MFEM attributes are 1-based; we store the STEP 0-based wire index + 1.
      const int wireId = (trimmingCurveWireIds.size() == numCurves) ? trimmingCurveWireIds[i] : 0;
      const int mfem_attribute = wireId + 1;
      const int v0 = 2 * i;
      const int v1 = 2 * i + 1;
      axom::fmt::format_to(std::back_inserter(content), "{} 1 {} {}\n", mfem_attribute, v0, v1);
    }

    axom::fmt::format_to(std::back_inserter(content), "\nboundary\n0\n\n");

    // Edge list provides the unique knotvector index and (optionally) encodes orientation.
    axom::fmt::format_to(std::back_inserter(content), "edges\n{}\n", numCurves);
    for(int i = 0; i < numCurves; ++i)
    {
      const int v0 = 2 * i;
      const int v1 = 2 * i + 1;
      axom::fmt::format_to(std::back_inserter(content), "{} {} {}\n", i, v0, v1);
    }

    axom::fmt::format_to(std::back_inserter(content), "\nvertices\n{}\n\n", 2 * numCurves);

    axom::fmt::format_to(std::back_inserter(content), "patches\n\n");
    for(int i = 0; i < numCurves; ++i)
    {
      const auto& curve = curves[i];
      SLIC_ASSERT(curve.isValidNURBS());
      const int wireId = (trimmingCurveWireIds.size() == numCurves) ? trimmingCurveWireIds[i] : 0;
      const int mfem_attribute = wireId + 1;

      const int degree = curve.getDegree();
      const int ncp = curve.getNumControlPoints();
      const auto& knots = curve.getKnots().getArray();
      SLIC_ASSERT(knots.size() == ncp + degree + 1);

      axom::fmt::format_to(std::back_inserter(content),
                           "# Curve wireId={} (mfem_attribute={}, output index {}): {} (degree {}, "
                           "{} control points)\n",
                           wireId,
                           mfem_attribute,
                           i,
                           curve.isRational() ? "rational" : "polynomial",
                           degree,
                           ncp);

      axom::fmt::format_to(std::back_inserter(content), "knotvectors\n1\n");
      axom::fmt::format_to(std::back_inserter(content), "{} {}", degree, ncp);
      for(const auto& kv : knots)
      {
        axom::fmt::format_to(std::back_inserter(content), "  {:.17g}", kv);
      }
      axom::fmt::format_to(std::back_inserter(content), "\n\n");

      axom::fmt::format_to(std::back_inserter(content), "dimension\n2\n\n");
      axom::fmt::format_to(std::back_inserter(content), "controlpoints_cartesian\n");

      const auto& cps = curve.getControlPoints();
      const auto& wts = curve.getWeights();
      for(int j = 0; j < ncp; ++j)
      {
        const double w = curve.isRational() ? wts[j] : 1.0;
        axom::fmt::format_to(std::back_inserter(content),
                             "{:.17g}  {:.17g}  {:.17g}\n",
                             cps[j][0],
                             cps[j][1],
                             w);
      }
      axom::fmt::format_to(std::back_inserter(content), "\n");
    }

    meshFile << axom::fmt::to_string(content);
    meshFile.close();

    SLIC_INFO_IF(m_verbose, axom::fmt::format("MFEM trim-curve mesh generated: '{}'", meshFilename));
  }

private:
  std::string m_outputDirectory;
  bool m_verbose {false};
  int m_numFillZeros {0};
};

/**
 * Class that writes a JSON stats file per patch summarizing trimming curves.
 *
 * The output includes:
 *  - number of wires (unique STEP wire indices)
 *  - number of trimming curves
 *  - min/max curve orders (order = degree + 1)
 *  - patch (u,v) knot-domain bounds
 */
class PatchTrimmingCurveStatsWriter
{
public:
  PatchTrimmingCurveStatsWriter() { }

  void setOutputDirectory(const std::string& dir) { m_outputDirectory = dir; }
  void setVerbosity(bool verbosityFlag) { m_verbose = verbosityFlag; }
  void setNumFillZeros(int num)
  {
    if(num >= 0)
    {
      m_numFillZeros = num;
    }
  }

  void writeStatsForPatch(int patchId,
                          const NURBSPatch& patch,
                          axom::ArrayView<const int> trimmingCurveWireIds,
                          bool was_originally_periodic_u,
                          bool was_originally_periodic_v) const
  {
    const auto& curves = patch.getTrimmingCurves();
    const int numCurves = curves.size();
    const bool is_trivially_trimmed = patch.isTriviallyTrimmed();

    const int num_knot_spans_u = static_cast<int>(patch.getKnots_u().getNumKnotSpans());
    const int num_knot_spans_v = static_cast<int>(patch.getKnots_v().getNumKnotSpans());

    int minOrder = 0;
    int maxOrder = 0;
    std::map<int, int> order_histogram;
    for(int i = 0; i < numCurves; ++i)
    {
      const int order = curves[i].getDegree() + 1;
      ++order_histogram[order];
      if(i == 0)
      {
        minOrder = order;
        maxOrder = order;
      }
      else
      {
        minOrder = std::min(minOrder, order);
        maxOrder = std::max(maxOrder, order);
      }
    }

    int numWires = 0;
    if(numCurves == 0)
    {
      numWires = 0;
    }
    else if(trimmingCurveWireIds.size() == numCurves)
    {
      std::set<int> uniqueWireIds;
      for(int i = 0; i < numCurves; ++i)
      {
        uniqueWireIds.insert(trimmingCurveWireIds[i]);
      }
      numWires = static_cast<int>(uniqueWireIds.size());
    }
    else
    {
      // If we can't trust the per-curve wire ids, conservatively report a single wire.
      numWires = 1;
    }

    using axom::utilities::filesystem::joinPath;

    const std::string outDir = joinPath(m_outputDirectory, "trim_curve_stats");
    if(!axom::utilities::filesystem::pathExists(outDir))
    {
      axom::utilities::filesystem::makeDirsForPath(outDir);
    }

    const std::string statsFilename =
      joinPath(outDir,
               axom::fmt::format("trim_curves_patch_{:0{}}.stats.json", patchId, m_numFillZeros));

    std::ofstream statsFile(statsFilename);
    if(!statsFile.is_open())
    {
      SLIC_WARNING(axom::fmt::format("Unable to open '{}' for writing.", statsFilename));
      return;
    }

    axom::fmt::memory_buffer content;
    axom::fmt::format_to(std::back_inserter(content), "{{\n");
    axom::fmt::format_to(std::back_inserter(content), "  \"patch_id\": {},\n", patchId);
    axom::fmt::format_to(std::back_inserter(content),
                         "  \"is_trivially_trimmed\": {},\n",
                         is_trivially_trimmed ? "true" : "false");
    axom::fmt::format_to(std::back_inserter(content),
                         "  \"was_originally_periodic_u\": {},\n",
                         was_originally_periodic_u ? "true" : "false");
    axom::fmt::format_to(std::back_inserter(content),
                         "  \"was_originally_periodic_v\": {},\n",
                         was_originally_periodic_v ? "true" : "false");
    axom::fmt::format_to(std::back_inserter(content),
                         "  \"num_knot_spans_u\": {},\n",
                         num_knot_spans_u);
    axom::fmt::format_to(std::back_inserter(content),
                         "  \"num_knot_spans_v\": {},\n",
                         num_knot_spans_v);
    axom::fmt::format_to(std::back_inserter(content), "  \"num_wires\": {},\n", numWires);
    axom::fmt::format_to(std::back_inserter(content), "  \"num_trimming_curves\": {},\n", numCurves);
    axom::fmt::format_to(std::back_inserter(content), "  \"min_curve_order\": {},\n", minOrder);
    axom::fmt::format_to(std::back_inserter(content), "  \"max_curve_order\": {},\n", maxOrder);
    axom::fmt::format_to(std::back_inserter(content), "  \"curves_by_order\": {{");
    {
      bool first = true;
      for(const auto& kv : order_histogram)
      {
        axom::fmt::format_to(std::back_inserter(content),
                             "{}\n    \"{}\": {}",
                             first ? "" : ",",
                             kv.first,
                             kv.second);
        first = false;
      }
      if(!order_histogram.empty())
      {
        axom::fmt::format_to(std::back_inserter(content), "\n  ");
      }
    }
    axom::fmt::format_to(std::back_inserter(content), "}},\n");
    axom::fmt::format_to(std::back_inserter(content),
                         "  \"uv_bbox\": {{\n"
                         "    \"min\": [{:.17g}, {:.17g}],\n"
                         "    \"max\": [{:.17g}, {:.17g}]\n"
                         "  }}\n",
                         patch.getMinKnot_u(),
                         patch.getMinKnot_v(),
                         patch.getMaxKnot_u(),
                         patch.getMaxKnot_v());
    axom::fmt::format_to(std::back_inserter(content), "}}\n");

    statsFile << axom::fmt::to_string(content);
    statsFile.close();

    SLIC_INFO_IF(m_verbose,
                 axom::fmt::format("Trim-curve stats JSON generated: '{}'", statsFilename));
  }

private:
  std::string m_outputDirectory;
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

  // Output options -----------------------------------------------------------
  auto* output_opts = app.add_option_group("Output", "Parameters associated with output");

  std::string output_dir = "step_output";
  output_opts->add_option("-o,--out,--output-dir", output_dir)
    ->description("Output directory for generated meshes")
    ->capture_default_str()
    ->check([](const std::string& dir) -> std::string {
      if(dir.find_first_of("\\:*?\"<>|") != std::string::npos)
      {
        return std::string("Output directory contains invalid characters.");
      }
      return std::string();
    });

  bool output_svg {false};
  output_opts->add_flag("--svg,--output-svg", output_svg)
    ->description("Generate SVG files for each NURBS patch")
    ->capture_default_str();

  bool output_mfem_trim_curves {false};
  output_opts->add_flag("--mfem-trim-curves,--output-mfem-trim-curves", output_mfem_trim_curves)
    ->description("Generate one MFEM NURBS mesh per trimmed patch containing its trimming curves")
    ->capture_default_str();

  bool output_trim_curve_stats_json {false};
  output_opts->add_flag("--stats-json,--output-trim-curve-stats-json", output_trim_curve_stats_json)
    ->description("Generate one JSON stats file per patch summarizing its trimming curves")
    ->capture_default_str();

  bool skip_trivial_trimmed_patches {false};
  output_opts
    ->add_flag("--skip-trivial,--skip-trivial-trimmed-patches", skip_trivial_trimmed_patches)
    ->description("Skip patch-wise SVG/MFEM outputs for trivially-trimmed patches")
    ->capture_default_str();

  // Triangulation options ----------------------------------------------------
  auto* tri_opts = app.add_option_group("Triangulation", "Parameters associated with triangulation");

  TriangleMeshOutputType output_trimmed {TriangleMeshOutputType::NONE};
  tri_opts->add_option("--tri.trimmed,--output-trimmed", output_trimmed)
    ->description("Output format for trimmed model triangulation: 'none', 'vtk', 'stl'")
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(validTriangleMeshOutputs));

  TriangleMeshOutputType output_untrimmed {TriangleMeshOutputType::NONE};
  tri_opts->add_option("--tri.untrimmed,--output-untrimmed", output_untrimmed)
    ->description("Output format for untrimmed model triangulation: 'none', 'vtk', 'stl'")
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(validTriangleMeshOutputs));

  double deflection {.1};
  tri_opts->add_option("--deflection", deflection)
    ->description("Max distance between actual geometry and triangulated geometry")
    ->capture_default_str();

  bool relative_deflection {false};
  tri_opts->add_flag("--relative", relative_deflection)
    ->description("Use relative deflection instead of absolute?")
    ->capture_default_str();

  double angular_deflection {0.5};
  tri_opts->add_option("--angular-deflection", angular_deflection)
    ->description("Angular deflection between adjacent normals when triangulating surfaces")
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

  // Ensure output directory exists iff we will write anything.
  const bool will_write_output = (output_trimmed != TriangleMeshOutputType::NONE) ||
    (output_untrimmed != TriangleMeshOutputType::NONE) || output_svg || output_mfem_trim_curves ||
    output_trim_curve_stats_json;
  if(is_root && will_write_output && !axom::utilities::filesystem::pathExists(output_dir))
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
    const int numFillZeros = (numPatches > 0) ? static_cast<int>(std::log10(numPatches)) + 1 : 1;

    PatchParametricSpaceProcessor patchProcessor;
    patchProcessor.setUnits(stepReader.getFileUnits());
    patchProcessor.setVerbosity(verbosity);
    patchProcessor.setOutputDirectory(output_dir);
    patchProcessor.setNumFillZeros(numFillZeros);

    for(int index = 0; index < numPatches; ++index)
    {
      if(skip_trivial_trimmed_patches && patches[index].isTriviallyTrimmed())
      {
        continue;
      }
      patchProcessor.generateSVGForPatch(index, patches[index]);
    }
  }

  //---------------------------------------------------------------------------
  // Optionally output an MFEM patch-wise NURBS mesh for each patch's trimming curves, only on root rank
  //---------------------------------------------------------------------------
  if(output_mfem_trim_curves && is_root)
  {
    const std::string outDir = joinPath(output_dir, "mfem_trim_curves");
    if(!axom::utilities::filesystem::pathExists(outDir))
    {
      axom::utilities::filesystem::makeDirsForPath(outDir);
    }

    SLIC_INFO(
      axom::fmt::format("Generating MFEM meshes for patch trimming curves in '{}' directory", outDir));

    const int numPatches = patches.size();
    const int numFillZeros = (numPatches > 0) ? static_cast<int>(std::log10(numPatches)) + 1 : 1;

    PatchMFEMTrimmingCurveWriter writer;
    writer.setVerbosity(verbosity);
    writer.setOutputDirectory(output_dir);
    writer.setNumFillZeros(numFillZeros);

    for(int index = 0; index < numPatches; ++index)
    {
      const int patch_id = stepReader.getPatchIds()[index];
      const auto wire_ids = stepReader.getTrimmingCurveWireIds(index);
      if(skip_trivial_trimmed_patches && patches[index].isTriviallyTrimmed())
      {
        continue;
      }

      writer.writeMFEMForPatch(patch_id, patches[index], wire_ids);
    }
  }

  //---------------------------------------------------------------------------
  // Optionally output trim-curve stats JSON per patch, only on root rank
  //---------------------------------------------------------------------------
  if(output_trim_curve_stats_json && is_root)
  {
    const std::string outDir = joinPath(output_dir, "trim_curve_stats");
    if(!axom::utilities::filesystem::pathExists(outDir))
    {
      axom::utilities::filesystem::makeDirsForPath(outDir);
    }

    SLIC_INFO(
      axom::fmt::format("Generating JSON trim-curve stats for patches in '{}' directory", outDir));

    const int numPatches = patches.size();
    const int numFillZeros = (numPatches > 0) ? static_cast<int>(std::log10(numPatches)) + 1 : 1;

    PatchTrimmingCurveStatsWriter stats_writer;
    stats_writer.setVerbosity(verbosity);
    stats_writer.setOutputDirectory(output_dir);
    stats_writer.setNumFillZeros(numFillZeros);

    for(int index = 0; index < numPatches; ++index)
    {
      const int patch_id = stepReader.getPatchIds()[index];
      const auto wire_ids = stepReader.getTrimmingCurveWireIds(index);
      stats_writer.writeStatsForPatch(patch_id,
                                      patches[index],
                                      wire_ids,
                                      stepReader.patchWasOriginallyPeriodic_u(index),
                                      stepReader.patchWasOriginallyPeriodic_v(index));
    }
  }

  return rc;
}
