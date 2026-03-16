// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_winding_number2d.cpp
 * \brief Example that computes the winding number of a grid of points
 * against a collection of 2D parametric rational curves.
 * Supports MFEM meshes in the cubic positive Bernstein basis or the (rational)
 * NURBS basis.
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/quest.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace primal = axom::primal;
namespace quest = axom::quest;

using Point2D = primal::Point<double, 2>;
using BoundingBox2D = primal::BoundingBox<double, 2>;

using NURBSCurve2D = primal::NURBSCurve<double, 2>;

namespace
{
/**
 * This helper class takes in an mfem mesh (potentially with variable order curves from mfem>4.9)
 * and writes out a version that is compatible with mfem@4.9 
 * In particular, this allows us to visualize it with current versions of VisIt which do not yet support this feature.
 * 
 * We can remove this class once downstream appications (such as VisIt) are updated to a version of mfem
 * that support the NURBS patches format.
 */
class MFEM49ElevatedNURBSMeshWriter
{
public:
  explicit MFEM49ElevatedNURBSMeshWriter(double tol = 1e-12) : m_tol(tol) { }

  bool writeElevatedMesh(const std::string& input_file, const std::string& output_file) const
  {
    mfem::Mesh mesh(input_file, /*generate_edges=*/1, /*refine=*/1);

    if(mesh.NURBSext == nullptr)
    {
      SLIC_WARNING(
        axom::fmt::format("Input mesh '{}' has no NURBS extension; skipping degree elevation",
                          input_file));
      return false;
    }

    // Note: NURBS curve meshes in mfem@4.9 must all have the same degree, and their knotvectors must have the same length
    const mfem::Array<int>& orders = mesh.NURBSext->GetOrders();
    int max_order = 0;
    for(int i = 0; i < orders.Size(); ++i)
    {
      max_order = std::max(max_order, orders[i]);
    }

    if(max_order <= 0)
    {
      SLIC_WARNING(
        axom::fmt::format("Input mesh '{}' has invalid NURBS orders; skipping degree elevation",
                          input_file));
      return false;
    }

    mesh.DegreeElevate(max_order, max_order);

    makeKnotVectorsUniform(mesh);

    // Write the modified mesh to output_file
    if(!writeMeshPrintFormat(mesh, output_file))
    {
      SLIC_WARNING(axom::fmt::format("Failed to write elevated mesh '{}'", output_file));
      return false;
    }

    SLIC_INFO(axom::fmt::format("Wrote elevated MFEM 4.9-compatible NURBS mesh '{}' (max order {})",
                                output_file,
                                max_order));
    return true;
  }

private:
  struct KnotBin
  {
    std::int64_t key {0};
    double value {0.0};
    int max_multiplicity {0};
  };

  // MFEM 4.9's 1D NURBS reader assumes all knotvectors have the same GetNE/GetNCP
  // and uses KnotVec(0) when computing offsets. To generate MFEM 4.9/VisIt
  // compatible files, insert knots so every knotvector matches the maximum
  // knot multiplicity pattern present in the mesh.
  void makeKnotVectorsUniform(mfem::Mesh& mesh) const
  {
    if(mesh.NURBSext == nullptr)
    {
      return;
    }
    if(mesh.Dimension() != 1 || mesh.NURBSext->Dimension() != 1)
    {
      return;
    }

    const int nkv = mesh.NURBSext->GetNKV();
    if(nkv <= 1)
    {
      return;
    }

    const double inv_tol = (m_tol > 0.0) ? (1.0 / m_tol) : 1.0e12;
    auto key_for = [inv_tol](double t) -> std::int64_t {
      return static_cast<std::int64_t>(std::llround(t * inv_tol));
    };

    std::vector<std::unordered_map<std::int64_t, int>> counts(nkv);
    std::unordered_map<std::int64_t, KnotBin> bins;

    for(int i = 0; i < nkv; ++i)
    {
      const mfem::KnotVector* kv = mesh.NURBSext->GetKnotVector(i);
      if(kv == nullptr)
      {
        continue;
      }

      auto& local = counts[i];
      for(int j = 0; j < kv->Size(); ++j)
      {
        const double value = (*kv)[j];
        const std::int64_t key = key_for(value);

        const int multiplicity = ++local[key];
        auto& bin = bins[key];
        bin.key = key;
        bin.value = value;
        bin.max_multiplicity = std::max(bin.max_multiplicity, multiplicity);
      }
    }

    std::vector<KnotBin> sorted_bins;
    sorted_bins.reserve(bins.size());
    for(const auto& it : bins)
    {
      sorted_bins.push_back(it.second);
    }
    std::sort(sorted_bins.begin(), sorted_bins.end(), [](const KnotBin& a, const KnotBin& b) {
      if(a.value != b.value)
      {
        return a.value < b.value;
      }
      return a.key < b.key;
    });

    int total_to_insert = 0;
    mfem::Array<mfem::Vector*> insertions(nkv);
    for(int i = 0; i < nkv; ++i)
    {
      int need_total = 0;
      for(const auto& bin : sorted_bins)
      {
        const auto found = counts[i].find(bin.key);
        const int have = (found == counts[i].end()) ? 0 : found->second;
        need_total += std::max(0, bin.max_multiplicity - have);
      }

      insertions[i] = new mfem::Vector(need_total);
      int pos = 0;
      for(const auto& bin : sorted_bins)
      {
        const auto found = counts[i].find(bin.key);
        const int have = (found == counts[i].end()) ? 0 : found->second;
        const int need = std::max(0, bin.max_multiplicity - have);
        for(int k = 0; k < need; ++k)
        {
          (*insertions[i])[pos++] = bin.value;
        }
      }
      total_to_insert += need_total;
    }

    if(total_to_insert > 0)
    {
      mesh.KnotInsert(insertions);
    }

    for(int i = 0; i < nkv; ++i)
    {
      delete insertions[i];
    }
  }

  bool writeMeshPrintFormat(const mfem::Mesh& mesh, const std::string& output_file) const
  {
    std::ostringstream oss;
    oss.precision(16);
    mesh.Print(oss);

    std::istringstream iss(oss.str());
    std::ofstream ofs(output_file);
    if(!ofs.is_open())
    {
      return false;
    }
    ofs.precision(16);

    // Note: mfem@4.9's mesh reader does not support the "NURBS2" finite element collection,
    // if it occurs, we replace it with the (essentially equivalent) "NURBS"
    std::string line;
    while(std::getline(iss, line))
    {
      constexpr const char* prefix = "FiniteElementCollection: NURBS";
      if(line.rfind(prefix, 0) == 0)
      {
        ofs << prefix << '\n';
      }
      else
      {
        ofs << line << '\n';
      }
    }

    return true;
  }

private:
  double m_tol {1e-12};
};

}  // namespace

//------------------------------------------------------------------------------
// CLI input
//------------------------------------------------------------------------------

class Input
{
public:
  std::string inputFile;
  std::string outputPrefix = {"winding2d"};

  bool verbose {false};
  std::string annotationMode {"none"};
  bool memoized {true};
  bool vis {true};
  bool stats {false};
  std::string elevatedMeshFile;

  const std::array<std::string, 2> valid_algorithms {"direct", "fast-approximation"};
  std::string algorithm {valid_algorithms[1]};  // fast-approximation

  bool linearize {false};
  int approximation_order {2};

  bool useUniformLinearization;
  int segmentsPerKnotSpan {10};
  double percentError {1.0};

  // Query mesh parameters
  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int queryOrder {1};

  primal::WindingTolerances tol;

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
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
    app.add_flag("--memoized,!--no-memoized", memoized, "Cache geometric data during query?")
      ->capture_default_str();
    app.add_flag("--vis,!--no-vis", vis, "Should we write out the results for visualization?")
      ->capture_default_str();
    app.add_flag("--stats,!--no-stats", stats, "Compute summary stats for query fields?")
      ->capture_default_str();

    app.add_option("--output-elevated-mesh", elevatedMeshFile)
      ->description(
        "Optional. Output MFEM mesh after elevating all NURBS curve orders to the maximum order of "
        "input")
      ->capture_default_str();

    // Options for query tolerances
    app.add_option("--edge-tol", tol.edge_tol)
      ->description("Relative edge tolerance for queries")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();
    app.add_option("--eps-tol", tol.EPS)
      ->description("Additional generic tolerance parameter")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();

#ifdef AXOM_USE_CALIPER
    app.add_option("--caliper", annotationMode)
      ->description(
        "caliper annotation mode. Valid options include 'none' and 'report'. "
        "Use 'help' to see full list.")
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    // Options for triangulation of the input STEP file
    auto* linearize_curves_subcommand =
      app.add_subcommand("linearize_curves")
        ->description("Options for linearizing NURBS curves. Default is ")
        ->fallthrough();

    auto* nsegments =
      linearize_curves_subcommand->add_option("--num-segments", segmentsPerKnotSpan)
        ->description(
          "Number of segments for each knot span of each input curve for a uniform linearization.")
        ->check(axom::CLI::PositiveNumber)
        ->capture_default_str();
    auto* perror =
      linearize_curves_subcommand->add_option("--percent-error", percentError)
        ->description(
          "The percent of error that is acceptable to stop refinement during non-uniform "
          "linearization.")
        ->check(axom::CLI::Range(0.0f, 100.0f))
        ->capture_default_str();
    linearize_curves_subcommand->add_option("--algorithm", algorithm)
      ->description(
        "Use direct evaluation instead of fast, heirarchical approximation? (significantly "
        "slower, slightly more precise)")
      ->capture_default_str()
      ->check(axom::CLI::IsMember(valid_algorithms));
    linearize_curves_subcommand
      ->add_option("--approximation-order",
                   approximation_order,
                   "The order of the Taylor expansion (lower is faster, less precise)")
      ->expected(0, 2)
      ->capture_default_str();

    auto* query_mesh_subcommand =
      app.add_subcommand("query_mesh")->description("Options for setting up a query mesh")->fallthrough();
    auto* minbb = query_mesh_subcommand->add_option("--min", boxMins)
                    ->description("Min bounds for box mesh (x,y)")
                    ->expected(2);
    auto* maxbb = query_mesh_subcommand->add_option("--max", boxMaxs)
                    ->description("Max bounds for box mesh (x,y)")
                    ->expected(2);
    query_mesh_subcommand->add_option("--res", boxResolution)
      ->description("Resolution of the box mesh (i,j)")
      ->expected(2)
      ->required();
    query_mesh_subcommand->add_option("--order", queryOrder)
      ->description("polynomial order of the query mesh")
      ->check(axom::CLI::PositiveNumber);

    // add some requirements -- if user provides minbb or maxbb, we need both
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    // If a user provides a number of refinements, we can't have a percent error
    nsegments->excludes(perror);
    useUniformLinearization = nsegments->count() > 0;

    app.parse(argc, argv);

    linearize = app.got_subcommand("linearize_curves");
  }
};

using GWNQueryType = std::variant<axom::quest::DirectGWN2D,
                                  axom::quest::PolylineGWN2D<0>,
                                  axom::quest::PolylineGWN2D<1>,
                                  axom::quest::PolylineGWN2D<2>>;

GWNQueryType make_gwn_query(bool linearize_curves, int approximation_order)
{
  if(linearize_curves)
  {
    if(approximation_order == 0)
    {
      return axom::quest::PolylineGWN2D<0> {};
    }
    else if(approximation_order == 1)
    {
      return axom::quest::PolylineGWN2D<1> {};
    }
    else
    {
      return axom::quest::PolylineGWN2D<2> {};
    }
  }

  return axom::quest::DirectGWN2D {};
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger raii_logger;

  // Parse command line arguments into input
  Input input;
  axom::CLI::App app {
    "Load mesh containing collection of curves"
    " and optionally generate a query mesh of winding numbers."};

  try
  {
    input.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  axom::utilities::raii::AnnotationsWrapper annotation_raii_wrapper(input.annotationMode);
  AXOM_ANNOTATE_SCOPE("winding number example");

  if(!input.elevatedMeshFile.empty())
  {
    AXOM_ANNOTATE_SCOPE("write_elevated_mesh");

    MFEM49ElevatedNURBSMeshWriter writer;
    writer.writeElevatedMesh(input.inputFile, input.elevatedMeshFile);
  }

  // Read curves from the MFEM mesh
  axom::Array<NURBSCurve2D> curves;
  {
    AXOM_ANNOTATE_SCOPE("read_mesh");

    axom::quest::MFEMReader mfem_reader;
    mfem_reader.setFileName(input.inputFile);

    const int ret = mfem_reader.read(curves);
    if(ret != axom::quest::MFEMReader::READ_SUCCESS)
    {
      SLIC_ERROR("Failed to read MFEM file.");
      return 1;
    }
  }

  // Linearize the input curves if asked for
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> poly_mesh(2, axom::mint::SEGMENT);
  if(input.linearize)
  {
    AXOM_ANNOTATE_SCOPE("linearization");

    axom::utilities::Timer timer(true);
    axom::quest::LinearizeCurves lc;
    if(input.useUniformLinearization)
    {
      lc.getLinearMeshUniform(curves.view(), &poly_mesh, input.segmentsPerKnotSpan);
    }
    else
    {
      lc.getLinearMeshNonUniform(curves.view(), &poly_mesh, input.percentError);
    }
    timer.stop();

    SLIC_INFO(axom::fmt::format(
      axom::utilities::locale(),
      "Discretized {} curves with {} segments in each knot span for a total of {:L} segments.",
      curves.size(),
      input.segmentsPerKnotSpan,
      poly_mesh.getNumberOfCells()));
    SLIC_INFO(
      axom::fmt::format("Preprocessing stage (linearization): {} s", timer.elapsedTimeInSec()));
  }

  // Early return if user didn't set up a query mesh
  if(input.boxResolution.empty())
  {
    return 0;
  }

  // Extract the curves and compute their bounding boxes along the way
  BoundingBox2D shape_bbox;
  for(const auto& cur : curves)
  {
    shape_bbox.addBox(cur.boundingBox());
  }
  SLIC_INFO(axom::fmt::format("Curve mesh bounding box: {}", shape_bbox));

  // Query grid setup;
  // if user did not provide a bounding box, user input bounding box scaled by 10%
  mfem::DataCollection dc("winding_query");
  {
    // Create the desired winding number query instance
    auto wn_query = make_gwn_query(app.got_subcommand("linearize_curves"), input.approximation_order);

    // Generate the query grid and fields
    quest::generate_gwn_query_mesh(dc,
                                   shape_bbox,
                                   input.boxMins,
                                   input.boxMaxs,
                                   input.boxResolution,
                                   input.queryOrder);

    // Run the preprocess
    std::visit(
      [&](auto& wn) {
        using T = std::decay_t<decltype(wn)>;
        if constexpr(quest::gwn_input_type_v<T> == quest::GWNInputType::Curve)
        {
          wn.preprocess(curves, input.memoized);
        }
        else if constexpr(quest::gwn_input_type_v<T> == quest::GWNInputType::Polyline)
        {
          wn.preprocess(&poly_mesh, input.algorithm == "direct");
        }
      },
      wn_query);

    // Run the query
    std::visit([&](auto& wn) { wn.query(dc, input.tol); }, wn_query);
  }

  // Postprocess query results: norms, ranges, and integral statistics
  if(input.stats)
  {
    AXOM_ANNOTATE_SCOPE("postprocess");

    auto& winding = *dc.GetField("winding");
    auto& inout = *dc.GetField("inout");

    const auto winding_stats = axom::quest::compute_field_stats(winding);
    const auto inout_stats = axom::quest::compute_field_stats(inout);
    const auto inout_integrals = axom::quest::compute_integrals(inout);

    std::int64_t pos_inout_dofs = 0;
    std::int64_t neg_inout_dofs = 0;
    for(int i = 0; i < inout.Size(); ++i)
    {
      const double v = inout[i];
      if(v > 0.0)
      {
        ++pos_inout_dofs;
      }
      else if(v < 0.0)
      {
        ++neg_inout_dofs;
      }
    }

    SLIC_INFO(
      axom::fmt::format("WN_STATS: dof_l2={:.6e} dof_linf={:.6e} l2={:.6e} min={:.6e} max={:.6e}",
                        winding_stats.dof_l2,
                        winding_stats.dof_linf,
                        winding_stats.l2,
                        winding_stats.min,
                        winding_stats.max));

    SLIC_INFO(axom::fmt::format(
      "INOUT_STATS: dof_l2={:.6e} dof_linf={:.6e} l2={:.6e} min={:.6e} max={:.6e} volume={:.6e} "
      "domain_volume={:.6e} vol_frac={:.6e} pos_dofs={} neg_dofs={}",
      inout_stats.dof_l2,
      inout_stats.dof_linf,
      inout_stats.l2,
      inout_stats.min,
      inout_stats.max,
      inout_integrals.integral,
      inout_integrals.domain_volume,
      (inout_integrals.domain_volume > 0.0 ? inout_integrals.integral / inout_integrals.domain_volume
                                           : 0.0),
      pos_inout_dofs,
      neg_inout_dofs));
  }

  // Save the query mesh and fields to disk using a format that can be viewed in VisIt
  if(input.vis)
  {
    AXOM_ANNOTATE_SCOPE("dump_mesh");

    mfem::VisItDataCollection windingDC(input.outputPrefix, dc.GetMesh());
    windingDC.RegisterField("winding", dc.GetField("winding"));
    windingDC.RegisterField("inout", dc.GetField("inout"));
    windingDC.Save();

    SLIC_INFO(axom::fmt::format("Outputting generated mesh '{}' to '{}'",
                                windingDC.GetCollectionName(),
                                axom::utilities::filesystem::getCWD()));
  }

  return 0;
}
