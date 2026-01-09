// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/STEPReader.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#ifndef AXOM_USE_OPENCASCADE
  #error STEPReader should only be included when Axom is configured with opencascade
#endif

#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include "opencascade/BRep_Tool.hxx"
#include "opencascade/BRepAdaptor_Curve.hxx"
#include "opencascade/BRepBuilderAPI_MakeFace.hxx"
#include "opencascade/BRepBuilderAPI_NurbsConvert.hxx"
#include "opencascade/BRepCheck_Analyzer.hxx"
#include "opencascade/BRepCheck_Edge.hxx"
#include "opencascade/BRepCheck_Face.hxx"
#include "opencascade/BRepCheck_Wire.hxx"
#include "opencascade/BRepCheck.hxx"
#include "opencascade/BRepLib.hxx"
#include "opencascade/BRepMesh_IncrementalMesh.hxx"
#include "opencascade/BRepTools.hxx"
#include "opencascade/Geom_BSplineSurface.hxx"
#include "opencascade/Geom_RectangularTrimmedSurface.hxx"
#include "opencascade/Geom_Surface.hxx"
#include "opencascade/Geom2d_BSplineCurve.hxx"
#include "opencascade/Geom2d_Curve.hxx"
#include "opencascade/Geom2dConvert.hxx"
#include "opencascade/GeomAbs_CurveType.hxx"
#include "opencascade/Interface_Static.hxx"
#include "opencascade/Poly_Triangulation.hxx"
#include "opencascade/Precision.hxx"
#include "opencascade/ShapeBuild_ReShape.hxx"
#include "opencascade/ShapeFix_Edge.hxx"
#include "opencascade/ShapeFix_Face.hxx"
#include "opencascade/ShapeFix_Wire.hxx"
#include "opencascade/STEPControl_Reader.hxx"
#include "opencascade/TColgp_Array2OfPnt.hxx"
#include "opencascade/TopExp_Explorer.hxx"
#include "opencascade/TopExp.hxx"
#include "opencascade/TopLoc_Location.hxx"
#include "opencascade/TopoDS_Edge.hxx"
#include "opencascade/TopoDS_Face.hxx"
#include "opencascade/TopoDS_Shape.hxx"
#include "opencascade/TopoDS_Wire.hxx"
#include "opencascade/TopoDS.hxx"

#include <iostream>

namespace axom
{
namespace quest
{
namespace internal
{
/// Struct to hold data associated with each surface patch of the mesh
struct PatchData
{
  int patchIndex {-1};
  bool wasOriginallyPeriodic_u {false};
  bool wasOriginallyPeriodic_v {false};
  axom::primal::BoundingBox<double, 2> parametricBBox;
  axom::primal::BoundingBox<double, 3> physicalBBox;
  axom::Array<bool> trimmingCurves_originallyPeriodic;
};

using PatchDataMap = std::map<int, PatchData>;

/**
 * Class to read in a STEP file representing trimmed NURBS meshes using Open Cascade 
 * and convert the patches and trimming curves to Axom's NURBSPatch and NURBSCurve primitives.
 * 
 * Implementation note: Since Axom's primitives do not support periodic knots, 
 * we must convert the Open Cascade analogues to an open/clamped representation, when necessary.
 */
class StepFileProcessor
{
public:
  enum class LoadStatus
  {
    UNINITIALIZED = 0,
    SUCCESS = 1 << 0,
    FAILED_TO_READ = 1 << 1,
    FAILED_NO_SHAPES = 1 << 2,
    FAILED_TO_CONVERT = 1 << 3,
    FAILED = FAILED_TO_READ | FAILED_TO_CONVERT
  };

  static constexpr int CurveDim = 2;
  static constexpr int SpaceDim = 3;
  using NCurve = axom::primal::NURBSCurve<double, CurveDim>;
  using NCurveArray = axom::Array<NCurve>;

  using NPatch = axom::primal::NURBSPatch<double, SpaceDim>;
  using NPatchArray = axom::Array<NPatch>;

  using PointType = axom::primal::Point<double, CurveDim>;
  using VectorType = axom::primal::Vector<double, CurveDim>;
  using PointType2D = axom::primal::Point<double, CurveDim>;
  using PointType3D = axom::primal::Point<double, SpaceDim>;
  using BBox2D = axom::primal::BoundingBox<double, CurveDim>;
  using BBox3D = axom::primal::BoundingBox<double, SpaceDim>;
  using PatchToTrimmingCurvesMap = std::map<int, NCurveArray>;

private:
  /// Returns a bounding box convering the patch's knot spans in 2D parametric space
  BBox2D faceBoundingBox(const TopoDS_Face& face) const
  {
    BBox2D bbox;

    opencascade::handle<Geom_Surface> surface = BRep_Tool::Surface(face);

    Standard_Real u1, u2, v1, v2;
    surface->Bounds(u1, u2, v1, v2);
    bbox.addPoint(PointType {u1, v1});
    bbox.addPoint(PointType {u2, v2});

    return bbox;
  }

  /// Helper class to convert the face geometry of CAD mesh to valid NURBSPatch instances
  /// The constructor converts the surface to a clamped (non-periodic) representation, if necessary
  class PatchProcessor
  {
  public:
    PatchProcessor() = delete;

    PatchProcessor(const opencascade::handle<Geom_BSplineSurface>& surface, bool verbose = false)
      : m_surface(surface)
      , m_verbose(verbose)
    {
      m_inputSurfaceWasPeriodic_u = m_surface->IsUPeriodic();
      m_inputSurfaceWasPeriodic_v = m_surface->IsVPeriodic();

      ensureClamped();
    }

    const opencascade::handle<Geom_BSplineSurface>& getSurface() const { return m_surface; }

    /// Returns a representation of the surface geometry as an axom::primal::NURBSPatch
    NPatch nurbsPatchGeometry() const
    {
      // Check if the surface is periodic in u or v
      const bool isUPeriodic = m_surface->IsUPeriodic();
      const bool isVPeriodic = m_surface->IsVPeriodic();
      SLIC_ERROR_IF(isUPeriodic || isVPeriodic,
                    "Axom's NURBSPatch only supports non-periodic patches");

      // Extract weights, if the surface is rational
      const bool isRational = m_surface->IsURational() || m_surface->IsVRational();

      // Create the NURBSPatch from control points, weights, and knots
      return isRational
        ? NPatch(extractControlPoints(),
                 extractWeights(),
                 extractCombinedKnots_u(),
                 extractCombinedKnots_v())
        : NPatch(extractControlPoints(), extractCombinedKnots_u(), extractCombinedKnots_v());
    }

    bool patchWasOriginallyPeriodic_u() const { return m_inputSurfaceWasPeriodic_u; }
    bool patchWasOriginallyPeriodic_v() const { return m_inputSurfaceWasPeriodic_v; }

    /**
     * Utility function to compare the represented surface to \a otherSurface
     * at \a numSamples uniformly sampled points in parameter space
     * Returns true when the sum of distances is less than the provided tolerance
     */
    bool compareToSurface(opencascade::handle<Geom_BSplineSurface> otherSurface,
                          int numSamples,
                          double sq_tol = 1e-8) const
    {
      if(numSamples <= 1)
      {
        return true;
      }

      auto knot_vals_u = extractKnotValues_u();
      axom::Array<double> params_u(numSamples);
      axom::numerics::linspace(knot_vals_u.front(), knot_vals_u.back(), params_u.data(), numSamples);

      auto knot_vals_v = extractKnotValues_v();
      axom::Array<double> params_v(numSamples);
      axom::numerics::linspace(knot_vals_v.front(), knot_vals_v.back(), params_v.data(), numSamples);

      auto evaluateSurface = [](auto surface, double u, double v) {
        gp_Pnt point;
        surface->D0(u, v, point);
        return PointType3D {point.X(), point.Y(), point.Z()};
      };

      double squared_sum = 0.;
      for(auto u : params_u)
      {
        for(auto v : params_v)
        {
          auto my_val = evaluateSurface(m_surface, u, v);
          auto other_val = evaluateSurface(otherSurface, u, v);

          const double sq_dist = axom::primal::squared_distance(my_val, other_val);
          squared_sum += sq_dist;
          SLIC_WARNING_IF(m_verbose && sq_dist > sq_tol,
                          axom::fmt::format("Distance between surfaces at evaluated param "
                                            "({},{}) exceeded tolerance {}.\n"
                                            "Point on my surface {}; Point on other surface "
                                            "{}; Squared distance {}  (running sum {})",
                                            u,
                                            v,
                                            sq_tol,
                                            my_val,
                                            other_val,
                                            sq_dist,
                                            squared_sum));
        }
      }

      return squared_sum <= sq_tol;
    }

    /// Logs some information about the patch
    void printSurfaceInfo() const
    {
      SLIC_INFO(axom::fmt::format("Patch is periodic in u: {}", m_surface->IsUPeriodic()));
      SLIC_INFO(axom::fmt::format("Patch is periodic in v: {}", m_surface->IsVPeriodic()));

      // Print control points
      {
        const auto patch_control_points = extractControlPoints();
        SLIC_INFO(axom::fmt::format("Patch control points ({} x {}): {}",
                                    patch_control_points.shape()[0],
                                    patch_control_points.shape()[1],
                                    patch_control_points));
      }

      // Print weights (if the surface is rational)
      if(const bool isRational = m_surface->IsURational() || m_surface->IsVRational(); isRational)
      {
        const auto patch_weights = extractWeights();
        SLIC_INFO(axom::fmt::format("Patch weights ({} x {}): {}",
                                    patch_weights.shape()[0],
                                    patch_weights.shape()[1],
                                    patch_weights));
      }
      else
      {
        SLIC_INFO("Patch is polynomial (uniform weights)");
      }
    }

  private:
    /// converts the surface from periodic knots to clamped knots, when necessary
    void ensureClamped()
    {
      if(!m_surface->IsUPeriodic() && !m_surface->IsVPeriodic())
      {
        return;  // nothing to do, return
      }

      // compute the period
      const Standard_Integer iu1 = m_surface->FirstUKnotIndex();
      const Standard_Integer iu2 = m_surface->LastUKnotIndex();
      const Standard_Integer iv1 = m_surface->FirstVKnotIndex();
      const Standard_Integer iv2 = m_surface->LastVKnotIndex();

      const Standard_Real Tu = m_surface->UKnot(iu2) - m_surface->UKnot(iu1);
      const Standard_Real Tv = m_surface->VKnot(iv2) - m_surface->VKnot(iv1);

      // Choose a start parameter and wrap to the patch's period
      Standard_Real U0 = m_surface->UKnot(iu1);
      Standard_Real V0 = m_surface->VKnot(iv1);
      m_surface->PeriodicNormalization(U0, V0);

      // Trim to one full period (updates control nets / knot vectors)
      m_surface->Segment(U0, U0 + Tu, V0, V0 + Tv);

      m_surface->SetUNotPeriodic();
      m_surface->SetVNotPeriodic();
    }

  public:
    /// extracts control points (poles) from the patch as a 2D axom::Array
    axom::Array<PointType3D, 2> extractControlPoints() const
    {
      axom::Array<PointType3D, 2> patch_control_points;

      TColgp_Array2OfPnt poles(1, m_surface->NbUPoles(), 1, m_surface->NbVPoles());
      m_surface->Poles(poles);

      patch_control_points.resize(axom::ArrayOptions::Uninitialized {},
                                  poles.ColLength(),
                                  poles.RowLength());
      for(int i = poles.LowerRow(); i <= poles.UpperRow(); ++i)
      {
        for(int j = poles.LowerCol(); j <= poles.UpperCol(); ++j)
        {
          gp_Pnt pole = poles(i, j);
          patch_control_points(i - 1, j - 1) = PointType3D {pole.X(), pole.Y(), pole.Z()};
        }
      }

      return patch_control_points;
    }

    /// extracts weights from the patch as a 2D axom::Array,
    /// each weight corresponds to a control point
    axom::Array<double, 2> extractWeights() const
    {
      axom::Array<double, 2> patch_weights;

      TColStd_Array2OfReal weights(1, m_surface->NbUPoles(), 1, m_surface->NbVPoles());
      m_surface->Weights(weights);

      patch_weights.resize(axom::ArrayOptions::Uninitialized {},
                           weights.ColLength(),
                           weights.RowLength());
      for(int i = weights.LowerRow(); i <= weights.UpperRow(); ++i)
      {
        for(int j = weights.LowerCol(); j <= weights.UpperCol(); ++j)
        {
          patch_weights(i - 1, j - 1) = weights(i, j);
        }
      }

      return patch_weights;
    }

    /// extracts the u knot vector from the patch, accounting for multiplicities
    axom::Array<double> extractCombinedKnots_u() const
    {
      const auto vals = extractKnotValues_u();
      const auto mults = extractKnotMultiplicities_u();
      const int total_knots = std::accumulate(mults.begin(), mults.end(), 0);

      axom::Array<double> knots_u(0, total_knots);
      for(int i = 0; i < vals.size(); ++i)
      {
        knots_u.insert(knots_u.end(), mults[i], vals[i]);
      }
      SLIC_ASSERT(knots_u.size() == total_knots);

      return knots_u;
    }

    /// extracts the v knot vector from the patch, accounting for multiplicities
    axom::Array<double> extractCombinedKnots_v() const
    {
      const auto vals = extractKnotValues_v();
      const auto mults = extractKnotMultiplicities_v();
      const int total_knots = std::accumulate(mults.begin(), mults.end(), 0);

      axom::Array<double> knots_v(0, total_knots);
      for(int i = 0; i < vals.size(); ++i)
      {
        knots_v.insert(knots_v.end(), mults[i], vals[i]);
      }
      SLIC_ASSERT(knots_v.size() == total_knots);

      return knots_v;
    }

    /// converts u knot values to axom Array w/o accounting for multiplicity
    /// /sa extractKnotMultiplicities_u
    axom::Array<double> extractKnotValues_u() const
    {
      const int num_knots = m_surface->NbUKnots();
      axom::Array<double> uKnotsArray(0, num_knots);

      TColStd_Array1OfReal uKnots(1, num_knots);
      m_surface->UKnots(uKnots);

      for(int i = uKnots.Lower(); i <= uKnots.Upper(); ++i)
      {
        uKnotsArray.push_back(uKnots(i));
      }

      return uKnotsArray;
    }

    /// converts u knot multiplicities to axom Array
    /// /sa extractKnotValues_u
    axom::Array<int> extractKnotMultiplicities_u() const
    {
      const int num_knots = m_surface->NbUKnots();
      axom::Array<int> uMultsArray(0, num_knots);

      TColStd_Array1OfInteger uMults(1, num_knots);
      m_surface->UMultiplicities(uMults);

      for(int i = uMults.Lower(); i <= uMults.Upper(); ++i)
      {
        uMultsArray.push_back(uMults(i));
      }

      return uMultsArray;
    }

    /// converts v knot values to axom Array w/o accounting for multiplicity
    /// /sa extractKnotMultiplicities_v
    axom::Array<double> extractKnotValues_v() const
    {
      const int num_knots = m_surface->NbVKnots();
      axom::Array<double> vKnotsArray(0, num_knots);

      TColStd_Array1OfReal vKnots(1, num_knots);
      m_surface->VKnots(vKnots);

      for(int i = vKnots.Lower(); i <= vKnots.Upper(); ++i)
      {
        vKnotsArray.push_back(vKnots(i));
      }

      return vKnotsArray;
    }

    /// converts v knot multiplicities to axom Array
    /// /sa extractKnotValues_v
    axom::Array<int> extractKnotMultiplicities_v() const
    {
      const int num_knots = m_surface->NbVKnots();
      axom::Array<int> vMultsArray(0, num_knots);

      TColStd_Array1OfInteger vMults(1, num_knots);
      m_surface->VMultiplicities(vMults);

      for(int i = vMults.Lower(); i <= vMults.Upper(); ++i)
      {
        vMultsArray.push_back(vMults(i));
      }

      return vMultsArray;
    }

  private:
    opencascade::handle<Geom_BSplineSurface> m_surface;

    bool m_verbose {false};
    bool m_inputSurfaceWasPeriodic_u {false};
    bool m_inputSurfaceWasPeriodic_v {false};
  };

  /// Helper class to convert trimming curves to valid NURBSCurve instances
  /// The constructor converts the curve to a clamped (non-periodic) representation, if necessary
  class CurveProcessor
  {
  public:
    CurveProcessor() = delete;

    CurveProcessor(const opencascade::handle<Geom2d_BSplineCurve>& curve,
                   bool verbose = false,
                   bool ensure_clamped = true)
      : m_curve(curve)
      , m_verbose(verbose)
    {
      m_inputSurfaceWasPeriodic = m_curve->IsPeriodic();

      if(ensure_clamped)
      {
        ensureClamped();
      }
    }

    /// Returns a representation of the trimming curve as an axom::primal::NURBSCurve
    NCurve nurbsCurve() const
    {
      const bool isPeriodic = m_curve->IsPeriodic();
      SLIC_ERROR_IF(isPeriodic, "Axom's NURBSCurve only supports non-periodic curves");

      return m_curve->IsRational()
        ? NCurve(extractControlPoints(), extractWeights(), extractCombinedKnots())
        : NCurve(extractControlPoints(), extractCombinedKnots());
    }

    bool curveWasOriginallyPeriodic() const { return m_inputSurfaceWasPeriodic; }

    /**
     * Utility function to compare the represented curve to \a otherCurve
     * at \a numSamples uniformly sampled points in parameter space
     * Returns true when the sum of distances is less than the provided tolerance
     */
    bool compareToCurve(opencascade::handle<Geom2d_BSplineCurve>& otherCurve,
                        int numSamples,
                        double sq_tol = 1e-8) const
    {
      if(numSamples <= 1)
      {
        return true;
      }

      auto knot_vals = extractKnotValues();

      axom::Array<double> params(numSamples);
      axom::numerics::linspace(knot_vals.front(), knot_vals.back(), params.data(), numSamples);

      auto evaluateCurve = [](auto curve, double t) {
        const gp_Pnt2d knot_point = curve->Value(t);
        return PointType {knot_point.X(), knot_point.Y()};
      };

      double squared_sum = 0.;
      for(auto val : params)
      {
        auto my_val = evaluateCurve(m_curve, val);
        auto other_val = evaluateCurve(otherCurve, val);

        const double sq_dist = axom::primal::squared_distance(my_val, other_val);
        squared_sum += sq_dist;
        SLIC_WARNING_IF(m_verbose && sq_dist > sq_tol,
                        axom::fmt::format("Distance between curves at evaluated param {} "
                                          "exceeded tolerance {}.\n"
                                          "Point on my curve {}; Point on other curve {}; "
                                          "Squared distance {}  (running sum {})",
                                          val,
                                          sq_tol,
                                          my_val,
                                          other_val,
                                          sq_dist,
                                          squared_sum));
      }

      return squared_sum <= sq_tol;
    }

  private:
    /// converts the curve from periodic knots to clamped knots, when necessary
    void ensureClamped()
    {
      if(!m_curve->IsPeriodic())
      {
        return;
      }

      // Compute the period
      const Standard_Integer i1 = m_curve->FirstUKnotIndex();
      const Standard_Integer i2 = m_curve->LastUKnotIndex();
      const Standard_Real T = m_curve->Knot(i2) - m_curve->Knot(i1);

      // Choose a start parameter, and normalize to curve's period
      Standard_Real U0 = m_curve->FirstParameter();
      m_curve->PeriodicNormalization(U0);  // no-op if not periodic

      // Trim to one full period (updates control points and knot vectors)
      m_curve->Segment(U0, U0 + T);
    }

  public:
    axom::Array<PointType> extractControlPoints() const
    {
      axom::Array<PointType> controlPoints;

      TColgp_Array1OfPnt2d paraPoints(1, m_curve->NbPoles());
      m_curve->Poles(paraPoints);

      for(Standard_Integer i = paraPoints.Lower(); i <= paraPoints.Upper(); ++i)
      {
        gp_Pnt2d paraPt = paraPoints(i);
        controlPoints.emplace_back(PointType2D {paraPt.X(), paraPt.Y()});
      }

      return controlPoints;
    }

    axom::Array<double> extractWeights() const
    {
      axom::Array<double> weights;
      if(m_curve->IsRational())
      {
        TColStd_Array1OfReal curveWeights(1, m_curve->NbPoles());
        m_curve->Weights(curveWeights);
        weights.reserve(curveWeights.Length());
        for(int i = 1; i <= curveWeights.Length(); ++i)
        {
          weights.push_back(curveWeights(i));
        }
      }
      return weights;
    }

    axom::Array<double> extractCombinedKnots() const
    {
      const auto vals = extractKnotValues();
      const auto mults = extractKnotMultiplicities();
      const int total_knots = std::accumulate(mults.begin(), mults.end(), 0);

      axom::Array<double> knots(0, total_knots);
      for(int i = 0; i < vals.size(); ++i)
      {
        knots.insert(knots.end(), mults[i], vals[i]);
      }
      SLIC_ASSERT(knots.size() == total_knots);

      return knots;
    }

    /// converts knot values to axom Array w/o accounting for multiplicity
    /// /sa extractKnotMultiplicities
    axom::Array<double> extractKnotValues() const
    {
      const int num_knots = m_curve->NbKnots();
      axom::Array<double> knots(0, num_knots);

      TColStd_Array1OfReal occ_knots(1, num_knots);
      m_curve->Knots(occ_knots);
      for(int i = 1; i <= occ_knots.Length(); ++i)
      {
        knots.push_back(occ_knots(i));
      }

      return knots;
    }

    /// converts knot multiplicities to axom Array
    /// /sa extractKnotValues
    axom::Array<int> extractKnotMultiplicities() const
    {
      const int num_knots = m_curve->NbKnots();
      axom::Array<int> mults(0, num_knots);

      TColStd_Array1OfInteger occ_mults(1, num_knots);
      m_curve->Multiplicities(occ_mults);
      for(int i = 1; i <= occ_mults.Length(); ++i)
      {
        mults.push_back(occ_mults(i));
      }

      return mults;
    }

  private:
    opencascade::handle<Geom2d_BSplineCurve> m_curve;
    bool m_verbose {false};
    bool m_inputSurfaceWasPeriodic {false};
  };

public:
  StepFileProcessor() = delete;

  StepFileProcessor(const std::string& filename, bool verbose = false) : m_verbose(verbose)
  {
    m_shape = loadStepFile(filename);
  }

  void setVerbosity(bool verbose) { m_verbose = verbose; }

  bool isLoaded() const { return m_loadStatus == LoadStatus::SUCCESS; }

  const TopoDS_Shape& getShape() const { return m_shape; }

  int numPatchesInFile() const
  {
    int count = 0;
    for(TopExp_Explorer ex(m_shape, TopAbs_FACE); ex.More(); ex.Next())
    {
      ++count;
    }
    return count;
  }

  /// Extracts data from the faces of the mesh and converts each patch to a NURBSPatch
  void extractPatches(NPatchArray& patches)
  {
    patches.resize(numPatchesInFile());

    int patchIndex = 0;
    for(TopExp_Explorer ex(m_shape, TopAbs_FACE); ex.More(); ex.Next(), ++patchIndex)
    {
      SLIC_DEBUG_IF(m_verbose, "*** Processing patch " << patchIndex);
      const TopoDS_Face& face = TopoDS::Face(ex.Current());

      opencascade::handle<Geom_Surface> surface = BRep_Tool::Surface(face);
      if(surface->IsKind(STANDARD_TYPE(Geom_BSplineSurface)))
      {
        opencascade::handle<Geom_BSplineSurface> bsplineSurface =
          opencascade::handle<Geom_BSplineSurface>::DownCast(surface);

        PatchProcessor patchProcessor(bsplineSurface, m_verbose);

        patches[patchIndex] = patchProcessor.nurbsPatchGeometry();

        // If the face is flipped in opencascade, we need to flip the primal primitive too
        if(face.Orientation() == TopAbs_REVERSED)
        {
          patches[patchIndex].reverseOrientation_u();
        }

        PatchData& patchData = m_patchData[patchIndex];
        patchData.patchIndex = patchIndex;
        patchData.wasOriginallyPeriodic_u = patchProcessor.patchWasOriginallyPeriodic_u();
        patchData.wasOriginallyPeriodic_v = patchProcessor.patchWasOriginallyPeriodic_v();
        patchData.parametricBBox = faceBoundingBox(face);
        patchData.physicalBBox = patches[patchIndex].boundingBox();

        if(patchData.wasOriginallyPeriodic_u || patchData.wasOriginallyPeriodic_v)
        {
          opencascade::handle<Geom_BSplineSurface> origSurface =
            opencascade::handle<Geom_BSplineSurface>::DownCast(surface);

          const bool withinThreshold = patchProcessor.compareToSurface(origSurface, 25);
          SLIC_WARNING_IF(!withinThreshold,
                          axom::fmt::format("[Patch {}] Patch geometry was not "
                                            "within threshold after clamping.",
                                            patchIndex,
                                            patches[patchIndex]));
        }
      }
      else
      {
        const std::string surfaceType = surface->DynamicType()->Name();
        SLIC_WARNING(fmt::format("Skipping patch {} with non-BSpline surface type: '{}'",
                                 patchIndex,
                                 surfaceType));
      }
    }
  }

  void validateBRep()
  {
    // lambda to return a string with all the (error) status problems in an entity
    auto printStatusList = [](const std::string& prefix, const BRepCheck_ListOfStatus& statusList) {
      auto statusToString = [](BRepCheck_Status st) -> std::string {
        // clang-format off
          switch(st)
          {
          case BRepCheck_NoError: return "NoError";
          case BRepCheck_InvalidPointOnCurve: return "InvalidPointOnCurve";
          case BRepCheck_InvalidPointOnCurveOnSurface: return "InvalidPointOnCurveOnSurface";
          case BRepCheck_InvalidPointOnSurface: return "InvalidPointOnSurface";
          case BRepCheck_No3DCurve: return "No3DCurve";
          case BRepCheck_Multiple3DCurve: return "Multiple3DCurve";
          case BRepCheck_Invalid3DCurve: return "Invalid3DCurve";
          case BRepCheck_NoCurveOnSurface: return "NoCurveOnSurface";
          case BRepCheck_InvalidCurveOnSurface: return "InvalidCurveOnSurface";
          case BRepCheck_InvalidCurveOnClosedSurface: return "InvalidCurveOnClosedSurface";
          case BRepCheck_InvalidSameRangeFlag: return "InvalidSameRangeFlag";
          case BRepCheck_InvalidSameParameterFlag: return "InvalidSameParameterFlag";
          case BRepCheck_InvalidDegeneratedFlag: return "InvalidDegeneratedFlag";
          case BRepCheck_FreeEdge: return "FreeEdge";
          case BRepCheck_InvalidMultiConnexity: return "InvalidMultiConnexity";
          case BRepCheck_InvalidRange: return "InvalidRange";
          case BRepCheck_EmptyWire: return "EmptyWire";
          case BRepCheck_RedundantEdge: return "RedundantEdge";
          case BRepCheck_SelfIntersectingWire: return "SelfIntersectingWire";
          case BRepCheck_NoSurface: return "NoFace";
          case BRepCheck_InvalidWire: return "InvalidWire";
          case BRepCheck_RedundantWire: return "RedundantWire";
          case BRepCheck_IntersectingWires: return "IntersectingWires";
          case BRepCheck_InvalidImbricationOfWires: return "InvalidImbricationOfWires";
          case BRepCheck_EmptyShell: return "EmptyShell";
          case BRepCheck_RedundantFace: return "RedundantFace";
          case BRepCheck_InvalidImbricationOfShells: return "InvalidImbricationOfShells";
          case BRepCheck_UnorientableShape: return "UnorientableShape";
          case BRepCheck_NotClosed: return "NotClosed";
          case BRepCheck_NotConnected: return "NotConnected";
          case BRepCheck_SubshapeNotInShape: return "SubshapeNotInShape";
          case BRepCheck_BadOrientation: return "BadOrientation";
          case BRepCheck_BadOrientationOfSubshape: return "BadOrientationOfSubshape";
          case BRepCheck_InvalidPolygonOnTriangulation: return "InvalidPolygonOnTriangulation";
          case BRepCheck_InvalidToleranceValue: return "InvalidToleranceValue";
          case BRepCheck_EnclosedRegion: return "EnclosedRegion";
          case BRepCheck_CheckFail: return "CheckFail";
          default: return "UnknownStatus";
          }
        // clang-format on
      };

      std::stringstream sstr;
      for(BRepCheck_ListIteratorOfListOfStatus it(statusList); it.More(); it.Next())
      {
        if(const BRepCheck_Status st = it.Value(); st != BRepCheck_NoError)
        {
          sstr << axom::fmt::format("  - {} ({})\n", static_cast<int>(st), statusToString(st));
        }
      }
      if(!sstr.str().empty())
      {
        SLIC_INFO(prefix << "\n" << sstr.str());
      }
    };

    // Diagnose each face
    int patchIndex = 0;
    for(TopExp_Explorer ex(m_shape, TopAbs_FACE); ex.More(); ex.Next(), ++patchIndex)
    {
      const TopoDS_Face& F = TopoDS::Face(ex.Current());

      // Use BRepCheck_Analyzer to validate the face
      BRepCheck_Analyzer analyzer(F, true);

      const bool face_is_valid = analyzer.IsValid(F);
      if(face_is_valid)
      {
        continue;
      }

      printStatusList(axom::fmt::format("[Patch {} analyzer status]", patchIndex),
                      analyzer.Result(F)->Status());

      BRepCheck_Face fc(F);
      fc.Blind();  // run full tests
      printStatusList(axom::fmt::format("[Patch {} face check]", patchIndex), fc.Status());

      // Validate wires and edges of face F
      int wireIdx = 0;
      for(TopExp_Explorer wireExp(F, TopAbs_WIRE); wireExp.More(); wireExp.Next(), ++wireIdx)
      {
        const TopoDS_Wire& W = TopoDS::Wire(wireExp.Current());
        printStatusList(axom::fmt::format("[Patch {} Wire{} wire check]", patchIndex, wireIdx),
                        analyzer.Result(W)->Status());

        int edgeIdx = 0;
        for(TopExp_Explorer edgeExp(W, TopAbs_EDGE); edgeExp.More(); edgeExp.Next(), ++edgeIdx)
        {
          const TopoDS_Edge& E = TopoDS::Edge(edgeExp.Current());
          if(!analyzer.IsValid(E))
          {
            printStatusList(axom::fmt::format("[Patch {} Wire {} Edge {} analyzer status]",
                                              patchIndex,
                                              wireIdx,
                                              edgeIdx),
                            analyzer.Result(E)->Status());
          }
        }
      }
    }
  }

  /// Extracts data from the trimming curves of each patch and converts the curves to a NURBSCurve representation
  void extractTrimmingCurves(NPatchArray& patches)
  {
    std::map<GeomAbs_CurveType, std::string> curveTypeMap = {{GeomAbs_Line, "Line"},
                                                             {GeomAbs_Circle, "Circle"},
                                                             {GeomAbs_Ellipse, "Ellipse"},
                                                             {GeomAbs_Hyperbola, "Hyperbola"},
                                                             {GeomAbs_Parabola, "Parabola"},
                                                             {GeomAbs_BezierCurve, "Bezier Curve"},
                                                             {GeomAbs_BSplineCurve, "BSpline Curve"},
                                                             {GeomAbs_OffsetCurve, "Offset Curve"},
                                                             {GeomAbs_OtherCurve, "Other Curve"}};

    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More(); faceExp.Next(), ++patchIndex)
    {
      PatchData& patchData = m_patchData[patchIndex];
      auto& patch = patches[patchIndex];

      // Get span of this patch in u and v directions
      BBox2D patchBbox = patchData.parametricBBox;
      auto expandedPatchBbox = patchBbox;
      expandedPatchBbox.scale(1. + 1e-3);

      SLIC_INFO_IF(m_verbose,
                   axom::fmt::format("[Patch {}]: BBox in parametric space: {}; expanded BBox {}",
                                     patchIndex,
                                     patchBbox,
                                     expandedPatchBbox));

      int wireIndex = 0;
      for(TopExp_Explorer wireExp(faceExp.Current(), TopAbs_WIRE); wireExp.More();
          wireExp.Next(), ++wireIndex)
      {
        const TopoDS_Wire& wire = TopoDS::Wire(wireExp.Current());

        int edgeIndex = 0;
        for(TopExp_Explorer edgeExp(wire, TopAbs_EDGE); edgeExp.More(); edgeExp.Next(), ++edgeIndex)
        {
          const TopoDS_Edge& edge = TopoDS::Edge(edgeExp.Current());
          const int curveIndex = patch.getNumTrimmingCurves();

          TopAbs_Orientation orientation = edge.Orientation();
          const bool isReversed = (orientation == TopAbs_REVERSED);

          if(m_verbose)
          {
            BRepAdaptor_Curve curveAdaptor(edge);
            SLIC_INFO(axom::fmt::format("[Patch {} Wire {} Edge {} Curve {}] Curve type: '{}'",
                                        patchIndex,
                                        wireIndex,
                                        edgeIndex,
                                        curveIndex,
                                        curveTypeMap[curveAdaptor.GetType()]));
          }

          Standard_Real first, last;
          opencascade::handle<Geom2d_Curve> parametricCurve =
            BRep_Tool::CurveOnSurface(edge, TopoDS::Face(faceExp.Current()), first, last);
          opencascade::handle<Geom2d_BSplineCurve> bsplineCurve =
            Geom2dConvert::CurveToBSplineCurve(parametricCurve);

          if(!parametricCurve.IsNull() && !bsplineCurve.IsNull())
          {
            const bool originalCurvePeriodic = bsplineCurve->IsPeriodic();

            CurveProcessor curveProcessor(bsplineCurve, m_verbose);
            auto curve = curveProcessor.nurbsCurve();
            patchData.trimmingCurves_originallyPeriodic.push_back(
              curveProcessor.curveWasOriginallyPeriodic());

            if(isReversed)  // Ensure consistency of curve w.r.t. patch
            {
              curve.reverseOrientation();
            }
            SLIC_ASSERT(curve.isValidNURBS());
            SLIC_ASSERT(curve.getDegree() == bsplineCurve->Degree());

            patch.addTrimmingCurve(curve);

            SLIC_INFO_IF(m_verbose,
                         axom::fmt::format("[Patch {} Wire {} Edge {} Curve {}] Added curve: {}",
                                           patchIndex,
                                           wireIndex,
                                           edgeIndex,
                                           curveIndex,
                                           curve));

            // Check to ensure that curve did not change geometrically after making non-periodic
            if(originalCurvePeriodic)
            {
              opencascade::handle<Geom2d_BSplineCurve> origCurve =
                Geom2dConvert::CurveToBSplineCurve(parametricCurve);
              const bool withinThreshold = curveProcessor.compareToCurve(origCurve, 25);
              SLIC_WARNING_IF(
                !withinThreshold,
                axom::fmt::format("[Patch {} Wire {} Edge {} Curve {}] Trimming curve was not "
                                  "within threshold after clamping.",
                                  patchIndex,
                                  wireIndex,
                                  edgeIndex,
                                  curveIndex));
            }

            // TODO: Check that curve control points are within UV patch after adjusting periodicity
          }
        }
      }

      // If the face is flipped, then the trimming curves all need to be reversed too
      if(patch.isTrimmed() &&
         TopoDS::Face(faceExp.Current()).Orientation() == TopAbs_Orientation::TopAbs_REVERSED)
      {
        patch.reverseTrimmingCurves();
      }
    }
  }

  const PatchDataMap& getPatchDataMap() const { return m_patchData; }

  std::string getFileUnits() const { return m_fileUnits; }

private:
  /// Returns the canonical representation of a unit string (e.g. "centimeter" -> "cm")
  std::string getCanonicalUnit(const std::string& unit) const
  {
    // we'll convert all units to lower case
    auto toLower = [](std::string str) {
      std::transform(str.begin(), str.end(), str.begin(), ::tolower);
      return str;
    };

    // start with imperial units
    std::map<std::string, std::string> unitCanonicalMap = {{"inch", "in"},
                                                           {"inches", "in"},
                                                           {"in", "in"},
                                                           {"foot", "ft"},
                                                           {"feet", "ft"},
                                                           {"ft", "ft"},
                                                           {"mile", "mi"},
                                                           {"miles", "mi"},
                                                           {"mi", "mi"}};

    // now add the SI units w/ several suffixes
    // we're going to reverse this for the map to canonical units
    std::map<std::string, std::string> prefixes = {
      {"am", "atto"},
      {"fm", "femto"},
      {"pm", "pico"},
      {"nm", "nano"},
      {"um", "micro"},
      {"mm", "milli"},
      {"cm", "centi"},
      {"dm", "deci"},
      {"m", ""},
      {"dam", "deca"},
      {"hm", "hecto"},
      {"km", "kilo"},
    };

    for(const auto& kv : prefixes)
    {
      const std::string& canonical = kv.first;
      const std::string& prefix = kv.second;
      unitCanonicalMap[canonical] = canonical;
      for(const std::string& suffix : {"meter", "meters", "metre", "metres"})
      {
        unitCanonicalMap[prefix + suffix] = canonical;
      }
    }

    return unitCanonicalMap[toLower(unit)];
  }

  /**
   *  Returns the conversion factor from an input unit to an output unit
   * 
   * \note Converts the units to their canonical form
   * \sa getCanonicalUnit
   */
  double getConversionFactor(const std::string& fileUnits, const std::string& defaultUnits = "mm") const
  {
    std::map<std::string, double> unitConversionMap = {{"am", 1e-15},
                                                       {"fm", 1e-12},
                                                       {"pm", 1e-9},
                                                       {"nm", 1e-6},
                                                       {"um", 1e-3},
                                                       {"mm", 1.0},
                                                       {"cm", 10.0},
                                                       {"dm", 100.0},
                                                       {"m", 1e3},
                                                       {"dam", 1e4},
                                                       {"hm", 1e5},
                                                       {"km", 1e6},
                                                       {"in", 25.4},
                                                       {"ft", 304.8},
                                                       {"mi", 1609344.0}};

    const double fileUnitFactor = unitConversionMap[getCanonicalUnit(fileUnits)];
    const double defaultUnitFactor = unitConversionMap[getCanonicalUnit(defaultUnits)];

    return fileUnitFactor / defaultUnitFactor;
  };

  /// Loads the step file \a filename from disk
  /// Uses the units from \a filename
  TopoDS_Shape loadStepFile(const std::string& filename)
  {
    STEPControl_Reader reader;

    IFSelect_ReturnStatus status = reader.ReadFile(filename.c_str());
    if(status != IFSelect_RetDone)
    {
      m_loadStatus = LoadStatus::FAILED_TO_READ;
      std::cerr << "Error: Cannot read the file." << std::endl;
      return TopoDS_Shape();
    }

    // adjust the units, as needed
    TColStd_SequenceOfAsciiString anUnitLengthNames;
    TColStd_SequenceOfAsciiString anUnitAngleNames;
    TColStd_SequenceOfAsciiString anUnitSolidAngleNames;
    reader.FileUnits(anUnitLengthNames, anUnitAngleNames, anUnitSolidAngleNames);
    if(anUnitLengthNames.Size() > 0)
    {
      m_fileUnits = getCanonicalUnit(anUnitLengthNames(1).ToCString());
      std::string defaultUnit = Interface_Static::CVal("xstep.cascade.unit");
      const double lengthUnit = getConversionFactor(m_fileUnits, defaultUnit);
      reader.SetSystemLengthUnit(lengthUnit);
    }

    Standard_Integer numRoots = reader.NbRootsForTransfer();
    reader.TransferRoots();
    TopoDS_Shape shape = reader.OneShape();
    if(shape.IsNull())
    {
      m_loadStatus = LoadStatus::FAILED_NO_SHAPES;
      std::cerr << "Error: No shape found in the file." << std::endl;
      return TopoDS_Shape();
    }

    // Convert to NURBS
    BRepBuilderAPI_NurbsConvert converter(shape);
    TopoDS_Shape nurbsShape = converter.Shape();

    if(nurbsShape.IsNull())
    {
      m_loadStatus = LoadStatus::FAILED_TO_CONVERT;
      std::cerr << "Error: Conversion to NURBS failed." << std::endl;
      return TopoDS_Shape();
    }

    m_loadStatus = LoadStatus::SUCCESS;
    SLIC_INFO_IF(m_verbose,
                 axom::fmt::format("Successfully read the STEP file with {} roots", numRoots));

    return nurbsShape;
  }

private:
  TopoDS_Shape m_shape;
  bool m_verbose {false};
  LoadStatus m_loadStatus {LoadStatus::UNINITIALIZED};

  std::string m_fileUnits {"mm"};

  PatchDataMap m_patchData;
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

  PatchTriangulator(const TopoDS_Shape& shape,
                    double deflection,
                    double angularDeflection,
                    bool deflectionIsRelative)
    : m_shape(shape)
    , m_deflection(deflection)
    , m_angularDeflection(angularDeflection)
    , m_deflectionIsRelative(deflectionIsRelative)
  {
    BRepTools::Clean(shape);
    BRepMesh_IncrementalMesh mesh(m_shape, m_deflection, m_deflectionIsRelative, m_angularDeflection);

    if(!mesh.IsDone())
    {
      throw std::runtime_error("Mesh generation failed.");
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

      const bool isReversed = (face.Orientation() == TopAbs_Orientation::TopAbs_REVERSED);

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

        if(isReversed)
        {
          std::swap(n1, n3);
        }

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

  /// Utility function to triangulate each patch, ignoring the trimming curves, and write it to disk as as STL mesh
  void triangulateUntrimmedPatches(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& output_mesh)
  {
    std::vector<int> patch_id;

    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More(); faceExp.Next(), ++patchIndex)
    {
      TopoDS_Face face = TopoDS::Face(faceExp.Current());

      const bool isReversed = (face.Orientation() == TopAbs_Orientation::TopAbs_REVERSED);

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
      BRepMesh_IncrementalMesh mesh(newFace, m_deflection, m_deflectionIsRelative, m_angularDeflection);

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

        if(isReversed)
        {
          std::swap(n1, n3);
        }

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
  bool m_deflectionIsRelative;
};

}  // end namespace internal

std::string STEPReader::getFileUnits() const { return m_stepProcessor->getFileUnits(); }

std::string STEPReader::getBRepStats() const
{
  // early return if the step file has not been loaded
  // (or, e.g. for the derived PSTEPReader class, if it's not loaded on this rank)
  if(!m_stepProcessor)
  {
    return "";
  }

  // Helper struct for simple stats over a collection of integers
  struct AccumStatistics
  {
    int min;
    int max;
    double mean;
    double stddev;
  };

  // Lambda to generate AccumStatistics for a list of integers
  auto computeStatistics = [](const std::vector<int>& data) -> AccumStatistics {
    AccumStatistics stats;
    stats.min = *std::min_element(data.begin(), data.end());
    stats.max = *std::max_element(data.begin(), data.end());

    const double sum = std::accumulate(data.begin(), data.end(), 0.0);
    const double sumSquared =
      std::accumulate(data.begin(), data.end(), 0.0, [](double a, double b) { return a + b * b; });
    stats.mean = sum / data.size();
    stats.stddev = std::sqrt(sumSquared / data.size() - stats.mean * stats.mean);
    return stats;
  };

  axom::fmt::memory_buffer out;

  axom::fmt::format_to(std::back_inserter(out), "Details about loaded mesh:\n");

  // summarize the number of patches and trimming curves
  {
    int totalTrimmingCurves = 0;
    for(const auto& patch : m_patches)
    {
      totalTrimmingCurves += patch.getNumTrimmingCurves();
    }

    axom::fmt::format_to(std::back_inserter(out),
                         " - Mesh has {} patches with a total of {} trimming curves\n",
                         m_patches.size(),
                         totalTrimmingCurves);
  }

  // compute and print the bounding box of the mesh in physical space
  {
    typename internal::StepFileProcessor::BBox3D meshBBox;
    for(const auto& [_, value] : m_stepProcessor->getPatchDataMap())
    {
      meshBBox.addBox(value.physicalBBox);
    }

    axom::fmt::format_to(std::back_inserter(out),
                         " - Bounding box of the mesh in physical space (in {}): {}\n",
                         this->getFileUnits(),
                         meshBBox);
  }

  // compute a histogram of the patch degrees w/ some additional info
  {
    struct Counts
    {
      int total {0};
      int rational {0};
      int periodic_u {0};
      int periodic_v {0};
    };

    std::map<std::pair<int, int>, Counts> patchDegrees;
    for(const auto& [_, value] : m_stepProcessor->getPatchDataMap())
    {
      const auto& patch = m_patches[value.patchIndex];
      auto& c = patchDegrees[{patch.getDegree_u(), patch.getDegree_v()}];
      c.total++;
      if(patch.isRational())
      {
        c.rational++;
      }
      if(value.wasOriginallyPeriodic_u)
      {
        c.periodic_u++;
      }
      if(value.wasOriginallyPeriodic_v)
      {
        c.periodic_v++;
      }
    }

    axom::fmt::format_to(std::back_inserter(out), " - Patch degree histogram:\n");
    for(const auto& [degs, counts] : patchDegrees)
    {
      axom::fmt::format_to(
        std::back_inserter(out),
        "   - Degree (u={}, v={}): {} patches ({} rational{}{})\n",
        degs.first,   // degree u
        degs.second,  // degree v
        counts.total,
        counts.rational,
        counts.periodic_u > 0 ? axom::fmt::format("; {} originally periodic in u", counts.periodic_u)
                              : "",
        counts.periodic_v > 0 ? axom::fmt::format("; {} originally periodic in v", counts.periodic_v)
                              : "");
    }
  }

  // Compute statistics on the number of spans per patch
  {
    std::vector<int> uSpansPerPatch;
    std::vector<int> vSpansPerPatch;
    for(const auto& patch : m_patches)
    {
      uSpansPerPatch.push_back(patch.getKnots_u().getNumKnotSpans());
      vSpansPerPatch.push_back(patch.getKnots_v().getNumKnotSpans());
    }

    AccumStatistics uSpansStats = computeStatistics(uSpansPerPatch);
    AccumStatistics vSpansStats = computeStatistics(vSpansPerPatch);

    axom::fmt::format_to(std::back_inserter(out),
                         " - Number of u-spans per patch:  min: {}, max: {}, "
                         "mean: {:.2f} (stdev: {:.2f})\n",
                         uSpansStats.min,
                         uSpansStats.max,
                         uSpansStats.mean,
                         uSpansStats.stddev);

    axom::fmt::format_to(std::back_inserter(out),
                         " - Number of v-spans per patch:  min: {}, max: {}, "
                         "mean: {:.2f} (stdev: {:.2f})\n",
                         vSpansStats.min,
                         vSpansStats.max,
                         vSpansStats.mean,
                         vSpansStats.stddev);
  }

  // Compute statistics on the number of trimming curves per patch
  {
    std::vector<int> trimmingCurvesPerPatch;
    for(const auto& patch : m_patches)
    {
      trimmingCurvesPerPatch.push_back(patch.getNumTrimmingCurves());
    }

    AccumStatistics trimmingCurvesStats = computeStatistics(trimmingCurvesPerPatch);

    axom::fmt::format_to(std::back_inserter(out),
                         " - Number of trimming curves per patch:  min: {}, "
                         "max: {}, mean: {:.2f} (stdev: {:.2f})\n",
                         trimmingCurvesStats.min,
                         trimmingCurvesStats.max,
                         trimmingCurvesStats.mean,
                         trimmingCurvesStats.stddev);
  }

  // Compute statistics on the degrees of the trimming curves in the mesh
  {
    struct Counts
    {
      int total {0};
      int rational {0};
      int periodic {0};
    };

    std::map<int, Counts> curveDegreeCounts;
    std::vector<int> curveDegreeList;
    for(const auto& [_, value] : m_stepProcessor->getPatchDataMap())
    {
      const auto& patch = m_patches[value.patchIndex];

      int idx = 0;
      for(const auto& curve : patch.getTrimmingCurves())
      {
        const int degree = curve.getDegree();
        auto& c = curveDegreeCounts[degree];
        c.total++;
        if(curve.isRational())
        {
          c.rational++;
        }
        if(value.trimmingCurves_originallyPeriodic[idx])
        {
          c.periodic++;
        }

        curveDegreeList.push_back(degree);
        ++idx;
      }
    }

    AccumStatistics curveDegreeStats = computeStatistics(curveDegreeList);

    // Output the results for curve orders
    axom::fmt::format_to(std::back_inserter(out), " - Mesh trimming curve degree histogram:\n");
    for(const auto& [deg, counts] : curveDegreeCounts)
    {
      axom::fmt::format_to(
        std::back_inserter(out),
        "   - Degree {}: {} curves ({} rational{})\n",
        deg,
        counts.total,
        counts.rational,
        counts.periodic > 0 ? axom::fmt::format("; {} originally periodic", counts.periodic) : "");
    }

    axom::fmt::format_to(std::back_inserter(out),
                         "   - Average trimming curve order: {:.2f} (stdev: {:.2f})\n",
                         curveDegreeStats.mean,
                         curveDegreeStats.stddev);
  }

  // Compute statistics on the number of spans in the trimming curves
  {
    std::vector<int> spansPerCurve;
    for(const auto& patch : m_patches)
    {
      for(const auto& curve : patch.getTrimmingCurves())
      {
        spansPerCurve.push_back(curve.getKnots().getNumKnotSpans());
      }
    }

    AccumStatistics spansStats = computeStatistics(spansPerCurve);

    axom::fmt::format_to(std::back_inserter(out),
                         " - Number of spans per trimming curve:  min: {}, "
                         "max: {}, mean: {:.2f} (stdev: {:.2f})\n",
                         spansStats.min,
                         spansStats.max,
                         spansStats.mean,
                         spansStats.stddev);
  }

  return axom::fmt::to_string(out);
}

STEPReader::~STEPReader()
{
  if(m_stepProcessor)
  {
    delete m_stepProcessor;
    m_stepProcessor = nullptr;
  }
}

int STEPReader::read(bool validate_model)
{
  m_stepProcessor = new internal::StepFileProcessor(m_fileName, m_verbosity);
  if(!m_stepProcessor->isLoaded())
  {
    return 1;
  }

  if(validate_model)
  {
    m_stepProcessor->validateBRep();
  }

  m_stepProcessor->extractPatches(m_patches);
  m_stepProcessor->extractTrimmingCurves(m_patches);

  return 0;
}

int STEPReader::getTriangleMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* mesh,
                                double linear_deflection,
                                double angular_deflection,
                                bool is_relative,
                                bool trimmed)
{
  if(!m_stepProcessor || !m_stepProcessor->isLoaded())
  {
    SLIC_WARNING("Cannot triangulate model until calling STEPReader::read()");
    return 1;
  }

  if(!mesh)
  {
    SLIC_WARNING("Passed in mesh instance was null. Skipping triangulation");
    return 1;
  }

  internal::PatchTriangulator patchTriangulator(m_stepProcessor->getShape(),
                                                linear_deflection,
                                                angular_deflection,
                                                is_relative);

  if(trimmed)
  {
    patchTriangulator.triangulateTrimmedPatches(*mesh);
  }
  else
  {
    patchTriangulator.triangulateUntrimmedPatches(*mesh);
  }

  return 0;
}

}  // end namespace quest
}  // end namespace axom
