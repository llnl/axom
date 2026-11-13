// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/STEPReader.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#ifndef AXOM_USE_OPENCASCADE
  #error STEPReader should only be included when Axom is configured with opencascade
#endif

#include "axom/slic.hpp"
#include "axom/fmt.hpp"
#include "axom/core/utilities/Units.hpp"

#include "opencascade/BRep_Tool.hxx"
#include "opencascade/BRepBuilderAPI_MakeFace.hxx"
#include "opencascade/BRepBuilderAPI_NurbsConvert.hxx"
#include "opencascade/BRepBuilderAPI_Copy.hxx"
#include "opencascade/BRepTools_NurbsConvertModification.hxx"
#include "opencascade/BndLib_Add2dCurve.hxx"
#include "opencascade/Bnd_Box2d.hxx"
#include "opencascade/BRepCheck_Analyzer.hxx"
#include "opencascade/BRepCheck_Face.hxx"
#include "opencascade/BRepMesh_IncrementalMesh.hxx"
#include "opencascade/BRepTools.hxx"
#include "opencascade/BRepTools_WireExplorer.hxx"
#include "opencascade/BRepBndLib.hxx"
#include "opencascade/Geom_BSplineSurface.hxx"
#include "opencascade/Geom_BSplineCurve.hxx"
#include "opencascade/Geom_Circle.hxx"
#include "opencascade/Geom_Ellipse.hxx"
#include "opencascade/Geom_Line.hxx"
#include "opencascade/GeomConvert.hxx"
#include "opencascade/GeomAdaptor_Curve.hxx"
#include "opencascade/GeomAdaptor_Surface.hxx"
#include "opencascade/ProjLib_ProjectedCurve.hxx"
#include "opencascade/ProjLib_ComputeApproxOnPolarSurface.hxx"
#include "opencascade/GeomAPI_ProjectPointOnCurve.hxx"
#include "opencascade/Geom_RectangularTrimmedSurface.hxx"
#include "opencascade/Geom_Surface.hxx"
#include "opencascade/Geom2d_BSplineCurve.hxx"
#include "opencascade/Geom2d_Curve.hxx"
#include "opencascade/Geom2d_TrimmedCurve.hxx"
#include "opencascade/Geom2dConvert.hxx"
#include "opencascade/Interface_Static.hxx"
#include "opencascade/Poly_Triangulation.hxx"
#include "opencascade/Precision.hxx"
#include "opencascade/Standard_ConstructionError.hxx"
#include "opencascade/STEPControl_Reader.hxx"
#include "opencascade/TopExp_Explorer.hxx"
#include "opencascade/TopLoc_Location.hxx"
#include "opencascade/TopoDS_Edge.hxx"
#include "opencascade/TopoDS_Face.hxx"
#include "opencascade/TopoDS_Shape.hxx"
#include "opencascade/TopoDS_Wire.hxx"
#include "opencascade/TopoDS.hxx"
#include "opencascade/TopTools_MapOfShape.hxx"
#include "opencascade/gp_Vec2d.hxx"
#include "opencascade/TColgp_Array1OfPnt2d.hxx"
#include "opencascade/TColgp_Array2OfPnt.hxx"
#include "opencascade/TColStd_Array1OfInteger.hxx"
#include "opencascade/TColStd_Array1OfReal.hxx"
#include "opencascade/TColStd_Array2OfReal.hxx"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

namespace axom::quest
{
namespace internal
{
namespace
{
constexpr int CurveDim = 2;
constexpr int SpaceDim = 3;
using NCurve = axom::primal::NURBSCurve<double, CurveDim>;

using NPatch = axom::primal::NURBSPatch<double, SpaceDim>;
using NPatchArray = axom::Array<NPatch>;

using PointType = axom::primal::Point<double, CurveDim>;
using PointType3D = axom::primal::Point<double, SpaceDim>;
using BBox2D = axom::primal::BoundingBox<double, CurveDim>;
using BBox3D = axom::primal::BoundingBox<double, SpaceDim>;

/// Accumulate bounds of the curve itself, rather than its control polygon.
void addPcurveBounds(Bnd_Box2d& bounds, const opencascade::handle<Geom2d_BSplineCurve>& curve)
{
  BndLib_Add2dCurve::AddOptimal(curve,
                                curve->FirstParameter(),
                                curve->LastParameter(),
                                Precision::PConfusion(),
                                bounds);
}

BBox2D uvBounds(Bnd_Box2d bounds)
{
  bounds.SetGap(0.);
  PointType lo, hi;
  bounds.Get(lo[0], lo[1], hi[0], hi[1]);
  return BBox2D(lo, hi);
}

// A control-point box is cheaper when only deciding whether to try interval projection.
BBox2D pcurveControlBounds(const opencascade::handle<Geom2d_BSplineCurve>& curve)
{
  BBox2D bounds;
  for(int i = 1; i <= curve->NbPoles(); ++i)
  {
    const auto p = curve->Pole(i);
    bounds.addPoint(PointType {p.X(), p.Y()});
  }
  return bounds;
}

bool spansMultiplePeriods(const opencascade::handle<Geom_BSplineSurface>& surface,
                          const BBox2D& bounds)
{
  const auto extent = bounds.range();
  return (surface->IsUPeriodic() && extent[0] > surface->UPeriod() + Precision::PConfusion()) ||
    (surface->IsVPeriodic() && extent[1] > surface->VPeriod() + Precision::PConfusion());
}

opencascade::handle<Geom2d_BSplineCurve> trimPcurve(const opencascade::handle<Geom2d_Curve>& curve,
                                                    double first,
                                                    double last)
{
  const opencascade::handle<Geom2d_TrimmedCurve> trimmed =
    new Geom2d_TrimmedCurve(curve, first, last, true, false);
  const auto spline = opencascade::handle<Geom2d_BSplineCurve>::DownCast(trimmed->BasisCurve());
  if(!spline.IsNull())
  {
    // The trimmed curve owns this copy. Retain Open Cascade's default tolerance for ordinary
    // ranges, but keep the endpoint tolerance neighborhoods disjoint on short ranges.
    const double tolerance = std::min(Precision::PConfusion(),
                                      0.25 * (trimmed->LastParameter() - trimmed->FirstParameter()));
    spline->Segment(trimmed->FirstParameter(), trimmed->LastParameter(), tolerance);
    return spline;
  }
  return Geom2dConvert::CurveToBSplineCurve(trimmed);
}

// Include samples in every knot span, even when a span is short.
axom::Array<double> pcurveSamples(const opencascade::handle<Geom2d_BSplineCurve>& curve)
{
  constexpr int numSamples = 1001;
  constexpr int samplesPerSpan = 21;
  axom::Array<double> samples(numSamples);
  const double first = curve->FirstParameter(), last = curve->LastParameter();
  axom::numerics::linspace(first, last, samples.data(), numSamples);
  for(int i = 1; i < curve->NbKnots(); ++i)
  {
    const double a = std::max(first, curve->Knot(i));
    const double b = std::min(last, curve->Knot(i + 1));
    if(b > a)
    {
      double spanSamples[samplesPerSpan];
      axom::numerics::linspace(a, b, spanSamples, samplesPerSpan);
      samples.append(axom::ArrayView<double>(spanSamples, samplesPerSpan));
    }
  }
  return samples;
}

/// Tolerance for accepting a pcurve of \a edge on \a face, in physical space.
double pcurveTolerance(const TopoDS_Edge& edge, const TopoDS_Face& face)
{
  return std::max({Precision::Confusion(), BRep_Tool::Tolerance(face), BRep_Tool::Tolerance(edge)});
}

/// Returns the whole-period translation that moves \a from nearest to \a to.
gp_Vec2d periodShift(const gp_Pnt2d& from,
                     const gp_Pnt2d& to,
                     const opencascade::handle<Geom_BSplineSurface>& surface)
{
  const double du = surface->IsUPeriodic()
    ? std::round((to.X() - from.X()) / surface->UPeriod()) * surface->UPeriod()
    : 0.;
  const double dv = surface->IsVPeriodic()
    ? std::round((to.Y() - from.Y()) / surface->VPeriod()) * surface->VPeriod()
    : 0.;
  return gp_Vec2d(du, dv);
}

void alignPeriodicPcurve(const opencascade::handle<Geom2d_BSplineCurve>& curve,
                         const gp_Pnt2d& start,
                         const opencascade::handle<Geom_BSplineSurface>& surface)
{
  curve->Translate(periodShift(curve->Value(curve->FirstParameter()), start, surface));
}

// Join nonrational C1 pieces without rescaling their original knot parameters.
opencascade::handle<Geom2d_BSplineCurve> joinPcurveIntervals(
  axom::ArrayView<const opencascade::handle<Geom2d_BSplineCurve>> pieces,
  int degree,
  double last)
{
  int numPoles = 1, numKnots = 1;
  for(auto& piece : pieces)
  {
    piece->IncreaseDegree(degree);
    numPoles += piece->NbPoles() - 1;
    numKnots += piece->NbKnots() - 1;
  }
  TColgp_Array1OfPnt2d poles(1, numPoles);
  TColStd_Array1OfReal knots(1, numKnots);
  TColStd_Array1OfInteger multiplicities(1, numKnots);
  int pole = 1, knot = 1;
  for(const auto& piece : pieces)
  {
    for(int i = pole == 1 ? 1 : 2; i <= piece->NbPoles(); ++i)
    {
      poles(pole++) = piece->Pole(i);
    }
    for(int i = 1; i < piece->NbKnots(); ++i)
    {
      knots(knot) = piece->Knot(i);
      multiplicities(knot) = knot > 1 && i == 1 ? degree : piece->Multiplicity(i);
      ++knot;
    }
  }
  knots(knot) = last;
  multiplicities(knot) = degree + 1;
  return new Geom2d_BSplineCurve(poles, knots, multiplicities, degree);
}

/// A 3D edge curve restricted to its topological parameter interval.
struct EdgeCurve
{
  explicit EdgeCurve(const TopoDS_Edge& edge) { curve = BRep_Tool::Curve(edge, first, last); }

  opencascade::handle<Geom_Curve> curve;
  double first {}, last {};
};

/// Project and validate pcurves against one source edge and one periodic surface.
/// Existing pcurves may use different parameters; new projections must match the source edge.
class PeriodicPcurveProjector
{
public:
  PeriodicPcurveProjector(const EdgeCurve& source,
                          const opencascade::handle<Geom_BSplineSurface>& surface,
                          double tolerance)
    : m_source(source)
    , m_surface(surface)
    , m_tolerance(tolerance)
  { }

  opencascade::handle<Geom2d_BSplineCurve> project(const gp_Pnt2d& oldStart,
                                                   bool retryIntervals = false) const
  {
    // Bound the projection search by the source tolerance used to accept the pcurve.
    // Open Cascade defaults to a much larger distance and can spend time fitting invalid guesses.
    ProjLib_ProjectedCurve projection;
    projection.Load(new GeomAdaptor_Surface(m_surface));
    projection.Load(std::max(1e-10, m_tolerance * .05));
    projection.SetMaxDist(m_tolerance);
    projection.Perform(new GeomAdaptor_Curve(m_source.curve, m_source.first, m_source.last));
    const auto projected = projection.BSpline();
    if(projected.IsNull())
    {
      return {};
    }
    auto candidate = trimPcurve(projected, m_source.first, m_source.last);
    if(alignAndCheck(candidate, oldStart))
    {
      return candidate;
    }
    if(retryIntervals && m_source.curve->IsKind(STANDARD_TYPE(Geom_BSplineCurve)))
    {
      candidate = projectIntervals();
      if(alignAndCheck(candidate, oldStart))
      {
        return candidate;
      }
    }
    return {};
  }

  // Keep the converted pcurve unless a replacement passes the source tolerance.
  opencascade::handle<Geom2d_BSplineCurve> refine(
    const opencascade::handle<Geom2d_BSplineCurve>& curve) const
  {
    try
    {
      if(!needsRefinement(curve))
      {
        return curve;
      }
      // A control-point box suffices to identify possible multi-period trims.
      // The interval fit still has to pass the source tolerance check.
      const bool retryIntervals = m_source.curve->IsKind(STANDARD_TYPE(Geom_BSplineCurve)) &&
        spansMultiplePeriods(m_surface, pcurveControlBounds(curve));
      const auto candidate = project(curve->Value(curve->FirstParameter()), retryIntervals);
      if(!candidate.IsNull())
      {
        return candidate;
      }
    }
    catch(const Standard_Failure&)
    {
      // Keep the converted pcurve if projection fails.
    }
    return curve;
  }

private:
  gp_Pnt pointOnSurface(const opencascade::handle<Geom2d_BSplineCurve>& curve, double t) const
  {
    const auto uv = curve->Value(t);
    return m_surface->Value(uv.X(), uv.Y());
  }

  bool needsRefinement(const opencascade::handle<Geom2d_BSplineCurve>& curve) const
  {
    double error = 0.;
    const auto start = m_source.curve->Value(m_source.first),
               end = m_source.curve->Value(m_source.last);
    for(double t : pcurveSamples(curve))
    {
      const auto xyz = pointOnSurface(curve, t);
      // A nearby point on the source edge bounds the nearest-point distance.
      if(t >= m_source.first && t <= m_source.last &&
         xyz.Distance(m_source.curve->Value(t)) <= m_tolerance)
      {
        continue;
      }
      GeomAPI_ProjectPointOnCurve projection(xyz, m_source.curve, m_source.first, m_source.last);
      if(projection.NbPoints() == 0)
      {
        // Without a reliable error estimate, retain the converted pcurve.
        return false;
      }
      error =
        std::max(error,
                 std::min({xyz.Distance(start), xyz.Distance(end), projection.LowerDistance()}));
      // Further samples cannot change the decision to try refining this pcurve.
      if(error > m_tolerance)
      {
        break;
      }
    }
    error = std::max({error,
                      pointOnSurface(curve, curve->FirstParameter()).Distance(start),
                      pointOnSurface(curve, curve->LastParameter()).Distance(end)});
    return !(error <= m_tolerance);
  }

  bool alignAndCheck(const opencascade::handle<Geom2d_BSplineCurve>& candidate,
                     const gp_Pnt2d& oldStart) const
  {
    if(candidate.IsNull())
    {
      return false;
    }
    alignPeriodicPcurve(candidate, oldStart, m_surface);
    // Projection preserves the source parameterization, so compare at matching parameters.
    for(double t : pcurveSamples(candidate))
    {
      const auto uv = candidate->Value(t);
      if(!(m_surface->Value(uv.X(), uv.Y()).Distance(m_source.curve->Value(t)) <= m_tolerance))
      {
        return false;
      }
    }
    return true;
  }

  // OCCT joins C1 intervals with at most a one-period shift. Multi-period trims
  // need the full integer shift at each join to preserve their winding.
  opencascade::handle<Geom2d_BSplineCurve> projectIntervals() const
  {
    const GeomAdaptor_Curve adaptor(m_source.curve, m_source.first, m_source.last);
    const int numIntervals = adaptor.NbIntervals(GeomAbs_C1);
    // A single interval has no periodic joins to repair.
    if(numIntervals <= 1)
    {
      return {};
    }
    TColStd_Array1OfReal intervals(1, numIntervals + 1);
    adaptor.Intervals(intervals, GeomAbs_C1);
    const opencascade::handle<GeomAdaptor_Surface> surfaceAdaptor =
      new GeomAdaptor_Surface(m_surface);
    axom::Array<opencascade::handle<Geom2d_BSplineCurve>> pieces;
    int degree = 0;
    for(int i = intervals.Lower(); i < intervals.Upper(); ++i)
    {
      ProjLib_ComputeApproxOnPolarSurface projection;
      // The higher-level projector floors its fitting tolerance at Precision::Confusion().
      // Use a tighter fit here, but retain the source tolerance for search and acceptance.
      projection.SetTolerance(std::max(1e-10, m_tolerance * .05));
      projection.SetMaxDist(m_tolerance);
      projection.Perform(new GeomAdaptor_Curve(m_source.curve, intervals(i), intervals(i + 1)),
                         surfaceAdaptor);
      const auto projected = projection.BSpline();
      if(projected.IsNull() || projected->IsRational())
      {
        return {};
      }
      auto piece = trimPcurve(projected, intervals(i), intervals(i + 1));
      if(!pieces.empty())
      {
        const auto& previous = pieces.back();
        alignPeriodicPcurve(piece, previous->Value(previous->LastParameter()), m_surface);
      }
      degree = std::max(degree, piece->Degree());
      pieces.push_back(piece);
    }

    return joinPcurveIntervals(pieces.view(), degree, m_source.last);
  }

  EdgeCurve m_source;
  opencascade::handle<Geom_BSplineSurface> m_surface;
  double m_tolerance;
};

opencascade::handle<Geom2d_BSplineCurve> refinePeriodicPcurve(
  const opencascade::handle<Geom2d_BSplineCurve>& curve,
  const opencascade::handle<Geom_BSplineSurface>& surface,
  const TopoDS_Edge& sourceEdge,
  const TopoDS_Face& sourceFace)
{
  if(BRep_Tool::Degenerated(sourceEdge))
  {
    return curve;
  }
  const EdgeCurve original(sourceEdge);
  if(original.curve.IsNull())
  {
    return curve;
  }
  return PeriodicPcurveProjector(original, surface, pcurveTolerance(sourceEdge, sourceFace)).refine(curve);
}

// Open Cascade can convert an analytic periodic surface to a closed, nonperiodic spline.
// Preserve its periodic basis so pcurve projection can cross the seam.
class PeriodicNurbsModification : public BRepTools_NurbsConvertModification
{
public:
  Standard_Boolean NewSurface(const TopoDS_Face& face,
                              opencascade::handle<Geom_Surface>& surface,
                              TopLoc_Location& location,
                              Standard_Real& tolerance,
                              Standard_Boolean& reverseWires,
                              Standard_Boolean& reverseFace) override
  {
    const bool changed = BRepTools_NurbsConvertModification::NewSurface(face,
                                                                        surface,
                                                                        location,
                                                                        tolerance,
                                                                        reverseWires,
                                                                        reverseFace);
    if(!changed)
    {
      return false;
    }
    auto source = BRep_Tool::Surface(face, location);
    // A trimmed wrapper can hide the periodic basis even for a full-period face.
    const auto trimmed = opencascade::handle<Geom_RectangularTrimmedSurface>::DownCast(source);
    if(!trimmed.IsNull())
    {
      source = trimmed->BasisSurface();
    }
    const auto spline = opencascade::handle<Geom_BSplineSurface>::DownCast(surface);
    if(spline.IsNull())
    {
      return changed;
    }
    if((source->IsUPeriodic() && spline->IsUClosed() && !spline->IsUPeriodic()) ||
       (source->IsVPeriodic() && spline->IsVClosed() && !spline->IsVPeriodic()))
    {
      const auto stockSurface = surface;
      try
      {
        surface = restorePeriodicity(face, source, spline, !trimmed.IsNull());
        m_periodicFaces.Add(face);
      }
      catch(const Standard_Failure&)
      {
        surface = stockSurface;
        SLIC_WARNING_ROOT(
          "Could not preserve a periodic STEP surface; using Open Cascade's default conversion.");
      }
    }
    return changed;
  }

  Standard_Boolean NewCurve2d(const TopoDS_Edge& edge,
                              const TopoDS_Face& face,
                              const TopoDS_Edge& newEdge,
                              const TopoDS_Face& newFace,
                              opencascade::handle<Geom2d_Curve>& curve,
                              Standard_Real& tolerance) override
  {
    double first, last;
    const auto source = BRep_Tool::Surface(face);
    const auto curve3d = BRep_Tool::Curve(edge, first, last);
    // Try direct projection for circles, ellipses and lines on wrapped periodic surfaces.
    // Open Cascade's seeded polar-surface approximation is expensive for these curves.
    if(m_periodicFaces.Contains(face) &&
       source->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface)) && !curve3d.IsNull() &&
       (curve3d->IsKind(STANDARD_TYPE(Geom_Circle)) ||
        curve3d->IsKind(STANDARD_TYPE(Geom_Ellipse)) || curve3d->IsKind(STANDARD_TYPE(Geom_Line))))
    {
      const double sourceTolerance = pcurveTolerance(edge, face);
      const auto projected = projectAnalyticPcurve(edge, face, newEdge, newFace, sourceTolerance);
      if(!projected.IsNull())
      {
        curve = projected;
        tolerance = sourceTolerance;
        return true;
      }
    }
    // Open Cascade leaves spline surfaces and 3D spline curves unchanged.
    // Reuse their pcurves to avoid projection errors at seams.
    if(source->IsKind(STANDARD_TYPE(Geom_BSplineSurface)) && !curve3d.IsNull() &&
       curve3d->IsKind(STANDARD_TYPE(Geom_BSplineCurve)))
    {
      try
      {
        const auto original = BRep_Tool::CurveOnSurface(edge, face, first, last);
        if(!original.IsNull())
        {
          curve = trimPcurve(original, first, last);
          tolerance = BRep_Tool::Tolerance(edge);
          return true;
        }
      }
      catch(const Standard_Failure&)
      {
        SLIC_WARNING_ROOT(
          "Could not preserve a STEP spline pcurve; using Open Cascade's default conversion.");
      }
    }
    return BRepTools_NurbsConvertModification::NewCurve2d(edge, face, newEdge, newFace, curve, tolerance);
  }

private:
  // Wrapped surfaces must retain the converted parameterization. Unwrapped
  // analytic surfaces can be converted from a basis trimmed only in nonperiodic directions.
  static opencascade::handle<Geom_BSplineSurface> restorePeriodicity(
    const TopoDS_Face& face,
    const opencascade::handle<Geom_Surface>& source,
    const opencascade::handle<Geom_BSplineSurface>& spline,
    bool wrapped)
  {
    if(wrapped)
    {
      // Keep the converted surface's parameterization when restoring periodicity.
      const auto periodic = opencascade::handle<Geom_BSplineSurface>::DownCast(spline->Copy());
      if(source->IsUPeriodic() && periodic->IsUClosed() && !periodic->IsUPeriodic())
      {
        periodic->SetUPeriodic();
      }
      if(source->IsVPeriodic() && periodic->IsVClosed() && !periodic->IsVPeriodic())
      {
        periodic->SetVPeriodic();
      }
      return periodic;
    }
    else
    {
      double u0, u1, v0, v1;
      BRepTools::UVBounds(face, u0, u1, v0, v1);
      opencascade::handle<Geom_Surface> basis = source;
      if(source->IsUPeriodic() && !source->IsVPeriodic())
      {
        basis = new Geom_RectangularTrimmedSurface(source, v0, v1, false);
      }
      if(source->IsVPeriodic() && !source->IsUPeriodic())
      {
        basis = new Geom_RectangularTrimmedSurface(source, u0, u1, true);
      }
      return GeomConvert::SurfaceToBSplineSurface(basis);
    }
  }

  static opencascade::handle<Geom2d_BSplineCurve> projectAnalyticPcurve(const TopoDS_Edge& edge,
                                                                        const TopoDS_Face& face,
                                                                        const TopoDS_Edge& newEdge,
                                                                        const TopoDS_Face& newFace,
                                                                        double sourceTolerance)
  {
    try
    {
      const auto convertedSurface =
        opencascade::handle<Geom_BSplineSurface>::DownCast(BRep_Tool::Surface(newFace));
      const EdgeCurve convertedCurve(newEdge);
      double pf, pl;
      const auto seed = BRep_Tool::CurveOnSurface(edge, face, pf, pl);
      if(!convertedSurface.IsNull() && !convertedCurve.curve.IsNull() && !seed.IsNull())
      {
        const PeriodicPcurveProjector projector(convertedCurve, convertedSurface, sourceTolerance);
        const auto projected = projector.project(seed->Value(pf));
        if(!projected.IsNull())
        {
          // A projection that crosses multiple periods needs Open Cascade's guided conversion.
          Bnd_Box2d box;
          addPcurveBounds(box, projected);
          box.SetGap(0.);
          double u0, v0, u1, v1;
          box.Get(u0, v0, u1, v1);
          if((!convertedSurface->IsUPeriodic() ||
              u1 - u0 <= convertedSurface->UPeriod() + Precision::PConfusion()) &&
             (!convertedSurface->IsVPeriodic() ||
              v1 - v0 <= convertedSurface->VPeriod() + Precision::PConfusion()))
          {
            return projected;
          }
        }
      }
    }
    catch(const Standard_Failure&)
    {
      // The caller falls back to Open Cascade's conversion.
    }
    return {};
  }

  TopTools_MapOfShape m_periodicFaces;
};

class PeriodicNurbsConverter : public BRepBuilderAPI_NurbsConvert
{
public:
  explicit PeriodicNurbsConverter(const TopoDS_Shape& shape)
  {
    myModification = new PeriodicNurbsModification;
    Perform(shape);
  }
};

/// Own both conversion histories so source faces and edge occurrences remain traceable.
/// The copy keeps Open Cascade's conversion from modifying the source geometry.
class NurbsConversion
{
public:
  explicit NurbsConversion(const TopoDS_Shape& source)
    : m_copy(source, true, false)
    , m_converter(m_copy.Shape())
  { }

  const TopoDS_Shape& shape() { return m_converter.Shape(); }

  TopoDS_Shape convertedShape(const TopoDS_Shape& source) const
  {
    return m_converter.ModifiedShape(m_copy.ModifiedShape(source));
  }

private:
  BRepBuilderAPI_Copy m_copy;
  PeriodicNurbsConverter m_converter;
};

/// Struct to hold data associated with each surface patch of the mesh
struct PatchData
{
  int patchIndex {-1};
  bool wasOriginallyPeriodic_u {false};
  bool wasOriginallyPeriodic_v {false};
  bool was_reversed_u {false};
  axom::primal::BoundingBox<double, 2> parametricBBox;
  axom::primal::BoundingBox<double, 3> physicalBBox;
  axom::Array<bool> trimmingCurves_originallyPeriodic;
  axom::Array<int> trimmingCurves_wireIds;
};

using PatchDataMap = std::map<int, PatchData>;

/// One oriented edge occurrence, restricted and repaired on the converted surface.
/// Seam edges can contribute two occurrences with distinct pcurves.
struct PreparedPcurve
{
  opencascade::handle<Geom2d_BSplineCurve> curve;
  int wireIndex;
  int edgeIndex;
  bool reversed;
  bool originallyPeriodic;
};

// Traverse source edge occurrences: a converted wire can hide one side of a seam.
axom::Array<PreparedPcurve> preparePcurves(int patchIndex,
                                           TopoDS_Face source_face,
                                           TopoDS_Face face,
                                           const opencascade::handle<Geom_BSplineSurface>& surface,
                                           const NurbsConversion& conversion)
{
  source_face.Orientation(TopAbs_FORWARD);
  face.Orientation(TopAbs_FORWARD);
  // Open Cascade may already produce a periodic spline
  // its converted pcurves still need checking.
  // Source spline surfaces keep their parameterization.
  TopLoc_Location sourceLocation;
  const bool refine = (surface->IsUPeriodic() || surface->IsVPeriodic()) &&
    !BRep_Tool::Surface(source_face, sourceLocation)->IsKind(STANDARD_TYPE(Geom_BSplineSurface));
  axom::Array<PreparedPcurve> curves;
  int wireIndex = 0;
  for(TopExp_Explorer wireExp(source_face, TopAbs_WIRE); wireExp.More(); wireExp.Next(), ++wireIndex)
  {
    const auto& wire = TopoDS::Wire(wireExp.Current());
    int expected_edges = 0;
    for(TopExp_Explorer edgeExp(wire, TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
    {
      ++expected_edges;
    }

    // A converted wire's UV chart can hide a seam occurrence from WireExplorer.
    // Traverse the source wire to retain both occurrences.
    int edgeIndex = 0;
    std::optional<gp_Pnt2d> wireEnd;
    BRepTools_WireExplorer edgeExp(wire, source_face);
    for(; edgeExp.More(); edgeExp.Next(), ++edgeIndex)
    {
      const auto& source_edge = edgeExp.Current();
      auto edge = TopoDS::Edge(conversion.convertedShape(source_edge));
      edge.Orientation(source_edge.Orientation());
      const bool reversed = edge.Orientation() == TopAbs_REVERSED;
      double first, last;
      const auto pcurve = BRep_Tool::CurveOnSurface(edge, face, first, last);
      if(pcurve.IsNull())
      {
        SLIC_WARNING_ROOT(axom::fmt::format(
          "[Patch {} Wire {} Edge {}] Missing pcurve; the exported boundary is incomplete.",
          patchIndex,
          wireIndex,
          edgeIndex));
        continue;
      }

      opencascade::handle<Geom2d_BSplineCurve> spline;
      try
      {
        // Trim to the edge's interval without shifting periodic parameters.
        spline = trimPcurve(pcurve, first, last);
        if(refine)
        {
          spline = refinePeriodicPcurve(spline, surface, source_edge, source_face);
          // Each refined pcurve is aligned with its own converted pcurve, which can sit
          // in another period. Chain it to the previous occurrence to keep the wire closed.
          if(wireEnd)
          {
            const auto start =
              spline->Value(reversed ? spline->LastParameter() : spline->FirstParameter());
            spline->Translate(periodShift(start, *wireEnd, surface));
          }
        }
      }
      catch(const Standard_Failure& failure)
      {
        SLIC_WARNING_ROOT(
          axom::fmt::format("[Patch {} Wire {} Edge {}] Could not convert pcurve ({}); "
                            "the exported boundary is incomplete.",
                            patchIndex,
                            wireIndex,
                            edgeIndex,
                            failure.GetMessageString()));
        continue;
      }
      wireEnd = spline->Value(reversed ? spline->FirstParameter() : spline->LastParameter());
      curves.push_back({spline, wireIndex, edgeIndex, reversed, pcurve->IsPeriodic()});
    }
    SLIC_WARNING_ROOT_IF(
      edgeIndex != expected_edges,
      axom::fmt::format("[Patch {} Wire {}] Source traversal visited {} of {} "
                        "edge occurrences; the exported boundary may be incomplete.",
                        patchIndex,
                        wireIndex,
                        edgeIndex,
                        expected_edges));
  }
  return curves;
}

// Copy Open Cascade knot sequences w/ repeated knots for multiplicity to an axom::Array
axom::Array<double> copyKnotSequence(const TColStd_Array1OfReal& knots)
{
  return axom::Array<double>(axom::ArrayView<const double>(&knots.First(), knots.Length()));
}

BBox2D surfaceBounds(const opencascade::handle<Geom_BSplineSurface>& surface)
{
  PointType lo, hi;
  surface->Bounds(lo[0], hi[0], lo[1], hi[1]);
  return BBox2D(lo, hi);
}

/// Choose fixed u/v seams enclosing the prepared trims. Keep at least one full
/// period, expanding to multiple periods when wrapping would split a wire.
BBox2D clampingBounds(const opencascade::handle<Geom_BSplineSurface>& surface,
                      axom::ArrayView<const PreparedPcurve> curves,
                      int patchIndex)
{
  const auto bounds = surfaceBounds(surface);
  if(!surface->IsUPeriodic() && !surface->IsVPeriodic())
  {
    return bounds;
  }
  Bnd_Box2d box;
  try
  {
    for(const auto& c : curves)
    {
      addPcurveBounds(box, c.curve);
    }
  }
  catch(const Standard_Failure&)
  {
    SLIC_WARNING_ROOT(
      axom::fmt::format("[Patch {}] Could not bound periodic trims; keeping the native chart. "
                        "Trims may extend outside the patch.",
                        patchIndex));
    return bounds;
  }
  if(box.IsVoid())
  {
    return bounds;
  }
  const auto trimBounds = uvBounds(box);
  const auto& trimMin = trimBounds.getMin();
  const auto& trimMax = trimBounds.getMax();
  auto lo = bounds.getMin(), hi = bounds.getMax();
  const double uvTolerance = 1e-7;
  for(int dim = 0; dim < CurveDim; ++dim)
  {
    if(!(dim == 0 ? surface->IsUPeriodic() : surface->IsVPeriodic()))
    {
      continue;
    }
    const double period = dim == 0 ? surface->UPeriod() : surface->VPeriod();
    if(!std::isfinite(trimMin[dim]) || !std::isfinite(trimMax[dim]))
    {
      SLIC_WARNING_ROOT(
        axom::fmt::format("[Patch {}] Trims have nonfinite bounds; "
                          "keeping the native domain. Trims may extend outside the patch.",
                          patchIndex));
      continue;
    }
    if(trimMax[dim] - trimMin[dim] > period + uvTolerance)
    {
      lo[dim] = trimMin[dim];
      hi[dim] = trimMax[dim];
      continue;
    }
    if(trimMin[dim] < lo[dim] - uvTolerance || trimMax[dim] > hi[dim] + uvTolerance)
    {
      lo[dim] = (trimMin[dim] + trimMax[dim] - period) * .5;
    }
    hi[dim] = lo[dim] + period;
  }
  return BBox2D(lo, hi);
}

// Number of native periods needed to cover an unwrapped chart.
int periodCopies(double span, double period, int poles, int knots)
{
  const double copies = std::max(1., std::ceil((span - Precision::PConfusion()) / period));
  if(copies > std::numeric_limits<int>::max() / std::max(poles, knots))
  {
    throw Standard_ConstructionError("Periodic STEP domain is too large to expand");
  }
  return static_cast<int>(copies);
}

/// Open Cascade stores distinct knots separately from their multiplicities.
struct SplineKnots
{
  TColStd_Array1OfReal values;
  TColStd_Array1OfInteger multiplicities;
};

SplineKnots repeatKnots(const TColStd_Array1OfReal& knots,
                        const TColStd_Array1OfInteger& mults,
                        int copies)
{
  const int size = knots.Length() - 1;
  const double period = knots.Last() - knots.First();
  TColStd_Array1OfReal repeated(1, copies * size + 1);
  TColStd_Array1OfInteger multiplicities(1, repeated.Length());
  for(int n = 0; n < copies; ++n)
  {
    for(int i = 1; i <= size; ++i)
    {
      repeated(n * size + i) = knots(i) + n * period;
      multiplicities(n * size + i) = mults(i);
    }
  }
  repeated.ChangeLast() = knots.Last() + (copies - 1) * period;
  multiplicities.ChangeLast() = mults.Last();
  return {repeated, multiplicities};
}

/// Helper class to convert the face geometry of CAD mesh to valid NURBSPatch instances
/// The constructor converts the surface to a clamped (non-periodic) representation, if necessary
class SurfaceConverter
{
public:
  SurfaceConverter(const opencascade::handle<Geom_BSplineSurface>& surface,
                   const BBox2D& bounds,
                   int patchIndex,
                   bool verbose = false)
    : m_surface(opencascade::handle<Geom_BSplineSurface>::DownCast(surface->Copy()))
    , m_verbose(verbose)
  {
    try
    {
      ensureClamped(bounds);
    }
    catch(const Standard_Failure&)
    {
      SLIC_WARNING_ROOT(axom::fmt::format(
        "[Patch {}] Could not clamp the periodic surface to the requested bounds; "
        "retrying its native domain. Trims may extend outside the patch.",
        patchIndex));
      m_surface = opencascade::handle<Geom_BSplineSurface>::DownCast(surface->Copy());
      ensureClamped(surfaceBounds(m_surface));
    }
  }

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
    return isRational ? NPatch(extractControlPoints(),
                               extractWeights(),
                               copyKnotSequence(m_surface->UKnotSequence()),
                               copyKnotSequence(m_surface->VKnotSequence()))
                      : NPatch(extractControlPoints(),
                               copyKnotSequence(m_surface->UKnotSequence()),
                               copyKnotSequence(m_surface->VKnotSequence()));
  }

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

    axom::Array<double> params_u(numSamples);
    axom::numerics::linspace(m_surface->UKnot(1),
                             m_surface->UKnot(m_surface->NbUKnots()),
                             params_u.data(),
                             numSamples);
    axom::Array<double> params_v(numSamples);
    axom::numerics::linspace(m_surface->VKnot(1),
                             m_surface->VKnot(m_surface->NbVKnots()),
                             params_v.data(),
                             numSamples);

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
        SLIC_WARNING_ROOT_IF(m_verbose && sq_dist > sq_tol,
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

private:
  // Segment accepts at most one period. Repeat the poles, weights and knots
  // so it can clamp the full trim domain without changing the parameterization.
  void repeatPeriods(const BBox2D& bounds)
  {
    const auto span = bounds.range();
    const int nu = m_surface->IsUPeriodic()
      ? periodCopies(span[0], m_surface->UPeriod(), m_surface->NbUPoles(), m_surface->NbUKnots())
      : 1;
    const int nv = m_surface->IsVPeriodic()
      ? periodCopies(span[1], m_surface->VPeriod(), m_surface->NbVPoles(), m_surface->NbVKnots())
      : 1;
    if(nu == 1 && nv == 1)
    {
      return;
    }

    const auto [uknots, umults] = repeatKnots(m_surface->UKnots(), m_surface->UMultiplicities(), nu);
    const auto [vknots, vmults] = repeatKnots(m_surface->VKnots(), m_surface->VMultiplicities(), nv);
    TColgp_Array2OfPnt poles(1, nu * m_surface->NbUPoles(), 1, nv * m_surface->NbVPoles());
    TColStd_Array2OfReal weights(1, poles.ColLength(), 1, poles.RowLength());
    for(int i = 1; i <= poles.ColLength(); ++i)
    {
      for(int j = 1; j <= poles.RowLength(); ++j)
      {
        const int u = (i - 1) % m_surface->NbUPoles() + 1;
        const int v = (j - 1) % m_surface->NbVPoles() + 1;
        poles(i, j) = m_surface->Pole(u, v);
        weights(i, j) = m_surface->Weight(u, v);
      }
    }
    m_surface = new Geom_BSplineSurface(poles,
                                        weights,
                                        uknots,
                                        vknots,
                                        umults,
                                        vmults,
                                        m_surface->UDegree(),
                                        m_surface->VDegree(),
                                        m_surface->IsUPeriodic(),
                                        m_surface->IsVPeriodic());
  }

  /// converts the surface from periodic knots to clamped knots, when necessary
  void ensureClamped(const BBox2D& bounds)
  {
    if(!m_surface->IsUPeriodic() && !m_surface->IsVPeriodic())
    {
      return;
    }
    const auto& lo = bounds.getMin();
    const auto& hi = bounds.getMax();
    repeatPeriods(bounds);
    m_surface->Segment(lo[0], hi[0], lo[1], hi[1]);
    m_surface->SetUNotPeriodic();
    m_surface->SetVNotPeriodic();
    restoreParameterOrigin(lo);
  }

  void restoreParameterOrigin(const PointType& origin)
  {
    // Segment can shift the domain by whole periods. Shift the knots back to
    // the requested domain so they use the same UV coordinates as the trims.
    const auto shift = origin - surfaceBounds(m_surface).getMin();
    TColStd_Array1OfReal uknots(1, m_surface->NbUKnots());
    TColStd_Array1OfReal vknots(1, m_surface->NbVKnots());
    m_surface->UKnots(uknots);
    m_surface->VKnots(vknots);
    for(int i = 1; i <= uknots.Length(); ++i)
    {
      uknots(i) += shift[0];
    }
    for(int i = 1; i <= vknots.Length(); ++i)
    {
      vknots(i) += shift[1];
    }
    m_surface->SetUKnots(uknots);
    m_surface->SetVKnots(vknots);
  }

  axom::Array<PointType3D, 2> extractControlPoints() const
  {
    axom::Array<PointType3D, 2> points(axom::ArrayOptions::Uninitialized {},
                                       m_surface->NbUPoles(),
                                       m_surface->NbVPoles());
    for(int i = 0; i < m_surface->NbUPoles(); ++i)
    {
      for(int j = 0; j < m_surface->NbVPoles(); ++j)
      {
        const auto& pole = m_surface->Pole(i + 1, j + 1);
        points(i, j) = PointType3D {pole.X(), pole.Y(), pole.Z()};
      }
    }
    return points;
  }

  axom::Array<double, 2> extractWeights() const
  {
    axom::Array<double, 2> weights(axom::ArrayOptions::Uninitialized {},
                                   m_surface->NbUPoles(),
                                   m_surface->NbVPoles());
    for(int i = 0; i < m_surface->NbUPoles(); ++i)
    {
      for(int j = 0; j < m_surface->NbVPoles(); ++j)
      {
        weights(i, j) = m_surface->Weight(i + 1, j + 1);
      }
    }
    return weights;
  }

private:
  opencascade::handle<Geom_BSplineSurface> m_surface;

  bool m_verbose {false};
};

/// Helper class to convert trimming curves to valid NURBSCurve instances
/// The constructor converts the curve to a clamped (non-periodic) representation, if necessary
class PcurveConverter
{
public:
  PcurveConverter(const opencascade::handle<Geom2d_BSplineCurve>& curve, bool verbose = false)
    : m_curve(opencascade::handle<Geom2d_BSplineCurve>::DownCast(curve->Copy()))
    , m_verbose(verbose)
  {
    ensureClamped();
  }

  /// Returns a representation of the trimming curve as an axom::primal::NURBSCurve
  NCurve nurbsCurve() const
  {
    const bool isPeriodic = m_curve->IsPeriodic();
    SLIC_ERROR_IF(isPeriodic, "Axom's NURBSCurve only supports non-periodic curves");

    return m_curve->IsRational()
      ? NCurve(extractControlPoints(), extractWeights(), copyKnotSequence(m_curve->KnotSequence()))
      : NCurve(extractControlPoints(), copyKnotSequence(m_curve->KnotSequence()));
  }

  /**
   * Utility function to compare the represented curve to \a otherCurve
   * at \a numSamples uniformly sampled points in parameter space
   * Returns true when the sum of distances is less than the provided tolerance
   */
  bool compareToCurve(const opencascade::handle<Geom2d_BSplineCurve>& otherCurve,
                      int numSamples,
                      double sq_tol = 1e-8) const
  {
    if(numSamples <= 1)
    {
      return true;
    }

    axom::Array<double> params(numSamples);
    axom::numerics::linspace(m_curve->Knot(1),
                             m_curve->Knot(m_curve->NbKnots()),
                             params.data(),
                             numSamples);

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
      SLIC_WARNING_ROOT_IF(
        m_verbose && sq_dist > sq_tol,
        axom::fmt::format(
          "Distance between curves at evaluated param {} exceeded tolerance {}.\n"
          "Point on my curve {}; Point on other curve {}; Squared distance {}  (running sum {})",
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

  axom::Array<PointType> extractControlPoints() const
  {
    axom::Array<PointType> points(axom::ArrayOptions::Uninitialized {}, m_curve->NbPoles());
    for(int i = 0; i < m_curve->NbPoles(); ++i)
    {
      const auto& pole = m_curve->Pole(i + 1);
      points[i] = PointType {pole.X(), pole.Y()};
    }
    return points;
  }

  axom::Array<double> extractWeights() const
  {
    axom::Array<double> weights(axom::ArrayOptions::Uninitialized {}, m_curve->NbPoles());
    for(int i = 0; i < m_curve->NbPoles(); ++i)
    {
      weights[i] = m_curve->Weight(i + 1);
    }
    return weights;
  }

private:
  opencascade::handle<Geom2d_BSplineCurve> m_curve;
  bool m_verbose {false};
};

}  // namespace

/**
 * Class to read in a STEP file representing trimmed NURBS meshes using Open Cascade
 * and convert the patches and trimming curves to Axom's NURBSPatch and NURBSCurve primitives.
 *
 * Implementation note: Since Axom's primitives do not support periodic knots,
 * we must convert the Open Cascade analogues to an open/clamped representation, when necessary.
 */
class StepFileProcessor
{
private:
  /// Add the prepared trims and apply the source face's orientation.
  void addTrimmingCurves(NPatch& patch,
                         PatchData& patchData,
                         axom::ArrayView<const PreparedPcurve> curves)
  {
    const int patchIndex = patchData.patchIndex;
    SLIC_INFO_ROOT_IF(m_verbose,
                      axom::fmt::format("[Patch {}]: BBox in parametric space: {}",
                                        patchIndex,
                                        patchData.parametricBBox));
    for(const auto& prepared : curves)
    {
      PcurveConverter processor(prepared.curve, m_verbose);
      auto curve = processor.nurbsCurve();
      SLIC_ASSERT(curve.isValidNURBS());
      SLIC_ASSERT(curve.getDegree() == prepared.curve->Degree());
      if(prepared.reversed)
      {
        curve.reverseOrientation();
      }
      patch.addTrimmingCurve(curve);
      patchData.trimmingCurves_originallyPeriodic.push_back(prepared.originallyPeriodic);
      patchData.trimmingCurves_wireIds.push_back(prepared.wireIndex);
      SLIC_INFO_ROOT_IF(m_verbose,
                        axom::fmt::format("[Patch {} Wire {} Edge {}] Added curve: {}",
                                          patchIndex,
                                          prepared.wireIndex,
                                          prepared.edgeIndex,
                                          curve));
      if(prepared.originallyPeriodic)
      {
        const bool withinThreshold = processor.compareToCurve(prepared.curve, 25);
        SLIC_WARNING_ROOT_IF(!withinThreshold,
                             axom::fmt::format("[Patch {} Wire {} Edge {}] Trimming curve was not "
                                               "within threshold after clamping.",
                                               patchIndex,
                                               prepared.wireIndex,
                                               prepared.edgeIndex));
      }
    }
    if(patchData.was_reversed_u)
    {
      // Reverse the surface, UV coordinates, and trim directions together.
      patch.reverseOrientation_u();
    }
  }

public:
  StepFileProcessor(const std::string& filename, bool verbose = false) : m_verbose(verbose)
  {
    m_shape = loadStepFile(filename);
  }

  bool isLoaded() const { return !m_shape.IsNull(); }

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

  /// Converts each source face and its trims to a NURBSPatch.
  void extractPatches(NPatchArray& patches)
  {
    patches.resize(m_source_faces.size());

    for(int patchIndex = 0; patchIndex < m_source_faces.size(); ++patchIndex)
    {
      SLIC_DEBUG_ROOT_IF(m_verbose, "*** Processing patch " << patchIndex);
      const auto& source_face = m_source_faces[patchIndex];
      const auto face = TopoDS::Face(m_conversion->convertedShape(source_face));

      opencascade::handle<Geom_Surface> surface = BRep_Tool::Surface(face);
      if(surface->IsKind(STANDARD_TYPE(Geom_BSplineSurface)))
      {
        opencascade::handle<Geom_BSplineSurface> bsplineSurface =
          opencascade::handle<Geom_BSplineSurface>::DownCast(surface);

        axom::Array<PreparedPcurve> curves;
        std::optional<SurfaceConverter> patchProcessor;
        try
        {
          curves = preparePcurves(patchIndex, source_face, face, bsplineSurface, *m_conversion);
          const auto bounds = clampingBounds(bsplineSurface, curves.view(), patchIndex);
          patchProcessor.emplace(bsplineSurface, bounds, patchIndex, m_verbose);
        }
        catch(const Standard_Failure& failure)
        {
          SLIC_WARNING_ROOT(
            axom::fmt::format("Skipping patch {}: could not convert its surface ({})",
                              patchIndex,
                              failure.GetMessageString()));
          continue;
        }

        patches[patchIndex] = patchProcessor->nurbsPatchGeometry();

        PatchData& patchData = m_patchData[patchIndex];
        patchData.patchIndex = patchIndex;
        patchData.wasOriginallyPeriodic_u = bsplineSurface->IsUPeriodic();
        patchData.wasOriginallyPeriodic_v = bsplineSurface->IsVPeriodic();
        // History maps may return a canonical orientation rather than this occurrence's.
        patchData.was_reversed_u = (source_face.Orientation() != TopAbs_FORWARD);
        patchData.parametricBBox =
          BBox2D(PointType {patches[patchIndex].getMinKnot_u(), patches[patchIndex].getMinKnot_v()},
                 PointType {patches[patchIndex].getMaxKnot_u(), patches[patchIndex].getMaxKnot_v()});
        patchData.physicalBBox = patches[patchIndex].boundingBox();

        if(patchData.wasOriginallyPeriodic_u || patchData.wasOriginallyPeriodic_v)
        {
          const bool withinThreshold = patchProcessor->compareToSurface(bsplineSurface, 25);
          SLIC_WARNING_ROOT_IF(!withinThreshold,
                               axom::fmt::format("[Patch {}] Patch geometry was not within "
                                                 "threshold after clamping.\n Patch data: {}.",
                                                 patchIndex,
                                                 patches[patchIndex]));
        }
        addTrimmingCurves(patches[patchIndex], patchData, curves.view());
      }
      else
      {
        const std::string surfaceType = surface->DynamicType()->Name();
        SLIC_WARNING_ROOT(fmt::format("Skipping patch {} with non-BSpline surface type: '{}'",
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
        SLIC_INFO_ROOT(prefix << "\n" << sstr.str());
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

  /// Check the exported trims in UV space
  void validatePatchCharts(const NPatchArray& patches) const
  {
    // Open Cascade's tolerances are in physical space. A small 3D gap can still be large in UV.
    constexpr double uvTolerance = 1e-7;
    constexpr int numSamples = 1001;
    constexpr int samplesPerSpan = 21;
    axom::Array<double> parameters(numSamples);
    for(const auto& entry : m_patchData)
    {
      const int patchIndex = entry.first;
      const auto& data = entry.second;
      const auto& curves = patches[patchIndex].getTrimmingCurves();
      auto bounds = data.parametricBBox;
      bounds.expand(uvTolerance);

      // A control point outside the domain doesn't necessarily put the curve outside it.
      // Sample the curve itself to avoid false positives.
      auto intervalIsInside = [&](const NCurve& curve, double first, double last, int count) {
        axom::numerics::linspace(first, last, parameters.data(), count);
        for(int i = 0; i < count; ++i)
        {
          if(!bounds.contains(curve.evaluate(parameters[i])))
          {
            return false;
          }
        }
        return true;
      };

      bool outside = false;
      int wireBegin = 0;
      // Extraction appends each wire's trimming curves contiguously.
      // Close each group separately so holes are not joined to the outer boundary.
      while(wireBegin < curves.size())
      {
        const int wireId = data.trimmingCurves_wireIds[wireBegin];
        int wireEnd = wireBegin + 1;
        while(wireEnd < curves.size() && data.trimmingCurves_wireIds[wireEnd] == wireId)
        {
          ++wireEnd;
        }

        double maxSquaredGap = 0.;
        for(int i = wireBegin; i < wireEnd; ++i)
        {
          const auto& curve = curves[i];
          outside |= !intervalIsInside(curve, curve.getMinKnot(), curve.getMaxKnot(), numSamples);
          // Include short knot spans that the whole-curve samples can miss.
          for(int k = 0; !outside && k + 1 < curve.getNumKnots(); ++k)
          {
            const double first = curve.getKnot(k), last = curve.getKnot(k + 1);
            if(last > first)
            {
              outside |= !intervalIsInside(curve, first, last, samplesPerSpan);
            }
          }

          const auto& next = curves[i + 1 == wireEnd ? wireBegin : i + 1];
          // Face reversal reverses each trim while preserving their array order.
          const auto end =
            curve.evaluate(data.was_reversed_u ? curve.getMinKnot() : curve.getMaxKnot());
          const auto start =
            next.evaluate(data.was_reversed_u ? next.getMaxKnot() : next.getMinKnot());
          maxSquaredGap = std::max(maxSquaredGap, axom::primal::squared_distance(end, start));
        }
        SLIC_WARNING_ROOT_IF(
          maxSquaredGap > uvTolerance * uvTolerance,
          axom::fmt::format(
            "[Patch {} Wire {}] Exported trimming wire has a UV endpoint gap of {}.",
            patchIndex,
            wireId,
            std::sqrt(maxSquaredGap)));
        wireBegin = wireEnd;
      }
      SLIC_WARNING_ROOT_IF(
        outside,
        axom::fmt::format("[Patch {}] Sampled trimming curves leave the exported UV domain.",
                          patchIndex));
    }
  }

  const PatchDataMap& getPatchDataMap() const { return m_patchData; }

  std::string getFileUnits() const { return m_fileUnits; }

private:
  /// Loads the step file \a filename from disk
  /// Uses the units from \a filename
  TopoDS_Shape loadStepFile(const std::string& filename)
  {
    STEPControl_Reader reader;

    IFSelect_ReturnStatus status = reader.ReadFile(filename.c_str());
    if(status != IFSelect_RetDone)
    {
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
      const auto fileUnits = axom::utilities::getLengthUnit(anUnitLengthNames(1).ToCString());
      m_fileUnits = axom::utilities::getLengthUnitName(fileUnits);
      const auto defaultUnit =
        axom::utilities::getLengthUnit(Interface_Static::CVal("xstep.cascade.unit"));
      const double lengthUnit = axom::utilities::getConversionFactor(fileUnits, defaultUnit);
      reader.SetSystemLengthUnit(lengthUnit);
    }

    Standard_Integer numRoots = reader.NbRootsForTransfer();
    reader.TransferRoots();
    TopoDS_Shape shape = reader.OneShape();
    if(shape.IsNull())
    {
      std::cerr << "Error: No shape found in the file." << std::endl;
      return TopoDS_Shape();
    }

    // Keep the source faces for traversal. Convert a copy to preserve their geometry.
    for(TopExp_Explorer faceExp(shape, TopAbs_FACE); faceExp.More(); faceExp.Next())
    {
      m_source_faces.push_back(TopoDS::Face(faceExp.Current()));
    }
    TopoDS_Shape nurbsShape;
    try
    {
      m_conversion = std::make_unique<NurbsConversion>(shape);
      nurbsShape = m_conversion->shape();
    }
    catch(const Standard_Failure& failure)
    {
      std::cerr << "Error: Conversion to NURBS failed: " << failure.GetMessageString() << std::endl;
      return TopoDS_Shape();
    }

    if(nurbsShape.IsNull())
    {
      std::cerr << "Error: Conversion to NURBS failed." << std::endl;
      return TopoDS_Shape();
    }

    SLIC_INFO_ROOT_IF(m_verbose,
                      axom::fmt::format("Successfully read the STEP file with {} roots", numRoots));

    return nurbsShape;
  }

private:
  TopoDS_Shape m_shape;
  bool m_verbose {false};

  std::string m_fileUnits {"mm"};

  PatchDataMap m_patchData;
  axom::Array<TopoDS_Face> m_source_faces;
  std::unique_ptr<NurbsConversion> m_conversion;
};

/// Triangulate the BRep or its untrimmed supports into a Mint mesh with source patch IDs.
class PatchTriangulator
{
public:
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

  void triangulatePatches(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& output_mesh,
                          bool trimmed)
  {
    axom::Array<int> patch_id;

    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More(); faceExp.Next(), ++patchIndex)
    {
      TopoDS_Face face = TopoDS::Face(faceExp.Current());
      const bool isReversed = (face.Orientation() != TopAbs_FORWARD);

      if(!trimmed)
      {
        const auto surface = BRep_Tool::Surface(face);
        double u1, u2, v1, v2;
        surface->Bounds(u1, u2, v1, v2);
        const opencascade::handle<Geom_RectangularTrimmedSurface> untrimmed_surface =
          new Geom_RectangularTrimmedSurface(surface, u1, u2, v1, v2);
        const auto orientation = face.Orientation();
        face = BRepBuilderAPI_MakeFace(untrimmed_surface, Precision::Confusion());
        face.Orientation(orientation);
        BRepMesh_IncrementalMesh mesh(face, m_deflection, m_deflectionIsRelative, m_angularDeflection);
      }

      TopLoc_Location loc;
      opencascade::handle<Poly_Triangulation> triangulation = BRep_Tool::Triangulation(face, loc);

      if(triangulation.IsNull())
      {
        SLIC_WARNING_ROOT(
          axom::fmt::format("Error: Triangulation could not be generated for {}patch {}",
                            trimmed ? "" : "untrimmed ",
                            patchIndex));
        if(trimmed)
        {
          continue;
        }
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

        // Ensure triangle orientation matches the oriented face normal.
        if(isReversed)
        {
          axom::utilities::swap(p2, p3);
        }

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

}  // namespace internal

std::string STEPReader::getFileUnits() const
{
  return m_stepProcessor && m_stepProcessor->isLoaded() ? m_stepProcessor->getFileUnits() : "";
}

axom::ArrayView<const int> STEPReader::getTrimmingCurveWireIds(int patchArrayIndex) const
{
  if(patchArrayIndex < 0 || patchArrayIndex >= static_cast<int>(m_trimmingCurveWireIds.size()))
  {
    return {};
  }

  return m_trimmingCurveWireIds[patchArrayIndex].view();
}

std::string STEPReader::getBRepStats() const
{
  // early return if the step file has not been loaded
  // (or, e.g. for the derived PSTEPReader class, if it's not loaded on this rank)
  if(!m_stepProcessor || !m_stepProcessor->isLoaded())
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
    if(data.empty())
    {
      return {};
    }

    const auto [min, max] = std::ranges::minmax(data);
    const double sum = std::accumulate(data.begin(), data.end(), 0.0);
    const double sumSquared =
      std::accumulate(data.begin(), data.end(), 0.0, [](double a, double b) { return a + b * b; });
    const double mean = sum / data.size();
    return {min, max, mean, std::sqrt(sumSquared / data.size() - mean * mean)};
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
    internal::BBox3D meshBBox;
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

axom::primal::BoundingBox<double, 3> STEPReader::getBRepBoundingBox(bool useTriangulation) const
{
  if(!m_stepProcessor || !m_stepProcessor->isLoaded())
  {
    SLIC_WARNING_ROOT("Cannot compute bounding box until calling STEPReader::read()");
    return axom::primal::BoundingBox<double, 3> {};
  }

  const auto shape = m_stepProcessor->getShape();

  Bnd_Box box;
  box.SetVoid();

  BRepBndLib::Add(shape, box, useTriangulation);
  box.SetGap(0.0);

  Standard_Real xmin, ymin, zmin, xmax, ymax, zmax;
  box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

  return axom::primal::BoundingBox<double, 3> {axom::primal::Point<double, 3> {xmin, ymin, zmin},
                                               axom::primal::Point<double, 3> {xmax, ymax, zmax}};
}

// Define these where StepFileProcessor is complete so unique_ptr can destroy it.
STEPReader::STEPReader() = default;
STEPReader::~STEPReader() = default;

void STEPReader::clear()
{
  m_stepProcessor.reset();
  m_patches.clear();
  m_patchIds.clear();
  m_trimmingCurveWireIds.clear();
  m_patchOriginallyPeriodic_u.clear();
  m_patchOriginallyPeriodic_v.clear();
}

bool STEPReader::patchWasOriginallyPeriodic_u(int patchArrayIndex) const
{
  SLIC_ASSERT(patchArrayIndex >= 0 && patchArrayIndex < m_patchOriginallyPeriodic_u.size());
  return m_patchOriginallyPeriodic_u[patchArrayIndex] != 0;
}

bool STEPReader::patchWasOriginallyPeriodic_v(int patchArrayIndex) const
{
  SLIC_ASSERT(patchArrayIndex >= 0 && patchArrayIndex < m_patchOriginallyPeriodic_v.size());
  return m_patchOriginallyPeriodic_v[patchArrayIndex] != 0;
}

int STEPReader::read(bool validate_model)
{
  clear();
  m_stepProcessor = std::make_unique<internal::StepFileProcessor>(m_fileName, m_verbosity);
  if(!m_stepProcessor->isLoaded())
  {
    return 1;
  }

  if(validate_model)
  {
    m_stepProcessor->validateBRep();
  }

  m_stepProcessor->extractPatches(m_patches);
  if(validate_model)
  {
    // Check the exported curves after chart shifts and face reversal.
    m_stepProcessor->validatePatchCharts(m_patches);
  }

  // Record stable ids that match the input STEP enumeration,
  // even if consumers later filter/skip patches or trimming curves.
  m_patchIds.resize(m_patches.size());
  m_trimmingCurveWireIds.resize(m_patches.size());
  m_patchOriginallyPeriodic_u.resize(m_patches.size());
  m_patchOriginallyPeriodic_u.fill(0);
  m_patchOriginallyPeriodic_v.resize(m_patches.size());
  m_patchOriginallyPeriodic_v.fill(0);

  const auto& patchDataMap = m_stepProcessor->getPatchDataMap();
  for(int patchArrayIndex = 0; patchArrayIndex < m_patches.size(); ++patchArrayIndex)
  {
    // Current implementation preserves patch array index == input face index.
    // Keep this explicit mapping in case future logic filters patches.
    m_patchIds[patchArrayIndex] = patchArrayIndex;

    auto it = patchDataMap.find(patchArrayIndex);
    if(it != patchDataMap.end())
    {
      m_trimmingCurveWireIds[patchArrayIndex] = it->second.trimmingCurves_wireIds;
      m_patchOriginallyPeriodic_u[patchArrayIndex] = it->second.wasOriginallyPeriodic_u ? 1 : 0;
      m_patchOriginallyPeriodic_v[patchArrayIndex] = it->second.wasOriginallyPeriodic_v ? 1 : 0;
    }
  }

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
    SLIC_WARNING_ROOT("Cannot triangulate model until calling STEPReader::read()");
    return 1;
  }

  if(!mesh)
  {
    SLIC_WARNING_ROOT("Passed in mesh instance was null. Skipping triangulation");
    return 1;
  }

  internal::PatchTriangulator patchTriangulator(m_stepProcessor->getShape(),
                                                linear_deflection,
                                                angular_deflection,
                                                is_relative);

  patchTriangulator.triangulatePatches(*mesh, trimmed);

  return 0;
}

}  // namespace axom::quest
