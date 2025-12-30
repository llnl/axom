// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/LinearizeCurves.hpp"

// These macros enable some debugging utilities for linearization.
// #define AXOM_DEBUG_LINEARIZE_VERBOSE
// #define AXOM_DEBUG_WRITE_LINES

using PointType = axom::primal::Point<double, 2>;

namespace axom
{
namespace quest
{
namespace internal
{

//---------------------------------------------------------------------------
/*!
 * \brief Transform a 2D point and return a 2D point.
 *
 * \param transform A 4x4 transformation matrix.
 * \param pt The input 2D point.
 *
 * \tparam InputPointType The input type needs to have a subscript operator
 *         which can be called on index 0 and 1
 * \return The transformed 2D point.
 */
template <typename InputPointType>
inline PointType transformPoint(const numerics::Matrix<double>& transform, const InputPointType& pt)
{
  // Turn the point2 into a vec4.
  double pt4[] = {pt[0], pt[1], 0., 1.};
  double newpt4[] = {0., 0., 0., 1.};

  // Transform the point.
  axom::numerics::matrix_vector_multiply(transform, pt4, newpt4);

  // Make a new 2D point from the 2D components.
  return PointType {newpt4[0], newpt4[1]};
}

//---------------------------------------------------------------------------
/*!
 * \brief This data type stores a segment, which is used in refining curves.
 */
struct Segment
{
  double u;         // The u parameter for the start of the segment.
  PointType point;  // The point at parameter u.
  double length;    // Distance from point to next segment's point.

  double next_u;          // Best place to split this segment.
  double next_length[2];  // left/right length values at next_u.
};

using SegmentArray = axom::Array<Segment>;

/*!
 * \brief Returns the index of the longest segment in the segments.
 * \param segments An array of Segments.
 * \return The index of the element with the longest next_length value.
 */
axom::IndexType nextLongest(const SegmentArray& segments)
{
  axom::IndexType maxIdx = 0;
  double max_length = 0.;
  const axom::IndexType n = segments.size();
  for(axom::IndexType i = 0; i < n; i++)
  {
    double next_length = segments[i].next_length[0] + segments[i].next_length[1];
    if(next_length > max_length)
    {
      max_length = next_length;
      maxIdx = i;
    }
  }
  return maxIdx;
}

//---------------------------------------------------------------------------
/*!
 * \brief Append a set of segments to a mint mesh, filtering out like points.
 * \param mesh The mint mesh to which we're adding points/cells.
 * \param S The segments to write.
 * \param EPS_SQ The squared distance for point matching.
 */
static void appendPoints(LinearizeCurves::SegmentMesh* mesh, SegmentArray& S, double EPS_SQ)
{
  // Check for simple vertex welding opportunities at endpoints of newly interpolated points
  {
    int numNodes = mesh->getNumberOfNodes();
    if(numNodes > 0)  // this is not the first Piece
    {
      PointType meshPt;
      // Fix start point if necessary; check against most recently added vertex in mesh
      mesh->getNode(numNodes - 1, meshPt.data());
      if(primal::squared_distance(S[0].point, meshPt) < EPS_SQ)
      {
        S[0].point = meshPt;
      }

      // Fix end point if necessary; check against 0th vertex in mesh
      const int endIdx = S.size() - 1;
      mesh->getNode(0, meshPt.data());
      if(primal::squared_distance(S[endIdx].point, meshPt) < EPS_SQ)
      {
        S[endIdx].point = meshPt;
      }
    }
    else  // This is the first, and possibly only span, check its endpoint, fix if necessary
    {
      int endIdx = S.size() - 1;
      if(primal::squared_distance(S[0].point, S[endIdx].point) < EPS_SQ)
      {
        S[endIdx].point = S[0].point;
      }
    }
  }

  // Add the new points and segments to the mesh, respecting welding checks from previous block
  {
    const int startNode = mesh->getNumberOfNodes();
    const int numNewNodes = S.size();
    mesh->reserveNodes(startNode + numNewNodes);

    for(int i = 0; i < numNewNodes; ++i)
    {
      mesh->appendNode(S[i].point[0], S[i].point[1]);
    }

    const int startCell = mesh->getNumberOfCells();
    const int numNewSegments = S.size() - 1;
    mesh->reserveCells(startCell + numNewSegments);
    for(int i = 0; i < numNewSegments; ++i)
    {
      IndexType seg[2] = {startNode + i, startNode + i + 1};
      mesh->appendCell(seg, mint::SEGMENT);
    }
  }
}

/// Compute the arc length of a curve by numerical quadrature
double computeArcLength(const LinearizeCurves::NURBSCurve& nurbs, int npts)
{
  double arcLength = 0.;
  for(const auto& bezier : nurbs.extractBezier())
  {
    arcLength +=
      primal::evaluate_line_integral(bezier, [](PointType /*x*/) -> double { return 1.; }, npts);
  }
  return arcLength;
}

//---------------------------------------------------------------------------
#ifdef AXOM_DEBUG_WRITE_LINES
/*!
 * \brief Write line segments to a VTK file for visualization.
 * \param filename The name of the file to write.
 * \param S The segments to write.
 */
static void writeLines(const std::string& filename, const SegmentArray& S)
{
  axom::fmt::memory_buffer out;
  axom::fmt::format_to(std::back_inserter(out),
                       R"(# vtk DataFile Version 4.2
vtk output
ASCII
DATASET POLYDATA
FIELD FieldData 2
CYCLE 1 1 int
1
TIME 1 1 double
1.0

)");

  // Write points
  // clang-format off
  const int npts = S.size();
  axom::fmt::format_to(std::back_inserter(out), "POINTS {} float\n", npts);
  for(int i = 0; i < npts; i += 3)
  {
    axom::fmt::format_to(std::back_inserter(out),"{}{}{}\n",
      fmt::format("{:.16f} {:.16f} 0.", S[i].point[0], S[i].point[1]),
      (i + 1) < npts ? fmt::format(" {:.16f} {:.16f} 0.", S[i + 1].point[0], S[i + 1].point[1]) : "",
      (i + 2) < npts ? fmt::format(" {:.16f} {:.16f} 0.", S[i + 2].point[0], S[i + 2].point[1]) : "");
  }
  axom::fmt::format_to(std::back_inserter(out), "\n");
  // clang-format on

  // Write ncells
  int nspans = npts - 1;
  axom::fmt::format_to(std::back_inserter(out), "LINES {} {}\n", nspans, 3 * nspans);
  for(int ispan = 0; ispan < nspans; ispan++)
  {
    axom::fmt::format_to(std::back_inserter(out), "2 {} {}\n", ispan, ispan + 1);
  }

  std::ofstream f(filename, std::ios::out | std::ios::trunc);
  f << axom::fmt::to_string(out);
}
#endif

}  // end namespace internal

//---------------------------------------------------------------------------
void LinearizeCurves::getLinearMeshUniform(LinearizeCurves::CurveArrayView curves,
                                           LinearizeCurves::SegmentMesh* mesh,
                                           int segmentsPerKnotSpan) const
{
  using axom::utilities::lerp;

  // Sanity checks
  SLIC_ERROR_IF(mesh == nullptr, "supplied mesh is null!");
  SLIC_ERROR_IF(mesh->getDimension() != 2, "C2C reader expects a 2D mesh!");
  SLIC_ERROR_IF(mesh->getCellType() != mint::SEGMENT, "C2C reader expects a segment mesh!");
  SLIC_ERROR_IF(segmentsPerKnotSpan < 1, "C2C reader: Need at least one segment per NURBs span");

  using PointsArray = axom::Array<PointType>;

  const double EPS_SQ = m_vertexWeldThreshold * m_vertexWeldThreshold;

  const double denom = static_cast<double>(segmentsPerKnotSpan);
  PointsArray pts;
  pts.reserve(segmentsPerKnotSpan + 1);
  for(const auto& nurbs : curves)
  {
    for(const auto& bezier : nurbs.extractBezier())
    {
      pts.clear();
      for(int i = 0; i <= segmentsPerKnotSpan; ++i)
      {
        pts.emplace_back(bezier.evaluate(lerp(0., 1., i / denom)));
      }

      // Check for simple vertex welding opportunities at endpoints of newly interpolated points
      {
        const int numNodes = mesh->getNumberOfNodes();
        if(numNodes > 0)  // this is not the first Piece
        {
          PointType meshPt;
          // Fix start point if necessary; check against most recently added vertex in mesh
          mesh->getNode(numNodes - 1, meshPt.data());
          if(primal::squared_distance(pts[0], meshPt) < EPS_SQ)
          {
            pts[0] = meshPt;
          }

          // Fix end point if necessary; check against 0th vertex in mesh
          const auto endIdx = pts.size() - 1;
          mesh->getNode(0, meshPt.data());
          if(primal::squared_distance(pts[endIdx], meshPt) < EPS_SQ)
          {
            pts[endIdx] = meshPt;
          }
        }
        else  // This is the first, and possibly only span, check its endpoint, fix if necessary
        {
          const auto endIdx = pts.size() - 1;
          if(primal::squared_distance(pts[0], pts[endIdx]) < EPS_SQ)
          {
            pts[endIdx] = pts[0];
          }
        }
      }

      // Add the new points and segments to the mesh, respecting welding checks from previous block
      {
        const auto startNode = mesh->getNumberOfNodes();
        const auto numNewNodes = pts.size();
        mesh->reserveNodes(startNode + numNewNodes);

        for(int i = 0; i < numNewNodes; ++i)
        {
          mesh->appendNode(pts[i][0], pts[i][1]);
        }

        const auto startCell = mesh->getNumberOfCells();
        const auto numNewSegments = pts.size() - 1;
        mesh->reserveCells(startCell + numNewSegments);
        for(int i = 0; i < numNewSegments; ++i)
        {
          IndexType seg[2] = {startNode + i, startNode + i + 1};
          mesh->appendCell(seg, mint::SEGMENT);
        }
      }
    }
  }
}

//---------------------------------------------------------------------------
void LinearizeCurves::getLinearMeshNonUniform(LinearizeCurves::CurveArrayView curves,
                                              LinearizeCurves::SegmentMesh* mesh,
                                              double percentError) const
{
  // Sanity checks
  SLIC_ERROR_IF(mesh == nullptr, "supplied mesh is null!");
  SLIC_ERROR_IF(mesh->getDimension() != 2, "C2C reader expects a 2D mesh!");
  SLIC_ERROR_IF(mesh->getCellType() != mint::SEGMENT, "C2C reader expects a segment mesh!");
  SLIC_ERROR_IF(percentError <= 0.,
                axom::fmt::format("C2C reader: percentError must be greater than 0. {} supplied.",
                                  percentError));
  SLIC_ERROR_IF(
    percentError >= 100.,
    axom::fmt::format("C2C reader: percentError must be less than 100. {} supplied.", percentError));

  /*!
   * \brief Checks curve lengths against tolerances and determines whether more
   *        refinement is needed.
   *
   * \param len The curve length that has been calculated.
   * \param maxlen The upper bound curve length that was calculated with a lot
   *               of segments.
   *
   * \note Assume that maxlen is longer than len because hi-res sampling should
   *       result in a longer arc length.
   */
  auto error_percent = [](double len, double maxlen) -> double {
    double errPct = 0.;
    if(maxlen > len)
    {
      errPct = 100. * (1. - (len / maxlen));
    }
    return errPct;
  };

  constexpr size_t INITIAL_GUESS_NPTS = 100;
  const double EPS_SQ = m_vertexWeldThreshold * m_vertexWeldThreshold;

  // Clamp the lower error bound so it does not get impractically small.
  percentError = axom::utilities::clampLower(percentError, 1.e-10);

  // clang complains about contourCount (-Wunused-but-set-variable); since we want keep it, let's mark it
  int contourCount {-1};
  AXOM_UNUSED_VAR(contourCount);

  // Iterate over the contours and linearize each of them.
  internal::SegmentArray S;
  S.reserve(INITIAL_GUESS_NPTS);
  for(const auto& nurbs : curves)
  {
    ++contourCount;

    const auto knots = nurbs.getKnots();

    // Get the contour start/end parameters.
    const double u0 = knots[0];
    const double u1 = knots[knots.getNumKnots() - 1];

    // Compute the start/end points of the whole curve.
    const PointType p0 = nurbs.evaluate(u0);
    const PointType p1 = nurbs.evaluate(u1);

    // This segment represents the whole [0,1] interval.
    internal::Segment first;
    first.u = u0;
    first.point = p0;
    first.length = sqrt(primal::squared_distance(p0, p1));
    first.next_u = (u0 + u1) / 2.;
    PointType next_point = nurbs.evaluate(first.next_u);
    first.next_length[0] = sqrt(primal::squared_distance(p0, next_point));
    first.next_length[1] = sqrt(primal::squared_distance(next_point, p1));

    // This segment just contains the end point of the interval.
    internal::Segment last;
    last.u = u1;
    last.point = p1;
    last.length = 0.;
    last.next_u = u1;
    last.next_length[0] = 0.;
    last.next_length[1] = 0.;

    // Store the initial segments.
    S.clear();
    S.push_back(first);
    S.push_back(last);

    // If there is a single span in the contour then we already added its points.
    if(knots.getNumKnotSpans() > 1)
    {
      constexpr int MAX_NUMBER_OF_SAMPLES = 2000;
      constexpr int quadrature_order = 30;
      const double hiCurveLen = internal::computeArcLength(nurbs, quadrature_order);

      // The initial curve length.
      double curveLength = first.length;

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
      // Print initial iteration.
      SLIC_INFO(fmt::format("hiCurveLen: {}", hiCurveLen));
      SLIC_INFO(fmt::format("curveLength: {}", curveLength));
      SLIC_INFO(fmt::format("percentError: {}", percentError));
#endif

      // Iterate until the difference between iterations is under the percent
      // error or we've reached the max number of samples. Each iteration, we
      // split the longest line segment so the line overall should get longer
      // until it reaches/approaches the upfront arc length. We get to put the
      // new points where they matter most early on so we should ideally get to
      // an acceptable arc length before we reach MAX_NUMBER_OF_SAMPLES.
      int iteration = 0;
      while((error_percent(curveLength, hiCurveLen) > percentError) &&
            (iteration < MAX_NUMBER_OF_SAMPLES))
      {
        // Get the index of the segment we'll split.
        size_t splitIndex = internal::nextLongest(S);
        size_t nextIndex = splitIndex + 1;

        // Partition segment.
        internal::Segment& left = S[splitIndex];
        internal::Segment right;

        // The old segment length pre-split.
        double oldSegLength = left.length;

        // The left and right segments will be the sum of the parent
        // segment's 2 pieces.
        double newSegLength = left.next_length[0] + left.next_length[1];

        // Split the left segment into left, right. We set next_u as the midpoint
        // of the subintervals. It was being set via a Newton solveu() function
        // that yielded the u value with the longest left+right line segments but
        // using the midpoint ended up with better behavior.
        right.u = left.next_u;
        right.point = nurbs.evaluate(right.u);
        right.length = left.next_length[1];

        right.next_u = (right.u + S[nextIndex].u) / 2.;
        PointType right_next = nurbs.evaluate(right.next_u);
        right.next_length[0] = sqrt(primal::squared_distance(right.point, right_next));
        right.next_length[1] = sqrt(primal::squared_distance(right_next, S[nextIndex].point));

        left.length = left.next_length[0];
        left.next_u = (left.u + right.u) / 2.;
        PointType left_next = nurbs.evaluate(left.next_u);
        left.next_length[0] = sqrt(primal::squared_distance(left.point, left_next));
        left.next_length[1] = sqrt(primal::squared_distance(left_next, right.point));

        // Insert the right segment after the left segment.
        S.insert(S.begin() + nextIndex, right);

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
        // Print initial iteration.
        SLIC_INFO(fmt::format("Iteration: {}", iteration));
        SLIC_INFO(fmt::format("\thiCurveLen: {}", hiCurveLen));
        SLIC_INFO(fmt::format("\tcurveLength: {}", curveLength));
#endif

        // Update the curve length.
        curveLength = curveLength - oldSegLength + newSegLength;

#ifdef AXOM_DEBUG_WRITE_LINES
        const auto filename = axom::fmt::format("lines{}_{:05}.vtk", contourCount, S.size());
        internal::writeLines(filename, S);
        SLIC_INFO(fmt::format("Wrote {}", filename));
#endif

        iteration++;
      }  // while(error_percent())

      SLIC_INFO(
        fmt::format("getLinearMesh: "
                    "percentError = {}"
                    ", hiCurveLength = {}"
                    ", curveLength = {}",
                    percentError,
                    hiCurveLen,
                    curveLength));
    }  // if(interpolator.numSpans() > 1)

    // Add the points to the mesh.
    internal::appendPoints(mesh, S, EPS_SQ);
  }
}

//---------------------------------------------------------------------------
double LinearizeCurves::revolvedVolume(const LinearizeCurves::NURBSCurve& nurbs,
                                       const numerics::Matrix<double>& transform) const
{
  using PointType = axom::primal::Point<double, 2>;
  using VectorType = axom::primal::Vector<double, 2>;

  // Use 5-point Gauss-Legendre Quadrature.
  numerics::QuadratureRule quad_rule = numerics::get_gauss_legendre(5);

  // Make a transform with no translation. We use this to transform
  // the derivative since we want to permit scaling and rotation but
  // translating it does not make sense.
  numerics::Matrix<double> transform2 = transform;
  transform2(0, 3) = 0.;
  transform2(1, 3) = 0.;
  transform2(2, 3) = 0.;

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
  SLIC_INFO(fmt::format("revolvedVolume"));
#endif

  // Break up the [0,1] interval and compute quadrature in the subintervals.
  double vol = 0.;
  for(const auto& bezier : nurbs.extractBezier())
  {
    constexpr double ad = 0;
    constexpr double bd = 1;
#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
    SLIC_INFO(fmt::format("interval ({}, {})", ad, bd));
#endif

    // Approximate the integral "Int pi*x'(u)*y(u)^2du" using quadrature.
    constexpr double scale = M_PI * (bd - ad);
    double sum = 0.;
    for(size_t i = 0; i < 5; i++)
    {
      // Map quad point x value [0,1] to [a,b].
      const double u = quad_rule.node(i) * (bd - ad) + ad;

      // Compute y(u) to get radius
      PointType eval;
      VectorType dt;
      bezier.evaluateFirstDerivative(u, eval, dt);

      const PointType p_uT = internal::transformPoint(transform, eval);
      const double r = p_uT[1];

      const PointType xprimeT = internal::transformPoint(transform2, dt);
      const double xp = xprimeT[0];

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
      SLIC_INFO(
        fmt::format("\ti={}, u={}, p(u)=({},{}), xp=({},{})", i, u, eval[0], eval[1], dt[0], dt[1]));
#endif

      // Accumulate weight times dx*r^2.
      sum += quad_rule.weight(i) * xp * (r * r);
    }
    // Guard against volumes being negative (if the curve went the wrong way)
    vol += axom::utilities::abs(scale * sum);
  }
#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
  SLIC_INFO(fmt::format("revolvedVolume={}", vol));
#endif
  return vol;
}

//---------------------------------------------------------------------------
double LinearizeCurves::getRevolvedVolume(const LinearizeCurves::CurveArrayView curves,
                                          const numerics::Matrix<double>& transform) const
{
  double vol = 0.;
  for(const auto& nurbs : curves)
  {
    vol += revolvedVolume(nurbs, transform);
  }
  return vol;
}

}  // end namespace quest
}  // end namespace axom
