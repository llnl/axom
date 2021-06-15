// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_curved_poly_impl.hpp
 *
 * \brief Consists of functions to test intersections among curved polygons
 */

#ifndef AXOM_PRIMAL_INTERSECTION_CURVED_POLYGON_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECTION_CURVED_POLYGON_IMPL_HPP_

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"

#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/operators/detail/intersect_bezier_impl.hpp"

#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
template <typename T>
class EdgeIntersectionInfo;

template <typename T>
bool orient(const BezierCurve<T, 2>& c1, const BezierCurve<T, 2>& c2, T s, T t);

template <typename T>
int isContained(const CurvedPolygon<T, 2>& p1,
                const CurvedPolygon<T, 2>& p2,
                double sq_tol = 1e-10);

/*!
 * Splits edges of each polygon based on the intersections with the other polygon
 * When the two polygons intersect, the split polygons are returned in \a psplit
 * The edges are labeled by the types of intersections in \a edgelabels
 * and \a orientation contains the relative orientation of the first intersection
 */

template <typename T, typename CurvedPolygonType>
int splitPolygonsAlongIntersections(const CurvedPolygonType& p1,
                                    const CurvedPolygonType& p2,
                                    double sq_tol,
                                    CurvedPolygonType psplit[2],
                                    std::vector<int> edgelabels[2],
                                    bool& orientation)
{
  // Object to store intersections
  std::vector<std::vector<EdgeIntersectionInfo<T>>> E1IntData(p1.numEdges());
  std::vector<std::vector<EdgeIntersectionInfo<T>>> E2IntData(p2.numEdges());
  EdgeIntersectionInfo<T> firstinter;  // Need to do orientation test on first intersection

  // Find all intersections and store
  int numinters = 0;
  for(int i = 0; i < p1.numEdges(); ++i)
  {
    for(int j = 0; j < p2.numEdges(); ++j)
    {
      std::vector<T> p1times;
      std::vector<T> p2times;
      intersect_bezier_curves(p1[i],
                              p2[j],
                              p1times,
                              p2times,
                              sq_tol,
                              p1[i].getOrder(),
                              p2[j].getOrder(),
                              0.,
                              1.,
                              0.,
                              1.);
      const int edgeIntersections = p1times.size();
      if(edgeIntersections > 0)
      {
        if(numinters == 0)
        {
          firstinter = {p1times[0], i, p2times[0], j, 1};
        }

        for(int k = 0; k < edgeIntersections; ++k, ++numinters)
        {
          E1IntData[i].push_back({p1times[k], i, p2times[k], j, numinters + 1});
          E2IntData[j].push_back({p2times[k], j, p1times[k], i, numinters + 1});

          SLIC_INFO(fmt::format(
            "Found intersection {} -- on edge {} of polygon1 at t={}"
            " and edge {} of polygon2 at t={}; intersection point {}",
            numinters + 1,
            i,
            p1times[k],
            j,
            p2times[k],
            p1[i].evaluate(p1times[k])));
        }
      }
    }
  }

  if(numinters > 0)
  {
    // Orient the first intersection point to be sure we get the intersection
    orientation = detail::orient(p1[firstinter.myEdge],
                                 p2[firstinter.otherEdge],
                                 firstinter.myTime,
                                 firstinter.otherTime);

    for(int i = 0; i < p1.numEdges(); ++i)
    {
      std::sort(E1IntData[i].begin(), E1IntData[i].end());
    }
    for(int i = 0; i < p2.numEdges(); ++i)
    {
      std::sort(E2IntData[i].begin(), E2IntData[i].end());
    }

    psplit[0] = p1;
    psplit[1] = p2;

    edgelabels[0].reserve(p1.numEdges() + numinters);
    edgelabels[1].reserve(p2.numEdges() + numinters);

    SLIC_INFO("Poly1 before split: " << psplit[0]);
    splitPolygon(psplit[0], E1IntData, edgelabels[0]);
    SLIC_INFO("Poly1 after split: " << psplit[0]);

    SLIC_INFO("Poly2 before split: " << psplit[1]);
    splitPolygon(psplit[1], E2IntData, edgelabels[1]);
    SLIC_INFO("Poly2 after split: " << psplit[1]);
  }

  return numinters;
}

/*!
 * \brief Test whether the regions bounded by CurvedPolygons \a p1 and \a p2 intersect
 * \return status true iff \a p1 intersects with \a p2, otherwise false.
 *
 * \param [in] p1, p2 CurvedPolygon objects to intersect
 * \param [in] sq_tol tolerance parameter for the base case of intersect_bezier_curve
 * \param [out] pnew vector of type CurvedPolygon holding CurvedPolygon objects representing boundaries of intersection regions. 
 */
template <typename T, int NDIMS>
bool intersect_polygon(const CurvedPolygon<T, NDIMS>& p1,
                       const CurvedPolygon<T, NDIMS>& p2,
                       std::vector<CurvedPolygon<T, NDIMS>>& pnew,
                       double sq_tol)
{
  // Objects to store completely split polygons (split at every intersection point) and vector with unique id for each intersection and zeros for corners of original polygons.
  std::vector<int> edgelabels[2];  // 0 for curves that end in original vertices, unique id for curves that end in intersection points
  CurvedPolygon<T, NDIMS> psplit[2];  // The two completely split polygons will be stored in this array
  bool orientation {false};

  int numinters = splitPolygonsAlongIntersections<T>(p1,
                                                     p2,
                                                     sq_tol,
                                                     psplit,
                                                     edgelabels,
                                                     orientation);

  if(numinters > 0)
  {
    // This performs the directional walking method using the completely split polygons
    std::vector<std::vector<int>::iterator> usedlabels;
    // When this is false, we are walking between intersection regions
    bool addingcurves = true;
    int startvertex = 1;  // Start at the vertex with "unique id" 1
    int nextvertex;       // The next vertex id is unknown at this time
    // This variable allows us to switch between the two elements
    bool currentelement = orientation;
    // This is the iterator pointing to the end vertex of the edge  of the completely split polygon we are on
    int currentit = std::find(edgelabels[currentelement].begin(),
                              edgelabels[currentelement].end(),
                              startvertex) -
      edgelabels[currentelement].begin();
    // This is the iterator to the end vertex of the starting edge on the starting polygon
    int startit = currentit;
    // This is the iterator to the end vertex of the next edge of whichever polygon we will be on next
    int nextit = (currentit + 1) % edgelabels[0].size();
    nextvertex = edgelabels[currentelement][nextit];  // This is the next vertex id
    while(numinters > 0)
    {
      CurvedPolygon<T, NDIMS> aPart;  // Object to store the current intersection polygon (could be multiple)

      // Once the end vertex of the current edge is the start vertex, we need to switch regions
      while(!(nextit == startit && currentelement == orientation) ||
            addingcurves == false)
      {
        if(nextit == currentit)
        {
          nextit = (currentit + 1) % edgelabels[0].size();
        }
        nextvertex = edgelabels[currentelement][nextit];
        while(nextvertex == 0)
        {
          currentit = nextit;
          if(addingcurves)
          {
            aPart.addEdge(psplit[currentelement][nextit]);
          }
          nextit = (currentit + 1) % edgelabels[0].size();
          nextvertex = edgelabels[currentelement][nextit];
        }
        if(edgelabels[currentelement][nextit] > 0)
        {
          if(addingcurves)
          {
            aPart.addEdge(psplit[currentelement][nextit]);
            edgelabels[currentelement][nextit] =
              -edgelabels[currentelement][nextit];
            currentelement = !currentelement;
            nextit = std::find(edgelabels[currentelement].begin(),
                               edgelabels[currentelement].end(),
                               nextvertex) -
              edgelabels[currentelement].begin();
            edgelabels[currentelement][nextit] =
              -edgelabels[currentelement][nextit];
            currentit = nextit;
            numinters -= 1;
          }
          else
          {
            addingcurves = true;
            startit = nextit;
            currentit = nextit;
            nextit = (currentit + 1) % edgelabels[0].size();
            orientation = currentelement;
          }
        }
        else
        {
          currentelement = !currentelement;
          nextit = std::find(edgelabels[currentelement].begin(),
                             edgelabels[currentelement].end(),
                             nextvertex) -
            edgelabels[currentelement].begin();
          edgelabels[currentelement][nextit] =
            -edgelabels[currentelement][nextit];
          currentit = nextit;
        }
      }
      pnew.push_back(aPart);
      currentelement = !currentelement;
      currentit = std::find(edgelabels[currentelement].begin(),
                            edgelabels[currentelement].end(),
                            -nextvertex) -
        edgelabels[currentelement].begin();
      nextit = (currentit + 1) % edgelabels[0].size();
      if(numinters > 0)
      {
        addingcurves = false;
      }
    }
    return true;
  }
  else  // If there are no intersection points, check for containment
  {
    int containment = isContained(p1, p2, sq_tol);
    switch(containment)
    {
    case 0:
      return false;
    case 1:
      pnew.push_back(p1);
      return true;
    case 2:
      pnew.push_back(p2);
      return true;
    }
    return false;  // Catch
  }
}

/*! \brief Checks if two polygons are mutually exclusive or if one includes the other, assuming that they have no intersection points
 *
 * \param [in] p1, p2 CurvedPolygons to be tested
 * \param [in] sq_tol tolerance parameter for the base case of intersect_bezier_curves
 * \return 0 if mutually exclusive, 1 if p1 is in p2, 2 if p2 is in p1
 */
template <typename T>
int isContained(const CurvedPolygon<T, 2>& p1,
                const CurvedPolygon<T, 2>& p2,
                double sq_tol)
{
  const int NDIMS = 2;
  using BCurve = BezierCurve<T, NDIMS>;
  using PointType = typename BCurve::PointType;

  int p1c = 0;
  int p2c = 0;
  T p1t = .5;
  T p2t = .5;
  PointType controlPoints[2] = {p1[p1c].evaluate(p1t), p2[p2c].evaluate(p2t)};
  BCurve LineGuess = BCurve(controlPoints, 1);
  T line1s = 0.0;
  T line2s = 1.0;
  for(int j = 0; j < p1.numEdges(); ++j)
  {
    std::vector<T> temps;
    std::vector<T> tempt;
    intersect_bezier_curves(LineGuess,
                            p1[j],
                            temps,
                            tempt,
                            sq_tol,
                            1,
                            p1[j].getOrder(),
                            0.,
                            1.,
                            0.,
                            1.);

    for(int i = 0; i < static_cast<int>(temps.size()); ++i)
    {
      if(temps[i] > line1s)
      {
        line1s = temps[i];
        p1c = j;
        p1t = tempt[i];
      }
    }
  }
  for(int j = 0; j < p2.numEdges(); ++j)
  {
    std::vector<T> temps;
    std::vector<T> tempt;
    intersect_bezier_curves(LineGuess,
                            p2[j],
                            temps,
                            tempt,
                            sq_tol,
                            1,
                            p2[j].getOrder(),
                            1.,
                            0.,
                            1.,
                            0.);
    intersect(LineGuess, p2[j], temps, tempt);
    for(int i = 0; i < static_cast<int>(temps.size()); ++i)
    {
      if(temps[i] < line2s && temps[i] > line1s)
      {
        line2s = temps[i];
        p2c = j;
        p2t = tempt[i];
      }
    }
  }

  using Vec3 = primal::Vector<T, 3>;
  const bool E1inE2 =
    Vec3::cross_product(p1[p1c].dt(p1t), LineGuess.dt(line1s))[2] < 0;
  const bool E2inE1 =
    Vec3::cross_product(p2[p2c].dt(p2t), LineGuess.dt(line2s))[2] < 0;
  if(E1inE2 && E2inE1)
  {
    return 1;
  }
  else if(!E1inE2 && !E2inE1)
  {
    return 2;
  }
  else
  {
    return 0;
  }
}

/// Splits a CurvedPolygon \a p1 at every intersection point stored in \a IntersectionData
/// \a edgeLabels stores intersection id for new vertices and 0 for original vertices
template <typename T, int NDIMS>
void splitPolygon(CurvedPolygon<T, NDIMS>& p1,
                  std::vector<std::vector<EdgeIntersectionInfo<T>>>& IntersectionData,
                  std::vector<int>& edgelabels)
{
  using axom::utilities::isNearlyEqual;

  int addedIntersections = 0;
  const int numEdges = p1.numEdges();
  for(int i = 0; i < numEdges; ++i)  // foreach edge
  {
    edgelabels.push_back(0);  //    mark start current vertex as 'original'
    const int nIntersect = IntersectionData[i].size();
    for(int j = 0; j < nIntersect; ++j)  //    foreach intersection on this edge
    {
      // split edge at parameter t_j
      const double t_j = IntersectionData[i][j].myTime;
      const int edgeIndex = i + addedIntersections;
      p1.splitEdge(edgeIndex, t_j);

      // update edge label
      const int label = IntersectionData[i][j].numinter;
      edgelabels.insert(edgelabels.begin() + edgeIndex, label);

      ++addedIntersections;

      // update remaining intersections on this edge; special case if already at end of curve
      if(!isNearlyEqual(1., t_j))
      {
        for(int k = j + 1; k < nIntersect; ++k)
        {
          const double t_k = IntersectionData[i][k].myTime;
          IntersectionData[i][k].myTime = (t_k - t_j) / (1.0 - t_j);
        }
      }
      else
      {
        for(int k = j + 1; k < nIntersect; ++k)
        {
          IntersectionData[i][k].myTime = 1.;
        }
      }
    }
  }
}

/*!
 * \brief Determines orientation of a bezier curve \a c1 with respect to another bezier curve \a c2, 
 *  given that they intersect at parameter values \a s and \a t, respectively
 *
 * \param [in] c1 the first bezier curve
 * \param [in] c2 the second bezier curve
 * \param [in] s the parameter value of intersection on \a c1
 * \param [in] t the parameter value of intersection on \a c2
 * \return True if \a c1's positive direction is counterclockwise from \a c2's positive direction
 */
template <typename T>
bool orient(const BezierCurve<T, 2>& c1, const BezierCurve<T, 2>& c2, T s, T t)
{
  const auto orientation =
    primal::Vector<T, 3>::cross_product(c1.dt(s), c2.dt(t))[2];
  return (orientation > 0);
}

/*! \class EdgeIntersectionInfo
 *
 * \brief For storing intersection points between edges of \a CurvedPolygon instances so they can be easily sorted by parameter value using std::sort
 */
template <typename T>
class EdgeIntersectionInfo
{
public:
  T myTime;    // parameter value of intersection on curve on first CurvePolygon
  int myEdge;  // index of curve on first CurvedPolygon
  T otherTime;  // parameter value of intersection on curve on second CurvedPolygon
  int otherEdge;  // index of curve on second CurvedPolygon
  int numinter;   // unique intersection point identifier

  /// \brief Comparison operator for sorting by parameter value
  bool operator<(const EdgeIntersectionInfo& other) const
  {
    return myTime < other.myTime;
  }
};

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECTION_CURVED_POLYGON_IMPL_HPP_
