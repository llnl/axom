// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_curved_poly_impl.hpp
 *
 * \brief Consists of functions to test intersection among geometric primitives.
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

namespace axom
{
namespace primal
{
namespace detail
{
template <typename T>
class IntersectionInfo;

template <typename T, int NDIMS>
bool orient(const BezierCurve<T, NDIMS> c1,
            const BezierCurve<T, NDIMS> c2,
            T s,
            T t);

template <typename T>
int isContained(const CurvedPolygon<T, 2> p1,
                const CurvedPolygon<T, 2> p2,
                double sq_tol = 1e-10);

/*!
 * \brief Test whether CurvedPolygons p1 and p2 intersect and find intersection points
 * \return status true iff p1 intersects with p2, otherwise false.
 *
 * \param p1, p2 CurvedPolygon objects to intersect
 * \param pnew vector of type CurvedPolygon holding intersection regions oriented as the original curves were. 
 */
template <typename T, int NDIMS>
bool intersect_polygon(CurvedPolygon<T, NDIMS>& p1,
                       CurvedPolygon<T, NDIMS>& p2,
                       std::vector<CurvedPolygon<T, NDIMS>>& pnew,
                       double sq_tol)
{
  // Object to store intersections
  std::vector<std::vector<IntersectionInfo<T>>> E1IntData(p1.numEdges());
  std::vector<std::vector<IntersectionInfo<T>>> E2IntData(p2.numEdges());
  IntersectionInfo<T> firstinter;  // Need to do orientation test on first intersection

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
                              1.,
                              0.,
                              1.,
                              0.);
      for(int k = 0; k < static_cast<int>(p1times.size()); ++k)
      {
        E1IntData[i].push_back({p1times[k], i, p2times[k], j, numinters + k + 1});
        E2IntData[j].push_back({p2times[k], j, p1times[k], i, numinters + k + 1});
        if(numinters == 0)
        {
          firstinter = {p1times[0], i, p2times[0], j, 1};
        }
      }
      numinters += p1times.size();
    }
  }
  if(numinters > 0)
  {
    for(int i = 0; i < p1.numEdges(); ++i)
    {
      std::sort(E1IntData[i].begin(), E1IntData[i].end());
      std::sort(E2IntData[i].begin(), E2IntData[i].end());
    }

    // Orient the first intersection point to be sure we get the intersection
    bool orientation = detail::orient(p1[firstinter.myEdge],
                                      p2[firstinter.otherEdge],
                                      firstinter.myTime,
                                      firstinter.otherTime);

    // Objects to store completely split polygons (split at every intersection point) and vector with unique id for each
    // intersection and zeros for corners of original polygons.
    std::vector<int> edgelabels[2];  // 0 for curves that end in original vertices, unique id for curves that end in intersection points
    CurvedPolygon<T, NDIMS> psplit[2];  // The two completely split polygons will be stored in this array
    psplit[0] = p1;
    psplit[1] = p2;

    //split polygon 1 at all the intersection points and store as psplit[0]
    int addedints = 0;
    for(int i = 0; i < p1.numEdges(); ++i)
    {
      edgelabels[0].push_back(0);
      for(int j = 0; j < static_cast<int>(E1IntData[i].size()); ++j)
      {
        psplit[0].splitEdge(i + addedints, E1IntData[i][j].myTime);
        edgelabels[0].insert(edgelabels[0].begin() + i + addedints,
                             E1IntData[i][j].numinter);
        addedints += 1;
        for(int k = j + 1; k < static_cast<int>(E1IntData[i].size()); ++k)
        {
          E1IntData[i][k].myTime =
            (E1IntData[i][k].myTime - E1IntData[i][j].myTime) /
            (1 - E1IntData[i][j].myTime);
        }
      }
    }

    //split polygon 2 at all the intersection points and store as psplit[1]
    addedints = 0;
    for(int i = 0; i < p2.numEdges(); ++i)
    {
      edgelabels[1].push_back(0);
      for(int j = 0; j < static_cast<int>(E2IntData[i].size()); ++j)
      {
        psplit[1].splitEdge(i + addedints, E2IntData[i][j].myTime);
        edgelabels[1].insert(edgelabels[1].begin() + i + addedints,
                             E2IntData[i][j].numinter);
        addedints += 1;
        for(int k = j + 1; k < static_cast<int>(E2IntData[i].size()); ++k)
        {
          E2IntData[i][k].myTime =
            (E2IntData[i][k].myTime - E2IntData[i][j].myTime) /
            (1 - E2IntData[i][j].myTime);
        }
      }
    }

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
  else
  {
    int containment = isContained(p1, p2, sq_tol);
    if(containment == 0)
    {
      return false;
    }
    else
    {
      if(containment == 1)
      {
        pnew.push_back(p1);
      }
      else
      {
        pnew.push_back(p2);
      }
      return true;
    }
  }
  // If there are no intersections
  return false;
}

/*! 
 * \brief Checks if two polygons are mutually exclusive or if one includes the other,
 * assuming that they have no intersection points
 *
 * \param [in] p1, p2 CurvedPolygons to be tested
 * \return 0 if mutually exclusive, 1 if p1 is in p2, 2 if p2 is in p1
 */
template <typename T>
int isContained(const CurvedPolygon<T, 2> p1,
                const CurvedPolygon<T, 2> p2,
                double sq_tol)
{
  const int NDIMS = 2;
  using PointType = primal::Point<T, NDIMS>;
  using BCurve = BezierCurve<T, NDIMS>;
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
                            1.,
                            0.,
                            1.,
                            0.);
    for(int i = 0; i < temps.size(); ++i)
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
    intersect(LineGuess, p2[j], temps, tempt);  // ?
    for(int i = 0; i < temps.size(); ++i)
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
  bool E1inE2 = Vec3::cross_product(p1[p1c].dt(p1t), LineGuess.dt(line1s))[2] < 0;
  bool E2inE1 = Vec3::cross_product(p2[p2c].dt(p2t), LineGuess.dt(line2s))[2] < 0;
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

// This determines with curve is "more" counterclockwise using the cross product of tangents
template <typename T, int NDIMS>
bool orient(const BezierCurve<T, NDIMS> c1, const BezierCurve<T, NDIMS> c2, T s, T t)
{
  const auto orientation =
    primal::Vector<T, 3>::cross_product(c1.dt(s), c2.dt(t))[2];
  return (orientation > 0);
}

// A class for storing intersection points so they can be easily sorted by parameter value
template <typename T>
class IntersectionInfo
{
public:
  T myTime;
  int myEdge;
  T otherTime;
  int otherEdge;
  int numinter;
  bool operator<(IntersectionInfo other) { return myTime < other.myTime; }
};

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECTION_CURVED_POLYGON_IMPL_HPP_
