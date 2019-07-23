// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_curved_poly.hpp
 *
 * \brief Consists of functions to test intersection among geometric primitives.
 */

#ifndef AXOM_PRIMAL_INTERSECTION_CURVED_POLYGON_HPP_
#define AXOM_PRIMAL_INTERSECTION_CURVED_POLYGON_HPP_

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
#include "axom/primal/operators/intersect.hpp"

namespace axom
{
namespace primal
{
template <typename T>
class IntersectionInfo;
/*
template <typename T, int NDIMS>
bool orient(BezierCurve<T,NDIMS> c1, BezierCurve<T,NDIMS> c2, T s, T t);
*/
/*!
 * \brief Tests if Bezier Curves c1 and c2 intersect.
 * \return status true iff c1 intersects with c2, otherwise false.
 *
 * \param c1, c2 BezierCurve objects to intersect
 * \param sp, tp vector of type T parameter space intersection points (t-values
 * and s-values) for c1 and c2, respectively
 */
template <typename T, int NDIMS>
bool intersect_polygon(CurvedPolygon<T, NDIMS>& p1,
                       CurvedPolygon<T, NDIMS>& p2,
                       std::vector<CurvedPolygon<T, NDIMS>>& pnew)
{
  // Object to store intersections
  std::vector<std::vector<IntersectionInfo<T>>> E1IntData(p1.numEdges());
  std::vector<std::vector<IntersectionInfo<T>>> E2IntData(p2.numEdges());

  IntersectionInfo<T> firstinter;

  // Find all intersections and store
  int numinters = 0;
  for(int i = 0; i < p1.numEdges(); ++i)
  {
    for(int j = 0; j < p2.numEdges(); ++j)
    {
      std::vector<T> p1times;
      std::vector<T> p2times;
      intersect(p1[i], p2[j], p1times, p2times);
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
  for(int i = 0; i < p1.numEdges(); ++i)
  {
    std::sort(E1IntData[i].begin(), E1IntData[i].end());
    std::sort(E2IntData[i].begin(), E2IntData[i].end());
  }

  bool orientation = orient(p1[firstinter.myEdge],
                            p2[firstinter.otherEdge],
                            firstinter.myTime,
                            firstinter.otherTime);
  std::vector<int> edgelabels[2];

  CurvedPolygon<T, NDIMS> psplit[2];
  psplit[0] = p1;
  psplit[1] = p2;
  int addedints = 0;
  for(int i = 0; i < static_cast<int>(p1.numEdges()); ++i)
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
  /*std::cout << psplit[0] << std::endl;
  for (int i=0; i < edgelabels[0].size(); ++i)
  {
    std::cout << edgelabels[0][i] << std::endl;
  }*/

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
  /*std::cout << psplit[1] << std::endl;
  for (int i=0; i < edgelabels[1].size(); ++i)
  {
    std::cout << edgelabels[1][i] << std::endl;
  }*/
  std::vector<std::vector<int>::iterator> usedlabels;
  if(numinters == 0)
  {
    return false;
  }  // If there are no intersections, return false
  else
  {
    bool addingcurves = true;
    int startinter = 1;  // Start at the first intersection
    int nextinter;
    bool currentelement = orientation;
    int currentit = std::find(edgelabels[currentelement].begin(),
                              edgelabels[currentelement].end(),
                              startinter) -
      edgelabels[currentelement].begin();
    int startit = currentit;
    int nextit = (currentit + 1) % edgelabels[0].size();
    nextinter = edgelabels[currentelement][nextit];
    while(numinters > 0)
    {
      CurvedPolygon<T, NDIMS> aPart;  // To store the current intersection polygon (could be multiple)
      while(!(nextit == startit && currentelement == orientation))
      {
        if(nextit == currentit)
        {
          currentit = nextit;
          nextit = (currentit + 1) % edgelabels[0].size();
        }
        nextinter = edgelabels[currentelement][nextit];
        while(nextinter == 0)
        {
          currentit = nextit;
          if(addingcurves)
          {
            aPart.addEdge(psplit[currentelement][nextit]);
          }
          nextit = (currentit + 1) % edgelabels[0].size();
          nextinter = edgelabels[currentelement][nextit];
        }
        if(addingcurves)
        {
          aPart.addEdge(psplit[currentelement][nextit]);
          currentelement = !currentelement;
          currentit = nextit;
          nextit = std::find(edgelabels[currentelement].begin(),
                             edgelabels[currentelement].end(),
                             nextinter) -
            edgelabels[currentelement].begin();
          numinters -= 1;
        }
        else
        {
          addingcurves = true;
          currentit = nextit;
          nextit = std::find(edgelabels[currentelement].begin(),
                             edgelabels[currentelement].end(),
                             nextinter) -
            edgelabels[currentelement].begin();
        }
        std::cout << aPart << std::endl;
      }
      pnew.push_back(aPart);
    }
    addingcurves = false;
  }
  std::cout << pnew[0] << std::endl;

  return true;
}

template <typename T, int NDIMS>
bool orient(const BezierCurve<T, NDIMS> c1, const BezierCurve<T, NDIMS> c2, T s, T t)
{
  Point<T, NDIMS> c1val1 = c1.evaluate(s + (1e-13));
  Point<T, NDIMS> c1val2 = c1.evaluate(s - (1e-13));
  Point<T, NDIMS> c2val = c2.evaluate(t + (1e-13));

  return ((-(c1val1[0] - c1val2[0]) * (c2val[1] - c1val2[1])) +
          (c1val1[1] - c1val2[1]) * (c2val[0] - c1val2[0])) > 0;
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

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECTION_CURVED_POLYGON_HPP_
