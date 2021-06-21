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
bool orient(const BezierCurve<T, 2>& c1, const BezierCurve<T, 2>& c2, T s, T t);

template <typename T>
int isContained(const CurvedPolygon<T, 2>& p1,
                const CurvedPolygon<T, 2>& p2,
                double sq_tol = 1e-10);

template <typename T, int NDIMS>
class DirectionalWalk
{
public:
  using CurvedPolygonType = CurvedPolygon<T, NDIMS>;
  using EdgeLabels = std::vector<int>;

public:
  /*! \class EdgeIntersectionInfo
  *
  * \brief For storing intersection points between edges of \a CurvedPolygon instances so they can be easily sorted by parameter value using std::sort
  */
  class EdgeIntersectionInfo
  {
  public:
    T myTime;  // parameter value of intersection on curve on first CurvePolygon
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

  /// Enum to represent different Junction states
  enum class JunctionState : int
  {
    UNINITIALIZED = -1,
    NON_JUNCTION,  /// Not a junction, e.g. a vertex of the original
    APPLIED,       /// Use after junction has been applied
    CROSS,         /// Typical case: edges of polygons intersect
    GRAZE          /// Atypical case: polygons intersecte at a terminal vertex
  };

  /*!
   * Helper class for encoding information about the incident edges at an interection
   */
  struct Junction
  {
    using JunctionIndex = int;
    using EdgeIndex = int;
    static const EdgeIndex INVALID_EDGE_INDEX = -1;
    static const int INVALID_JUNCTION_INDEX = -1;

    // indices of polygon edges leading into junction
    EdgeIndex edgeIndex[2] {INVALID_EDGE_INDEX, INVALID_EDGE_INDEX};
    // describes the junction type and status
    JunctionState junctionState {JunctionState::UNINITIALIZED};
    JunctionIndex index {INVALID_JUNCTION_INDEX};

    bool isActive() const { return junctionState > JunctionState::APPLIED; }

    EdgeIndex currentEdgeIndex(bool active) const { return edgeIndex[active]; }
    EdgeIndex nextEdgeIndex(bool active) const { return edgeIndex[active] + 1; }

    bool operator==(const Junction& other) const
    {
      return index == other.index;
    }
    bool operator!=(const Junction& other) const { return !(*this == other); }
  };

  /*!
   * Helper class for indexing an edge of a polygon
   * Assumption is that the polygon is not empty 
   */
  class PolygonEdge
  {
  public:
    using VertexIndex = int;
    using EdgeIndex = int;

    PolygonEdge(int polygonID, const EdgeLabels& labels)
      : m_polygonID(polygonID)
      , m_labels(labels)
      , m_numEdges(labels.size())
    {
      SLIC_ASSERT_MSG(m_numEdges > 0, "Polygon " << polygonID << " was empty");
    }

    EdgeIndex getIndex() const { return m_index; }
    void setIndex(EdgeIndex idx) { m_index = (idx % m_numEdges); }

    int getStartLabel() const { return m_labels[prevIndex()]; }
    int getEndLabel() const { return m_labels[m_index]; }
    void advance() { m_index = nextIndex(); }

    bool isJunction() const { return getEndLabel() > 0; }

    bool operator==(const PolygonEdge& other)
    {
      return m_index == other.m_index && m_polygonID == other.m_polygonID;
    }

    bool operator!=(const PolygonEdge& other) { return !(*this == other); }

  private:
    EdgeIndex lastIndex() const { return m_numEdges - 1; }
    EdgeIndex nextIndex() const
    {
      return m_index < lastIndex() ? m_index + 1 : 0;
    }
    EdgeIndex prevIndex() const
    {
      return m_index > 0 ? m_index - 1 : lastIndex();
    }

  private:
    EdgeIndex m_index;

    const int m_polygonID;
    const EdgeLabels& m_labels;
    const int m_numEdges;
  };

public:
  explicit DirectionalWalk(bool verbose = false) : m_verbose(verbose) { }

  /*!
  * Splits edges of each polygon based on the intersections with the other polygon
  * When the two polygons intersect, the split polygons are returned in \a psplit
  * The edges are labeled by the types of intersections in \a edgelabels
  * and \a orientation contains the relative orientation of the first intersection
  */
  int splitPolygonsAlongIntersections(const CurvedPolygonType& p1,
                                      const CurvedPolygonType& p2,
                                      double sq_tol)
  {
    // We store intersection information for each edge in EdgeIntersectionInfo structures
    std::vector<std::vector<EdgeIntersectionInfo>> E1IntData(p1.numEdges());
    std::vector<std::vector<EdgeIntersectionInfo>> E2IntData(p2.numEdges());
    EdgeIntersectionInfo firstinter;  // Need to do orientation test on first intersection

    // Find all intersections and store
    numinters = 0;
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

            if(m_verbose)
            {
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

      junctions = std::vector<Junction>(numinters + 1);

      if(m_verbose)
      {
        SLIC_INFO("Poly1 before split: " << psplit[0]);
        SLIC_INFO("Poly2 before split: " << psplit[1]);
      }

      splitPolygon(psplit[0], E1IntData, edgelabels[0], 0);
      splitPolygon(psplit[1], E2IntData, edgelabels[1], 1);

      if(m_verbose)
      {
        SLIC_INFO("Poly1 after split: " << psplit[0]);
        SLIC_INFO("Poly2 after split: " << psplit[1]);
      }
    }

    return numinters;
  }

public:
  /*!
   * Finds all intersection regions of the two input polygons
   *  
   * This function must be called after splitPolygonsAlongIntersections
   * which inserts the intersection points into each polygon and generates 
   * the \a edgelabels structure.
   * 
   * The results will be a vector of polygons inserted into OUT parameter \a pnew
   */
  void findIntersectionRegions(std::vector<CurvedPolygonType>& pnew)
  {
    PolygonEdge currentEdge[2] = {PolygonEdge(0, edgelabels[0]),
                                  PolygonEdge(1, edgelabels[1])};

    // We use junctionIndex to loop through the active junctions starting with index 1
    int junctionIndex = 1;

    // Find all connected regions of the intersection (we're done when we've found all junctions)
    const int numJunctions = junctions.size();
    while(junctionIndex < numJunctions)
    {
      // Attempt to get an initial "active" junction for this polygon
      Junction* startJunction = nullptr;
      do
      {
        startJunction = (junctionIndex < numJunctions)
          ? &(junctions[junctionIndex++])
          : nullptr;
      } while(startJunction != nullptr && !startJunction->isActive());

      // If we've found all the active junctions, we're done
      if(startJunction == nullptr)
      {
        break;
      }
      else
      {
        startJunction->junctionState = JunctionState::APPLIED;
      }

      // This variable allows us to switch between the two polygons
      bool active = orientation;

      // Set the index of the active edge to that of the start index
      PolygonEdge* activeEdge = &(currentEdge[active]);
      const auto startEdgeIndex = startJunction->nextEdgeIndex(active);
      activeEdge->setIndex(startEdgeIndex);

      if(m_verbose)
      {
        SLIC_INFO(""
                  << "Starting with junction " << startJunction->index
                  << "\n\t edges[0] {in: " << startJunction->currentEdgeIndex(0)
                  << " -- " << psplit[0][startJunction->currentEdgeIndex(0)]
                  << "; out: " << startJunction->nextEdgeIndex(0) << " -- "
                  << psplit[0][startJunction->nextEdgeIndex(0)] << "}"
                  << "\n\t edges[1] {in: " << startJunction->currentEdgeIndex(1)
                  << " -- " << psplit[1][startJunction->currentEdgeIndex(1)]
                  << "; out: " << startJunction->nextEdgeIndex(1) << " -- "
                  << psplit[1][startJunction->nextEdgeIndex(1)] << "}"
                  << "\n\t active: " << (active ? 1 : 0));
      }

      CurvedPolygonType aPart;  // Tracks the edges of the current intersection polygon

      // Each polygon iterates until it returns to the startJunctions
      Junction* junction = nullptr;
      while(junction == nullptr || *junction != *startJunction)
      {
        // Add all edges until we find a junction edge
        while(!activeEdge->isJunction())
        {
          const auto edgeIndex = activeEdge->getIndex();
          if(m_verbose)
          {
            SLIC_INFO("Adding edge (non-junction): " << psplit[active][edgeIndex]);
          }
          aPart.addEdge(psplit[active][edgeIndex]);
          activeEdge->advance();
        }

        // Add last leg of previous segment
        {
          const auto edgeIndex = activeEdge->getIndex();
          if(m_verbose)
          {
            SLIC_INFO("Adding edge (end of last): " << psplit[active][edgeIndex]);
          }
          aPart.addEdge(psplit[active][edgeIndex]);
        }

        // Handle junction
        junction = &(junctions[activeEdge->getEndLabel()]);
        SLIC_ASSERT(junction != nullptr);

        if(m_verbose)
        {
          SLIC_INFO(""
                    << "Swapped to junction " << junction->index
                    << "\n\t edges[0] {in: " << junction->currentEdgeIndex(0)
                    << " -- " << psplit[0][junction->currentEdgeIndex(0)]
                    << "; out: " << junction->nextEdgeIndex(0) << " -- "
                    << psplit[0][junction->nextEdgeIndex(0)] << "}"
                    << "\n\t edges[1] {in: " << junction->currentEdgeIndex(1)
                    << " -- " << psplit[1][junction->currentEdgeIndex(1)]
                    << "; out: " << junction->nextEdgeIndex(1) << " -- "
                    << psplit[1][junction->nextEdgeIndex(1)] << "}"
                    << "\n\t active: " << (active ? 1 : 0));
        }

        if(junction->isActive())
        {
          // swap active edge using junction data
          active = !active;
          const auto nextEdgeIndex = junction->nextEdgeIndex(active);
          activeEdge = &(currentEdge[active]);
          activeEdge->setIndex(nextEdgeIndex);

          junction->junctionState = JunctionState::APPLIED;

          if(m_verbose)
          {
            SLIC_INFO("Swapped to other polygon, edge index: "
                      << nextEdgeIndex
                      << "; edge: " << psplit[active][activeEdge->getIndex()]);
          }
        }
      }
      // Finalize polygon
      pnew.push_back(aPart);
    }
  }

private:
  /// Splits a CurvedPolygon \a p1 at every intersection point stored in \a IntersectionData
  /// \a edgeLabels stores intersection id for new vertices and 0 for original vertices
  void splitPolygon(CurvedPolygon<T, NDIMS>& polygon,
                    std::vector<std::vector<EdgeIntersectionInfo>>& IntersectionData,
                    std::vector<int>& edgelabels,
                    int polygonID)
  {
    using axom::utilities::isNearlyEqual;

    int addedIntersections = 0;
    const int numEdges = polygon.numEdges();
    for(int i = 0; i < numEdges; ++i)  // foreach edge
    {
      edgelabels.push_back(0);  //    mark start current vertex as 'original'
      const int nIntersect = IntersectionData[i].size();
      for(int j = 0; j < nIntersect; ++j)  //    foreach intersection on this edge
      {
        // split edge at parameter t_j
        const double t_j = IntersectionData[i][j].myTime;
        const int edgeIndex = i + addedIntersections;

        // TODO: Handle case where t_j is nearly 0.

        polygon.splitEdge(edgeIndex, t_j);

        // update edge label
        const int label = IntersectionData[i][j].numinter;
        edgelabels.insert(edgelabels.begin() + edgeIndex, label);
        junctions[label].edgeIndex[polygonID] = edgeIndex;
        junctions[label].junctionState =
          JunctionState::CROSS;  // TODO: Figure out what junction type needs to be
        junctions[label].index = label;

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

public:
  // Objects to store completely split polygons (split at every intersection point) and vector with unique id for each intersection and zeros for corners of original polygons.
  CurvedPolygonType psplit[2];  // The two completely split polygons will be stored in this array
  EdgeLabels edgelabels[2];  // 0 for curves that end in original vertices, unique id for curves that end in intersection points
  std::vector<Junction> junctions;
  int numinters {0};
  bool orientation {false};

  bool m_verbose {false};
};

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
  // Intersection is empty if either of the two polygons are empty
  if(p1.empty() || p2.empty())
  {
    return false;
  }

  DirectionalWalk<T, NDIMS> walk;

  int numinters = walk.splitPolygonsAlongIntersections(p1, p2, sq_tol);

  if(numinters > 0)
  {
    // This performs the directional walking method using the completely split polygons
    walk.findIntersectionRegions(pnew);
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

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECTION_CURVED_POLYGON_IMPL_HPP_
