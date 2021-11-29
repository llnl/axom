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

#include "axom/core.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"

#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/operators/orientation.hpp"
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

template <typename T>
class DirectionalWalk
{
public:
  static constexpr int NDIMS = 2;
  using CurvedPolygonType = CurvedPolygon<T, NDIMS>;
  using IndexArray = std::vector<int>;

public:
  /*! \class EdgeIntersectionInfo
  *
  * \brief For storing intersection points between edges of \a CurvedPolygon instances so they can be easily sorted by parameter value using std::sort
  */
  struct EdgeIntersectionInfo
  {
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
    GRAZE          /// Atypical case: polygons intersect at a terminal vertex
  };

  /*!
   * Helper class for encoding information about the incident edges at an interection
   */
  struct Junction
  {
    using JunctionIndex = int;
    using EdgeIndex = int;
    static constexpr EdgeIndex INVALID_EDGE_INDEX = -1;
    static constexpr int INVALID_JUNCTION_INDEX = -1;
    static constexpr int NON_JUNCTION_INDEX = 0;

    bool isActive() const { return junctionState > JunctionState::APPLIED; }

    EdgeIndex currentEdgeIndex(bool active) const
    {
      return incomingEdgeIndex[active];
    }
    EdgeIndex nextEdgeIndex(bool active) const
    {
      return incomingEdgeIndex[active] + 1;
    }

    bool operator==(const Junction& other) const
    {
      return index == other.index;
    }
    bool operator!=(const Junction& other) const { return !(*this == other); }

  public:
    // indices of polygon edges leading into junction
    EdgeIndex incomingEdgeIndex[2] {INVALID_EDGE_INDEX, INVALID_EDGE_INDEX};
    // describes the junction type and status
    JunctionState junctionState {JunctionState::UNINITIALIZED};
    JunctionIndex index {INVALID_JUNCTION_INDEX};
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

    PolygonEdge(int polygonID, const IndexArray& endpointJunctionIndices)
      : m_polygonID(polygonID)
      , m_endpointJunctionIndices(endpointJunctionIndices)
      , m_numEdges(endpointJunctionIndices.size())
    {
      SLIC_ASSERT_MSG(m_numEdges > 0, "Polygon " << polygonID << " was empty");
    }

    EdgeIndex getIndex() const { return m_index; }
    void setIndex(EdgeIndex idx) { m_index = (idx % m_numEdges); }

    int getStartLabel() const { return m_endpointJunctionIndices[prevIndex()]; }
    int getEndLabel() const { return m_endpointJunctionIndices[m_index]; }
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
    const IndexArray& m_endpointJunctionIndices;
    const int m_numEdges;
  };

public:
  explicit DirectionalWalk(bool verbose = true) : m_verbose(verbose) { }

  /*!
  * Splits edges of each polygon based on the intersections with the other polygon
  * When the two polygons intersect, the split polygons are returned in \a psplit
  * The edges are labeled by the types of intersections in \a endpointJunctionIndices
  * and \a orientation contains the relative orientation of the first intersection
  */
  int splitPolygonsAlongIntersections(const CurvedPolygonType& p1,
                                      const CurvedPolygonType& p2,
                                      double sq_tol)
  {
    // We store intersection information for each edge in EdgeIntersectionInfo structures
    using EdgeIntersectionInfoArray =
      std::vector<std::vector<EdgeIntersectionInfo>>;
    EdgeIntersectionInfoArray p1IntersectionData(p1.numEdges());
    EdgeIntersectionInfoArray p2IntersectionData(p2.numEdges());
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
            p1IntersectionData[i].push_back(
              {p1times[k], i, p2times[k], j, numinters + 1});
            p2IntersectionData[j].push_back(
              {p2times[k], j, p1times[k], i, numinters + 1});

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
        std::sort(p1IntersectionData[i].begin(), p1IntersectionData[i].end());
      }
      for(int i = 0; i < p2.numEdges(); ++i)
      {
        std::sort(p2IntersectionData[i].begin(), p2IntersectionData[i].end());
      }

      psplit[0] = p1;
      psplit[1] = p2;

      endpointJunctionIndices[0].reserve(p1.numEdges() + numinters);
      endpointJunctionIndices[1].reserve(p2.numEdges() + numinters);

      junctions = std::vector<Junction>(numinters + 1);

      if(m_verbose)
      {
        SLIC_INFO("Poly1 before split: " << psplit[0]);
        SLIC_INFO("Poly2 before split: " << psplit[1]);
      }

      this->splitPolygon(0, p1IntersectionData);
      this->splitPolygon(1, p2IntersectionData);

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
   * the \a endpointJunctionIndices structure.
   *
   * The results will be a vector of polygons inserted into OUT parameter \a pnew
   */
  void findIntersectionRegions(std::vector<CurvedPolygonType>& pnew)
  {
    PolygonEdge currentEdge[2] = {PolygonEdge(0, endpointJunctionIndices[0]),
                                  PolygonEdge(1, endpointJunctionIndices[1])};

    if(m_verbose)
    {
      const int nJunctions = junctions.size();
      SLIC_INFO("At start of 'findIntersectionRegions', there are "
                << nJunctions << " junctions:");
      for(int i = 0; i < nJunctions; ++i)
      {
        SLIC_INFO(
          axom::fmt::format("\tJunction {}: incomingEdge[0]: {}, "
                            "incomingEdge[1]: {}, state: {}, index: {} ",
                            i,
                            junctions[i].incomingEdgeIndex[0],
                            junctions[i].incomingEdgeIndex[1],
                            junctions[i].junctionState,
                            junctions[i].index));
      }
    }

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
      // switch(startJunction->junctionState)
      // {
      // case JunctionState::CROSS:
      //   break;
      // case JunctionState::GRAZE:
      //   break;
      // default:
      //   break;
      // }
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
        const int p0Edges = psplit[0].numEdges();
        const int p1Edges = psplit[1].numEdges();
        const int p0CurrIdx = startJunction->currentEdgeIndex(0);
        const int p0NextIdx = startJunction->nextEdgeIndex(0) % p0Edges;
        const int p1CurrIdx = startJunction->currentEdgeIndex(1);
        const int p1NextIdx = startJunction->nextEdgeIndex(1) % p1Edges;

        SLIC_INFO(
          axom::fmt::format("Starting with junction {}"
                            "\n\t edges[0] {} -> {} {{in: {}; out: {} }}"
                            "\n\t edges[1] {} -> {} {{in: {}; out: {} }}"
                            "\n\t active: {}",
                            startJunction->index,
                            p0CurrIdx,
                            p0NextIdx,
                            psplit[0][p0CurrIdx],
                            psplit[0][p0NextIdx],
                            p1CurrIdx,
                            p1NextIdx,
                            psplit[1][p1CurrIdx],
                            psplit[1][p1NextIdx],
                            (active ? 1 : 0)));
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
  /// Splits the polygon with id \a polygonID at every intersection point in \a edgeIntersections.
  /// Uses internal array \a edgeLabels to store ids for vertices (junction indices or 0 for original vertices)
  void splitPolygon(
    int polygonID,
    std::vector<std::vector<EdgeIntersectionInfo>>& allEdgeIntersections)
  {
    using axom::utilities::isNearlyEqual;

    CurvedPolygonType& polygon = psplit[polygonID];
    IndexArray& junctionIndices = endpointJunctionIndices[polygonID];

    bool fixupBeginning = false;
    int fixupBeginningJunctionIdx = Junction::INVALID_JUNCTION_INDEX;

    int addedIntersections = 0;
    const int numOrigEdges = polygon.numEdges();
    for(int i = 0; i < numOrigEdges; ++i)  // foreach edge
    {
      // mark this edge's endpoint as 'original'
      junctionIndices.push_back(Junction::NON_JUNCTION_INDEX);
      double previous_tj = 0.;
      auto& curEdgeIntersections = allEdgeIntersections[i];
      const int nIntersect = curEdgeIntersections.size();
      for(int j = 0; j < nIntersect; ++j)  //    foreach intersection on this edge
      {
        // split edge at parameter t_j
        const double t_j = curEdgeIntersections[j].myTime;
        const int edgeIndex = i + addedIntersections;

        if(m_verbose)
        {
          SLIC_INFO(
            fmt::format("i {}, j {}, added {}, t_j {}, previous t_j {}, "
                        "edgeIndex {} (sz: {}), polygon size {}",
                        i,
                        j,
                        addedIntersections,
                        t_j,
                        previous_tj,
                        edgeIndex,
                        junctionIndices.size(),
                        polygon.numEdges()));
        }

        const bool nearlyZero = isNearlyEqual(t_j, 0.);

        // TODO: Handle case where t_j is nearly 0.
        if(!nearlyZero)
        {
          polygon.splitEdge(edgeIndex, t_j);

          // update edge label
          const int idx = curEdgeIntersections[j].numinter;
          junctionIndices.insert(junctionIndices.begin() + edgeIndex, idx);
          junctions[idx].incomingEdgeIndex[polygonID] = edgeIndex;
          //if(junctions[idx].junctionState == JunctionState::UNINITIALIZED)
          //{
          junctions[idx].junctionState = JunctionState::CROSS;
          //}
          junctions[idx].index = idx;

          ++addedIntersections;
          previous_tj = t_j;
        }
        else
        {
          const int idx = curEdgeIntersections[j].numinter;

          if(edgeIndex > 0)
          {
            junctionIndices[edgeIndex] = idx;
            junctions[idx].incomingEdgeIndex[polygonID] = edgeIndex;
          }
          else
          {
            fixupBeginning = true;
            fixupBeginningJunctionIdx = idx;
          }

          junctions[idx].junctionState =
            JunctionState::CROSS;  //JunctionState::GRAZE;
          junctions[idx].index = idx;
        }

        if(m_verbose)
        {
          SLIC_INFO(axom::fmt::format("junctionIndices: {}\n junctions:",
                                      axom::fmt::join(junctionIndices, " ")));
          for(const auto& jj : junctions)
          {
            SLIC_INFO(fmt::format(
              "\t{{index: {}, state: {}, edgeIndex[0]: {}, edgeIndex[1]: {} }}",
              jj.index,
              jj.junctionState,
              jj.incomingEdgeIndex[0],
              jj.incomingEdgeIndex[1]));
          }
        }

        // update remaining intersections on this edge; special case if already at end of curve
        if(!isNearlyEqual(1., t_j))
        {
          for(int k = j + 1; k < nIntersect; ++k)
          {
            const double t_k = curEdgeIntersections[k].myTime;
            curEdgeIntersections[k].myTime = (t_k - t_j) / (1.0 - t_j);
          }
        }
        else
        {
          for(int k = j + 1; k < nIntersect; ++k)
          {
            curEdgeIntersections[k].myTime = 1.;
          }
        }
      }
    }

    // If the first vertex of the first edge was a junction, we need to fix it up using the current last edge
    if(fixupBeginning)
    {
      const int edgeIndex = junctionIndices.size() - 1;
      junctionIndices[edgeIndex] = fixupBeginningJunctionIdx;
      junctions[fixupBeginningJunctionIdx].incomingEdgeIndex[polygonID] =
        edgeIndex;
    }
  }

public:
  // Objects to store completely split polygons (split at every intersection point) and vector with unique id for each intersection and zeros for corners of original polygons.
  CurvedPolygonType psplit[2];  // The two completely split polygons will be stored in this array
  IndexArray endpointJunctionIndices[2];  // 0 for curves that end in original vertices, unique id for curves that end in intersection points
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
template <typename T>
bool intersect_polygon(const CurvedPolygon<T, 2>& p1,
                       const CurvedPolygon<T, 2>& p2,
                       std::vector<CurvedPolygon<T, 2>>& pnew,
                       double sq_tol)
{
  // Intersection is empty if either of the two polygons are empty
  if(p1.empty() || p2.empty())
  {
    return false;
  }

  DirectionalWalk<T> walk;

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

enum class JunctionIntersectionType : int
{
  Non_Intersection = -1,  //
  Cross,
  Type1,
  Type2,
  Type3,
  Type4
};

// when interesection is at a vertex of the original polygon(s),
// determine if the in/out edges of second polygon is on the same
// side of the first polygon
template <typename T>
JunctionIntersectionType getJunctionIntersectionType(
  const BezierCurve<T, 2>& p0In,
  const BezierCurve<T, 2>& p0Out,
  const BezierCurve<T, 2>& p1In,
  const BezierCurve<T, 2>& p1Out)
{
  using SegmentType = primal::Segment<T, 2>;
  using VectorType = primal::Vector<T, 2>;
  using PointType = primal::Point<T, 2>;

  // TODO: Error checking to ensure the four curves have a common intersection point
  // P0In(1) == P0Out(0) == P1In(1) == P1Out(0)
  // If not, return JunctionIntersectionType::Non_Intersection

  const VectorType p0Tangents[2] = {-p0In.dt(1.).unitVector(),
                                    p0Out.dt(0.).unitVector()};
  const SegmentType p0Seg(PointType {p0Tangents[0][0], p0Tangents[0][1]},
                          PointType {p0Tangents[1][0], p0Tangents[1][1]});

  const VectorType p1Tangents[2] = {-p1In.dt(1.).unitVector(),
                                    p1Out.dt(0.).unitVector()};
  const SegmentType p1Seg(PointType {p1Tangents[0][0], p1Tangents[0][1]},
                          PointType {p1Tangents[1][0], p1Tangents[1][1]});

  const int orientIn = primal::orientation(p1Seg[0], p0Seg);
  const int orientOut = primal::orientation(p1Seg[1], p0Seg);

  if(orientIn == orientOut)
  {
    const int orient_p1_p0In = primal::orientation(p0Seg[0], p1Seg);
    if(orientIn == ON_POSITIVE_SIDE)
    {
      return orient_p1_p0In == ON_POSITIVE_SIDE
        ? JunctionIntersectionType::Type3
        : JunctionIntersectionType::Type1;
    }
    else
    {
      {
        return orient_p1_p0In == ON_POSITIVE_SIDE
          ? JunctionIntersectionType::Type4
          : JunctionIntersectionType::Type2;
      }
    }
  }
  else
  {
    return JunctionIntersectionType::Cross;
  }
}

template <typename T>
constexpr int DirectionalWalk<T>::Junction::NON_JUNCTION_INDEX;

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECTION_CURVED_POLYGON_IMPL_HPP_
