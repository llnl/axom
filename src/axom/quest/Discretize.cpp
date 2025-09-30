// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/Discretize.hpp"

#include "axom/core.hpp"

#include "axom/primal/constants.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/operators/split.hpp"

#include "axom/mint/config.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/mint/mesh/CellTypes.hpp"

#include <cmath>

namespace axom
{
namespace quest
{
/* ------------------------------------------------------------ */
/* Project a Point onto a sphere.
 */
Point3D project_to_shape(const Point3D& p, const SphereType& sphere)
{
  const auto& ctr = sphere.getCenter();
  const double dist2 = primal::squared_distance(ctr, p);
  const double dist = sqrt(dist2);
  const double drat = sphere.getRadius() * dist / (dist2 + primal::PRIMAL_TINY);

  return Point3D::lerp(ctr, p, drat);
}

/* Return an octahedron whose six points lie on the given sphere.
 */
OctType from_sphere(const SphereType& sphere)
{
  const NAType& center = sphere.getCenter().array();
  const NAType ihat({1., 0., 0.});
  const NAType jhat({0., 1., 0.});
  const NAType khat({0., 0., 1.});

  const Point3D P = project_to_shape(Point3D(center + ihat), sphere);
  const Point3D Q = project_to_shape(Point3D(center + jhat), sphere);
  const Point3D R = project_to_shape(Point3D(center + khat), sphere);
  const Point3D S = project_to_shape(Point3D(center - ihat), sphere);
  const Point3D T = project_to_shape(Point3D(center - jhat), sphere);
  const Point3D U = project_to_shape(Point3D(center - khat), sphere);

  return OctType(P, Q, R, S, T, U);
}

/* How many octahedra will we generate?  One oct for the central octahedron
 * (level 0), eight more for its faces (level 1), and four more to refine
 * each exposed face of the last generation.
 *
 * Total = 1 + 8*sum[i=0 to levels-1](4^i)
 */
int count_sphere_octahedra(int levels)
{
  int octcount = 1;
  for(int level = levels; level > 0; --level)
  {
    octcount *= 4;
    if(level == 1)
    {
      octcount *= 2;
    }
    octcount += 1;
  }

  return octcount;
}

/* Given a sphere, a parent octahedron with vertices lying on the
 * sphere, and vertex indices s, t, u defining a face on that
 * octahedron, return a new oct sharing the face (s,t,u) and all other
 * faces looking "outward" toward the sphere.
 */
OctType new_inscribed_oct(const SphereType& sphere, OctType& o, int s, int t, int u)
{
  Point3D P = Point3D::midpoint(o[t], o[u]);
  Point3D Q = Point3D::midpoint(o[s], o[u]);
  Point3D R = Point3D::midpoint(o[s], o[t]);

  P = project_to_shape(P, sphere);
  Q = project_to_shape(Q, sphere);
  R = project_to_shape(R, sphere);

  return OctType(P, Q, R, o[s], o[t], o[u]);
}

/* Given a sphere and level of refinement, place the list of octahedra that
 * approximate that sphere at that LOR into the output argument.  Return true for
 * valid input and no errors, false for bad input or an error.  If we return
 * false, no octahedra are placed in the output argument.
 *
 * As noted in the doxygen, this function chops a sphere into O(4^levels)
 * octahedra.  That is exponential growth.  Use appropriate caution.
 *
 * This routine allocates an array pointed to by \a out.  The caller is responsible
 * to free the array.
 */
bool discretize(const SphereType& sphere, int levels, axom::Array<OctType>& out, int& octcount)
{
  // Check input.  Negative radius: return false.
  if(sphere.getRadius() < 0)
  {
    return false;
  }
  // Zero radius: return true without generating octahedra.
  if(sphere.getRadius() < primal::PRIMAL_TINY)
  {
    octcount = 0;
    return true;
  }

  octcount = count_sphere_octahedra(levels);

  out = axom::Array<OctType>(octcount, octcount);

  // index points to an octahedron of the last generation.  We'll generate
  // new octahedra based on out[index].
  int index = 0;

  // outindex points to where the next new octahedron should be stored.
  int outindex = 0;

  // Establish an octahedron with all of its points lying on the sphere (level 0)
  out[outindex++] = from_sphere(sphere);

  // max_last_gen indexes to the last oct of the last generation.
  int max_last_gen = 0;

  // Refine: add an octahedron to each exposed face.  Perform "levels"
  // refinements beyond the level-0 octahedron.
  for(int level = 0; level < levels; ++level)
  {
    max_last_gen = outindex - 1;
    while(index <= max_last_gen)
    {
      // Octahedra are defined by points P, Q, R, S, T, U (indexes 0--5; see
      // enum at file beginning).  Point oct[i] is opposite point oct[(i+3)%6].
      // Convention for new octahedra: P, Q, R are new points from midpoints,
      // S, T, U are inherited from the last gen.
      /* newoct[0] uses (P,Q,R)-old. */
      /* newoct[1] uses (T,P,R)-old. */
      /* newoct[2] uses (P,U,Q)-old. */
      /* newoct[3] uses (R,Q,S)-old. */
      out[outindex++] = new_inscribed_oct(sphere, out[index], P, Q, R);
      out[outindex++] = new_inscribed_oct(sphere, out[index], T, P, R);
      out[outindex++] = new_inscribed_oct(sphere, out[index], P, U, Q);
      out[outindex++] = new_inscribed_oct(sphere, out[index], R, Q, S);
      if(index == 0)
      {
        out[outindex++] = new_inscribed_oct(sphere, out[index], P, T, U);
        out[outindex++] = new_inscribed_oct(sphere, out[index], Q, U, S);
        out[outindex++] = new_inscribed_oct(sphere, out[index], T, R, S);
        out[outindex++] = new_inscribed_oct(sphere, out[index], U, T, S);
      }

      index += 1;
    }
  }

  return true;
}

namespace
{
constexpr mint::CellType CELL_TYPE = mint::TET;
using TetType = primal::Tetrahedron<double, 3>;
using PolyhedronType = primal::Polyhedron<double, 3>;

double octPolyVolume(const OctType& o)
{
  // Convert Octahedron into Polyhedrom
  PolyhedronType octPoly;

  octPoly.addVertex(o[0]);
  octPoly.addVertex(o[1]);
  octPoly.addVertex(o[2]);
  octPoly.addVertex(o[3]);
  octPoly.addVertex(o[4]);
  octPoly.addVertex(o[5]);

  octPoly.addNeighbors(0, {1, 5, 4, 2});
  octPoly.addNeighbors(1, {0, 2, 3, 5});
  octPoly.addNeighbors(2, {0, 4, 3, 1});
  octPoly.addNeighbors(3, {1, 2, 4, 5});
  octPoly.addNeighbors(4, {0, 5, 3, 2});
  octPoly.addNeighbors(5, {0, 1, 3, 4});

  return octPoly.volume();
}
}  // namespace

//------------------------------------------------------------------------------
int mesh_from_discretized_polyline(const axom::ArrayView<OctType>& octs,
                                   int octcount,
                                   int segcount,
                                   mint::Mesh*& mesh)
{
  SLIC_ASSERT(octs.data() != nullptr);

  const int tetcount = 8 * octcount;
  const int vertcount = 4 * tetcount;
  int octPerSeg = octcount / segcount;
  SLIC_ASSERT_MSG((octcount % segcount) == 0,  // remainderOcts
                  "Total octahedron count is not evenly divisible by segment count");

  // Step 0: create the UnstructuredMesh
  mint::UnstructuredMesh<mint::SINGLE_SHAPE>* um =
    new mint::UnstructuredMesh<mint::SINGLE_SHAPE>(3, CELL_TYPE, vertcount, tetcount);

  // Step 1: Add fields
  int* octlevel = um->createField<int>("level_of_refinement", mint::CELL_CENTERED);
  int* octidx = um->createField<int>("octahedron_index", mint::CELL_CENTERED);
  int* segidx = um->createField<int>("segment_index", mint::CELL_CENTERED);
  double* vol = um->createField<double>("octahedron_volume", mint::CELL_CENTERED);
  double* pvol = um->createField<double>("oct_as_polyhedron_volume", mint::CELL_CENTERED);

  // Step 2: for each oct,
  //    - split it into tets
  //    - add the nodes
  //    - add the cells
  //    - add the fields
  constexpr int TETS_PER_OCT = 8;
  constexpr int NODES_PER_TET = 4;
  int level = 0, octInLevel = 1, maxOctIdxInLevel = 0;
  Array<TetType> tets;
  for(int o = 0; o < octcount; ++o)
  {
    tets.clear();
    int segmentIndex = o / octPerSeg;  // integer division
    int octInSegment = o % octPerSeg;
    if(octInSegment == 0)
    {
      level = 0;
      octInLevel = 1;
      maxOctIdxInLevel = 0;
    }
    if(octInSegment > maxOctIdxInLevel)
    {
      level += 1;
      if(level == 1)
      {
        octInLevel *= 3;
      }
      else
      {
        octInLevel *= 2;
      }
      maxOctIdxInLevel += octInLevel;
    }

    primal::split(octs[o], tets);

    double octvol = 0.;
    for(int t = 0; t < 8; ++t)
    {
      TetType& tet = tets[t];
      for(int n = 0; n < 4; ++n)
      {
        um->appendNode(tet[n][0], tet[n][1], tet[n][2]);
      }
      axom::IndexType nidx = o * (TETS_PER_OCT * NODES_PER_TET) + t * NODES_PER_TET;
      axom::IndexType cell[NODES_PER_TET] = {nidx + 0, nidx + 2, nidx + 1, nidx + 3};
      um->appendCell(cell);

      octvol += tets[t].signedVolume();
    }

    double pvolume = octPolyVolume(octs[o]);

    axom::IndexType idx = o * TETS_PER_OCT;
    for(int t = 0; t < 8; ++t)
    {
      vol[idx + t] = octvol;
      pvol[idx + t] = pvolume;
      octidx[idx + t] = octInSegment;
      segidx[idx + t] = segmentIndex;
      octlevel[idx + t] = level;
    }
  }

  mesh = um;
  return 0;
}

}  // end namespace quest
}  // end namespace axom
