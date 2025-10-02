// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_DISCRETIZE_DETAIL_
#define AXOM_QUEST_DISCRETIZE_DETAIL_

#include "axom/primal/constants.hpp"

namespace
{
enum
{
  P = 0,
  Q,
  R,
  S,
  T,
  U
};

using SphereType = axom::quest::SphereType;
using OctType = axom::quest::OctType;
using Point2D = axom::quest::Point2D;
using Point3D = axom::primal::Point<double, 3>;
using NAType = axom::NumericArray<double, 3>;

/* Return an octahedron whose six points lie on the truncated cone
 * described by rotating the line segment ab around the positive X-axis
 * in a right-handed way (thumb points out the axis, fingers spin segment
 * toward wrist).
 */
inline OctType from_segment(const Point2D &a, const Point2D &b)
{
  const double SQ_3_2 = sqrt(3.) / 2.;

  Point3D p {a[0], 0., a[1]};
  Point3D q {b[0], -SQ_3_2 * b[1], -0.5 * b[1]};
  Point3D r {a[0], -SQ_3_2 * a[1], -0.5 * a[1]};
  Point3D s {b[0], SQ_3_2 * b[1], -0.5 * b[1]};
  Point3D t {a[0], SQ_3_2 * a[1], -0.5 * a[1]};
  Point3D u {b[0], 0., b[1]};

  return OctType(p, q, r, s, t, u);
}

/* How many octahedra will we generate in each segment of the polyline,
 * summed over all levels of refinement?  We generate the octahedra in
 * discretizeSegment() (which see for further details).
 *  - One oct for the central octahedron (level 0),
 *  - three more for its side faces (level 1),
 *  - and for level i > 1, two times the octahedron count in level i-1,
 *    to refine each exposed face.
 *
 * Total = 1 + 3*sum[i=0 to levels-1](2^i)
 *
 * For levels < 0, we return 0.  This lets us call the routine to find
 * the offset in an array to store prisms for level (levels + 1).
 */
inline int count_segment_prisms(int levels)
{
  int octcount = 1;
  for(int level = levels; level > 0; --level)
  {
    if(level == 1)
    {
      octcount *= 3;
    }
    else
    {
      octcount *= 2;
    }
    octcount += 1;
  }
  if(levels < 0)
  {
    octcount = 0;
  }

  return octcount;
}

AXOM_HOST_DEVICE
Point3D rescale_YZ(const Point3D &p, double new_dst)
{
  const double cur_dst =
    axom::utilities::clampLower(sqrt(p[1] * p[1] + p[2] * p[2]), axom::primal::PRIMAL_TINY);

  Point3D retval;
  retval[0] = p[0];
  retval[1] = p[1] * new_dst / cur_dst;
  retval[2] = p[2] * new_dst / cur_dst;
  return retval;
}

AXOM_HOST_DEVICE
inline OctType new_inscribed_prism(OctType &old_oct,
                                   int p_off,
                                   int s_off,
                                   int t_off,
                                   int u_off,
                                   Point2D pa,
                                   Point2D pb)
{
  OctType retval;
  retval[P] = old_oct[p_off];
  Point3D new_q = Point3D::lerp(old_oct[u_off], old_oct[s_off], 0.5);
  retval[Q] = rescale_YZ(new_q, pb[1]);
  Point3D new_r = Point3D::lerp(old_oct[p_off], old_oct[t_off], 0.5);
  retval[R] = rescale_YZ(new_r, pa[1]);
  retval[S] = old_oct[s_off];
  retval[T] = old_oct[t_off];
  retval[U] = old_oct[u_off];
  return retval;
}

/* ------------------------------------------------------------ */
/* Discretize a cylinder (or truncated cone) into a hierarchy of
 * triangular prisms (or truncated tetrahedra) stored as octahedra.
 * Each level of refinement places a new prism/tet on all exposed faces.
 *
 * Input:  2D points a and b, number of levels of refinement, index into
 * the output array.  The routine assumes a.x <= b.x, a.y >= 0, b.y >= 0.
 * The routine also assumes that the output array is resized to hold
 * the generated prisms at the specified index.
 *
 * Output: the output array of prisms (placed in the output array
 * starting at the index).
 *
 * Return value: the number of prisms generated: zero for degenerate segments,
 * or a value of order 2^(levels) calculated within the routine.
 *
 * Conceptually, end points a and b are revolved around the X-axis, describing
 * circles that are the truncated cone's end-caps.  The segment ab revolved
 * about the X-axis becomes the side-wall.
 *
 * The level-zero prism is constructed by inscribing a triangle in each of the
 * end-cap circles.  In one circle, put the points PRT; in the other circle, put
 * the points UQS.  The edges UP, QR, and ST lie in the side-wall described by
 * rotating the segment ab about the x-axis and are each co-planar with the
 * X-axis.  (This five-sided prism, with two end-caps, three side-walls, and six
 * vertices, is stored in an Octahedron data record.  The edges PQ, RS, and TU
 * split the quadrilateral prism side-walls into pairs of coplanar triangles.)
 *
 * Each subsequent level of refinement adds a prism to each exposed
 * quadrilateral side-wall.
 */
template <typename ExecSpace>
int discrSeg(const Point2D &a, const Point2D &b, int levels, axom::ArrayView<OctType> &out, int idx)
{
  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Assert input assumptions
  SLIC_ASSERT(b[0] - a[0] >= 0);
  SLIC_ASSERT(a[1] >= 0);
  SLIC_ASSERT(b[1] >= 0);

  // Deal with degenerate segments
  if(b[0] - a[0] < axom::primal::PRIMAL_TINY)
  {
    return 0;
  }
  if(a[1] < axom::primal::PRIMAL_TINY && b[1] < axom::primal::PRIMAL_TINY)
  {
    return 0;
  }

  int total_count = count_segment_prisms(levels);

  // Establish a prism (in an octahedron record) with one triangular
  // end lying on the circle described by rotating point a around the
  // x-axis and the other lying on circle from rotating b.
  OctType *oct_from_seg = axom::allocate<OctType>(1, hostAllocID);
  oct_from_seg[0] = from_segment(a, b);

  axom::copy(out.data() + idx + 0, oct_from_seg, sizeof(OctType));
  axom::deallocate(oct_from_seg);

  // curr_lvl indexes to the first prism in the level we're currently refining
  int curr_lvl = idx;
  int curr_lvl_count = 1;
  int next_lvl = curr_lvl + curr_lvl_count;

  Point2D pa, pb;

  // Refine: add an octahedron to each exposed face.  Perform "levels"
  // refinements beyond the level-0 octahedron.
  for(int level = 0; level < levels; ++level)
  {
    // Each level of refinement generates a prism for each exposed
    // face, so twice the prisms in the preceding level.  Refining the
    // initial prism is the only different step, since all three of its
    // side-faces are exposed.
    int lvl_factor = 2;
    if(level == 0)
    {
      lvl_factor = 3;
    }
    //int index = count_segment_prisms(level - 1);

    // The ends of the prisms switch each level.
    if(level & 1)
    {
      pa = a;
      pb = b;
    }
    else
    {
      pa = b;
      pb = a;
    }

    // This loop generates the prisms of the next level of refinement.
    // The specified vertices ensure that the new prisms always have
    // triangular end-caps QSU and RTP, and that side-face PTSU faces
    // the parent-level prism.
    // Of note, the child-level end-cap QSU is coplanar with parent-
    // level cap RTP, and vice versa.  Hence the preceding if-statement
    // with comment "the ends switch each level."
    axom::for_all<ExecSpace>(
      curr_lvl_count,
      AXOM_LAMBDA(axom::IndexType i) mutable {
        out[next_lvl + i * lvl_factor + 0] =
          new_inscribed_prism(out[curr_lvl + i], Q, T, S, R, pa, pb);
        out[next_lvl + i * lvl_factor + 1] =
          new_inscribed_prism(out[curr_lvl + i], U, R, Q, P, pa, pb);
        if(level == 0)
        {
          out[next_lvl + i * lvl_factor + 2] =
            new_inscribed_prism(out[curr_lvl + i], S, P, U, T, pa, pb);
        }
      });

    curr_lvl = next_lvl;
    curr_lvl_count *= lvl_factor;
    next_lvl = curr_lvl + curr_lvl_count;
  }

  return total_count;
}

}  // end anonymous namespace

namespace axom
{
namespace quest
{
/* Given a surface of revolution and level of refinement, place the list of
 * octahedra that approximate that shape at that LOR in an output argument.
 * Return true for valid input and lack of errors, false otherwise.  If we
 * return false, put nothing in the output argument.
 *
 * As noted in the doxygen, this function chops a surface of revolution into
 * n*O(2^levels) octahedra, where n is the number of segments in the SoR (one
 * less than the polyline's length).  That is exponential growth.  Use
 * appropriate caution.
 *
 * This routine initializes an Array pointed to by \a out.
 */
template <typename ExecSpace>
bool discretize(const axom::ArrayView<Point2D> &polyline,
                int pointcount,
                int levels,
                axom::Array<OctType> &out,
                int &octcount)
{
  int allocId = axom::execution_space<ExecSpace>::allocatorID();
  // Check for invalid input.  If any segment is invalid, exit returning false.
  bool stillValid = true;
  int segmentcount = pointcount - 1;
  for(int seg = 0; seg < segmentcount && stillValid; ++seg)
  {
    const Point2D &a = polyline[seg];
    const Point2D &b = polyline[seg + 1];
    // invalid if a.x > b.x
    if(a[0] > b[0])
    {
      stillValid = false;
    }
    if(a[1] < 0 || b[1] < 0)
    {
      stillValid = false;
    }
  }
  if(!stillValid)
  {
    return false;
  }

  int segoctcount = count_segment_prisms(levels);

  // That was the octahedron count for one segment.  Multiply by the number
  // of segments we will compute.
  int totaloctcount = segoctcount * segmentcount;
  out = axom::Array<OctType>(totaloctcount, totaloctcount, allocId);
  axom::ArrayView<OctType> out_view = out.view();
  octcount = 0;

  for(int seg = 0; seg < segmentcount; ++seg)
  {
    int segment_prism_count =
      discrSeg<ExecSpace>(polyline[seg], polyline[seg + 1], levels, out_view, octcount);
    octcount += segment_prism_count;
  }

  // TODO check for errors in each segment's computation
  return true;
}

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_DISCRETIZE_DETAIL_
