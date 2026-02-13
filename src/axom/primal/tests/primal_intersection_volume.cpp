// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace Primal3D
{
using PointType = axom::primal::Point<double, 3>;
using VectorType = axom::primal::Vector<double, 3>;
using BoundingBoxType = axom::primal::BoundingBox<double, 3>;
using HexahedronType = axom::primal::Hexahedron<double, 3>;
using TriangleType = axom::primal::Triangle<double, 3>;
using TetrahedronType = axom::primal::Tetrahedron<double, 3>;
using OctahedronType = axom::primal::Octahedron<double, 3>;
using PolyhedronType = axom::primal::Polyhedron<double, 3>;
using PolygonType = axom::primal::Polygon<double, 3>;
using PlaneType = axom::primal::Plane<double, 3>;
using PolyhedronType = axom::primal::Polyhedron<double, 3>;
}  // namespace Primal3D

namespace
{
using Primal3D::HexahedronType;
using Primal3D::TetrahedronType;

double intersection_volume_via_hex_triangulation(const HexahedronType& hex,
                                                 const TetrahedronType& tet,
                                                 double eps)
{
  axom::StackArray<TetrahedronType, 24> tets;
  hex.triangulate(tets);

  double vol = 0.;
  constexpr bool TRY_FIX_ORIENTATION = true;
  for(int i = 0; i < 24; ++i)
  {
    vol += axom::primal::intersection_volume<double>(tets[i], tet, eps, TRY_FIX_ORIENTATION);
  }
  return vol;
}

/**
 * This utility class computes the intersection volume of a pair of tetrahedra, or a tetrahedron and a hexahedron
 * by finding the convex hull of their component planes. We're using it in this test as a comparison to the
 * primal::intersection_volume routines, which uses clipping planes and the primal::Polyhedron class
 */
class HalfspaceIntersectionVolumeLD
{
public:
  using Real = long double;
  using Point3 = axom::primal::Point<Real, 3>;
  using Vec3 = axom::primal::Vector<Real, 3>;
  using Vec2 = axom::primal::Vector<Real, 2>;
  using Plane = axom::primal::Plane<Real, 3>;

  /// Finds the intersection volume of a pair of tets by via the convex hull of their intersections
  Real intersection_volume(const TetrahedronType& a, const TetrahedronType& b) const
  {
    // add the planes
    std::vector<Plane> planes;
    for(const auto& p : tet_halfspaces(a))
    {
      planes.push_back(p);
    }
    for(const auto& p : tet_halfspaces(b))
    {
      planes.push_back(p);
    }

    // get scale for relative tolerance
    const Real scale = [&a, &b]() {
      Real max_abs = 0.;
      for(int i = 0; i < 4; ++i)
      {
        max_abs = std::max(max_abs, abs(to_point_ld(a[i]).array()).max());
        max_abs = std::max(max_abs, abs(to_point_ld(b[i]).array()).max());
      }
      return max_abs;
    }();

    constexpr Real abs_tol = 1e-14L;
    const Real tol = abs_tol * std::max(Real(1), scale);

    // remove duplicate planes (up to tolerance)
    planes = dedup_planes(planes, abs_tol);

    // find planes defining the convex hull
    std::vector<Point3> candidates;
    candidates.reserve(64);

    for(std::size_t i = 0; i < planes.size(); ++i)
    {
      for(std::size_t j = i + 1; j < planes.size(); ++j)
      {
        for(std::size_t k = j + 1; k < planes.size(); ++k)
        {
          Point3 x;
          if(!intersect_three_planes(planes[i], planes[j], planes[k], x, abs_tol))
          {
            continue;
          }

          bool inside = true;
          for(const auto& p : planes)
          {
            if(p.signedDistance(x) > tol)
            {
              inside = false;
              break;
            }
          }
          if(inside)
          {
            candidates.push_back(x);
          }
        }
      }
    }

    const std::vector<Point3> verts = unique_points(candidates, tol);
    return volume_from_halfspaces(planes, verts);
  }

  Real intersection_volume(const HexahedronType& hex, const TetrahedronType& tet) const
  {
    axom::StackArray<TetrahedronType, 24> tets;
    hex.triangulate(tets);

    Real vol = 0.;
    for(int i = 0; i < 24; ++i)
    {
      vol += intersection_volume(tets[i], tet);
    }
    return vol;
  }

private:
  static Point3 to_point_ld(const Primal3D::PointType& p)
  {
    return Point3 {static_cast<Real>(p[0]), static_cast<Real>(p[1]), static_cast<Real>(p[2])};
  }

  static Plane make_oriented_plane(const Point3& a,
                                   const Point3& b,
                                   const Point3& c,
                                   const Point3& interior)
  {
    Plane p(Vec3::cross_product(b - a, c - a), a);
    if(p.signedDistance(interior) > 0.L)
    {
      p.flip();
    }
    return p;
  }

  static std::vector<Plane> tet_halfspaces(const TetrahedronType& tet)
  {
    const Point3 p0 = to_point_ld(tet[0]);
    const Point3 p1 = to_point_ld(tet[1]);
    const Point3 p2 = to_point_ld(tet[2]);
    const Point3 p3 = to_point_ld(tet[3]);
    const Point3 interior((p0.array() + p1.array() + p2.array() + p3.array()) / 4.L);

    std::vector<Plane> planes;
    planes.reserve(4);
    planes.push_back(make_oriented_plane(p0, p1, p2, interior));
    planes.push_back(make_oriented_plane(p0, p3, p1, interior));
    planes.push_back(make_oriented_plane(p0, p2, p3, interior));
    planes.push_back(make_oriented_plane(p1, p3, p2, interior));
    return planes;
  }

  static std::vector<Plane> dedup_planes(const std::vector<Plane>& planes, Real eps)
  {
    std::vector<Plane> out;
    out.reserve(planes.size());

    for(const auto& p : planes)
    {
      bool is_dup = false;
      for(const auto& q : out)
      {
        const Real dn = (p.getNormal() - q.getNormal()).norm();
        const Real dd = std::abs(p.getOffset() - q.getOffset());
        if(dn <= eps && dd <= eps)
        {
          is_dup = true;
          break;
        }
      }
      if(!is_dup)
      {
        out.push_back(p);
      }
    }

    return out;
  }

  static bool intersect_three_planes(const Plane& p1,
                                     const Plane& p2,
                                     const Plane& p3,
                                     Point3& out,
                                     Real det_eps)
  {
    const Vec3 n1 = p1.getNormal();
    const Vec3 n2 = p2.getNormal();
    const Vec3 n3 = p3.getNormal();

    const Real det = Vec3::scalar_triple_product(n1, n2, n3);
    if(std::abs(det) <= det_eps)
    {
      return false;
    }

    const Vec3 x = (1 / det) *
      (p1.getOffset() * Vec3::cross_product(n2, n3) +  //
       p2.getOffset() * Vec3::cross_product(n3, n1) +  //
       p3.getOffset() * Vec3::cross_product(n1, n2));
    out = Point3(x.array());
    return true;
  }

  static std::vector<Point3> unique_points(const std::vector<Point3>& points, Real eps)
  {
    std::vector<Point3> out;
    out.reserve(points.size());

    const Real eps2 = eps * eps;
    for(const auto& p : points)
    {
      bool is_dup = false;
      for(const auto& q : out)
      {
        const Vec3 d = p - q;
        if(d.dot(d) <= eps2)
        {
          is_dup = true;
          break;
        }
      }
      if(!is_dup)
      {
        out.push_back(p);
      }
    }

    return out;
  }

  struct HullPt2
  {
    Vec2 p;
    int idx {-1};
  };

  static Real cross2(const HullPt2& a, const HullPt2& b, const HullPt2& c)
  {
    const Real abx = b.p[0] - a.p[0];
    const Real aby = b.p[1] - a.p[1];
    const Real acx = c.p[0] - a.p[0];
    const Real acy = c.p[1] - a.p[1];
    return abx * acy - aby * acx;
  }

  static std::vector<int> convex_hull_2d_indices(std::vector<HullPt2> pts, Real eps)
  {
    if(pts.size() < 3)
    {
      return {};
    }

    std::sort(pts.begin(), pts.end(), [](const HullPt2& a, const HullPt2& b) {
      return (a.p[0] < b.p[0]) || (a.p[0] == b.p[0] && a.p[1] < b.p[1]);
    });

    std::vector<HullPt2> lower;
    for(const auto& p : pts)
    {
      while(lower.size() >= 2 && cross2(lower[lower.size() - 2], lower[lower.size() - 1], p) <= eps)
      {
        lower.pop_back();
      }
      lower.push_back(p);
    }

    std::vector<HullPt2> upper;
    for(int i = static_cast<int>(pts.size()) - 1; i >= 0; --i)
    {
      const auto& p = pts[i];
      while(upper.size() >= 2 && cross2(upper[upper.size() - 2], upper[upper.size() - 1], p) <= eps)
      {
        upper.pop_back();
      }
      upper.push_back(p);
    }

    lower.pop_back();
    upper.pop_back();

    std::vector<int> hull;
    hull.reserve(lower.size() + upper.size());
    for(const auto& p : lower)
    {
      hull.push_back(p.idx);
    }
    for(const auto& p : upper)
    {
      hull.push_back(p.idx);
    }

    return hull;
  }

  static Real volume_from_halfspaces(const std::vector<Plane>& planes,
                                     const std::vector<Point3>& verts)
  {
    if(verts.size() < 4)
    {
      return 0.;
    }

    const Real scale = [&verts]() {
      Real max_abs = 0.;
      for(const auto& v : verts)
      {
        max_abs = std::max(max_abs, std::abs(v[0]));
        max_abs = std::max(max_abs, std::abs(v[1]));
        max_abs = std::max(max_abs, std::abs(v[2]));
      }
      return max_abs;
    }();

    const Real on_plane_tol = 1e-14L * (1.L + scale);
    const Real hull_eps = 1e-18L;

    Real signed_vol = 0.;

    for(const auto& plane : planes)
    {
      Real max_dist = -std::numeric_limits<Real>::infinity();
      for(const auto& v : verts)
      {
        max_dist = std::max(max_dist, plane.signedDistance(v));
      }
      if(max_dist < -on_plane_tol)
      {
        continue;
      }

      std::vector<int> on_plane;
      on_plane.reserve(verts.size());
      for(int i = 0; i < static_cast<int>(verts.size()); ++i)
      {
        if(std::abs(plane.signedDistance(verts[i])) <= on_plane_tol)
        {
          on_plane.push_back(i);
        }
      }
      if(on_plane.size() < 3)
      {
        continue;
      }

      const Vec3 n = plane.getNormal();
      const Vec3 a0 = (std::abs(n[0]) < 0.9L) ? Vec3 {1.L, 0.L, 0.L} : Vec3 {0.L, 1.L, 0.L};
      const Vec3 u = Vec3::cross_product(a0, n).unitVector();
      const Vec3 v = Vec3::cross_product(n, u);

      const Point3 origin = verts[on_plane[0]];
      std::vector<HullPt2> pts2;
      pts2.reserve(on_plane.size());
      for(const int idx : on_plane)
      {
        const Vec3 r = verts[idx] - origin;
        pts2.push_back(HullPt2 {Vec2 {r.dot(u), r.dot(v)}, idx});
      }

      std::vector<int> hull = convex_hull_2d_indices(std::move(pts2), hull_eps);
      if(hull.size() < 3)
      {
        continue;
      }

      const Point3 p0 = verts[hull[0]];
      const Point3 p1 = verts[hull[1]];
      const Point3 p2 = verts[hull[2]];
      const Vec3 face_n = Vec3::cross_product(p1 - p0, p2 - p0);
      if(face_n.dot(n) < 0.L)
      {
        std::reverse(hull.begin(), hull.end());
      }

      const auto f0 = Vec3(verts[hull[0]]);
      for(std::size_t i = 1; i + 1 < hull.size(); ++i)
      {
        signed_vol += f0.dot(Vec3::cross_product(Vec3(verts[hull[i]]), Vec3(verts[hull[i + 1]])));
      }
    }

    return std::abs(signed_vol / 6.L);
  }
};

}  // namespace

TEST(primal_intersection_volume, hex_tet_user_regression_cases)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr double sixth = 1. / 6.;

  HexahedronType hex(PointType {42. + sixth, -66, -178.5},
                     PointType {52. - sixth, -66, -178.5},
                     PointType {52. - sixth, -55, -178.5},
                     PointType {42. + sixth, -55, -178.5},
                     PointType {42. + sixth, -66, -170},
                     PointType {52. - sixth, -66, -170},
                     PointType {52. - sixth, -55, -170},
                     PointType {42. + sixth, -55, -170});

  const TetrahedronType cases[] = {
    TetrahedronType {PointType {27.5859, -19.5363, -148.01},
                     PointType {44.2539, -58.1624, -171.152},
                     PointType {43.7957, -57.9539, -146.494},
                     PointType {15.8564, -32.9522, -147.246}},
    TetrahedronType {PointType {76.8265, -45.3561, -146.396},
                     PointType {79.5055, -43.6324, -171.152},
                     PointType {77.7366, -63.3987, -171.152},
                     PointType {44.2539, -58.1624, -171.152}},
    TetrahedronType {PointType {77.4836, -63.4271, -145.519},
                     PointType {76.8265, -45.3561, -146.396},
                     PointType {44.2539, -58.1624, -171.152},
                     PointType {43.7957, -57.9539, -146.494}},
    TetrahedronType {PointType {77.4836, -63.4271, -145.519},
                     PointType {76.8265, -45.3561, -146.396},
                     PointType {77.7366, -63.3987, -171.152},
                     PointType {44.2539, -58.1624, -171.152}},
  };

  const double expected_volumes[] = {0.30774371225316693,
                                     15.302276033131164,
                                     0.2431953788553831,
                                     3.1668745904354338};

  for(int i = 0; i < 4; ++i)
  {
    const auto& tet = cases[i];
    const TetrahedronType tet_flipped(tet[1], tet[0], tet[2], tet[3]);

    constexpr bool fix_orient = true;
    const double direct_hex_tet =
      axom::primal::intersection_volume<double>(hex, tet, EPS, fix_orient);
    const double ref_subdiv_hex_tet = intersection_volume_via_hex_triangulation(hex, tet, EPS);
    const double ref_halfspace =
      static_cast<double>(HalfspaceIntersectionVolumeLD {}.intersection_volume(hex, tet));

    // check that intersection volumes are non-negative
    EXPECT_GT(direct_hex_tet, 0.) << "case " << i;
    EXPECT_GT(ref_subdiv_hex_tet, 0.) << "case " << i;
    EXPECT_GT(ref_halfspace, 0.) << "case " << i;

    // check that primal::intersection_volume results matches expectations and other ways to compute it
    const double expected = expected_volumes[i];
    EXPECT_NEAR(direct_hex_tet, expected, EPS) << "case " << i;
    EXPECT_NEAR(direct_hex_tet, ref_subdiv_hex_tet, EPS) << "case " << i;
    EXPECT_NEAR(direct_hex_tet, ref_halfspace, EPS) << "case " << i;
  }
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
