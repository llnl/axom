// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/quest/Discretize.hpp"

using SphereType = axom::primal::Sphere<double, 3>;
using OctType = axom::primal::Octahedron<double, 3>;
using Point3D = axom::primal::Point<double, 3>;
using Point2D = axom::primal::Point<double, 2>;
using NAType = axom::NumericArray<double, 3>;

// gtest includes
#include "gtest/gtest.h"

#include <vector>

//------------------------------------------------------------------------------
bool check_generation(axom::Array<OctType>& standard,
                      axom::Array<OctType>& test,
                      int generation,
                      int offset,
                      int count)
{
  std::vector<int> matched(count, 0);

  for(int i = 0; i < count; ++i)
  {
    int has_matched = -1;
    for(int j = 0; j < count; ++j)
    {
      if(!matched[j] && standard[offset + i].equals(test[offset + j]))
      {
        matched[j] = 1;
        has_matched = offset + j;
      }
    }
    if(has_matched < 0)
    {
      std::cout << "Gen " << generation << " standard oct " << (offset + i)
                << " didn't match: " << standard[offset + i] << std::endl;
    }
  }

  int matchcount = 0;
  for(int i = 0; i < count; ++i)
  {
    matchcount += matched[i];
  }

  bool matches = (matchcount == count);

  if(!matches)
  {
    // a little more reporting here
    std::cout << "Generation " << generation << " had " << (count - matchcount)
              << " test octahedra not matched to standard octahedra:" << std::endl;
    for(int i = 0; i < count; ++i)
    {
      if(!matched[i])
      {
        std::cout << "Test oct " << (offset + i) << " not matched:" << test[offset + i] << std::endl;
      }
    }
  }

  return matches;
}

/* Return a handwritten list of the octahedra discretizing a one-segment polyline.
 *
 * The routine allocates and returns an Array of octahedrons, out.
 */
void discretized_segment(Point2D a, Point2D b, axom::Array<OctType>& out)
{
  // We're going to return three generations in the out-vector:
  // one in the first generation, three in the second (covering each of the
  // side-walls of the first generation), eight in the third (covering the
  // two exposed side-walls in each of the second-gen octs).
  //
  // These octahedra are actually triangular prisms (if the polyline
  // is parallel to the X-axis) or truncated tetrahedra.  They have two parallel
  // triangle end-caps and three quadrilateral side-walls.  Because we store
  // them in Octahedron objects, those side-walls are actually two triangles
  // that are coplanar.
  constexpr int ZEROTH_GEN_COUNT = 1;
  constexpr int FIRST_GEN_COUNT = 3;
  constexpr int SECOND_GEN_COUNT = 6;
  out = axom::Array<OctType>(ZEROTH_GEN_COUNT + FIRST_GEN_COUNT + SECOND_GEN_COUNT);

  // The first generation puts a triangle in the end-discs of the truncated
  // cones with vertices at 12, 4, and 8 o'clock.
  // The second generation puts three triangles in the end-discs;
  //     the first has vertices at 12, 2, and 4 o'clock,
  //     the second at 4, 6, and 8 o'clock,
  //     the third at 8, 10, and 12 o'clock.
  // The third generation puts six more triangles in the ends;
  //     first at 12, 1, 2 o'clock,
  //     second at 2, 3, 4 o'clock,
  //     third at 4, 5, 6, o'clock,
  //     fourth at 6, 7, 8 o'clock,
  //     fifth at 8, 9, 10 o'clock,
  //     sixth at 10, 11, 12 o'clock.
  // These are constructed with 1, 1/2, and sqrt(3)/2 in discs perpendicular
  // to the X-axis at a[0] and b[0].  We're looking out the X-axis so a
  // clockwise spin corresponds to the right-hand spin in the function
  // from_segment() in Discretize.cpp.

  const double SQ_3_2 = sqrt(3.) / 2.;

  // clang-format off
  Point3D Ia   {a[0],    -0.5*a[1],  SQ_3_2*a[1]};
  Point3D IIa  {a[0], -SQ_3_2*a[1],     0.5*a[1]};
  Point3D IIIa {a[0],      -1*a[1],       0*a[1]};
  Point3D IVa  {a[0], -SQ_3_2*a[1],    -0.5*a[1]};
  Point3D Va   {a[0],    -0.5*a[1], -SQ_3_2*a[1]};
  Point3D VIa  {a[0],       0*a[1],      -1*a[1]};
  Point3D VIIa {a[0],     0.5*a[1], -SQ_3_2*a[1]};
  Point3D VIIIa{a[0],  SQ_3_2*a[1],    -0.5*a[1]};
  Point3D IXa  {a[0],       1*a[1],       0*a[1]};
  Point3D Xa   {a[0],  SQ_3_2*a[1],     0.5*a[1]};
  Point3D XIa  {a[0],     0.5*a[1],  SQ_3_2*a[1]};
  Point3D XIIa {a[0],       0*a[1],       1*a[1]};

  Point3D Ib   {b[0],    -0.5*b[1],  SQ_3_2*b[1]};
  Point3D IIb  {b[0], -SQ_3_2*b[1],     0.5*b[1]};
  Point3D IIIb {b[0],      -1*b[1],       0*b[1]};
  Point3D IVb  {b[0], -SQ_3_2*b[1],    -0.5*b[1]};
  Point3D Vb   {b[0],    -0.5*b[1], -SQ_3_2*b[1]};
  Point3D VIb  {b[0],       0*b[1],      -1*b[1]};
  Point3D VIIb {b[0],     0.5*b[1], -SQ_3_2*b[1]};
  Point3D VIIIb{b[0],  SQ_3_2*b[1],    -0.5*b[1]};
  Point3D IXb  {b[0],       1*b[1],       0*b[1]};
  Point3D Xb   {b[0],  SQ_3_2*b[1],     0.5*b[1]};
  Point3D XIb  {b[0],     0.5*b[1],  SQ_3_2*b[1]};
  Point3D XIIb {b[0],       0*b[1],       1*b[1]};
  // clang-format on

  // First generation: a single prism
  out[0] = OctType(XIIa, IVb, IVa, VIIIb, VIIIa, XIIb);

  // Second generation: a prism on the outside of each side wall.
  out[1] = OctType(IVb, VIa, VIb, VIIIa, VIIIb, IVa);
  out[2] = OctType(XIIb, IIa, IIb, IVa, IVb, XIIa);
  out[3] = OctType(VIIIb, Xa, Xb, XIIa, XIIb, VIIIa);

  // Third generation: a prism on the outside of each exterior side wall.
  out[4] = OctType(XIIa, Ib, Ia, IIb, IIa, XIIb);
  out[5] = OctType(IIa, IIIb, IIIa, IVb, IVa, IIb);
  out[6] = OctType(IVa, Vb, Va, VIb, VIa, IVb);
  out[7] = OctType(VIa, VIIb, VIIa, VIIIb, VIIIa, VIb);
  out[8] = OctType(VIIIa, IXb, IXa, Xb, Xa, VIIIb);
  out[9] = OctType(Xa, XIb, XIa, XIIb, XIIa, Xb);
}

template <typename ExecPolicy>
void degenerate_segment_test(const char* label, axom::Array<Point2D>& polyline, int len, bool expsuccess)
{
  constexpr int gens = 2;
  int octcount = 0;

  SCOPED_TRACE(label);

  axom::Array<OctType> generated;
  if(expsuccess)
  {
    EXPECT_TRUE(axom::quest::discretize<ExecPolicy>(polyline, len, gens, generated, octcount));
  }
  else
  {
    EXPECT_FALSE(axom::quest::discretize<ExecPolicy>(polyline, len, gens, generated, octcount));
  }
  EXPECT_EQ(0, octcount);
}

//------------------------------------------------------------------------------
template <typename ExecPolicy>
void run_degen_segment_tests()
{
  // Test each of the three generations.
  // We don't know what order they'll be in, but we do know how many octahedra
  // will be in each generation.

  int allocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Use unified memory for frequent movement between device operations
  // and value checking on host
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecPolicy>::onDevice())
  {
    allocID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  }
#endif

  axom::Array<Point2D> polyline(2, 2, allocID);

  polyline[0] = {0., 0.};
  polyline[1] = {0., 0.};
  degenerate_segment_test<ExecPolicy>("a.x == b.x, a.y == b.y", polyline, 2, true);

  polyline[0] = {1., 0.};
  polyline[1] = {1., 1.};
  degenerate_segment_test<ExecPolicy>("a.x == b.x, a.y != b.y", polyline, 2, true);

  polyline[0] = {1., -0.1};
  polyline[1] = {1.5, 1.};
  degenerate_segment_test<ExecPolicy>("a.y < 0", polyline, 2, false);

  polyline[0] = {1., -0.1};
  polyline[1] = {1.5, 1.};
  degenerate_segment_test<ExecPolicy>("a.y < 0", polyline, 2, false);

  polyline[0] = {1., 1.};
  polyline[1] = {1.5, -.1};
  degenerate_segment_test<ExecPolicy>("b.y < 0", polyline, 2, false);

  polyline[0] = {.5, 1.};
  polyline[1] = {0., 1.};
  degenerate_segment_test<ExecPolicy>("a.x > b.x", polyline, 2, false);
}

template <typename ExecPolicy>
void segment_test(const char* label, axom::Array<Point2D>& polyline, int len)
{
  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  SCOPED_TRACE(label);

  // Test each of the three generations.
  constexpr int generations = 2;

  // We don't know what order they'll be in, but we do know how many octahedra
  // will be in each generation.
  constexpr int ZEROTH_GEN_COUNT = 1;
  constexpr int FIRST_GEN_COUNT = 3;
  constexpr int SECOND_GEN_COUNT = 6;

  int generation = 0;

  // The discretized_segment() routine produces a list of 10 hand-calculated
  // octahedra (three generations) that discretize the surface of revolution
  // (SoR) produced by revolving a one-segment polyline around the positive
  // X-axis.
  axom::Array<OctType> handcut;
  discretized_segment(polyline[0], polyline[1], handcut);

  axom::Array<OctType> generatedDevice;
  int octcount = 0;
  axom::quest::discretize<ExecPolicy>(polyline, len, generations, generatedDevice, octcount);

  // Copy generated back to host
  axom::Array<OctType> generated = axom::Array<OctType>(generatedDevice, hostAllocID);

  EXPECT_TRUE(check_generation(handcut, generated, generation, 0, ZEROTH_GEN_COUNT));
  generation += 1;
  EXPECT_TRUE(check_generation(handcut, generated, generation, ZEROTH_GEN_COUNT, FIRST_GEN_COUNT));
  generation += 1;
  EXPECT_TRUE(check_generation(handcut,
                               generated,
                               generation,
                               ZEROTH_GEN_COUNT + FIRST_GEN_COUNT,
                               SECOND_GEN_COUNT));
}

//------------------------------------------------------------------------------
template <typename ExecPolicy>
void run_single_segment_tests()
{
  int allocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Use unified memory for frequent movement between device operations
  // and value checking on host
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecPolicy>::onDevice())
  {
    allocID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  }
#endif

  axom::Array<Point2D> polyline(2, 2, allocID);

  polyline[0] = Point2D {0.5, 0.};
  polyline[1] = Point2D {1.8, 0.8};
  segment_test<ExecPolicy>("Cone (pointing left)", polyline, 2);

  polyline[0] = Point2D {1.0, 1.0};
  polyline[1] = Point2D {1.8, 0.};
  segment_test<ExecPolicy>("Cone (pointing right)", polyline, 2);

  polyline[0] = Point2D {1.0, 1.0};
  polyline[1] = Point2D {1.8, 0.8};
  segment_test<ExecPolicy>("Truncated cone", polyline, 2);

  polyline[0] = Point2D {1., 1.};
  polyline[1] = Point2D {3., 1.};
  segment_test<ExecPolicy>("Cylinder", polyline, 2);

  polyline[0] = Point2D {-.4, 1.2};
  polyline[1] = Point2D {1.2, 1.};
  segment_test<ExecPolicy>("a.x < 0, b.x > 0", polyline, 2);
}

template <typename ExecPolicy>
void multi_segment_test(const char* label, axom::Array<Point2D>& polyline, int len)
{
  SCOPED_TRACE(label);

  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Test each of the three generations.
  constexpr int generations = 2;

  // We don't know what order they'll be in, but we do know how many octahedra
  // will be in each generation.
  constexpr int ZEROTH_GEN_COUNT = 1;
  constexpr int FIRST_GEN_COUNT = 3;
  constexpr int SECOND_GEN_COUNT = 6;
  constexpr int TOTAL_COUNT = ZEROTH_GEN_COUNT + FIRST_GEN_COUNT + SECOND_GEN_COUNT;

  int generation = 0;

  axom::Array<OctType> generatedDevice;
  int octcount = 0;
  axom::quest::discretize<ExecPolicy>(polyline, len, generations, generatedDevice, octcount);

  // Copy generated back to host
  axom::Array<OctType> generated = axom::Array<OctType>(generatedDevice, hostAllocID);

  int segcount = len - 1;
  axom::Array<OctType> handcut(segcount * TOTAL_COUNT);

  for(int segidx = 0; segidx < segcount; ++segidx)
  {
    std::stringstream seglabel;
    seglabel << "Segment " << segidx;
    SCOPED_TRACE(seglabel.str());

    // The discretized_segment() routine produces a list of 10 hand-calculated
    // octahedra (three generations) that discretize the surface of revolution
    // (SoR) produced by revolving a one-segment polyline around the positive
    // X-axis.
    axom::Array<OctType> octpointer;
    discretized_segment(polyline[segidx], polyline[segidx + 1], octpointer);

    EXPECT_TRUE(
      check_generation(octpointer, generated, generation, segidx * TOTAL_COUNT + 0, ZEROTH_GEN_COUNT));
    generation += 1;
    EXPECT_TRUE(check_generation(octpointer,
                                 generated,
                                 generation,
                                 segidx * TOTAL_COUNT + ZEROTH_GEN_COUNT,
                                 FIRST_GEN_COUNT));
    generation += 1;
    EXPECT_TRUE(check_generation(octpointer,
                                 generated,
                                 generation,
                                 segidx * TOTAL_COUNT + ZEROTH_GEN_COUNT + FIRST_GEN_COUNT,
                                 SECOND_GEN_COUNT));
  }
}

//------------------------------------------------------------------------------
template <typename ExecPolicy>
void run_multi_segment_tests()
{
  int allocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  //Use unified memory for frequent movement between device operations
  // and value checking on host
#if defined(AXOM_USE_GPU) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecPolicy>::onDevice())
  {
    allocID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  }
#endif

  constexpr int pointcount = 5;
  axom::Array<Point2D> polyline(pointcount, pointcount, allocID);

  polyline[0] = Point2D {1.0, 0.5};
  polyline[1] = Point2D {1.6, 0.3};
  polyline[2] = Point2D {1.7, 2.0};
  polyline[3] = Point2D {3.1, 1.5};
  polyline[4] = Point2D {4.0, 2.0};
  segment_test<ExecPolicy>("multi-segment", polyline, pointcount);
}

//------------------------------------------------------------------------------
enum ReflectDimension
{
  X = 0,
  Y,
  Z
};

/* Reflect an octahedron in a given dimension.
 */
OctType reflect(ReflectDimension d, OctType o)
{
  OctType out(o);
  for(int i = 0; i < OctType::NUM_VERTS; ++i)
  {
    Point3D& pt = out[i];
    pt[d] *= -1;
  }
  return out;
}

/* Return a handwritten list of the octahedra discretizing the unit sphere.
 *
 * The routine allocates and returns an Array of octahedrons, out.
 */
void discretized_sphere(axom::Array<OctType>& out)
{
  // We're going to return three generations in the out-vector:
  // one in the first generation, eight in the second (covering each
  // of the faces of the first generation), 32 in the third (covering
  // the four exposed faces in each of the second-gen octs).
  constexpr int ZEROTH_GEN_COUNT = 1;
  constexpr int FIRST_GEN_COUNT = 8;
  constexpr int SECOND_GEN_COUNT = 32;
  out = axom::Array<OctType>(ZEROTH_GEN_COUNT + FIRST_GEN_COUNT + SECOND_GEN_COUNT);

  // First generation: one octahedron, with vertices on the unit vectors.
  NAType ihat({1., 0., 0.});
  NAType jhat({0., 1., 0.});
  NAType khat({0., 0., 1.});

  OctType center(Point3D(ihat),
                 Point3D(jhat),
                 Point3D(khat),
                 Point3D(-1 * ihat),
                 Point3D(-1 * jhat),
                 Point3D(-1 * khat));
  out[0] = center;

  // Second generation: eight octs, one for each face of the unit oct.
  // Point ij is halfway between (1, 0, 0) and (0, 1, 0),
  // point jk is halfway between (0, 1, 0) and (0, 0, 1),
  // point ki is halfway between (0, 0, 1) and (1, 0, 0).
  Point3D ij(NAType({M_SQRT1_2, M_SQRT1_2, 0.}));
  Point3D jk(NAType({0., M_SQRT1_2, M_SQRT1_2}));
  Point3D ki(NAType({M_SQRT1_2, 0., M_SQRT1_2}));

  OctType second_gen(Point3D(ihat), Point3D(jhat), Point3D(khat), jk, ki, ij);
  out[1] = second_gen;
  out[2] = reflect(X, second_gen);
  out[3] = reflect(Y, second_gen);
  out[4] = reflect(Z, second_gen);
  out[5] = reflect(Y, reflect(X, second_gen));
  out[6] = reflect(Z, reflect(Y, second_gen));
  out[7] = reflect(X, reflect(Z, second_gen));
  out[8] = reflect(X, reflect(Y, reflect(Z, second_gen)));

  // Third generation: 32 new octs, one for each exposed face of the previous
  // generation.
  double SQRT1_6 = 1. / sqrt(6.);
  // There are three interior points, derived from ij, jk, and ki.
  // Point a is halfway between ij and ki, at (1/sqrt(6))(2, 1, 1).
  Point3D a(SQRT1_6 * NAType({2, 1, 1}));
  // Point b is halfway between ij and jk, at (1/sqrt(6))(1, 2, 1).
  Point3D b(SQRT1_6 * NAType({1, 2, 1}));
  // Point c is halfway between jk and ki, at (1/sqrt(6))(1, 1, 2).
  Point3D c(SQRT1_6 * NAType({1, 1, 2}));

  // There are six edge points, derived from the original corner points and
  // ij, jk, and ki.
  // Point d is halfway between ihat and ij, at
  // (1/sqrt(4 + 2 sqrt(2)))(1+sqrt(2), 1, 0)
  double FACTOR_3G = 1. / sqrt(2. * M_SQRT2 + 4);
  Point3D d(FACTOR_3G * NAType({1. + M_SQRT2, 1., 0}));
  // Point e splits jhat and ij, at (1/sqrt(2 sqrt(2) + 2))(1, 1+sqrt(2), 0)
  Point3D e(FACTOR_3G * NAType({1., 1. + M_SQRT2, 0}));
  // Point f splits jhat and jk, at (1/sqrt(2 sqrt(2) + 2))(0, 1+sqrt(2), 1)
  Point3D f(FACTOR_3G * NAType({0, 1. + M_SQRT2, 1.}));
  // Point g splits khat and jk, at (1/sqrt(2 sqrt(2) + 2))(0, 1, 1+sqrt(2))
  Point3D g(FACTOR_3G * NAType({0, 1., 1. + M_SQRT2}));
  // Point m splits khat and ki, at (1/sqrt(2 sqrt(2) + 2))(1, 0, 1+sqrt(2))
  Point3D m(FACTOR_3G * NAType({1., 0, 1. + M_SQRT2}));
  // Point n splits ihat and ki, at (1/sqrt(2 sqrt(2) + 2))(1+sqrt(2), 0, 1)
  Point3D n(FACTOR_3G * NAType({1. + M_SQRT2, 0, 1.}));

  int offset = ZEROTH_GEN_COUNT + FIRST_GEN_COUNT;
  // Here's the first octant of third-generation octahedra (octant 0).
  int octant = 0;
  // First, the interior oct
  out[offset + 0 + octant * 4] = OctType(ij, jk, ki, c, a, b);
  // The one next to ihat
  out[offset + 1 + octant * 4] = OctType(Point3D(ihat), ij, ki, a, n, d);
  // The one next to jhat
  out[offset + 2 + octant * 4] = OctType(Point3D(jhat), jk, ij, b, e, f);
  // The one next to khat
  out[offset + 3 + octant * 4] = OctType(Point3D(khat), ki, jk, c, g, m);

  // Now we get to transform these into all seven remaining octants.
  // Reflect in X
  octant = 1;
  ReflectDimension rd0 = X;
  out[offset + 0 + octant * 4] = reflect(rd0, out[offset + 0]);
  out[offset + 1 + octant * 4] = reflect(rd0, out[offset + 1]);
  out[offset + 2 + octant * 4] = reflect(rd0, out[offset + 2]);
  out[offset + 3 + octant * 4] = reflect(rd0, out[offset + 3]);
  // Reflect in Y
  octant = 2;
  rd0 = Y;
  out[offset + 0 + octant * 4] = reflect(rd0, out[offset + 0]);
  out[offset + 1 + octant * 4] = reflect(rd0, out[offset + 1]);
  out[offset + 2 + octant * 4] = reflect(rd0, out[offset + 2]);
  out[offset + 3 + octant * 4] = reflect(rd0, out[offset + 3]);
  // Reflect in Z
  octant = 3;
  rd0 = Z;
  out[offset + 0 + octant * 4] = reflect(rd0, out[offset + 0]);
  out[offset + 1 + octant * 4] = reflect(rd0, out[offset + 1]);
  out[offset + 2 + octant * 4] = reflect(rd0, out[offset + 2]);
  out[offset + 3 + octant * 4] = reflect(rd0, out[offset + 3]);
  // Reflect in X, then Y
  octant = 4;
  rd0 = X;
  ReflectDimension rd1 = Y;
  out[offset + 0 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 0]));
  out[offset + 1 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 1]));
  out[offset + 2 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 2]));
  out[offset + 3 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 3]));
  // Reflect in Y, then Z
  octant = 5;
  rd0 = Y;
  rd1 = Z;
  out[offset + 0 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 0]));
  out[offset + 1 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 1]));
  out[offset + 2 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 2]));
  out[offset + 3 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 3]));
  // Reflect in Z, then X
  octant = 6;
  rd0 = Z;
  rd1 = X;
  out[offset + 0 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 0]));
  out[offset + 1 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 1]));
  out[offset + 2 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 2]));
  out[offset + 3 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 3]));
  // And the grand finale: reflect in Z, then Y, then X
  octant = 7;
  rd0 = Z;
  rd1 = Y;
  ReflectDimension rd2 = X;
  out[offset + 0 + octant * 4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset + 0])));
  out[offset + 1 + octant * 4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset + 1])));
  out[offset + 2 + octant * 4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset + 2])));
  out[offset + 3 + octant * 4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset + 3])));
}

//------------------------------------------------------------------------------
TEST(quest_discretize, sphere_test)
{
  // The discretize() routine chops up a given sphere into the specified
  // number of generations of octahedra.  Here, we'll discretize the unit
  // sphere in three generations, to match the hand-calculated octs from
  // discretized_sphere().
  SphereType sph;  // Unit sphere at the origin
  constexpr int generations = 2;
  axom::Array<OctType> generated;
  int octcount = 0;
  axom::quest::discretize(sph, generations, generated, octcount);

  // The discretized_sphere() routine produces a list of 41 hand-calculated
  // octahedra (three generations) that discretize the unit sphere.
  axom::Array<OctType> handcut;
  discretized_sphere(handcut);

  // Test each of the three generations.
  // We don't know what order they'll be in, but we do know how many octahedra
  // will be in each generation.
  constexpr int ZEROTH_GEN_COUNT = 1;
  constexpr int FIRST_GEN_COUNT = 8;
  constexpr int SECOND_GEN_COUNT = 32;
  int generation = 0;

  EXPECT_TRUE(check_generation(handcut, generated, generation, 0, ZEROTH_GEN_COUNT));
  generation += 1;
  EXPECT_TRUE(check_generation(handcut, generated, generation, ZEROTH_GEN_COUNT, FIRST_GEN_COUNT));
  generation += 1;
  EXPECT_TRUE(check_generation(handcut,
                               generated,
                               generation,
                               ZEROTH_GEN_COUNT + FIRST_GEN_COUNT,
                               SECOND_GEN_COUNT));
}

//------------------------------------------------------------------------------
TEST(quest_discretize, degenerate_sphere_test)
{
  constexpr int generations = 2;
  axom::Array<OctType> generated;
  int octcount = 0;

  {
    SCOPED_TRACE("Negative sphere radius");
    SphereType sph(-.2);  // BAD sphere at the origin with negative radius
    EXPECT_FALSE(axom::quest::discretize(sph, generations, generated, octcount));
    EXPECT_EQ(0, octcount);
  }

  {
    SCOPED_TRACE("Zero sphere radius");
    SphereType sph(0.);  // Degenerate sphere at the origin with zero radius
    EXPECT_TRUE(axom::quest::discretize(sph, generations, generated, octcount));
    EXPECT_EQ(0, octcount);
  }
}

//------------------------------------------------------------------------------
TEST(quest_discretize, degenerate_segment_test)
{
  SLIC_INFO("Discretizing sequentially");
  {
    SCOPED_TRACE("sequential execution");
    run_degen_segment_tests<axom::SEQ_EXEC>();
  }

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  SLIC_INFO("Discretizing with OpenMP");
  {
    SCOPED_TRACE("OpenMP execution");
    run_degen_segment_tests<axom::OMP_EXEC>();
  }
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__CUDACC__)
  SLIC_INFO("Discretizing with CUDA");
  {
    SCOPED_TRACE("32-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<32>>();
  }

  {
    SCOPED_TRACE("64-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<64>>();
  }

  {
    SCOPED_TRACE("128-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<128>>();
  }

  {
    SCOPED_TRACE("256-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<256>>();
  }

  {
    SCOPED_TRACE("512-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<512>>();
  }
#endif

#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__HIPCC__)
  SLIC_INFO("Discretizing with HIP");
  {
    SCOPED_TRACE("32-wide HIP execution");
    run_degen_segment_tests<axom::HIP_EXEC<32>>();
  }

  {
    SCOPED_TRACE("64-wide HIP execution");
    run_degen_segment_tests<axom::HIP_EXEC<64>>();
  }

  {
    SCOPED_TRACE("128-wide HIP execution");
    run_degen_segment_tests<axom::HIP_EXEC<128>>();
  }

  {
    SCOPED_TRACE("256-wide HIP execution");
    run_degen_segment_tests<axom::HIP_EXEC<256>>();
  }

  {
    SCOPED_TRACE("512-wide HIP execution");
    run_degen_segment_tests<axom::HIP_EXEC<512>>();
  }
#endif
}

//------------------------------------------------------------------------------
TEST(quest_discretize, segment_test)
{
  SLIC_INFO("Discretizing sequentially");
  {
    SCOPED_TRACE("sequential execution");
    run_single_segment_tests<axom::SEQ_EXEC>();
  }

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  SLIC_INFO("Discretizing with OpenMP");
  {
    SCOPED_TRACE("OpenMP execution");
    run_single_segment_tests<axom::OMP_EXEC>();
  }
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__CUDACC__)
  SLIC_INFO("Discretizing with CUDA");
  {
    SCOPED_TRACE("32-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<32>>();
  }

  {
    SCOPED_TRACE("64-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<64>>();
  }

  {
    SCOPED_TRACE("128-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<128>>();
  }

  {
    SCOPED_TRACE("256-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<256>>();
  }

  {
    SCOPED_TRACE("512-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<512>>();
  }
#endif

#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__HIPCC__)
  SLIC_INFO("Discretizing with HIP");
  {
    SCOPED_TRACE("32-wide HIP execution");
    run_single_segment_tests<axom::HIP_EXEC<32>>();
  }

  {
    SCOPED_TRACE("64-wide HIP execution");
    run_single_segment_tests<axom::HIP_EXEC<64>>();
  }

  {
    SCOPED_TRACE("128-wide HIP execution");
    run_single_segment_tests<axom::HIP_EXEC<128>>();
  }

  {
    SCOPED_TRACE("256-wide HIP execution");
    run_single_segment_tests<axom::HIP_EXEC<256>>();
  }

  {
    SCOPED_TRACE("512-wide HIP execution");
    run_single_segment_tests<axom::HIP_EXEC<512>>();
  }
#endif
}

//------------------------------------------------------------------------------
TEST(quest_discretize, multi_segment_test)
{
  SLIC_INFO("Discretizing sequentially");
  {
    SCOPED_TRACE("sequential execution");
    run_multi_segment_tests<axom::SEQ_EXEC>();
  }

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  SLIC_INFO("Discretizing with OpenMP");
  {
    SCOPED_TRACE("OpenMP execution");
    run_multi_segment_tests<axom::OMP_EXEC>();
  }
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__CUDACC__)
  SLIC_INFO("Discretizing with CUDA");
  {
    SCOPED_TRACE("32-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<32>>();
  }

  {
    SCOPED_TRACE("64-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<64>>();
  }

  {
    SCOPED_TRACE("128-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<128>>();
  }

  {
    SCOPED_TRACE("256-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<256>>();
  }

  {
    SCOPED_TRACE("512-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<512>>();
  }
#endif

#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__HIPCC__)
  SLIC_INFO("Discretizing with HIP");
  {
    SCOPED_TRACE("32-wide HIP execution");
    run_multi_segment_tests<axom::HIP_EXEC<32>>();
  }

  {
    SCOPED_TRACE("64-wide HIP execution");
    run_multi_segment_tests<axom::HIP_EXEC<64>>();
  }

  {
    SCOPED_TRACE("128-wide HIP execution");
    run_multi_segment_tests<axom::HIP_EXEC<128>>();
  }

  {
    SCOPED_TRACE("256-wide HIP execution");
    run_multi_segment_tests<axom::HIP_EXEC<256>>();
  }

  {
    SCOPED_TRACE("512-wide HIP execution");
    run_multi_segment_tests<axom::HIP_EXEC<512>>();
  }
#endif
}

//------------------------------------------------------------------------------
TEST(quest_discretize, to_tet_mesh)
{
  constexpr int pointcount = 5;
  constexpr int segcount = pointcount - 1;
  constexpr int generations = 2;
  axom::Array<Point2D> polyline(pointcount);

  polyline[0] = Point2D {1.0, 0.5};
  polyline[1] = Point2D {1.6, 0.3};
  polyline[2] = Point2D {1.7, 2.0};
  polyline[3] = Point2D {3.1, 1.5};
  polyline[4] = Point2D {4.0, 2.0};

  axom::Array<OctType> generated;
  int octcount = 0;
  axom::quest::discretize<axom::SEQ_EXEC>(polyline, pointcount, generations, generated, octcount);

  axom::ArrayView<OctType> generated_view = generated.view();
  axom::mint::Mesh* mesh;
  axom::quest::mesh_from_discretized_polyline(generated_view, octcount, segcount, mesh);
  axom::mint::write_vtk(mesh, "tet_mesh.vtk");

  delete mesh;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
