// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file analytic_test_surfaces.hpp
 *
 * \brief Consists of methods that generate arrays of NURBSPatch objects which depict
 *         objects whose containment decisions can be made analytically, primarily
 *         for the purpose of testing GWN containment methods
 *
 * \sa primal_solid_angle.cpp
 */

#ifndef AXOM_PRIMAL_ANALYTIC_TEST_SURFACES_HPP
#define AXOM_PRIMAL_ANALYTIC_TEST_SURFACES_HPP

#include "axom/config.hpp"
#include "axom/primal.hpp"

namespace axom
{
namespace primal
{
namespace detail
{

axom::Array<primal::NURBSPatch<double, 3>> make_sphere_biquartic()
{
  using Point3D = primal::Point<double, 3>;
  using NPatch = primal::NURBSPatch<double, 3>;

  double rt2 = sqrt(2), rt3 = sqrt(3), rt6 = sqrt(6);

  // Define the nodes and weights for one of six rational, biquartic Bezier patches
  //  that compose the unit sphere. These will be rotated to form the other 5.
  // Nodes and weights obtained from the technical report
  // "Tiling the Sphere with Rational Bezier Patches",
  //  James E. Cobb, University of Utah, 1988

  // clang-format off
    axom::Array<Point3D> node_data = {
      Point3D {4*(1-rt3),     4*(1-rt3),     4*(1-rt3)}, Point3D {rt2*(rt3-4),            -rt2, rt2*(rt3-4)}, Point3D {4*(1-2*rt3)/3,   0, 4*(1-2*rt3)/3}, Point3D {rt2*(rt3-4),           rt2,   rt2*(rt3-4)}, Point3D {4*(1-rt3),     4*(rt3-1),     4*(1-rt3)},
      Point3D {     -rt2, rt2*(rt3 - 4), rt2*(rt3 - 4)}, Point3D {(2-3*rt3)/2,     (2-3*rt3)/2,  -(rt3+6)/2}, Point3D {rt2*(2*rt3-7)/3, 0,      -5*rt6/3}, Point3D {(2-3*rt3)/2,   (3*rt3-2)/2,    -(rt3+6)/2}, Point3D {     -rt2,   rt2*(4-rt3),   rt2*(rt3-4)},
      Point3D {        0, 4*(1-2*rt3)/3, 4*(1-2*rt3)/3}, Point3D {          0, rt2*(2*rt3-7)/3,    -5*rt6/3}, Point3D {0,               0,   4*(rt3-5)/3}, Point3D {          0, rt2*(7-2*rt3)/3,    -5*rt6/3}, Point3D {        0, 4*(2*rt3-1)/3, 4*(1-2*rt3)/3},
      Point3D {      rt2, rt2*(rt3 - 4), rt2*(rt3 - 4)}, Point3D {(3*rt3-2)/2,     (2-3*rt3)/2,  -(rt3+6)/2}, Point3D {rt2*(7-2*rt3)/3, 0,      -5*rt6/3}, Point3D {(3*rt3-2)/2,   (3*rt3-2)/2,    -(rt3+6)/2}, Point3D {      rt2,   rt2*(4-rt3),   rt2*(rt3-4)},
      Point3D {4*(rt3-1),     4*(1-rt3),     4*(1-rt3)}, Point3D {rt2*(4-rt3),            -rt2, rt2*(rt3-4)}, Point3D {4*(2*rt3-1)/3,   0, 4*(1-2*rt3)/3}, Point3D {rt2*(4-rt3),           rt2,   rt2*(rt3-4)}, Point3D {4*(rt3-1),     4*(rt3-1),     4*(1-rt3)}};
  
    axom::Array<double> weight_data = {
           4*(3-rt3), rt2*(3*rt3-2),   4*(5-rt3)/3, rt2*(3*rt3-2),     4*(3-rt3),
       rt2*(3*rt3-2),     (rt3+6)/2, rt2*(rt3+6)/3,     (rt3+6)/2, rt2*(3*rt3-2),
         4*(5-rt3)/3, rt2*(rt3+6)/3, 4*(5*rt3-1)/9, rt2*(rt3+6)/3,   4*(5-rt3)/3,
       rt2*(3*rt3-2),     (rt3+6)/2, rt2*(rt3+6)/3,     (rt3+6)/2, rt2*(3*rt3-2),
           4*(3-rt3), rt2*(3*rt3-2),   4*(5-rt3)/3, rt2*(3*rt3-2),     4*(3-rt3)};
  // clang-format on

  axom::Array<NPatch> sphere_faces(6);
  for(int n = 0; n < 6; ++n)
  {
    sphere_faces[n].setParameters(5, 5, 4, 4);
    sphere_faces[n].makeRational();
  }

  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      const int idx = 5 * i + j;
      for(int n = 0; n < 6; ++n)
      {
        sphere_faces[n].setWeight(i, j, weight_data[idx]);
      }

      // Set up each face by rotating one of the patch faces
      sphere_faces[0](i, j)[0] = node_data[idx][1];
      sphere_faces[0](i, j)[1] = node_data[idx][0];
      sphere_faces[0](i, j)[2] = node_data[idx][2];
      sphere_faces[0](i, j).array() /= weight_data[idx];

      sphere_faces[1](i, j)[0] = -node_data[idx][0];
      sphere_faces[1](i, j)[1] = -node_data[idx][1];
      sphere_faces[1](i, j)[2] = -node_data[idx][2];
      sphere_faces[1](i, j).array() /= weight_data[idx];

      sphere_faces[2](i, j)[0] = node_data[idx][2];
      sphere_faces[2](i, j)[1] = node_data[idx][1];
      sphere_faces[2](i, j)[2] = node_data[idx][0];
      sphere_faces[2](i, j).array() /= weight_data[idx];

      sphere_faces[3](i, j)[0] = -node_data[idx][1];
      sphere_faces[3](i, j)[1] = -node_data[idx][2];
      sphere_faces[3](i, j)[2] = -node_data[idx][0];
      sphere_faces[3](i, j).array() /= weight_data[idx];

      sphere_faces[4](i, j)[0] = node_data[idx][0];
      sphere_faces[4](i, j)[1] = node_data[idx][2];
      sphere_faces[4](i, j)[2] = node_data[idx][1];
      sphere_faces[4](i, j).array() /= weight_data[idx];

      sphere_faces[5](i, j)[0] = -node_data[idx][2];
      sphere_faces[5](i, j)[1] = -node_data[idx][0];
      sphere_faces[5](i, j)[2] = -node_data[idx][1];
      sphere_faces[5](i, j).array() /= weight_data[idx];
    }
  }

  return sphere_faces;
}

axom::Array<primal::NURBSPatch<double, 3>> make_sphere_bicubic()
{
  // Generate a sphere using (degenerate) bicubic patches

  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using BCurve = primal::BezierCurve<double, 2>;
  using NPatch = primal::NURBSPatch<double, 3>;

  double rt2 = sqrt(2);

  // Define BezierCurves which, when rotated around the z-axis,
  //  form the top and bottom halves of a sphere.
  BCurve top_curve(2), bottom_curve(2);
  top_curve[0] = Point2D {0.0, 1.0};
  top_curve[1] = Point2D {1.0, 1.0};
  top_curve[2] = Point2D {1.0, 0.0};

  bottom_curve[0] = Point2D {1.0, 0.0};
  bottom_curve[1] = Point2D {1.0, -1.0};
  bottom_curve[2] = Point2D {0.0, -1.0};

  top_curve.makeRational();
  top_curve.setWeight(1, 1.0 / rt2);

  bottom_curve.makeRational();
  bottom_curve.setWeight(1, 1.0 / rt2);

  axom::Array<NPatch> sphere_faces(8);
  for(int n = 0; n < 8; ++n)
  {
    sphere_faces[n].setParameters(3, 3, 2, 2);
    sphere_faces[n].makeRational();
  }

  for(int n = 0; n < 2; ++n)
  {
    auto& curve = (n == 0) ? top_curve : bottom_curve;

    for(int i = 0; i <= 2; ++i)
    {
      // clang-format off
      sphere_faces[4*n + 0](i, 0) = Point3D {curve[i][0], 0.0, curve[i][1]};
      sphere_faces[4*n + 0](i, 1) = Point3D {curve[i][0], curve[i][0], curve[i][1]};
      sphere_faces[4*n + 0](i, 2) = Point3D {0.0, curve[i][0], curve[i][1]};

      sphere_faces[4*n + 1](i, 0) = Point3D {0.0, curve[i][0], curve[i][1]};
      sphere_faces[4*n + 1](i, 1) = Point3D {-curve[i][0], curve[i][0], curve[i][1]};
      sphere_faces[4*n + 1](i, 2) = Point3D {-curve[i][0], 0.0, curve[i][1]};

      sphere_faces[4*n + 2](i, 0) = Point3D {-curve[i][0], 0.0, curve[i][1]};
      sphere_faces[4*n + 2](i, 1) = Point3D {-curve[i][0], -curve[i][0], curve[i][1]};
      sphere_faces[4*n + 2](i, 2) = Point3D {0.0, -curve[i][0], curve[i][1]};

      sphere_faces[4*n + 3](i, 0) = Point3D {0.0, -curve[i][0], curve[i][1]};
      sphere_faces[4*n + 3](i, 1) = Point3D {curve[i][0], -curve[i][0], curve[i][1]};
      sphere_faces[4*n + 3](i, 2) = Point3D {curve[i][0], 0.0, curve[i][1]};

      const double the_weight = curve.getWeight(i);

      sphere_faces[4*n + 0].setWeight(i, 0, the_weight);
      sphere_faces[4*n + 1].setWeight(i, 0, the_weight);
      sphere_faces[4*n + 2].setWeight(i, 0, the_weight);
      sphere_faces[4*n + 3].setWeight(i, 0, the_weight);

      sphere_faces[4*n + 0].setWeight(i, 1, the_weight / rt2);
      sphere_faces[4*n + 1].setWeight(i, 1, the_weight / rt2);
      sphere_faces[4*n + 2].setWeight(i, 1, the_weight / rt2);
      sphere_faces[4*n + 3].setWeight(i, 1, the_weight / rt2);

      sphere_faces[4*n + 0].setWeight(i, 2, the_weight);
      sphere_faces[4*n + 1].setWeight(i, 2, the_weight);
      sphere_faces[4*n + 2].setWeight(i, 2, the_weight);
      sphere_faces[4*n + 3].setWeight(i, 2, the_weight);
      // clang-format on
    }
  }

  return sphere_faces;
}

axom::Array<primal::NURBSPatch<double, 3>> make_teardrop()
{
  // Generate a teardrop shape using (degenerate) bicubic patches.
  //  The bottom portion is a bicubic sphere,
  //  The top portion is defined by the solid of revolution of a cubic Bezier curve

  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using BCurve = primal::BezierCurve<double, 2>;
  using NPatch = primal::NURBSPatch<double, 3>;

  double rt2 = sqrt(2);

  // Define BezierCurves which, when rotated around the z-axis,
  //  form the top and bottom halves of a sphere.
  BCurve top_curve(3), bottom_curve(2);
  top_curve[0] = Point2D {0.0, 1.0};
  top_curve[1] = Point2D {0.0, 0.0};
  top_curve[2] = Point2D {1.0, 0.0};
  top_curve[3] = Point2D {1.0, -1.0};

  bottom_curve[0] = Point2D {1.0, -1.0};
  bottom_curve[1] = Point2D {1.0, -2.0};
  bottom_curve[2] = Point2D {0.0, -2.0};

  bottom_curve.makeRational();
  bottom_curve.setWeight(1, 1.0 / rt2);

  axom::Array<NPatch> teardrop_faces(8);
  for(int n = 0; n < 4; ++n)
  {
    teardrop_faces[n].setParameters(4, 3, 3, 2);
    teardrop_faces[n].makeRational();

    teardrop_faces[4 + n].setParameters(3, 3, 2, 2);
    teardrop_faces[4 + n].makeRational();
  }

  // Top faces
  for(int i = 0; i <= 3; ++i)
  {
    // clang-format off
    teardrop_faces[0](i, 0) = Point3D {top_curve[i][0], 0.0, top_curve[i][1]};
    teardrop_faces[0](i, 1) = Point3D {top_curve[i][0], top_curve[i][0], top_curve[i][1]};
    teardrop_faces[0](i, 2) = Point3D {0.0, top_curve[i][0], top_curve[i][1]};

    teardrop_faces[1](i, 0) = Point3D {0.0, top_curve[i][0], top_curve[i][1]};
    teardrop_faces[1](i, 1) = Point3D {-top_curve[i][0], top_curve[i][0], top_curve[i][1]};
    teardrop_faces[1](i, 2) = Point3D {-top_curve[i][0], 0.0, top_curve[i][1]};

    teardrop_faces[2](i, 0) = Point3D {-top_curve[i][0], 0.0, top_curve[i][1]};
    teardrop_faces[2](i, 1) = Point3D {-top_curve[i][0], -top_curve[i][0], top_curve[i][1]};
    teardrop_faces[2](i, 2) = Point3D {0.0, -top_curve[i][0], top_curve[i][1]};

    teardrop_faces[3](i, 0) = Point3D {0.0, -top_curve[i][0], top_curve[i][1]};
    teardrop_faces[3](i, 1) = Point3D {top_curve[i][0], -top_curve[i][0], top_curve[i][1]};
    teardrop_faces[3](i, 2) = Point3D {top_curve[i][0], 0.0, top_curve[i][1]};

    teardrop_faces[0].setWeight(i, 1, 1.0 / rt2);
    teardrop_faces[1].setWeight(i, 1, 1.0 / rt2);
    teardrop_faces[2].setWeight(i, 1, 1.0 / rt2);
    teardrop_faces[3].setWeight(i, 1, 1.0 / rt2);
    // clang-format on
  }

  // Bottom faces
  for(int i = 0; i <= 2; ++i)
  {
    // clang-format off
    teardrop_faces[4](i, 0) = Point3D {bottom_curve[i][0], 0.0, bottom_curve[i][1]};
    teardrop_faces[4](i, 1) = Point3D {bottom_curve[i][0], bottom_curve[i][0], bottom_curve[i][1]};
    teardrop_faces[4](i, 2) = Point3D {0.0, bottom_curve[i][0], bottom_curve[i][1]};

    teardrop_faces[5](i, 0) = Point3D {0.0, bottom_curve[i][0], bottom_curve[i][1]};
    teardrop_faces[5](i, 1) = Point3D {-bottom_curve[i][0], bottom_curve[i][0], bottom_curve[i][1]};
    teardrop_faces[5](i, 2) = Point3D {-bottom_curve[i][0], 0.0, bottom_curve[i][1]};

    teardrop_faces[6](i, 0) = Point3D {-bottom_curve[i][0], 0.0, bottom_curve[i][1]};
    teardrop_faces[6](i, 1) = Point3D {-bottom_curve[i][0], -bottom_curve[i][0], bottom_curve[i][1]};
    teardrop_faces[6](i, 2) = Point3D {0.0, -bottom_curve[i][0], bottom_curve[i][1]};

    teardrop_faces[7](i, 0) = Point3D {0.0, -bottom_curve[i][0], bottom_curve[i][1]};
    teardrop_faces[7](i, 1) = Point3D {bottom_curve[i][0], -bottom_curve[i][0], bottom_curve[i][1]};
    teardrop_faces[7](i, 2) = Point3D {bottom_curve[i][0], 0.0, bottom_curve[i][1]};

    const double the_weight = bottom_curve.getWeight(i);

    teardrop_faces[4].setWeight(i, 0, the_weight);
    teardrop_faces[5].setWeight(i, 0, the_weight);
    teardrop_faces[6].setWeight(i, 0, the_weight);
    teardrop_faces[7].setWeight(i, 0, the_weight);

    teardrop_faces[4].setWeight(i, 1, the_weight / rt2);
    teardrop_faces[5].setWeight(i, 1, the_weight / rt2);
    teardrop_faces[6].setWeight(i, 1, the_weight / rt2);
    teardrop_faces[7].setWeight(i, 1, the_weight / rt2);

    teardrop_faces[4].setWeight(i, 2, the_weight);
    teardrop_faces[5].setWeight(i, 2, the_weight);
    teardrop_faces[6].setWeight(i, 2, the_weight);
    teardrop_faces[7].setWeight(i, 2, the_weight);
    // clang-format on
  }

  return teardrop_faces;
}

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_ANALYTIC_TEST_SURFACES_HPP
