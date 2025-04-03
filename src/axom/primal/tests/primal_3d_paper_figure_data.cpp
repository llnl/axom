// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_bezier_curve.cpp
 * \brief This file tests primal's Bezier curve functionality
 */

#include "gtest/gtest.h"

#include "axom/core/FlatMap.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/operators/winding_number_2d.hpp"
#include "axom/primal/operators/winding_number_3d.hpp"
#include "axom/primal/operators/detail/winding_number_2d_impl.hpp"
#include "axom/primal/operators/detail/winding_number_3d_impl.hpp"

#include "axom/primal/operators/printers.hpp"

#include <chrono>
#include <string>
#include <iomanip>

namespace primal = axom::primal;
using Point2D = primal::Point<double, 2>;
using Point3D = primal::Point<double, 3>;
using Vector3D = primal::Vector<double, 3>;
using BCurve2D = primal::BezierCurve<double, 2>;
using BPatch = primal::BezierPatch<double, 3>;
using BBox = primal::BoundingBox<double, 3>;
using Tri = primal::Triangle<double, 3>;

TEST(primal_3d_paper_figure_data, plotting_demo)
{
  return;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\plotting_demo";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  double rt2 = sqrt(2), rt3 = sqrt(3), rt6 = sqrt(6);

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

  BPatch sphere_faces[6];
  for(int n = 0; n < 6; ++n)
  {
    sphere_faces[n].setOrder(4, 4);
    sphere_faces[n].makeRational();
  }

  sphere_faces[0].setOrder(4, 4);
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

  double the_wn = winding_number_casting(Point3D {0.1, 0.1, 0.0},
                                         sphere_faces[1],
                                         edge_tol,
                                         quad_tol,
                                         EPS);
  std::cout << the_wn << std::endl;

  return;

  BBox bbox;
  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      for(int n = 0; n < 6; ++n)
      {
        bbox.addPoint(sphere_faces[n](i, j));
      }
    }
  }

  axom::Array<BPatch> patches;
  patches.push_back(sphere_faces[0]);

  auto wn_accurate = [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
    double wn = 0.0;
    for(const auto& patch : patches)
      wn += winding_number_casting(query, patch, edge_tol, quad_tol, EPS);
    return wn;
  };

  // clang-format off
  // exportScalarFieldToVTK(data_dir + "/sphere_field.vtk", wn_accurate, bbox, 50, 50, 50);
  // exportScalarFieldToVTK(data_dir + "/sphere_field_stokes.vtk", wn_stokes, bbox, 50, 50, 50);
  // clang-format on

  // exportSurfaceToSTL(data_dir + "/sphere_face.stl", patches);
}

TEST(primal_3d_paper_figure_data, patch_subdivision_test)
{
  return;

  double rt2 = sqrt(2), rt3 = sqrt(3), rt6 = sqrt(6);

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

  BPatch sphere_face;
  sphere_face.setOrder(4, 4);
  sphere_face.makeRational();

  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      const int idx = 5 * i + j;
      sphere_face.setWeight(i, j, weight_data[idx]);

      // Set up each face by rotating one of the patch faces
      sphere_face(i, j)[0] = -node_data[idx][0];
      sphere_face(i, j)[1] = -node_data[idx][1];
      sphere_face(i, j)[2] = -node_data[idx][2];

      sphere_face(i, j).array() /= weight_data[idx];
    }
  }

  // Define a recursive lambda function to check each subpatch for isPlanar
  std::function<bool(const BPatch&, int n)> check_planar =
    [&](const BPatch& patch, int n) {
      if(patch.isPlanar(1e-4))
      {
        std::cout << "Subpatch " << n << " is planar" << std::endl;
        return true;
      }
      else
      {
        BPatch p1, p2, p3, p4;
        patch.split(0.5, 0.5, p1, p2, p3, p4);
        return check_planar(p1, n + 1) && check_planar(p2, n + 1) &&
          check_planar(p3, n + 1) && check_planar(p4, n + 1);
      }
    };

  check_planar(sphere_face, 0);
}

TEST(primal_3d_paper_figure_data, super_cusp_example)
{
  return;

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\very_degenerate_"
    "case";

  BPatch tent;
  tent.setOrder(4, 4);

  double D = 16.0 / 7.0;

  // clang-format off
  tent(0, 0)[0] = 1.0;  tent(0, 1)[0] = 2.0;    tent(0, 2)[0] = 3.0;  tent(0, 3)[0] = 4.0;    tent(0, 4)[0] = 5.0; 
  tent(1, 0)[0] = 1.0;  tent(1, 1)[0] = 2.0+D;  tent(1, 2)[0] = 3.0;  tent(1, 3)[0] = 4.0-D;  tent(1, 4)[0] = 5.0; 
  tent(2, 0)[0] = 1.0;  tent(2, 1)[0] = 2.0+D;  tent(2, 2)[0] = 3.0;  tent(2, 3)[0] = 4.0-D;  tent(2, 4)[0] = 5.0; 
  tent(3, 0)[0] = 1.0;  tent(3, 1)[0] = 2.0+D;  tent(3, 2)[0] = 3.0;  tent(3, 3)[0] = 4.0-D;  tent(3, 4)[0] = 5.0; 
  tent(4, 0)[0] = 1.0;  tent(4, 1)[0] = 2.0;    tent(4, 2)[0] = 3.0;  tent(4, 3)[0] = 4.0;    tent(4, 4)[0] = 5.0; 

  tent(0, 0)[1] = 1.0;  tent(0, 1)[1] = 1.0;    tent(0, 2)[1] = 1.0;    tent(0, 3)[1] = 1.0;    tent(0, 4)[1] = 1.0; 
  tent(1, 0)[1] = 2.0;  tent(1, 1)[1] = 2.0+D;  tent(1, 2)[1] = 2.0+D;  tent(1, 3)[1] = 2.0+D;  tent(1, 4)[1] = 2.0; 
  tent(2, 0)[1] = 3.0;  tent(2, 1)[1] = 3.0;    tent(2, 2)[1] = 3.0;    tent(2, 3)[1] = 3.0;    tent(2, 4)[1] = 3.0; 
  tent(3, 0)[1] = 4.0;  tent(3, 1)[1] = 4.0-D;  tent(3, 2)[1] = 4.0-D;  tent(3, 3)[1] = 4.0-D;  tent(3, 4)[1] = 4.0; 
  tent(4, 0)[1] = 5.0;  tent(4, 1)[1] = 5.0;    tent(4, 2)[1] = 5.0;    tent(4, 3)[1] = 5.0;    tent(4, 4)[1] = 5.0; 

  tent(0, 0)[2] = 1.0;  tent(0, 1)[2] = 1.0;  tent(0, 2)[2] =  1.0;  tent(0, 3)[2] = 1.0;  tent(0, 4)[2] = 1.0; 
  tent(1, 0)[2] = 1.0;  tent(1, 1)[2] = 1.0;  tent(1, 2)[2] =  1.0;  tent(1, 3)[2] = 1.0;  tent(1, 4)[2] = 1.0; 
  tent(2, 0)[2] = 1.0;  tent(2, 1)[2] = 1.0;  tent(2, 2)[2] = -2.0;  tent(2, 3)[2] = 1.0;  tent(2, 4)[2] = 1.0; 
  tent(3, 0)[2] = 1.0;  tent(3, 1)[2] = 1.0;  tent(3, 2)[2] =  1.0;  tent(3, 3)[2] = 1.0;  tent(3, 4)[2] = 1.0; 
  tent(4, 0)[2] = 1.0;  tent(4, 1)[2] = 1.0;  tent(4, 2)[2] =  1.0;  tent(4, 3)[2] = 1.0;  tent(4, 4)[2] = 1.0;

  tent.reverseOrientation_u();
  // clang-format on

  primal::Polygon<double, 3> closures[6];
  closures[0].addVertex(Point3D {1.0, 1.0, 1.0});
  closures[0].addVertex(Point3D {1.0, 5.0, 1.0});
  closures[0].addVertex(Point3D {1.0, 5.0, -1.0});
  closures[0].addVertex(Point3D {1.0, 1.0, -1.0});

  closures[1].addVertex(Point3D {5.0, 1.0, 1.0});
  closures[1].addVertex(Point3D {1.0, 1.0, 1.0});
  closures[1].addVertex(Point3D {1.0, 1.0, -1.0});
  closures[1].addVertex(Point3D {5.0, 1.0, -1.0});

  closures[2].addVertex(Point3D {5.0, 5.0, 1.0});
  closures[2].addVertex(Point3D {5.0, 1.0, 1.0});
  closures[2].addVertex(Point3D {5.0, 1.0, -1.0});
  closures[2].addVertex(Point3D {5.0, 5.0, -1.0});

  closures[3].addVertex(Point3D {1.0, 5.0, 1.0});
  closures[3].addVertex(Point3D {5.0, 5.0, 1.0});
  closures[3].addVertex(Point3D {5.0, 5.0, -1.0});
  closures[3].addVertex(Point3D {1.0, 5.0, -1.0});

  closures[4].addVertex(Point3D {1.0, 1.0, -1.0});
  closures[4].addVertex(Point3D {1.0, 5.0, -1.0});
  closures[4].addVertex(Point3D {5.0, 5.0, -1.0});
  closures[4].addVertex(Point3D {5.0, 1.0, -1.0});

  closures[5].addVertex(Point3D {5.0, 1.0, 1.0});
  closures[5].addVertex(Point3D {5.0, 5.0, 1.0});
  closures[5].addVertex(Point3D {1.0, 5.0, 1.0});
  closures[5].addVertex(Point3D {1.0, 1.0, 1.0});

  Point3D the_query = tent.evaluate(0.5, 0.5);

  double closure_wn = 0.0;
  for(int i = 0; i < 5; ++i)
  {
    closure_wn += winding_number(the_query, closures[i]);
  }

  std::cout << "Closure winding number: " << closure_wn << std::endl;

  // exportSurfaceToSTL(data_dir + "/super_neato_tent.stl", tent);

  std::pair<double, double> wn;
  double wn_ref, wn_direct;

  // wn =
  // winding_number_casting_split(the_query, tent, edge_tol, quad_tol, EPS);
  // wn_ref = winding_number(the_query, tent, edge_tol, quad_tol, EPS);
  // wn_direct =
  // winding_number_direct(the_query, tent, edge_tol, quad_tol, EPS);

  // std::cout << wn.first << " + " << wn.second << " -> " << wn.first + wn.second
  // << " (" << wn.first + wn.second + closure_wn << ")" << std::endl;
  // std::cout << wn_ref << " (" << wn_ref + closure_wn << ")" << std::endl;
  // std::cout << wn_direct << " (" << wn_direct + closure_wn << ")" << std::endl
  // << std::endl
  // << std::endl;

  BPatch tents[4];
  tent.split(0.5, 0.5, tents[0], tents[1], tents[2], tents[3]);

  std::pair<double, double> wns[4];
  double wn_refs[4], wn_directs[4];

  // std::cout << tents[0] << std::endl;

  for(int i = 0; i < 4; ++i)
  {
    wns[i] =
      winding_number_casting_split(the_query, tents[i], edge_tol, quad_tol, EPS);
    wn_refs[i] = winding_number(the_query, tents[i], edge_tol, quad_tol, EPS);
    wn_directs[i] =
      winding_number_direct(the_query, tents[i], edge_tol, quad_tol, EPS);

    std::cout << wns[i].first << " + " << wns[i].second << " -> "
              << wns[i].first + wns[i].second << " ("
              << wns[i].first + wns[i].second << ")" << std::endl;
    std::cout << wn_refs[i] << " (" << wn_refs[i] << ")" << std::endl;
    std::cout << wn_directs[i] << " (" << wn_directs[i] << ")" << std::endl
              << std::endl
              << std::endl;
  }

  std::cout << "Total winding number: "
            << wns[0].first + wns[0].second + wns[1].first + wns[1].second +
      wns[2].first + wns[2].second + wns[3].first + wns[3].second + closure_wn
            << std::endl;
}

TEST(primal_3d_paper_figure_data, cusp_example)
{
  return;

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\very_degenerate_"
    "case";

  BPatch tent;
  tent.setOrder(3, 3);

  // clang-format off
  tent(0, 0)[0] = -1.0;  tent(0, 1)[0] = 1.0;  tent(0, 2)[0] = -1.0;  tent(0, 3)[0] = 1.0; 
  tent(1, 0)[0] = -1.0;  tent(1, 1)[0] = 1.0;  tent(1, 2)[0] = -1.0;  tent(1, 3)[0] = 1.0; 
  tent(2, 0)[0] = -1.0;  tent(2, 1)[0] = 1.0;  tent(2, 2)[0] = -1.0;  tent(2, 3)[0] = 1.0; 
  tent(3, 0)[0] = -1.0;  tent(3, 1)[0] = 1.0;  tent(3, 2)[0] = -1.0;  tent(3, 3)[0] = 1.0; 

  tent(0, 0)[1] = -1.0;  tent(0, 1)[1] = 1.0;  tent(0, 2)[1] = -1.0;  tent(0, 3)[1] = 1.0; 
  tent(1, 0)[1] =  0.0;  tent(1, 1)[1] = 0.0;  tent(1, 2)[1] =  0.0;  tent(1, 3)[1] = 0.0; 
  tent(2, 0)[1] =  0.0;  tent(2, 1)[1] = 0.0;  tent(2, 2)[1] =  0.0;  tent(2, 3)[1] = 0.0; 
  tent(3, 0)[1] =  1.0;  tent(3, 1)[1] = 1.0;  tent(3, 2)[1] =  1.0;  tent(3, 3)[1] = 1.0; 

  tent(0, 0)[2] =  0.0;  tent(0, 1)[2] = 1.0;  tent(0, 2)[2] =  1.0;  tent(0, 3)[2] = 0.0; 
  tent(1, 0)[2] =  0.0;  tent(1, 1)[2] = 1.0;  tent(1, 2)[2] =  1.0;  tent(1, 3)[2] = 0.0; 
  tent(2, 0)[2] =  0.0;  tent(2, 1)[2] = 1.0;  tent(2, 2)[2] =  1.0;  tent(2, 3)[2] = 0.0; 
  tent(3, 0)[2] =  0.0;  tent(3, 1)[2] = 1.0;  tent(3, 2)[2] =  1.0;  tent(3, 3)[2] = 0.0;
  // clang-format on

  // exportSurfaceToSTL(data_dir + "/neato_tent.stl", tent);

  Point3D the_query = tent.evaluate(0.5, 0.5);

  auto wn =
    winding_number_casting_split(the_query, tent, edge_tol, quad_tol, EPS);
  auto wn_ref = winding_number(the_query, tent, edge_tol, quad_tol, EPS);
  auto wn_direct =
    winding_number_direct(the_query, tent, edge_tol, quad_tol, EPS);

  std::cout << wn.first << " + " << wn.second << " -> " << wn.first + wn.second
            << std::endl;
  std::cout << wn_ref << std::endl;
  std::cout << wn_direct << std::endl;
}

TEST(primal_3d_paper_figure_data, very_degenerate_example)
{
  return;

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\very_degenerate_"
    "case";

  double D = 8.0 / 9.0;

  BPatch flower;
  flower.setOrder(3, 3);

  flower(0, 0) = Point3D {-1.0, -1.0, 1.0};
  flower(1, 0) = Point3D {1.0, -1.0, 0.0};
  flower(2, 0) = Point3D {1.0, 1.0, 0.0};
  flower(3, 0) = Point3D {-1.0, 1.0, 1.0};

  flower(0, 1) = Point3D {-D, D, 0.0};
  flower(1, 1) = Point3D {0.0, 0.0, 0.0};
  flower(2, 1) = Point3D {0.0, 0.0, 0.0};
  flower(3, 1) = Point3D {-D, -D, 0.0};

  flower(0, 2) = Point3D {D, D, 0.0};
  flower(1, 2) = Point3D {0.0, 0.0, 0.0};
  flower(2, 2) = Point3D {0.0, 0.0, 0.0};
  flower(3, 2) = Point3D {D, -D, 0.0};

  flower(0, 3) = Point3D {1.0, -1.0, 1.0};
  flower(1, 3) = Point3D {-1.0, -1.0, 0.0};
  flower(2, 3) = Point3D {-1.0, 1.0, 0.0};
  flower(3, 3) = Point3D {1.0, 1.0, 1.0};

  // exportSurfaceToSTL(data_dir + "/wacko_flower.stl", flower);

  std::cout << flower.evaluate(0.25, 0.85) << std::endl;
  std::cout << flower.du(0.5, 0.5) << std::endl;
  std::cout << flower.dv(0.5, 0.5) << std::endl;
  std::cout << flower.normal(0.5, 0.5) << std::endl;

  Point3D the_query = {0.0, 0.0, 0.0625};

  auto wn =
    winding_number_casting_split(the_query, flower, edge_tol, quad_tol, EPS);
  auto wn_ref = winding_number(the_query, flower, edge_tol, quad_tol, EPS);
  auto wn_direct =
    winding_number_direct(the_query, flower, edge_tol, quad_tol, EPS);

  std::cout << wn.first << " + " << wn.second << " -> " << wn.first + wn.second
            << std::endl;
  std::cout << wn_ref << std::endl;
  std::cout << wn_direct << std::endl;
}

TEST(primal_3d_paper_figure_data, edge_case_explorer)
{
  return;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\full_patch_test";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // Make a patch that is a unit square
  BPatch square;
  square.setOrder(3, 3);

  square(0, 0) = Point3D {0.0, 0.0, 0.0};
  square(1, 0) = Point3D {1.0, -1.0, 0.0};
  square(2, 0) = Point3D {1.0, -1.0, 1.0};
  square(3, 0) = Point3D {0.0, 0.0, 2.0};

  square(0, 1) = Point3D {0.0, 0.0, 0.0};
  square(1, 1) = Point3D {1.0, 0.0, 0.0};
  square(2, 1) = Point3D {1.0, 0.0, 1.0};
  square(3, 1) = Point3D {0.0, 0.0, 2.0};

  square(0, 2) = Point3D {0.0, 0.5, 0.0};
  square(1, 2) = Point3D {1.0, 0.5, 0.0};
  square(2, 2) = Point3D {1.0, 0.5, 1.0};
  square(3, 2) = Point3D {0.0, 0.5, 2.0};

  square(0, 3) = Point3D {0.0, 1.0, 0.0};
  square(1, 3) = Point3D {1.0, 1.0, 0.0};
  square(2, 3) = Point3D {1.0, 1.0, 1.0};
  square(3, 3) = Point3D {0.0, 1.0, 2.0};

  // exportSurfaceToSTL(data_dir + "/edge_case.stl", square);

  Point3D the_query = {0.0, 0.0, 0.0};

  auto wn =
    winding_number_casting_split(the_query, square, edge_tol, quad_tol, EPS);
  auto wn_ref = winding_number(the_query, square, edge_tol, quad_tol, EPS);
  auto wn_direct =
    winding_number_direct(the_query, square, edge_tol, quad_tol, EPS);

  std::cout << wn.first << " + " << wn.second << " -> " << wn.first + wn.second
            << std::endl;
  std::cout << wn_ref << std::endl;
  std::cout << wn_direct << std::endl;
}

TEST(primal_3d_paper_figure_data, full_patch_test)
{
  //   return;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\full_patch_test";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  double rt2 = sqrt(2), rt3 = sqrt(3), rt6 = sqrt(6);

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

  BPatch sphere_faces[6];
  for(int n = 0; n < 6; ++n)
  {
    sphere_faces[n].setOrder(4, 4);
    sphere_faces[n].makeRational();
  }

  sphere_faces[0].setOrder(4, 4);
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

  // Lambda to generate an entirely random unit vector
  auto random_unit = []() -> axom::primal::Vector<double, 3> {
    double theta = axom::utilities::random_real(0.0, 2 * M_PI);
    double u = axom::utilities::random_real(-1.0, 1.0);
    return axom::primal::Vector<double, 3> {sin(theta) * sqrt(1 - u * u),
                         cos(theta) * sqrt(1 - u * u),
                         u};
  };

  // Evaluate the winding number for random points in the box [-2.028, 2.028]^3
  for(int N = 0; N < 10; ++N)
  {
    auto the_direction = random_unit() * ( 1 - 1e-4 );

    Point3D in_query = {the_direction[0], the_direction[1], the_direction[2]};

    double wn_in = 0.0;

    int indices[] = {0, 1, 2, 3, 4, 5};
    for(int n : indices)
    {
      wn_in +=
        winding_number_casting(in_query, sphere_faces[n], edge_tol, quad_tol, EPS);
    }

    std::cout << std::setprecision(15);

    // Check the inside point
    if( !axom::utilities::isNearlyEqual( wn_in, 1.0, 6 * quad_tol ) )
    {
      std::cout << '\t' << "Inside point: " << in_query << " -> " << wn_in << std::endl;
    }
  }
}

TEST(primal_3d_paper_figure_data, rotating_patch)
{
  return;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\rotating_patch";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  double rt2 = sqrt(2), rt3 = sqrt(3), rt6 = sqrt(6);

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

  BPatch sphere_face;
  sphere_face.setOrder(4, 4);
  sphere_face.makeRational();

  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      const int idx = 5 * i + j;
      sphere_face.setWeight(i, j, weight_data[idx]);

      // Set up each face by rotating one of the patch faces
      sphere_face(i, j)[0] = -node_data[idx][0];
      sphere_face(i, j)[1] = -node_data[idx][1];
      sphere_face(i, j)[2] = -node_data[idx][2];

      sphere_face(i, j).array() /= weight_data[idx];
    }
  }

  auto the_query = sphere_face.evaluate(0.0, 0.0);

  auto wn_true = winding_number(the_query, sphere_face, edge_tol, quad_tol, EPS);
  auto wn_split =
    winding_number_casting_split(the_query, sphere_face, edge_tol, quad_tol, EPS);

  std::cout << wn_true << std::endl;
  std::cout << wn_split.first << " " << wn_split.second << std::endl;

  if(axom::utilities::isNearlyEqual(wn_true, wn_split.first + wn_split.second, 1e-5))
  {
    std::cout << "Winding numbers match" << std::endl;
  }
  else
  {
    std::cout << "Winding numbers do not match" << std::endl;
  }

  //   return;

  BBox bbox;
  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      bbox.addPoint(sphere_face(i, j));
    }
  }

  // Make the bounding box a cube with the same centroid
  // int max_dim = bbox.getLongestDimension();
  // double max_len = bbox.getMax()[max_dim] - bbox.getMin()[max_dim];

  // primal::Point<double, 3> centroid = bbox.getCentroid();

  // primal::Point<double, 3> new_min {centroid[0] - max_len / 2.0,
  //                                   centroid[1] - max_len / 2.0,
  //                                   centroid[2] - max_len / 2.0};
  // primal::Point<double, 3> new_max {centroid[0] + max_len / 2.0,
  //                                   centroid[1] + max_len / 2.0,
  //                                   centroid[2] + max_len / 2.0};

  // bbox.addPoint(new_min);
  // bbox.addPoint(new_max);
  // bbox.scale(1.1);

  axom::Array<BPatch> patches;
  patches.push_back(sphere_face);

  auto wn_ground_truth =
    [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& patch : patches)
        wn += winding_number(query, patch, edge_tol, quad_tol, EPS);
      return wn;
    };

  auto wn_casting = [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
    std::pair<double, double> wn_split = {0.0, 0.0};
    for(const auto& patch : patches)
    {
      auto val =
        winding_number_casting_split(query, patch, edge_tol, quad_tol, EPS);
      wn_split.first += val.first;
      wn_split.second += val.second;
    }
    return wn_split;
  };

  auto wn_comparison = [&patches, &edge_tol, &quad_tol, &EPS](
                         const Point3D& query) {
    std::pair<double, double> wn_split = {0.0, 0.0};
    for(const auto& patch : patches)
    {
      auto val =
        winding_number_casting_split(query, patch, edge_tol, quad_tol, EPS);
      wn_split.first += val.first;
      wn_split.second += val.second;
    }

    double wn = 0.0;
    for(const auto& patch : patches)
      wn += winding_number(query, patch, edge_tol, quad_tol, EPS);

    if(!axom::utilities::isNearlyEqual(wn, wn_split.first + wn_split.second, 1e-5))
    {
      std::cout << "Winding numbers do not match!" << std::endl;
      std::cout << '\t' << query << ": " << wn << " vs "
                << wn_split.first + wn_split.second << std::endl;
    }

    return wn_split;
  };

  // clang-format off
  // exportScalarFieldToVTK(data_dir + "/sphere_field.vtk", wn_ground_truth, bbox, 100, 100, 100);
  // exportSplitScalarFieldToVTK(data_dir + "/sphere_field_casting_complete.vtk", wn_comparison, bbox, 100, 100, 100);
  // clang-format on

  // exportSurfaceToSTL(data_dir + "/sphere_face.stl", patches);
}

TEST(primal_3d_paper_figure_data, trimmed_sphere)
{
  return;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\trimmed_sphere";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  double rt2 = sqrt(2), rt3 = sqrt(3), rt6 = sqrt(6);

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

  BPatch sphere_faces[6];
  for(int n = 0; n < 6; ++n)
  {
    sphere_faces[n].setOrder(4, 4);
    sphere_faces[n].makeRational();
  }

  sphere_faces[0].setOrder(4, 4);
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

  BBox bbox;
  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      for(int n = 0; n < 6; ++n)
      {
        bbox.addPoint(sphere_faces[n](i, j));
      }
    }
  }

  // Create a patch that has the inner 90% of sphere_faces[5]
  BPatch p1, p2, p3, p4, p5;
  sphere_faces[5].split(0.9, 0.9, p1, p2, p3, p4);
  p1.split(0.09, 0.09, p2, p3, p4, p5);
  sphere_faces[5] = p5;

  axom::Array<BPatch> patches;
  for(int n = 0; n < 6; ++n) patches.push_back(sphere_faces[n]);

  auto wn_accurate = [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
    double wn = 0.0;
    for(const auto& patch : patches)
      wn += winding_number(query, patch, edge_tol, quad_tol, EPS);
    return wn;
  };

  // exportSurfaceToSTL(data_dir + "/sphere_face.stl", patches, 50, 50);

  // clang-format off
  // exportScalarFieldToVTK(data_dir + "/sphere_field.vtk", wn_accurate, bbox, 300, 300, 300);
  // clang-format on
}

TEST(primal_3d_paper_figure_data, vase_shape)
{
  return;
  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\vase_shape";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // lambda function that rotates a given BezierCurve around the z-axis
  auto rotate_curve = [](const BCurve2D& curve) -> axom::Array<BPatch> {
    const int ord = curve.getOrder();
    axom::Array<BPatch> rs(4);
    for(int i = 0; i < 4; ++i)
    {
      rs[i].setOrder(ord, 2);
      rs[i].makeRational();
    }

    for(int i = 0; i <= ord; ++i)
    {
      rs[0](i, 0) = Point3D {curve[i][0], 0.0, curve[i][1]};
      rs[0](i, 1) = Point3D {curve[i][0], curve[i][0], curve[i][1]};
      rs[0](i, 2) = Point3D {0.0, curve[i][0], curve[i][1]};

      rs[1](i, 0) = Point3D {0.0, curve[i][0], curve[i][1]};
      rs[1](i, 1) = Point3D {-curve[i][0], curve[i][0], curve[i][1]};
      rs[1](i, 2) = Point3D {-curve[i][0], 0.0, curve[i][1]};

      rs[2](i, 0) = Point3D {-curve[i][0], 0.0, curve[i][1]};
      rs[2](i, 1) = Point3D {-curve[i][0], -curve[i][0], curve[i][1]};
      rs[2](i, 2) = Point3D {0.0, -curve[i][0], curve[i][1]};

      rs[3](i, 0) = Point3D {0.0, -curve[i][0], curve[i][1]};
      rs[3](i, 1) = Point3D {curve[i][0], -curve[i][0], curve[i][1]};
      rs[3](i, 2) = Point3D {curve[i][0], 0.0, curve[i][1]};

      for(int j = 0; j < 4; ++j)
      {
        rs[j].setWeight(i, 1, 1.0 / std::sqrt(2));
      }
    }

    return rs;
  };

  BCurve2D vase1(5);
  vase1[0] = Point2D {0.0, -1.0};
  vase1[1] = Point2D {1.0, -1.0};
  vase1[2] = Point2D {0.5, 0.5};
  vase1[3] = Point2D {0.0, 0.8};
  vase1[4] = Point2D {0.8, 0.5};
  vase1[5] = Point2D {0.8, 1.0};

  BCurve2D vase2(5);
  vase2[0] = Point2D {0.8, 1.0};
  vase2[1] = Point2D {0.8, 1.5};
  vase2[2] = Point2D {0.0, 0.8};
  vase2[3] = Point2D {0.3, 1.5};
  vase2[4] = Point2D {0.5, -0.8};
  vase2[5] = Point2D {0.0, -0.8};

  auto vase1_patches = rotate_curve(vase1);
  auto vase2_patches = rotate_curve(vase2);

  BCurve2D trim(2);
  trim[0] = Point2D {0.55, 0.2};
  trim[1] = Point2D {1.0, 0.4};
  trim[2] = Point2D {0.4, 0.4};

  auto trim_patches = rotate_curve(trim);

  BBox bbox;
  bbox.addPoint(Point3D {1.25, 1.25, -1.25});
  bbox.addPoint(Point3D {-1.25, -1.25, 1.75});

  axom::Array<BPatch> patches;
  for(int n = 0; n < 4; ++n)
  {
    patches.push_back(vase1_patches[n]);
    patches.push_back(vase2_patches[n]);
  }

  auto vase_wn_accurate =
    [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& patch : patches)
        wn += winding_number(query, patch, edge_tol, quad_tol, EPS);
      return wn;
    };

  auto trim_wn_accurate =
    [&trim_patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& trim_patches : trim_patches)
        wn += winding_number(query, trim_patches, edge_tol, quad_tol, EPS);
      return wn;
    };

  auto bad_query = Point3D {-0.38793103448275867429,
                            -0.56034482758620685061,
                            -0.93965517241379314939};

  for(int i = 0; i < patches.size(); ++i)
  {
    axom::Array<BPatch> single_patch;
    single_patch.push_back(patches[i]);
    // exportSurfaceToSTL(data_dir + "/vase_shape_" + std::to_string(i) + ".stl",
                      //  single_patch,
                      //  50,
                      //  50);
  }

  for(int i = 0; i < trim_patches.size(); ++i)
  {
    axom::Array<BPatch> single_patch;
    single_patch.push_back(trim_patches[i]);
    // exportSurfaceToSTL(data_dir + "/trim_shape_" + std::to_string(i) + ".stl",
                      //  single_patch,
                      //  50,
                      //  50);
  }

  // exportSurfaceToSTL(data_dir + "/vase_shape.stl", patches, 50, 50);
  // exportSurfaceToSTL(data_dir + "/trim_shape.stl", trim_patches, 50, 50);

  // clang-format off
  //exportScalarFieldToVTK(data_dir + "/vase_field.vtk", vase_wn_accurate, bbox, 300, 300, 300);
  //exportScalarFieldToVTK(data_dir + "/trim_field.vtk", trim_wn_accurate, bbox, 300, 300, 300);
  // clang-format on
}

TEST(primal_3d_paper_figure_data, two_teapots)
{
  return;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\teapot_shape\\";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  std::ifstream file(data_dir + "teapot.stl");

  BBox bbox;
  bbox.addPoint(Point3D {0, 0, -55});
  bbox.addPoint(Point3D {160, 90, 55});

  if(!file.is_open())
  {
    std::cerr << "Unable to open file" << std::endl;
    return;
  }

  std::string line;
  Point3D p1, p2;
  Point3D p3;
  int vertexIndex = 0;

  axom::Array<Tri> half_patches;
  axom::Array<Tri> mirror_patches;

  while(std::getline(file, line))
  {
    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword;

    if(keyword == "vertex")
    {
      Point3D the_point;
      iss >> the_point[0] >> the_point[1] >> the_point[2];

      switch(vertexIndex)
      {
      case 0:
        p1 = the_point;
        break;
      case 1:
        p2 = the_point;
        break;
      case 2:
        p3 = the_point;
        break;
      }

      vertexIndex++;
    }
    else if(keyword == "endloop")
    {
      half_patches.push_back(Tri {p1, p2, p3});
      p1[2] = -p1[2];
      p2[2] = -p2[2];
      p3[2] = -p3[2];
      mirror_patches.push_back(Tri {p1, p2, p3});
      vertexIndex = 0;
    }
  }

  auto half_teapot_wn =
    [&half_patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& patch : half_patches) wn += winding_number(query, patch);
      return wn;
    };

  auto mirror_teapot_wn =
    [&mirror_patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& patch : mirror_patches)
        wn += winding_number(query, patch);
      return wn;
    };

  // exportScalarFieldToVTK(data_dir + "/teapot_half_field_high.vtk",
                        //  half_teapot_wn,
                        //  bbox,
                        //  400,
                        //  200,
                        //  200);
  // exportScalarFieldToVTK(data_dir + "/teapot_mirror_field_high.vtk",
                        //  mirror_teapot_wn,
                        //  bbox,
                        //  400,
                        //  200,
                        //  200);
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
