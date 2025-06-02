// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file printers.hpp
 *
 * \brief Consists of jacob's neat little debugging print statements
 * 
 * \note If you see this in a repo, it reflects a failure on jacob's part.
 *       Reprimand him for his error.
 */

#ifndef PRIMAL_PRINTERS_HPP_
#define PRIMAL_PRINTERS_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/winding_number.hpp"
#include "axom/core/numerics/matvecops.hpp"

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <regex>
#include <vector>

namespace axom
{
namespace primal
{

template <typename T>
void python_print(const Point<T, 3>& point, const char* plot_lines = nullptr)
{
  printf("plot_query( fig, ax, (%.17f, %.17f, %.17f), ", point[0], point[1], point[2]);
  if(plot_lines) printf("plot_lines=\"%s\", ", plot_lines);
  printf("size=111)\n\n");
}

template <typename T>
void python_print(const Vector<T, 3>& vec, const Point<T, 3>& point, const char* color = "blue")
{
  printf("plot_vector( fig, ax, (%.17f, %.17f, %.17f), ", vec[0], vec[1], vec[2]);
  printf("origin=(%.17f, %.17f, %.17f),", point[0], point[1], point[2]);
  printf("color=\"%s\")\n", color);
}

template <typename T>
void desmos_print(const Point<T, 3>& point)
{
  printf("Q = [%.17f, %.17f, %.17f]\n", point[0], point[1], point[2]);
}

template <typename T>
void desmos_print(const Point<T, 2>& point)
{
  printf("Q = (%.17f, %.17f)\n", point[0], point[1]);
}

template <typename T>
void python_print(const BezierCurve<T, 3>& curve,
                  bool plot_tangent = false,
                  bool plot_endpoints = false,
                  char color = 'k')
{
  const int ord = curve.getOrder();

  printf("my_BezierCurve( [");
  printf("(%.17f,%.17f,%.17f)", curve[0][0], curve[0][1], curve[0][2]);
  for(int i = 1; i <= ord; ++i)
    printf(",(%.17f,%.17f,%.17f)", curve[i][0], curve[i][1], curve[i][2]);
  if(curve.isRational())
  {
    printf(" ], [%.17f", curve.getWeight(0));
    if(curve.isRational())
      for(int i = 1; i <= ord; ++i) printf(",%.17f", curve.getWeight(i));
  }
  printf("] ).plot( fig, ax, None, '%c'", color);
  if(plot_tangent) printf(", plot_tangent=True");
  if(plot_endpoints) printf(", plot_endpoints=True");
  printf(" )\n\n");
}

template <typename T>
void desmos_print(const BezierCurve<T, 2>& curve, int num = 0)
{
  printf("C(P_{%d}, W_{%d})\n", num, num);
  printf("P_{%d} = [(%.17f,%.17f)", num, curve[0][0], curve[0][1]);
  for(int p = 1; p <= curve.getOrder(); ++p) printf(",(%.17f,%.17f)", curve[p][0], curve[p][1]);
  printf("]\n");

  if(curve.isRational())
  {
    printf("W_{%d} = [%.17f", num, curve.getWeight(0));
    for(int p = 1; p <= curve.getOrder(); ++p) printf(",%.17f", curve.getWeight(p));
    printf("]\n");
  }
  else
  {
    printf("W_{%d} = [1", num);
    for(int p = 1; p <= curve.getOrder(); ++p) printf(",1");
    printf("]\n");
  }
}

template <typename T>
void desmos_print(const BezierCurve<T, 3>& curve)
{
  const int ord = curve.getOrder();
  for(int i = 0; i <= ord; ++i)
  {
    printf("P_{%d} = [%.17f, %.17f, %.17f]\n", i + 1, curve[i][0], curve[i][1], curve[i][2]);
  }

  if(curve.isRational())
  {
    printf("W = [%.17f", curve.getWeight(0));
    for(int i = 1; i <= ord; ++i)
    {
      printf(", %.17f", curve.getWeight(i));
    }
    printf("]\n");
  }
}

template <typename T>
void python_print(const BezierPatch<T, 3>& patch,
                  const char* cmap = "cm.Reds",
                  char cpcolor = '\0',
                  bool plot_normal = false)
{
  const int ord_u = patch.getOrder_u();
  const int ord_v = patch.getOrder_v();

  printf("my_BezierPatch(%d, %d, [", ord_u, ord_v);
  printf("(%.17f,%.17f,%.17f)", patch(0, 0)[0], patch(0, 0)[1], patch(0, 0)[2]);
  for(int i = 0; i <= ord_u; ++i)
    for(int j = (i == 0 ? 1 : 0); j <= ord_v; ++j)
      printf(",(%.17f,%.17f,%.17f)", patch(i, j)[0], patch(i, j)[1], patch(i, j)[2]);
  if(patch.isRational())
  {
    printf(" ], [%.17f", patch.getWeight(0, 0));
    if(patch.isRational())
      for(int i = 0; i <= ord_u; ++i)
        for(int j = (i == 0 ? 1 : 0); j <= ord_v; ++j) printf(",%.17f", patch.getWeight(i, j));
  }
  printf("] ).plot( fig, ax, %s", cmap);
  if(plot_normal) printf(", plot_normal=True");
  if(cpcolor != '\0') printf(", cpcolor='%c'", cpcolor);
  printf(")\n\n");
}

template <typename T>
void desmos_print(const BezierPatch<T, 3>& patch)
{
  // clang-format off
  const int ord_u = patch.getOrder_u();
  const int ord_v = patch.getOrder_v();

  printf("M = %d\n", ord_u);
  printf("N = %d\n", ord_v);

  // Print u basis functions
  for(int u = 0; u <= ord_u; ++u)
    printf("B_{%dM}(u) = nCr\\left(M, %d\\right)(1 - u)^{(M - %d)}u^{%d}\n", u + 1, u, u, u);

  // Print v basis functions
  for(int v = 0; v <= ord_v; ++v)
    printf("B_{%dN}(v) = nCr\\left(N, %d\\right)(1 - v)^{(N - %d)}v^{%d}\n", v + 1, v, v, v);

  // ===================== Print weight ======================
  printf("w(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}(u)B_{%dN}(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nw_{u}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}'(u)B_{%dN}(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nw_{v}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}(u)B_{%dN}'(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nw_{uu}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}''(u)B_{%dN}(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
    
  printf("\nw_{vv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}(u)B_{%dN}''(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nw_{uv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}'(u)B_{%dN}'(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  // ===================== Print X ======================
  printf("\nX(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}(v)",  v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
  
  printf("\nX_{u}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nX_{v}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nX_{uu}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}''(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
    
  printf("\nX_{vv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}''(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nX_{uv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
 
  // ===================== Print Y ======================
  printf("\nY(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
  
  printf("\nY_{u}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nY_{v}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nY_{uu}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}''(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
    
  printf("\nY_{vv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}''(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nY_{uv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  // ===================== Print Z ======================
  printf("\nZ(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
  
  printf("\nZ_{u}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nZ_{v}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nZ_{uu}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}''(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
    
  printf("\nZ_{vv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}''(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nZ_{uv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nP(u, v) = [X(u, v), Y(u, v), Z(u, v)]");
  printf("\nP_{u}(u, v) = [X_{u}(u, v), Y_{u}(u, v), Z_{u}(u, v)]");
  printf("\nP_{v}(u, v) = [X_{v}(u, v), Y_{v}(u, v), Z_{v}(u, v)]");
  printf("\nP_{uu}(u, v) = [X_{uu}(u, v), Y_{uu}(u, v), Z_{uu}(u, v)]");
  printf("\nP_{vv}(u, v) = [X_{vv}(u, v), Y_{vv}(u, v), Z_{vv}(u, v)]");
  printf("\nP_{uv}(u, v) = [X_{uv}(u, v), Y_{uv}(u, v), Z_{uv}(u, v)]");

  printf("\nS(u, v) = \\frac{P(u, v)}{w(u, v)}");
  printf("\nS_u(u, v) = \\frac{P_{u}(u, v) - S(u, v)w_{u}(u, v)}{w(u, v)}");
  printf("\nS_v(u, v) = \\frac{P_{v}(u, v) - S(u, v)w_{v}(u, v)}{w(u, v)}");
  printf("\nS_{uu}(u, v) = \\frac{P_{uu}(u, v) - 2S_u(u, v)w_u(u, v) - S(u, v)w_{uu}(u, v)}{w(u, v)}");
  printf("\nS_{vv}(u, v) = \\frac{P_{vv}(u, v) - 2S_v(u, v)w_v(u, v) - S(u, v)w_{vv}(u, v)}{w(u, v)}");
  printf("\nS_{uv}(u, v) = \\frac{P_{uv}(u, v) - S_u(u, v)w_v(u, v) - S_v(u, v)w_u(u, v) - S(u, v)w_{uv}(u, v)}{w(u, v)}");

  printf("\ns = 0.0000000001");
  printf("\nu_0 = ");
  printf("\nv_0 = ");
  
  printf("\nS_{uu}(u_0, v_0)[1]");
  printf("\nS_{uu}(u_0, v_0)[2]");
  printf("\nS_{uu}(u_0, v_0)[3]");
  
  printf("\nS_{vv}(u_0, v_0)[1]");
  printf("\nS_{vv}(u_0, v_0)[2]");
  printf("\nS_{vv}(u_0, v_0)[3]");

  printf("\nS_{uv}(u_0, v_0)[1]");
  printf("\nS_{uv}(u_0, v_0)[2]");
  printf("\nS_{uv}(u_0, v_0)[3]");
  printf("\n\n");
  for(int v = 0; v <= ord_v; ++v)
  {
    printf("W_{%d} = [%f", v + 1, patch.getWeight(0, v));
    for(int u = 1; u <= ord_u; ++u)
    {
      printf(", %f", patch.getWeight(u, v));
    }
    printf("]\nX_{%d} = [%f", v + 1, patch(0, v)[0]);
    for(int u = 1; u <= ord_u; ++u)
    {
      printf(", %f", patch(u, v)[0]);
    }
    printf("]\nY_{%d} = [%f", v + 1, patch(0, v)[1]);
    for(int u = 1; u <= ord_u; ++u)
    {
      printf(", %f", patch(u, v)[1]);
    }
    printf("]\nZ_{%d} = [%f", v + 1, patch(0, v)[2]);
    for(int u = 1; u <= ord_u; ++u)
    {
      printf(", %f", patch(u, v)[2]);
    }
    printf("]\n");
  }
  // clang-format on
}

template <typename T>
void python_print(const OrientedBoundingBox<T, 3> oBox, const char* color = "blue")
{
  auto verts = oBox.vertices();
  printf("plot_box( fig, ax, [");
  printf("(%.17f,%.17f,%.17f)", verts[0][0], verts[0][1], verts[0][2]);
  for(int i = 1; i < verts.size(); ++i)
    printf(",(%.17f,%.17f,%.17f)", verts[i][0], verts[i][1], verts[i][2]);
  printf("], color='%s' )\n\n", color);
}

template <class Lambda, typename T>
void python_print(const OrientedBoundingBox<T, 3> oBox,
                  Lambda rotate_point,
                  numerics::Matrix<T> rotator,
                  const char* color = "blue")
{
  auto verts = oBox.vertices();
  printf("plot_box( fig, ax, [");
  auto rot_point = rotate_point(rotator, verts[0]);
  printf("(%.17f,%.17f,%.17f)", rot_point[0], rot_point[1], rot_point[2]);
  for(int i = 1; i < verts.size(); ++i)
  {
    rot_point = rotate_point(rotator, verts[i]);
    printf(",(%.17f,%.17f,%.17f)", rot_point[0], rot_point[1], rot_point[2]);
  }
  printf("], color='%s' )\n\n", color);
}

template <typename T>
void python_print(const BoundingBox<T, 3>& bBox, const char* color = "blue")
{
  std::vector<Point<T, 3>> verts;
  BoundingBox<T, 3>::getPoints(bBox, verts);
  printf("plot_bounding_box( fig, ax, [");
  printf("(%.17f,%.17f,%.17f)", verts[0][0], verts[0][1], verts[0][2]);
  for(int i = 1; i < verts.size(); ++i)
    printf(",(%.17f,%.17f,%.17f)", verts[i][0], verts[i][1], verts[i][2]);
  printf("], color='%s' )\n", color);
}

template <typename T>
void python_print(const Polygon<T, 3>& poly, bool closed = true)
{
  printf("plot_polygon(fig, ax, [");
  printf("(%.17f,%.17f,%.17f)", poly[0][0], poly[0][1], poly[0][2]);
  for(int i = 1; i < poly.numVertices(); ++i)
    printf(",(%.17f,%.17f,%.17f)", poly[i][0], poly[i][1], poly[i][2]);
  printf("]");
  if(!closed) printf(", closed=False");
  printf(")\n");
}

template <typename T>
void desmos_print(const Polygon<T, 2>& poly)
{
  printf("polygon(");
  printf("(%.17f,%.17f)", poly[0][0], poly[0][1]);
  for(int i = 1; i < poly.numVertices(); ++i) printf(",(%.17f,%.17f)", poly[i][0], poly[i][1]);
  printf(")\n");
}

template <typename T>
void convert_from_svg(std::string filename, axom::Array<BezierCurve<T, 2>>& curves)
{
  // Get regex pattern for all path elements
  std::regex pattern("<\\s*path[^>]*\\bd\\s*=\\s*[\"']([^\"']+)\"[^>]*>");

  // Read in file
  std::ifstream svg_file(filename);

  if(!svg_file.is_open())
  {
    std::cerr << "Failed to open the SVG file." << std::endl;
    return;  // Exit with an error code
  }

  std::string input((std::istreambuf_iterator<char>(svg_file)), std::istreambuf_iterator<char>());

  std::sregex_iterator iter(input.begin(), input.end(), pattern);
  std::sregex_iterator end;

  while(iter != end)
  {
    std::smatch match = *iter;
    std::string dAttribute = match[1].str();

    // Replace all ',' with ' '
    dAttribute = std::regex_replace(dAttribute, std::regex(","), std::string(" "));

    // Put a space before and after every letter
    dAttribute = std::regex_replace(dAttribute, std::regex("([a-zA-Z])"), std::string(" $1 "));

    char command = '\0';
    double init_x = 0.0, init_y = 0.0;
    double curr_x = 0.0, curr_y = 0.0;
    double x1, y1, x2, y2, x, y;
    double dx1, dy1, dx2, dy2, dx, dy;

    // Iterate over words in the string
    std::istringstream iss(dAttribute);
    std::string word;

    while(iss >> word)
    {
      // Check if the word is a command
      if(isalpha(word[0]))
      {
        command = word[0];
        if(command == 'Z' || command == 'z')
        {
          axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                            Point<T, 2> {init_x, init_y}};

          // Add a line segment
          curves.push_back(BezierCurve<T, 2>(nodes, 1));

          // Update the current position
          curr_x = init_x;
          curr_y = init_y;
        }
      }
      else
      {
        // Check if the command is a move command
        if(command == 'M' || command == 'm')
        {
          // If the command is relative, add the value to the current position
          if(command == 'm')
          {
            dx = std::stod(word);
            iss >> dy;
            curr_x += dx;
            curr_y += dy;
            command = 'l';
          }
          else
          {
            x = std::stod(word);
            iss >> y;
            curr_x = x;
            curr_y = y;
            command = 'L';
          }

          // Set the initial position
          init_x = curr_x;
          init_y = curr_y;
        }
        else if(command == 'L' || command == 'l')
        {
          // If the command is relative, add the value to the current position
          if(command == 'l')
          {
            dx = std::stod(word);
            iss >> dy;
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {curr_x + dx, curr_y + dy}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_x += dx;
            curr_y += dy;
          }
          else
          {
            x = std::stod(word);
            iss >> y;
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y}, Point<T, 2> {x, y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_x = x;
            curr_y = y;
          }
        }
        else if(command == 'C' || command == 'c')
        {
          // If the command is relative, add the value to the current position
          if(command == 'c')
          {
            dx1 = std::stod(word);
            iss >> dy1 >> dx2 >> dy2 >> dx >> dy;
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {curr_x + dx1, curr_y + dy1},
                                              Point<T, 2> {curr_x + dx2, curr_y + dy2},
                                              Point<T, 2> {curr_x + dx, curr_y + dy}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 3));

            // Update the current position
            curr_x += dx;
            curr_y += dy;
          }
          else
          {
            x1 = std::stod(word);
            iss >> y1 >> x2 >> y2 >> x >> y;
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {x1, y1},
                                              Point<T, 2> {x2, y2},
                                              Point<T, 2> {x, y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 3));

            // Update the current position
            curr_x = x;
            curr_y = y;
          }
        }
        else if(command == 'H' || command == 'h')
        {
          // If the command is relative, add the value to the current position
          if(command == 'h')
          {
            dx = std::stod(word);
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {curr_x + dx, curr_y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_x += dx;
          }
          else
          {
            x = std::stod(word);
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y}, Point<T, 2> {x, curr_y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_x = x;
          }
        }
        else if(command == 'V' || command == 'v')
        {
          // If the command is relative, add the value to the current position
          if(command == 'v')
          {
            dy = std::stod(word);
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {curr_x, curr_y + dy}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_y += dy;
          }
          else
          {
            y = std::stod(word);
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y}, Point<T, 2> {curr_x, y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_y = y;
          }
        }
        else
        {
          std::cout << word << " is an unrecognized command. Whoopsie!" << std::endl;
        }
      }
    }
    ++iter;
  }

  svg_file.close();
}

template <typename T>
BoundingBox<T, 2> curves_bbox(axom::Array<BezierCurve<T, 2>>& curves,
                              double scale_factor = 1.0,
                              bool make_square = false)
{
  // Find maximum and minimum x and y values
  double min_x = curves[0][0][0], min_y = curves[0][0][1];
  double max_x = curves[0][0][0], max_y = curves[0][0][1];

  for(auto& curve : curves)
  {
    for(int p = 0; p <= curve.getOrder(); ++p)
    {
      min_x = axom::utilities::min(min_x, curve[p][0]);
      max_x = axom::utilities::max(max_x, curve[p][0]);

      min_y = axom::utilities::min(min_y, curve[p][1]);
      max_y = axom::utilities::max(max_y, curve[p][1]);
    }
  }

  // Put those in a bounding box, then expand it
  primal::Point<double, 2> bounds[2] = {primal::Point<double, 2> {min_x, min_y},
                                        primal::Point<double, 2> {max_x, max_y}};
  primal::BoundingBox<double, 2> bb(bounds, 2);
  bb.scale(scale_factor);

  if(make_square)
  {
    int max_dim = bb.getLongestDimension();
    double max_len = bb.getMax()[max_dim] - bb.getMin()[max_dim];

    primal::Point<double, 2> centroid = bb.getCentroid();

    primal::Point<double, 2> new_min {centroid[0] - max_len / 2.0, centroid[1] - max_len / 2.0};
    primal::Point<double, 2> new_max {centroid[0] + max_len / 2.0, centroid[1] + max_len / 2.0};

    bb.addPoint(new_min);
    bb.addPoint(new_max);
  }

  return bb;
}

template <typename T>
void simple_grid_test(axom::Array<BezierCurve<T, 2>>& curves,
                      const BoundingBox<T, 2>& bb,
                      int npts_x,
                      int npts_y,
                      std::ofstream& wn_out)
{
  for(int xi = 0; xi < npts_x; ++xi)
  {
    double x = bb.getMin()[0] + xi * (bb.getMax()[0] - bb.getMin()[0]) / npts_x;

    //std::cout << x << std::endl;
    printLoadingBar((x - bb.getMin()[0]) / (bb.getMax()[0] - bb.getMin()[0]) * 100, 100);
    for(int yi = 0; yi < npts_y; ++yi)
    {
      double y = bb.getMin()[1] + yi * (bb.getMax()[1] - bb.getMin()[1]) / npts_y;
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn = 0.0;
      for(int i = 0; i < curves.size(); ++i)
      {
        wn += primal::winding_number(query, curves[i], 1e-16, 1e-16);
      }

      wn_out << x << "," << y << "," << wn << std::endl;
    }
  }
}

//template <typename T>
//void simple_grid_test_memoized(
//  axom::Array<std::unordered_map<std::pair<int, int>,
//                                       primal::detail::BezierCurveMemo<double>,
//                                       primal::detail::PairHash>>& marray,
//  const BoundingBox<T, 2>& bb,
//  int npts_x,
//  int npts_y,
//  std::ofstream& wn_out)
//{
//  primal::Polygon<double, 2> temp_approxogon(20);
//
//  for(double x = bb.getMin()[0]; x <= bb.getMax()[0];
//      x += (bb.getMax()[0] - bb.getMin()[0]) / npts_x)
//  {
//    //std::cout << x << std::endl;
//    printLoadingBar((x - bb.getMin()[0]) / (bb.getMax()[0] - bb.getMin()[0]) * 100,
//                    100);
//    for(double y = bb.getMin()[1]; y <= bb.getMax()[1];
//        y += (bb.getMax()[1] - bb.getMin()[1]) / npts_y)
//    {
//      primal::Point<double, 2> query({x, y});
//
//      int nevals = 0;
//
//      double wn = 0.0;
//
//      wn += primal::winding_number_approxogon_memoized(query,
//                                                       marray,
//                                                       temp_approxogon,
//                                                       1e-16);
//
//      wn_out << x << "," << y << "," << wn << std::endl;
//    }
//  }
//}

template <typename T>
void simple_timing_test(axom::Array<BezierCurve<T, 2>>& curves,
                        const BoundingBox<T, 2>& bb,
                        int npts_x,
                        int npts_y,
                        std::ofstream& wn_out)
{
  for(double x = bb.getMin()[0]; x <= bb.getMax()[0]; x += (bb.getMax()[0] - bb.getMin()[0]) / npts_x)
  {
    //std::cout << x << std::endl;
    printLoadingBar((x - bb.getMin()[0]) / (bb.getMax()[0] - bb.getMin()[0]) * 100, 100);
    for(double y = bb.getMin()[1]; y <= bb.getMax()[1];
        y += (bb.getMax()[1] - bb.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn = 0.0;
      for(int i = 0; i < curves.size(); ++i)
      {
        wn += primal::winding_number(query, curves[i], nevals, 1e-16, 1e-16);
      }

      wn_out << x << "," << y << "," << wn << std::endl;
    }
  }
}

template <typename T>
inline void printLoadingBar(T progress, int total, int barWidth = 40)
{
  float percentage = static_cast<float>(progress) / total;
  int progressWidth = static_cast<int>(percentage * barWidth);

  std::cout << "[";
  for(int i = 0; i < progressWidth; ++i) std::cout << "=";
  for(int i = progressWidth; i < barWidth; ++i) std::cout << "-";
  std::cout << "] " << std::setw(3) << static_cast<int>(percentage * 100.0) << "%\r";
  std::cout.flush();
}

template <typename T>
void exportScalarFieldToVTK(const std::string& filename,
                            std::function<T(Point3D)> scalarField,
                            primal::BoundingBox<T, 3> bbox,
                            int xSteps,
                            int ySteps,
                            int zSteps)
{
  std::ofstream file(filename);

  if(!file.is_open())
  {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  T dx = (xSteps > 1) ? (bbox.getMax()[0] - bbox.getMin()[0]) / (xSteps - 1) : 0.0;
  T dy = (ySteps > 1) ? (bbox.getMax()[1] - bbox.getMin()[1]) / (ySteps - 1) : 0.0;
  T dz = (zSteps > 1) ? (bbox.getMax()[2] - bbox.getMin()[2]) / (zSteps - 1) : 0.0;

  // Write VTK header
  file << "# vtk DataFile Version 3.0\n";
  file << "Scalar field data\n";
  file << "ASCII\n";
  file << "DATASET STRUCTURED_POINTS\n";
  file << "DIMENSIONS " << xSteps << " " << ySteps << " " << zSteps << "\n";
  file << "ORIGIN " << bbox.getMin()[0] << " " << bbox.getMin()[1] << " " << bbox.getMin()[2] << "\n";
  file << "SPACING " << dx << " " << dy << " " << dz << "\n";
  file << "POINT_DATA " << xSteps * ySteps * zSteps << "\n";
  file << "SCALARS scalars float\n";
  file << "LOOKUP_TABLE default\n";

  file << std::setprecision(20);

  // Write scalar field data
  for(int k = 0; k < zSteps; ++k)
  {
    printLoadingBar(k, zSteps);

    for(int j = 0; j < ySteps; ++j)
    {
      for(int i = 0; i < xSteps; ++i)
      {
        T x = bbox.getMin()[0] + i * dx;
        T y = bbox.getMin()[1] + j * dy;
        T z = bbox.getMin()[2] + k * dz;

        auto query = Point3D({x, y, z});

        T scalarValue = scalarField(query);

        // Print the query if scalarValue is nan
        if(scalarValue != scalarValue) std::cout << std::setprecision(20) << query << std::endl;

        //for(int n = 0; n < patches.size(); ++n)
        //scalarValue +=
        //winding_number(query, patches[n], edge_tol, quad_tol, EPS);

        file << scalarValue << "\n";
      }
    }
  }

  file.close();
}

template <typename T>
void exportScalarFieldToVTK(const std::string& filename,
                            const std::string& timing_filename,
                            std::function<T(Point3D)> scalarField,
                            primal::BoundingBox<T, 3> bbox,
                            int xSteps,
                            int ySteps,
                            int zSteps)
{
  std::ofstream file(filename);
  std::ofstream timing_file(timing_filename);

  if(!file.is_open())
  {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  T dx = (xSteps > 1) ? (bbox.getMax()[0] - bbox.getMin()[0]) / (xSteps - 1) : 0.0;
  T dy = (ySteps > 1) ? (bbox.getMax()[1] - bbox.getMin()[1]) / (ySteps - 1) : 0.0;
  T dz = (zSteps > 1) ? (bbox.getMax()[2] - bbox.getMin()[2]) / (zSteps - 1) : 0.0;

  // Write VTK header
  file << "# vtk DataFile Version 3.0\n";
  file << "Scalar field data\n";
  file << "ASCII\n";
  file << "DATASET STRUCTURED_POINTS\n";
  file << "DIMENSIONS " << xSteps << " " << ySteps << " " << zSteps << "\n";
  file << "ORIGIN " << bbox.getMin()[0] << " " << bbox.getMin()[1] << " " << bbox.getMin()[2] << "\n";
  file << "SPACING " << dx << " " << dy << " " << dz << "\n";
  file << "POINT_DATA " << xSteps * ySteps * zSteps << "\n";
  file << "SCALARS scalars float\n";
  file << "LOOKUP_TABLE default\n";

  int case_count[4] = {0, 0, 0, 0};
  double case_time[4] = {0.0, 0.0, 0.0, 0.0};

  // Write scalar field data
  for(int k = 0; k < zSteps; ++k)
  {
    printLoadingBar(k, zSteps);

    for(int j = 0; j < ySteps; ++j)
    {
      for(int i = 0; i < xSteps; ++i)
      {
        T x = bbox.getMin()[0] + i * dx;
        T y = bbox.getMin()[1] + j * dy;
        T z = bbox.getMin()[2] + k * dz;

        auto query = Point3D({x, y, z});

        T scalarValue = scalarField(query);

        // Print the query if scalarValue is nan
        if(scalarValue != scalarValue) std::cout << std::setprecision(20) << query << std::endl;

        //for(int n = 0; n < patches.size(); ++n)
        //scalarValue +=
        //winding_number(query, patches[n], edge_tol, quad_tol, EPS);

        file << scalarValue << "\n";
      }
    }
  }

  file.close();
}

template <typename T>
void exportSumSplitScalarFieldToVTK(const std::string& filename,
                                    std::function<std::pair<double, double>(Point3D)> scalarField,
                                    primal::BoundingBox<T, 3> bbox,
                                    int xSteps,
                                    int ySteps,
                                    int zSteps)
{
  // Arrays to store the scalar fields
  std::vector<double> scalarField1;
  std::vector<double> scalarField2;
  std::vector<double> scalarField3;

  // Reserve space for the arrays
  int totalPoints = xSteps * ySteps * zSteps;
  scalarField1.reserve(totalPoints);
  scalarField2.reserve(totalPoints);
  scalarField3.reserve(totalPoints);

  for(int k = 0; k < zSteps; ++k)
  {
    printLoadingBar(k, zSteps);

    for(int j = 0; j < ySteps; ++j)
    {
      for(int i = 0; i < xSteps; ++i)
      {
        double x = bbox.getMin()[0] + i * (bbox.getMax()[0] - bbox.getMin()[0]) / (xSteps - 1);
        double y = bbox.getMin()[1] + j * (bbox.getMax()[1] - bbox.getMin()[1]) / (ySteps - 1);
        double z = bbox.getMin()[2] + k * (bbox.getMax()[2] - bbox.getMin()[2]) / (zSteps - 1);

        auto query = Point3D({x, y, z});

        auto splitScalarValue = scalarField(query);

        scalarField1.push_back(splitScalarValue.first);
        scalarField2.push_back(splitScalarValue.second);
        scalarField3.push_back(splitScalarValue.first + splitScalarValue.second);
      }
    }
  }

  std::ofstream file(filename);

  if(!file.is_open())
  {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  // Write VTK header
  file << "# vtk DataFile Version 3.0\n";
  file << "Scalar field data\n";
  file << "ASCII\n";
  file << "DATASET STRUCTURED_POINTS\n";
  file << "DIMENSIONS " << xSteps << " " << ySteps << " " << zSteps << "\n";
  file << "ORIGIN " << bbox.getMin()[0] << " " << bbox.getMin()[1] << " " << bbox.getMin()[2] << "\n";
  file << "SPACING " << (bbox.getMax()[0] - bbox.getMin()[0]) / (xSteps - 1) << " "
       << (bbox.getMax()[1] - bbox.getMin()[1]) / (ySteps - 1) << " "
       << (bbox.getMax()[2] - bbox.getMin()[2]) / (zSteps - 1) << "\n";
  file << "POINT_DATA " << xSteps * ySteps * zSteps << "\n";

  // Write scalar field data
  file << "SCALARS Stokes double 1\n";
  file << "LOOKUP_TABLE default\n";
  for(int i = 0; i < totalPoints; ++i)
  {
    file << scalarField1[i] << "\n";
  }

  file << "SCALARS JumpCondition double 1\n";
  file << "LOOKUP_TABLE default\n";
  for(int i = 0; i < totalPoints; ++i)
  {
    file << scalarField2[i] << "\n";
  }

  file << "SCALARS GWN double 1\n";
  file << "LOOKUP_TABLE default\n";
  for(int i = 0; i < totalPoints; ++i)
  {
    file << scalarField3[i] << "\n";
  }

  file.close();
}

template <typename T>
void exportSliceScalarFieldToVTK(const std::string& filename,
                                 const std::function<T(Point3D)>& fieldFunc,
                                 const Point3D& origin,
                                 const Vector<T, 3>& normal,
                                 double planeWidth,
                                 double planeHeight,
                                 int uSteps,
                                 int vSteps)
{
  // Define two orthogonal vectors to normal
  Vector<T, 3> u;
  if(axom::utilities::isNearlyEqual(normal[0], normal[1]))
  {
    u = Vector<T, 3>({normal[2], normal[2], -normal[0] - normal[1]}).unitVector();
    // Pick a u that is orthogonal to this one
    // u = Vector<T, 3>({normal[1], -normal[0], 0}).unitVector();
    // u = Vector<T, 3>({0, normal[2], -normal[1]}).unitVector();
  }
  else
  {
    u = Vector<T, 3>({-normal[1] - normal[2], normal[0], normal[0]}).unitVector();
  }
  // u = Vector<T, 3>( {-0.6536123100031038, -0.10738450531201763, -0.7491725543766935});
  Vector<T, 3> v = Vector<T, 3>::cross_product(normal, u).unitVector();

  exportSliceScalarFieldToVTK<T>(filename, fieldFunc, origin, u, v, planeWidth, planeHeight, uSteps, vSteps);
}

template <typename T>
void exportSliceScalarFieldToVTK(const std::string& filename,
                                 const std::function<T(Point3D)>& fieldFunc,
                                 const Point3D& origin,
                                 Vector<T, 3> u,
                                 Vector<T, 3> v,
                                 double planeWidth,
                                 double planeHeight,
                                 int uSteps,
                                 int vSteps)
{
  u = u.unitVector();
  v = v.unitVector();

  std::ofstream file(filename);
  if(!file.is_open())
  {
    std::cerr << "Failed to open file for writing: " << filename << std::endl;
    return;
  }

  file << "# vtk DataFile Version 3.0\n";
  file << "Scalar field on a plane\n";
  file << "ASCII\n";
  file << "DATASET STRUCTURED_GRID\n";
  file << "DIMENSIONS " << uSteps << " " << vSteps << " 1\n";

  // Write points
  file << "POINTS " << (uSteps * vSteps) << " float\n";
  for(int j = 0; j < vSteps; ++j)
  {
    for(int i = 0; i < uSteps; ++i)
    {
      double s = planeWidth * (i / static_cast<double>(uSteps - 1) - 0.5);
      double t = planeHeight * (j / static_cast<double>(vSteps - 1) - 0.5);
      double x = origin[0] + s * u[0] + t * v[0];
      double y = origin[1] + s * u[1] + t * v[1];
      double z = origin[2] + s * u[2] + t * v[2];
      file << x << " " << y << " " << z << "\n";
    }
  }

  file << std::setprecision(20);

  // Write scalar field
  file << "POINT_DATA " << (uSteps * vSteps) << "\n";
  file << "SCALARS scalar_field float 1\n";
  file << "LOOKUP_TABLE default\n";
  for(int j = 0; j < vSteps; ++j)
  {
    printLoadingBar(j, vSteps);
    for(int i = 0; i < uSteps; ++i)
    {
      double s = planeWidth * (i / static_cast<double>(uSteps - 1) - 0.5);
      double t = planeHeight * (j / static_cast<double>(vSteps - 1) - 0.5);
      double x = origin[0] + s * u[0] + t * v[0];
      double y = origin[1] + s * u[1] + t * v[1];
      double z = origin[2] + s * u[2] + t * v[2];

      auto val = fieldFunc(Point3D({x, y, z}));
      if(val != val)
      {
        std::cout << std::setprecision(20);
        std::cout << "NAN recorded at (" << x << " " << y << " " << z << ")" << std::endl;
        std::cout << std::setprecision(6);
        val = 0.0;
      }
      file << val << "\n";
    }
  }

  file.close();
}

template <typename T>
void exportSplitScalarSliceFieldToVTK(const std::string& filename,
                                      std::function<std::pair<double, double>(Point3D)> scalarField,
                                      const Point3D& origin,
                                      Vector<T, 3>& u,
                                      Vector<T, 3>& v,
                                      double planeWidth,
                                      double planeHeight,
                                      int uSteps,
                                      int vSteps)
{
  u = u.unitVector();
  v = v.unitVector();

  std::ofstream file(filename);

  // Arrays to store the scalar fields
  std::vector<double> scalarField1, scalarField2;

  // Reserve space for the arrays
  int totalPoints = uSteps * vSteps;
  scalarField1.reserve(totalPoints);
  scalarField2.reserve(totalPoints);

  for(int j = 0; j < vSteps; ++j)
  {
    printLoadingBar(j, vSteps);
    for(int i = 0; i < uSteps; ++i)
    {
      double s = planeWidth * (i / static_cast<double>(uSteps - 1) - 0.5);
      double t = planeHeight * (j / static_cast<double>(vSteps - 1) - 0.5);
      double x = origin[0] + s * u[0] + t * v[0];
      double y = origin[1] + s * u[1] + t * v[1];
      double z = origin[2] + s * u[2] + t * v[2];

      auto val = scalarField(Point3D({x, y, z}));
      if(val.first != val.first)
      {
        std::cout << std::setprecision(20);
        std::cout << "NAN recorded at (" << x << " " << y << " " << z << ")" << std::endl;
        std::cout << std::setprecision(6);
        val.first = 0.0;
      }
      scalarField1.push_back(val.first);
      scalarField2.push_back(val.second);
    }
  }

  if(!file.is_open())
  {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  // Write VTK header
  file << "# vtk DataFile Version 3.0\n";
  file << "Scalar field on a plane\n";
  file << "ASCII\n";
  file << "DATASET STRUCTURED_GRID\n";
  file << "DIMENSIONS " << uSteps << " " << vSteps << " 1\n";

  // Write points
  file << "POINTS " << (uSteps * vSteps) << " float\n";
  for(int j = 0; j < vSteps; ++j)
  {
    for(int i = 0; i < uSteps; ++i)
    {
      double s = planeWidth * (i / static_cast<double>(uSteps - 1) - 0.5);
      double t = planeHeight * (j / static_cast<double>(vSteps - 1) - 0.5);
      double x = origin[0] + s * u[0] + t * v[0];
      double y = origin[1] + s * u[1] + t * v[1];
      double z = origin[2] + s * u[2] + t * v[2];
      file << x << " " << y << " " << z << "\n";
    }
  }

  // Write scalar field data
  file << "POINT_DATA " << (uSteps * vSteps) << "\n";
  file << "SCALARS gwn double 1\n";
  file << "LOOKUP_TABLE default\n";
  for(int i = 0; i < totalPoints; ++i)
  {
    file << scalarField1[i] << "\n";
  }

  file << "SCALARS timing double 1\n";
  file << "LOOKUP_TABLE default\n";
  for(int i = 0; i < totalPoints; ++i)
  {
    file << scalarField2[i] << "\n";
  }

  file.close();
}

template <typename T>
void exportSurfaceToSTL(const std::string& filename,
                        const axom::Array<primal::BezierPatch<T, 3>>& patches,
                        int uSteps = 17,
                        int vSteps = 17)
{
  std::ofstream file(filename);

  if(!file.is_open())
  {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  // Write STL header
  file << "solid parametric_surface\n";
  file << std::setprecision(20);

  // Generate surface data
  for(const auto& patch : patches)
  {
    for(int j = 0; j < vSteps - 1; ++j)
    {
      for(int i = 0; i < uSteps - 1; ++i)
      {
        double u0 = 1.0 * i / (uSteps - 1);
        double v0 = 1.0 * j / (vSteps - 1);
        double u1 = 1.0 * (i + 1) / (uSteps - 1);
        double v1 = 1.0 * (j + 1) / (vSteps - 1);

        auto pt0 = patch.evaluate(u0, v0);
        auto pt1 = patch.evaluate(u1, v0);
        auto pt2 = patch.evaluate(u0, v1);
        auto pt3 = patch.evaluate(u1, v1);

        // First triangle
        file << "  facet normal 0 0 0\n";
        file << "    outer loop\n";
        file << "      vertex " << pt0[0] << " " << pt0[1] << " " << pt0[2] << "\n";
        file << "      vertex " << pt1[0] << " " << pt1[1] << " " << pt1[2] << "\n";
        file << "      vertex " << pt2[0] << " " << pt2[1] << " " << pt2[2] << "\n";
        file << "    endloop\n";
        file << "  endfacet\n";

        // Second triangle
        file << "  facet normal 0 0 0\n";
        file << "    outer loop\n";
        file << "      vertex " << pt2[0] << " " << pt2[1] << " " << pt2[2] << "\n";
        file << "      vertex " << pt1[0] << " " << pt1[1] << " " << pt1[2] << "\n";
        file << "      vertex " << pt3[0] << " " << pt3[1] << " " << pt3[2] << "\n";
        file << "    endloop\n";
        file << "  endfacet\n";
      }
    }
  }

  file << "endsolid parametric_surface\n";
  file.close();
}

template <typename T>
void exportSurfaceToSTL(const std::string& filename,
                        const primal::BezierPatch<T, 3>& patches,
                        int uSteps = 17,
                        int vSteps = 17)
{
  axom::Array<primal::BezierPatch<T, 3>> patchArray;
  patchArray.push_back(patches);
  exportSurfaceToSTL(filename, patchArray, uSteps, vSteps);
}

template <typename T>
void exportSurfaceToSTL(const std::string& filename,
                        const primal::NURBSPatch<T, 3>& patch,
                        int u_steps = 17,
                        int v_steps = 17)
{
  auto beziers = patch.extractBezier();
  exportSurfaceToSTL(filename, beziers, u_steps, v_steps);
}

template <typename T>
void lamePrinter(const primal::NURBSPatch<T, 3>& patch)
{
  const int npts = 20;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(patch.getMinKnot_u(), patch.getMaxKnot_u(), u_pts, npts);
  axom::numerics::linspace(patch.getMinKnot_v(), patch.getMaxKnot_v(), v_pts, npts);

  std::cout << std::fixed;
  std::cout << "==========================" << std::endl;
  for(int j = 0; j < npts; ++j)
  {
    for(int i = 0; i < npts; ++i)
    {
      if(patch.isVisible(u_pts[i], v_pts[j]))
      {
        auto pt = patch.evaluate(u_pts[i], v_pts[j]);
        std::cout << "(" << pt[0] << "," << pt[1] << "," << pt[2] << "), ";
      }
    }
  }
  std::cout << "==========================" << std::endl;
}

template <typename T>
void printPatchBoundaries(const primal::BezierPatch<T, 3>& patch)
{
  printPatchBoundaries<T>(primal::NURBSPatch<T, 3>(patch));
}

template <typename T>
void printPatchBoundaries(const primal::NURBSPatch<T, 3>& patch, bool printUntrimmedBoundary = false)
{
  const int n_pts = 20;

  auto min_u = patch.getMinKnot_u();
  auto max_u = patch.getMaxKnot_u();
  auto min_v = patch.getMinKnot_v();
  auto max_v = patch.getMaxKnot_v();

  std::cout << "=============" << std::endl;
  std::cout << std::fixed << std::endl;

  if(!patch.isTrimmed())
  {
    double u_pts[n_pts], v_pts[n_pts];
    axom::numerics::linspace(min_u, max_u, u_pts, n_pts);
    axom::numerics::linspace(min_v, max_v, v_pts, n_pts);

    // Print the u/v = 0/1 isocurves
    for(int i = 0; i < n_pts; ++i)
    {
      auto pt1 = patch.evaluate(u_pts[i], min_v);
      auto pt2 = patch.evaluate(u_pts[i], max_v);
      auto pt3 = patch.evaluate(min_u, v_pts[i]);
      auto pt4 = patch.evaluate(max_u, v_pts[i]);

      std::cout << "(" << pt1[0] << ", " << pt1[1] << ", " << pt1[2] << "), ";
      std::cout << "(" << pt2[0] << ", " << pt2[1] << ", " << pt2[2] << "), ";
      std::cout << "(" << pt3[0] << ", " << pt3[1] << ", " << pt3[2] << "), ";
      std::cout << "(" << pt4[0] << ", " << pt4[1] << ", " << pt4[2] << "), ";
    }

    std::cout << std::endl;
  }
  else
  {
    for(auto& curve : patch.getTrimmingCurves())
    {
      double t_pts[n_pts];
      axom::numerics::linspace(curve.getMinKnot(), curve.getMaxKnot(), t_pts, n_pts);

      for(int i = 0; i < n_pts; ++i)
      {
        auto ppt = curve.evaluate(t_pts[i]);
        auto spt = patch.evaluate(ppt[0], ppt[1]);
        std::cout << "(" << spt[0] << ", " << spt[1] << ", " << spt[2] << "), ";
      }
    }

    std::cout << std::endl;
  }

  std::cout << "=============" << std::endl;

  // Print the untrimmed patch boundaries
  if(printUntrimmedBoundary)
  {
    double u_pts[n_pts], v_pts[n_pts];
    axom::numerics::linspace(min_u, max_u, u_pts, n_pts);
    axom::numerics::linspace(min_v, max_v, v_pts, n_pts);

    // Print the u/v = 0/1 isocurves
    for(int i = 0; i < n_pts; ++i)
    {
      auto pt1 = patch.evaluate(u_pts[i], min_v);
      auto pt2 = patch.evaluate(u_pts[i], max_v);
      auto pt3 = patch.evaluate(min_u, v_pts[i]);
      auto pt4 = patch.evaluate(max_u, v_pts[i]);

      std::cout << "(" << pt1[0] << ", " << pt1[1] << ", " << pt1[2] << "), ";
      std::cout << "(" << pt2[0] << ", " << pt2[1] << ", " << pt2[2] << "), ";
      std::cout << "(" << pt3[0] << ", " << pt3[1] << ", " << pt3[2] << "), ";
      std::cout << "(" << pt4[0] << ", " << pt4[1] << ", " << pt4[2] << "), ";
    }

    std::cout << "=============" << std::endl;
  }
}

}  // namespace primal

}  // namespace axom

#endif
