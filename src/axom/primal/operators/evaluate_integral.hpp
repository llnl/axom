// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file evaluate_integral.hpp
 *
 * \brief Consists of methods that evaluate integrals over regions defined
 * by Bezier and NURBS curves, such as 2D area integrals and scalar/vector field line integrals
 *
 * Line integrals are evaluated numerically with Gauss-Legendre quadrature
 *
 * 2D area integrals computed with "Spectral Mesh-Free Quadrature for Planar 
 * Regions Bounded by Rational Parametric Curves" by David Gunderman et al.
 */

#ifndef PRIMAL_EVAL_INTEGRAL_HPP_
#define PRIMAL_EVAL_INTEGRAL_HPP_

// Axom includes
#include "axom/config.hpp"

#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"

#include "axom/primal/operators/detail/evaluate_integral_impl.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
/*!
 * \brief Evaluate a line integral along the boundary of a CurvedPolygon object 
 *        on a scalar field.
 *
 * The line integral is evaluated on each curve in the CurvedPolygon, and added
 * together to represent the total integral. The Polygon need not be connected.
 * 
 * Evaluate the scalar field line integral with Gauss-Legendre quadrature
 * 
 * \param [in] cpoly the CurvedPolygon object
 * \param [in] scalar_integrand the lambda function representing the integrand. 
 * Must accept a Point<T, NDIM> as input and return a double
 * \param [in] npts the number of quadrature points to evaluate the line integral
 *                  on each edge of the CurvedPolygon
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
double evaluate_scalar_line_integral(const primal::CurvedPolygon<T, NDIMS> cpoly,
                                     Lambda&& scalar_integrand,
                                     int npts)
{
  double total_integral = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    // Compute the line integral along each component.
    total_integral +=
      detail::evaluate_scalar_line_integral_component(cpoly[i], scalar_integrand, npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a line integral on a single Bezier curve on a scalar field
 *
 * \param [in] c the Bezier curve object
 * \param [in] scalar_integrand the lambda function representing the integrand. 
 * Must accept a Point<T, NDIMS> as input, and return a double.
 * \param [in] npts the number of quadrature nodes
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
double evaluate_scalar_line_integral(const primal::BezierCurve<T, NDIMS>& c,
                                     Lambda&& scalar_integrand,
                                     int npts)
{
  return detail::evaluate_scalar_line_integral_component(c, scalar_integrand, npts);
}

/*!
 * \brief Evaluate a line integral on a single NURBS curve on a scalar field
 *
 * \param [in] n the NURBS curve object
 * \param [in] scalar_integrand the lambda function representing the integrand. 
 * Must accept a Point<T, NDIMS> as input, and return a double.
 * \param [in] npts the number of quadrature nodes per knot span
 * 
 * \note The NURBS curve is decomposed into Bezier segments, and the Gaussian quadrature
 *   is computed using npts on each segment. 
 *   
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
double evaluate_scalar_line_integral(const primal::NURBSCurve<T, NDIMS>& n,
                                     Lambda&& scalar_integrand,
                                     int npts)
{
  // Compute the line integral along each component.
  auto beziers = n.extractBezier();
  double total_integral = 0.0;
  for(int i = 0; i < beziers.size(); i++)
  {
    total_integral +=
      detail::evaluate_scalar_line_integral_component(beziers[i], scalar_integrand, npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a line integral on an array of NURBS curves on a scalar field
 *
 * \param [in] narray the array of NURBS curve object
 * \param [in] scalar_integrand the lambda function representing the integrand. 
 * Must accept a Point<T, NDIMS> as input, and return a double.
 * \param [in] npts the number of quadrature nodes per curve per knot span
 * 
 * \note Each NURBS curve is decomposed into Bezier segments, and the Gaussian quadrature
 *   is computed using npts on each segment
 * 
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
double evaluate_scalar_line_integral(const axom::Array<primal::NURBSCurve<T, NDIMS>>& narray,
                                     Lambda&& scalar_integrand,
                                     int npts)
{
  double total_integral = 0.0;
  for(int i = 0; i < narray.size(); i++)
  {
    total_integral += evaluate_scalar_line_integral(narray[i], scalar_integrand, npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a line integral along the boundary of a CurvedPolygon object 
 *        on a vector field.
 *
 * The line integral is evaluated on each curve in the CurvedPolygon, and added
 * together to represent the total integral. The Polygon need not be connected.
 *
 * Evaluate the vector field line integral with Gauss-Legendre quadrature
 * 
 * \param [in] cpoly the CurvedPolygon object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * Must accept a Point<T, NDIM> as input and return a Vector<double, NDIM>
 * \param [in] npts the number of quadrature points to evaluate the line integral
 *                  on each edge of the CurvedPolygon
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
double evaluate_vector_line_integral(const primal::CurvedPolygon<T, NDIMS> cpoly,
                                     Lambda&& vector_integrand,
                                     int npts)
{
  double total_integral = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    // Compute the line integral along each component.
    total_integral +=
      detail::evaluate_vector_line_integral_component(cpoly[i], vector_integrand, npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a line integral on a single Bezier curve on a vector field
 *
 * \param [in] c the Bezier curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * Must accept a Point<T, NDIMS> as input, and return a Vector<double, NDIMS>.
 * \param [in] npts the number of quadrature nodes
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
double evaluate_vector_line_integral(const primal::BezierCurve<T, NDIMS>& c,
                                     Lambda&& vector_integrand,
                                     int npts)
{
  return detail::evaluate_vector_line_integral_component(c, vector_integrand, npts);
}

/*!
 * \brief Evaluate a line integral on a single NURBS curve on a vector field
 *
 * \param [in] n the NURBS curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * Must accept a Point<T, NDIMS> as input, and return a Vector<double, NDIMS>.
 * \param [in] npts the number of quadrature nodes per knot span
 * 
 * \note The NURBS curve is decomposed into Bezier segments, and the Gaussian quadrature
 *   is computed using npts on each segment
 * 
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
double evaluate_vector_line_integral(const primal::NURBSCurve<T, NDIMS>& n,
                                     Lambda&& vector_integrand,
                                     int npts)
{
  // Compute the line integral along each component.
  auto beziers = n.extractBezier();
  double total_integral = 0.0;
  for(int i = 0; i < beziers.size(); i++)
  {
    total_integral +=
      detail::evaluate_vector_line_integral_component(beziers[i], vector_integrand, npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a line integral on an array of NURBS curves on a vector field
 *
 * \param [in] narray the array of NURBS curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * Must accept a Point<T, NDIMS> as input and return a Vector<double, NDIMS>.
 * \param [in] npts the number of quadrature nodes per curve per knot span
 * 
 * \note Each NURBS curve is decomposed into Bezier segments, and the Gaussian quadrature
 *   is computed using npts on each segment
 *
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
double evaluate_vector_line_integral(const axom::Array<primal::NURBSCurve<T, NDIMS>>& narray,
                                     Lambda&& vector_integrand,
                                     int npts)
{
  double total_integral = 0.0;
  for(int i = 0; i < narray.size(); i++)
  {
    total_integral += evaluate_vector_line_integral(narray[i], vector_integrand, npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate an integral on the interior of a CurvedPolygon object.
 *
 * See above definition for details.
 * 
 * \param [in] cs the array of Bezier curve objects that bound the region
 * \param [in] integrand the lambda function representing the integrand. 
 * Must accept a 2D point as input and return a double
 * \param [in] npts_Q the number of quadrature points to evaluate the line integral
 * \param [in] npts_P the number of quadrature points to evaluate the antiderivative
 * \return the value of the integral
 */
template <class Lambda, typename T>
double evaluate_area_integral(const primal::CurvedPolygon<T, 2> cpoly,
                              Lambda&& integrand,
                              int npts_Q,
                              int npts_P = 0)
{
  if(npts_P <= 0)
  {
    npts_P = npts_Q;
  }

  // Use minimum y-coord of control nodes as lower bound for integration
  double int_lb = cpoly[0][0][1];
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    for(int j = 1; j < cpoly[i].getOrder() + 1; j++)
    {
      int_lb = std::min(int_lb, cpoly[i][j][1]);
    }
  }

  // Evaluate the antiderivative line integral along each component
  double total_integral = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    total_integral +=
      detail::evaluate_area_integral_component(cpoly[i], integrand, int_lb, npts_Q, npts_P);
  }

  return total_integral;
}

/*!
 * \brief Evaluate an integral on the interior of an array of NURBS curves.
 *
 * See above definition for details.
 * 
 * \param [in] narray the array of NURBS curve objects that bound the region
 * \param [in] integrand the lambda function representing the integrand. 
 * Must accept a 2D point as input and return a double
 * \param [in] npts_Q the number of quadrature points to evaluate the line integral
 * \param [in] npts_P the number of quadrature points to evaluate the antiderivative
 *  
 * \note Each NURBS curve is decomposed into Bezier segments, and the Gaussian quadrature
 *   is computed using npts_Q * npts_P on each segment
 * 
 * \note The numerical result is only meaningful if the curves enclose a region
 * 
 * \return the value of the integral
 */
template <class Lambda, typename T>
double evaluate_area_integral(const axom::Array<primal::NURBSCurve<T, 2>> narray,
                              Lambda&& integrand,
                              int npts_Q,
                              int npts_P = 0)
{
  if(npts_P <= 0)
  {
    npts_P = npts_Q;
  }

  if(narray.empty())
  {
    return 0.0;
  }

  // Use minimum y-coord of control nodes as lower bound for integration
  double int_lb = narray[0][0][1];
  for(int i = 0; i < narray.size(); i++)
  {
    for(int j = 1; j < narray[i].getNumControlPoints(); j++)
    {
      int_lb = std::min(int_lb, narray[i][j][1]);
    }
  }

  // Evaluate the antiderivative line integral along each component
  double total_integral = 0.0;
  for(int i = 0; i < narray.size(); i++)
  {
    auto beziers = narray[i].extractBezier();

    for(int j = 0; j < beziers.size(); j++)
    {
      total_integral +=
        detail::evaluate_area_integral_component(beziers[j], integrand, int_lb, npts_Q, npts_P);
    }
  }

  return total_integral;
}

}  // namespace primal
}  // end namespace axom

#endif
