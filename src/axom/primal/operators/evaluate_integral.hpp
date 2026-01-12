// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file evaluate_integral.hpp
 *
 * \brief Consists of methods that evaluate scalar-field integrals on curves and 
 *  regions defined by 2D curves, and vector-field integrals on curves
 *
 * All integrals are evaluated numerically with Gauss-Legendre quadrature
 * 
 * Scalar-field line integrals and scalar-field area integrals are of form 
 * int_D f(x) dr, with f : R^n -> R^m, D is a curve or a 2D region bound by curves
 * 
 * Vector-field line integrals are of form int_C f(x) \cdot d\vec{r}, 
 *  with f : R^n -> R^n, C is a curve
 * 
 * 2D area integrals computed with "Spectral Mesh-Free Quadrature for Planar 
 * Regions Bounded by Rational Parametric Curves" by David Gunderman et al.
 */

#ifndef PRIMAL_EVAL_INTEGRAL_HPP_
#define PRIMAL_EVAL_INTEGRAL_HPP_

// Axom includes
#include "axom/config.hpp"

#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/detail/evaluate_integral_impl.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
///@{
/// \name Evaluates scalar-field line integrals for functions f : R^n -> R^m

/*!
 * \brief Evaluate a line integral along the boundary of a CurvedPolygon object 
 *        for a function with an arbitrary return type.
 *
 * The line integral is evaluated on each curve in the CurvedPolygon, and added
 * together to represent the total integral. The curved polygon need not be connected.
 * 
 * Evaluate the line integral with Gauss-Legendre quadrature
 * 
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \tparam Lambda A callable type taking a CurveType's PointType and returning an integrable type
 * \tparam LambdaRetType A type which supports addition and scalar multiplication
 * \param [in] cpoly the CurvedPolygon object
 * \param [in] integrand the lambda function representing the integrand. 
 * \param [in] npts the number of quadrature points to evaluate the line integral
 *                  on each edge of the CurvedPolygon
 * \return the value of the integral
 */
template <typename CurveType,
          typename Lambda,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename CurveType::PointType>>
LambdaRetType evaluate_line_integral(const primal::CurvedPolygon<CurveType> cpoly,
                                     Lambda&& integrand,
                                     int npts)
{
  static_assert(
    detail::internal::is_integrable_v<typename CurveType::NumericType, LambdaRetType>,
    "evaluate_integral methods require addition and scalar multiplication for lambda function "
    "return type");

  LambdaRetType total_integral = LambdaRetType {};
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    // Compute the line integral along each component.
    total_integral +=
      detail::evaluate_line_integral_component(cpoly[i], std::forward<Lambda>(integrand), npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a line integral along the boundary of a generic curve
 *        for a function with an arbitrary return type.
 *
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \tparam Lambda A callable type taking a CurveType's PointType and returning an integrable type
 * \tparam LambdaRetType A type which supports addition and scalar multiplication
 * \param [in] c the generic curve object
 * \param [in] integrand the lambda function representing the integrand. 
 * \param [in] npts the number of quadrature nodes
 * \return the value of the integral
 */
template <typename CurveType,
          typename Lambda,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename CurveType::PointType>>
LambdaRetType evaluate_line_integral(const CurveType& c, Lambda&& integrand, int npts)
{
  static_assert(
    detail::internal::is_integrable_v<typename CurveType::NumericType, LambdaRetType>,
    "evaluate_integral methods require addition and scalar multiplication for lambda function "
    "return type");

  return detail::evaluate_line_integral_component(c, std::forward<Lambda>(integrand), npts);
}

/*!
 * \brief Evaluate a line integral on an array of NURBS curves on a scalar field
 *        for a function with an arbitrary return type.
 *
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \tparam Lambda A callable type taking a CurveType's PointType and returning an integrable type
 * \tparam LambdaRetType A type which supports addition and scalar multiplication
 * \param [in] carray The array of generic curve objects
 * \param [in] integrand the lambda function representing the integrand. 
 * \param [in] npts the number of quadrature nodes per curve per knot span
 * 
 * \note Each NURBS curve is decomposed into Bezier segments, and the Gaussian quadrature
 *   is computed using npts on each segment
 * 
 * \return the value of the integral
 */
template <typename CurveType,
          typename Lambda,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename CurveType::PointType>>
LambdaRetType evaluate_line_integral(const axom::Array<CurveType>& carray, Lambda&& integrand, int npts)
{
  static_assert(
    detail::internal::is_integrable_v<typename CurveType::NumericType, LambdaRetType>,
    "evaluate_integral methods require addition and scalar multiplication for lambda function "
    "return type");

  LambdaRetType total_integral = LambdaRetType {};
  for(int i = 0; i < carray.size(); i++)
  {
    total_integral +=
      detail::evaluate_line_integral_component(carray[i], std::forward<Lambda>(integrand), npts);
  }

  return total_integral;
}
//@}

///@{
/// \name Evaluates vector-field line integrals for functions f : R^n -> R^n

/*!
 * \brief Evaluate a vector-field line integral along the boundary of a CurvedPolygon object 
 *
 * The line integral is evaluated on each curve in the CurvedPolygon, and added
 * together to represent the total integral. The Polygon need not be connected.
 *
 * Evaluate the vector field line integral with Gauss-Legendre quadrature
 *
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \tparam Lambda A callable type taking a CurveType's PointType and returning its numeric type
 * \tparam FuncRetType The CurveType's numeric type
 * \param [in] cpoly the CurvedPolygon object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * \param [in] npts the number of quadrature points to evaluate the line integral
 *                  on each edge of the CurvedPolygon
 * \pre Lambda must return the CurveTypes's vector type
 * \return the value of the integral
 */
template <typename CurveType, typename Lambda, typename FuncRetType = typename CurveType::NumericType>
FuncRetType evaluate_vector_line_integral(const CurvedPolygon<CurveType> cpoly,
                                          Lambda&& vector_integrand,
                                          int npts)
{
  FuncRetType total_integral = FuncRetType {};
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    // Compute the line integral along each component.
    total_integral +=
      detail::evaluate_vector_line_integral_component(cpoly[i],
                                                      std::forward<Lambda>(vector_integrand),
                                                      npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a vector-field line integral on a single generic curve
 *
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \tparam Lambda A callable type taking a CurveType's PointType and returning its numeric type
 * \tparam FuncRetType The CurveType's numeric type
 * \param [in] c the generic curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * \param [in] npts the number of quadrature nodes
 * 
 * \pre Lambda must return the CurveTypes's vector type
 * \return the value of the integral
 */
template <typename CurveType, typename Lambda, typename FuncRetType = typename CurveType::NumericType>
FuncRetType evaluate_vector_line_integral(const CurveType& c, Lambda&& vector_integrand, int npts)
{
  return detail::evaluate_vector_line_integral_component(c,
                                                         std::forward<Lambda>(vector_integrand),
                                                         npts);
}

/*!
 * \brief Evaluate a line integral on an array of generic curves on a vector field
 *
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \tparam Lambda A callable type taking a CurveType's PointType and returning its numeric type
 * \tparam FuncRetType The CurveType's numeric type
 * \param [in] carray The array of generic curve objects
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * \param [in] npts the number of quadrature nodes per curve per knot span
 * 
 * \note Each NURBS curve is decomposed into Bezier segments, and the Gaussian quadrature
 *   is computed using npts on each segment
 *
 * \return the value of the integral
 */
template <typename CurveType, typename Lambda, typename FuncRetType = typename CurveType::NumericType>
FuncRetType evaluate_vector_line_integral(const axom::Array<CurveType>& carray,
                                          Lambda&& vector_integrand,
                                          int npts)
{
  FuncRetType total_integral = FuncRetType {};
  for(int i = 0; i < carray.size(); i++)
  {
    total_integral +=
      detail::evaluate_vector_line_integral_component(carray[i],
                                                      std::forward<Lambda>(vector_integrand),
                                                      npts);
  }

  return total_integral;
}
//@}

///@{
/// \name Evaluates scalar-field 2D area integrals for functions f : R^2 -> R^m

/*!
 * \brief Evaluate an integral on the interior of a CurvedPolygon object.
 *
 * Evaluates the integral using a Spectral Mesh-Free Quadrature derived from 
 * Green's theorem, evaluating the area integral as a line integral of the 
 * antiderivative over each component curve.
 * 
 * For algorithm details, see "Spectral Mesh-Free Quadrature for Planar 
 * Regions Bounded by Rational Parametric Curves" by David Gunderman et al.
 * 
 * \tparam Lambda A callable type taking a CurveType's PointType and returning an integrable type
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the geometry
 * \tparam LambdaRetType A type which supports addition and scalar multiplication
 * \param [in] cpoly the CurvedPolygon object
 * \param [in] integrand the lambda function representing the integrand. 
 * \param [in] npts_Q the number of quadrature points to evaluate the line integral
 * \param [in] npts_P the number of quadrature points to evaluate the antiderivative
 * \return the value of the integral
 */
template <typename CurveType,
          typename Lambda,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename CurveType::PointType>>
LambdaRetType evaluate_area_integral(const primal::CurvedPolygon<CurveType>& cpoly,
                                     Lambda&& integrand,
                                     int npts_Q,
                                     int npts_P = 0)
{
  using T = typename CurveType::NumericType;

  static_assert(
    detail::internal::is_integrable_v<T, LambdaRetType>,
    "evaluate_integral methods require addition and scalar multiplication for lambda function "
    "return type");

  if(npts_P <= 0)
  {
    npts_P = npts_Q;
  }

  // Use minimum y-coord of control nodes as lower bound for integration
  T lower_bound_y = cpoly[0][0][1];
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    for(int j = 1; j < cpoly[i].getOrder() + 1; j++)
    {
      lower_bound_y = std::min(lower_bound_y, cpoly[i][j][1]);
    }
  }

  // Evaluate the antiderivative line integral along each component
  LambdaRetType total_integral = LambdaRetType {};
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    total_integral += detail::evaluate_area_integral_component(cpoly[i],
                                                               std::forward<Lambda>(integrand),
                                                               lower_bound_y,
                                                               npts_Q,
                                                               npts_P);
  }

  return total_integral;
}

/*!
 * \brief Evaluate an integral on the interior of a region bound by 2D curves
 *
 * See above definition for details.
 * 
 * \tparam Lambda A callable type taking a CurveType's PointType and returning an integrable type
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the geometry
 * \tparam LambdaRetType A type which supports addition and scalar multiplication
 * \param [in] carray the array of generic curve objects that bound the region
 * \param [in] integrand the lambda function representing the integrand. 
 * \param [in] npts_Q the number of quadrature points to evaluate the line integral
 * \param [in] npts_P the number of quadrature points to evaluate the antiderivative
 * 
 * \note The numerical result is only meaningful if the curves enclose a region
 * 
 * \return the value of the integral
 */
template <typename CurveType,
          typename Lambda,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename CurveType::PointType>>
LambdaRetType evaluate_area_integral(const axom::Array<CurveType>& carray,
                                     Lambda&& integrand,
                                     int npts_Q,
                                     int npts_P = 0)
{
  using T = typename CurveType::NumericType;

  static_assert(
    detail::internal::is_integrable_v<T, LambdaRetType>,
    "evaluate_integral methods require addition and scalar multiplication for lambda function "
    "return type");

  if(npts_P <= 0)
  {
    npts_P = npts_Q;
  }

  if(carray.empty())
  {
    return LambdaRetType {};
  }

  // Use minimum y-coord of control nodes as lower bound for integration
  T lower_bound_y = carray[0][0][1];
  for(int i = 0; i < carray.size(); i++)
  {
    for(int j = 1; j < carray[i].getNumControlPoints(); j++)
    {
      lower_bound_y = std::min(lower_bound_y, carray[i][j][1]);
    }
  }

  // Evaluate the antiderivative line integral along each component
  LambdaRetType total_integral = LambdaRetType {};
  for(int i = 0; i < carray.size(); i++)
  {
    for(const auto& bez : carray[i].extractBezier())
    {
      total_integral += detail::evaluate_area_integral_component(bez,
                                                                 std::forward<Lambda>(integrand),
                                                                 lower_bound_y,
                                                                 npts_Q,
                                                                 npts_P);
    }
  }

  return total_integral;
}
//@}

}  // namespace primal
}  // end namespace axom

#endif