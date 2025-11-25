// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_EVAL_INTEGRAL_IMPL_HPP_
#define PRIMAL_EVAL_INTEGRAL_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile-time configuration options
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/operators/detail/winding_number_2d_memoization.hpp"

#include "axom/core/numerics/quadrature.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
namespace detail
{

/*!
 * \brief Evaluate a scalar field line integral on a single Bezier curve.
 *
 * Evaluate the scalar field line integral with Gauss-Legendre quadrature
 *
 * \tparam Lambda A callable type taking an axom::Point type returning a numeric type
 * \param [in] c the Bezier curve object
 * \param [in] scalar_integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
inline double evaluate_scalar_line_integral_component(const primal::BezierCurve<T, NDIMS>& c,
                                                      Lambda&& scalar_integrand,
                                                      const int npts)
{
  const axom::numerics::QuadratureRule& quad = axom::numerics::get_gauss_legendre(npts);

  // Store/compute quadrature result
  double full_quadrature = 0.0;
  for(int q = 0; q < npts; q++)
  {
    // Get intermediate quadrature point
    //  at which to evaluate tangent vector
    auto x_q = c.evaluate(quad.node(q));
    auto dx_q = c.dt(quad.node(q));

    full_quadrature += quad.weight(q) * scalar_integrand(x_q) * dx_q.norm();
  }

  return full_quadrature;
}

/*!
 * \brief Evaluate a scalar field line integral on a single NURBS curve.
 *
 * Decompose the NURBS curve into Bezier segments, then sum the integral on each
 *
 * \tparam Lambda A callable type taking an axom::Point type returning a numeric type
 * \param [in] n The NURBS curve object
 * \param [in] scalar_integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
inline double evaluate_scalar_line_integral_component(const primal::NURBSCurve<T, NDIMS>& n,
                                                      Lambda&& scalar_integrand,
                                                      const int npts)
{
  const auto beziers = n.extractBezier();
  double total_integral = 0.0;
  for(int i = 0; i < beziers.size(); ++i)
  {
    total_integral +=
      detail::evaluate_scalar_line_integral_component(beziers[i],
                                                      std::forward<Lambda>(scalar_integrand),
                                                      npts);
  }
  return total_integral;
}

/*!
 * \brief Evaluate the scalar integral on a single NURBS curve with cached data for GWN evaluation.
 *
 * The cache object has already decomposed the NURBS curve into Bezier segments, 
 *  which are used to evaluate the integral over each
 *
 * \tparam Lambda A callable type taking an axom::Point type returning a numeric type
 * \param [in] n The NURBS curve object
 * \param [in] scalar_integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
inline double evaluate_scalar_line_integral_component(const primal::detail::NURBSCurveGWNCache<T>& nc,
                                                      Lambda&& scalar_integrand,
                                                      const int npts)
{
  double total_integral = 0.0;
  for(int i = 0; i < nc.getNumKnotSpans(); ++i)
  {
    // Assuming the cache is properly initialized, this operation will never add to the cache
    const auto& this_bezier_data = nc.getSubdivisionData(i, 0, 0);
    total_integral +=
      detail::evaluate_scalar_line_integral_component(this_bezier_data.getCurve(),
                                                      std::forward<Lambda>(scalar_integrand),
                                                      npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a vector field line integral on a single Bezier curve.
 *
 * Evaluate the scalar field line integral with Gauss-Legendre quadrature
 *
 * \tparam Lambda A callable type taking an axom::Point type returning an axom::Vector type
 * \param [in] c the Bezier curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
inline double evaluate_vector_line_integral_component(const primal::BezierCurve<T, NDIMS>& c,
                                                      Lambda&& vector_integrand,
                                                      const int npts)
{
  const axom::numerics::QuadratureRule& quad = axom::numerics::get_gauss_legendre(npts);

  // Store/compute quadrature result
  double full_quadrature = 0.0;
  for(int q = 0; q < npts; q++)
  {
    // Get intermediate quadrature point
    //  on which to evaluate dot product
    auto x_q = c.evaluate(quad.node(q));
    auto dx_q = c.dt(quad.node(q));
    auto func_val = vector_integrand(x_q);

    full_quadrature += quad.weight(q) * Vector<T, NDIMS>::dot_product(func_val, dx_q);
  }
  return full_quadrature;
}

/*!
 * \brief Evaluate a vector field line integral on a single NURBS curve.
 *
 * Decompose the NURBS curve into Bezier segments, then sum the integral on each
 *
 * \tparam Lambda A callable type taking an axom::Point type returning an axom::Vector type
 * \param [in] n The NURBS curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
inline double evaluate_vector_line_integral_component(const primal::NURBSCurve<T, NDIMS>& n,
                                                      Lambda&& vector_integrand,
                                                      const int npts)
{
  const auto beziers = n.extractBezier();
  double total_integral = 0.0;
  for(int i = 0; i < beziers.size(); ++i)
  {
    total_integral +=
      detail::evaluate_vector_line_integral_component(beziers[i],
                                                      std::forward<Lambda>(vector_integrand),
                                                      npts);
  }
  return total_integral;
}

/*!
 * \brief Evaluate the vector integral on a single NURBS curve with cached data for GWN evaluation.
 *
 * The cache object has already decomposed the NURBS curve into Bezier segments, 
 *  which are used to evaluate the integral over each
 *
 * \tparam Lambda A callable type taking an axom::Point type returning an axom::Vector type
 * \param [in] n The NURBS curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
inline double evaluate_vector_line_integral_component(const primal::detail::NURBSCurveGWNCache<T>& nc,
                                                      Lambda&& vector_integrand,
                                                      const int npts)
{
  double total_integral = 0.0;
  for(int i = 0; i < nc.getNumKnotSpans(); ++i)
  {
    // Assuming the cache is properly initialized, this operation will never add to the cache
    const auto& this_bezier_data = nc.getSubdivisionData(i, 0, 0);
    total_integral +=
      detail::evaluate_vector_line_integral_component(this_bezier_data.getCurve(),
                                                      std::forward<Lambda>(vector_integrand),
                                                      npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate the area integral across one component of the curved polygon.
 *
 * Intended to be called for each BezierCurve object in a curved polygon.
 * Uses a Spectral Mesh-Free Quadrature derived from Green's theorem, evaluating
 * the area integral as a line integral of the antiderivative over the curve.
 * For algorithm details, see "Spectral Mesh-Free Quadrature for Planar 
 * Regions Bounded by Rational Parametric Curves" by David Gunderman et al.
 * 
 * \tparam Lambda A callable type taking a 2D axom::Point type returning a numeric type
 * \param [in] c The component Bezier curve 
 * \param [in] scalar_integrand The lambda function representing the scalar integrand. 
 * \param [in] The lower bound of integration for the antiderivatives
 * \param [in] npts_Q The number of quadrature points for the line integral
 * \param [in] npts_P The number of quadrature points for the antiderivative
 * \return the value of the integral, which is mathematically meaningless.
 */
template <typename Lambda, typename T>
double evaluate_area_integral_component(const primal::BezierCurve<T, 2>& c,
                                        Lambda&& scalar_integrand,
                                        double int_lb,
                                        const int npts_Q,
                                        const int npts_P)
{
  const axom::numerics::QuadratureRule& quad_Q = axom::numerics::get_gauss_legendre(npts_Q);
  const axom::numerics::QuadratureRule& quad_P = axom::numerics::get_gauss_legendre(npts_P);

  // Store some intermediate values
  double antiderivative = 0.0;

  // Store/compute quadrature result
  double full_quadrature = 0.0;
  for(int q = 0; q < npts_Q; q++)
  {
    // Get intermediate quadrature point
    //  on which to evaluate antiderivative
    auto x_q = c.evaluate(quad_Q.node(q));

    // Evaluate the antiderivative at x_q, add it to full quadrature
    for(int xi = 0; xi < npts_P; xi++)
    {
      // Define interior quadrature points
      auto x_qxi = Point<T, 2>({x_q[0], (x_q[1] - int_lb) * quad_P.node(xi) + int_lb});

      antiderivative = quad_P.weight(xi) * (x_q[1] - int_lb) * scalar_integrand(x_qxi);

      full_quadrature += quad_Q.weight(q) * c.dt(quad_Q.node(q))[0] * -antiderivative;
    }
  }

  return full_quadrature;
}

/*!
 * \brief Evaluate the area integral across one component NURBS of the region.
 *
 * Intended to be called for each NURBSCurve object bounding a closed 2D region.
 *
 * \tparam Lambda A callable type taking a 2D axom::Point type returning a numeric type
 * \param [in] n The component NURBSCurve object
 * \param [in] scalar_integrand The lambda function representing the scalar integrand. 
 * \param [in] The lower bound of integration for the antiderivatives
 * \param [in] npts_Q The number of quadrature points for the line integral
 * \param [in] npts_P The number of quadrature points for the antiderivative
 * \return the value of the integral, which is mathematically meaningless.
 */
template <typename Lambda, typename T>
double evaluate_area_integral_component(const primal::NURBSCurve<T, 2>& n,
                                        Lambda&& scalar_integrand,
                                        double int_lb,
                                        const int npts_Q,
                                        const int npts_P)
{
  auto beziers = n.extractBezier();
  double total_integral = 0.0;
  for(int i = 0; i < beziers.size(); ++i)
  {
    total_integral += detail::evaluate_area_integral_component(beziers[i],
                                                               std::forward<Lambda>(scalar_integrand),
                                                               int_lb,
                                                               npts_Q,
                                                               npts_P);
  }
  return total_integral;
}

/*!
 * \brief Evaluate the area integral across one component NURBSCurveGWNCache of the region.
 *
 * Intended to be called for each NURBSCurveGWNCache object bounding a closed 2D region.
 *
 * \tparam Lambda A callable type taking a 2D axom::Point type returning a numeric type
 * \param [in] nc The component NURBSCurveGWNCache object
 * \param [in] scalar_integrand The lambda function representing the scalar integrand. 
 * \param [in] The lower bound of integration for the antiderivatives
 * \param [in] npts_Q The number of quadrature points for the line integral
 * \param [in] npts_P The number of quadrature points for the antiderivative
 * \return the value of the integral, which is mathematically meaningless.
 */
template <typename Lambda, typename T, int NDIMS>
inline double evaluate_area_integral_component(const primal::detail::NURBSCurveGWNCache<T>& nc,
                                               Lambda&& scalar_integrand,
                                               double int_lb,
                                               const int npts_Q,
                                               const int npts_P)
{
  double total_integral = 0.0;
  for(int i = 0; i < nc.getNumKnotSpans(); ++i)
  {
    // Assuming the cache is properly initialized, this operation will never add to the cache
    const auto& this_bezier_data = nc.getSubdivisionData(i, 0, 0);
    total_integral += detail::evaluate_area_integral_component(this_bezier_data.getCurve(),
                                                               std::forward<Lambda>(scalar_integrand),
                                                               int_lb,
                                                               npts_Q,
                                                               npts_P);
  }

  return total_integral;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif