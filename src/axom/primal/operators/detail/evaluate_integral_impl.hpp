// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_EVAL_INTEGRAL_IMPL_HPP_
#define PRIMAL_EVAL_INTEGRAL_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile-time configuration options
#include "axom/primal.hpp"

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#else
  #error "Primal's integral evaluation functions require mfem library."
#endif

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
 * Evaluate the scalar field line integral with MFEM integration rule
 *
 * \param [in] c the Bezier curve object
 * \param [in] integrand the lambda function representing the integrand. 
 * Must accept a 2D point as input and return a double
 * \param [in] quad the mfem integration rule containing nodes and weights
 * \return the value of the integral
 */
template <class Lambda, typename T, int NDIMS>
inline double evaluate_scalar_line_integral_component(const primal::BezierCurve<T, NDIMS>& c,
                                                      Lambda&& scalar_integrand,
                                                      const mfem::IntegrationRule& quad)
{
  // Store/compute quadrature result
  double full_quadrature = 0.0;
  for(int q = 0; q < quad.GetNPoints(); q++)
  {
    // Get intermediate quadrature point
    //  at which to evaluate tangent vector
    auto x_q = c.evaluate(quad.IntPoint(q).x);
    auto dx_q = c.dt(quad.IntPoint(q).x);

    full_quadrature += quad.IntPoint(q).weight * scalar_integrand(x_q) * dx_q.norm();
  }

  return full_quadrature;
}

/*!
 * \brief Evaluate a vector field line integral on a single Bezier curve.
 *
 * Evaluate the vector field line integral with MFEM integration rule
 *
 * \param [in] c the Bezier curve object
 * \param [in] integrand the lambda function representing the integrand. 
 * Must accept a 2D point as input and return a 2D axom vector object
 * \param [in] quad the mfem integration rule containing nodes and weights
 * \return the value of the integral
 */
template <class Lambda, typename T, int NDIMS>
inline double evaluate_vector_line_integral_component(const primal::BezierCurve<T, NDIMS>& c,
                                                      Lambda&& vec_field,
                                                      const mfem::IntegrationRule& quad)
{
  // Store/compute quadrature result
  double full_quadrature = 0.0;
  for(int q = 0; q < quad.GetNPoints(); q++)
  {
    // Get intermediate quadrature point
    //  on which to evaluate dot product
    auto x_q = c.evaluate(quad.IntPoint(q).x);
    auto dx_q = c.dt(quad.IntPoint(q).x);
    auto func_val = vec_field(x_q);

    full_quadrature += quad.IntPoint(q).weight * Vector<T, NDIMS>::dot_product(func_val, dx_q);
  }
  return full_quadrature;
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
 * \param [in] cs the array of Bezier curve objects that bound the region
 * \param [in] integrand the lambda function representing the integrand. 
 * Must accept a 2D point as input and return a double
 * \param [in] The lower bound of integration for the antiderivatives
 * \param [in] quad_Q the quadrature rule for the line integral
 * \param [in] quad_P the quadrature rule for the antiderivative
 * \return the value of the integral, which is mathematically meaningless.
 */
template <class Lambda, typename T>
double evaluate_area_integral_component(const primal::BezierCurve<T, 2>& c,
                                        Lambda&& integrand,
                                        double int_lb,
                                        const mfem::IntegrationRule& quad_Q,
                                        const mfem::IntegrationRule& quad_P)
{
  // Store some intermediate values
  double antiderivative = 0.0;

  // Store/compute quadrature result
  double full_quadrature = 0.0;
  for(int q = 0; q < quad_Q.GetNPoints(); q++)
  {
    // Get intermediate quadrature point
    //  on which to evaluate antiderivative
    auto x_q = c.evaluate(quad_Q.IntPoint(q).x);

    // Evaluate the antiderivative at x_q, add it to full quadrature
    for(int xi = 0; xi < quad_P.GetNPoints(); xi++)
    {
      // Define interior quadrature points
      auto x_qxi = Point2D({x_q[0], (x_q[1] - int_lb) * quad_P.IntPoint(xi).x + int_lb});

      antiderivative = quad_P.IntPoint(xi).weight * (x_q[1] - int_lb) * integrand(x_qxi);

      full_quadrature += quad_Q.IntPoint(q).weight * c.dt(quad_Q.IntPoint(q).x)[0] * -antiderivative;
    }
  }

  return full_quadrature;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif