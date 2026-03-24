// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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

#ifndef PRIMAL_EVAL_INTEGRAL_IMPL_HPP_
#define PRIMAL_EVAL_INTEGRAL_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile-time configuration options
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"
#include "axom/primal/operators/detail/winding_number_2d_memoization.hpp"

#include "axom/core/numerics/quadrature.hpp"

// C++ includes
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <utility>

namespace axom
{
namespace primal
{
namespace detail
{
namespace internal
{
///@{
/// \name Type traits to support integrals of functions with general return types,
///   provided it supports addition and scalar multiplication
template <typename U, typename = void>
struct has_addition : std::false_type
{ };

template <typename U>
struct has_addition<U, std::void_t<decltype(std::declval<U>() + std::declval<U>())>> : std::true_type
{ };

template <typename T, typename U, typename = void>
struct has_scalar_multiplication : std::false_type
{ };

template <typename T, typename U>
struct has_scalar_multiplication<T, U, std::void_t<decltype(std::declval<T>() * std::declval<U>())>>
  : std::true_type
{ };

template <typename T, typename U>
using is_integrable = std::conjunction<has_addition<U>, has_scalar_multiplication<T, U>>;

template <typename T, typename U>
constexpr bool is_integrable_v = is_integrable<T, U>::value;
///@}
}  // namespace internal

///@{
/// \name Evaluates scalar-field line integrals for functions f : R^n -> R^m

/*!
 * \brief Evaluate a line integral on a single Bezier curve.
 *
 * Evaluate the line integral with Gauss-Legendre quadrature
 *
 * \tparam Lambda A callable type taking an CurveType's PointType and returning an integrable type
 * \tparam LambdaRetType A type which supports addition and scalar multiplication
 * \param [in] c the Bezier curve object
 * \param [in] integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda,
          typename T,
          int NDIMS,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierCurve<T, NDIMS>::PointType>>
inline LambdaRetType evaluate_line_integral_component(const BezierCurve<T, NDIMS>& c,
                                                      Lambda&& integrand,
                                                      const int npts)
{
  const axom::numerics::QuadratureRule& quad = axom::numerics::get_gauss_legendre(npts);

  // Store/compute quadrature result
  LambdaRetType full_quadrature = LambdaRetType {};
  for(int q = 0; q < npts; q++)
  {
    // Get intermediate quadrature point
    //  at which to evaluate tangent vector
    auto x_q = c.evaluate(quad.node(q));
    auto dx_q = c.dt(quad.node(q));

    full_quadrature += quad.weight(q) * integrand(x_q) * dx_q.norm();
  }

  return full_quadrature;
}

/*!
 * \brief Evaluate a line integral on a single NURBS curve.
 *
 * Decompose the NURBS curve into Bezier segments, then sum the integral on each
 *
 * \tparam Lambda A callable type taking an CurveType's PointType and returning an integrable type
 * \tparam LambdaRetType The return type of Lambda, which must support addition and scalar multiplication
 * \param [in] n The NURBS curve object
 * \param [in] integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda,
          typename T,
          int NDIMS,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSCurve<T, NDIMS>::PointType>>
inline LambdaRetType evaluate_line_integral_component(const NURBSCurve<T, NDIMS>& n,
                                                      Lambda&& integrand,
                                                      const int npts)
{
  LambdaRetType total_integral = LambdaRetType {};
  for(const auto& bez : n.extractBezier())
  {
    total_integral +=
      detail::evaluate_line_integral_component(bez, std::forward<Lambda>(integrand), npts);
  }
  return total_integral;
}

/*!
 * \brief Evaluate an integral on a single NURBS curve with cached data for GWN evaluation.
 *
 * The cache object has already decomposed the NURBS curve into Bezier segments, 
 *  which are used to evaluate the integral over each
 *
 * \tparam Lambda A callable type taking an CurveType's PointType and returning an integrable type
 * \tparam LambdaRetType A type which supports addition and scalar multiplication
 * \param [in] n The NURBS curve object
 * \param [in] integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSCurveGWNCache<T>::PointType>>
inline LambdaRetType evaluate_line_integral_component(const NURBSCurveGWNCache<T>& nc,
                                                      Lambda&& integrand,
                                                      const int npts)
{
  LambdaRetType total_integral = LambdaRetType {};
  for(int i = 0; i < nc.getNumKnotSpans(); ++i)
  {
    // Assuming the cache is properly initialized, this operation will never add to the cache
    const auto& this_bezier_data = nc.getSubdivisionData(i, 0, 0);
    total_integral += detail::evaluate_line_integral_component(this_bezier_data.getCurve(),
                                                               std::forward<Lambda>(integrand),
                                                               npts);
  }

  return total_integral;
}
//@}

///@{
/// \name Evaluates vector-field line integrals for functions f : R^n -> R^n

/*!
 * \brief Evaluate a vector field line integral on a single Bezier curve.
 *
 * Evaluate the vector field line integral with Gauss-Legendre quadrature
 *
 * \tparam Lambda A callable type taking an CurveType's PointType and returning its numeric type
 * \param [in] c the Bezier curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
inline T evaluate_vector_line_integral_component(const primal::BezierCurve<T, NDIMS>& c,
                                                 Lambda&& vector_integrand,
                                                 const int npts)
{
  const axom::numerics::QuadratureRule& quad = axom::numerics::get_gauss_legendre(npts);

  // Store/compute quadrature result
  T full_quadrature = T {};
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
 * \tparam Lambda A callable type taking an CurveType's PointType and returning its numeric type
 * \param [in] n The NURBS curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda, typename T, int NDIMS>
inline T evaluate_vector_line_integral_component(const primal::NURBSCurve<T, NDIMS>& n,
                                                 Lambda&& vector_integrand,
                                                 const int npts)
{
  T total_integral = T {};
  for(const auto& bez : n.extractBezier())
  {
    total_integral +=
      detail::evaluate_vector_line_integral_component(bez,
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
 * \tparam Lambda A callable type taking an CurveType's PointType and returning its numeric type
 * \param [in] n The NURBS curve object
 * \param [in] vector_integrand the lambda function representing the integrand. 
 * \param [in] npts The number of quadrature points in the rule
 * \return the value of the integral
 */
template <typename Lambda, typename T>
inline T evaluate_vector_line_integral_component(const NURBSCurveGWNCache<T>& nc,
                                                 Lambda&& vector_integrand,
                                                 const int npts)
{
  T total_integral = T {};
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
//@}

///@{
/// \name Evaluates scalar-field 2D area integrals for functions f : R^2 -> R^m

/*!
 * \brief Evaluate the area integral across one component of the curved polygon.
 *
 * Intended to be called for each BezierCurve object in a curved polygon.
 * Uses a Spectral Mesh-Free Quadrature derived from Green's theorem, evaluating
 * the area integral as a line integral of the antiderivative over the curve.
 * For algorithm details, see "Spectral Mesh-Free Quadrature for Planar 
 * Regions Bounded by Rational Parametric Curves" by David Gunderman et al.
 * 
 * \tparam Lambda A callable type taking an CurveType's PointType and returning an integrable type
 * \tparam LambdaRetType A type which supports addition and scalar multiplication
 * \param [in] c The component Bezier curve 
 * \param [in] integrand The lambda function representing the scalar integrand. 
 * \param [in] The lower bound of integration for the antiderivatives
 * \param [in] npts_Q The number of quadrature points for the line integral
 * \param [in] npts_P The number of quadrature points for the antiderivative
 * \return the value of the integral, which is mathematically meaningless.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierCurve<T, 2>::PointType>>
LambdaRetType evaluate_area_integral_component(const primal::BezierCurve<T, 2>& c,
                                               Lambda&& integrand,
                                               double int_lb,
                                               const int npts_Q,
                                               const int npts_P)
{
  const axom::numerics::QuadratureRule& quad_Q = axom::numerics::get_gauss_legendre(npts_Q);
  const axom::numerics::QuadratureRule& quad_P = axom::numerics::get_gauss_legendre(npts_P);

  // Store/compute quadrature result
  LambdaRetType full_quadrature = LambdaRetType {};
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

      auto antiderivative = quad_P.weight(xi) * (x_q[1] - int_lb) * integrand(x_qxi);

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
 * \tparam Lambda A callable type taking an CurveType's PointType and returning an integrable type
 * \tparam LambdaRetType A type which supports addition and scalar multiplication
 * \param [in] n The component NURBSCurve object
 * \param [in] integrand The lambda function representing the scalar integrand. 
 * \param [in] int_lb The lower bound of integration for the antiderivatives
 * \param [in] npts_Q The number of quadrature points for the line integral
 * \param [in] npts_P The number of quadrature points for the antiderivative
 * \return the value of the integral, which is mathematically meaningless.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSCurve<T, 2>::PointType>>
LambdaRetType evaluate_area_integral_component(const primal::NURBSCurve<T, 2>& n,
                                               Lambda&& integrand,
                                               double int_lb,
                                               const int npts_Q,
                                               const int npts_P)
{
  LambdaRetType total_integral = LambdaRetType {};
  for(const auto& bez : n.extractBezier())
  {
    total_integral += detail::evaluate_area_integral_component(bez,
                                                               std::forward<Lambda>(integrand),
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
 * \tparam Lambda A callable type taking an CurveType's PointType and returning an integrable type
 * \tparam RetType The return type of Lambda, which must support addition and scalar multiplication
 * \param [in] nc The component NURBSCurveGWNCache object
 * \param [in] integrand The lambda function representing the scalar integrand. 
 * \param [in] int_lb The lower bound of integration for the antiderivatives
 * \param [in] npts_Q The number of quadrature points for the line integral
 * \param [in] npts_P The number of quadrature points for the antiderivative
 * \return the value of the integral, which is mathematically meaningless.
 */
template <typename Lambda,
          typename T,
          typename RetType = std::invoke_result_t<Lambda, typename NURBSCurveGWNCache<T>::PointType>>
inline RetType evaluate_area_integral_component(const NURBSCurveGWNCache<T>& nc,
                                                Lambda&& integrand,
                                                double int_lb,
                                                const int npts_Q,
                                                const int npts_P)
{
  RetType total_integral = RetType {};
  for(int i = 0; i < nc.getNumKnotSpans(); ++i)
  {
    // Assuming the cache is properly initialized, this operation will never add to the cache
    const auto& this_bezier_data = nc.getSubdivisionData(i, 0, 0);
    total_integral += detail::evaluate_area_integral_component(this_bezier_data.getCurve(),
                                                               std::forward<Lambda>(integrand),
                                                               int_lb,
                                                               npts_Q,
                                                               npts_P);
  }

  return total_integral;
}
///@{
/// \name Helper routines for patch-based integrals in 3D

template <typename CurveType>
inline typename CurveType::NumericType curve_array_lower_bound_y(const axom::Array<CurveType>& carray)
{
  using T = typename CurveType::NumericType;

  SLIC_ASSERT(!carray.empty() && carray[0].getNumControlPoints() > 0);

  T lower_bound_y = carray[0][0][1];
  for(int i = 0; i < carray.size(); ++i)
  {
    for(int j = 0; j < carray[i].getNumControlPoints(); ++j)
    {
      lower_bound_y = axom::utilities::min(lower_bound_y, carray[i][j][1]);
    }
  }

  return lower_bound_y;
}

template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierPatch<T, 3>::PointType>>
inline LambdaRetType evaluate_surface_integral_component(const primal::BezierPatch<T, 3>& b,
                                                         Lambda&& integrand,
                                                         const int npts)
{
  const axom::numerics::QuadratureRule& quad = axom::numerics::get_gauss_legendre(npts);

  LambdaRetType full_quadrature = LambdaRetType {};
  for(int qu = 0; qu < npts; ++qu)
  {
    for(int qv = 0; qv < npts; ++qv)
    {
      const auto x_q = b.evaluate(quad.node(qu), quad.node(qv));
      const auto n_q = b.normal(quad.node(qu), quad.node(qv));

      full_quadrature += quad.weight(qu) * quad.weight(qv) * integrand(x_q) * n_q.norm();
    }
  }

  return full_quadrature;
}

template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSPatch<T, 3>::PointType>>
inline LambdaRetType evaluate_surface_integral_component(const primal::NURBSPatch<T, 3>& n,
                                                         Lambda&& integrand,
                                                         const int npts_Q,
                                                         const int npts_P)
{
  if(!n.isTrimmed())
  {
    LambdaRetType total_integral = LambdaRetType {};
    for(const auto& bez : n.extractBezier())
    {
      total_integral += detail::evaluate_surface_integral_component(bez, integrand, npts_Q);
    }

    return total_integral;
  }

  LambdaRetType total_integral = LambdaRetType {};
  for(const auto& split_patch : n.extractTrimmedBezier())
  {
    const auto& curves = split_patch.getTrimmingCurves();
    if(curves.empty())
    {
      continue;
    }

    // clamp the lower bound to lie within the parametric space of the patch
    const auto lower_bound_y = axom::utilities::clampVal(detail::curve_array_lower_bound_y(curves),
                                                         split_patch.getMinKnot_v(),
                                                         split_patch.getMaxKnot_v());
    for(int i = 0; i < curves.size(); ++i)
    {
      total_integral += detail::evaluate_area_integral_component(
        curves[i],
        [&split_patch, &integrand](Point<T, 2> uv) -> LambdaRetType {
          const auto x_q = split_patch.evaluate(uv[0], uv[1]);
          const auto n_q = split_patch.normal(uv[0], uv[1]);
          return integrand(x_q) * n_q.norm();
        },
        lower_bound_y,
        npts_Q,
        npts_P);
    }
  }

  return total_integral;
}

template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierPatch<T, 3>::PointType>>
inline LambdaRetType evaluate_volume_integral_component(const primal::BezierPatch<T, 3>& b,
                                                        Lambda&& integrand,
                                                        double int_lb,
                                                        const int npts_uv,
                                                        const int npts_z)
{
  const axom::numerics::QuadratureRule& quad_uv = axom::numerics::get_gauss_legendre(npts_uv);
  const axom::numerics::QuadratureRule& quad_z = axom::numerics::get_gauss_legendre(npts_z);

  LambdaRetType full_quadrature = LambdaRetType {};
  for(int qu = 0; qu < npts_uv; ++qu)
  {
    for(int qv = 0; qv < npts_uv; ++qv)
    {
      const auto x_q = b.evaluate(quad_uv.node(qu), quad_uv.node(qv));
      const auto n_q = b.normal(quad_uv.node(qu), quad_uv.node(qv));

      LambdaRetType antiderivative = LambdaRetType {};
      const T z_scale = x_q[2] - int_lb;
      for(int qz = 0; qz < npts_z; ++qz)
      {
        const auto x_qz = Point<T, 3>({x_q[0], x_q[1], z_scale * quad_z.node(qz) + int_lb});
        antiderivative += quad_z.weight(qz) * z_scale * integrand(x_qz);
      }

      full_quadrature += quad_uv.weight(qu) * quad_uv.weight(qv) * antiderivative * n_q[2];
    }
  }

  return full_quadrature;
}

template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSPatch<T, 3>::PointType>>
inline LambdaRetType evaluate_volume_integral_component(const primal::NURBSPatch<T, 3>& n,
                                                        Lambda&& integrand,
                                                        double int_lb,
                                                        const int npts_Q,
                                                        const int npts_P,
                                                        const int npts_Z)
{
  if(!n.isTrimmed())
  {
    LambdaRetType total_integral = LambdaRetType {};
    for(const auto& bez : n.extractBezier())
    {
      total_integral +=
        detail::evaluate_volume_integral_component(bez, integrand, int_lb, npts_Q, npts_Z);
    }

    return total_integral;
  }

  const axom::numerics::QuadratureRule& quad_z = axom::numerics::get_gauss_legendre(npts_Z);

  LambdaRetType total_integral = LambdaRetType {};
  for(const auto& split_patch : n.extractTrimmedBezier())
  {
    const auto& curves = split_patch.getTrimmingCurves();
    if(curves.empty())
    {
      continue;
    }

    // clamp the lower bound to lie within the parametric space of the patch
    const auto lower_bound_y = axom::utilities::clampVal(detail::curve_array_lower_bound_y(curves),
                                                         split_patch.getMinKnot_v(),
                                                         split_patch.getMaxKnot_v());
    for(int i = 0; i < curves.size(); ++i)
    {
      total_integral += detail::evaluate_area_integral_component(
        curves[i],
        [&split_patch, &integrand, &int_lb, &quad_z, &npts_Z](Point<T, 2> uv) -> LambdaRetType {
          const auto x_q = split_patch.evaluate(uv[0], uv[1]);
          const auto n_q = split_patch.normal(uv[0], uv[1]);

          LambdaRetType antiderivative = LambdaRetType {};
          const T z_scale = x_q[2] - int_lb;
          for(int qz = 0; qz < npts_Z; ++qz)
          {
            const auto x_qz = Point<T, 3>({x_q[0], x_q[1], z_scale * quad_z.node(qz) + int_lb});
            antiderivative += quad_z.weight(qz) * z_scale * integrand(x_qz);
          }

          return antiderivative * n_q[2];
        },
        lower_bound_y,
        npts_Q,
        npts_P);
    }
  }

  return total_integral;
}

//@}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
