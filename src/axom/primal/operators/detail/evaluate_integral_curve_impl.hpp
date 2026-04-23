// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file evaluate_integral_curve_impl.hpp
 *
 * \brief Implementation helpers for curve/area integral evaluation.
 *
 * \note This header intentionally avoids including surface/volume (patch)
 * dependencies to prevent circular include chains (e.g. with NURBSPatch).
 */

#ifndef PRIMAL_EVAL_INTEGRAL_CURVE_IMPL_HPP_
#define PRIMAL_EVAL_INTEGRAL_CURVE_IMPL_HPP_

// Axom includes
#include "axom/core.hpp"
#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/operators/detail/winding_number_2d_memoization.hpp"

#include "axom/core/numerics/quadrature.hpp"

// C++ includes
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
}  // namespace internal

///@{
/// \name Scalar-field line integrals

template <typename Lambda,
          typename T,
          int NDIMS,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierCurve<T, NDIMS>::PointType>>
inline LambdaRetType evaluate_line_integral_component(const BezierCurve<T, NDIMS>& c,
                                                      Lambda&& integrand,
                                                      const int npts)
{
  const axom::numerics::QuadratureRule& quad = axom::numerics::get_gauss_legendre(npts);

  LambdaRetType full_quadrature = LambdaRetType {};
  for(int q = 0; q < npts; q++)
  {
    auto x_q = c.evaluate(quad.node(q));
    auto dx_q = c.dt(quad.node(q));

    full_quadrature += quad.weight(q) * integrand(x_q) * dx_q.norm();
  }

  return full_quadrature;
}

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
    const auto& this_bezier_data = nc.getSubdivisionData(i, 0, 0);
    total_integral += detail::evaluate_line_integral_component(this_bezier_data.getCurve(),
                                                               std::forward<Lambda>(integrand),
                                                               npts);
  }

  return total_integral;
}
//@}

///@{
/// \name Vector-field line integrals

template <typename Lambda, typename T, int NDIMS>
inline T evaluate_vector_line_integral_component(const primal::BezierCurve<T, NDIMS>& c,
                                                 Lambda&& vector_integrand,
                                                 const int npts)
{
  const axom::numerics::QuadratureRule& quad = axom::numerics::get_gauss_legendre(npts);

  T full_quadrature = T {};
  for(int q = 0; q < npts; q++)
  {
    auto x_q = c.evaluate(quad.node(q));
    auto dx_q = c.dt(quad.node(q));
    auto func_val = vector_integrand(x_q);

    full_quadrature += quad.weight(q) * Vector<T, NDIMS>::dot_product(func_val, dx_q);
  }

  return full_quadrature;
}

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

template <typename Lambda, typename T>
inline T evaluate_vector_line_integral_component(const NURBSCurveGWNCache<T>& nc,
                                                 Lambda&& vector_integrand,
                                                 const int npts)
{
  T total_integral = T {};
  for(int i = 0; i < nc.getNumKnotSpans(); ++i)
  {
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
/// \name Scalar-field 2D area integrals

template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierCurve<T, 2>::PointType>>
inline LambdaRetType evaluate_area_integral_component(const primal::BezierCurve<T, 2>& c,
                                                      Lambda&& integrand,
                                                      double int_lb,
                                                      const int npts_Q,
                                                      const int npts_P)
{
  const axom::numerics::QuadratureRule& quad_Q = axom::numerics::get_gauss_legendre(npts_Q);
  const axom::numerics::QuadratureRule& quad_P = axom::numerics::get_gauss_legendre(npts_P);

  LambdaRetType full_quadrature = LambdaRetType {};
  for(int q = 0; q < npts_Q; q++)
  {
    auto x_q = c.evaluate(quad_Q.node(q));

    for(int xi = 0; xi < npts_P; xi++)
    {
      auto x_qxi = Point<T, 2>({x_q[0], (x_q[1] - int_lb) * quad_P.node(xi) + int_lb});

      auto antiderivative = quad_P.weight(xi) * (x_q[1] - int_lb) * integrand(x_qxi);

      full_quadrature += quad_Q.weight(q) * c.dt(quad_Q.node(q))[0] * -antiderivative;
    }
  }

  return full_quadrature;
}

template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSCurve<T, 2>::PointType>>
inline LambdaRetType evaluate_area_integral_component(const primal::NURBSCurve<T, 2>& n,
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
    const auto& this_bezier_data = nc.getSubdivisionData(i, 0, 0);
    total_integral += detail::evaluate_area_integral_component(this_bezier_data.getCurve(),
                                                               std::forward<Lambda>(integrand),
                                                               int_lb,
                                                               npts_Q,
                                                               npts_P);
  }

  return total_integral;
}
//@}

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

//@}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
