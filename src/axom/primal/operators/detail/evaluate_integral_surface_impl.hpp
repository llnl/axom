// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file evaluate_integral_surface_impl.hpp
 *
 * \brief Implementation helpers for surface/volume integral evaluation.
 */

#ifndef PRIMAL_EVAL_INTEGRAL_SURFACE_IMPL_HPP_
#define PRIMAL_EVAL_INTEGRAL_SURFACE_IMPL_HPP_

// Axom includes
#include "axom/core.hpp"
#include "axom/config.hpp"

#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"
#include "axom/primal/operators/detail/evaluate_integral_curve_impl.hpp"

// C++ includes
#include <type_traits>
#include <utility>

namespace axom
{
namespace primal
{
namespace detail
{
///@{
/// \name Scalar-field surface/volume integral components

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
