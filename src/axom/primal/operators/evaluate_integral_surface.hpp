// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file evaluate_integral_surface.hpp
 *
 * \brief Consists of methods that evaluate surface and volume integrals on 
 * surface patches and regions defined by collections of patches.
 *
 * All integrals are evaluated numerically with Gauss-Legendre quadrature
 * 
 * 3D integrals computed with "High-accuracy mesh-free quadrature 
 * for trimmed parametric surfaces and volumes" by D. Gunderman et al.
 * https://doi.org/10.1016/j.cad.2021.103093
 */

#ifndef PRIMAL_EVAL_INTEGRAL_SURFACE_HPP_
#define PRIMAL_EVAL_INTEGRAL_SURFACE_HPP_

// Axom includes
#include "axom/core.hpp"
#include "axom/config.hpp"

#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"
#include "axom/primal/operators/detail/evaluate_integral_surface_impl.hpp"

namespace axom
{
namespace primal
{
///@{
/// \name Evaluates scalar-field surface integrals for functions f : R^3 -> R^m

/*!
 * \brief Evaluate a scalar surface integral on a single Bezier patch.
 *
 * Uses tensor-product Gauss-Legendre quadrature in the patch parameter space.
 *
 * \param [in] patch the Bezier patch
 * \param [in] integrand callable representing the integrand
 * \param [in] npts the number of quadrature points in each parametric direction
 *
 * \pre The patch parameterization must be valid on its full parameter domain.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierPatch<T, 3>::PointType>>
LambdaRetType evaluate_surface_integral(const primal::BezierPatch<T, 3>& patch,
                                        Lambda&& integrand,
                                        int npts)
{
  static_assert(detail::internal::is_integrable_v<T, LambdaRetType>,
                "evaluate_integral methods require addition and scalar multiplication for lambda "
                "function return type");

  return detail::evaluate_surface_integral_component(patch, std::forward<Lambda>(integrand), npts);
}

/*!
 * \brief Evaluate a scalar surface integral on a single NURBS patch.
 *
 * Untrimmed patches are integrated by Bezier extraction followed by tensor-product
 * Gauss-Legendre quadrature. Trimmed patches are integrated by reducing the
 * parameter-space area integral to line integrals over the trimming curves.
 *
 * \param [in] patch the NURBS patch
 * \param [in] integrand callable representing the integrand
 * \param [in] npts_Q the number of quadrature points on each trimming curve or
 *                    in each parametric direction for untrimmed Bezier pieces
 * \param [in] npts_P the number of quadrature points used for numerical
 *                    antidifferentiation in parameter space
 *
 * \pre The patch parameterization must be valid on its full parameter domain.
 * \pre If the patch is trimmed, its trimming curves must bound the intended
 *      interior region in parameter space.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSPatch<T, 3>::PointType>>
LambdaRetType evaluate_surface_integral(const primal::NURBSPatch<T, 3>& patch,
                                        Lambda&& integrand,
                                        int npts_Q,
                                        int npts_P = 0)
{
  static_assert(detail::internal::is_integrable_v<T, LambdaRetType>,
                "evaluate_integral methods require addition and scalar multiplication for lambda "
                "function return type");

  if(npts_P <= 0)
  {
    npts_P = npts_Q;
  }

  return detail::evaluate_surface_integral_component(patch,
                                                     std::forward<Lambda>(integrand),
                                                     npts_Q,
                                                     npts_P);
}

/*!
 * \brief Evaluate a scalar surface integral on a collection of Bezier patches.
 *
 * The result is the sum of the surface integrals over each patch in the array.
 *
 * \param [in] patches the patch collection
 * \param [in] integrand callable representing the integrand
 * \param [in] npts the number of quadrature points in each parametric direction
 *
 * \pre Each patch parameterization must be valid on its full parameter domain.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierPatch<T, 3>::PointType>>
LambdaRetType evaluate_surface_integral(const axom::Array<BezierPatch<T, 3>>& patches,
                                        Lambda&& integrand,
                                        int npts)
{
  static_assert(detail::internal::is_integrable_v<T, LambdaRetType>,
                "evaluate_integral methods require addition and scalar multiplication for lambda "
                "function return type");

  LambdaRetType total_integral = LambdaRetType {};
  for(int i = 0; i < patches.size(); ++i)
  {
    total_integral += detail::evaluate_surface_integral_component(patches[i], integrand, npts);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a scalar surface integral on a collection of NURBS patches.
 *
 * The result is the sum of the surface integrals over each patch in the array.
 *
 * \param [in] patches the patch collection
 * \param [in] integrand callable representing the integrand
 * \param [in] npts_Q the number of quadrature points on each trimming curve or
 *                    in each parametric direction for untrimmed Bezier pieces
 * \param [in] npts_P the number of quadrature points used for numerical
 *                    antidifferentiation in parameter space
 *
 * \pre Each patch parameterization must be valid on its full parameter domain.
 * \pre Any trimmed patch in the array must have trimming curves that bound its
 *      intended interior region in parameter space.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSPatch<T, 3>::PointType>>
LambdaRetType evaluate_surface_integral(const axom::Array<NURBSPatch<T, 3>>& patches,
                                        Lambda&& integrand,
                                        int npts_Q,
                                        int npts_P = 0)
{
  static_assert(detail::internal::is_integrable_v<T, LambdaRetType>,
                "evaluate_integral methods require addition and scalar multiplication for lambda "
                "function return type");

  if(npts_P <= 0)
  {
    npts_P = npts_Q;
  }

  LambdaRetType total_integral = LambdaRetType {};
  for(int i = 0; i < patches.size(); ++i)
  {
    total_integral +=
      detail::evaluate_surface_integral_component(patches[i], integrand, npts_Q, npts_P);
  }

  return total_integral;
}
//@}

///@{
/// \name Evaluates scalar-field volume integrals for functions f : R^3 -> R^m

/*!
 * \brief Evaluate a scalar volume-integral contribution from a single Bezier patch.
 *
 * This applies the Stokes-based reduction used for the full volume algorithm to
 * one patch using a z-directed numerical antiderivative.
 *
 * \param [in] patch the Bezier patch
 * \param [in] integrand callable representing the integrand
 * \param [in] lower_bound_z the shared lower integration bound used for the
 *             z-directed antiderivative across the full boundary
 * \param [in] npts_uv the number of quadrature points in each patch parameter direction
 * \param [in] npts_z the number of quadrature points used for numerical
 *                    antidifferentiation in z
 *
 * \pre The patch parameterization must be valid on its full parameter domain.
 * \pre The returned value is geometrically meaningful as a volume contribution only
 *      when this patch is interpreted as part of a closed, consistently oriented
 *      boundary that uses the same lower integration bound.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierPatch<T, 3>::PointType>>
LambdaRetType evaluate_volume_integral(const primal::BezierPatch<T, 3>& patch,
                                       Lambda&& integrand,
                                       T lower_bound_z,
                                       int npts_uv,
                                       int npts_z = 0)
{
  static_assert(detail::internal::is_integrable_v<T, LambdaRetType>,
                "evaluate_integral methods require addition and scalar multiplication for lambda "
                "function return type");

  if(npts_z <= 0)
  {
    npts_z = npts_uv;
  }

  return detail::evaluate_volume_integral_component(patch,
                                                    std::forward<Lambda>(integrand),
                                                    lower_bound_z,
                                                    npts_uv,
                                                    npts_z);
}

/*!
 * \brief Evaluate a scalar volume-integral contribution from a single NURBS patch.
 *
 * Trimmed patches use the same Green/Stokes reduction as the surface-integral
 * algorithm, combined with a z-directed numerical antiderivative for the volume
 * reduction.
 *
 * \param [in] patch the NURBS patch
 * \param [in] integrand callable representing the integrand
 * \param [in] lower_bound_z the shared lower integration bound used for the
 *             z-directed antiderivative across the full boundary
 * \param [in] npts_Q the number of quadrature points on each trimming curve or
 *                    in each parametric direction for untrimmed Bezier pieces
 * \param [in] npts_P the number of quadrature points used for numerical
 *                    antidifferentiation in parameter space
 * \param [in] npts_Z the number of quadrature points used for numerical
 *                    antidifferentiation in z
 *
 * \pre The patch parameterization must be valid on its full parameter domain.
 * \pre If the patch is trimmed, its trimming curves must bound the intended
 *      interior region in parameter space.
 * \pre The returned value is geometrically meaningful as a volume contribution only
 *      when this patch is interpreted as part of a closed, consistently oriented
 *      boundary that uses the same lower integration bound.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSPatch<T, 3>::PointType>>
LambdaRetType evaluate_volume_integral(const primal::NURBSPatch<T, 3>& patch,
                                       Lambda&& integrand,
                                       T lower_bound_z,
                                       int npts_Q,
                                       int npts_P = 0,
                                       int npts_Z = 0)
{
  static_assert(detail::internal::is_integrable_v<T, LambdaRetType>,
                "evaluate_integral methods require addition and scalar multiplication for lambda "
                "function return type");

  if(npts_P <= 0)
  {
    npts_P = npts_Q;
  }
  if(npts_Z <= 0)
  {
    npts_Z = npts_Q;
  }

  return detail::evaluate_volume_integral_component(patch,
                                                    std::forward<Lambda>(integrand),
                                                    lower_bound_z,
                                                    npts_Q,
                                                    npts_P,
                                                    npts_Z);
}

/*!
 * \brief Evaluate a scalar volume integral over a collection of Bezier patches.
 *
 * The result is obtained by summing the Stokes-based contribution from each
 * patch in the collection.
 *
 * \param [in] patches the patch collection
 * \param [in] integrand callable representing the integrand
 * \param [in] npts_uv the number of quadrature points in each patch parameter direction
 * \param [in] npts_z the number of quadrature points used for numerical
 *                    antidifferentiation in z
 *
 * \pre Each patch parameterization must be valid on its full parameter domain.
 * \pre The patch collection must represent a closed, consistently oriented
 *      boundary of the target volume.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename BezierPatch<T, 3>::PointType>>
LambdaRetType evaluate_volume_integral(const axom::Array<BezierPatch<T, 3>>& patches,
                                       Lambda&& integrand,
                                       int npts_uv,
                                       int npts_z = 0)
{
  static_assert(detail::internal::is_integrable_v<T, LambdaRetType>,
                "evaluate_integral methods require addition and scalar multiplication for lambda "
                "function return type");

  if(npts_z <= 0)
  {
    npts_z = npts_uv;
  }

  if(patches.empty())
  {
    return LambdaRetType {};
  }

  T lower_bound_z = patches[0].boundingBox().getMin()[2];
  for(int i = 1; i < patches.size(); ++i)
  {
    lower_bound_z = axom::utilities::min(lower_bound_z, patches[i].boundingBox().getMin()[2]);
  }

  LambdaRetType total_integral = LambdaRetType {};
  for(int i = 0; i < patches.size(); ++i)
  {
    total_integral +=
      detail::evaluate_volume_integral_component(patches[i], integrand, lower_bound_z, npts_uv, npts_z);
  }

  return total_integral;
}

/*!
 * \brief Evaluate a scalar volume integral over a collection of NURBS patches.
 *
 * The result is obtained by summing the Stokes-based contribution from each
 * patch in the collection.
 *
 * \param [in] patches the patch collection
 * \param [in] integrand callable representing the integrand
 * \param [in] npts_Q the number of quadrature points on each trimming curve or
 *                    in each parametric direction for untrimmed Bezier pieces
 * \param [in] npts_P the number of quadrature points used for numerical
 *                    antidifferentiation in parameter space
 * \param [in] npts_Z the number of quadrature points used for numerical
 *                    antidifferentiation in z
 *
 * \pre Each patch parameterization must be valid on its full parameter domain.
 * \pre Any trimmed patch in the array must have trimming curves that bound its
 *      intended interior region in parameter space.
 * \pre The patch collection must represent a closed, consistently oriented
 *      boundary of the target volume.
 */
template <typename Lambda,
          typename T,
          typename LambdaRetType = std::invoke_result_t<Lambda, typename NURBSPatch<T, 3>::PointType>>
LambdaRetType evaluate_volume_integral(const axom::Array<NURBSPatch<T, 3>>& patches,
                                       Lambda&& integrand,
                                       int npts_Q,
                                       int npts_P = 0,
                                       int npts_Z = 0)
{
  static_assert(detail::internal::is_integrable_v<T, LambdaRetType>,
                "evaluate_integral methods require addition and scalar multiplication for lambda "
                "function return type");

  if(npts_P <= 0)
  {
    npts_P = npts_Q;
  }
  if(npts_Z <= 0)
  {
    npts_Z = npts_Q;
  }

  if(patches.empty())
  {
    return LambdaRetType {};
  }

  T lower_bound_z = patches[0].boundingBox().getMin()[2];
  for(int i = 1; i < patches.size(); ++i)
  {
    lower_bound_z = axom::utilities::min(lower_bound_z, patches[i].boundingBox().getMin()[2]);
  }

  LambdaRetType total_integral = LambdaRetType {};
  for(int i = 0; i < patches.size(); ++i)
  {
    total_integral += detail::evaluate_volume_integral_component(patches[i],
                                                                 integrand,
                                                                 lower_bound_z,
                                                                 npts_Q,
                                                                 npts_P,
                                                                 npts_Z);
  }

  return total_integral;
}
//@}

}  // namespace primal
}  // end namespace axom

#endif
