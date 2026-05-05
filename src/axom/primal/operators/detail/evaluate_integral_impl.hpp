// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file evaluate_integral_impl.hpp
 *
 * \brief Header to include implementation files for evaluating integrals
 *
 * Splitting the implementation avoids circular include dependencies when patches
 * internally evaluate curve-based integrals (e.g. via trimming curves).
 */

#ifndef PRIMAL_EVAL_INTEGRAL_IMPL_HPP_
#define PRIMAL_EVAL_INTEGRAL_IMPL_HPP_

#include "axom/primal/operators/detail/evaluate_integral_curve_impl.hpp"
#include "axom/primal/operators/detail/evaluate_integral_surface_impl.hpp"

#endif
