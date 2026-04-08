// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file evaluate_integral.hpp
 *
 * \brief Header for including integral evaluation functions
 *
 * The API is split into:
 * - evaluate_integral_curve.hpp (line/area)
 * - evaluate_integral_surface.hpp (surface/volume)
 *
 * This split avoids circular include dependencies when patches internally
 * evaluate curve-based integrals (e.g. via trimming curves).
 */

#ifndef PRIMAL_EVAL_INTEGRAL_HPP_
#define PRIMAL_EVAL_INTEGRAL_HPP_

#include "axom/primal/operators/evaluate_integral_curve.hpp"
#include "axom/primal/operators/evaluate_integral_surface.hpp"

#endif
