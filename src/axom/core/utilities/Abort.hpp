// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Abort.hpp
 *
 * \brief Declarations for process-abort utilities.
 *
 * \note This header intentionally keeps dependencies minimal and is safe to include
 * in low-level facilities that cannot depend on higher-level utilities headers.
 */

#ifndef AXOM_CORE_UTILITIES_ABORT_HPP_
#define AXOM_CORE_UTILITIES_ABORT_HPP_

#include "axom/config.hpp"

namespace axom::utilities
{
/*!
 * \brief Gracefully aborts the application
 */
[[noreturn]] void processAbort();
}  // namespace axom::utilities

#endif  // AXOM_CORE_UTILITIES_ABORT_HPP_
