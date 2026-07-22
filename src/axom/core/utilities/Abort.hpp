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

#pragma once

#include "axom/config.hpp"

namespace axom::utilities
{
/*!
 * \brief Gracefully aborts the application
 */
[[noreturn]] void processAbort();
}  // namespace axom::utilities
