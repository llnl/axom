// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// This file exists only to enable documentation of accelerated code within Axom. Do not
// delete this file, nor add any lines of code.

/** @page accelerated_list Acceleration in Axom
 *  @brief Acceleration page
 * 
 *  Axom supports acceleration in two ways.  First, the \ref for_all function, \ref execution_space struct,
 *  and \ref synchronize function work with the core memory management routines to allow users to
 *  launch and control code execution on accelerators.  Second, many of the routines throughout Axom
 *  can be called on either the host or the device; that is, they are callable within Axom for_all loops or even
 *  RAJA loops.
*/
