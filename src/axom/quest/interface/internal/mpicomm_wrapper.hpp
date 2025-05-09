// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_MPICOMM_WRAPPER_HPP_
#define QUEST_MPICOMM_WRAPPER_HPP_

#include "axom/config.hpp"  // for Axom compile-time definitions

/*!
 * \file
 *
 * \brief Provides aliases for MPI_Comm and MPI_COMM_SELF, which are used in
 *  the interface, when Axom is compiled without MPI. This avoids complicated
 *  ifdef logic at the user-facing API, which was otherwise, necessary.
 *
 * \note When Axom is compiled with MPI, this file just includes the MPI header.
 */

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#else
using MPI_Comm = int;
constexpr int MPI_COMM_SELF = -1;
#endif

#endif /* QUEST_MPICOMM_WRAPPER_HPP_ */
