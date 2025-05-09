// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Annotations.hpp
 *
 * \brief Defines functions to help with performance annotations
 * 
 * The annotations API and macros are always available but they are effectively no-ops
 * unless axom is built with caliper and adiak support
 */

#ifndef AXOM_CORE_ANNOTATIONS_HPP_
#define AXOM_CORE_ANNOTATIONS_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

#ifdef AXOM_USE_ADIAK
  #include "adiak.hpp"
#endif

#ifdef AXOM_USE_CALIPER
  #include "caliper/cali.h"
#endif

#include <string>
#include <map>

namespace axom
{
namespace utilities
{
namespace annotations
{
namespace detail
{
/// \note Intended to be called from within the axom::utilities::annotations API
#ifdef AXOM_USE_MPI
void initialize_adiak(MPI_Comm comm = MPI_COMM_WORLD);
#else
void initialize_adiak();
#endif

/// \note Intended to be called from within the axom::utilities::annotations API
void initialize_caliper(const std::string& mode);

// Checks if the provided annotation mode is valid
bool is_mode_valid(const std::string& mode);

// Returns a help string for valid annotation modes
std::string mode_help_string();

}  // namespace detail

/// Initializes the Annotation API
#ifdef AXOM_USE_MPI
void initialize(MPI_Comm comm, const std::string& mode);
#endif
void initialize(const std::string& mode);

/// Finalizes and flushes the Annotation API
void finalize();

/// Begins an annotation within a region
void begin(const std::string& name);

/// Ends an annotation within a region
void end(const std::string& name);

/// Declares metadata for this run
template <typename T>
void declare_metadata(const std::string& name, const T& value, std::string category = "")
{
#ifdef AXOM_USE_ADIAK
  detail::initialize_adiak();
  adiak::value(name, value, adiak_general, category);
#else
  AXOM_UNUSED_VAR(name);
  AXOM_UNUSED_VAR(value);
  AXOM_UNUSED_VAR(category);
#endif
}

/**
 * \brief Access registered metadata from adiak (when available)
 * and returns the result as a map of key-value pairs of strings
 * 
 * \note Returns an empty map in configurations without adiak
 */
std::map<std::string, std::string> retrieve_metadata();

}  // namespace annotations
}  // namespace utilities
}  // namespace axom

#endif  // AXOM_CORE_ANNOTATIONS_HPP_
