// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file RAII.hpp
 *
 * \brief Defines several RAII (Resource Acquisition Is Initialization) utility classes to help with 
 * initializing and finalizing libraries and/or components with \a initialize and \a finalize calls
 *
 * For more information about RAII, see: https://en.cppreference.com/w/cpp/language/raii
 */

#ifndef AXOM_CORE_RAII_HPP_
#define AXOM_CORE_RAII_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Annotations.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

namespace axom
{
namespace utilities
{
namespace raii
{
/**
 * RAII wrapper class to initialize and finalize MPI
 * Can also be used when Axom is not configured with MPI
 */
class MPIWrapper
{
public:
  /// Initialize MPI when Axom is configured with MPI, else no-op
  MPIWrapper(int argc, char** argv)
  {
#ifdef AXOM_USE_MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);

    if(mpi_initialized == 0)
    {
      MPI_Init(&argc, &argv);
      m_should_finalize = true;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_numranks);
#else
    AXOM_UNUSED_VAR(argc);
    AXOM_UNUSED_VAR(argv);
#endif
  }

  /// Finalize MPI when Axom is configured with MPI, else no-op
  ~MPIWrapper()
  {
#ifdef AXOM_USE_MPI
    if(m_should_finalize)
    {
      MPI_Finalize();
    }
#endif
  }

  /// Returns number of MPI ranks (1 when Axom not configured with MPI)
  int num_ranks() const { return m_numranks; }

  /// Returns id of current of MPI rank (0 when Axom not configured with MPI)
  int my_rank() const { return m_rank; }

  // Returns true is MPI has been initialized (false when Axom not configured with MPI)
  bool is_initialized() const
  {
#ifdef AXOM_USE_MPI
    int mpi = 0;
    MPI_Initialized(&mpi);
    return (mpi != 0);
#else
    return false;
#endif
  }

private:
  int m_rank {0};
  int m_numranks {1};
  bool m_should_finalize {false};
};

/**
 * Basic RAII wrapper class for initializing and finalizing axom's annotations API.
 * Calls \a annotations::initialize() in constructor and \a annotations::finalize() in destructor
 * 
 * \note Assumes MPI_COMM_WORLD is the MPI communicator. Applications using a custom communicator
 * should directly call the annotations
 * \note In Axom configurations with MPI, this must be called after initializing MPI
 */
class AnnotationsWrapper
{
public:
  AnnotationsWrapper(const std::string& mode)
  {
#ifdef AXOM_USE_MPI
    axom::utilities::annotations::initialize(MPI_COMM_WORLD, mode);
#else
    axom::utilities::annotations::initialize(mode);
#endif
  }

  ~AnnotationsWrapper() { axom::utilities::annotations::finalize(); }
};

}  // namespace raii
}  // namespace utilities
}  // namespace axom

#endif  // AXOM_CORE_RAII_HPP_
