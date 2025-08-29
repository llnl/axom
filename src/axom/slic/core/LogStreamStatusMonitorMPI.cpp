// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file LogStreamStatusMonitorMPI.cpp
 *
 * \brief Implementation of the LogStreamStatusMonitorMPI class.
 *
 ******************************************************************************
 */

#include <mpi.h>

#include "LogStreamStatusMonitorMPI.hpp"
#include "axom/lumberjack/MPIUtility.hpp"
#include "axom/slic/core/LogStreamStatusMonitor.hpp"

#include <iostream>
#include <algorithm>

namespace axom
{
namespace slic
{

//------------------------------------------------------------------------------
LogStreamStatusMonitorMPI::LogStreamStatusMonitorMPI()
  : LogStreamStatusMonitor(),
    m_useMPI(false),
    m_mpiComm(MPI_COMM_NULL)
{
}

//------------------------------------------------------------------------------
void LogStreamStatusMonitorMPI::addStream(LogStream* ls)
{
  LogStreamStatusMonitor::addStream(ls);
  if (ls->isUsingMPI() == true)
  {
    m_useMPI = true;
    if (m_mpiComm == MPI_COMM_NULL && 
      ls->comm() != MPI_COMM_NULL)
    {
      m_mpiComm = ls->comm();
    }
    else if (m_mpiComm != MPI_COMM_NULL && m_mpiComm != ls->comm()) {
      std::cerr << "ERROR: multiple MPI communicators passed to LogStreamStatusMonitor" << std::endl;
    }
  }
}

//------------------------------------------------------------------------------
bool LogStreamStatusMonitorMPI::hasPendingMessages() const
{
  int has_pending_messages = LogStreamStatusMonitor::hasPendingMessages();

  if (m_useMPI)
  {
    int rank, size;
    MPI_Comm_rank(m_mpiComm, &rank);
    MPI_Comm_size(m_mpiComm, &size);

    int local_has_pending_messages = has_pending_messages;

    MPI_Allreduce(&local_has_pending_messages, &has_pending_messages, 1, MPI_INT, MPI_MAX, m_mpiComm);
  }

  return has_pending_messages > 0;
}

}
}