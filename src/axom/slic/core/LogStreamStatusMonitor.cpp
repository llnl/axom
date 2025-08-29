// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file LogStreamStatusMonitor.cpp
 *
 * \brief Implementation of the LogStreamStatusMonitor class.
 *
 ******************************************************************************
 */

#include "LogStreamStatusMonitor.hpp"

#if defined(AXOM_USE_MPI)
#include "axom/lumberjack/MPIUtility.hpp"
#include "axom/slic/core/LogStreamStatusMonitor.hpp"
#endif

#include <iostream>
#include <algorithm>

namespace axom
{
namespace slic
{

//------------------------------------------------------------------------------
LogStreamStatusMonitor::LogStreamStatusMonitor()
  : m_streamVec(),
#if defined(AXOM_USE_MPI)
    m_useMPI(false),
    m_mpiComm(MPI_COMM_NULL)
#endif
{
}

//------------------------------------------------------------------------------
void LogStreamStatusMonitor::addStream(LogStream* ls)
{
  m_streamVec.push_back(ls);
#if defined(AXOM_USE_MPI)
  if (ls->isUsingMPI() == true)
  {
    m_useMPI = true;
    if (m_mpiComm == MPI_COMM_NULL && 
        ls->comm() != MPI_COMM_NULL)
    {
      m_mpiComm = ls->comm();
    }
    else if (m_mpiComm != MPI_COMM_NULL && m_mpiComm != ls->comm()) {
      std::cerr << "ERROR: attempting to register a logstream with an incompatible "
                   "MPI communicator to LogStreamStatusMonitor's existing communicator" << std::endl;
    }
  }
#endif
}

//------------------------------------------------------------------------------
bool LogStreamStatusMonitor::hasPendingMessages() const
{
  int has_pending_messages = 0;

  for (auto& stream : m_streamVec)
  {
    has_pending_messages += static_cast<int>(stream->hasPendingMessages());
  }

#if defined(AXOM_USE_MPI)
  if (m_useMPI)
  {
    int local_has_pending_messages = has_pending_messages;
    MPI_Allreduce(&local_has_pending_messages, &has_pending_messages, 1, MPI_INT, MPI_MAX, m_mpiComm);
  }
#endif

  return has_pending_messages > 0;
}

//------------------------------------------------------------------------------
void LogStreamStatusMonitor::finalize()
{
  m_streamVec.clear();
  m_useMPI = false;
}

}
}
