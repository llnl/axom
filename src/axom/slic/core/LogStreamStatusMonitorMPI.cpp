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

        MPI_Group mpi_group;
        MPI_Comm_group(ls->comm(), &mpi_group);

        lumberjack::addGroupToComm(&m_mpiComm, mpi_group);
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

//------------------------------------------------------------------------------
void LogStreamStatusMonitorMPI::finalize()
{
    if (m_mpiComm != MPI_COMM_WORLD && m_mpiComm != MPI_COMM_NULL)
    {
        MPI_Comm_free(&m_mpiComm);
    }
}

}
}