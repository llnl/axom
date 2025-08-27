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
#include "axom/lumberjack/MPIUtility.hpp"

namespace axom
{
namespace slic
{

//------------------------------------------------------------------------------
LogStreamStatusMonitor::LogStreamStatusMonitor()
    : m_streamVec(),
      m_useMPI(false),
      m_mpiComm(MPI_COMM_NULL)
{
}

//------------------------------------------------------------------------------
void LogStreamStatusMonitor::addStream(LogStream* ls)
{
    m_streamVec.push_back(ls);
    if (ls->isUsingMPI() == true)
    {
        m_useMPI = true;

        MPI_Group mpi_group;
        MPI_Comm_group(ls->comm(), &mpi_group);

        lumberjack::addGroupToComm(&m_mpiComm, mpi_group);

        MPI_Group_free(&mpi_group);
    }
}

//------------------------------------------------------------------------------
bool LogStreamStatusMonitor::hasPendingMessages() const
{
    int has_pending_messages = 0;
    for (auto& stream : m_streamVec) {
        has_pending_messages += static_cast<int>(stream->hasPendingMessages());
    }

    if (m_useMPI) {
        int rank, size;
        MPI_Comm_rank(m_mpiComm, &rank);
        MPI_Comm_size(m_mpiComm, &size);

        const int local_has_pending_messages = has_pending_messages;

        MPI_Allreduce(&local_has_pending_messages, &has_pending_messages, 1, MPI_INT, MPI_MAX, m_mpiComm);
    }

    return has_pending_messages > 0;
}

}
}