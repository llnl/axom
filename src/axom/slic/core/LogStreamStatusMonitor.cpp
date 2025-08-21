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

namespace axom
{
namespace slic
{

LogStreamStatusMonitor::LogStreamStatusMonitor()
    : m_streamVec(),
      m_useMPI(false),
      m_mpiComms()
      //m_mpiComm(MPI_COMM_NULL)
{
}

void LogStreamStatusMonitor::addStream(LogStream* ls)
{
    m_streamVec.push_back(ls);
    if (ls->useMpi() == true)
    {
        m_useMPI = true;
        m_mpiComms.push_back(ls->communicator()->comm());
    }
    // also at some point merge communicators
}

}
}