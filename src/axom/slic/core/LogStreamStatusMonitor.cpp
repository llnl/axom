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
  #include "axom/slic/core/LogStreamStatusMonitor.hpp"
#endif

#include <algorithm>

namespace axom
{
namespace slic
{

//------------------------------------------------------------------------------
LogStreamStatusMonitor::LogStreamStatusMonitor()
  : m_streamVec()
#if defined(AXOM_USE_MPI)
  , m_useMPI(false)
  , m_mpiComms()
#endif
{ }

//------------------------------------------------------------------------------
void LogStreamStatusMonitor::addStream(LogStream* ls)
{
  m_streamVec.push_back(ls);
#if defined(AXOM_USE_MPI)
  if(ls->isUsingMPI() == true)
  {
    m_useMPI = true;

    /* 
      Lambda returns true if ranks participating in input MPI communicator are 
      identical to those of the current log stream's MPI commmunicator
    */
    auto commCompareFunc = [&](MPI_Comm& comm) {
      MPI_Group groupCommA, groupCommB;
      MPI_Comm_group(comm, &groupCommA);
      MPI_Comm_group(ls->comm(), &groupCommB);

      int result;
      MPI_Group_compare(groupCommA, groupCommB, &result);

      MPI_Group_free(&groupCommA);
      MPI_Group_free(&groupCommB);
      return result == MPI_IDENT;
    };

    auto it = std::find_if(m_mpiComms.begin(), m_mpiComms.end(), commCompareFunc);

    /* 
      Add the stream's MPI communicator only if its participating ranks differ
      from those of all existing communicators in m_mpiComms.
    */
    if(it == m_mpiComms.end())
    {
      m_mpiComms.push_back(ls->comm());
    }
  }
#endif
}

//------------------------------------------------------------------------------
bool LogStreamStatusMonitor::hasPendingMessages() const
{
#if defined(AXOM_USE_MPI)

  if(!m_useMPI)
  {
    return false;
  }

  int has_pending_messages = 0;

  for(auto& stream : m_streamVec)
  {
    if(stream->isUsingMPI())
    {
      has_pending_messages += static_cast<int>(stream->hasPendingMessages());
    }
  }

  int local_has_pending_messages = has_pending_messages;
  for(auto& comm : m_mpiComms)
  {
    if(comm != MPI_COMM_NULL)
    {
      MPI_Allreduce(&local_has_pending_messages, &has_pending_messages, 1, MPI_INT, MPI_MAX, comm);
    }
  }

  return has_pending_messages > 0;
#else
  return false;
#endif
}

//------------------------------------------------------------------------------
void LogStreamStatusMonitor::finalize()
{
  m_streamVec.clear();
#if defined(AXOM_USE_MPI)
  m_useMPI = false;
  m_mpiComms.clear();
#endif
}

}  // namespace slic
}  // namespace axom
