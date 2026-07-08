// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file NonCollectiveRootCommunicator.cpp
 *
 * \brief Implementation of the NonCollectiveRootCommunicator class.
 *
 ******************************************************************************
 */

#include <iostream>
#include <cstring>
#include <limits>
#include <utility>

#include "axom/lumberjack/NonCollectiveRootCommunicator.hpp"
#include "axom/lumberjack/MPIUtility.hpp"

namespace axom
{
namespace lumberjack
{
void NonCollectiveRootCommunicator::initialize(MPI_Comm comm, int ranksLimit)
{
  m_pendingSends.clear();
  if(ranksLimit < 1)
  {
    std::cerr << "Error: Ranks limit passed to NonCollectiveRootCommunicator "
              << "is not positive" << std::endl;
  }
  MPI_Comm_dup(comm, &m_mpiComm);

  MPI_Barrier(m_mpiComm);
  m_startTime = MPI_Wtime();

  MPI_Comm_rank(m_mpiComm, &m_mpiCommRank);
  MPI_Comm_size(m_mpiComm, &m_mpiCommSize);
  m_ranksLimit = ranksLimit;
}

void NonCollectiveRootCommunicator::finalize()
{
  releasePendingSends();
  MPI_Comm_free(&m_mpiComm);
}

MPI_Comm NonCollectiveRootCommunicator::comm() { return m_mpiComm; }

int NonCollectiveRootCommunicator::rank() { return m_mpiCommRank; }

void NonCollectiveRootCommunicator::ranksLimit(int value) { m_ranksLimit = value; }

int NonCollectiveRootCommunicator::ranksLimit() { return m_ranksLimit; }

int NonCollectiveRootCommunicator::numPushesToFlush() { return 1; }

void NonCollectiveRootCommunicator::push(const char* packedMessagesToBeSent,
                                         std::vector<const char*>& receivedPackedMessages)
{
  if(m_mpiCommRank == 0)
  {
    const char* currPackedMessages = nullptr;
    bool receive_messages = true;
    while(receive_messages)
    {
      currPackedMessages = mpiBlockingReceiveIfMessagesExist(m_mpiComm);

      if(isPackedMessagesEmpty(currPackedMessages))
      {
        if(currPackedMessages == nullptr)
        {
          receive_messages = false;
        }
        else
        {
          delete[] currPackedMessages;
        }
      }
      else
      {
        receivedPackedMessages.push_back(currPackedMessages);
      }

      currPackedMessages = nullptr;
    }
  }
  else
  {
    drainCompletedSends();
    if(isPackedMessagesEmpty(packedMessagesToBeSent) == false)
    {
      const int messageSize = static_cast<int>(std::strlen(packedMessagesToBeSent));
      PendingSend pendingSend;
      pendingSend.request = MPI_REQUEST_NULL;
      pendingSend.buffer.reset(new char[messageSize + 1]);
      std::memcpy(pendingSend.buffer.get(), packedMessagesToBeSent, messageSize + 1);

      pendingSend.request =
        mpiNonBlockingSendMessagesWithRequest(m_mpiComm, 0, pendingSend.buffer.get());
      m_pendingSends.push_back(std::move(pendingSend));
    }
    drainCompletedSends();
  }
}

void NonCollectiveRootCommunicator::drainCompletedSends()
{
  for(auto it = m_pendingSends.begin(); it != m_pendingSends.end();)
  {
    int complete = 0;
    MPI_Test(&it->request, &complete, MPI_STATUS_IGNORE);
    if(complete)
    {
      it = m_pendingSends.erase(it);
    }
    else
    {
      ++it;
    }
  }
}

void NonCollectiveRootCommunicator::releasePendingSends()
{
  for(auto& pendingSend : m_pendingSends)
  {
    if(pendingSend.request == MPI_REQUEST_NULL)
    {
      continue;
    }

    int complete = 0;
    MPI_Test(&pendingSend.request, &complete, MPI_STATUS_IGNORE);
    if(!complete)
    {
      // Keep the send buffer valid without blocking finalization. This
      // communicator is used for best-effort non-collective error reporting;
      // waiting here could hang if root is no longer receiving.
      MPI_Request_free(&pendingSend.request);
      pendingSend.buffer.release();
    }
  }
  m_pendingSends.clear();
}

bool NonCollectiveRootCommunicator::isOutputNode()
{
  if(m_mpiCommRank == 0)
  {
    return true;
  }
  return false;
}

double NonCollectiveRootCommunicator::startTime() { return m_startTime; }

}  // end namespace lumberjack
}  // end namespace axom
