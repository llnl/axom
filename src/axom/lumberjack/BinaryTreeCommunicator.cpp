// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file BinaryTreeCommunicator.cpp
 *
 * \brief Implementation of the BinaryTreeCommunicator class.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/BinaryTreeCommunicator.hpp"

#include "axom/core/utilities/Utilities.hpp"

#include "mpi.h"

#include "axom/lumberjack/MPIUtility.hpp"
#include "axom/lumberjack/Message.hpp"

namespace axom
{
namespace lumberjack
{
void BinaryTreeCommunicator::initialize(MPI_Comm comm, int ranksLimit)
{
  m_mpiComm = comm;

  MPI_Barrier(m_mpiComm);
  m_startTime = MPI_Wtime();

  MPI_Comm_rank(m_mpiComm, &m_mpiCommRank);
  MPI_Comm_size(m_mpiComm, &m_mpiCommSize);
  m_ranksLimit = ranksLimit;

  // Calculate tree information about this rank
  m_treeHeight = int(axom::utilities::log2(m_mpiCommSize) + 1);
  m_parentRank = (m_mpiCommRank - 1) >> 1;
  m_leftChildRank = (m_mpiCommRank * 2) + 1;
  m_rightChildRank = (m_mpiCommRank * 2) + 2;
  m_childCount = 0;
  if(m_leftChildRank < m_mpiCommSize)
  {
    ++m_childCount;
  }
  else
  {
    m_leftChildRank = -1;
  }
  if(m_rightChildRank < m_mpiCommSize)
  {
    ++m_childCount;
  }
  else
  {
    m_rightChildRank = -1;
  }
}

void BinaryTreeCommunicator::finalize() { }

MPI_Comm BinaryTreeCommunicator::comm() { return m_mpiComm; }

int BinaryTreeCommunicator::rank() { return m_mpiCommRank; }

void BinaryTreeCommunicator::ranksLimit(int value) { m_ranksLimit = value; }

int BinaryTreeCommunicator::ranksLimit() { return m_ranksLimit; }

int BinaryTreeCommunicator::numPushesToFlush() { return m_treeHeight - 1; }

void BinaryTreeCommunicator::push(const char* packedMessagesToBeSent,
                                  std::vector<const char*>& receivedPackedMessages)
{
  MPI_Barrier(m_mpiComm);
  if(m_mpiCommRank != 0)
  {
    if(isPackedMessagesEmpty(packedMessagesToBeSent))
    {
      mpiNonBlockingSendMessages(m_mpiComm, m_parentRank, zeroMessage);
    }
    else
    {
      mpiNonBlockingSendMessages(m_mpiComm, m_parentRank, packedMessagesToBeSent);
    }
  }

  int childrenDoneCount = 0;
  const char* currPackedMessages;
  while(childrenDoneCount < m_childCount)
  {
    currPackedMessages = mpiBlockingReceiveMessages(m_mpiComm);
    if(isPackedMessagesEmpty(currPackedMessages))
    {
      if((currPackedMessages != nullptr) || (currPackedMessages[0] == '\0'))
      {
        delete[] currPackedMessages;
      }
    }
    else
    {
      receivedPackedMessages.push_back(currPackedMessages);
    }
    ++childrenDoneCount;
  }

  MPI_Barrier(m_mpiComm);
}

bool BinaryTreeCommunicator::isOutputNode()
{
  if(m_mpiCommRank == 0)
  {
    return true;
  }
  return false;
}

double BinaryTreeCommunicator::startTime() { return m_startTime; }

}  // end namespace lumberjack
}  // end namespace axom
