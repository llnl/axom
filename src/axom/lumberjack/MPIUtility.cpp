// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file MPIUtility.cpp
 *
 * \brief Implementation of the MPI utility methods.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/MPIUtility.hpp"
#include <cstring>

namespace axom
{
namespace lumberjack
{
constexpr int LJ_TAG = 32766;

const char* mpiBlockingReceiveMessages(MPI_Comm comm)
{
  char* charArray = nullptr;
  int messageSize = -1;
  MPI_Status mpiStatus;

  // Get size and source of MPI message
  MPI_Probe(MPI_ANY_SOURCE, LJ_TAG, comm, &mpiStatus);
  MPI_Get_count(&mpiStatus, MPI_CHAR, &messageSize);

  // Setup where to receive the char array
  charArray = new char[messageSize + 1];
  charArray[messageSize] = '\0';

  // Receive packed Message
  MPI_Recv(charArray, messageSize, MPI_CHAR, mpiStatus.MPI_SOURCE, LJ_TAG, comm, &mpiStatus);

  return charArray;
}

const char* mpiBlockingReceiveIfMessagesExist(MPI_Comm comm)
{
  char* charArray = nullptr;
  int messageSize = -1;
  MPI_Status mpiStatus;

  // Get size and source of MPI message
  int mpiFlag = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, LJ_TAG, comm, &mpiFlag, &mpiStatus);

  if(mpiFlag == 1)
  {
    MPI_Get_count(&mpiStatus, MPI_CHAR, &messageSize);

    // Setup where to receive the char array
    charArray = new char[messageSize + 1];
    charArray[messageSize] = '\0';

    // Receive packed Message
    MPI_Recv(charArray, messageSize, MPI_CHAR, mpiStatus.MPI_SOURCE, LJ_TAG, comm, &mpiStatus);
  }

  return charArray;
}

void mpiNonBlockingSendMessages(MPI_Comm comm, int destinationRank, const char* packedMessagesToBeSent)
{
  MPI_Request mpiRequest;
  MPI_Isend(const_cast<char*>(packedMessagesToBeSent),
            strlen(packedMessagesToBeSent),
            MPI_CHAR,
            destinationRank,
            LJ_TAG,
            comm,
            &mpiRequest);
  MPI_Request_free(&mpiRequest);
}

}  // end namespace lumberjack
}  // end namespace axom
