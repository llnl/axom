// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 *******************************************************************************
 * \file NonCollectiveRootCommunicator.hpp
 *
 * \brief This file contains the class definition of the 
 *  NonCollectiveRootCommunicator.
 *******************************************************************************
 */

#include "axom/lumberjack/Lumberjack.hpp"
#include "axom/lumberjack/Communicator.hpp"

#include <memory>
#include <vector>

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class NonCollectiveRootCommunicator
 *
 * \brief Pushes messages from any rank to root non-collectively.
 *
 * This communicator is intended for best-effort logging in abort scenarios.
 * Non-root sends are non-blocking. Send buffers stay valid until MPI completes
 * or finalization releases ownership to avoid waiting on root.
 *******************************************************************************
 */
class NonCollectiveRootCommunicator : public axom::lumberjack::Communicator
{
public:
  /*!
   *****************************************************************************
   * \brief Called to initialize the Communicator.
   *
   * This performs any setup work the Communicator needs before doing any work.
   * It is required that this is called before using the Communicator.
   *
   * \param [in] comm The MPI Communicator
   * \param [in] ranksLimit Upper Limit on number of ranks that are tracked per
   *  Message.
   *****************************************************************************
   */
  void initialize(MPI_Comm comm, int ranksLimit);

  /*!
   *****************************************************************************
   * \brief Called to finalize the Communicator.
   *
   * This performs any cleanup work the Communicator needs to do before going
   * away.It is required that this is the last function called by the
   * Communicator.
   *****************************************************************************
   */
  void finalize();

  /*!
   *****************************************************************************
   * \brief Returns the MPI communicator
   *****************************************************************************
   */
  MPI_Comm comm();

  /*!
   *****************************************************************************
   * \brief Returns the MPI rank of this node
   *****************************************************************************
   */
  int rank();

  /*!
   *****************************************************************************
   * \brief Sets the rank limit.
   *
   * This is the limit on how many ranks generated a given message are
   * individually tracked per Message.  After the limit has been reached, only
   * the Message::rankCount is incremented.
   *
   * \param [in] value Limits how many ranks are tracked per Message.
   *****************************************************************************
   */
  void ranksLimit(int value);

  /*!
   *****************************************************************************
   * \brief Returns the rank limit.
   *
   * This is the limit on how many ranks generated a given message are
   * individually tracked per Message.  After the limit has been reached, only
   * the Message::rankCount is incremented.
   *****************************************************************************
   */
  int ranksLimit();

  /*!
   *****************************************************************************
   * \brief Function used by the Lumberjack class to indicate how many
   * individual pushes fully flush all currently held Message classes to the
   * root node. The Communicator class's tree structure dictates this.
   *****************************************************************************
   */
  int numPushesToFlush();

  /*!
   *****************************************************************************
   * \brief This pushes all messages to the root node.
   *
   * All messages are pushed to the root node. This is the same as
   * RootCommunicator::pushMessagesFully for this Communicator.
   *
   * \param [in] packedMessagesToBeSent All of this rank's Message classes
   *  packed into a single buffer.
   * \param [in,out] receivedPackedMessages Received packed message buffers from
   *  this nodes children.
   *****************************************************************************
   */
  void push(const char* packedMessagesToBeSent, std::vector<const char*>& receivedPackedMessages);

  /*!
   *****************************************************************************
   * \brief Function used by the Lumberjack to indicate whether this node should
   *  be outputting messages. Only the root node outputs messages.
   *****************************************************************************
   */
  bool isOutputNode();

  /*!
   *****************************************************************************
   * \brief This function returns a start time that is consistent across ranks.
   * This time corresponds to the time that the Communicator is initialized
   *
   * \return Double value that corresponds to a global start time
   *****************************************************************************
   */
  double startTime();

private:
  struct PendingSend
  {
    MPI_Request request;
    std::unique_ptr<char[]> buffer;
  };

  /*!
   *****************************************************************************
   * \brief Polls pending sends and releases completed send buffers.
   *****************************************************************************
   */
  void drainCompletedSends();

  /*!
   *****************************************************************************
   * \brief Releases pending send state without waiting for incomplete sends.
   *
   * Incomplete send buffers are intentionally leaked to keep them valid for MPI.
   *****************************************************************************
   */
  void releasePendingSends();

  MPI_Comm m_mpiComm;
  int m_mpiCommRank;
  int m_mpiCommSize;
  int m_ranksLimit;
  double m_startTime;
  std::vector<PendingSend> m_pendingSends;
};

}  // end namespace lumberjack
}  // end namespace axom
