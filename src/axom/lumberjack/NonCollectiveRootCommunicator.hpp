// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file NonCollectiveRootCommunicator.hpp
 *
 * \brief This file contains the class definition of the 
 *  NonCollectiveRootCommunicator.
 *******************************************************************************
 */

#ifndef NONCOLLECTIVEROOTCOMMUNICATOR_HPP
#define NONCOLLECTIVEROOTCOMMUNICATOR_HPP

#include "axom/lumberjack/Lumberjack.hpp"
#include "axom/lumberjack/Communicator.hpp"

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class NonCollectiveRootCommunicator
 *
 * \brief Based off of RootCommunicator. This communicator pushes
    messages from any rank to root non-collectively, if any messages are sent.
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
  MPI_Comm m_mpiComm;
  int m_mpiCommRank;
  int m_mpiCommSize;
  int m_ranksLimit;
  double m_startTime;
};

}  // end namespace lumberjack
}  // end namespace axom

#endif
