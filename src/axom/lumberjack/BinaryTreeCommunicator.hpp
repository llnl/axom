// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file BinaryTreeCommunicator.hpp
 *
 * \brief This file contains the class definition of the BinaryTreeCommunicator.
 *******************************************************************************
 */

#ifndef BINARYTREECOMMUNICATOR_HPP
#define BINARYTREECOMMUNICATOR_HPP

#include "mpi.h"

#include "axom/lumberjack/Communicator.hpp"

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class BinaryTreeCommunicator
 *
 * \brief A Communicator that utilizes a binary tree stucture of MPI nodes
 * to scalably pass Message classes.
 *
 *  You will need to add your Communicator using Lumberjack::initialize.
 *
 * \see Communicator Lumberjack
 *******************************************************************************
 */
class BinaryTreeCommunicator : public Communicator
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
   * \param [in] ranksLimit Limit on how many ranks are individually tracked per
   *  Message.
   *****************************************************************************
   */
  void initialize(MPI_Comm comm, int ranksLimit);

  /*!
   *****************************************************************************
   * \brief Called to finalize the Communicator.
   *
   * This performs any cleanup work the Communicator needs to do before going
   * away. It is required that this is the last function called by the
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
   * \param [in] value limits how many ranks are tracked per Message.
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
   *  individual pushes fully flush all currently held Message classes to the
   *  root node. The Communicator class's tree structure dictates this.
   *****************************************************************************
   */
  int numPushesToFlush();

  /*!
   *****************************************************************************
   * \brief All children push their Message classes pushed to their parent and
   *  then combined.
   *
   * \param [in] packedMessagesToBeSent All of this rank's Message classes
   *  packed into a single buffer.
   * \param [in,out] receivedPackedMessages Received packed message buffers from
   *  this nodes children.
   *
   * This is helpful when you want to spread your Lumberjack work over a set of
   * your work. Instead of doing all of Lumberjack's work at one point in your
   * program. This does not guarantee your Message will be ready to output.
   * After this call only Message classes at the root node will be returned from
   * Lumberjack::getMessages.
   *
   *****************************************************************************
   */
  void push(const char* packedMessagesToBeSent, std::vector<const char*>& receivedPackedMessages);

  /*!
   *****************************************************************************
   * \brief Function used by the Lumberjack to indicate whether this node should
   *  be outputting messages. Only the root node outputs messages.
   *
   * \return Boolean indicates whether you should output messages
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
  int m_treeHeight;
  int m_parentRank;
  int m_leftChildRank;
  int m_rightChildRank;
  int m_childCount;
  double m_startTime;
};

}  // end namespace lumberjack
}  // end namespace axom

#endif
