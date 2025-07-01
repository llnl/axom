// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Lumberjack.hpp
 *
 * \brief This file contains the class definition of the Lumberjack class,
 * which is the main class users will interact with.
 *******************************************************************************
 */

#ifndef LUMBERJACK_HPP
#define LUMBERJACK_HPP

#include <string>

#include "mpi.h"

#include "axom/lumberjack/Combiner.hpp"
#include "axom/lumberjack/Communicator.hpp"
#include "axom/lumberjack/Message.hpp"
#include "axom/lumberjack/TextTagCombiner.hpp"

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class Lumberjack
 *
 * \brief Class that all user interactions with Lumberjack are through.
 *
 *  This class performs all the high level functionality of Lumberjack, such as,
 *  holding and combining Message classes and telling the given Communicator
 *  to push Message classes.
 *
 * \see BinaryTreeCommunicator RootCommunicator Message Combiner
 *******************************************************************************
 */
class Lumberjack
{
public:
  /*!
   *****************************************************************************
   * \brief Called to initialize the Lumberjack.
   *
   * This performs any setup work the Lumberjack needs before doing any work.
   * It is required that this is called before using the Lumberjack.
   *
   * \param [in] communicator The Lumberjack Communicator that will send/receive
   * messages
   * \param [in] ranksLimit Limit on how many ranks are individually tracker per
   * Message.
   * \param [in] isCommunicatorOwned When set to true, Lumberjack will be 
   * responsible for freeing communication object when finalize() is called.
   *****************************************************************************
   */
  void initialize(Communicator* communicator, int ranksLimit, bool isCommunicatorOwned = false);

  /*!
   *****************************************************************************
   * \brief Called to finalize the Lumberjack.
   *
   * This performs any cleanup work the Lumberjack needs to do before going
   * away. It is required that this is the last function called by Lumberjack.
   *****************************************************************************
   */
  void finalize();

  /*!
   *****************************************************************************
   * \brief Adds a Combiner to the Lumberjack.
   *
   * Lumberjack can have multiple Combiner classes.  This is helpful when
   * you have different combining criteria for different Message classes. If
   * you try to add the same Combiner more than once, the second Combiner will
   * not be added.  This is determined solely on Combiner::id. Combiner classes
   * that are added first have precedence. Lumberjack will call delete on the
   * Combiner when finalize is called.
   *
   * \param [in] combiner The Combiner that will be added.
   *****************************************************************************
   */
  void addCombiner(Combiner* combiner);

  /*!
   *****************************************************************************
   * \brief Removes a Combiner from the Lumberjack.
   *
   * \param [in] combinerIdentifier The Combiner identifier that will be
   *  removed.
   *
   * This removes and calls delete on a Combiner held by the Lumberjack. If no
   * Combiner::id matches the given identifier than nothing is removed.
   *****************************************************************************
   */
  void removeCombiner(const std::string& combinerIdentifier);

  /*!
   *****************************************************************************
   * \brief Clears all Combiner classes from the Lumberjack.
   *
   * This removes and calls delete on all Combiner classes held by the
   * Lumberjack.
   *****************************************************************************
   */
  void clearCombiners();

  /*!
   *****************************************************************************
   * \brief Sets the rank limit.
   *
   * This is the limit on how many ranks generated a given message are
   * individually tracked per Message.  After the limit has been reached, only
   * the Message::rankCount is incremented.
   *
   * \param [in] value Limits how many ranks are tracked per  Message.
   *****************************************************************************
   */
  void ranksLimit(int value);

  /*!
   *****************************************************************************
   * \brief Returns the limit on tracked ranks.
   *
   * This is the limit on how many ranks generated a given message are
   * individually tracked per Message.  After the limit has been reached, only
   * the Message::rankCount is incremented.
   *
   * \return The limit on tracked ranks
   *****************************************************************************
   */
  int ranksLimit();

  /*!
   *****************************************************************************
   * \brief Clears all Message classes from the Lumberjack.
   *
   * This removes and calls delete on all Messages classes held by Lumberjack.
   *****************************************************************************
   */
  void clearMessages();

  /*!
   *****************************************************************************
   * \brief Returns a const reference vector with all currently held Message
   * classes.
   *
   * This returns a const reference to the vector that holds the Message classes
   * held by this node. You should check isOutputNode() to indicate if you
   * should output messages depending on your communication scheme.
   *****************************************************************************
   */
  const std::vector<Message*>& getMessages() const;

  /*!
   *****************************************************************************
   * \brief Queues a message to be sent and combined
   *
   * \param [in] text Text of the Message.
   *
   * This creates a Message and queues it to be sent through the Communicator
   * to the root node.  This message may be combined with others depending on
   * the given criteria by the already defined Combiner classes. Depending on
   * the behavior of the Communicator, the message will not be outputted
   * immediately.
   *
   *****************************************************************************
   */
  void queueMessage(const std::string& text, double creationTime = 0.0);

  /*!
   *****************************************************************************
   * \brief Queues a message to be sent and combined
   *
   * This creates a Message and queues it to be sent through the Communicator
   * to the root node.  This message may be combined with others depending on
   * the given criteria by the already defined Combiner classes. Depending on
   * the behavior of the Communicator, the message will not be outputted
   * immediately.
   *
   * \param [in] text Text of the Message
   * \param [in] fileName File name of Message
   * \param [in] lineNumber Line number of Message
   * \param [in] level The level of the severity of the Message.
   * \param [in] tag The tag of where the Message originated.
   *****************************************************************************
   */
  void queueMessage(const std::string& text,
                    const std::string& fileName,
                    const int lineNumber,
                    int level,
                    double creationTime,
                    const std::string& tag);

  /*!
   *****************************************************************************
   * \brief Queues a message to be sent and combined
   *
   * This creates a Message and queues it to be sent through the Communicator
   * to the root node with a pre-existing ranks list.  This message may be 
   * combined with others depending on the given criteria by the already defined 
   * Combiner classes. Depending on the behavior of the Communicator, the message 
   * will not be outputted immediately.
   *
   * \param [in] text Text of the Message
   * \param [in] text Ranks list
   * \param [in] fileName File name of Message
   * \param [in] lineNumber Line number of Message
   * \param [in] level The level of the severity of the Message.
   * \param [in] tag The tag of where the Message originated.
   *****************************************************************************
   */
  void queueMessage(const std::string& text,
                    const std::vector<int>& ranks,
                    const int count,
                    const std::string& fileName,
                    const int lineNumber,
                    int level,
                    double creationTime,
                    const std::string& tag);

  /*!
   *****************************************************************************
   * \brief This pushes all messages once up the Communicator class's tree
   * structure.
   *
   * All of the children push their Message classes to their parent node. This
   * is helpful if you want to spread the work load of Lumberjack over a large
   * set of work. Before and after Message classes are pushed, they are
   * combined.
   *****************************************************************************
   */
  void pushMessagesOnce();

  /*!
   *****************************************************************************
   * \brief This pushes all messages fully up the Communicator class's tree
   *  structure.
   *
   * All messages are continually pushed until all messages are pushed to the
   * root node.After Message classes are pushed once they are combined before
   * being pushed again.
   *****************************************************************************
   */
  void pushMessagesFully();

  /*!
   *****************************************************************************
   * \brief Function indicates whether this node should be outputting messages.
   *  The Communicator class's communication structure dictates this.
   *
   * \return Boolean indicates whether you should output messages
   *****************************************************************************
   */
  bool isOutputNode();

  /*!
   *****************************************************************************
   * \brief set communicator pointer stored in object as well as ownership
   *
   *****************************************************************************
   */
  void setCommunicator(Communicator* communicator, bool isCommunicatorOwned);

  /*!
   *****************************************************************************
   * \brief get communicator pointer stored in object
   *
   * \return Communicator pointer
   *****************************************************************************
   */
  Communicator* getCommunicator();

  /*!
   *****************************************************************************
   * \brief Returns Boolean flag that controls whether communicator instance is 
   * owned by Lumberjack
   *
   * \return Communicator pointer
   *****************************************************************************
   */
  bool isCommunicatorOwned();

private:
  /*!
   *****************************************************************************
   * \brief All Message classes are combined by the currently held Combiner
   *  classes.
   *****************************************************************************
   */
  void combineMessages();

  bool m_isInitialized = false;

  Communicator* m_communicator;
  bool m_isCommunicatorOwned;
  int m_ranksLimit;
  std::vector<Combiner*> m_combiners;
  std::vector<Message*> m_messages;
};

}  // end namespace lumberjack
}  // end namespace axom

#endif
