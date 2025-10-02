// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file BinaryTreeCommunicator.cpp
 *
 * \brief Implementation of the Lumberjack class.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/Lumberjack.hpp"
#include <algorithm>
#include <iostream>

namespace axom
{
namespace lumberjack
{
void Lumberjack::initialize(Communicator* communicator, int ranksLimit, bool isCommunicatorOwned)
{
  if(m_isInitialized == false)
  {
    m_isInitialized = true;
    m_communicator = communicator;
    m_isCommunicatorOwned = isCommunicatorOwned;
    m_ranksLimit = ranksLimit;
    m_combiners.push_back(new TextTagCombiner);
  }
  else
  {
    std::cerr << "Lumberjack::initialize called more than once" << std::endl;
  }
}

void Lumberjack::finalize()
{
  if(m_isCommunicatorOwned && m_communicator != nullptr)
  {
    m_communicator->finalize();
    delete m_communicator;
  }
  m_communicator = nullptr;
  clearCombiners();
  clearMessages();
}

void Lumberjack::addCombiner(Combiner* combiner)
{
  bool identifierFound = false;
  for(int i = 0; i < (int)m_combiners.size(); ++i)
  {
    if(m_combiners[i]->id() == combiner->id())
    {
      identifierFound = true;
      break;
    }
  }
  if(!identifierFound)
  {
    m_combiners.push_back(combiner);
  }
}

void Lumberjack::removeCombiner(const std::string& combinerIdentifier)
{
  int combinerToBeRemoved = -1;
  for(int i = 0; i < (int)m_combiners.size(); ++i)
  {
    if(m_combiners[i]->id() == combinerIdentifier)
    {
      delete m_combiners[i];
      combinerToBeRemoved = i;
      break;
    }
  }
  if(combinerToBeRemoved != -1)
  {
    m_combiners.erase(m_combiners.begin() + combinerToBeRemoved);
  }
}

void Lumberjack::clearCombiners()
{
  for(int i = 0; i < (int)m_combiners.size(); ++i)
  {
    delete m_combiners[i];
  }
  m_combiners.clear();
}

const std::vector<Message*>& Lumberjack::getMessages() const { return m_messages; }

void Lumberjack::ranksLimit(int value)
{
  m_ranksLimit = value;
  m_communicator->ranksLimit(value);
}

int Lumberjack::ranksLimit() { return m_ranksLimit; }

void Lumberjack::clearMessages()
{
  for(int i = 0; i < (int)m_messages.size(); ++i)
  {
    delete m_messages[i];
  }
  m_messages.clear();
}

void Lumberjack::queueMessage(const std::string& text, double creationTime)
{
  queueMessage(text, "", -1, 0, creationTime, "");
}

void Lumberjack::queueMessage(const std::string& text,
                              const std::string& fileName,
                              const int lineNumber,
                              int level,
                              double creationTime,
                              const std::string& tag)
{
  const double elapsedTime = creationTime - m_communicator->startTime();
  Message* mi =
    new Message(text, m_communicator->rank(), fileName, lineNumber, level, elapsedTime, tag);
  m_messages.push_back(mi);
}

void Lumberjack::queueMessage(const std::string& text,
                              const std::vector<int>& ranks,
                              const int count,
                              const std::string& fileName,
                              const int lineNumber,
                              int level,
                              double creationTime,
                              const std::string& tag)
{
  const double elapsedTime = creationTime - m_communicator->startTime();
  Message* mi =
    new Message(text, ranks, count, m_ranksLimit, fileName, lineNumber, level, elapsedTime, tag);
  m_messages.push_back(mi);
}

void Lumberjack::pushMessagesOnce()
{
  const char* packedMessagesToBeSent = "";
  if(!m_communicator->isOutputNode())
  {
    combineMessages();
    packedMessagesToBeSent = packMessages(m_messages);
    clearMessages();
  }
  std::vector<const char*> receivedPackedMessages;

  m_communicator->push(packedMessagesToBeSent, receivedPackedMessages);

  if(!m_communicator->isOutputNode() && !isPackedMessagesEmpty(packedMessagesToBeSent))
  {
    delete[] packedMessagesToBeSent;
  }

  for(int i = 0; i < (int)receivedPackedMessages.size(); ++i)
  {
    unpackMessages(m_messages, receivedPackedMessages[i], m_ranksLimit);
    delete[] receivedPackedMessages[i];
  }
  receivedPackedMessages.clear();

  combineMessages();

  std::sort(m_messages.begin(), m_messages.end(), [](Message* const a, Message* const b) {
    return a->creationTime() < b->creationTime();
  });
}

void Lumberjack::pushMessagesFully()
{
  const char* packedMessagesToBeSent = "";
  std::vector<const char*> receivedPackedMessages;
  int numPushesToFlush = m_communicator->numPushesToFlush();
  for(int i = 0; i < numPushesToFlush; ++i)
  {
    if(!m_communicator->isOutputNode())
    {
      combineMessages();
      packedMessagesToBeSent = packMessages(m_messages);
      clearMessages();
    }

    m_communicator->push(packedMessagesToBeSent, receivedPackedMessages);

    if(!m_communicator->isOutputNode() && !isPackedMessagesEmpty(packedMessagesToBeSent))
    {
      delete[] packedMessagesToBeSent;
    }

    for(int i = 0; i < (int)receivedPackedMessages.size(); ++i)
    {
      unpackMessages(m_messages, receivedPackedMessages[i], m_ranksLimit);
      delete[] receivedPackedMessages[i];
    }
    receivedPackedMessages.clear();
  }

  combineMessages();

  std::sort(m_messages.begin(), m_messages.end(), [](Message* const a, Message* const b) {
    return a->creationTime() < b->creationTime();
  });
}

bool Lumberjack::isOutputNode() { return m_communicator->isOutputNode(); }

void Lumberjack::setCommunicator(Communicator* communicator, bool isCommunicatorOwned)
{
  if(m_isCommunicatorOwned && m_communicator != nullptr)
  {
    m_communicator->finalize();
    delete m_communicator;
  }
  m_isCommunicatorOwned = isCommunicatorOwned;
  m_communicator = communicator;
}

Communicator* Lumberjack::getCommunicator() { return m_communicator; }

bool Lumberjack::isCommunicatorOwned() { return m_isCommunicatorOwned; }

void Lumberjack::combineMessages()
{
  int messagesSize = (int)m_messages.size();
  if(messagesSize < 2)
  {
    return;
  }

  std::vector<Message*> finalMessages;
  std::vector<int> indexesToBeDeleted;
  int combinersSize = (int)m_combiners.size();

  if(combinersSize == 0)
  {
    return;
  }

  bool combinedMessage = false;
  finalMessages.push_back(m_messages[0]);
  for(int allIndex = 1; allIndex < messagesSize; ++allIndex)
  {
    combinedMessage = false;
    for(int finalIndex = 0; finalIndex < (int)finalMessages.size(); ++finalIndex)
    {
      for(int combinerIndex = 0; combinerIndex < combinersSize; ++combinerIndex)
      {
        if(m_combiners[combinerIndex]->shouldMessagesBeCombined(*finalMessages[finalIndex],
                                                                *m_messages[allIndex]))
        {
          m_combiners[combinerIndex]->combine(*finalMessages[finalIndex],
                                              *m_messages[allIndex],
                                              m_ranksLimit);
          indexesToBeDeleted.push_back(allIndex);
          combinedMessage = true;
          break;
        }
      }
      if(combinedMessage)
      {
        break;
      }
    }
    if(!combinedMessage)
    {
      finalMessages.push_back(m_messages[allIndex]);
    }
  }

  for(int i = 0; i < (int)indexesToBeDeleted.size(); ++i)
  {
    delete m_messages[indexesToBeDeleted[i]];
  }
  m_messages.swap(finalMessages);
}
}  // end namespace lumberjack
}  // end namespace axom
