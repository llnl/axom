// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Message.cpp
 *
 * \brief Implementation of the Message class.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/Message.hpp"

#include <algorithm>
#include <iostream>
#include <string>

namespace axom
{
namespace lumberjack
{

// Global workaround: https://rzlc.llnl.gov/jira/browse/ELCAP-851
char zeroMessageStorage[] = "0";
const char* const zeroMessage = zeroMessageStorage;

//Getters

std::string Message::text() const { return m_text; }

std::vector<int> Message::ranks() const { return m_ranks; }

int Message::count() const { return m_count; }

std::string Message::fileName() const { return m_fileName; }

int Message::lineNumber() const { return m_lineNumber; }

int Message::level() const { return m_level; }

double Message::creationTime() const { return m_creationTime; }

std::string Message::tag() const { return m_tag; }

std::string Message::stringOfRanks(std::string delimiter) const
{
  std::string returnString = "";
  int ranksSize = m_ranks.size();
  for(int i = 0; i < ranksSize; ++i)
  {
    returnString += std::to_string(m_ranks[i]);
    if(i < (ranksSize - 1))
    {
      returnString += delimiter;
    }
  }

  // Indicate we might have more ranks we aren't individually tracking
  if(m_ranksLimitReached)
  {
    returnString += "...";
  }
  return returnString;
}

//Setters

void Message::text(const std::string& newText) { m_text = newText; }

void Message::fileName(const std::string& newFileName) { m_fileName = newFileName; }

void Message::lineNumber(int newLineNumber) { m_lineNumber = newLineNumber; }

void Message::level(int newLevel) { m_level = newLevel; }

void Message::creationTime(double newCreationTime) { m_creationTime = newCreationTime; }

void Message::tag(const std::string& newTag) { m_tag = newTag; }

void Message::addRank(int newRank, int ranksLimit)
{
  // If ranksLimit has already been reached don't add newRank to m_ranks
  if(m_ranks.size() < (std::vector<int>::size_type)ranksLimit)
  {
    // If newRank is already in m_ranks then don't add it
    std::vector<int>::iterator iter = std::find(m_ranks.begin(), m_ranks.end(), newRank);
    if((m_ranks.size() == 0) || (iter == m_ranks.end()))
    {
      m_ranks.push_back(newRank);
    }
  }

  if(!m_ranksLimitReached && (m_ranks.size() == (std::vector<int>::size_type)ranksLimit))
  {
    m_ranksLimitReached = true;
  }

  // Always increment message count
  m_count++;
}

void Message::addRanks(const std::vector<int>& newRanks, int count, int ranksLimit)
{
  int newRanksSize = newRanks.size();
  for(int i = 0; i < newRanksSize; ++i)
  {
    // If ranksLimit has already been reached don't add newRank to m_ranks
    if(m_ranks.size() >= (std::vector<int>::size_type)ranksLimit)
    {
      break;
    }
    // If newRank is already in m_ranks then don't add it
    std::vector<int>::iterator iter = std::find(m_ranks.begin(), m_ranks.end(), newRanks[i]);
    if((m_ranks.size() == 0) || (iter == m_ranks.end()))
    {
      m_ranks.push_back(newRanks[i]);
    }
  }

  if(!m_ranksLimitReached && (m_ranks.size() == (std::vector<int>::size_type)ranksLimit))
  {
    m_ranksLimitReached = true;
  }

  // Always increment message count
  m_count += count;
}

// Utilities

std::string Message::pack()
{
  std::string packedMessage;

  int ranksSize = (int)m_ranks.size();
  for(int i = 0; i < ranksSize; ++i)
  {
    packedMessage += std::to_string(m_ranks[i]);
    if(i < (ranksSize - 1))
    {
      packedMessage += rankDelimiter;
    }
  }
  packedMessage += memberDelimiter;

  packedMessage += std::to_string(m_count) + memberDelimiter;

  packedMessage += m_fileName + memberDelimiter;

  if(m_lineNumber > 0)
  {
    packedMessage += std::to_string(m_lineNumber);
  }
  packedMessage += memberDelimiter;

  packedMessage += std::to_string(m_level) + memberDelimiter;

  packedMessage += std::to_string(m_creationTime) + memberDelimiter;

  packedMessage += m_tag + memberDelimiter;

  packedMessage += m_text;

  return packedMessage;
}

void Message::unpack(const std::string& packedMessage, int ranksLimit)
{
  std::size_t start, end;
  std::string section;

  // Grab ranks
  end = packedMessage.find(memberDelimiter);
  if(end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack received a truncated message "
              << "that ended in the ranks section." << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  unpackRanks(packedMessage.substr(0, end), ranksLimit);
  start = end + 1;

  //Grab message count since it can differ from list that is sent
  end = packedMessage.find(memberDelimiter, start);
  if(end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack received a truncated message "
              << "that ended in the rank count section." << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  m_count = std::stoi(packedMessage.substr(start, end - start));
  start = end + 1;

  //Grab file name
  end = packedMessage.find(memberDelimiter, start);
  if(end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack received a truncated message "
              << "that ended in the file name section." << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  m_fileName = packedMessage.substr(start, end - start);
  start = end + 1;

  //Grab line number
  end = packedMessage.find(memberDelimiter, start);
  if(end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack received a truncated message "
              << "that ended in the line number section." << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  if(end - start > 0)
  {
    m_lineNumber = std::stoi(packedMessage.substr(start, end - start));
  }
  start = end + 1;

  //Grab level
  end = packedMessage.find(memberDelimiter, start);
  if(end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack received a truncated message "
              << "that ended in the level section." << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  m_level = std::stoi(packedMessage.substr(start, end - start));
  start = end + 1;

  //Grab creation time
  end = packedMessage.find(memberDelimiter, start);
  if(end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack received a truncated message "
              << "that ended in the level section." << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  m_creationTime = std::stod(packedMessage.substr(start, end - start));
  start = end + 1;

  //Grab tag
  end = packedMessage.find(memberDelimiter, start);
  if(end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack received a truncated message "
              << "that ended in the tag section." << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  m_tag = packedMessage.substr(start, end - start);
  start = end + 1;

  //Grab message
  m_text = packedMessage.substr(start);
}

void Message::unpackRanks(const std::string& ranksString, int ranksLimit)
{
  m_ranks.clear();
  if(ranksString.empty())
  {
    std::cerr << "Error: Lumberjack received an empty rank section." << std::endl;
    return;
  }

  std::size_t start, end = ranksString.find(rankDelimiter);
  start = 0;
  while(true)
  {
    if(end == std::string::npos)
    {
      addRank(std::stoi(ranksString.substr(start)), ranksLimit);
      break;
    }
    else
    {
      addRank(std::stoi(ranksString.substr(start, end - start)), ranksLimit);
    }
    start = end + 1;
    end = ranksString.find(rankDelimiter, start);
  }
}

const char* packMessages(const std::vector<Message*>& messages)
{
  if(messages.size() == 0)
  {
    return zeroMessage;
  }

  int totalSize = 1;  // include size for null terminator

  //Calculate total size of char array after all messages are combined.
  std::vector<std::string> packedMessages;
  std::vector<std::string> sizeStrings;
  int currSize = 0;
  int messageCount = (int)messages.size();
  for(int i = 0; i < messageCount; ++i)
  {
    packedMessages.push_back(messages[i]->pack());
    currSize = packedMessages[i].size();
    sizeStrings.push_back(std::to_string(currSize));
    //           message size + size string size + memberDelimiter size
    totalSize += currSize + sizeStrings[i].size() + 1;
  }

  // Create and calculate size of message count
  std::string messageCountString = std::to_string(messageCount) + memberDelimiter;
  totalSize += messageCountString.size();

  const char* packedMessagesString = new char[totalSize];
  char* packedMessagesIndex = (char*)packedMessagesString;

  // Copy message count to start of packed message
  std::memcpy(packedMessagesIndex, messageCountString.c_str(), messageCountString.size());
  packedMessagesIndex += messageCountString.size();

  for(int i = 0; i < messageCount; ++i)
  {
    // Copy current message size
    std::memcpy(packedMessagesIndex, sizeStrings[i].c_str(), sizeStrings[i].size());
    packedMessagesIndex += sizeStrings[i].size();
    // Copy memberDelimiter
    // ToDo: better way to copy this I'm sure
    std::memcpy(packedMessagesIndex, &memberDelimiter, sizeof(char));
    packedMessagesIndex += 1;
    // Copy packed message
    std::memcpy(packedMessagesIndex, packedMessages[i].c_str(), packedMessages[i].size());
    packedMessagesIndex += packedMessages[i].size();
  }

  packedMessagesIndex[0] = '\0';
  return packedMessagesString;
}

void unpackMessages(std::vector<Message*>& messages, const char* packedMessages, const int ranksLimit)
{
  std::string packedMessagesString = std::string(packedMessages);
  std::size_t start, end;
  std::string tempSubString = "";

  // Get message count
  end = packedMessagesString.find(memberDelimiter);
  tempSubString = packedMessagesString.substr(0, end);
  int messageCount = std::stoi(tempSubString);
  start = end + 1;

  // Grab each message
  Message* message;
  int messageSize;
  for(int j = 0; j < messageCount; ++j)
  {
    //Get current message size
    end = packedMessagesString.find(memberDelimiter, start);
    tempSubString = packedMessagesString.substr(start, end - start);
    messageSize = std::stoi(tempSubString);
    start = end + 1;

    //Create current message and save
    message = new Message();
    tempSubString = packedMessagesString.substr(start, messageSize);
    message->unpack(tempSubString, ranksLimit);
    messages.push_back(message);
    start += messageSize;
  }
}

}  // end namespace lumberjack
}  // end namespace axom
