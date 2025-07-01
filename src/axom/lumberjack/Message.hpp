// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Message.hpp
 *
 * \brief This file contains the class definition of the Message.
 *******************************************************************************
 */

#ifndef MESSAGE_HPP
#define MESSAGE_HPP

#include <cstring>
#include <string>
#include <vector>

namespace axom
{
namespace lumberjack
{
/*!
 *****************************************************************************
 * \brief Message indicating no messages need to be sent from child node.
 *****************************************************************************
 */
extern const char* const zeroMessage;

/*!
 *****************************************************************************
 * \brief Delimiter used for packing messages and separating their members.
 *****************************************************************************
 */
const char memberDelimiter = '*';

/*!
 *****************************************************************************
 * \brief Delimiter used for packing messages and separating their individual
 *  tracked ranks.
 *****************************************************************************
 */
const char rankDelimiter = ',';

/*!
 *******************************************************************************
 * \class Message
 *
 * \brief Holds information about messages and where they came from.
 *
 *  This class holds all the information about an individual message and where
 *  it came from, such as rank, file name, and line number.
 *
 * \see Combiner Lumberjack
 *******************************************************************************
 */
class Message
{
public:
  // Constructors
  /*!
   *****************************************************************************
   * \brief Basic constructor where everything defaults to nothing.
   *****************************************************************************
   */
  Message()
    : m_text("")
    , m_ranks()
    , m_ranksLimitReached(false)
    , m_count(0)
    , m_fileName("")
    , m_lineNumber(0)
    , m_level(0)
    , m_creationTime(0)
    , m_tag("")
  { }

  /*!
   *****************************************************************************
   * \brief Constructor where you can specify all values for a Message that
   *  originated from a specific rank.
   *
   * \param [in] text Actual text of the Message.
   * \param [in] rank The rank where the Message originated.
   * \param [in] fileName The file name where the Message originated.
   * \param [in] lineNumber The line number where the Message originated.
   * \param [in] level The level of the severity of the Message.
   * \param [in] tag The tag of where the Message originated.
   *****************************************************************************
   */
  Message(const std::string& text,
          int rank,
          const std::string& fileName,
          int lineNumber,
          int level,
          double creationTime,
          const std::string& tag)
    : m_text(text)
    , m_ranks(1, rank)
    , m_ranksLimitReached(false)
    , m_count(1)
    , m_fileName(fileName)
    , m_lineNumber(lineNumber)
    , m_level(level)
    , m_creationTime(creationTime)
    , m_tag(tag)
  { }

  /*!
   *****************************************************************************
   * \brief Constructor where you can specify all values for a Message that
   *  originated from a multiple ranks.
   *
   * \param [in] text Actual text of the Message.
   * \param [in] ranks The rank where the Message originated.
   * \param [in] count Total number of instances of this Message.
   * \param [in] ranksLimit Limit on how many ranks are individually tracked per
   *  Message.
   * \param [in] fileName The file name where the Message originated.
   * \param [in] lineNumber The line number where the Message originated.
   * \param [in] level The level of the severity of the Message.
   * \param [in] tag The tag of where the Message originated.
   *****************************************************************************
   */
  Message(const std::string& text,
          const std::vector<int>& ranks,
          int count,
          int ranksLimit,
          const std::string& fileName,
          int lineNumber,
          int level,
          double creationTime,
          const std::string& tag)
    : m_text(text)
    , m_ranks()
    , m_ranksLimitReached(false)
    , m_count(0)
    , m_fileName(fileName)
    , m_lineNumber(lineNumber)
    , m_level(level)
    , m_creationTime(creationTime)
    , m_tag(tag)
  {
    addRanks(ranks, count, ranksLimit);
  }

  // Getters
  /*!
   *****************************************************************************
   * \brief Returns the text of the Message.
   *****************************************************************************
   */
  std::string text() const;

  /*!
   *****************************************************************************
   * \brief Returns the vector of the ranks where this Message originated.
   *****************************************************************************
   */
  std::vector<int> ranks() const;

  /*!
   *****************************************************************************
   * \brief Returns the total count of this Message.
   *****************************************************************************
   */
  int count() const;

  /*!
   *****************************************************************************
   * \brief Returns a string of ranks delimited by ',' unless otherwise
   *  specified.
   *
   * \param [in] delimiter The delimiter used to separate the ranks in returned
   *  string.
   *****************************************************************************
   */
  std::string stringOfRanks(std::string delimiter = ",") const;

  /*!
   *****************************************************************************
   * \brief Returns the file name of where this Message originated.
   *****************************************************************************
   */
  std::string fileName() const;

  /*!
   *****************************************************************************
   * \brief Returns the line number of where this Message originated.
   *****************************************************************************
   */
  int lineNumber() const;

  /*!
   *****************************************************************************
   * \brief Returns the level of the severity of the Message.
   *****************************************************************************
   */
  int level() const;

  /*!
   *****************************************************************************
   * \brief Returns the creation time of the Message.
   *****************************************************************************
   */
  double creationTime() const;

  /*!
   *****************************************************************************
   * \brief Returns the tag of where the Message originated.
   *****************************************************************************
   */
  std::string tag() const;

  // Setters

  /*!
   *****************************************************************************
   * \brief Sets a new text for this Message.
   *
   * \param [in] newText The new text to be set for this Message.
   *****************************************************************************
   */
  void text(const std::string& newText);

  /*!
   *****************************************************************************
   * \brief Sets a new file name for this Message.
   *
   * \param [in] newFileName The new file name to be set for this Message.
   *****************************************************************************
   */
  void fileName(const std::string& newFileName);

  /*!
   *****************************************************************************
   * \brief Sets a new line number for this Message.
   *
   * \param [in] newLineNumber The new line number to be set for this message.
   *****************************************************************************
   */
  void lineNumber(int newLineNumber);

  /*!
   *****************************************************************************
   * \brief Sets a new severity level for this Message.
   *
   * \param [in] newLevel The new severity level for this Message.
   *****************************************************************************
   */
  void level(int newLevel);

  /*!
   *****************************************************************************
   * \brief Sets the creation time of the Message.
   *****************************************************************************
   */
  void creationTime(double newCreationTime);

  /*!
   *****************************************************************************
   * \brief Sets a new tag of where the Message originated.
   *
   * \param [in] newTag The new tag of where the Message originated.
   *****************************************************************************
   */
  void tag(const std::string& newTag);

  /*!
   *****************************************************************************
   * \brief Adds a rank to this Message.
   *
   * \param [in] newRank The new rank to be added.
   * \param [in] ranksLimit Limits how many ranks are tracked per  Message.
   *****************************************************************************
   */
  void addRank(int newRank, int ranksLimit);

  /*!
   *****************************************************************************
   * \brief Adds multiple ranks to this Message.  count tracks how many times
   * this message has occurred since duplicates are being filtered.
   *
   * \param [in] newRanks The new ranks to be added.
   * \param [in] count Count to add to Message::count
   * \param [in] ranksLimit Limits how many ranks are tracked per Message.
   *****************************************************************************
   */
  void addRanks(const std::vector<int>& newRanks, int count, int ranksLimit);

  // utilities

  /*!
   *****************************************************************************
   * \brief Returns a string of all information about this Message packed into a
   *  string.
   *
   * The Message is packed into a string utilizing the following format:
   *  <ranks delimited by ,>*<rank count>*<file name>*<line number>
   *  *<message level>*<tag>*<text>
   *
   *****************************************************************************
   */
  std::string pack();

  /*!
   *****************************************************************************
   * \brief Overrides the information in this Message with the given packed
   *  string.
   *
   * \param [in] packedMessage Packed Message containing the new information.
   * \param [in] ranksLimit Limits how many ranks are tracked per Message.
   *
   * The Message is unpacked from a string utilizing the format defined by
   * pack()
   *****************************************************************************
   */
  void unpack(const std::string& packedMessage, int ranksLimit);

private:
  void unpackRanks(const std::string& ranksString, int ranksLimit);

  std::string m_text;
  std::vector<int> m_ranks;
  bool m_ranksLimitReached;
  int m_count;
  std::string m_fileName;
  int m_lineNumber;
  int m_level;
  double m_creationTime;
  std::string m_tag;
};

/*!
 *****************************************************************************
 * \brief This packs all given Message classes into one const char buffer
 *
 * The messages are packed into the following format:
 *  <message count>[*<packed message size>*<packed message>]...
 * This function does not alter the messages vector.
 *
 * \param [in] messages Message classes to be packed for sending
 *
 * \return Packed char array of all given messages
 * \note It is the caller's responsibility to deallocate the returned buffer
 *****************************************************************************
 */
const char* packMessages(const std::vector<Message*>& messages);

/*!
 *****************************************************************************
 * \brief This unpacks the given const char buffer and adds the created Messages
 *  classes to the given vector.
 *
 * The messages are packed into the following format:
 *  <message count>[*<packed message size>*<packed message>]...
 * This function only adds to the messages vector and does not alter the
 * packagedMessages parameter.
 *
 * \note It is the caller's responsibility to deallocate the Message added 
 * to the \a messages vector
 * 
 * \param [in,out] messages Vector to append created messages to
 * \param [in]  packedMessages Packed messages to be unpacked
 * \param [in]  ranksLimit Limits how many ranks are tracked per Message.
 *****************************************************************************
 */
void unpackMessages(std::vector<Message*>& messages, const char* packedMessages, const int ranksLimit);

/*!
 *****************************************************************************
 * \brief This checks if a given set of packed messages is empty.
 *
 * \param [in]  packedMessages Packed messages to be checked for empty.
 *
 * \return Whether or not the given packed messages is empty.
 *****************************************************************************
 */
inline bool isPackedMessagesEmpty(const char* packedMessages)
{
  return (packedMessages == nullptr) || (packedMessages[0] == '\0') ||
    (strcmp(packedMessages, zeroMessage) == 0);
}

}  // end namespace lumberjack
}  // end namespace axom

#endif
