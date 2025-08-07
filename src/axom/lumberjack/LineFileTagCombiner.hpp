// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file LineFileTagCombiner.hpp
 *
 * \brief This file contains the class implementation of the
 * LineFileTagCombiner.
 *******************************************************************************
 */

#ifndef LINEFILETAGCOMBINER_HPP
#define LINEFILETAGCOMBINER_HPP

#include "axom/lumberjack/Combiner.hpp"
#include "axom/lumberjack/Message.hpp"

#include <string>

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class LineFileTagCombiner
 *
 * \brief Combines Message classes if their Message::fileName, 
 *        Message::lineNumer, and Message::tag are equal.
 *
 *******************************************************************************
 */
class LineFileTagCombiner : public axom::lumberjack::Combiner
{
public:
  LineFileTagCombiner() { }

  /*!
   *****************************************************************************
   * \brief Returns the unique string identifier for this combiner. Used by
   *  Lumberjack to differentiate between other combiners.
   *****************************************************************************
   */
  const std::string id() { return m_id; }

  /*!
   *****************************************************************************
   * \brief Function used by Lumberjack to indicate whether two messages should
   *  be combined.
   *
   * They are not actually combined by this function. Message classes are
   * triggered for combination if Message::fileName, Message::lineNumer, 
   * and Message::tag are equal.
   *
   * \param [in] leftMessage One of the Messages to be compared.
   * \param [in] rightMessage One of the Messages to be compared.
   *****************************************************************************
   */
  bool shouldMessagesBeCombined(const axom::lumberjack::Message& leftMessage,
                                const axom::lumberjack::Message& rightMessage)
  {
    return ((leftMessage.lineNumber() == rightMessage.lineNumber()) &&
            leftMessage.fileName().compare(rightMessage.fileName()) == 0 &&
            leftMessage.tag().compare(rightMessage.tag()) == 0);
  }

  /*!
   *****************************************************************************
   * \brief Combines the combinee into the combined Message.
   *
   * The only things truly combined in this Combiner are the ranks from combinee
   * to combined and the creation time.  The text will not be combined, even if 
   * it is not equal.  Only text from the first message will be saved.
   * The creation time of the first message (lowest creation time) will be saved
   * to the combined message.
   *
   * \param [in,out] combined the Message that will be modified.
   * \param [in] combinee the Message that is combined into the other.
   * \param [in] ranksLimit The limit on how many individual ranks are tracked
   *  in the combined Message. Message::rankCount is always incremented.
   *
   * \pre shouldMessagesBeCombined(combined, combinee) must be true
   *****************************************************************************
   */
  void combine(axom::lumberjack::Message& combined,
               const axom::lumberjack::Message& combinee,
               const int ranksLimit)
  {
    combined.addRanks(combinee.ranks(), combinee.count(), ranksLimit);
    if(combinee.creationTime() < combined.creationTime())
    {
      combined.creationTime(combinee.creationTime());
    }
  }

private:
  const std::string m_id = "LineFileTagCombiner";
};

}  // end namespace lumberjack
}  // end namespace axom

#endif