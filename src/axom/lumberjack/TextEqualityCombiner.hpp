// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file TextEqualityCombiner.hpp
 *
 * \brief This file contains the class implementation of the
 * TextEqualityCombiner.
 *******************************************************************************
 */

#ifndef TEXTEQUALITYCOMBINER_HPP
#define TEXTEQUALITYCOMBINER_HPP

#include "axom/lumberjack/Combiner.hpp"
#include "axom/lumberjack/Message.hpp"

#include <string>

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class TextEqualityCombiner
 *
 * \brief Combines Message classes if their Message::text are equal.
 *
 *  This class instance is automatically added to Lumberjack's Lumberjack for
 *  you. If you want it removed call Lumberjack::removeCombiner with the string
 * "TextEqualityCombiner" as it's parameter.
 *
 * \warning Using the TextEqualityCombiner with Message::tag has undefined
 *          behavior.
 *
 * \see Combiner Lumberjack
 *******************************************************************************
 */
class TextEqualityCombiner : public Combiner
{
public:
  TextEqualityCombiner() { }

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
   * triggered for combination if both Message::text are equal.
   *
   * \param [in] leftMessage One of the Messages to be compared.
   * \param [in] rightMessage One of the Messages to be compared.
   *****************************************************************************
   */
  bool shouldMessagesBeCombined(const Message& leftMessage, const Message& rightMessage)
  {
    return (leftMessage.text().compare(rightMessage.text()) == 0);
  }

  /*!
   *****************************************************************************
   * \brief Combines the combinee into the combined Message.
   *
   * The only things truly combined in this Combiner are the ranks from combinee
   * to combined and the creation time, since text is already equal.
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
  void combine(Message& combined, const Message& combinee, const int ranksLimit)
  {
    combined.addRanks(combinee.ranks(), combinee.count(), ranksLimit);
    if(combinee.creationTime() < combined.creationTime())
    {
      combined.creationTime(combinee.creationTime());
    }
  }

private:
  const std::string m_id = "TextEqualityCombiner";
};

}  // end namespace lumberjack
}  // end namespace axom

#endif
