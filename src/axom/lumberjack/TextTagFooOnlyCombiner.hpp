// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file TextTagFooOnlyCombiner.hpp
 *
 * \brief This file contains the class implementation of the
 * TextTagFooOnlyCombiner.
 *******************************************************************************
 */

#ifndef TextTagFOOONLYCombiner_HPP
#define TextTagFOOONLYCombiner_HPP

#include "axom/lumberjack/TextTagCombiner.hpp"
#include "axom/lumberjack/Message.hpp"

#include <string>

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class TextTagFooOnlyCombiner
 *
 * \brief Combines Message classes if they have a tag of Foo and their 
 *        Message::text and Message::tag are equal.
 *
 *  This class can be added to Lumberjack's Lumberjack by calling
 *  Lumberjack::addCombiner with a
 *  TextTagFooOnlyCombiner instance as its parameter.
 *
 * \see Combiner Lumberjack
 *******************************************************************************
 */
class TextTagFooOnlyCombiner : public TextTagCombiner
{
public:
  TextTagFooOnlyCombiner() { }

  /*!
   *****************************************************************************
   * \brief Function used by Lumberjack to indicate whether a Message class
   * should be considered for this combiner, must have tag "foo".
   *
   * \param [in] m The Message to be considered
   *****************************************************************************
   */
  bool isMessageCandidateForCombiner(const Message& m) { return m.tag() == "foo"; }

private:
  const std::string m_id = "TextTagFooOnlyCombiner";
};

}  // end namespace lumberjack
}  // end namespace axom

#endif
