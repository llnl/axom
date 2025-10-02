// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Combiner.hpp
 *
 * \brief This file contains the abstract base class defining the interface of
 *  all Combiners.
 *******************************************************************************
 */

#ifndef COMBINER_HPP
#define COMBINER_HPP

#include "axom/lumberjack/Message.hpp"

namespace axom
{
namespace lumberjack
{
/*!
 *******************************************************************************
 * \class Combiner
 *
 * \brief Abstract base class defining the interface of all Combiner classes.
 *
 *  Concrete instances need to inherit from this class and implement these
 *  functions. You will need to add your Combiner using Lumberjack::addCombiner
 *
 * \see MessageEqualityCombiner Lumberjack
 *******************************************************************************
 */
class Combiner
{
public:
  /*!
   *****************************************************************************
   * \brief Virtual destructor.
   *****************************************************************************
   */
  virtual ~Combiner() { };

  /*!
   *****************************************************************************
   * \brief Returns the unique string identifier for this combiner. Used by
   *  Lumberjack to differentiate between other combiners.
   *****************************************************************************
   */
  virtual const std::string id() = 0;

  /*!
   *****************************************************************************
   * \brief Function used by Lumberjack to indicate whether two Message classes
   * should be combined.  They are not actually combined by this function.
   *
   * \param [in] leftMessage The left Message to be compared.
   * \param [in] rightMessage The right Message to be compared.
   *****************************************************************************
   */
  virtual bool shouldMessagesBeCombined(const Message& leftMessage, const Message& rightMessage) = 0;

  /*!
   *****************************************************************************
   * \brief Combines the combinee into the combined Message.
   *
   * \param [in,out] combined the Message that will be modified.
   * \param [in] combinee the Message that is combined into the other.
   * \param [in] ranksLimit limit on how many individual ranks are tracked
   *  in combined Messages. Message::rankCount is always incremented.
   *****************************************************************************
   */
  virtual void combine(Message& combined, const Message& combinee, const int ranksLimit) = 0;
};

}  // end namespace lumberjack
}  // end namespace axom

#endif
