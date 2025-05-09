// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MessageLevel.h
 *
 */

#ifndef MESSAGELEVEL_H_
#define MESSAGELEVEL_H_

#include <string>

namespace axom
{
namespace slic
{
namespace message
{
/*!
 * \enum MessageType
 *
 * \brief Enumerates the different types of messaging.
 *
 * \note The ordering of the enumerators reflects the level of severity of
 *  a message.
 *
 * \see Logger
 */
enum Level
{
  Error,    //!< ERROR log an error that *may* be recoverable.
  Warning,  //!< WARNING log a warning.
  Info,     //!< INFO log information that is useful for users & developers.
  Debug,    //!< DEBUG log information that is useful for developers.

  Num_Levels  //!< Num_Levels
};

/*!
 * \brief Returns the string name representation of the given level.
 * \param [in] l the level in query.
 * \return name a string corresponding to the name of the given leve.
 * \pre l >= 0 && l < Num_levels
 * \post name="UNKNOWN-LEVEL" \iff (l < 0) || (l >= Num_Levels)
 */
static inline std::string getLevelAsString(Level l)
{
  const std::string MessageLevelName[Num_Levels] = {"ERROR", "WARNING", "INFO", "DEBUG"};
  if(l < 0 || l >= Num_Levels)
  {
    return ("UNKNOWN-LEVEL");
  }

  return (MessageLevelName[l]);
}

} /* namespace message */

namespace inherit
{
/*!
 * \enum flags
 *
 * \brief Holds the bit flags associated with each level.
 */
enum flags
{
  nothing = 0,       //!< nothing, no bit is set.
  error = 1 << 0,    //!< error,   zeroth bit is set.
  warning = 1 << 1,  //!< warning, 1st bit is set.
  info = 1 << 2,     //!< info,    2nd bit is set.
  debug = 1 << 3     //!< debug,   3rd bit is set.
};

/*!
 * \brief Convenience bit mask that turns on all four bits.
 */
static const char everything = error | warning | info | debug;

/*!
 * \brief Convenience bit mask that turns on fatal,error and warning bits.
 */
static const char errors_and_warnings = error | warning;

/*!
 * \brief Array of bit masks corresponding to each level.
 * \note Used to loop through the bit mask associated with each level.
 */
static const flags masks[message::Num_Levels] = {error, warning, info, debug};

} /* namespace inherit */

} /* namespace slic */

} /* namespace axom */

#endif /* MESSAGELEVEL_H_ */
