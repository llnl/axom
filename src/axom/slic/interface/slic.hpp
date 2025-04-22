// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file slic.hpp
 */

#ifndef SLIC_HPP_
#define SLIC_HPP_

#include "axom/config.hpp"
#include "axom/slic/core/Logger.hpp"
#include "axom/slic/core/LogStream.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"
#include "axom/slic/core/MessageLevel.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/export/slic.h"

// C/C++ includes
#include <iostream>
#include <sstream>

namespace axom
{
namespace slic
{
struct debug
{
  AXOM_SLIC_EXPORT static bool checksAreErrors;
};

/*!
 * \brief Initializes the SLIC logging environment and allows setting which
 * ranks write messages (root nodes). This defaults to all ranks writing unless
 * selectively set by the user and can be set afterwards via slic::setIsRoot().
 * 
 * \sa slic::setIsRoot()
 * 
 * \param [in] is_root Enables selective logging macros based on root.
 */
void initialize(bool is_root = true);

/*!
 * \brief Checks if the SLIC logging environment is initialized.
 *
 * \return status true if initialized, else, false.
 */
bool isInitialized();

/*!
 * \brief Checks if we are on the root rank. Used for selective logging on root nodes.
 * 
 * \sa slic::initialize()
 * \sa slic::setIsRoot()
 * 
 * \return status true if on root, else, false.
 */
bool isRoot();

/*!
 * \brief Marks this node as root and enables the selective root logging filtering
 * used in the SLIC_*_ROOT_* macros.
 *
 * \param [in] is_root Enables selective logging macros based on root.
 */
void setIsRoot(bool is_root);

/*!
 * \brief Ensures the SLIC logging environment is initialized.
 *
 * If SLIC is not initialized when this method is called, initialize SLIC
 * to print all messages to std::cout and log a warning that prompts
 * developers to properly call slic::initialize() and slic::finalize().
 */
void ensureInitialized();

/*!
 * \brief Creates a new logger associated with the given name.
 *
 * \param [in] name the name to associate with the new logger.
 * \param [in] imask inheritance mask, indicates the log level streams that
 *
 *  will be inherited from the "root" logger. By default, nothing is inherited.
 */
void createLogger(const std::string& name, char imask = inherit::nothing);

/*!
 * \brief Activates the logger associated with the given name.
 *
 * \param [in] name the name of the logger to activate.
 *
 * \return True if the named logger was activated, False otherwise
 */
bool activateLogger(const std::string& name);

/*!
 * \brief Returns the name of the active logger.
 *
 * \return s a string corresponding to the name of the active logger.
 */
std::string getActiveLoggerName();

/*!
 * \brief Sets desired logging level.
 *
 * \param [in] level user-supplied level to log.
 */
void setLoggingMsgLevel(message::Level level);

/*!
 * \brief Gets the current logging level.
 */
message::Level getLoggingMsgLevel();

/*!
 * \brief Toggles the abort behavior for errors on the current active logger.
 *
 * \param [in] status user-supplied flag.
 */
void setAbortOnError(bool status);

/*!
 * \brief Enables aborts on error messages for the current active logger.
 *
 * \note This is equivalent to calling slic::setAbortOnError( true )
 * \post slic::isAbortOnErrorsEnabled() == true.
 * \pre slic::isInitialized() == true
 */
void enableAbortOnError();

/*!
 * \brief Disables aborts on error messages for the current active logger.
 *
 * \note this is equivalent to calling slic::setAbortOnError( false )
 * \post slic::isAbortOnErrorsEnabled() == false.
 * \pre slic::isInitialized() == true
 */
void disableAbortOnError();

/*!
 * \brief Checks whether aborts on errors are enabled for the current logger.
 *
 * \return status true if aborts on errors are enabled, otherwise, false.
 *
 * \pre slic::isInitialized() == true.
 */
bool isAbortOnErrorsEnabled();

/*!
 * \brief Toggles the abort behavior for warnings on the current active logger.
 *
 * \param [in] status user-supplied flag.
 */
void setAbortOnWarning(bool status);

/*!
 * \brief Enables aborts on warnings messages for the current active logger.
 *
 * \note This is equivalent to calling slic::setAbortOnWarning( true )
 * \post slic::isAbortOnWarningsEnabled() == true.
 * \pre slic::isInitialized() == true.
 */
void enableAbortOnWarning();

/*!
 * \brief Disables aborts on warnings messages for the current active logger.
 *
 * \note This is equivalent to calling slic::setAbortOnWarnings( false )
 * \post slic::isAbortOnWarnigsEnabled() == true.
 * \pre slic::isInitialized() == true.
 */
void disableAbortOnWarning();

/*!
 * \brief Checks whether aborts on warnings are enabled for the current logger.
 *
 * \return status true if aborts on warnings are enabled, otherwise, false.
 * \pre slic::isInitialized() == true.
 */
bool isAbortOnWarningsEnabled();

/*!
 * \brief Sets the function to call when program abort is requested
 *
 * \param [in] abort_func The user-specified function to call
 * \pre slic::isInitialized() == true.
 *
 * \warning No collective calls should be made in the given function.
 *          Collective calls may cause the program to hang on abort.
 */
void setAbortFunction(AbortFunctionPtr abort_func);

/*!
 * \brief Adds the given stream to the given level.
 *
 * \param [in] ls pointer to the log stream.
 * \param [in] level the level to log.
 * \pre ls != nullptr
 */
void addStreamToMsgLevel(LogStream* ls, message::Level level);

/*!
 * \brief Adds the given GenericOutputStream to the given level.
 *
 * \param [in] ls pointer to the GenericOutputStream.
 * \param [in] level the level to log.
 * \pre ls != nullptr
 */
void addStreamToMsgLevel(GenericOutputStream* ls, message::Level level);

/*!
 * \brief Adds the given stream to all levels.
 *
 * \param [in] ls pointer to the log stream.
 * \pre ls != nullptr.
 */
void addStreamToAllMsgLevels(LogStream* ls);

/*!
 * \brief Adds the given GenericOutputStream to all levels.
 *
 * \param [in] ls pointer to the GenericOutputStream.
 * \pre ls != nullptr.
 */
void addStreamToAllMsgLevels(GenericOutputStream* ls);

/*!
* \brief Binds the given stream to the given tag.
*
* \param [in] ls pointer to the user-supplied LogStream object.
* \param [in] tag the tag that this stream will be associated with.
*
* \pre ls != nullptr.
*/
void addStreamToTag(LogStream* ls, const std::string& tag);

/*!
* \brief Binds the given GenericOutputStream to the given tag.
*
* \param [in] ls pointer to the user-supplied GenericOutputStream.
* \param [in] tag the tag that this stream will be associated with.
*
* \pre ls != nullptr.
*/
void addStreamToTag(GenericOutputStream* ls, const std::string& tag);

/*!
* \brief Binds the given stream to all the tags.
*
* \param [in] ls pointer to the user-supplied LogStream object.
*
* \pre ls != nullptr.
*/
void addStreamToAllTags(LogStream* ls);

/*!
* \brief Binds the given GenericOutputStream to all the tags.
*
* \param [in] ls pointer to the user-supplied GenericOutputStream.
*
* \pre ls != nullptr.
*/
void addStreamToAllTags(GenericOutputStream* ls);

/*!
* \brief Returns the number of streams for a given tag.
*        Returns 0 if the tag does not exist.
*
* \param [in] tag the tag in query.
*
* \return N the number of streams for the given tag.
* \post N >= 0
*/
int getNumStreamsWithTag(const std::string& tag);

/*!
 * \brief Logs the given message to all registered streams.
 *
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] filter_duplicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 */
void logMessage(message::Level level, const std::string& message, bool filter_duplicates = false);

/*!
 * \brief Logs the given message to all registered streams.
 *
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] tag user-supplied associated with this message.
 * \param [in] filter_duplicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 * /param [in] tag_stream_only optional parameter that indicates whether the
 * message will go only to streams bound to tagName. Default is false.
 */
void logMessage(message::Level level,
                const std::string& message,
                const std::string& tag,
                bool filter_duplicates = false,
                bool tag_stream_only = false);

/*!
 * \brief Logs the given message to all registered streams.
 *
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] fileName the name of the file this message is logged from.
 * \param [in] line the line number within the file this message is logged.
 * \param [in] filter_duplicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 */
void logMessage(message::Level level,
                const std::string& message,
                const std::string& fileName,
                int line,
                bool filter_duplicates = false);

/*!
 * \brief Logs the given message to all registered streams.
 *
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] tag user-supplied tag associated with the message.
 * \param [in] fileName the name of the file this message is logged from.
 * \param [in] line the line number within the file this message is logged.
 * \param [in] filter_duplicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 * /param [in] tag_stream_only optional parameter that indicates whether the
 * message will go only to streams bound to tagName. Default is false.
 */
void logMessage(message::Level level,
                const std::string& message,
                const std::string& tag,
                const std::string& fileName,
                int line,
                bool filter_duplicates = false,
                bool tag_stream_only = false);

/*!
 * \brief Convenience method to log an error message.
 *
 * \param [in] message user-supplied message.
 * \param [in] fileName the name of the file this message is logged from.
 * \param [in] line the line number within the file that the message is logged.
 */
void logErrorMessage(const std::string& message, const std::string& fileName, int line);

/*!
 * \brief Convenience method to log warning messages.
 *
 * \param [in] message user-supplied message.
 * \param [in] fileName the name of the file this message is logged from.
 * \param [in] line the line number within the file that the message is logged.
 */
void logWarningMessage(const std::string& message, const std::string& fileName, int line);

/*!
 * \brief For the current rank, outputs messages from all streams to the
 *        console
 *
 * \warning outputLocalMessages() is used before a rank aborts.
 *          flushStreams() is preferred over this function,
 *          as outputLocalMessages() may put LogStreams in an undesirable
 *          state. This call is not collective.
 */
void outputLocalMessages();

///@{
//! \name Collective Methods
//!
//! \attention These methods are collective operations.
//! All ranks in the user-supplied communicator must call the method
//! when used within an MPI distributed environment.
//!

/*!
 * \brief Calls abort function. Default is abort() or MPI_Abort() in a
 *        MPI distributed environment.
 *
 * \collective
 *
 */
void abort();

/*!
 * \brief Flushes all streams for all ranks.
 *
 * \collective
 * \see Logger::flushStreams.
 * \note When used within an MPI distributed environment, flushStreams is
 *  a collective operation. All ranks in the
 *  user-supplied communicator must call this method.
 */
void flushStreams();

/*!
 * \brief Pushes all streams.
 *
 * \collective
 * \see Logger::pushStreams.
 * \note When used within an MPI distributed environment, pushStreams is
 *  a collective operation. All ranks in the user-supplied communicator must
 *  call this method.
 */
void pushStreams();

/*!
 * \brief Finalizes the slic logging environment.
 *
 * \collective
 */
void finalize();

///@}

/*!
 * \brief Uses glibc's backtrace() functionality to return a stacktrace
 *
 * \returns a string corresponding to the stacktrace.
 */
std::string stacktrace();

} /* namespace slic */

} /* namespace axom */

#endif /* SLIC_HPP_ */
