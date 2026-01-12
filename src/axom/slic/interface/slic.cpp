// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/internal/stacktrace.hpp"

#include <sstream>  // for std::ostringstream

namespace axom
{
namespace slic
{
//------------------------------------------------------------------------------
// Initialize static variables for controlling runtime behavior of asserts and
// error macros.
//------------------------------------------------------------------------------
bool debug::checksAreErrors = false;

static bool s_is_root = true;

void initialize(bool is_root)
{
  axom::slic::s_is_root = is_root;
  Logger::initialize();
}

//------------------------------------------------------------------------------
bool isInitialized() { return (Logger::getActiveLogger() != nullptr); }

//------------------------------------------------------------------------------
bool isRoot() { return axom::slic::s_is_root; }

//------------------------------------------------------------------------------
void setIsRoot(bool is_root) { axom::slic::s_is_root = is_root; }

//------------------------------------------------------------------------------
void ensureInitialized()
{
  if(!isInitialized())
  {
    initialize();
    disableAbortOnError();
    disableAbortOnWarning();
    setLoggingMsgLevel(slic::message::Debug);
    std::string format = std::string("<TIMESTAMP>\n") + std::string("[<LEVEL>]: <MESSAGE> \n") +
      std::string("FILE=<FILE>\n") + std::string("LINE=<LINE>\n\n");
    addStreamToAllMsgLevels(new slic::GenericOutputStream(&std::cout, format));

    logMessage(message::Warning,
               "slic::initialize() must be called before any other calls to SLIC\n."
               "The SLIC library called slic::initialize() for you and set up a minimal "
               "configuration\nto allow log messages to print.\n"
               "Please call slic::initialize() near the beginning of the code\n"
               "to fix this error and get rid of this message.\n"
               "Please call slic::finalize() after all other calls to SLIC.\n");
  }
}

//------------------------------------------------------------------------------
void createLogger(const std::string& name, char imask)
{
  ensureInitialized();
  Logger::createLogger(name, imask);
}

//------------------------------------------------------------------------------
bool activateLogger(const std::string& name)
{
  ensureInitialized();
  return Logger::activateLogger(name);
}

//------------------------------------------------------------------------------
std::string getActiveLoggerName()
{
  ensureInitialized();
  return (Logger::getActiveLoggerName());
}

//------------------------------------------------------------------------------
message::Level getLoggingMsgLevel()
{
  ensureInitialized();
  return Logger::getActiveLogger()->getLoggingMsgLevel();
}

//------------------------------------------------------------------------------
void abort()
{
  ensureInitialized();
  Logger::getActiveLogger()->abort();
}

//------------------------------------------------------------------------------
void setLoggingMsgLevel(message::Level level)
{
  ensureInitialized();
  Logger::getActiveLogger()->setLoggingMsgLevel(level);
}

//------------------------------------------------------------------------------
void setAbortOnError(bool status)
{
  ensureInitialized();

  Logger::getActiveLogger()->setAbortOnError(status);
}

//------------------------------------------------------------------------------
void enableAbortOnError() { setAbortOnError(true); }

//------------------------------------------------------------------------------
void disableAbortOnError() { setAbortOnError(false); }

//------------------------------------------------------------------------------
bool isAbortOnErrorsEnabled()
{
  ensureInitialized();
  return (Logger::getActiveLogger()->isAbortOnErrorsEnabled());
}

//------------------------------------------------------------------------------
void setAbortOnWarning(bool status)
{
  ensureInitialized();
  Logger::getActiveLogger()->setAbortOnWarning(status);
}

//------------------------------------------------------------------------------
void enableAbortOnWarning() { setAbortOnWarning(true); }

//------------------------------------------------------------------------------
void disableAbortOnWarning() { setAbortOnWarning(false); }

//------------------------------------------------------------------------------
bool isAbortOnWarningsEnabled()
{
  ensureInitialized();
  return (Logger::getActiveLogger()->isAbortOnWarningsEnabled());
}

//------------------------------------------------------------------------------
void setAbortFunction(AbortFunctionPtr abort_func)
{
  ensureInitialized();
  Logger::getActiveLogger()->setAbortFunction(abort_func);
}

//------------------------------------------------------------------------------
AbortFunctionPtr getAbortFunction()
{
  ensureInitialized();
  return Logger::getActiveLogger()->getAbortFunction();
}

//------------------------------------------------------------------------------
void addStreamToMsgLevel(LogStream* ls, message::Level level)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToMsgLevel(ls, level);
}

//------------------------------------------------------------------------------
void addStreamToMsgLevel(GenericOutputStream* ls, message::Level level)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToMsgLevel(ls, level);
}

//------------------------------------------------------------------------------
void addStreamToAllMsgLevels(LogStream* ls)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToAllMsgLevels(ls);
}

//------------------------------------------------------------------------------
void addStreamToAllMsgLevels(GenericOutputStream* ls)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToAllMsgLevels(ls);
}

//------------------------------------------------------------------------------
void addStreamToTag(LogStream* ls, const std::string& tag)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToTag(ls, tag);
}

//------------------------------------------------------------------------------
void addStreamToTag(GenericOutputStream* ls, const std::string& tag)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToTag(ls, tag);
}

//------------------------------------------------------------------------------
void addStreamToAllTags(LogStream* ls)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToAllTags(ls);
}

//------------------------------------------------------------------------------
void addStreamToAllTags(GenericOutputStream* ls)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToAllTags(ls);
}

//------------------------------------------------------------------------------
int getNumStreamsWithTag(const std::string& tag)
{
  ensureInitialized();
  return Logger::getActiveLogger()->getNumStreamsWithTag(tag);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level, const std::string& message, bool filter_duplicates)
{
  ensureInitialized();
  Logger::getActiveLogger()->logMessage(level, message, filter_duplicates);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                const std::string& tag,
                bool filter_duplicates,
                bool tag_stream_only)
{
  ensureInitialized();
  Logger::getActiveLogger()->logMessage(level, message, tag, filter_duplicates, tag_stream_only);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                const std::string& fileName,
                int line,
                bool filter_duplicates)
{
  ensureInitialized();
  Logger::getActiveLogger()->logMessage(level, message, fileName, line, filter_duplicates);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                const std::string& tag,
                const std::string& fileName,
                int line,
                bool filter_duplicates,
                bool tag_stream_only)
{
  ensureInitialized();
  Logger::getActiveLogger()
    ->logMessage(level, message, tag, fileName, line, filter_duplicates, tag_stream_only);
}

//------------------------------------------------------------------------------
void logErrorMessage(const std::string& message, const std::string& fileName, int line)
{
  std::ostringstream oss;
  oss << message << slic::internal::stacktrace();

  slic::logMessage(message::Error, oss.str(), fileName, line);
}

//------------------------------------------------------------------------------
void logWarningMessage(const std::string& message, const std::string& fileName, int line)
{
  slic::logMessage(message::Warning, message, fileName, line);
}

//------------------------------------------------------------------------------
void outputLocalMessages()
{
  ensureInitialized();
  Logger::getActiveLogger()->outputLocalMessages();
}

//------------------------------------------------------------------------------
void flushStreams()
{
  ensureInitialized();
  Logger::getActiveLogger()->flushStreams();
}

//------------------------------------------------------------------------------
void pushStreams()
{
  ensureInitialized();
  Logger::getActiveLogger()->pushStreams();
}

//------------------------------------------------------------------------------
bool hasPendingMessages() { return Logger::getActiveLogger()->hasPendingMessages(); }

//------------------------------------------------------------------------------
void finalize() { Logger::finalize(); }

} /* namespace slic */

} /* namespace axom */
