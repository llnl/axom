// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/core/Logger.hpp"
#include "axom/slic/core/LogStream.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/slic/core/LogStreamStatusMonitor.hpp"

// C/C++ includes
#include <iostream>

namespace axom
{
namespace slic
{
using Loggermap = std::map<std::string, Logger*>;

//------------------------------------------------------------------------------
// This is a singleton, scope-limited to this file.
Loggermap& getLoggers()
{
  static Loggermap s_loggers;
  return s_loggers;
}

//------------------------------------------------------------------------------
// This is a singleton, scope-limited to this file.
Logger*& getLogger()
{
  static Logger* s_Logger = nullptr;
  return s_Logger;
}

//------------------------------------------------------------------------------
// This is a singleton, scope-limited to this file.
LogStreamStatusMonitor& getLogStreamStatusMonitor()
{
  static LogStreamStatusMonitor s_logStreamStatusMonitor;
  return s_logStreamStatusMonitor;
}

//------------------------------------------------------------------------------
Logger::Logger()
  : m_abortOnError(true)
  , m_abortOnWarning(false)
  , m_abortFunction(axom::utilities::processAbort)
{
  // by default, all message streams are disabled
  for(int i = 0; i < message::Num_Levels; ++i)
  {
    m_isEnabled[i] = false;
  }
}

//------------------------------------------------------------------------------
Logger::Logger(const std::string& name)
  : m_name(name)
  , m_abortOnError(true)
  , m_abortOnWarning(false)
  , m_abortFunction(axom::utilities::processAbort)
{
  // by default, all message streams are disabled
  for(int i = 0; i < message::Num_Levels; ++i)
  {
    m_isEnabled[i] = false;
  }
}

//------------------------------------------------------------------------------
Logger::~Logger()
{
  for(auto& kv : m_streamObjectsManager)
  {
    delete kv.second;
  }

  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    m_logStreams[level].clear();
  }

  m_taggedStreams.clear();
}

//------------------------------------------------------------------------------
void Logger::abort()
{
  outputLocalMessages();
  m_abortFunction();
}

//------------------------------------------------------------------------------
void Logger::setAbortFunction(AbortFunctionPtr abort_func)
{
  if(abort_func == nullptr)
  {
    std::cerr << "WARNING: slic::Logger::setAbortFunction() -- supplied abort "
                 "function is invalid!\n";
    return;
  }

  m_abortFunction = abort_func;
}

//------------------------------------------------------------------------------
message::Level Logger::getLoggingMsgLevel()
{
  int lev = 0;
  const int maxLevels = static_cast<int>(message::Num_Levels);

  // find first disabled level
  while(lev < maxLevels && m_isEnabled[lev])
  {
    ++lev;
  }

  // Subtract one level (with clamp) since we overshot by one
  lev = axom::utilities::clampVal(lev - 1, 0, maxLevels - 1);

  return static_cast<message::Level>(lev);
}

//------------------------------------------------------------------------------
void Logger::setLoggingMsgLevel(message::Level level)
{
  for(int i = 0; i < message::Num_Levels; ++i)
  {
    m_isEnabled[i] = (i <= level) ? true : false;
  }
}

//------------------------------------------------------------------------------
void Logger::addStreamToMsgLevel(LogStream* ls, message::Level level, bool pass_ownership)
{
  if(ls == nullptr)
  {
    std::cerr << "WARNING: supplied log stream is NULL!\n";
    return;
  }

  m_logStreams[level].push_back(ls);

  if(pass_ownership)
  {
    m_streamObjectsManager[ls] = ls;
  }

  getLogStreamStatusMonitor().addStream(ls);
}

//------------------------------------------------------------------------------
void Logger::addStreamToTag(LogStream* ls, const std::string& tag, bool pass_ownership)
{
  if(ls == nullptr)
  {
    std::cerr << "WARNING: supplied log stream is NULL!\n";
    return;
  }

  if(m_taggedStreams.find(tag) == m_taggedStreams.end())
  {
    m_taggedStreams[tag] = std::vector<LogStream*> {ls};
  }
  else
  {
    m_taggedStreams[tag].push_back(ls);
  }

  if(pass_ownership)
  {
    m_streamObjectsManager[ls] = ls;
  }

  getLogStreamStatusMonitor().addStream(ls);
}

//------------------------------------------------------------------------------
void Logger::addStreamToAllMsgLevels(LogStream* ls, bool pass_ownership)
{
  if(ls == nullptr)
  {
    std::cerr << "WARNING: supplied log stream is NULL!\n";
    return;
  }

  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    this->addStreamToMsgLevel(ls, static_cast<message::Level>(level), pass_ownership);
  }
}

//------------------------------------------------------------------------------
void Logger::addStreamToAllTags(LogStream* ls, bool pass_ownership)
{
  if(ls == nullptr)
  {
    std::cerr << "WARNING: supplied log stream is NULL!\n";
    return;
  }

  if(m_taggedStreams.empty())
  {
    std::cerr << "WARNING: no tags are available!\n";
    if(pass_ownership)
    {
      m_streamObjectsManager[ls] = ls;
    }

    return;
  }

  for(auto& kv : m_taggedStreams)
  {
    this->addStreamToTag(ls, kv.first, pass_ownership);
  }
}

//------------------------------------------------------------------------------
int Logger::getNumStreamsAtMsgLevel(message::Level level)
{
  return static_cast<int>(m_logStreams[level].size());
}

//------------------------------------------------------------------------------
int Logger::getNumStreamsWithTag(const std::string& tag)
{
  if(m_taggedStreams.find(tag) == m_taggedStreams.end())
  {
    return 0;
  }

  return static_cast<int>(m_taggedStreams[tag].size());
}

//------------------------------------------------------------------------------
LogStream* Logger::getStream(message::Level level, int i)
{
  if(i < 0 || i >= static_cast<int>(m_logStreams[level].size()))
  {
    std::cerr << "ERROR: stream index is out-of-bounds!\n";
    return nullptr;
  }

  return m_logStreams[level][i];
}

//------------------------------------------------------------------------------
LogStream* Logger::getStream(const std::string& tag, int i)
{
  if(m_taggedStreams.find(tag) == m_taggedStreams.end())
  {
    std::cerr << "ERROR: tag does not exist!\n";
    return nullptr;
  }

  if(i < 0 || i >= static_cast<int>(m_taggedStreams[tag].size()))
  {
    std::cerr << "ERROR: stream index is out-of-bounds!\n";
    return nullptr;
  }

  return m_taggedStreams[tag][i];
}

//------------------------------------------------------------------------------
void Logger::logMessage(message::Level level, const std::string& message, bool filter_duplicates)
{
  this->logMessage(level,
                   message,
                   MSG_IGNORE_TAG,
                   MSG_IGNORE_FILE,
                   MSG_IGNORE_LINE,
                   filter_duplicates,
                   false);
}

//------------------------------------------------------------------------------
void Logger::logMessage(message::Level level,
                        const std::string& message,
                        const std::string& tagName,
                        bool filter_duplicates,
                        bool tag_stream_only)
{
  this->logMessage(level,
                   message,
                   tagName,
                   MSG_IGNORE_FILE,
                   MSG_IGNORE_LINE,
                   filter_duplicates,
                   tag_stream_only);
}

//------------------------------------------------------------------------------
void Logger::logMessage(message::Level level,
                        const std::string& message,
                        const std::string& fileName,
                        int line,
                        bool filter_duplicates)
{
  this->logMessage(level, message, MSG_IGNORE_TAG, fileName, line, filter_duplicates, false);
}

//------------------------------------------------------------------------------
void Logger::logMessage(message::Level level,
                        const std::string& message,
                        const std::string& tagName,
                        const std::string& fileName,
                        int line,
                        bool filter_duplicates,
                        bool tag_stream_only)
{
  if(m_isEnabled[level] == false && tag_stream_only == false)
  {
    return;
  }

  if(tag_stream_only == true && tagName == MSG_IGNORE_TAG)
  {
    std::cerr << "ERROR: message for tagged streams does not have a tag!\n";
    return;
  }

  if(tag_stream_only == true && m_taggedStreams.find(tagName) == m_taggedStreams.end())
  {
    std::cerr << "ERROR: tag does not exist!\n";
    return;
  }

  // Message for message levels
  if(tag_stream_only == false)
  {
    unsigned nstreams = static_cast<unsigned>(m_logStreams[level].size());
    for(unsigned istream = 0; istream < nstreams; ++istream)
    {
      m_logStreams[level][istream]
        ->append(level, message, tagName, fileName, line, filter_duplicates, tag_stream_only);
    }
  }

  // Message for tagged streams
  else
  {
    for(unsigned int i = 0; i < m_taggedStreams[tagName].size(); i++)
    {
      m_taggedStreams[tagName][i]
        ->append(level, message, tagName, fileName, line, filter_duplicates, tag_stream_only);
    }
  }
}

//------------------------------------------------------------------------------
void Logger::outputLocalMessages()
{
  //Output for all message levels
  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    unsigned nstreams = static_cast<unsigned>(m_logStreams[level].size());
    for(unsigned istream = 0; istream < nstreams; ++istream)
    {
      m_logStreams[level][istream]->outputLocal();
    }  // END for all streams

  }  // END for all levels

  // Output for all tagged streams
  std::map<std::string, std::vector<LogStream*>>::iterator it;

  for(it = m_taggedStreams.begin(); it != m_taggedStreams.end(); it++)
  {
    for(unsigned int i = 0; i < it->second.size(); i++)
    {
      it->second[i]->outputLocal();
    }
  }
}

//------------------------------------------------------------------------------
void Logger::flushStreams()
{
  /*
    check if any MPI-based stream has pending messages on any rank.
    This is needed to avoid unecessary pushes or flushes when there are no
    pending messages.
   */
  const bool pendingMessages = hasPendingMessages();

  //Flush for all message levels
  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    unsigned nstreams = static_cast<unsigned>(m_logStreams[level].size());
    for(unsigned istream = 0; istream < nstreams; ++istream)
    {
      const bool streamUsesMPI = m_logStreams[level][istream]->isUsingMPI();
      if(shouldPushMessages(pendingMessages, streamUsesMPI))
      {
        m_logStreams[level][istream]->flush();
      }

    }  // END for all streams

  }  // END for all levels

  // Flush for all tagged streams
  std::map<std::string, std::vector<LogStream*>>::iterator it;

  for(it = m_taggedStreams.begin(); it != m_taggedStreams.end(); it++)
  {
    for(unsigned int i = 0; i < it->second.size(); i++)
    {
      const bool streamUsesMPI = it->second[i]->isUsingMPI();
      if(shouldPushMessages(pendingMessages, streamUsesMPI))
      {
        it->second[i]->flush();
      }
    }
  }
}

//------------------------------------------------------------------------------
void Logger::pushStreams()
{
  /*
    check if any MPI-based stream has pending messages on any rank.
    This is needed to avoid unecessary pushes or flushes when there are no
    pending messages.
   */
  const bool pendingMessages = hasPendingMessages();
  //Push for all message levels
  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    unsigned nstreams = static_cast<unsigned>(m_logStreams[level].size());
    for(unsigned istream = 0; istream < nstreams; ++istream)
    {
      const bool streamUsesMPI = m_logStreams[level][istream]->isUsingMPI();
      if(shouldPushMessages(pendingMessages, streamUsesMPI))
      {
        m_logStreams[level][istream]->push();
      }

    }  // END for all streams

  }  // END for all levels

  // Push for all tagged streams
  std::map<std::string, std::vector<LogStream*>>::iterator it;

  for(it = m_taggedStreams.begin(); it != m_taggedStreams.end(); it++)
  {
    for(unsigned int i = 0; i < it->second.size(); i++)
    {
      const bool streamUsesMPI = it->second[i]->isUsingMPI();
      if(shouldPushMessages(pendingMessages, streamUsesMPI))
      {
        it->second[i]->push();
      }
    }
  }
}

//------------------------------------------------------------------------------
bool Logger::hasPendingMessages() { return getLogStreamStatusMonitor().hasPendingMessages(); }

//------------------------------------------------------------------------------
//                Static Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void Logger::initialize()
{
  Logger::createLogger("root");
  Logger::activateLogger("root");
}

//------------------------------------------------------------------------------
bool Logger::createLogger(const std::string& name, char imask)
{
  Loggermap& loggers = getLoggers();
  if(loggers.find(name) != loggers.end())
  {
    std::cerr << "ERROR: " << name << " logger is duplicated!\n";
    return false;
  }

  loggers[name] = new Logger(name);

  if(imask == inherit::nothing)
  {
    /* short-circuit */
    return true;
  }  // END if inherit nothing

  Logger* rootLogger = Logger::getRootLogger();
  if(rootLogger == nullptr)
  {
    std::cerr << "ERROR: no root logger found!\n";
    return false;
  }

  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    message::Level current_level = static_cast<message::Level>(level);

    int nstreams = rootLogger->getNumStreamsAtMsgLevel(current_level);
    if(nstreams == 0)
    {
      continue;
    }

    if(imask & inherit::masks[level])
    {
      for(int istream = 0; istream < nstreams; ++istream)
      {
        loggers[name]->addStreamToMsgLevel(rootLogger->getStream(current_level, istream),
                                           current_level,
                                           /* pass_ownership */ false);

      }  // END for all streams at this level

    }  // END if inherit streams at the given level

  }  // END for all message levels

  return true;
}

//------------------------------------------------------------------------------
bool Logger::activateLogger(const std::string& name)
{
  Loggermap& loggers = getLoggers();
  if(loggers.find(name) != loggers.end())
  {
    getLogger() = loggers[name];
    return true;
  }

  return false;
}

//------------------------------------------------------------------------------
void Logger::finalize()
{
  Loggermap& loggers = getLoggers();
  for(auto& kv : loggers)
  {
    kv.second->flushStreams();
  }

  for(auto& kv : loggers)
  {
    delete kv.second;
  }
  loggers.clear();

  getLogger() = nullptr;

  LogStreamStatusMonitor& logStreamStatusMonitor = getLogStreamStatusMonitor();

  logStreamStatusMonitor.finalize();
}

//------------------------------------------------------------------------------
std::string Logger::getActiveLoggerName() { return getLogger()->getName(); }

//------------------------------------------------------------------------------
Logger* Logger::getActiveLogger() { return getLogger(); }

//------------------------------------------------------------------------------
Logger* Logger::getRootLogger()
{
  Loggermap& loggers = getLoggers();
  if(loggers.find("root") == loggers.end())
  {
    // no root logger
    return nullptr;
  }

  return (loggers["root"]);
}

bool Logger::shouldPushMessages(const bool hasPendingMessages, const bool streamUsesMPI) const
{
  return (!streamUsesMPI || hasPendingMessages);
}

} /* namespace slic */

} /* namespace axom */
