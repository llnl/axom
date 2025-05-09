// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file LogStream.cpp
 *
 */

#include "axom/slic/core/LogStream.hpp"

// C/C++ includes
#include <ctime>
#include <sstream>

namespace axom
{
namespace slic
{
//------------------------------------------------------------------------------
LogStream::LogStream()
  : m_formatString("*****\n[<LEVEL>]\n\n <MESSAGE> \n\n <FILE>\n<LINE>\n****\n")
{ }

//------------------------------------------------------------------------------
LogStream::~LogStream() { }

//------------------------------------------------------------------------------
void LogStream::replaceKey(std::string& msg,
                           const std::string& key,
                           const std::string& value,
                           std::size_t pos)
{
  if(pos == std::string::npos)
  {
    pos = msg.find(key);
  }

  if(pos != std::string::npos)
  {
    msg = msg.substr(0, pos) + value + msg.substr(pos + key.length(), msg.length() - 1);

  }  // END if
}

//------------------------------------------------------------------------------
std::string LogStream::getTimeStamp()
{
#ifdef WIN32
  #pragma warning(disable : 4996)  // _CRT_SECURE_NO_WARNINGS
#endif

  std::time_t t;
  std::time(&t);
  std::string timestamp(std::asctime(std::localtime(&t)));

  // Remove trailing newline added by previous line
  if(timestamp[timestamp.size() - 1] == '\n')
  {
    timestamp.erase(timestamp.size() - 1);
  }
  return timestamp;
}

//------------------------------------------------------------------------------
std::string LogStream::getFormatedMessage(const std::string& msgLevel,
                                          const std::string& message,
                                          const std::string& tagName,
                                          const std::string& rank,
                                          const std::string& rank_count,
                                          const std::string& fileName,
                                          int line)
{
  std::string msg = m_formatString;

  this->replaceKey(msg, "<LEVEL>", msgLevel);
  this->replaceKey(msg, "<MESSAGE>", message);
  this->replaceKey(msg, "<TAG>", tagName);
  this->replaceKey(msg, "<FILE>", fileName);
  this->replaceKey(msg, "<RANK>", rank);
  this->replaceKey(msg, "<RANK_COUNT>", rank_count);

  if(line != MSG_IGNORE_LINE)
  {
    std::ostringstream oss;
    oss << line;

    this->replaceKey(msg, "<LINE>", oss.str());
  }
  else
  {
    this->replaceKey(msg, "<LINE>", "");
  }

  std::size_t pos = msg.find("<TIMESTAMP>");
  if(pos != std::string::npos)
  {
    this->replaceKey(msg, "<TIMESTAMP>", this->getTimeStamp(), pos);
  }

  return (msg);
}

} /* namespace slic */

} /* namespace axom */
