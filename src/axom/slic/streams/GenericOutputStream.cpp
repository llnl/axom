// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/streams/GenericOutputStream.hpp"

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

namespace axom
{
namespace slic
{
GenericOutputStream::GenericOutputStream(std::ostream* os)
  : m_stream(os)
  , m_file_name()
  , m_opened(true)
  , m_isOstreamOwnedBySLIC(false)
{ }

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(const std::string& stream)
{
  if(stream == "cout")
  {
    m_stream = &std::cout;
    m_file_name = std::string();
    m_opened = true;
    m_isOstreamOwnedBySLIC = false;
  }
  else if(stream == "cerr")
  {
    m_stream = &std::cerr;
    m_file_name = std::string();
    m_opened = true;
    m_isOstreamOwnedBySLIC = false;
  }
  else
  {
    m_stream = new std::ostringstream();
    m_file_name = stream;
    m_opened = false;
    m_isOstreamOwnedBySLIC = true;
  }
}

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(std::ostream* os, const std::string& format)
  : m_stream(os)
  , m_file_name()
  , m_opened(true)
  , m_isOstreamOwnedBySLIC(false)
{
  this->setFormatString(format);
}

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(const std::string& stream, const std::string& format)
  : GenericOutputStream::GenericOutputStream(stream)
{
  // Fix newline and tab characters if needed
  std::string format_fixed = axom::utilities::string::replaceAllInstances(
    axom::utilities::string::replaceAllInstances(format, "\\n", "\n"),
    "\\t",
    "\t");
  this->setFormatString(format_fixed);
}

//------------------------------------------------------------------------------
GenericOutputStream::~GenericOutputStream()
{
  if(m_isOstreamOwnedBySLIC)
  {
    delete m_stream;
    m_stream = static_cast<std::ostream*>(nullptr);
  }
}

//------------------------------------------------------------------------------
void GenericOutputStream::append(message::Level msgLevel,
                                 const std::string& message,
                                 const std::string& tagName,
                                 const std::string& fileName,
                                 int line,
                                 bool AXOM_UNUSED_PARAM(filtered_duplicates),
                                 bool AXOM_UNUSED_PARAM(tag_stream_only))
{
  if(m_stream == nullptr)
  {
    std::cerr << "ERROR: NULL stream!\n";
    return;
  }

  (*m_stream) << this->getFormatedMessage(message::getLevelAsString(msgLevel),
                                          message,
                                          tagName,
                                          "",
                                          "",
                                          fileName,
                                          line);
}

//------------------------------------------------------------------------------
void GenericOutputStream::openBeforeFlush()
{
  if(m_isOstreamOwnedBySLIC && !m_opened)
  {
    std::ostringstream* oss = dynamic_cast<std::ostringstream*>(m_stream);
    if(oss != nullptr)
    {
      // Converting stream from ostringstream to ofstream and
      // writing ostringstream's string buffer to ofstream
      std::string buffer = oss->str();
      if(!buffer.empty())
      {
        delete m_stream;
        m_stream = new std::ofstream(m_file_name);
        (*m_stream) << buffer;
        m_opened = true;
      }
    }
  }
}

//------------------------------------------------------------------------------
void GenericOutputStream::outputLocal()
{
  openBeforeFlush();
  m_stream->flush();
}

//------------------------------------------------------------------------------
void GenericOutputStream::flush()
{
  openBeforeFlush();
  m_stream->flush();
}

} /* namespace slic */

} /* namespace axom */
