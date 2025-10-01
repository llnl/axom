// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/streams/LumberjackStream.hpp"
#include <mpi.h>

#include <vector>

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

#include "axom/lumberjack/BinaryTreeCommunicator.hpp"
#include "axom/lumberjack/Lumberjack.hpp"
#include "axom/lumberjack/TextTagCombiner.hpp"

namespace axom
{
namespace slic
{
//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(std::ostream* stream, MPI_Comm comm, int ranksLimit)
  : m_isLJOwnedBySLIC(false)
  , m_isOstreamOwnedBySLIC(false)
  , m_stream(stream)
  , m_file_name()
  , m_opened(true)
{
  this->initializeLumberjack(comm, ranksLimit);
}

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(std::ostream* stream,
                                   MPI_Comm comm,
                                   int ranksLimit,
                                   const std::string& format)
  : m_isLJOwnedBySLIC(false)
  , m_isOstreamOwnedBySLIC(false)
  , m_stream(stream)
  , m_file_name()
  , m_opened(true)
{
  this->initializeLumberjack(comm, ranksLimit);
  this->setFormatString(format);
}

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(std::ostream* stream, axom::lumberjack::Lumberjack* lj)
  : m_lj(lj)
  , m_isLJOwnedBySLIC(false)
  , m_isOstreamOwnedBySLIC(false)
  , m_stream(stream)
  , m_file_name()
  , m_opened(true)
{ }

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(std::ostream* stream,
                                   axom::lumberjack::Lumberjack* lj,
                                   const std::string& format)
  : m_lj(lj)
  , m_isLJOwnedBySLIC(false)
  , m_isOstreamOwnedBySLIC(false)
  , m_stream(stream)
  , m_file_name()
  , m_opened(true)
{
  this->setFormatString(format);
}

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(const std::string stream, MPI_Comm comm, int ranksLimit)
{
  this->initializeLumberjack(comm, ranksLimit);

  if(stream == "cout")
  {
    m_isOstreamOwnedBySLIC = false;
    m_stream = &std::cout;
    m_file_name = std::string();
    m_opened = true;
  }
  else if(stream == "cerr")
  {
    m_isOstreamOwnedBySLIC = false;
    m_stream = &std::cerr;
    m_file_name = std::string();
    m_opened = true;
  }
  else
  {
    m_isOstreamOwnedBySLIC = true;
    m_stream = new std::ofstream();
    m_file_name = stream;
    m_opened = false;
  }
}

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(const std::string stream,
                                   MPI_Comm comm,
                                   int ranksLimit,
                                   const std::string& format)
  : LumberjackStream::LumberjackStream(stream, comm, ranksLimit)
{
  // Fix newline and tab characters if needed
  std::string format_fixed = axom::utilities::string::replaceAllInstances(
    axom::utilities::string::replaceAllInstances(format, "\\n", "\n"),
    "\\t",
    "\t");
  this->setFormatString(format_fixed);
}

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(const std::string stream, axom::lumberjack::Lumberjack* lj)
{
  m_lj = lj;
  m_isLJOwnedBySLIC = false;

  if(stream == "cout")
  {
    m_isOstreamOwnedBySLIC = false;
    m_stream = &std::cout;
    m_file_name = std::string();
    m_opened = true;
  }
  else if(stream == "cerr")
  {
    m_isOstreamOwnedBySLIC = false;
    m_stream = &std::cerr;
    m_file_name = std::string();
    m_opened = true;
  }
  else
  {
    m_isOstreamOwnedBySLIC = true;
    m_stream = new std::ofstream();
    m_file_name = stream;
    m_opened = false;
  }
}

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(const std::string stream,
                                   axom::lumberjack::Lumberjack* lj,
                                   const std::string& format)
  : LumberjackStream::LumberjackStream(stream, lj)
{
  // Fix newline and tab characters if needed
  std::string format_fixed = axom::utilities::string::replaceAllInstances(
    axom::utilities::string::replaceAllInstances(format, "\\n", "\n"),
    "\\t",
    "\t");
  this->setFormatString(format_fixed);
}

//------------------------------------------------------------------------------
LumberjackStream::~LumberjackStream()
{
  if(m_isLJOwnedBySLIC)
  {
    this->finalizeLumberjack();
  }

  if(m_isOstreamOwnedBySLIC)
  {
    delete m_stream;
    m_stream = static_cast<std::ostream*>(nullptr);
  }
}

//------------------------------------------------------------------------------
void LumberjackStream::append(message::Level msgLevel,
                              const std::string& message,
                              const std::string& tagName,
                              const std::string& fileName,
                              int line,
                              bool AXOM_UNUSED_PARAM(filter_duplicates),
                              bool AXOM_UNUSED_PARAM(tag_stream_only))
{
  if(m_lj == nullptr)
  {
    std::cerr << "ERROR: NULL Lumberjack instance in LumberjackStream::append!\n";
    return;
  }

  m_lj->queueMessage(message, fileName, line, msgLevel, MPI_Wtime(), tagName);
}

//------------------------------------------------------------------------------
void LumberjackStream::outputLocal()
{
  if(m_lj == nullptr)
  {
    std::cerr << "ERROR: NULL Lumberjack instance in LumberjackStream::flush!\n";
    return;
  }

  //Non-collective write to console
  this->write(true);
}

//------------------------------------------------------------------------------
void LumberjackStream::flush()
{
  if(m_lj == nullptr)
  {
    std::cerr << "ERROR: NULL Lumberjack instance in LumberjackStream::flush!\n";
    return;
  }

  // Collective push of messages to output node followed by write to console
  m_lj->pushMessagesFully();
  this->write();
}

//------------------------------------------------------------------------------
void LumberjackStream::push()
{
  if(m_lj == nullptr)
  {
    std::cerr << "ERROR: NULL Lumberjack instance in LumberjackStream::push!\n";
    return;
  }

  m_lj->pushMessagesOnce();
}

//------------------------------------------------------------------------------
void LumberjackStream::write(bool local)
{
  if(m_lj == nullptr)
  {
    std::cerr << "ERROR: NULL Lumberjack instance in LumberjackStream::write!\n";
    return;
  }

  if(m_lj->isOutputNode() || local)
  {
    for(const auto* curr_message : m_lj->getMessages())
    {
      if(curr_message == nullptr)
      {
        continue;
      }

      if(m_isOstreamOwnedBySLIC && !m_opened)
      {
        std::ofstream* ofs = dynamic_cast<std::ofstream*>(m_stream);
        if(ofs != nullptr)
        {
          ofs->open(m_file_name);
          m_opened = true;
        }
      }

      (*m_stream) << this->getFormatedMessage(
        message::getLevelAsString(static_cast<message::Level>(curr_message->level())),
        curr_message->text(),
        curr_message->tag(),
        curr_message->stringOfRanks(),
        std::to_string(curr_message->count()),
        curr_message->fileName(),
        curr_message->lineNumber());
    }

    m_stream->flush();
    m_lj->clearMessages();
  }
}

//------------------------------------------------------------------------------
bool LumberjackStream::hasPendingMessages() { return m_lj->getMessages().size() > 0; }

//------------------------------------------------------------------------------
bool LumberjackStream::isUsingMPI() { return true; }

MPI_Comm LumberjackStream::comm() { return m_lj->getCommunicator()->comm(); }

//------------------------------------------------------------------------------
void LumberjackStream::initializeLumberjack(MPI_Comm comm, int ranksLimit)
{
  m_ljComm = new axom::lumberjack::BinaryTreeCommunicator;
  m_ljComm->initialize(comm, ranksLimit);
  m_lj = new axom::lumberjack::Lumberjack;
  m_lj->initialize(m_ljComm, ranksLimit);
  m_isLJOwnedBySLIC = true;
}

//------------------------------------------------------------------------------
void LumberjackStream::finalizeLumberjack()
{
  m_lj->finalize();
  m_ljComm->finalize();
  delete m_lj;
  delete m_ljComm;
  m_isLJOwnedBySLIC = false;
}

} /* namespace slic */
} /* namespace axom */
