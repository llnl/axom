// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/streams/SynchronizedStream.hpp"

#include <vector>

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

namespace axom
{
namespace slic
{
struct SynchronizedStream::MessageCache
{
  std::vector<std::string> messages;

  void printMessages(std::ostream* stream)
  {
    if(stream == nullptr)
    {
      std::cerr << "ERROR: cannot write to NULL stream!\n";
      return;
    }

    const unsigned N = messages.size();

    if(N == 0)
    {
      /* short-circuit */
      return;
    }

    for(unsigned i = 0; i < N; ++i)
    {
      (*stream) << messages[i];
    }  // END for all messages

    stream->flush();
    messages.clear();
  }
};

//------------------------------------------------------------------------------
SynchronizedStream::SynchronizedStream(std::ostream* stream, MPI_Comm comm)
  : m_comm(comm)
  , m_cache(new MessageCache())
  , m_stream(stream)
  , m_file_name()
  , m_isOstreamOwnedBySLIC(false)
  , m_opened(false)
{ }

//------------------------------------------------------------------------------
SynchronizedStream::SynchronizedStream(std::ostream* stream, MPI_Comm comm, const std::string& format)
  : m_comm(comm)
  , m_cache(new MessageCache)
  , m_stream(stream)
  , m_file_name()
  , m_isOstreamOwnedBySLIC(false)
  , m_opened(false)
{
  this->setFormatString(format);
}

//------------------------------------------------------------------------------
SynchronizedStream::SynchronizedStream(const std::string stream, MPI_Comm comm)
  : m_comm(comm)
  , m_cache(new MessageCache)
{
  if(stream == "cout")
  {
    m_stream = &std::cout;
    m_file_name = std::string();
    m_isOstreamOwnedBySLIC = false;
    m_opened = true;
  }
  else if(stream == "cerr")
  {
    m_stream = &std::cerr;
    m_file_name = std::string();
    m_isOstreamOwnedBySLIC = false;
    m_opened = true;
  }
  else
  {
    m_stream = new std::ofstream();
    m_file_name = stream;
    m_isOstreamOwnedBySLIC = true;
    m_opened = false;
  }
}

//------------------------------------------------------------------------------
SynchronizedStream::SynchronizedStream(const std::string stream,
                                       MPI_Comm comm,
                                       const std::string& format)
  : SynchronizedStream::SynchronizedStream(stream, comm)
{
  // Fix newline and tab characters if needed
  std::string format_fixed = axom::utilities::string::replaceAllInstances(
    axom::utilities::string::replaceAllInstances(format, "\\n", "\n"),
    "\\t",
    "\t");
  this->setFormatString(format_fixed);
}

//------------------------------------------------------------------------------
SynchronizedStream::~SynchronizedStream()
{
  delete m_cache;
  m_cache = static_cast<MessageCache*>(nullptr);

  if(m_isOstreamOwnedBySLIC)
  {
    delete m_stream;
    m_stream = static_cast<std::ostream*>(nullptr);
  }
}

//------------------------------------------------------------------------------
void SynchronizedStream::append(message::Level msgLevel,
                                const std::string& message,
                                const std::string& tagName,
                                const std::string& fileName,
                                int line,
                                bool AXOM_UNUSED_PARAM(filter_duplicates),
                                bool AXOM_UNUSED_PARAM(tag_stream_only))
{
  if(m_cache == nullptr)
  {
    std::cerr << "ERROR: NULL cache!\n";
    return;
  }

  int rank = -1;
  MPI_Comm_rank(m_comm, &rank);

  // STEP 1: cache formatted message
  m_cache->messages.push_back(this->getFormatedMessage(message::getLevelAsString(msgLevel),
                                                       message,
                                                       tagName,
                                                       std::to_string(rank),
                                                       "1",
                                                       fileName,
                                                       line));
}

//------------------------------------------------------------------------------
void SynchronizedStream::openBeforeFlush()
{
  if(m_isOstreamOwnedBySLIC && !m_opened && !m_cache->messages.empty())
  {
    std::ofstream* ofs = dynamic_cast<std::ofstream*>(m_stream);
    if(ofs != nullptr)
    {
      ofs->open(m_file_name);
      m_opened = true;
    }
  }
}

//------------------------------------------------------------------------------
void SynchronizedStream::outputLocal()
{
  if(m_cache == nullptr)
  {
    std::cerr << "ERROR: NULL cache!\n";
    return;
  }

  if(m_comm == MPI_COMM_NULL)
  {
    std::cerr << "ERROR: NULL communicator!\n";
    return;
  }

  openBeforeFlush();

  // print messages for this rank
  m_cache->printMessages(m_stream);
}

//------------------------------------------------------------------------------
void SynchronizedStream::flush()
{
  if(m_cache == nullptr)
  {
    std::cerr << "ERROR: NULL cache!\n";
    return;
  }

  if(m_comm == MPI_COMM_NULL)
  {
    std::cerr << "ERROR: NULL communicator!\n";
    return;
  }

  openBeforeFlush();

  // Collective flush
  int rank = -1;
  int nranks = 0;
  MPI_Comm_rank(m_comm, &rank);
  MPI_Comm_size(m_comm, &nranks);

  const int prevrank = rank - 1;
  const int nextrank = rank + 1;

  if(rank > 0)
  {
    // wait for signal from previous rank
    MPI_Recv(nullptr, 0, MPI_INT, prevrank, MPI_ANY_TAG, m_comm, MPI_STATUSES_IGNORE);
  }

  // print messages for this rank
  m_cache->printMessages(m_stream);

  // signal next rank
  if(nranks > 1 && nextrank < nranks)
  {
    MPI_Request null_request = MPI_REQUEST_NULL;
    MPI_Isend(nullptr, 0, MPI_INT, nextrank, 0, m_comm, &null_request);
    MPI_Request_free(&null_request);
  }
}

//------------------------------------------------------------------------------
bool SynchronizedStream::hasPendingMessages() { return m_cache->messages.size() > 0; }

//------------------------------------------------------------------------------
bool SynchronizedStream::isUsingMPI() { return true; }

//------------------------------------------------------------------------------
MPI_Comm SynchronizedStream::comm() { return m_comm; }

} /* namespace slic */
} /* namespace axom */
