// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file SynchronizedStream.hpp
 *
 */

#ifndef SYNCHRONIZEDSTREAM_HPP_
#define SYNCHRONIZEDSTREAM_HPP_

#include "axom/slic/core/LogStream.hpp"

#include "axom/core/Macros.hpp"

// C/C++ includes
#include <iostream>  // for std::ostream
#include <fstream>   // for ofstream

// MPI
#include <mpi.h>  // For MPI

namespace axom
{
namespace slic
{
/*!
 * \class SynchronizedStream
 *
 * \brief A concrete instance of LogStream that dumps messages to a C++
 *  std::ostream object.
 *
 * \note The intent of this class is to illustrate how to using the Logging
 *  facility within an MPI distributed environment and provide a utility that
 *  could be useful for debugging problems at small scales.
 *
 * \warning Do not use this for large-scale production runs.
 * \warning The intent of this class is to be used primarily with std::cout,
 *  std::cerr, etc. It is suggested that applications do not use this class
 *  with an std::ofstream object.
 */
class SynchronizedStream : public LogStream
{
public:
  /*!
   * \brief Constructs a SynchronizedStream instance with the given stream
   *  and MPI communicator.
   * \param [in] os pointer to a user-supplied ostream instance.
   * \param [in] comm MPI communicator
   * \pre stream != NULL
   */
  SynchronizedStream(std::ostream* stream, MPI_Comm comm);

  /*!
   * \brief Constructs a SynchronizedStream instance with the given stream,
   *  MPI communicator, and message formatting.
   * \param [in] os pointer to a user-supplied ostream instance.
   * \param [in] comm MPI communicator
   * \param [in] format the format string.
   * \pre stream != NULL
   * \see LogStream::setFormatString for the format string.
   */
  SynchronizedStream(std::ostream* stream, MPI_Comm comm, const std::string& format);

  /*!
   * \brief Constructs a SynchronizedStream instance specified by the given
   *  string, and MPI communicator.
   *  The string input determines the stream as follows:
   *   - "cout" makes std::cout the output stream
   *   - "cerr" makes std::cerr the output stream
   *   - Any other input will construct a std::ofstream associated with input
   * \param [in] stream the string to control type of stream created
   * \param [in] comm MPI communicator
   * \pre stream != NULL
   * \see LogStream::setFormatString for the format string.
   *
   * \note This constructor avoids creating an empty file if this
   *       SynchronizedStream never flushes a message.
   */
  SynchronizedStream(std::string stream, MPI_Comm comm);

  /*!
   * \brief Constructs a SynchronizedStream instance specified by the given
   *  string, MPI communicator, and message formatting.
   *  The string input determines the stream as follows:
   *   - "cout" makes std::cout the output stream
   *   - "cerr" makes std::cerr the output stream
   *   - Any other input will construct a std::ofstream associated with input
   * \param [in] stream the string to control type of stream created
   * \param [in] comm MPI communicator
   * \param [in] format the format string.
   * \pre stream != NULL
   * \see LogStream::setFormatString for the format string.
   *
   * \note This constructor avoids creating an empty file if this
   *       SynchronizedStream never flushes a message.
   */
  SynchronizedStream(std::string stream, MPI_Comm comm, const std::string& format);

  virtual ~SynchronizedStream();

  /*!
   * \brief Appends the given message to the stream.
   *
   * \param [in] msgLevel the level of the message.
   * \param [in] message the user-supplied message.
   * \param [in] tagName user-supplied tag to associate with the given message.
   * \param [in] fileName the file where this message is appended
   * \param [in] line the line within the file at which the message is appended.
   * \param [in] filter_duplicates optional parameter that indicates whether
   * duplicate messages resulting from running in parallel will be filtered out.
   * /param [in] tag_stream_only optional parameter that indicates whether the
   * message will go only to streams bound to tagName.
   *
   * \note This method doesn't put anything to the console. Instead the
   *  messages are cached locally to each ranks and are dumped to the console
   *  in rank order when flush is called.
   */
  virtual void append(message::Level msgLevel,
                      const std::string& message,
                      const std::string& tagName,
                      const std::string& fileName,
                      int line,
                      bool filter_duplicates,
                      bool tag_stream_only) override;

  /*!
   * \brief Dumps the messages from the current rank directly to the
   *        console (non-collectively).
   *
   * \warning This method is being called before slic aborts.
   */
  virtual void outputLocal() override;

  /*!
   * \brief Dumps the messages to the console in rank-order for all ranks.
   *
   * \collective
   * \note This method is a collective operation
   *  intended for a synchronization checkpoint.
   */
  virtual void flush() override;

  /*!
   * \brief Tests whether there are any pending messages that need to be flushed
   *
   * \return Returns true if there are pending messages that need to be flushed
   */
  virtual bool hasPendingMessages() override;

  /*!
   * \brief Tests whether this class relies on MPI
   *
   * \return Returns true if this class relies on MPI
   */
  virtual bool isUsingMPI() override;

  /*!
   * \brief Get the communicator
   *
   * \return Returns the communicator if it exists, or MPI_COMM_NULL otherwise
   */
  virtual MPI_Comm comm() override;

private:
  /// Forward declarations
  struct MessageCache;

  /// \name Private Members
  /// @{

  MPI_Comm m_comm;
  MessageCache* m_cache;
  std::ostream* m_stream;
  std::string m_file_name;
  bool m_isOstreamOwnedBySLIC;
  bool m_opened;
  /// @}

  /*!
   * \brief Opens a file before flushing stream when SynchronizedStream
   *        has ownership of ostream to a file (std::string constructor)
   */
  void openBeforeFlush();

  /*!
   * \brief Default constructor. Made private to prevent applications from
   *  using it. Instead the constructor that passes the underlying MPI comm
   *  should be used.
   */
  SynchronizedStream()
    : m_comm(MPI_COMM_NULL)
    , m_cache(static_cast<MessageCache*>(nullptr))
    , m_stream(static_cast<std::ostream*>(nullptr))
    , m_file_name()
    , m_opened(false) { };

  DISABLE_COPY_AND_ASSIGNMENT(SynchronizedStream);
  DISABLE_MOVE_AND_ASSIGNMENT(SynchronizedStream);
};

} /* namespace slic */
} /* namespace axom */

#endif /* SYNCHRONIZEDSTREAM_HPP_ */
