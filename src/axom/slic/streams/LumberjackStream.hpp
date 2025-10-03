// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file LumberjackStream.hpp
 *
 */

#ifndef LUMBERJACKSTREAM_HPP_
#define LUMBERJACKSTREAM_HPP_

#include "axom/slic/core/LogStream.hpp"

#include "axom/core/Macros.hpp"

// C/C++ includes
#include <iostream>  // for std::ostream
#include <fstream>   // for ofstream

// MPI
#include <mpi.h>  // For MPI

// Forward declarations
namespace axom
{
namespace lumberjack
{
class Lumberjack;
class Communicator;

}  // namespace lumberjack
}  // namespace axom

namespace axom
{
namespace slic
{
/*!
 * \class LumberjackStream
 *
 * \brief A concrete instance of LogStream that utilizes Lumberjack to
 *  filter and pass messages between MPI nodes.
 *
 */
class LumberjackStream : public LogStream
{
public:
  /*!
   * \brief Constructs a LumberjackStream instance with the given stream,
   *  MPI communicator, and rank limit.
   * \param [in] stream pointer to a user-supplied ostream instance.
   * \param [in] comm MPI communicator
   * \param [in] ranksLimit limit on how many ranks are individually tracked per
   *             message
   * \pre stream != NULL
   */
  LumberjackStream(std::ostream* stream, MPI_Comm comm, int ranksLimit);

  /*!
   * \brief Constructs a LumberjackStream instance with the given stream,
   *  MPI communicator, rank limit, and message formatting.
   * \param [in] stream pointer to a user-supplied ostream instance.
   * \param [in] comm MPI communicator
   * \param [in] ranksLimit limit on how many ranks are individually tracked per
   *             message
   * \param [in] format the format string.
   * \pre stream != NULL
   * \see LogStream::setFormatString for the format string.
   */
  LumberjackStream(std::ostream* stream, MPI_Comm comm, int ranksLimit, const std::string& format);

  /*!
   * \brief Constructs a LumberjackStream instance with the given stream
   *        and Lumberjack communicator.
   * \param [in] stream pointer to a user-supplied ostream instance.
   * \param [in] lj Lumberjack communicator
   * \pre stream != NULL
   */
  LumberjackStream(std::ostream* stream, axom::lumberjack::Lumberjack* lj);

  /*!
   * \brief Constructs a LumberjackStream instance with the given stream
   *        Lumberjack communicator, and message formatting.
   * \param [in] stream pointer to a user-supplied ostream instance.
   * \param [in] lj Lumberjack communicator
   * \param [in] format the format string.
   * \pre stream != NULL
   */
  LumberjackStream(std::ostream* stream, axom::lumberjack::Lumberjack* lj, const std::string& format);

  /*!
   * \brief Constructs a LumberjackStream instance specified by the given
   *  string, MPI communicator, and rank limit.
   *  The string input determines the stream as follows:
   *   - "cout" makes std::cout the output stream
   *   - "cerr" makes std::cerr the output stream
   *   - Any other input will construct a std::ofstream associated with input
   * \param [in] stream the string to control type of stream created
   * \param [in] comm MPI communicator
   * \param [in] ranksLimit limit on how many ranks are individually tracked per
   *             message
   * \pre stream != NULL
   *
   * \note This constructor avoids creating an empty file if this
   *       LumberjackStream never flushes a message.
   */
  LumberjackStream(std::string stream, MPI_Comm comm, int ranksLimit);

  /*!
   * \brief Constructs a LumberjackStream instance specified by the given
   *  string, MPI communicator, rank limit, and message formatting.
   *  The string input determines the stream as follows:
   *   - "cout" makes std::cout the output stream
   *   - "cerr" makes std::cerr the output stream
   *   - Any other input will construct a std::ofstream associated with input
   * \param [in] stream the string to control type of stream created
   * \param [in] comm MPI communicator
   * \param [in] ranksLimit limit on how many ranks are individually tracked per
   *             message
   * \param [in] format the format string.
   * \pre stream != NULL
   * \see LogStream::setFormatString for the format string.
   *
   * \note This constructor avoids creating an empty file if this
   *       LumberjackStream never flushes a message.
   */
  LumberjackStream(std::string stream, MPI_Comm comm, int ranksLimit, const std::string& format);

  /*!
   * \brief Constructs a LumberjackStream instance specified by the given
   *  string and Lumberjack communicator.
   *  The string input determines the stream as follows:
   *   - "cout" makes std::cout the output stream
   *   - "cerr" makes std::cerr the output stream
   *   - Any other input will construct a std::ofstream associated with input
   * \param [in] stream the string to control type of stream created
   * \param [in] lj Lumberjack communicator
   * \pre stream != NULL
   *
   * \note This constructor avoids creating an empty file if this
   *       LumberjackStream never flushes a message.
   */
  LumberjackStream(std::string stream, axom::lumberjack::Lumberjack* lj);

  /*!
   * \brief Constructs a LumberjackStream instance specified by the given
   *  string, Lumberjack communicator, and message formatting.
   *  The string input determines the stream as follows:
   *   - "cout" makes std::cout the output stream
   *   - "cerr" makes std::cerr the output stream
   *   - Any other input will construct a std::ofstream associated with input
   * \param [in] stream the string to control type of stream created
   * \param [in] lj Lumberjack communicator
   * \param [in] format the format string.
   * \pre stream != NULL
   * \see LogStream::setFormatString for the format string.
   *
   * \note This constructor avoids creating an empty file if this
   *       LumberjackStream never flushes a message.
   */
  LumberjackStream(std::string stream, axom::lumberjack::Lumberjack* lj, const std::string& format);

  virtual ~LumberjackStream();

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
   * \brief Pushes the messages from the current rank directly to the
   *        console (non-collectively).
   *
   * \warning This method is being called before slic aborts.
   */
  virtual void outputLocal() override;

  /*!
   * \brief Pushes all messages to the output node according to Lumberjack's
   *  Communication scheme. Then writes it to the console.
   *
   * \collective
   * \note This method is a collective operation
   *  intended for a synchronization checkpoint.
   */
  virtual void flush() override;

  /*!
   * \brief Pushes all messages once to their parent node according to
   *  Lumberjack's Communication scheme.
   *
   * \collective
   * \note This method is a collective operation
   *  intended for a synchronization checkpoint.
   * \note This does not guarantee all messages have reached the output node.
   * \note This does not write out to the given stream.
   */
  virtual void push() override;

  /*!
   * \brief Writes the messages to the given stream that are at the output node
   *  or at the current node if local is true
   *
   *  param [in] local If true, writes out messages at the current node.
   *             If false, only writes out messages at the output node.
   *             Default is false.
   *
   *  It does not flush any messages and not all messages are guaranteed to be
   *  at the output node.
   */
  void write(bool local = false);

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

protected:
  void initializeLumberjack(MPI_Comm comm, int ranksLimit);
  void finalizeLumberjack();

  /// \name Protected Members
  /// @{

  axom::lumberjack::Lumberjack* m_lj;
  axom::lumberjack::Communicator* m_ljComm;
  bool m_isLJOwnedBySLIC;
  bool m_isOstreamOwnedBySLIC;
  std::ostream* m_stream;
  std::string m_file_name;
  bool m_opened;

  /// @}

private:
  /*!
   * \brief Default constructor. Made private to prevent applications from
   *  using it. Instead the constructor that passes the underlying Lumberjack
   * instance
   *  should be used.
   */
  LumberjackStream()
    : m_lj(static_cast<axom::lumberjack::Lumberjack*>(nullptr))
    , m_stream(static_cast<std::ostream*>(nullptr)) { };

  DISABLE_COPY_AND_ASSIGNMENT(LumberjackStream);
  DISABLE_MOVE_AND_ASSIGNMENT(LumberjackStream);
};

} /* namespace slic */
} /* namespace axom */

#endif /* LUMBERJACKSTREAM_HPP_ */
