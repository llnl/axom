// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file GenericOutputStream.hpp
 *
 */

#ifndef GENERICOUTPUTSTREAM_HPP_
#define GENERICOUTPUTSTREAM_HPP_

#include "axom/slic/core/LogStream.hpp"

#include "axom/core/Macros.hpp"

// C/C++ includes
#include <iostream>  // for ostream
#include <fstream>   // for ofstream

namespace axom
{
namespace slic
{
/*!
 * \class GenericOutputStream
 *
 * \brief Concrete instance of LogStream, which implements functionality for
 *  logging messages to a C++ ostream object, e.g., std::cout, std::cerr, a
 *  std::ofstream, std::ostringstream, etc.
 *
 * \see LogStream Logger
 */
class GenericOutputStream : public LogStream
{
public:
  /*!
   * \brief Constructs a GenericOutpuStream instance with the given stream.
   * \param [in] os pointer to a user-supplied ostream instance.
   * \pre os != NULL
   */
  GenericOutputStream(std::ostream* os);

  /*!
   * \brief Constructs a GenericOutputStream instance specified by the given
   *  string. The string input determines the stream as follows:
   *   - "cout" makes std::cout the output stream
   *   - "cerr" makes std::cerr the output stream
   *   - Any other input will construct a std::ofstream associated with input
   * \param [in] stream the string to control type of stream created
   *
   * \note This constructor avoids creating an empty file if this
   *       GenericOutputStream never flushes a message.
   */
  GenericOutputStream(const std::string& stream);

  /*!
   * \brief Constructs a GenericOutputStream instance with the given stream and
   *  message formatting.
   * \param [in] os pointer to a user-supplied ostream instance.
   * \param [in] format the format string.
   * \pre os != NULL
   * \see LogStream::setFormatString for the format string.
   */
  GenericOutputStream(std::ostream* os, const std::string& format);

  /*!
   * \brief Constructs a GenericOutputStream instance specified by the given
   *  string and message formatting.
   *  The string input determines the stream as follows:
   *   - "cout" makes std::cout the output stream
   *   - "cerr" makes std::cerr the output stream
   *   - Any other input will construct a std::ofstream associated with input
   * \param [in] stream the string to control type of stream created
   * \param [in] format the format string.
   * \see LogStream::setFormatString for the format string.
   *
   * \note This constructor avoids creating an empty file if this
   *       GenericOutputStream never flushes a message.
   */
  GenericOutputStream(const std::string& stream, const std::string& format);

  /*!
   * \brief Destructor.
   */
  virtual ~GenericOutputStream();

  /// \see LogStream::append
  virtual void append(message::Level msgLevel,
                      const std::string& message,
                      const std::string& tagName,
                      const std::string& fileName,
                      int line,
                      bool filter_duplicates,
                      bool tag_stream_only) override;

  /*!
   * \brief Outputs the log stream to the console.
   */
  virtual void outputLocal() override;

  /*!
   * \brief Flushes the log stream.
   */
  virtual void flush() override;

private:
  std::ostream* m_stream;
  std::string m_file_name;
  bool m_opened;
  bool m_isOstreamOwnedBySLIC;

  /*!
   * \brief Default constructor.
   * \note Made private to prevent applications from using it.
   */
  GenericOutputStream()
    : m_stream(static_cast<std::ostream*>(nullptr))
    , m_file_name()
    , m_opened(false)
    , m_isOstreamOwnedBySLIC(false) { };

  /*!
   * \brief Opens a file before flushing stream when GenericOutputStream
   *        has ownership of ostream to a file (std::string constructor)
   */
  void openBeforeFlush();

  DISABLE_COPY_AND_ASSIGNMENT(GenericOutputStream);
  DISABLE_MOVE_AND_ASSIGNMENT(GenericOutputStream);
};

} /* namespace slic */

} /* namespace axom */

#endif /* GENERICOUTPUTSTREAM_HPP_ */
