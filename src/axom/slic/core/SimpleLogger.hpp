// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *  \file SimpleLogger.hpp
 *
 *  \brief Header file containing definition of SimpleLogger class.
 *
 */

#ifndef SLIC_SIMPLELOGGER_HPP_
#define SLIC_SIMPLELOGGER_HPP_

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

// slic component headers
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

namespace axom
{
namespace slic
{
/*!
 * \class SimpleLogger
 *
 * \brief SimpleLogger is a simple wrapper around the initialization and
 * finalize operations of the slic::Logger class that is helpful for
 * unit tests and simple applications in axom.
 *
 * To use, create an instance of this class before tests are run. This
 * initializes the slic logger. When the object is destroyed (e.g., goes out
 * of scope), the slic logger is finalized. For example, when using gtest,
 * a simple main program can be used in each test source file to do this:
 *
 * \verbatim
 *
 *  int main(int argc, char* argv[])
 *  {
 *     int result = 0;
 *
 *     ::testing::InitGoogleTest(&argc, argv);
 *
 *     // create & initialize test logger, finalized when exiting main scope
 *     axom::slic::SimpleLogger logger;
 *
 *     result = RUN_ALL_TESTS();
 *
 *     return 0;
 *  }
 *
 * \endverbatim
 */
class SimpleLogger
{
public:
  /*!
   * \brief Constructor initializes slic logging environment.
   */
  explicit SimpleLogger(message::Level level = message::Debug)
  {
    initialize();
    setLoggingMsgLevel(level);

    // Formatting for warning, errors and fatal message
    std::string wefFormatStr = std::string("\n***********************************\n") +
      std::string("[<LEVEL> in line <LINE> of file <FILE>]\n") +
      std::string("MESSAGE=<MESSAGE>\n") + std::string("***********************************\n");

    // Simple formatting for debug and info messages
    std::string diFormatStr = "[<LEVEL>] <MESSAGE> \n";

    GenericOutputStream* wefStream = new GenericOutputStream(&std::cout, wefFormatStr);
    GenericOutputStream* diStream = new GenericOutputStream(&std::cout, diFormatStr);

    addStreamToMsgLevel(wefStream, message::Error);
    addStreamToMsgLevel(wefStream, message::Warning);
    addStreamToMsgLevel(diStream, message::Info);
    addStreamToMsgLevel(diStream, message::Debug);
  }

  /*!
   * \brief Destructor finalizes slic loging environment.
   */
  ~SimpleLogger()
  {
    if(isInitialized())
    {
      flushStreams();
      finalize();
    }
  }

private:
  DISABLE_COPY_AND_ASSIGNMENT(SimpleLogger);
  DISABLE_MOVE_AND_ASSIGNMENT(SimpleLogger);
};

}  // namespace slic
}  // namespace axom

#endif  // SLIC_SIMPLELOGGER_HPP_
