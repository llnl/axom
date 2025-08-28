// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file LogStreamStatusMonitor.hpp
 *
 */

#ifndef LOGSTREAMSTATUS_MONITOR_HPP_
#define LOGSTREAMSTATUS_MONITOR_HPP_


#include <vector>
#include "axom/slic/core/LogStream.hpp"

namespace axom
{
namespace slic
{
/*!
 * \class LogStreamStatusMonitor
 *
 * \brief Add a description later
 *
 */
class LogStreamStatusMonitor
{
public:
  LogStreamStatusMonitor();
  virtual ~LogStreamStatusMonitor();

  /*!
   * \brief ADD DESCRIPTION HERE.  This is where we add to m_streamVec and combine its communicator with the current m_mpiComm.
   *
   * \param [in] ls pointer to the user-supplied LogStream object.
   *
   * \note Something something something add later
   */
  virtual void addStream(LogStream* ls);

  /*!
   * \brief ADD DESCRIPTION HERE.  This is the main point of this class!
   *
   * \note Something something something add later
   */
  virtual bool hasPendingMessages() const;

  /*!
   * \brief ADD DESCRIPTION HERE. Free MPI comm and that's it.
   *
   * \note Something something something add later
   */
  virtual void finalize();

protected:

  std::vector<LogStream*> m_streamVec;

};

} /* namespace slic */

} /* namespace axom */

#endif /* LOGSTREAMStatusMonitor_HPP_ */
