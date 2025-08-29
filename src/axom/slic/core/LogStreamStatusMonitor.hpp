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
 * \brief Monitor log streams to see if there are any pending messages
 *
 */
class LogStreamStatusMonitor
{
public:
  LogStreamStatusMonitor();
  virtual ~LogStreamStatusMonitor();

  /*!
   * \brief Add LogStream pointer to vector stored in LogStreamStatusMonitor.
   * \param [in] ls pointer to the user-supplied LogStream object.
   *
   */
  virtual void addStream(LogStream* ls);

  /*!
   * \brief Checks to see if any pending messages exist on the current MPI communicator
   */
  virtual bool hasPendingMessages() const;

protected:

  std::vector<LogStream*> m_streamVec;

};

} /* namespace slic */

} /* namespace axom */

#endif /* LOGSTREAMSTATUSMONITOR_HPP_ */
