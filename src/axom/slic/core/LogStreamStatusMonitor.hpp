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

#if defined(AXOM_USE_MPI)
  #include <mpi.h>
#endif

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

  /*!
   * \brief Add LogStream pointer to vector stored in LogStreamStatusMonitor.
   * \param [in] ls pointer to the user-supplied LogStream object.
   *
   * \note All ranks must call this function with the same set of MPI communicators 
   * in the same order.
   */
  void addStream(LogStream* ls);

  /*!
   * \brief Checks to see if any pending messages exist on the current MPI communicator.
   * 
   * \note This call is collective.  All ranks in m_mpiComm must call this function.
   */
  bool hasPendingMessages() const;

  /*!
   * \brief Finalize/clear member data
   */
  void finalize();

protected:
  std::vector<LogStream*> m_streamVec;
#if defined(AXOM_USE_MPI)
  bool m_useMPI;
  std::vector<MPI_Comm> m_mpiComms;
#endif
};

} /* namespace slic */

} /* namespace axom */

#endif /* LOGSTREAMSTATUSMONITOR_HPP_ */
