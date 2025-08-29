// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file LogStreamStatusMonitorMPI.hpp
 *
 */

#ifndef LOGSTREAMSTATUSMPI_MONITOR_HPP_
#define LOGSTREAMSTATUSMPI_MONITOR_HPP_


#include <mpi.h>
#include "axom/slic/core/LogStreamStatusMonitor.hpp"

namespace axom
{
namespace slic
{
/*!
 * \class LogStreamStatusMonitorMPI
 *
 * \brief Monitor log streams to see if there are any pending messages using MPI
 *
 */
class LogStreamStatusMonitorMPI : public LogStreamStatusMonitor
{
public:
  LogStreamStatusMonitorMPI();

  /*!
   * \brief Add LogStream pointer to vector stored in LogStreamStatusMonitor
      and set m_useMPI if MPI is used.  Also add MPI communicator.
   * \param [in] ls pointer to the user-supplied LogStream object.
   *
   * \note It is assumed that the same MPI communicator is used for all LogStream objects
   *  added to this class.
   */
  void addStream(LogStream* ls) override;

  /*!
   * \brief Checks to see if any pending messages exist on the current MPI communicator
   */
  bool hasPendingMessages() const override;

private:

  bool m_useMPI;
  MPI_Comm m_mpiComm;

};

} /* namespace slic */

} /* namespace axom */

#endif /* LOGSTREAMStatusMonitorMPI_HPP_ */
