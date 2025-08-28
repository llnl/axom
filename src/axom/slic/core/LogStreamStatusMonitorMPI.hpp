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
 * \brief Add a description later
 *
 */
class LogStreamStatusMonitorMPI : public LogStreamStatusMonitor
{
public:
  LogStreamStatusMonitorMPI();

  /*!
   * \brief ADD DESCRIPTION HERE.  This is where we add to m_streamVec and combine its communicator with the current m_mpiComm.
   *
   * \param [in] ls pointer to the user-supplied LogStream object.
   *
   * \note Something something something add later
   */
  void addStream(LogStream* ls) override;

  /*!
   * \brief ADD DESCRIPTION HERE.  This is the main point of this class!
   *
   * \note Something something something add later
   */
  bool hasPendingMessages() const override;

  /*!
   * \brief ADD DESCRIPTION HERE. Free MPI comm and that's it.
   *
   * \note Something something something add later
   */
  void finalize() override;

private:

  bool m_useMPI;
  MPI_Comm m_mpiComm;
  //std::vector<MPI_Comm> m_mpiComm;

};

} /* namespace slic */

} /* namespace axom */

#endif /* LOGSTREAMStatusMonitorMPI_HPP_ */
