// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file IOBaton.hpp
 *
 * \brief   Header file containing definition of IOBaton class.
 *
 ******************************************************************************
 */

#ifndef SIDRE_IOBATON_HPP_
#define SIDRE_IOBATON_HPP_

#include "mpi.h"

// Other axom headers
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/slic/interface/slic.hpp"
#include "axom/export/sidre.h"

namespace axom
{
namespace sidre
{
/*!
 * \class IOBaton
 *
 * \brief IOBaton ensures that during I/O operations, only one rank will
 * interact with a particular file at one time.
 *
 * Each rank is placed into a set of ranks, with the number of sets being
 * equal to the number of I/O files, and then the ranks use the wait and
 * pass methods to pass control of I/O operations from one rank to the next
 * within the set.
 */
class IOBaton
{
public:
  /*!
   * \brief Constructor
   *
   * \param com        MPI communicator
   * \param num_files  Number of files involved in an I/O operation
   * \param num_groups Number of groups stored in the files (usually
   *                   equivalent to the number of MPI ranks providing data).
   */
  IOBaton(MPI_Comm com, int num_files, int num_groups);

  /*!
   * \brief Destructor
   */
  ~IOBaton();

  /*!
   * \brief Wait for previous rank to pass control to the local rank.
   *
   * \return An integer id for the set of which this rank is a member.
   */
  int wait();

  /*!
   * \brief Pass control to the next rank.
   *
   * \return  0 if sucessful, -1 if not.
   */
  int pass();

  /*!
   * \brief Size of local rank's set.
   *
   * \return Number of ranks in the set.
   */
  int setSize() const { return m_my_rank < m_first_regular_set_rank ? m_set_size + 1 : m_set_size; }

  /*!
   * \brief Tells if the local rank is the first (lowest) in its set.
   */
  bool isFirstInGroup() const { return (m_rank_within_set == 0); }

  /*!
   * \brief Tells if the local rank is the last (highest) in its set.
   */
  bool isLastInGroup() const { return (m_rank_after_me == s_invalid_rank_id); }

  /*!
   * \brief Get the number of files involved in the I/O operation.
   */
  int getNumFiles() const { return m_num_files; }

private:
  DISABLE_COPY_AND_ASSIGNMENT(IOBaton);

  void setupReducedRanks();

  AXOM_SIDRE_EXPORT static const int s_invalid_rank_id;

  MPI_Comm m_mpi_comm;

  int m_comm_size;        // num procs in the MPI communicator
  int m_my_rank;          // rank of this proc
  int m_num_files;        // number of files
  int m_num_groups;       // number of groups (input ranks)
  int m_num_larger_sets;  // some sets have one extra
  int m_set_size;         // regular set size (m_comm_size / m_num_files) w/o
                          // remainder
  int m_set_id;
  int m_first_regular_set_rank;
  int m_rank_within_set;
  int m_rank_before_me;
  int m_rank_after_me;

  int m_mpi_tag;
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_IOBATON_HPP_ */
