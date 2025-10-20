// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Associated header file
#include "IOBaton.hpp"

const int axom::sidre::IOBaton::s_invalid_rank_id = -1;

namespace axom
{
namespace sidre
{
/*
 *************************************************************************
 *
 * IOBaton constructor
 *
 *************************************************************************
 */
IOBaton::IOBaton(MPI_Comm comm, int num_files, int num_groups)
  : m_mpi_comm(comm)
  , m_comm_size(1)
  , m_my_rank(0)
  , m_rank_before_me(s_invalid_rank_id)
  , m_rank_after_me(s_invalid_rank_id)
  , m_mpi_tag(MPI_ANY_TAG)
{
  SLIC_ASSERT(num_files > 0);
  SLIC_ASSERT(num_groups > 0);

  MPI_Comm_size(comm, &m_comm_size);
  MPI_Comm_rank(comm, &m_my_rank);
  m_num_files = num_files;
  m_num_groups = num_groups;

  int active_comm_size = m_comm_size;
  if(m_comm_size > m_num_groups)
  {
    active_comm_size = m_num_groups;
  }
  m_num_larger_sets = active_comm_size % num_files;
  if(m_my_rank < active_comm_size)
  {
    m_set_size = active_comm_size / m_num_files;
  }
  else
  {
    m_set_size = 1;
  }
  m_first_regular_set_rank = (m_set_size + 1) * m_num_larger_sets;
  if(m_my_rank < m_first_regular_set_rank)
  {
    m_set_id = m_my_rank / (m_set_size + 1);
    m_rank_within_set = m_my_rank % (m_set_size + 1);
    if(m_rank_within_set < m_set_size)
    {
      m_rank_after_me = m_my_rank + 1;
    }
  }
  else if(m_my_rank < active_comm_size)
  {
    m_set_id = m_num_larger_sets + (m_my_rank - m_first_regular_set_rank) / m_set_size;
    m_rank_within_set = (m_my_rank - m_first_regular_set_rank) % m_set_size;
    if(m_rank_within_set < m_set_size - 1)
    {
      m_rank_after_me = m_my_rank + 1;
    }
  }
  else
  {
    m_set_id = m_my_rank;
    m_rank_within_set = 0;
  }
  if(m_rank_within_set > 0)
  {
    m_rank_before_me = m_my_rank - 1;
  }
}

/*
 *************************************************************************
 *
 * IOBaton dtor destroys all contents.
 *
 *************************************************************************
 */
IOBaton::~IOBaton() { }

int IOBaton::wait()
{
  int return_val = -1;
  if(m_rank_before_me != s_invalid_rank_id)
  {
    MPI_Status mpi_stat;
    int baton;
    int mpi_err = MPI_Recv(&baton, 1, MPI_INT, m_rank_before_me, m_mpi_tag, m_mpi_comm, &mpi_stat);
    if(mpi_err == MPI_SUCCESS)
    {
      return_val = m_set_id;
    }
  }
  else
  {
    return_val = m_set_id;
  }
  return return_val;
}

int IOBaton::pass()
{
  int return_val = 0;
  if(m_rank_after_me != s_invalid_rank_id)
  {
    int baton = 0;
    int mpi_err = MPI_Ssend(&baton, 1, MPI_INT, m_rank_after_me, 0, m_mpi_comm);
    if(mpi_err != MPI_SUCCESS)
    {
      return_val = -1;
    }
  }
  return return_val;
}

} /* end namespace sidre */
} /* end namespace axom */
