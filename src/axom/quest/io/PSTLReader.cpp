// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/PSTLReader.hpp"

namespace axom::quest
{
namespace
{
constexpr int READER_SUCCESS = 0;
constexpr int READER_FAILED = -1;
}  // namespace

//------------------------------------------------------------------------------
PSTLReader::PSTLReader(MPI_Comm comm) : m_comm(comm) { MPI_Comm_rank(m_comm, &m_my_rank); }

//------------------------------------------------------------------------------
PSTLReader::~PSTLReader() = default;

//------------------------------------------------------------------------------
int PSTLReader::read()
{
  SLIC_ASSERT(m_comm != MPI_COMM_NULL);

  // Clear internal data-structures
  this->clear();

  int rc = -1;  // return code

  switch(m_my_rank)
  {
  case 0:

    rc = STLReader::read();
    if(rc == READER_SUCCESS)
    {
      MPI_Bcast(&m_num_nodes, 1, axom::mpi_traits<axom::IndexType>::type, 0, m_comm);
      MPI_Bcast(&m_nodes[0], m_num_nodes * 3, MPI_DOUBLE, 0, m_comm);
    }
    else
    {
      MPI_Bcast(&rc, 1, MPI_INT, 0, m_comm);
    }
    break;

  default:

    // Rank 0 broadcasts the number of nodes, a positive integer, if the
    // STL file is read successfully, or send a READER_FAILED flag, indicating
    // that the read was not successful.
    MPI_Bcast(&m_num_nodes, 1, axom::mpi_traits<axom::IndexType>::type, 0, m_comm);

    if(m_num_nodes != READER_FAILED)
    {
      rc = READER_SUCCESS;
      m_num_faces = m_num_nodes / 3;
      m_nodes.resize(m_num_nodes * 3);
      MPI_Bcast(&m_nodes[0], m_num_nodes * 3, MPI_DOUBLE, 0, m_comm);
    }
  }

  return (rc);
}

}  // namespace axom::quest
