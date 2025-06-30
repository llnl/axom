// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/PProEReader.hpp"

namespace axom
{
namespace quest
{
namespace
{
constexpr int READER_SUCCESS = 0;
constexpr int READER_FAILED = -1;
}  // namespace

//------------------------------------------------------------------------------
PProEReader::PProEReader(MPI_Comm comm) : m_comm(comm) { MPI_Comm_rank(m_comm, &m_my_rank); }

//------------------------------------------------------------------------------
PProEReader::~PProEReader()
{
  // Auto-generated destructor stub
}

//------------------------------------------------------------------------------
int PProEReader::read()
{
  SLIC_ASSERT(m_comm != MPI_COMM_NULL);

  // Clear internal data-structures
  this->clear();

  axom::IndexType rc = -1;  // return code

  switch(m_my_rank)
  {
  case 0:

    rc = ProEReader::read();
    if(rc == READER_SUCCESS)
    {
      MPI_Bcast(&m_num_nodes, 1, axom::mpi_traits<axom::IndexType>::type, 0, m_comm);
      MPI_Bcast(&m_num_tets, 1, axom::mpi_traits<axom::IndexType>::type, 0, m_comm);
      MPI_Bcast(&m_nodes[0], m_num_nodes * 3, MPI_DOUBLE, 0, m_comm);
      MPI_Bcast(&m_tets[0], m_num_tets * 4, MPI_INT, 0, m_comm);
    }  // END if
    else
    {
      MPI_Bcast(&rc, 1, axom::mpi_traits<axom::IndexType>::type, 0, m_comm);
    }  // END else
    break;

  default:

    // Rank 0 broadcasts the number of nodes, a positive integer, if the
    // Pro/E file is read successfully, or send a READER_FAILED flag, indicating
    // that the read was not successful.
    MPI_Bcast(&m_num_nodes, 1, axom::mpi_traits<axom::IndexType>::type, 0, m_comm);

    if(m_num_nodes != READER_FAILED)
    {
      rc = READER_SUCCESS;
      m_nodes.resize(m_num_nodes * 3);
      MPI_Bcast(&m_num_tets, 1, axom::mpi_traits<axom::IndexType>::type, 0, m_comm);
      m_tets.resize(m_num_tets * 4);
      MPI_Bcast(&m_nodes[0], m_num_nodes * 3, MPI_DOUBLE, 0, m_comm);
      MPI_Bcast(&m_tets[0], m_num_tets * 4, MPI_INT, 0, m_comm);
    }

  }  // END switch

  return (rc);
}

}  // end namespace quest
}  // end namespace axom
