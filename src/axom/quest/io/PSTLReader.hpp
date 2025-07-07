// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_PSTLREADER_HPP_
#define QUEST_PSTLREADER_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/quest/io/STLReader.hpp"  // base class

#include "mpi.h"

namespace axom
{
namespace quest
{
class PSTLReader : public STLReader
{
public:
  PSTLReader() = delete;
  PSTLReader(MPI_Comm comm);

  virtual ~PSTLReader();

  /*!
   * \brief Reads in an STL file to all ranks in the associated communicator.
   * 
   * \note Rank 0 reads in the STL mesh file and broadcasts to the other ranks.
   * \return status set to zero on success; set to a non-zero value otherwise.
   */
  int read() final override;

private:
  MPI_Comm m_comm {MPI_COMM_NULL};
  int m_my_rank {0};

  DISABLE_COPY_AND_ASSIGNMENT(PSTLReader);
  DISABLE_MOVE_AND_ASSIGNMENT(PSTLReader);
};

}  // namespace quest
}  // namespace axom

#endif /* QUEST_PSTLREADER_HPP_ */
