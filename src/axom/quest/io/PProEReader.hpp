// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_PPROEREADER_HPP_
#define QUEST_PPROEREADER_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/quest/io/ProEReader.hpp"  // base class

#include "mpi.h"

namespace axom
{
namespace quest
{
class PProEReader : public ProEReader
{
public:
  PProEReader() = delete;
  PProEReader(MPI_Comm comm);

  virtual ~PProEReader();

  /*!
   * \brief Reads in a Pro/E file to all ranks in the associated communicator.
   * 
   * \note Rank 0 reads in the Pro/E mesh file and broadcasts to the other ranks.
   * \return status set to zero on success; set to a non-zero value otherwise.
   */
  int read() final override;

private:
  MPI_Comm m_comm {MPI_COMM_NULL};
  int m_my_rank {0};

  DISABLE_COPY_AND_ASSIGNMENT(PProEReader);
  DISABLE_MOVE_AND_ASSIGNMENT(PProEReader);
};

}  // namespace quest
}  // namespace axom

#endif /* QUEST_PPROEREADER_HPP_ */
