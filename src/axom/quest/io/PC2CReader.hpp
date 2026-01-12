// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_PC2CREADER_HPP_
#define QUEST_PC2CREADER_HPP_

#include "axom/config.hpp"

#ifndef AXOM_USE_C2C
  #error PC2CReader should only be included when Axom is configured with C2C
#endif

#include "axom/core/Macros.hpp"
#include "axom/quest/io/C2CReader.hpp"  // base class

#include "mpi.h"

namespace axom
{
namespace quest
{

class PC2CReader : public C2CReader
{
public:
  PC2CReader() = delete;
  PC2CReader(MPI_Comm comm);

  virtual ~PC2CReader() = default;

  /*!
   * \brief Reads in a C2C file to all ranks in the associated communicator
   *
   * \note Rank 0 reads in the C2C mesh file and broadcasts to the other ranks
   * \return status set to zero on success; non-zero otherwise
   */
  int read() final override;

private:
  /// MPI broadcasts an integer from rank 0 and returns the value to all ranks
  /// (other ranks do not need to supply a value)
  int bcast_int(int value = 0);

  /// MPI broadcasts a bool from rank 0 and returns the value to all ranks
  /// (other ranks do not need to supply a value)
  bool bcast_bool(bool value = false);

  /// MPI broadcasts an Array<double> or Array<PointType>
  template <typename ArrayType>
  void bcast_array(ArrayType& arr)
  {
    // Check that the value_type of the array is either double or Point<double, DIM>
    using value_type = typename ArrayType::value_type;
    static_assert(std::is_same_v<value_type, double> || primal::detail::is_point_v<value_type>);

    const bool is_root = (m_my_rank == 0);

    // first, send/receive the size(s)
    if(is_root)
    {
      bcast_int(arr.size());
    }
    else
    {
      arr.clear();
      const int size = bcast_int();
      arr.resize(ArrayOptions::Uninitialized {}, size);
    }

    // then, send/receive the data
    if constexpr(std::is_same_v<value_type, double>)
    {
      MPI_Bcast(arr.data(), arr.size(), axom::mpi_traits<double>::type, 0, m_comm);
    }
    else if constexpr(primal::detail::is_point_v<value_type>)
    {
      // since the data is contiguous, we can cast the start of the data to a double*
      const int numVals = arr.size() * value_type::DIMENSION;
      double* dataPtr = (numVals > 0) ? arr[0].data() : nullptr;
      MPI_Bcast(dataPtr, numVals, axom::mpi_traits<double>::type, 0, m_comm);
    }
  }

private:
  MPI_Comm m_comm {MPI_COMM_NULL};
  int m_my_rank {0};
  int m_num_ranks {-1};

  DISABLE_COPY_AND_ASSIGNMENT(PC2CReader);
  DISABLE_MOVE_AND_ASSIGNMENT(PC2CReader);
};

}  // namespace quest
}  // namespace axom

#endif  // QUEST_PC2CREADER_HPP_
