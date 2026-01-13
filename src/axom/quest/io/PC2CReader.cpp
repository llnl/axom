// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/PC2CReader.hpp"

#ifndef AXOM_USE_C2C
  #error PC2CReader should only be included when Axom is configured with C2C
#endif

#include "axom/core.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace quest
{
namespace
{
constexpr int READER_SUCCESS = 0;
constexpr int READER_FAILED = -1;

}  // end anonymous namespace

//------------------------------------------------------------------------------
PC2CReader::PC2CReader(MPI_Comm comm) : m_comm(comm)
{
  MPI_Comm_rank(m_comm, &m_my_rank);
  MPI_Comm_size(m_comm, &m_num_ranks);
}

//------------------------------------------------------------------------------
int PC2CReader::read()
{
  SLIC_ASSERT(m_comm != MPI_COMM_NULL);

  // Clear internal data-structures
  this->clear();

  int rc = READER_FAILED;  // return code

  switch(m_my_rank)
  {
  // handle rank 0
  case 0:
    rc = C2CReader::read();
    if(m_num_ranks <= 1)
    {
      return rc;
    }

    bcast_int(rc);
    if(rc == READER_SUCCESS)
    {
      // broadcast number of curves, followed by curve data
      bcast_int(m_nurbsData.size());
      for(auto& curve : m_nurbsData)
      {
        // broadcast knot vector
        bcast_array(curve.getKnots().getArray());

        // broadcast control points
        bcast_array(curve.getControlPoints());

        // broadcast rational flag and weights
        const bool isRational = bcast_bool(curve.isRational());
        if(isRational)
        {
          bcast_array(curve.getWeights());
        }
      }
    }
    break;
  // handle other ranks
  default:
    rc = bcast_int();
    if(rc == READER_SUCCESS)
    {
      // Receive and reconstruct each NURBSCurve
      const int numNURBS = bcast_int();
      m_nurbsData.clear();
      m_nurbsData.reserve(numNURBS);
      for(int i = 0; i < numNURBS; ++i)
      {
        // receive knot vector
        axom::Array<double> knotsArr;
        bcast_array(knotsArr);

        // receive control points
        NURBSCurve::CoordsVec ctrlPts;
        bcast_array(ctrlPts);

        // receive rationality flag and weights
        const bool isRational = bcast_bool();
        if(isRational)
        {
          axom::Array<double> weights;
          bcast_array(weights);

          m_nurbsData.emplace_back(NURBSCurve {ctrlPts, weights, knotsArr});
        }
        else
        {
          m_nurbsData.emplace_back(NURBSCurve {ctrlPts, knotsArr});
        }
      }
    }
    break;
  }

  return rc;
}

/// MPI broadcasts an integer from rank 0
int PC2CReader::bcast_int(int value)
{
  MPI_Bcast(&value, 1, axom::mpi_traits<int>::type, 0, m_comm);
  return value;
}

/// MPI broadcasts a bool from rank 0
bool PC2CReader::bcast_bool(bool value)
{
  int intValue = value ? 1 : 0;
  MPI_Bcast(&intValue, 1, axom::mpi_traits<int>::type, 0, m_comm);
  return static_cast<bool>(intValue);
}

}  // namespace quest
}  // namespace axom
