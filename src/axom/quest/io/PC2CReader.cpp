// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/readers/PC2CReader.hpp"

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
PC2CReader::PC2CReader(MPI_Comm comm) : m_comm(comm) { MPI_Comm_rank(m_comm, &m_my_rank); }

//------------------------------------------------------------------------------
int PC2CReader::read()
{
  using CoordsVec = typename NURBSCurve::CoordsVec;
  constexpr int DIM = 2;

  SLIC_ASSERT(m_comm != MPI_COMM_NULL);

  // Clear internal data-structures
  this->clear();

  int rc = -1;  // return code

  switch(m_my_rank)
  {
  case 0:
    rc = C2CReader::read();
    if(rc == READER_SUCCESS)
    {
      bcast_int(m_nurbsData.size());

      for(const auto& curve : m_nurbsData)
      {
        // broadcast sizes
        const int numCP = bcast_int(curve.getNumControlPoints());
        const int numKnots = bcast_int(curve.getNumKnots());
        const bool isRational = bcast_bool(curve.isRational());

        AXOM_UNUSED_VAR(numKnots);

        // broadcast control points as an array of doubles
        auto pt_v = numCP > 0
          ? axom::ArrayView<double>(const_cast<double*>(curve[0].data()), numCP * DIM)
          : axom::ArrayView<double>(nullptr, 0);
        bcast_data(pt_v);

        // broadcast knot vector as an array of doubles
        bcast_data(curve.getKnotsArray().view());

        // optionally broadcast control point weights
        if(isRational)
        {
          auto wts_v = numCP > 0
            ? axom::ArrayView<double>(const_cast<double*>(&curve.getWeight(0)), numCP)
            : axom::ArrayView<double>(nullptr, 0);
          bcast_data(wts_v);
        }
      }
    }
    else
    {
      bcast_int(rc);
    }
    break;

  default:

    // Rank 0 broadcasts the number of NURBSCurve entities, a positive integer, if the
    // C2C file is read successfully, or sends a READER_FAILED flag, indicating
    // that the read was not successful.
    rc = bcast_int();
    if(rc != READER_FAILED)
    {
      const int numNURBS = rc;
      rc = READER_SUCCESS;

      // Receive and reconstruct each NURBSCurve
      m_nurbsData.reserve(numNURBS);
      for(int i = 0; i < numNURBS; ++i)
      {
        // receive sizes
        const int numCP = bcast_int();
        const int numKnots = bcast_int();
        const bool isRational = bcast_bool();

        // receive control points
        // since they're contiguous, we can recieve them as array of doubles
        auto pts = CoordsVec(axom::ArrayOptions::Uninitialized {}, numCP, numCP);
        auto pts_v = numCP > 0 ? axom::ArrayView<double>(pts[0].data(), numCP * DIM)
                               : axom::ArrayView<double>(nullptr, 0);
        bcast_data(pts_v);

        // receive knot vector
        auto knots = axom::Array<double>(axom::ArrayOptions::Uninitialized {}, numKnots, numKnots);
        bcast_data(knots.view());

        // construct NURBSCurve (w/ optional weights)
        if(isRational)
        {
          // receive weights
          auto weights = axom::Array<double>(axom::ArrayOptions::Uninitialized {}, numCP, numCP);
          bcast_data(weights.view());

          m_nurbsData.emplace_back(NURBSCurve {pts, weights, knots});
        }
        else
        {
          m_nurbsData.emplace_back(NURBSCurve {pts, knots});
        }
      }
    }
  }

  return (rc);
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
  int intValue = static_cast<int>(value);
  MPI_Bcast(&intValue, 1, axom::mpi_traits<int>::type, 0, m_comm);
  return static_cast<bool>(intValue);
}

/// MPI broadcasts an array of doubles in-place in an axom::ArrayView
// assume all ranks already have the correct size
void PC2CReader::bcast_data(axom::ArrayView<double> arr)
{
  const int sz = arr.size();
  MPI_Bcast(arr.data(), sz, axom::mpi_traits<double>::type, 0, m_comm);
}

}  // namespace quest
}  // namespace axom
