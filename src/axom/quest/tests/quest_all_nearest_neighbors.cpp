// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

// _quest_allnear_include_start
#include "axom/quest/AllNearestNeighbors.hpp"
// _quest_allnear_include_end
#include "axom/quest/detail/AllNearestNeighbors_detail.hpp"

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <fstream>
#include <sstream>

char* fname;
char* outfname;

template <typename T>
void verify_array(T* standard, T* expt, int n)
{
  int mismatches = 0;

  for(int i = 0; i < n; ++i)
  {
    if(!axom::utilities::isNearlyEqual(static_cast<double>(standard[i]), static_cast<double>(expt[i])))
    {
      ++mismatches;
      SLIC_INFO("i " << i << " std " << standard[i] << " expt " << expt[i]);
    }
  }

  if(mismatches > 0)
  {
    ADD_FAILURE() << " with " << mismatches << " mismatches.";
  }
  else
  {
    SUCCEED();
  }
}

/*!
 * \brief Find the closest point (in another region) to each given point
 * \param [in] x X-coordinates of input points
 * \param [in] y Y-coordinates of input points
 * \param [in] z Z-coordinates of input points
 * \param [in] region Region of each point
 * \param [in] n Number of points
 * \param [in] limit Max distance for all-nearest-neighbors query
 * \param [out] neighbor Index of nearest neighbor not in the same class
 * \param [out] sqdistance Squared distance to nearest neighbor
 * \pre x, y, z, and region have n entries
 * \pre neighbor is allocated with room for n entries
 *
 * This method compares each point to each other point, taking O(n^2) time.
 */
void all_nearest_neighbors_bruteforce(const double* x,
                                      const double* y,
                                      const double* z,
                                      const int* region,
                                      int n,
                                      double limit,
                                      int* neighbor,
                                      double* sqdistance)
{
  // O(n^2) brute force approach.  For each point i, test distance to all other
  // points and report result.

  double sqlimit = limit * limit;

  for(int i = 0; i < n; ++i)
  {
    sqdistance[i] = DBL_MAX;
    neighbor[i] = -1;
  }

  for(int i = 0; i < n; ++i)
  {
    // The j-loop has to go from 0, not i+1, because "closest" is not
    // symmetrical.
    for(int j = 0; j < n; ++j)
    {
      if(region[i] != region[j])
      {
        double sqdist = axom::quest::detail::squared_distance(x[i], y[i], z[i], x[j], y[j], z[j]);
        if(sqdist < sqdistance[i] && sqdist < sqlimit)
        {
          sqdistance[i] = sqdist;
          neighbor[i] = j;
        }
      }
    }
  }
}

//----------------------------------------------------------------------
TEST(quest_all_nearnbr, simple_2D_query)
{
  SLIC_INFO("*** Simple 2D all-nearest-neighbors query.");

  // _quest_allnear_input_start
  double x[] = {-1.2, -1.0, -0.8, -1.0, 0.8, 1.0, 1.2, 1.0};
  double y[] = {0.0, -0.2, 0.0, 0.2, 0.0, -0.2, 0.0, 0.2};
  double z[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int region[] = {0, 0, 0, 0, 1, 1, 1, 1};
  int n = 8;
  double limit = 1.9;
  // _quest_allnear_input_end
  int neighbor[] = {-1, -1, -1, -1, -1, -1, -1, -1};
  int expneighbor[] = {-1, 4, 4, 4, 2, 2, -1, 2};
  double dsq[8];
  double expdsq[] = {DBL_MAX, 3.28, 2.56, 3.28, 2.56, 3.28, DBL_MAX, 3.28};

  {
    SCOPED_TRACE("brute force limit 1.9");
    all_nearest_neighbors_bruteforce(x, y, z, region, n, limit, neighbor, dsq);
    verify_array(expneighbor, neighbor, n);
    verify_array(expdsq, dsq, n);
  }
  {
    SCOPED_TRACE("indexed limit 1.9");
    // _quest_allnear_query_start
    axom::quest::all_nearest_neighbors(x, y, z, region, n, limit, neighbor, dsq);
    // _quest_allnear_query_end
    verify_array(expneighbor, neighbor, n);
    verify_array(expdsq, dsq, n);
  }
}

//----------------------------------------------------------------------
TEST(quest_all_nearnbr, simple_3D_query)
{
  SLIC_INFO("*** Simple 3D all-nearest-neighbors query.");

  double x[] = {-1.2, -1.0, -0.8, -1.0, 0.8, 1.0, 1.2, 1.0};
  double y[] = {0.0, -0.2, 0.0, -0.1, 0.0, 0.2, 0.0, 0.1};
  double z[] = {0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.2};
  int region[] = {0, 0, 0, 0, 1, 1, 1, 1};
  int n = 8;
  double limit = 1.9;
  int neighbor[] = {-1, -1, -1, -1, -1, -1, -1, -1};
  int expneighbor[] = {-1, 4, 4, 4, 2, 2, -1, 2};
  double dsq[8];
  double expdsq[] = {DBL_MAX, 3.28, 2.56, 3.29, 2.56, 3.28, DBL_MAX, 3.29};

  {
    SCOPED_TRACE("brute force limit 1.9");
    all_nearest_neighbors_bruteforce(x, y, z, region, n, limit, neighbor, dsq);
    verify_array(expneighbor, neighbor, n);
    verify_array(expdsq, dsq, n);
  }
  {
    SCOPED_TRACE("indexed limit 1.9");
    axom::quest::all_nearest_neighbors(x, y, z, region, n, limit, neighbor, dsq);
    verify_array(expneighbor, neighbor, n);
    verify_array(expdsq, dsq, n);
  }
}

//----------------------------------------------------------------------
TEST(quest_all_nearnbr, cplx_13region_query)
{
  SLIC_INFO("*** 13-region closely-packed query.");

  // clang-format off
  double x[] = { -2.7, -2.3, -1.5, -1.2, -0.8, -0.9, -1.8,
                 -0.8, -0.3,  0.4,  1.4,  1.5,  0.9,
                 -2.6, -2.5, -2.0, -1.7, -1.4, -1.7, -2.0, -1.3,
                 -1.6, -1.3, -0.9, -0.8, -0.9, -1.1, -1.3, -1.4,
                 -0.9, -0.3,  0.2,  0.9,  0.9,  0.5, -0.5, -0.7,
                 1.0,  1.1,  1.3,  1.6,  2.0,  2.0,  2.3,  1.6,
                 -2.5, -1.9, -1.3, -1.2, -0.9, -1.1, -1.5, -2.1, -2.3,
                 -1.0, -0.9, -0.4,  0.0, -0.1, -0.6, -1.0,
                 0.1,  0.5,  1.1,  1.3,  0.9,  0.4,
                 1.3,  1.4,  2.1,  2.4,  2.3,  1.9,
                 -1.0, -0.3,  0.0, -0.3, -0.8,
                 -0.1,  0.0,  0.3,  0.5,  0.4,  0.3,  0.1,
                 0.7,  1.1,  1.8,  2.0,  1.8,  1.4,  1.0,  0.7,
                 3.5,  3.7,  4.0,  3.6 };

  double y[] = {  1.3,  1.2,  1.2,  1.2,  1.7,  2.3,  1.8,
                  2.0,  1.3,  0.9,  1.3,  2.1,  2.9,
                  0.8,  0.2,  0.4,  0.5,  0.8,  1.1,  1.0,  1.1,
                  0.5,  0.2,  0.3,  0.9,  1.4,  1.2,  1.0,  0.6,
                  0.0, -0.2, -0.4, -0.1,  0.7,  0.9,  1.3,  0.6,
                  0.7,  0.3, -0.1,  0.4,  0.0,  0.6,  0.8,  1.4,
                  -0.6, -0.9, -1.7, -0.8, -0.2,  0.1,  0.2,  0.3,  0.0,
                  -1.3, -1.9, -1.8, -1.0, -0.7, -0.2, -0.7,
                  -0.8, -1.3, -1.3, -0.5, -0.4, -0.5,
                  -0.3, -1.1, -1.5, -1.0,  0.4, -0.3,
                  -2.5, -2.8, -2.4, -2.0, -2.0,
                  -1.5, -2.2, -2.6, -1.9, -1.6, -1.3, -1.1,
                  -2.2, -2.6, -2.3, -1.7, -1.4, -1.3, -1.5,  0.7,
                  1.0,  0.8,  0.9,  1.5 };

  double z[97];

  int region[] = {  1,  1,  1,  1,  1,  1,  1,
                    2,  2,  2,  2,  2,  2,
                    3,  3,  3,  3,  3,  3,  3,  3,
                    4,  4,  4,  4,  4,  4,  4,  4,
                    5,  5,  5,  5,  5,  5,  5,  5,
                    6,  6,  6,  6,  6,  6,  6,  6,
                    7,  7,  7,  7,  7,  7,  7,  7,  7,
                    8,  8,  8,  8,  8,  8,  8,
                    9,  9,  9,  9,  9,  9,
                    10, 10, 10, 10, 10, 10,
                    11, 11, 11, 11, 11,
                    12, 12, 12, 12, 12, 12, 12,
                    13, 13, 13, 13, 13, 13, 13, 13,
                    14, 14, 14, 14 };
  // clang-format on

  int n = 97;
  double limit = 1.4;
  int bfneighbor[97];
  int idxneighbor[97];
  double bfsqdst[97];
  double idxsqdst[97];

  for(int i = 0; i < n; ++i)
  {
    z[i] = 0.;
    bfneighbor[i] = -1;
    idxneighbor[i] = -1;
  }

  {
    SCOPED_TRACE("Comparing brute force with indexed, limit 1.4");
    all_nearest_neighbors_bruteforce(x, y, z, region, n, limit, bfneighbor, bfsqdst);
    axom::quest::all_nearest_neighbors(x, y, z, region, n, limit, idxneighbor, idxsqdst);
    verify_array(bfneighbor, idxneighbor, n);
    verify_array(bfsqdst, idxsqdst, n);
  }
}

void readPointsFile(char* filename,
                    std::vector<double>& x,
                    std::vector<double>& y,
                    std::vector<double>& z,
                    std::vector<int>& region)
{
  std::ifstream infile(filename);
  std::string theline;
  double xi, yi, zi;
  int regioni;

  if(infile.good())
  {
    std::getline(infile, theline);
    do
    {
      if(theline.size() > 0)
      {
        std::stringstream line(theline);
        line >> xi >> yi >> zi >> regioni;
        x.push_back(xi);
        y.push_back(yi);
        z.push_back(zi);
        region.push_back(regioni);
      }
      std::getline(infile, theline);
    } while(!infile.eof() && infile.good());
  }
}

void writeNeigborsFile(char* filename, int* neighbors, int n)
{
  std::ofstream outfile(filename);
  if(outfile.good())
  {
    for(int i = 0; i < n; ++i)
    {
      outfile << neighbors[i] << std::endl;
    }
  }
}

TEST(quest_all_nearnbr, file_query)
{
  if(fname != nullptr)
  {
    SLIC_INFO("About to read file " << fname);

    std::vector<double> x, y, z;
    std::vector<int> region;
    int* bfneighbor;
    int* idxneighbor;
    double* bfsqdst;
    double* idxsqdst;

    readPointsFile(fname, x, y, z, region);

    const int n = static_cast<int>(region.size());

    SLIC_INFO("n is " << n);

    if((n > 0) && (n == static_cast<int>(x.size())) && (n == static_cast<int>(y.size())) &&
       (n == static_cast<int>(z.size())))
    {
      double limit = 2.1;
      bfneighbor = new int[n];
      idxneighbor = new int[n];
      bfsqdst = new double[n];
      idxsqdst = new double[n];

      for(int i = 0; i < n; ++i)
      {
        bfneighbor[i] = -1;
        idxneighbor[i] = -1;
      }

      {
        SCOPED_TRACE("Read file, compare brute force with indexed, limit 2.1");
        all_nearest_neighbors_bruteforce(&x[0], &y[0], &z[0], &region[0], n, limit, bfneighbor, bfsqdst);
        axom::quest::all_nearest_neighbors(&x[0], &y[0], &z[0], &region[0], n, limit, idxneighbor, idxsqdst);
        verify_array(bfneighbor, idxneighbor, n);
        verify_array(bfsqdst, idxsqdst, n);
      }

      if(outfname != nullptr)
      {
        writeNeigborsFile(outfname, idxneighbor, n);
      }

      delete[] bfneighbor;
      delete[] idxneighbor;
      delete[] bfsqdst;
      delete[] idxsqdst;
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;
  fname = nullptr;
  outfname = nullptr;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  if(argc > 1)
  {
    fname = argv[1];
    if(argc > 2)
    {
      outfname = argv[2];
    }
  }

  result = RUN_ALL_TESTS();

  return result;
}
