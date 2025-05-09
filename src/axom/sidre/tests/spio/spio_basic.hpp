// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sidre/spio/IOManager.hpp"
#include "axom/sidre/spio/IOBaton.hpp"
#include "axom/fmt.hpp"

#include "mpi.h"

#include <map>
#include <string>

using axom::sidre::IOManager;

//------------------------------------------------------------------------------

// Check that the mapping between sidre and relay protocols is accurate
TEST(spio_basic, root_name)
{
  // Set up map from sidre protocols to expected relay protocols
  std::map<std::string, std::string> protocolMap;
  protocolMap["sidre_hdf5"] = "hdf5";
  protocolMap["conduit_hdf5"] = "hdf5";
  protocolMap["sidre_json"] = "json";
  protocolMap["conduit_bin"] = "json";
  protocolMap["json"] = "json";
  protocolMap["sidre_conduit_json"] = "conduit_json";

  for(const auto& kv : protocolMap)
  {
    const std::string& sidreProtocol = kv.first;
    const std::string& expRelayProtocol = kv.second;

    const std::string relayProtocol = IOManager::correspondingRelayProtocol(sidreProtocol);

    EXPECT_EQ(expRelayProtocol, relayProtocol);
  }
}

// Helper function to check spio's baton interface
void checkBaton(int numFiles, int numRanks, int myRank)
{
  // Sanity checks -- one or more files; at least one rank per file
  EXPECT_GT(numFiles, 0);
  EXPECT_GE(numRanks, numFiles);

  // Create a simple struct to hold expected values
  struct
  {
    int m_size;
    int m_firstIdx;
    int m_lastIdx;
    int m_id;
  } myGroup;

  // Find the floor of the group size
  myGroup.m_size = numRanks / numFiles;

  // The first few groups might be larger
  int numLarger = numRanks % numFiles;

  // Find the index that splits the larger and smaller ranks
  int splitIndex = (myGroup.m_size + 1) * numLarger;

  // Find the index of the group and its first element
  // And adjust group size, if necessary
  if(myRank < splitIndex)
  {
    myGroup.m_size += 1;

    myGroup.m_id = myRank / myGroup.m_size;
    myGroup.m_firstIdx = myGroup.m_size * myGroup.m_id;
  }
  else
  {
    myGroup.m_id = (myRank - splitIndex) / myGroup.m_size;
    myGroup.m_firstIdx = splitIndex + myGroup.m_size * myGroup.m_id;
    myGroup.m_id += numLarger;
  }
  myGroup.m_lastIdx = myGroup.m_firstIdx + myGroup.m_size - 1;

  //// Test the baton interface
  axom::sidre::IOBaton baton(MPI_COMM_WORLD, numFiles, numRanks);

  // Num files and groups
  EXPECT_EQ(numFiles, baton.getNumFiles());
  EXPECT_EQ(myGroup.m_size, baton.setSize());

  // First and last index in the group
  const bool expFirst = (myRank == myGroup.m_firstIdx);
  EXPECT_EQ(expFirst, baton.isFirstInGroup());

  const bool expLast = (myRank == myGroup.m_lastIdx);
  EXPECT_EQ(expLast, baton.isLastInGroup());

  // wait() returns group index
  int batonGroup = baton.wait();
  EXPECT_EQ(myGroup.m_id, batonGroup);

  // pass() returns 0 for success, 1 for fail
  int success = baton.pass();
  EXPECT_EQ(0, success);
}

// Test basic baton operations
TEST(spio_basic, baton)
{
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  // Test baton for different numbers of files
  for(int nFiles = 1; nFiles <= num_ranks; ++nFiles)
  {
    SCOPED_TRACE(axom::fmt::format("Checking baton for {} files with {} ranks", nFiles, num_ranks));
    checkBaton(nFiles, num_ranks, my_rank);
  }
}
