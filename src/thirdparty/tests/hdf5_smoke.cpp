// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: hdf5_smoke.cpp
///
/// Simple smoke test for hdf5 library
///
/// The test assumes an Axom configuration with HDF5 enabled
/// and does not guard against its availability.
//-----------------------------------------------------------------------------

#include "axom/config.hpp"  // Brings in AXOM_USE_HDF5

#include "hdf5.h"

#include <iostream>
#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
TEST(hdf5_smoke, check_axom_define)
{
  // This test checks that axom has the correct defines for hdf5
  // This is brought in through axom's config.hpp

#ifdef AXOM_USE_HDF5
  SUCCEED();
#else
  FAIL() << "Bad configuration of Axom. "
         << "Using HDF5 smoke test, but AXOM_USE_HDF5 not defined";
#endif
}

TEST(hdf5_smoke, print_version)
{
  // This test prints the version of hdf5. It uses only the hdf5 headers
  std::cout << "HDF5 version is " << H5_VERSION << std::endl;
}

TEST(hdf5_smoke, create_dset)
{
  // This test uses the hdf5 dset API from a compiled library.
  // It is a slightly modified version of an hdf5 example obtained from:
  //   https://support.hdfgroup.org/HDF5/examples/intro.html#c
  // Valid return codes are non-negative

  hid_t file_id, dataset_id, dataspace_id;  // identifiers
  hsize_t dims[2];
  herr_t status;
  const char* FILE = "dset.h5";

  // Create a new file using default properties.
  file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  EXPECT_GE(file_id, 0);

  // Create the data space for the dataset.
  dims[0] = 4;
  dims[1] = 6;
  dataspace_id = H5Screate_simple(2, dims, nullptr);
  EXPECT_GE(dataspace_id, 0);

  // Create the dataset.
  dataset_id =
    H5Dcreate2(file_id, "/dset", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  EXPECT_GE(dataset_id, 0);

  // End access to the dataset and release resources used by it.
  status = H5Dclose(dataset_id);
  EXPECT_GE(status, 0);

  // Terminate access to the data space.
  status = H5Sclose(dataspace_id);
  EXPECT_GE(status, 0);

  // Close the file.
  status = H5Fclose(file_id);
  EXPECT_GE(status, 0);
}
