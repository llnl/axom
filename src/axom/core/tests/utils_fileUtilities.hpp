// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include <fstream>
#include <exception>

#include "axom/config.hpp"
#include "axom/core/utilities/FileUtilities.hpp"

namespace fs = axom::utilities::filesystem;

TEST(utils_fileUtilities, getCWD_smoke)
{
  // This test just checks that we can call the getCWD function
  // It does not in any way confirm the results

  std::cout << "Checking that we can call getCWD()" << std::endl;

  std::string cwd = fs::getCWD();

  std::cout << " CWD is: " << cwd << std::endl;

  SUCCEED();
}

TEST(utils_fileUtilities, joinPath)
{
  std::cout << "Testing joinPath() function" << std::endl;

  std::string fdir = "abc";
  std::string fdirWithSlash = "abc/";
  std::string fname = "def";

  std::string fullfile = "abc/def";

  EXPECT_EQ(fullfile, fs::joinPath(fdir, fname));

  EXPECT_EQ(fullfile, fs::joinPath(fdirWithSlash, fname));

  std::string fnameWithSubdir = "def/ghi";
  EXPECT_EQ("abc/def/ghi", fs::joinPath(fdir, fnameWithSubdir));
}

TEST(utils_fileUtilities, pathExists)
{
  std::cout << "Testing pathExists() function" << std::endl;

  std::string cwd = fs::getCWD();

  // Test file that we know is present (i.e. the working directory)
  {
    EXPECT_TRUE(fs::pathExists(cwd));
  }

  // Test on file that we know is not present
  {
    const std::string missing = "m_i_s_s_i_n_g__f_i_l_e";
    EXPECT_FALSE(fs::pathExists(fs::joinPath(cwd, missing)));
  }

  // Test a file relative to AXOM_SRC_DIR
  {
    EXPECT_TRUE(fs::pathExists(AXOM_SRC_DIR));

    std::string dataDir = fs::joinPath(AXOM_SRC_DIR, "axom");
    EXPECT_TRUE(fs::pathExists(dataDir));

    std::string fileName = fs::joinPath(dataDir, "config.hpp.in");
    EXPECT_TRUE(fs::pathExists(fileName));
  }
}

TEST(utils_fileUtilities, changeCWD_smoke)
{
  std::cout << "Testing 'changeCWD()'" << std::endl;

  // Save copy of cwd at start of test
  std::string origCWD = fs::getCWD();
  std::cout << "[Original cwd]: '" << origCWD << "'" << std::endl;

  // Update cwd to new directory
  std::string newCWD = AXOM_SRC_DIR;
  EXPECT_TRUE(fs::pathExists(newCWD));
  std::cout << "Changing directory to: '" << newCWD << "'" << std::endl;

  int rc = fs::changeCWD(newCWD);
  EXPECT_EQ(0, rc);
  std::cout << "[Updated cwd]: '" << fs::getCWD() << "'" << std::endl;

  // Note: newCWD might contain symbolic links,
  // so don't directly compare newCWD and getCWD()
  if(origCWD != newCWD)
  {
    EXPECT_NE(origCWD, fs::getCWD());
  }

  // Change back to original directory
  rc = fs::changeCWD(origCWD);
  EXPECT_EQ(0, rc);

  EXPECT_EQ(origCWD, fs::getCWD());
  std::cout << "[cwd after change]: '" << fs::getCWD() << "'" << std::endl;
}

TEST(utils_fileUtilities, prefixRelativePath)
{
  using namespace axom::utilities::filesystem;
  EXPECT_EQ(prefixRelativePath("rel/path", "/pre/fix"), "/pre/fix/rel/path");
  EXPECT_EQ(prefixRelativePath("rel/path", ""), "rel/path");
  EXPECT_THROW(prefixRelativePath("", "/pre/fix"), std::invalid_argument);

  // These full paths should not change.
  EXPECT_EQ(prefixRelativePath("/full/path", "/pre/fix"), "/full/path");
  EXPECT_EQ(prefixRelativePath("/full/path", ""), "/full/path");
}

TEST(utils_fileUtilities, getParentPath)
{
  using namespace axom::utilities::filesystem;
  EXPECT_EQ(getParentPath("/full/multi/level/path"), "/full/multi/level");
  EXPECT_EQ(getParentPath("/full/multi/level"), "/full/multi");
  EXPECT_EQ(getParentPath("rel/multi/level/path"), "rel/multi/level");
  EXPECT_EQ(getParentPath("rel/multi/level"), "rel/multi");
  EXPECT_EQ(getParentPath("level"), "");
  EXPECT_EQ(getParentPath("/level0/level1"), "/level0");
  EXPECT_EQ(getParentPath("/level0"), "/");
  EXPECT_EQ(getParentPath("/"), "");
}
