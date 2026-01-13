// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include <fstream>
#include <exception>
#include <iterator>

#include "axom/config.hpp"
#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

#include "axom/fmt.hpp"

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

  // test with empty dir or name
  {
    EXPECT_EQ("def", fs::joinPath("", "def"));
    EXPECT_EQ("abc", fs::joinPath("abc", ""));
  }

  // test with default separator
  {
    EXPECT_EQ("abc/def", fs::joinPath("abc", "def"));
    EXPECT_EQ("abc/def", fs::joinPath("abc/", "def"));
    EXPECT_EQ("abc/def", fs::joinPath("abc", "/def"));
    EXPECT_EQ("abc/def", fs::joinPath("abc/", "/def"));
  }

  // test with a different separator
  {
    EXPECT_EQ("abc.def", fs::joinPath("abc", "def", "."));
    EXPECT_EQ("abc.def", fs::joinPath("abc.", "def", "."));
    EXPECT_EQ("abc.def", fs::joinPath("abc", ".def", "."));
    EXPECT_EQ("abc.def", fs::joinPath("abc.", ".def", "."));
  }

  // test with backslash separator (w/ and w/o raw string literals)
  {
    EXPECT_EQ("abc\\def", fs::joinPath("abc", "def", "\\"));
    EXPECT_EQ("abc\\def", fs::joinPath("abc\\", "def", "\\"));
    EXPECT_EQ("abc\\def", fs::joinPath("abc", "\\def", "\\"));
    EXPECT_EQ("abc\\def", fs::joinPath("abc\\", "\\def", "\\"));

    EXPECT_EQ(R"(abc\def)", fs::joinPath(R"(abc)", R"(def)", R"(\)"));
    EXPECT_EQ(R"(abc\def)", fs::joinPath(R"(abc\)", R"(def)", R"(\)"));
    EXPECT_EQ(R"(abc\def)", fs::joinPath(R"(abc)", R"(\def)", R"(\)"));
    EXPECT_EQ(R"(abc\def)", fs::joinPath(R"(abc\)", R"(\def)", R"(\)"));
  }

  // test a string that has the separator in other positions
  {
    EXPECT_EQ("abc/def/ghi", fs::joinPath("abc", "def/ghi"));
    EXPECT_EQ("abc/def/ghi", fs::joinPath("abc/", "def/ghi"));
    EXPECT_EQ("abc/def/ghi", fs::joinPath("abc", "/def/ghi"));
    EXPECT_EQ("abc/def/ghi", fs::joinPath("abc/", "/def/ghi"));

    EXPECT_EQ("abc/def/ghi", fs::joinPath("abc/def", "ghi"));
    EXPECT_EQ("abc/def/ghi", fs::joinPath("abc/def/", "ghi"));
    EXPECT_EQ("abc/def/ghi", fs::joinPath("abc/def", "/ghi"));
    EXPECT_EQ("abc/def/ghi", fs::joinPath("abc/def/", "/ghi"));
  }
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
  EXPECT_EQ(fs::prefixRelativePath("rel/path", "/pre/fix"), "/pre/fix/rel/path");
  EXPECT_EQ(fs::prefixRelativePath("rel/path", ""), "rel/path");
  EXPECT_THROW(fs::prefixRelativePath("", "/pre/fix"), std::invalid_argument);

  // These full paths should not change.
  EXPECT_EQ(fs::prefixRelativePath("/full/path", "/pre/fix"), "/full/path");
  EXPECT_EQ(fs::prefixRelativePath("/full/path", ""), "/full/path");
}

TEST(utils_fileUtilities, getParentPath)
{
  EXPECT_EQ(fs::getParentPath("/full/multi/level/path"), "/full/multi/level");
  EXPECT_EQ(fs::getParentPath("/full/multi/level"), "/full/multi");
  EXPECT_EQ(fs::getParentPath("rel/multi/level/path"), "rel/multi/level");
  EXPECT_EQ(fs::getParentPath("rel/multi/level"), "rel/multi");
  EXPECT_EQ(fs::getParentPath("level"), "");
  EXPECT_EQ(fs::getParentPath("/level0/level1"), "/level0");
  EXPECT_EQ(fs::getParentPath("/level0"), "/");
  EXPECT_EQ(fs::getParentPath("/"), "");
}

TEST(utils_fileUtilities, TempFile_basic)
{
  const std::string fname = "foo";
  std::string actual_path;
  {
    fs::TempFile fooFile(fname);
    actual_path = fooFile.getPath();
    std::cout << "Created temp file: " << actual_path << std::endl;
    EXPECT_TRUE(fs::pathExists(actual_path));
  }
  // file is removed once it is out of scope
  EXPECT_FALSE(fs::pathExists(actual_path));
}

TEST(utils_fileUtilities, TempFile_delete_during_destruction)
{
  const std::string fname = "foo";
  const std::string file_contents = "some string";

  for(bool should_retain : {true, false})
  {
    std::string actual_path;
    {
      fs::TempFile fooFile(fname);
      fooFile.retain(should_retain);
      actual_path = fooFile.getPath();

      fooFile.write(file_contents);
      EXPECT_EQ(file_contents, fooFile.getFileContents());

      EXPECT_TRUE(fs::pathExists(actual_path));
    }

    if(should_retain)
    {
      EXPECT_TRUE(fs::pathExists(actual_path));

      // check that the file contents match what we wrote above
      std::ifstream infile(actual_path);
      std::string contents((std::istreambuf_iterator<char>(infile)),
                           std::istreambuf_iterator<char>());
      EXPECT_EQ(contents, file_contents);

      fs::removeFile(actual_path);  // clean up
    }
    else
    {
      EXPECT_FALSE(fs::pathExists(actual_path));
    }
  }
}

TEST(utils_fileUtilities, TempFile_write)
{
  const std::string fname = "foo";
  const std::string file_contents = "some string";

  for(bool start_open : {true, false})
  {
    for(bool preserve_state : {true, false})
    {
      fs::TempFile fooFile(fname);

      if(start_open)
      {
        fooFile.open();
      }
      EXPECT_EQ(fooFile.is_open(), start_open);

      fooFile.write(file_contents, std::ios::out, preserve_state);

      if(preserve_state)
      {
        EXPECT_EQ(fooFile.is_open(), start_open)
          << axom::fmt::format("start_open {} -- preserve_state {}", start_open, preserve_state);
      }
      else
      {
        EXPECT_TRUE(fooFile.is_open())
          << axom::fmt::format("start_open {} -- preserve_state {}", start_open, preserve_state);
      }
    }
  }
}

TEST(utils_fileUtilities, TempFile_two)
{
  // check that we can create two TempFiles with the same prefix
  const std::string fname = "foo";

  fs::TempFile fooFile1(fname);
  fs::TempFile fooFile2(fname);

  const std::string actual_path1 = fooFile1.getPath();
  const std::string actual_path2 = fooFile2.getPath();
  EXPECT_NE(actual_path1, actual_path2);

  std::cout << "Created temp file 1: " << actual_path1 << std::endl;
  std::cout << "Created temp file 2: " << actual_path2 << std::endl;
}

TEST(utils_fileUtilities, TempFile_extension)
{
  for(const auto nm : {"", "foo"})
  {
    for(const auto ext : {"", ".json", "json"})
    {
      fs::TempFile temp(nm, ext);

      std::cout << "Creating temp file: '" << temp.getPath() << "' with extension '" << ext << "'\n";
      EXPECT_TRUE(axom::utilities::string::endsWith(temp.getPath(), ext));
    }
  }
}
