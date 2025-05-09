// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/FileUtilities.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <cerrno>

#include <cstdio>  // defines FILENAME_MAX

#ifdef WIN32
  #include <direct.h>
  #include <sys/stat.h>

  #define GetCurrentDir _getcwd
  #define ChangeCurrentDir _chdir
  #define Stat _stat
  #define Unlink _unlink
#else
  #include <unistd.h>    // for getcwd
  #include <sys/stat.h>  // for stat

  #define GetCurrentDir getcwd
  #define ChangeCurrentDir chdir
  #define Stat stat
  #define Unlink unlink
#endif

// Note: The hard-wired path separator character in this file
// should be set to the backslash when on Windows.

namespace axom
{
namespace utilities
{
namespace filesystem
{
std::string getCWD()
{
  char cCurrentPath[FILENAME_MAX];

  if(!GetCurrentDir(cCurrentPath, FILENAME_MAX))
  {
    //Note: Cannot use logging in COMMON component -- topic of JIRA issue
    // ATK-463
    //SLIC_WARNING("Common::Could not find cwd");

    return std::string("./");
  }

  return std::string(cCurrentPath);
}

int changeCWD(const std::string& dirName) { return ChangeCurrentDir(dirName.c_str()); }

//-----------------------------------------------------------------------------
bool pathExists(const std::string& fileName)
{
  // Uses system's stat() function.
  // Return code 0 indicates file exists
  struct Stat buffer;
  return (Stat(fileName.c_str(), &buffer) == 0);
}

//-----------------------------------------------------------------------------
std::string joinPath(const std::string& fileDir,
                     const std::string& fileName,
                     const std::string& separator)
{
  // Check if we need to add a separator
  bool pathNeedsSep = !fileDir.empty() && (fileDir[fileDir.size() - 1] != separator[0]);

  // Concatenate the path with the fileName to create the full path
  std::stringstream fullFileNameStream;
  fullFileNameStream << fileDir << (pathNeedsSep ? separator : "") << fileName;

  return fullFileNameStream.str();
}

//-----------------------------------------------------------------------------
int makeDirsForPath(const std::string& path)
{
  char separator = '/';
  std::string::size_type pos = 0;
  int err = 0;

  do
  {
    pos = path.find(separator, pos + 1);
    std::string dir_name = path.substr(0, pos);
#ifdef WIN32
    err = _mkdir(dir_name.c_str());
#else
    mode_t mode = 0777;  // rwx permissions for everyone
    err = mkdir(dir_name.c_str(), mode);
#endif
    err = (err && (errno != EEXIST)) ? 1 : 0;

  } while(pos != std::string::npos);

  return err;
}

//-----------------------------------------------------------------------------
std::string prefixRelativePath(const std::string& path, const std::string& prefix)
{
  if(path.empty())
  {
    throw std::invalid_argument("path must not be empty");
  };
  if(path[0] == '/' || prefix.empty())
  {
    return path;
  }
  return utilities::filesystem::joinPath(prefix, path);
}

//-----------------------------------------------------------------------------
std::string getParentPath(const std::string& path)
{
  if(path.empty())
  {
    throw std::invalid_argument("path must not be empty");
  };

  char separator = '/';

  std::string parent;

  if(path.size() == 1 && path[0] == separator)
  {
    // path is root, so parent is blank.
  }
  else
  {
    std::size_t found = path.rfind(separator);

    if(found != std::string::npos)
    {
      if(found == 0)
      {
        ++found;
      }
      parent = path.substr(0, found);
    }
  }

  return parent;
}

//-----------------------------------------------------------------------------
void getDirName(std::string& dir, const std::string& path)
{
  char separator = '/';

  std::size_t found = path.rfind(separator);
  if(found != std::string::npos)
  {
    dir = path.substr(0, found);
  }
  else
  {
    dir = "";
  }
}

//-----------------------------------------------------------------------------
int removeFile(const std::string& filename) { return Unlink(filename.c_str()); }

}  // end namespace filesystem
}  // end namespace utilities
}  // end namespace axom
