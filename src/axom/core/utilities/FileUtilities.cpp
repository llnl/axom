// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/fmt.hpp"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cerrno>
#include <cstdio>  // defines FILENAME_MAX
#include <cstdlib>

#ifdef WIN32
  #include <windows.h>
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
    //Note: Cannot use slic logging in core component
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
  namespace sutil = axom::utilities::string;

  // Check if we need to add or remove a separator
  const bool has_empties = fileDir.empty() || fileName.empty();
  const int sep_count = (!fileDir.empty() && sutil::endsWith(fileDir, separator) ? 1 : 0) +
    (!fileName.empty() && sutil::startsWith(fileName, separator) ? 1 : 0);

  // Concatenate the path with the fileName, adding or removing a separator, as needed
  return axom::fmt::format("{}{}{}",
                           (sep_count == 2) ? sutil::removeSuffix(fileDir, separator) : fileDir,
                           !has_empties && sep_count == 0 ? separator : "",
                           fileName);
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

//-----------------------------------------------------------------------------
TempFile::TempFile(const std::string& file_name, const std::string& ext)
{
#ifdef WIN32
  char temp_dir[MAX_PATH];
  char temp_file_name[MAX_PATH];
  GetTempPathA(MAX_PATH, temp_dir);

  // Combine file_name and ext for the prefix (max 3 chars for prefix in GetTempFileNameA)
  GetTempFileNameA(temp_dir, file_name.c_str(), 0, temp_file_name);

  if(!ext.empty())
  {
    // remove ".tmp" (if present), add the requested extension and rename the file
    const std::string new_path =
      joinPath(axom::utilities::string::removeSuffix(temp_file_name, ".tmp"), ext, ".");
    if(std::rename(temp_file_name, new_path.c_str()) != 0)
    {
      throw std::ios::failure {"Failed to rename temp file to include extension"};
    }
    m_path = new_path;
  }
  else
  {
    m_path = std::string(temp_file_name);
  }
#else
  // create a tmp file with the requested prefix
  // note: mkstemp requires the last six chars to be "XXXXXX"
  const char* tmpdir = getenv("TMPDIR");
  const std::string dir = tmpdir ? tmpdir : "/tmp";
  const std::string tmp_file_name = joinPath(dir, file_name + "XXXXXX");
  std::vector<char> buf(tmp_file_name.begin(), tmp_file_name.end());
  buf.push_back('\0');

  const int fd = mkstemp(buf.data());
  if(fd == -1)
  {
    throw std::ios::failure {"Failed to create temp file"};
  }

  // rename to use the provided extension
  if(!ext.empty())
  {
    const std::string new_path = joinPath(buf.data(), ext, ".");
    if(std::rename(buf.data(), new_path.c_str()) != 0)
    {
      ::close(fd);

      throw std::ios::failure {"Failed to rename temp file to include extension"};
    }
    m_path = new_path;
  }
  else
  {
    m_path = buf.data();
  }
  ::close(fd);
#endif
}

TempFile::~TempFile()
{
  this->close();

  if(!m_retain_file)
  {
    removeFile(m_path);
  }
}

std::string TempFile::getFileContents() const
{
  std::stringstream buffer;

  std::ifstream ifs(m_path.c_str(), std::ios::in);
  if(ifs.is_open())
  {
    buffer << ifs.rdbuf();
  }

  return buffer.str();
}

}  // end namespace filesystem
}  // end namespace utilities
}  // end namespace axom
