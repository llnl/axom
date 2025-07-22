// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef COMMON_FILE_UTILITIES_H_
#define COMMON_FILE_UTILITIES_H_

#include <string>
#include <fstream>
#include <utility>

namespace axom
{
namespace utilities
{
namespace filesystem
{
/*!
 * \brief Returns the current working directory within the file system
 * \return path of cwd if found, else, returns "./"
 */
std::string getCWD();

/*!
 * \brief Changes the current working directory within the file system
 * \param [in] dirName an absolute or relative directory path
 * \return Status code 0 for success, non-zero otherwise
 */
int changeCWD(const std::string& dirName);

/*!
 * \brief Checks if the path pointed to by fileName exists
 * \param [in] fileName name of the file to check
 * \return true if file system contains a file named fileName, false otherwise
 */
bool pathExists(const std::string& fileName);

/*!
 * \brief Joins two strings with the given separator char
 *
 * \param [in] fileDir The directory of the file
 * \param [in] fileName The name of the file
 * \param [in] separator a single character to seperate the two strings
 *
 * \note fileDir can be the empty string
 * \note fileName can include directories (e.g. a/b/c.txt), but should not be an
 * absolute path f a non-empty fileDir is supplied.
 *
 * \returns The concatenated string: fileDir + fileName, with separator in
 *  between, if necessary
 * \note Example1:  joinPath("abc", "def") -> "abc/def"
 * \note Example2:  joinPath("abc/", "def") -> "abc/def"
 * \note Example3:  joinPath("abc/", "def/ghi") -> "abc/def/ghi"
 */
std::string joinPath(const std::string& fileDir,
                     const std::string& fileName,
                     const std::string& separator = "/");

/*!
 * \brief Make directories for a given path string
 *
 * \param [in] path  string representing an absolute or relative directory path
 *
 * Everything in the path is assumed to be intended to be a directory.  If
 * a directory in the path already exists, nothing is done.  If a directory
 * doesn't exist, it is created.
 */
int makeDirsForPath(const std::string& path);

/*!
 * \brief Add a prefix to a path if it is relative.
 *
 * \param [in] path non-empty string representing an absolute or relative path
 * \param [in] prefix string representing a directory path
 *
 * \return \c prefix + \c path (with delimiter) if the path is
 * relative or \c prefix is empty, or just \c path otherwise.
 */
std::string prefixRelativePath(const std::string& path, const std::string& prefix);

/*!
 * \brief Get parent path name from a filesystem path.
 *
 * \param [in] path a non-empty absolute or relative filesystem path
 *
 * \return a directory path formed by removing the last part of the
 *    input path
 *
 * If path is the system root, its parent is the empty string.
 */
std::string getParentPath(const std::string& path);

/*!
 * \brief Get directory name from a path that contains a file name
 *
 * \param [out] dir  a directory path formed by removing the file name from
 *                   the input path
 * \param [in] path  an absolute or relative directory/file path
 *
 * This function assumes that the input path has a file name at the end, and
 * it removes that file name, leaving a string containing only a directory
 * path.
 *
 * For example, if the path string is "abc/def/ghi/file.txt", the output dir
 * string will be "abc/def/ghi".  If the path string is "file.txt", the
 * output dir string will be "" (the empty string).
 */
void getDirName(std::string& dir, const std::string& path);

/*!
 * \brief Remove the specified file.
 * \param filename The name of the file.
 * \return 0 on success, -1 on failure. errno can obtain more information
 *         about the failure.
 *
 * \note On Windows, this function calls _unlink() and will fail if there are
 *       any open file handles to the specified file.
 *       On Linux, this function calls unlink() and will succeed even if
 *       there are open file handles to the specified file.
 */
int removeFile(const std::string& filename);

/**
 * /brief Utility class for managing temporary files.
 *
 * The TempFile class provides a convenient way to create and manage temporary files.
 * It ensures that the file is deleted upon destruction, unless the user calls \a retain(true)
 * The class is non-copyable to prevent accidental duplication of file handles.
 *
 * \note The path to the temp file is likely different from the supplied \a file_name
 */
class TempFile
{
public:
  /**
   * \brief Creates a temp file using \a file_name and \a file_ext
   * 
   * \param file_name The name of the temporary file to create without extension
   * \param file_ext An optional extension for the temp file
   * \note When creating the temp file, the name will likely be changed. You can get the actual file name
   * using the \a getPath() function after the file is created.
   */
  explicit TempFile(const std::string& file_name, const std::string& file_ext = "");

  ~TempFile();

  // Non-copyable
  TempFile(const TempFile&) = delete;
  TempFile(TempFile&&) = delete;
  TempFile& operator=(const TempFile&) = delete;
  TempFile& operator=(TempFile&&) = delete;

  /// If set to true, we will retain the file after the instance is destroyed
  void retain(bool should_retain) { m_retain_file = should_retain; }

  /// Returns true if the file will be retained after the instance is destroyed, false otherwise (default)
  bool retain() const { return m_retain_file; }

  /**
   * \brief Opens the temporary file for writing
   * 
   * \param mode ios flags for opening the file (not used if file is already open)
   * \return true if the file was successfully opened, false otherwise.
   */
  bool open(std::ios_base::openmode mode = std::ios::out | std::ios::trunc)
  {
    if(!m_ofs.is_open())
    {
      m_ofs.open(m_path, mode);
    }
    return m_ofs.is_open();
  }

  /// \brief Closes the temporary file if it is open
  void close()
  {
    if(m_ofs.is_open())
    {
      m_ofs.close();
    }
  }

  /**
   * Writes \a contents to the file
   * 
   * Opens the file if it is not already open using \a mode and optionally closes the file
   * \param content A string to write to the file
   * \param mode The ios mode to open the file
   * \param preserve_file_state Controls whether the file should be closed if we had to open it
   *   If false, the file will be left open. If true, it will be left in its initial state
   */
  template <typename StringType>
  void write(StringType&& content,
             std::ios::openmode mode = std::ios::out | std::ios::trunc,
             bool preserve_file_state = true)
  {
    bool opened_here = false;

    if(!is_open())
    {
      open(mode);
      opened_here = true;
    }

    m_ofs << std::forward<StringType>(content);

    if(opened_here && preserve_file_state)
    {
      close();
    }
  }

  /**
   * \brief Checks if the temporary file is currently open.
   * \return true if the file is open, false otherwise.
   */
  bool is_open() const { return m_ofs.is_open(); }

  /// Returns the path to the temporary file
  const std::string& getPath() const { return m_path; }

  /// Loads the contents of the file into a string and returns it
  std::string getFileContents() const;

  /// Overload the << operator to write data to the file (for general types)
  template <typename T>
  TempFile& operator<<(const T& data)
  {
    m_ofs << data;
    return *this;
  }

  /// Overload the << operator to support manipulators (e.g., std::endl)
  TempFile& operator<<(std::ostream& (*manip)(std::ostream&))
  {
    m_ofs << manip;
    return *this;
  }

private:
  std::string m_path;
  std::ofstream m_ofs;
  bool m_retain_file {false};  // should the temp file persist after the instance goes out of scope?
};

}  // end namespace filesystem
}  // end namespace utilities
}  // end namespace axom

#endif  //  COMMON_FILE_UTILITIES_H_
