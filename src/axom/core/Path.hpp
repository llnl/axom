// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_PATH_H_
#define AXOM_CORE_PATH_H_

#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

namespace axom
{
/*!
 * \brief Path class for performing basic path operations with a user-selectable
 * delimiter character.
 */
class Path
{
public:
  /*!
     * \brief Constructs a path from a string
     *
     * \param [in] path The string representing the path
     * \param [in] delim The character that delimits the path
     *
     * Empty parts in \p path are removed
     * 
     * FIXME: Should this be marked \p explicit ?
     */
  Path(const std::string& path, const char delim = '/');

  /*!
     * \brief Default constructor for an empty path
     */
  Path() = default;

  /*!
     * \brief Constructs a path by concatenating a set of paths
     *
     * \param [in] paths The paths to concatenate/join
     * \param [in] delim The character to join the paths with
     * TODO: Should this be something easier-to-use, like the delimiter used
     * by the first Path object in the list?
     * 
     * Empty paths ins \p paths will be ignored.
     * 
     * \return The joined paths
     */
  static Path join(std::initializer_list<Path> paths, const char delim = '/');

  /*!
     * \brief Converts a path to a string
     * 
     * FIXME: Should this be marked \p explicit ?
     */
  operator std::string() const;

  /*!
     * \brief Returns a path representing the parent of the calling path
     * If the calling path is the root, the root will be returned
     */
  Path parent() const;

  /*!
     * \brief Returns the basename of the path (the last component)
     * \see basename(3)
     */
  std::string baseName() const;

  /*!
     * \brief Returns dirname of the path (all but the last component)
     * \see dirname(3)
     */
  std::string dirName() const;

  /*!
     * \brief Splits the path into the "parent" (dirname) and the last component (basename)
     * 
     * \return A pair {dirname, basename}
     */
  std::pair<std::string, std::string> split() const;

  /*!
   * \brief Returns a const reference to the parts of the path, to be iterated over
   */
  const std::vector<std::string>& parts() const { return m_components; }

private:
  // The components (tokens) that make up the path
  std::vector<std::string> m_components;

  // The delimiter character
  char m_delim = '/';

  // Indicates whether the path has a leading delimiter character
  bool m_leading_delim = false;
};

/*!
 * \brief Equality operator (equals)
 * \param [in] lhs The first Path to compare
 * \param [in] rhs The second Path to compare
 * Two paths are equal if they use the same delimiter and contain the same components
 */
bool operator==(const Path& lhs, const Path& rhs);

/*!
 * \brief Equality operator (not equals)
 * \param [in] lhs The first Path to compare
 * \param [in] rhs The second Path to compare
 * Two paths are equal if they use the same delimiter and contain the same components
 */
inline bool operator!=(const Path& lhs, const Path& rhs) { return !(lhs == rhs); }

}  // end namespace axom

#endif  // AXOM_CORE_PATH_H_
