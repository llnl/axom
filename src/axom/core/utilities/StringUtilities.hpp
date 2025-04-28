// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef CORE_STRING_UTILITIES_H_
#define CORE_STRING_UTILITIES_H_

#include <sstream>
#include <string>
#include <vector>

namespace axom
{
namespace utilities
{
namespace string
{

/*!
 * \brief An invalid name string used in axom components
 *
 * This is an empty string that is available for axom components to compare
 * as an invalid name.
 */
static const std::string InvalidName;

/*!
 * \brief Tests whether a string ends with another string
 * https://stackoverflow.com/questions/20446201/how-to-check-if-string-ends-with-txt/20446257
 * \param [in] str string to be searched
 * \param [in] suffix string to be checked for
 * \return     boolean value whether suffix was found at the end of str
 */
inline bool endsWith(const std::string& str, const std::string& suffix)
{
  return str.size() >= suffix.size() &&
    0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

/*!
 * \brief Tests whether a string ends with a char
 * \param [in] str string to be searched
 * \param [in] suffix char to be checked for
 * \return     boolean value whether suffix was found at the end of str
 */
inline bool endsWith(const std::string& str, const char suffix)
{
  return str.size() >= 1 && str[str.size() - 1] == suffix;
}

/*!
 * \brief Tests whether a string starts with another string
 * \param [in] str string to be searched
 * \param [in] prefix string to be checked for
 * \return     boolean value whether prefix was found at start of str
 */
inline bool startsWith(const std::string& str, const std::string& prefix)
{
  return str.size() >= prefix.size() && 0 == str.compare(0, prefix.size(), prefix);
}

/*!
 * \brief Tests whether a string starts with a char
 * \param [in] str string to be searched
 * \param [in] prefix char to be checked for
 * \return     boolean value whether prefix was found at start of str
 */
inline bool startsWith(const std::string& str, const char prefix)
{
  return str.size() >= 1 && str[0] == prefix;
}

/*!
 * \brief Removes \a suffix from the end of \a str, if present
 * \param [in] str string to be searched
 * \param [in] suffix string to check for
 * \return     String with (one copy of) \a suffix removed from end, if it was present  
 */
inline std::string removeSuffix(const std::string& str, const std::string& suffix)
{
  if(!suffix.empty() && endsWith(str, suffix))
  {
    return str.substr(0, str.size() - suffix.size());
  }
  return str;
}

/*!
 * \brief Splits the given string based on the given delimiter
 * \param [in]  str       string to be tokenized
 * \param [in]  delimiter char to split string on
 * 
 * \return vector containing the found tokens
 */
std::vector<std::string> split(const std::string& str, const char delimiter);

/*!
 * \brief Converts a string to lowercase
 * \param [inout] str    string to be converted
 */
void toLower(std::string& str);

/*!
 * \brief Converts a string to uppercase
 * \param [inout] str    string to be converted
 */
void toUpper(std::string& str);

/**
 * \brief Retrieves the last "n" tokens of a string split with the specified delimiter
 * \param [in] input The string to split
 * \param [in] n The number of tokens to retrieve
 * \param [in] delim The delimiter to split with
 * 
 * \return A list of tokens (of size \a n )
 * \note If @a delim is not found in the input string, the full input string is returned
 * 
 * Splits a string starting from the end of the string into a maximum of \a n tokens
 */
std::vector<std::string> rsplitN(const std::string& input, const std::size_t n, const char delim);

/*!
*****************************************************************************
* \brief This function appends the prefix name to the ending name.
*
* \param [in] prefix The prefix string name.
* \param [in] name The ending string name.
* \param [in] delim The delimiting character inserted between the prefix and the name, if both are non-empty
*
* \return The appended string.  If @a prefix is empty, then just @a name is returned.
*****************************************************************************
*/
std::string appendPrefix(const std::string& prefix, const std::string& name, const char delim = '/');

/*!
*****************************************************************************
* \brief This function returns the tail of a string following a given prefix plus a delimiter
*
* \param [in] prefix The prefix of the name, to be removed.
* \param [in] name The full name.
* \param [in] delim The delimiting character "expected" to separate the prefix and the return value
*
* \return The extracted string.
*
* For example, removePrefix("foo", "foo/bar") will return "bar".
* If the prefix is empty or not found the full @a name string is returned.
*****************************************************************************
*/
std::string removePrefix(const std::string& prefix, const std::string& name, const char delim = '/');

/*!
*****************************************************************************
* \brief This function extracts the substring following the last instance
* of the delimiting character
*
* \param [in] str The string to extract from
* \param [in] delim The delimiting character
*
* \return The extracted string.
*****************************************************************************
*/
std::string removeBeforeDelimiter(const std::string& str, const char delim = '/');

/*!
*****************************************************************************
* \brief This function removes all instances of the substring from the target
* string
*
* \param [in] target The string to operate on
* \param [in] substr The string to remove
*
* \return The filtered string.
*****************************************************************************
*/
std::string removeAllInstances(const std::string& target, const std::string& substr);

/*!
*****************************************************************************
* \brief This function finds all instances of a substring and replaces it
* with another string for the target string
*
* \param [in] target The string to operate on
* \param [in] findstr The string to find
* \param [in] replacestr The string to replace with
*
* \return The modified string. If @a findstr and @a replacestr are the same or
* if @a findstr is an empty string, @a target is returned.
*****************************************************************************
*/
std::string replaceAllInstances(const std::string& target,
                                const std::string& findstr,
                                const std::string& replacestr);

}  // end namespace string
}  // end namespace utilities
}  // end namespace axom

#endif  //  CORE_STRING_UTILITIES_H_
