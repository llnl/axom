// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <memory>
#include <utility>

#include "axom/sidre.hpp"
#include "axom/fmt.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/core/Path.hpp"

#ifndef INLET_UTILS_HPP
  #define INLET_UTILS_HPP

namespace axom
{
namespace inlet
{
enum class ReaderResult
{
  Success,         // Found with no issue
  NotFound,        // Path does not exist in the input file
  NotHomogeneous,  // Found, but elements of other type exist
  WrongType        // Found, but item at specified path was not of requested type
};

/*!
 *****************************************************************************
 * \brief Information on an Inlet verification error
 *****************************************************************************
 */
struct VerificationError
{
  /// \brief The path to the container/field/function with the error
  const axom::Path path;
  /// \brief The error message
  const std::string message;
  /// \brief Returns whether a given substring is present in the error message
  bool messageContains(const std::string substr) const
  {
    return message.find(substr) != std::string::npos;
  }
};

/*!
 *****************************************************************************
 * \brief Utility macro for selecting between logging to SLIC and logging
 * to a list of errors
 * \param path The path within the input file to warn on
 * \param msg The warning message
 * \param errs The list of errors, must be of type \p std::vector<VerificationError>*
 *****************************************************************************
 */
  #define INLET_VERIFICATION_WARNING(path, msg, errs) \
    if(errs)                                          \
    {                                                 \
      errs->push_back({axom::Path {path}, msg});      \
    }                                                 \
    else                                              \
    {                                                 \
      SLIC_WARNING(msg);                              \
    }

/*!
*****************************************************************************
* \brief This function is used to mark if anything went wrong during the 
* defining phase of inlet so verify() will properly fail.
*
* \param [in] root Pointer to the Sidre Root Group where the warning flag 
* will be set.
*****************************************************************************
*/
void setWarningFlag(axom::sidre::Group* root);

/*!
*****************************************************************************
* \brief This function is used to add a flag to the Inlet object
* corresponding to the provided Sidre group
*
* \param [in] target Reference to the Sidre group to set the required 
* status of
* \param [in] root Reference to the Sidre Root Group where the warning flag 
* will be set on failure
* \param [in] flag The name of the flag to set
* \param [in] value The value of the flag
*****************************************************************************
*/
void setFlag(axom::sidre::Group& target, axom::sidre::Group& root, const std::string& flag, bool value);

/*!
*****************************************************************************
* \brief This function is used to determine the value of a flag for the
* Inlet object corresponding to the provided Sidre group
*
* \param [in] target Reference to the Sidre group to check the required 
* status of
* \param [in] root Reference to the Sidre Root Group where the warning flag 
* will be set on failure
* \param [in] flag The name of the flag to check
* \return The value of the flag
*****************************************************************************
*/
bool checkFlag(const axom::sidre::Group& target, axom::sidre::Group& root, const std::string& flag);

/*!
*****************************************************************************
* \brief This function is used to verify the required-ness of the Inlet object
* corresponding to the provided Sidre group
*
* \param [in] target Reference to the Sidre group to verify the required-ness of
* \param [in] condition The condition that must be true if the object is required
* \param [in] type The type of the object as a string, for use in the warning message
* \param [in] errors An optional vector of errors to append to in the case
* of verification failure
* 
* \return False if the object was required but \p condition was false, True otherwise
* \post If the function returns False, a warning message will be emitted
*****************************************************************************
*/
bool verifyRequired(const axom::sidre::Group& target,
                    const bool condition,
                    const std::string& type,
                    std::vector<VerificationError>* errors = nullptr);

namespace detail
{
/*!
  *******************************************************************************
  * Names of the internal collection data and collection index groups/fields
  * used for managing arrays/dictionaries
  *******************************************************************************
  */
const std::string COLLECTION_GROUP_NAME = "_inlet_collection";
const std::string COLLECTION_INDICES_NAME = "_inlet_collection_indices";
const std::string STRUCT_COLLECTION_FLAG = "_inlet_struct_collection";
const std::string REQUIRED_FLAG = "required";
const std::string STRICT_FLAG = "strict";
}  // namespace detail

/*!
*****************************************************************************
* \brief Determines whether a Container is a collection group
*
* \param [in] name The name of the container
*****************************************************************************
*/
inline bool isCollectionGroup(const std::string& name)
{
  return axom::utilities::string::endsWith(name, detail::COLLECTION_GROUP_NAME);
}

/*!
*****************************************************************************
* \brief Marks the sidre::Group as a "struct collection" by adding a
* corresponding flag to the group
*
* \param [inout] target The group to tag
*****************************************************************************
*/
void markAsStructCollection(axom::sidre::Group& target);

/*!
*****************************************************************************
* \brief Adds a ReaderResult to a sidre::Group corresponding to an inlet
* object
*
* \param [inout] target The group to tag
* \param [in] result The retrieval result
*****************************************************************************
*/
void markRetrievalStatus(axom::sidre::Group& target, const ReaderResult result);

/*!
*****************************************************************************
* \brief Returns the corresponding retrieval result for a collection depending
* on whether the collection contained any elements of the requested or of
* other type
*
* \param [in] contains_other_type Whether any collection elements were of type
* other than the requested type
* \param [in] contains_requested_type Whether the collection of requested type
* was not empty, i.e., if any elements of the requested type were present
*****************************************************************************
*/
ReaderResult collectionRetrievalResult(const bool contains_other_type,
                                       const bool contains_requested_type);

}  // namespace inlet
}  // namespace axom

#endif
