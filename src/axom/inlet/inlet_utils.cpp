// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/inlet_utils.hpp"

namespace axom
{
namespace inlet
{
void setWarningFlag(axom::sidre::Group* root)
{
  if(!root->hasView("warningFlag"))
  {
    root->createViewScalar("warningFlag", 1);
  }
}

void setFlag(axom::sidre::Group& target, axom::sidre::Group& root, const std::string& flag, bool value)
{
  const std::int8_t bval = value ? 1 : 0;
  if(target.hasView(flag))
  {
    auto flagView = target.getView(flag);
    if(flagView->getData<std::int8_t>() != bval)
    {
      const std::string msg =
        fmt::format("[Inlet] '{0}' value has already been defined for: {1}", flag, target.getName());

      SLIC_WARNING(msg);
      setWarningFlag(&root);
    }
  }
  else
  {
    if(value)
    {
      target.createViewScalar(flag, bval);
    }
    else
    {
      target.createViewScalar(flag, bval);
    }
  }
}

bool checkFlag(const axom::sidre::Group& target, axom::sidre::Group& root, const std::string& flag)
{
  if(!target.hasView(flag))
  {
    return false;
  }
  const axom::sidre::View* valueView = target.getView(flag);
  const std::int8_t intValue = valueView->getScalar();
  if(intValue < 0 || intValue > 1)
  {
    const std::string msg = fmt::format(
      "[Inlet] Invalid integer value stored in "
      " boolean value named {0} for flag '{1}'",
      target.getName(),
      flag);
    SLIC_WARNING(msg);
    setWarningFlag(&root);
    return static_cast<bool>(intValue);
  }

  return static_cast<bool>(intValue);
}

bool verifyRequired(const axom::sidre::Group& target,
                    const bool condition,
                    const std::string& type,
                    std::vector<VerificationError>* errors)
{
  // Assume that it wasn't found
  ReaderResult status = ReaderResult::NotFound;
  if(target.hasView("retrieval_status"))
  {
    status =
      static_cast<ReaderResult>(static_cast<int>(target.getView("retrieval_status")->getData()));
  }

  if(target.hasView("required"))
  {
    std::int8_t required = target.getView("required")->getData();
    // If it wasn't found at all, it's only an error if the object was required and not provided
    // The retrieval_status will typically (but not always) be NotFound in these cases, but that
    // information isn't needed here unless it's a collection group - empty collections are permissible,
    // so they shouldn't impede verification if they existed but were empty (hence Success check)
    if(required && !condition &&
       (!isCollectionGroup(target.getPathName()) || status != ReaderResult::Success))
    {
      const std::string msg = fmt::format(
        "[Inlet] Required {0} not "
        "specified: {1}",
        type,
        target.getPathName());
      INLET_VERIFICATION_WARNING(target.getPathName(), msg, errors);
      return false;
    }
  }

  // If it was the wrong type or part of a non-homogeneous array, it's always an error,
  // even if the object wasn't marked as required
  if(status == ReaderResult::WrongType || status == ReaderResult::NotHomogeneous)
  {
    const std::string reason =
      (status == ReaderResult::WrongType) ? "of the wrong type" : "not homogeneous";
    const std::string msg =
      fmt::format("[Inlet] {0} '{1}' was {2}", type, target.getPathName(), reason);
    INLET_VERIFICATION_WARNING(target.getPathName(), msg, errors);
    return false;
  }
  return true;
}

void markAsStructCollection(axom::sidre::Group& target)
{
  if(target.hasView(detail::STRUCT_COLLECTION_FLAG))
  {
    // This flag should only ever be one, so we verify that and error otherwise
    const sidre::View* flag = target.getView(detail::STRUCT_COLLECTION_FLAG);
    SLIC_ERROR_IF(!flag->isScalar(),
                  fmt::format("[Inlet] Struct collection flag of group '{0}' was not a scalar",
                              target.getName()));
    const std::int8_t value = flag->getScalar();
    SLIC_ERROR_IF(value != 1,
                  fmt::format("[Inlet] Struct collection flag of group '{0}' "
                              "had a value other than 1",
                              target.getName()));
  }
  else
  {
    target.createViewScalar(detail::STRUCT_COLLECTION_FLAG, static_cast<std::int8_t>(1));
  }
}

void markRetrievalStatus(axom::sidre::Group& target, const ReaderResult result)
{
  if(!target.hasView("retrieval_status"))
  {
    target.createViewScalar("retrieval_status", static_cast<int>(result));
  }
}

ReaderResult collectionRetrievalResult(const bool contains_other_type,
                                       const bool contains_requested_type)
{
  // First check if the collection was entirely the wrong type
  if(contains_other_type && !contains_requested_type)
  {
    return ReaderResult::WrongType;
  }
  // Then check if some values were of the correct type, but others weren't
  else if(contains_other_type)
  {
    return ReaderResult::NotHomogeneous;
  }
  // Otherwise we mark it as successful - having just an empty collection
  // counts as success
  return ReaderResult::Success;
}

}  // namespace inlet
}  // namespace axom
