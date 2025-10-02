// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file File.cpp
 *
 * \brief   Implementation file for Sina File class
 *
 ******************************************************************************
 */

#include "axom/sina/core/File.hpp"
#include "axom/sina/core/ConduitUtil.hpp"

#include <stdexcept>
#include <utility>
#include <sstream>

#include "conduit.hpp"

namespace axom
{
namespace sina
{

namespace
{
char const MIMETYPE_KEY[] = "mimetype";
char const FILE_TYPE_NAME[] = "File";
char const TAGS_KEY[] = "tags";
}  // namespace

File::File(std::string uri_) : uri {std::move(uri_)} { }

File::File(std::string uri_, conduit::Node const &asNode)
  : uri {std::move(uri_)}
  , mimeType {getOptionalString(MIMETYPE_KEY, asNode, FILE_TYPE_NAME)}
{
  if(asNode.has_child(TAGS_KEY))
  {
    auto tagsIter = asNode[TAGS_KEY].children();
    while(tagsIter.has_next())
    {
      auto &tag = tagsIter.next();
      if(tag.dtype().is_string())
        tags.emplace_back(tag.as_string());
      else
      {
        std::ostringstream message;
        message << "The optional field '" << TAGS_KEY << "' must be an array of strings. Found '"
                << tag.dtype().name() << "' instead.";
        throw std::invalid_argument(message.str());
      }
    }
  }
}

void File::setMimeType(std::string mimeType_) { File::mimeType = std::move(mimeType_); }

void File::setTags(std::vector<std::string> tags_) { File::tags = std::move(tags_); }

conduit::Node File::toNode() const
{
  conduit::Node asNode;
  if(!mimeType.empty())
  {
    asNode[MIMETYPE_KEY] = mimeType;
  }
  if(tags.size() > 0)
  {
    std::vector<std::string> tags_copy(tags);
    addStringsToNode(asNode, TAGS_KEY, tags_copy);
  }
  return asNode;
}

}  // namespace sina
}  // namespace axom
