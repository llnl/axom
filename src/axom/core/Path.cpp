// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/Path.hpp"

#include <iostream>
#include <iterator>

#include "axom/fmt.hpp"

#include "axom/core/utilities/StringUtilities.hpp"

namespace axom
{
Path::Path(const std::string& path, const char delim) : m_delim(delim)
{
  size_t first_position = 0;  // position of first non-delimiter char in path

  if(axom::utilities::string::startsWith(path, m_delim))
  {
    m_leading_delim = true;
    first_position = 1;
  }

  // Check if the path has more than one component
  if(path.find(delim, first_position) != std::string::npos)
  {
    m_components = utilities::string::split(path, delim);

    // Remove empty parts
    m_components.erase(std::remove_if(m_components.begin(),
                                      m_components.end(),
                                      [](std::string const& s) { return s.empty(); }),
                       m_components.end());
  }
  else if(!path.empty())
  {
    if(m_leading_delim)
    {
      if(path.size() > 1)
      {
        m_components.push_back(path.substr(1));
      }
    }
    else
    {
      m_components.push_back(path);
    }
  }
}

Path Path::join(std::initializer_list<Path> paths, const char delim)
{
  Path result;
  result.m_delim = delim;
  for(const auto& path : paths)
  {
    // Insert everything in order
    std::move(path.m_components.begin(),
              path.m_components.end(),
              std::back_inserter(result.m_components));
  }

  if(paths.size() != 0 && (*paths.begin()).m_leading_delim)
  {
    result.m_leading_delim = true;
  }
  else
  {
    result.m_leading_delim = false;
  }

  return result;
}

Path::operator std::string() const
{
  return fmt::format("{0}{1}",
                     m_leading_delim ? std::string(1, m_delim) : "",
                     fmt::join(m_components, std::string(1, m_delim)));
}

Path Path::parent() const
{
  Path result(*this);
  if(!result.m_components.empty())  // [[likely]] - uncomment when C++20 available
  {
    result.m_components.pop_back();
  }
  return result;
}

std::string Path::baseName() const
{
  if(m_components.empty())
  {
    return {};
  }
  return m_components.back();
}

std::string Path::dirName() const { return static_cast<std::string>(parent()); }

std::pair<std::string, std::string> Path::split() const
{
  return std::make_pair(dirName(), baseName());
}

bool operator==(const Path& lhs, const Path& rhs)
{
  return static_cast<std::string>(lhs) == static_cast<std::string>(rhs);
}

}  // end namespace axom
