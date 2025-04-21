// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file DataHolder.cpp
 *
 * \brief   Implementation file for Sina DataHolder class
 *
 ******************************************************************************
 */

#include "axom/sina/core/DataHolder.hpp"

#include "axom/sina/core/ConduitUtil.hpp"
#include "axom/sina/core/Datum.hpp"

#include <stdexcept>

namespace
{

char const DATA_FIELD[] = "data";
char const CURVE_SETS_FIELD[] = "curve_sets";
char const LIBRARY_DATA_FIELD[] = "library_data";
char const USER_DEFINED_FIELD[] = "user_defined";

}  // namespace

namespace axom
{
namespace sina
{

void DataHolder::add(std::string name, Datum datum)
{
  auto existing = data.find(name);
  if(existing == data.end())
  {
    data.emplace(std::make_pair(std::move(name), datum));
  }
  else
  {
    existing->second = datum;
  }
}

void DataHolder::add(CurveSet curveSet)
{
  auto name = curveSet.getName();
  auto existing = curveSets.find(name);
  if(existing == curveSets.end())
  {
    curveSets.emplace(name, std::move(curveSet));
  }
  else
  {
    existing->second = std::move(curveSet);
  }
}

std::shared_ptr<DataHolder> DataHolder::addLibraryData(std::string const &name)
{
  auto existing = libraryData.find(name);
  if(existing == libraryData.end())
  {
    libraryData.emplace(name, std::make_shared<DataHolder>());
  }
  else
  {
    existing->second = std::make_shared<DataHolder>();
  }
  return libraryData.at(name);
}

std::shared_ptr<DataHolder> DataHolder::addLibraryData(std::string const &name,
                                                       conduit::Node existingLibraryData)
{
  auto existing = libraryData.find(name);
  if(existing == libraryData.end())
  {
    libraryData.emplace(name, std::make_shared<DataHolder>(existingLibraryData));
  }
  else
  {
    existing->second = std::make_shared<DataHolder>(existingLibraryData);
  }
  return libraryData.at(name);
}

void DataHolder::setUserDefinedContent(conduit::Node userDefined_)
{
  userDefined = std::move(userDefined_);
}

conduit::Node DataHolder::toNode(CurveSet::CurveOrder curveOrder) const
{
  conduit::Node asNode;
  asNode.set(conduit::DataType::object());
  if(!libraryData.empty())
  {
    //Loop through vector of data and append Json
    conduit::Node libRef;
    for(auto &lib : libraryData)
    {
      libRef.add_child(lib.first) = lib.second->toNode(curveOrder);
    }
    asNode[LIBRARY_DATA_FIELD] = libRef;
  }
  if(!curveSets.empty())
  {
    conduit::Node curveSetsNode;
    for(auto &entry : curveSets)
    {
      curveSetsNode.add_child(entry.first) = entry.second.toNode(curveOrder);
    }
    asNode[CURVE_SETS_FIELD] = curveSetsNode;
  }
  if(!data.empty())
  {
    //Loop through vector of data and append Json
    conduit::Node datumRef;
    for(auto &datum : data)
    {
      datumRef.add_child(datum.first) = datum.second.toNode();
    }
    asNode[DATA_FIELD] = datumRef;
  }
  if(!userDefined.dtype().is_empty())
  {
    asNode[USER_DEFINED_FIELD] = userDefined;
  }
  return asNode;
}

DataHolder::DataHolder(conduit::Node const &asNode)
{
  if(asNode.has_child(DATA_FIELD))
  {
    auto dataIter = asNode[DATA_FIELD].children();
    //Loop through DATA_FIELD objects and add them to data:
    while(dataIter.has_next())
    {
      auto &namedDatum = dataIter.next();
      data.emplace(std::make_pair(dataIter.name(), Datum(namedDatum)));
    }
  }
  if(asNode.has_child(CURVE_SETS_FIELD))
  {
    auto curveSetsIter = asNode[CURVE_SETS_FIELD].children();
    while(curveSetsIter.has_next())
    {
      auto &curveSetNode = curveSetsIter.next();
      std::string name = curveSetsIter.name();
      CurveSet cs {name, curveSetNode};
      curveSets.emplace(std::make_pair(std::move(name), std::move(cs)));
    }
  }
  if(asNode.has_child(LIBRARY_DATA_FIELD))
  {
    auto libraryIter = asNode[LIBRARY_DATA_FIELD].children();
    while(libraryIter.has_next())
    {
      auto &libraryDataNode = libraryIter.next();
      std::string name = libraryIter.name();
      libraryData.emplace(
        std::make_pair(std::move(name), std::make_shared<DataHolder>(libraryDataNode)));
    }
  }
  if(asNode.has_child(USER_DEFINED_FIELD))
  {
    userDefined = asNode[USER_DEFINED_FIELD];
    if(!userDefined.dtype().is_object())
    {
      throw std::invalid_argument("user_defined must be an object Node");
    }
  }
}

}  // namespace sina
}  // namespace axom
