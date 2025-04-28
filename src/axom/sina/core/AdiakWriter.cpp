// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file AdiakWriter.cpp
 *
 * \brief   Implementation file for the Adiak Sina callback function.
 *
 ******************************************************************************
 */

#include "axom/sina/core/AdiakWriter.hpp"

#ifdef AXOM_USE_ADIAK

  #include <stdexcept>
  #include <utility>
  #include <iostream>
  #include <fstream>
  #include <stdexcept>

extern "C" {
  #include "adiak_tool.h"
}

  #include "axom/sina/core/ConduitUtil.hpp"
  #include "axom/sina/core/Record.hpp"
  #include "axom/sina/core/Datum.hpp"
  #include "axom/sina/core/Document.hpp"

namespace axom
{
namespace sina
{
namespace
{

/**
* Adiak has a much wider array of supported types than Sina. We will convert
* Adiak types to ones Sina understands; SinaType holds the possibilities.
**/
enum SinaType
{
  sina_scalar,
  sina_string,
  sina_list,
  sina_file,
  sina_unknown
};

/**
* Add a axom::sina::Datum object to a Record. These are the sina equivalent
* of an Adiak datapoint. Since we track slightly different info, this function
* harvests what it can and hands it off to the Record.
**/
template <typename T>
void addDatum(const std::string &name,
              T sina_safe_val,
              const std::vector<std::string> &tags,
              axom::sina::Record *record)
{
  axom::sina::Datum datum {sina_safe_val};
  datum.setTags(std::move(tags));
  record->add(name, datum);
}

/**
* Add a axom::sina::File object to our current Record. Adiak stores paths,
* which are essentially the same as Sina's idea of storing files.
**/
void addFile(const std::string &name, const std::string &uri, axom::sina::Record *record)
{
  // We don't care about type here, there's only one adiak type that acts as a file
  axom::sina::File file {uri};
  file.setTags(std::vector<std::string> {name});
  record->add(std::move(file));
}

/**
* Given an Adiak type, return its corresponding Sina type.
**/
SinaType findSinaType(adiak_datatype_t *t)
{
  switch(t->dtype)
  {
  case adiak_long:
  case adiak_ulong:
  case adiak_int:
  case adiak_uint:
  case adiak_double:
  case adiak_timeval:
    return sina_scalar;
  case adiak_date:
  case adiak_version:
  case adiak_string:
  case adiak_catstring:
    return sina_string;
  case adiak_path:
    return sina_file;
  case adiak_set:
  case adiak_tuple:
  case adiak_range:
  case adiak_list:
    return sina_list;
  case adiak_type_unset:
    return sina_unknown;
  default:
    return sina_unknown;
  }
}

/**
* Several Adiak types become what Sina views as a "scalar" (a double).
* Manage the conversions from various Adiak types to the final double
* representation
**/
double toScalar(adiak_value_t *val, adiak_datatype_t *adiak_type)
{
  switch(adiak_type->dtype)
  {
  case adiak_long:
  case adiak_ulong:
    return static_cast<double>(val->v_long);
  case adiak_int:
  case adiak_uint:
    return static_cast<double>(val->v_int);
  case adiak_double:
    return val->v_double;
  case adiak_timeval:
  {
    struct timeval *tval = static_cast<struct timeval *>(val->v_ptr);
    return static_cast<double>(tval->tv_sec) + (static_cast<double>(tval->tv_usec) / 1000000.0);
  }
  // None of the rest of these should ever be reachable, so special error message
  case adiak_date:
  case adiak_version:
  case adiak_string:
  case adiak_catstring:
  case adiak_path:
  case adiak_set:
  case adiak_tuple:
  case adiak_range:
  case adiak_list:
  case adiak_type_unset:
  {
    std::string msg("Logic error, contact maintainer: Adiak-to-Sina double converter given ");
    char *s = adiak_type_to_string(adiak_type, 1);
    msg += s;
    free(s);
    throw std::runtime_error(msg);
  }
  default:
    throw std::runtime_error(
      "Adiak-to-Sina double converter given something not convertible to "
      "double");
  }
}

/**
* Some Adiak types become what Sina views as a string.
* Manage the conversions from various Adiak types to said string.
**/
std::string toString(adiak_value_t *val, adiak_datatype_t *adiak_type)
{
  switch(adiak_type->dtype)
  {
  case adiak_date:
  {
    char datestr[512];
    signed long seconds_since_epoch = static_cast<signed long>(val->v_long);
    struct tm *loc = localtime(&seconds_since_epoch);
    strftime(datestr, sizeof(datestr), "%a, %d %b %Y %T %z", loc);
    return static_cast<std::string>(datestr);
  }
  case adiak_catstring:
  case adiak_version:
  case adiak_string:
  case adiak_path:
    return std::string(static_cast<char *>(val->v_ptr));
  case adiak_long:
  case adiak_ulong:
  case adiak_int:
  case adiak_uint:
  case adiak_double:
  case adiak_timeval:
  case adiak_set:
  case adiak_tuple:
  case adiak_range:
  case adiak_list:
  case adiak_type_unset:
  {
    std::string msg("Logic error, contact maintainer: Adiak-to-Sina string converter given ");
    char *s = adiak_type_to_string(adiak_type, 1);
    msg += s;
    free(s);
    throw std::runtime_error(msg);
  }
  default:
    throw std::runtime_error(
      "Adiak-to-Sina string converter given something not convertible to "
      "string");
  }
}

/**
* Some Adiak types become a list of some form. Sina, being concerned
* with queries and visualization, only handles lists that are all scalars
* or all strings. Manage conversions from various Adiak list types that
* contain scalars to a simple list (vector) of scalars.
**/
std::vector<double> toScalarList(adiak_value_t *subvals, adiak_datatype_t *t)
{
  std::vector<double> sina_safe_list;
  for(int i = 0; i < t->num_elements; i++)
  {
    sina_safe_list.emplace_back(toScalar(subvals + i, t->subtype[0]));
  }
  return sina_safe_list;
}

/**
* Partner method to toScalarList, invoked when the children of an adiak list
* type are strings (according to Sina).
**/
std::vector<std::string> toStringList(adiak_value_t *subvals, adiak_datatype_t *t)
{
  std::vector<std::string> sina_safe_list;
  for(int i = 0; i < t->num_elements; i++)
  {
    sina_safe_list.emplace_back(toString(subvals + i, t->subtype[0]));
  }
  return sina_safe_list;
}

}  // namespace

void adiakSinaCallback(const char *name,
                       adiak_category_t,
                       const char *subcategory,
                       adiak_value_t *val,
                       adiak_datatype_t *adiak_type,
                       void *void_record)
{
  const SinaType sina_type = findSinaType(adiak_type);
  axom::sina::Record *record = static_cast<axom::sina::Record *>(void_record);
  std::vector<std::string> tags;
  if(subcategory && subcategory[0] != '\0')
  {
    tags.emplace_back(subcategory);
  }
  switch(sina_type)
  {
  case sina_unknown:
    // If we don't know what it is, we can't store it, so as above...
    throw std::runtime_error("Unknown Adiak type cannot be added to Sina record.");
  case sina_scalar:
  {
    char *s = adiak_type_to_string(adiak_type, 1);
    tags.emplace_back(s);
    free(s);
    addDatum(name, toScalar(val, adiak_type), tags, record);
    break;
  }
  case sina_string:
  {
    char *s = adiak_type_to_string(adiak_type, 1);
    tags.emplace_back(s);
    free(s);
    addDatum(name, toString(val, adiak_type), tags, record);
    break;
  }
  case sina_file:
    addFile(name, toString(val, adiak_type), record);
    break;
  case sina_list:
  {
    // Sina doesn't really know/care the difference between list, tuple, set
    // Further simplification: everything has to be the same type
    // Even further simplification: nothing nested. In the future, depth>1 lists
    // should be sent to user_defined
    adiak_value_t *subvals = static_cast<adiak_value_t *>(val->v_ptr);
    SinaType list_type = findSinaType(adiak_type->subtype[0]);
    char *s = adiak_type_to_string(adiak_type->subtype[0], 1);
    tags.emplace_back(s);
    free(s);
    switch(list_type)
    {
    case sina_string:
      addDatum(name, toStringList(subvals, adiak_type), tags, record);
      break;
    // Weird case wherein we're given a list of filenames, which we can somewhat manage
    case sina_file:
      for(int i = 0; i < adiak_type->num_elements; i++)
      {
        addFile(name, toString(subvals + i, adiak_type->subtype[0]), record);
      }
      break;
    case sina_scalar:
      addDatum(name, toScalarList(subvals, adiak_type), tags, record);
      break;
    case sina_unknown:
      throw std::runtime_error(
        "Type must not be unknown for list entries to be added to a Sina "
        "record");
    case sina_list:
      throw std::runtime_error(
        "Lists must not be nested for list entries to be added to a Sina "
        "record");
    default:
      throw std::runtime_error("Type must be set for list entries to be added to a Sina record");
    }
  }
  }
}

}  // namespace sina
}  // namespace axom

#endif  // AXOM_USE_ADIAK
