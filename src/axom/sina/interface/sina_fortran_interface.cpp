// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <string.h>

#include "axom/core.hpp"
#include "axom/sina/interface/sina_fortran_interface.h"
#include "axom/sina/core/Document.hpp"
#include "axom/sina/core/Record.hpp"
#include "axom/sina/core/CurveSet.hpp"
#include <cstring>

std::vector<std::unique_ptr<axom::sina::Record>> sinaRecordsList;
axom::sina::Document *sina_document;
char default_record_type[25] = "fortran_code_output";

// Helper function to check if modifications are allowed
inline bool can_modify_records()
{
  if(sina_document != nullptr)
  {
    std::cerr << "ERROR: Cannot modify records after document has been created. "
              << "Call sina_write_document with preserve=0 first" << std::endl;
    return false;
  }
  return true;
}
// Helper function to clean Fortran strings
inline std::string fortran_to_cpp_string(const char *str, int str_len)
{
  if(!str || str_len <= 0) return "";

  // Find null terminator within the given length
  int actual_len = 0;
  for(int i = 0; i < str_len; ++i)
  {
    if(str[i] == '\0')
    {
      actual_len = i;
      break;
    }
  }
  if(actual_len == 0 && str_len > 0) actual_len = str_len;

  return std::string(str, actual_len);
}

extern "C" void sina_set_default_record_type_(char *record_type)
{
  strcpy(default_record_type, record_type);
}

extern "C" char *Get_File_Extension(char *input_fn)
{
  char *ext = strrchr(input_fn, '.');
  if(!ext)
  {
    return (new char[1] {'\0'});
  }
  return (ext + 1);
}

extern "C" void sina_create_record_(char *recID, char *recType, int recId_length, int recType_length)
{
  if(!can_modify_records()) return;

  // Clean up recID
  std::string id_str = fortran_to_cpp_string(recID, recId_length);

  axom::sina::ID id {id_str, axom::sina::IDType::Global};

  // Clean up recType
  std::string type_str = fortran_to_cpp_string(recType, recType_length);

  if(type_str.empty())
  {
    type_str = default_record_type;  // Use default if empty
  }

  // Pass the string object, not a pointer - Record should copy it
  sinaRecordsList.emplace_back(std::make_unique<axom::sina::Record>(id, type_str));
}

extern "C" axom::sina::Record *Sina_Get_Record(char *recId = NULL)
{
  if(recId == NULL || recId[0] == '\0')
  {
    std::unique_ptr<axom::sina::Record> const &myRecord = sinaRecordsList.front();
    return myRecord.get();
  }
  else
  {
    axom::sina::ID id {recId, axom::sina::IDType::Global};
    for(const std::unique_ptr<axom::sina::Record> &myRecord : sinaRecordsList)
    {
      const char *current_id_str = myRecord->getId().getId().c_str();
      // Compare the input C-string (recId) with the current record's C-string ID
      if(strcmp(recId, current_id_str) == 0)
      {
        return myRecord.get();
      }
    }
    // Didn't match we will return a new record
    char empty_str[] = "";
    sina_create_record_(recId, empty_str, strlen(recId), 0);
    return sinaRecordsList.back().get();
  }
  return nullptr;
}

extern "C" void sina_add_logical_(char *key,
                                  bool *value,
                                  char *units,
                                  char *tags,
                                  char *recId,
                                  int key_len,
                                  int units_len,
                                  int tags_len,
                                  int recId_len)
{
  if(!can_modify_records()) return;

  std::string recId_str = fortran_to_cpp_string(recId, recId_len);
  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));

  std::string key_name = fortran_to_cpp_string(key, key_len);

  // Convert Fortran logical properly (Fortran .true. might be -1)
  axom::sina::Datum datum {*value ? 1.0 : 0.0};

  std::string key_units = fortran_to_cpp_string(units, units_len);
  if(!key_units.empty())
  {
    datum.setUnits(key_units);
  }

  std::string tags_str = fortran_to_cpp_string(tags, tags_len);
  if(!tags_str.empty())
  {
    datum.setTags({tags_str});
  }

  sina_record->add(key_name, datum);
}

extern "C" void sina_add_long_(char *key,
                               long long int *value,
                               char *units,
                               char *tags,
                               char *recId,
                               int key_len,
                               int units_len,
                               int tags_len,
                               int recId_len)
{
  if(!can_modify_records()) return;

  std::string recId_str = fortran_to_cpp_string(recId, recId_len);
  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));

  std::string key_name = fortran_to_cpp_string(key, key_len);
  axom::sina::Datum datum {static_cast<double>(*value)};

  std::string key_units = fortran_to_cpp_string(units, units_len);
  if(!key_units.empty())
  {
    datum.setUnits(key_units);
  }

  std::string tags_str = fortran_to_cpp_string(tags, tags_len);
  if(!tags_str.empty())
  {
    datum.setTags({tags_str});
  }

  sina_record->add(key_name, datum);
}

extern "C" void sina_add_int_(char *key,
                              int *value,
                              char *units,
                              char *tags,
                              char *recId,
                              int key_len,
                              int units_len,
                              int tags_len,
                              int recId_len)
{
  if(!can_modify_records()) return;

  std::string recId_str = fortran_to_cpp_string(recId, recId_len);
  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));

  std::string key_name = fortran_to_cpp_string(key, key_len);
  axom::sina::Datum datum {static_cast<double>(*value)};

  std::string key_units = fortran_to_cpp_string(units, units_len);
  if(!key_units.empty())
  {
    datum.setUnits(key_units);
  }

  std::string tags_str = fortran_to_cpp_string(tags, tags_len);
  if(!tags_str.empty())
  {
    datum.setTags({tags_str});
  }

  sina_record->add(key_name, datum);
}

extern "C" void sina_add_double_(char *key,
                                 double *value,
                                 char *units,
                                 char *tags,
                                 char *recId,
                                 int key_len,
                                 int units_len,
                                 int tags_len,
                                 int recId_len)
{
  if(!can_modify_records()) return;

  std::string recId_str = fortran_to_cpp_string(recId, recId_len);
  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));

  std::string key_name = fortran_to_cpp_string(key, key_len);
  axom::sina::Datum datum {*value};

  std::string key_units = fortran_to_cpp_string(units, units_len);
  if(!key_units.empty())
  {
    datum.setUnits(key_units);
  }

  std::string tags_str = fortran_to_cpp_string(tags, tags_len);
  if(!tags_str.empty())
  {
    datum.setTags({tags_str});
  }

  sina_record->add(key_name, datum);
}

extern "C" void sina_add_float_(char *key,
                                float *value,
                                char *units,
                                char *tags,
                                char *recId,
                                int key_len,
                                int units_len,
                                int tags_len,
                                int recId_len)
{
  if(!can_modify_records()) return;

  std::string recId_str = fortran_to_cpp_string(recId, recId_len);
  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));

  std::string key_name = fortran_to_cpp_string(key, key_len);
  axom::sina::Datum datum {*value};

  std::string key_units = fortran_to_cpp_string(units, units_len);
  if(!key_units.empty())
  {
    datum.setUnits(key_units);
  }

  std::string tags_str = fortran_to_cpp_string(tags, tags_len);
  if(!tags_str.empty())
  {
    datum.setTags({tags_str});
  }

  sina_record->add(key_name, datum);
}

// Fix for sina_add_string_ - remove value_len parameter since it's not in the header
extern "C" void sina_add_string_(char *key,
                                 char *value,
                                 char *units,
                                 char *tags,
                                 char *recId,
                                 int key_len,
                                 int units_len,
                                 int tags_len,
                                 int recId_len)
{
  AXOM_UNUSED_VAR(units_len);
  AXOM_UNUSED_VAR(tags_len);

  if(!can_modify_records()) return;

  std::string recId_str = fortran_to_cpp_string(recId, recId_len);
  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));

  std::string key_name = fortran_to_cpp_string(key, key_len);
  // Since we don't have value_len, we need to figure out the length
  // Assuming value follows same pattern as other string params with hidden length
  // We'll use strlen as fallback since we don't have the length parameter
  std::string key_value(value);  // This will use the null terminator

  axom::sina::Datum datum {key_value};

  std::string key_units(units);
  if(!key_units.empty())
  {
    datum.setUnits(key_units);
  }

  std::string tags_str(tags);
  if(!tags_str.empty())
  {
    datum.setTags({tags_str});
  }

  sina_record->add(key_name, datum);
}

// Fix for sina_add_file_ - cast to non-const char*
extern "C" void sina_add_file_(char *filename,
                               char *mime_type,
                               char *recId,
                               int file_len,
                               int mime_len,
                               int recId_len)
{
  if(!can_modify_records()) return;

  std::string recId_str = fortran_to_cpp_string(recId, recId_len);
  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));

  std::string filename_str = fortran_to_cpp_string(filename, file_len);
  std::string mime_type_str = fortran_to_cpp_string(mime_type, mime_len);

  axom::sina::File my_file {filename_str};

  bool is_effectively_empty = mime_type_str.find_first_not_of('\0') == std::string::npos;

  if(!is_effectively_empty)
  {
    my_file.setMimeType(mime_type_str);
  }
  else
  {
    std::string ext = Get_File_Extension(const_cast<char *>(filename_str.c_str()));
    my_file.setMimeType(ext);
  }

  if(sina_record)
  {
    sina_record->add(my_file);
  }
}

extern "C" void sina_write_document_all_args_(char *input_fn,
                                              int *protocol,
                                              int *preserve,
                                              int *mergeProtocol)
{
  // Create the document if needed
  if(sina_document == nullptr)
  {
    sina_document = new axom::sina::Document();

    // Move all records into the document
    for(auto &uniquePtr : sinaRecordsList)
    {
      if(uniquePtr)
      {
        sina_document->add(std::move(uniquePtr));
      }
    }
    sinaRecordsList.clear();
  }

  // input_fn should be null-terminated from Fortran's make_cstring
  std::string filename(input_fn);
  axom::sina::Protocol proto = static_cast<axom::sina::Protocol>(*protocol);

  // Save everything
  axom::sina::appendDocument(*sina_document, filename.c_str(), *mergeProtocol, proto);

  // Do we want to bring it back?
  if(*preserve == 0)
  {
    delete sina_document;
    sina_document = nullptr;
    sinaRecordsList.clear();
  }
}

extern "C" void sina_write_document_noprotocol_nopreserve_nomerge_(char *input_fn)
{
  int default_protocol = static_cast<int>(axom::sina::Protocol::AUTO_DETECT);
  int default_merge_protocol = 0;
  int default_preserve = 0;

  sina_write_document_all_args_(input_fn, &default_protocol, &default_preserve, &default_merge_protocol);
}

extern "C" void sina_write_document_protocol_nopreserve_nomerge_(char *input_fn, int *protocol)
{
  int default_merge_protocol = 0;
  int default_preserve = 0;

  sina_write_document_all_args_(input_fn, protocol, &default_preserve, &default_merge_protocol);
}

extern "C" void sina_write_document_protocol_preserve_nomerge_(char *input_fn,
                                                               int *protocol,
                                                               int *preserve)
{
  int default_merge_protocol = 0;

  sina_write_document_all_args_(input_fn, protocol, preserve, &default_merge_protocol);
}

extern "C" void sina_add_curveset_(char *name, char *recId, int name_len, int recId_len)
{
  if(!can_modify_records()) return;

  std::string recId_str = fortran_to_cpp_string(recId, recId_len);
  std::string name_str = fortran_to_cpp_string(name, name_len);

  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));
  if(sina_record)
  {
    axom::sina::CurveSet cs {name_str};
    sina_record->add(cs);
  }
}

extern "C" void sina_add_curve_double_(char *curveset_name,
                                       char *curve_name,
                                       double *values,
                                       int *n,
                                       int *independent,
                                       char *recId)
{
  if(!can_modify_records()) return;

  std::string recId_str(recId);
  std::string curveset_str(curveset_name);
  std::string curvename_str(curve_name);
  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));

  if(sina_record)
  {
    axom::sina::Curve curve {curvename_str, values, static_cast<size_t>(*n)};
    auto &curvesets = sina_record->getCurveSets();
    if(curvesets.find(curveset_str) == curvesets.end())
    {
      // Create curveset directly
      axom::sina::CurveSet cs {curveset_str};
      if(*independent != 0)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
    }
    else
    {
      // Curveset already exists
      axom::sina::CurveSet cs = curvesets.at(curveset_str);
      if(*independent != 0)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
    }
  }
}

extern "C" void sina_add_curve_float_(char *curveset_name,
                                      char *curve_name,
                                      float *values,
                                      int *n,
                                      int *independent,
                                      char *recId)
{
  if(!can_modify_records()) return;

  std::string recId_str(recId);
  std::string curveset_str(curveset_name);
  std::string curvename_str(curve_name);

  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));
  if(sina_record)
  {
    std::vector<double> y(*n);
    for(int i = 0; i < *n; i++)
    {
      y[i] = static_cast<double>(values[i]);
    }
    axom::sina::Curve curve {curvename_str, y};

    auto &curvesets = sina_record->getCurveSets();
    if(curvesets.find(curveset_str) == curvesets.end())
    {
      // Create curveset directly
      axom::sina::CurveSet cs {curveset_str};
      if(*independent != 0)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
    }
    else
    {
      // Curveset already exists
      axom::sina::CurveSet cs = curvesets.at(curveset_str);
      if(*independent != 0)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
    }
  }
}

extern "C" void sina_add_curve_int_(char *curveset_name,
                                    char *curve_name,
                                    int *values,
                                    int *n,
                                    int *independent,
                                    char *recId)
{
  if(!can_modify_records()) return;

  std::string recId_str(recId);
  std::string curveset_str(curveset_name);
  std::string curvename_str(curve_name);

  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));
  if(sina_record)
  {
    std::vector<double> y(*n);
    for(int i = 0; i < *n; i++)
    {
      y[i] = static_cast<double>(values[i]);
    }
    axom::sina::Curve curve {curvename_str, y};

    auto &curvesets = sina_record->getCurveSets();
    if(curvesets.find(curveset_str) == curvesets.end())
    {
      // Create curveset directly
      axom::sina::CurveSet cs {curveset_str};
      if(*independent != 0)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
    }
    else
    {
      // Curveset already exists
      axom::sina::CurveSet cs = curvesets.at(curveset_str);
      if(*independent != 0)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
    }
  }
}

extern "C" void sina_add_curve_long_(char *curveset_name,
                                     char *curve_name,
                                     long long int *values,
                                     int *n,
                                     int *independent,
                                     char *recId)
{
  if(!can_modify_records()) return;

  std::string recId_str(recId);
  std::string curveset_str(curveset_name);
  std::string curvename_str(curve_name);

  axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(recId_str.c_str()));
  if(sina_record)
  {
    std::vector<double> y(*n);
    for(int i = 0; i < *n; i++)
    {
      y[i] = static_cast<double>(values[i]);
    }
    axom::sina::Curve curve {curvename_str, y};

    auto &curvesets = sina_record->getCurveSets();
    if(curvesets.find(curveset_str) == curvesets.end())
    {
      // Create curveset directly
      axom::sina::CurveSet cs {curveset_str};
      if(*independent != 0)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
    }
    else
    {
      // Curveset already exists
      axom::sina::CurveSet cs = curvesets.at(curveset_str);
      if(*independent != 0)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
    }
  }
}

//=============================================================================
// Curve Ordering Functions
//=============================================================================

extern "C" void sina_set_curves_order_(int *curve_order)
{
  axom::sina::CurveSet::CurveOrder order;
  switch(*curve_order)
  {
  case 0:
    order = axom::sina::CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST;
    break;
  case 1:
    order = axom::sina::CurveSet::CurveOrder::REGISTRATION_NEWEST_FIRST;
    break;
  case 2:
    order = axom::sina::CurveSet::CurveOrder::ALPHABETIC;
    break;
  case 3:
    order = axom::sina::CurveSet::CurveOrder::REVERSE_ALPHABETIC;
    break;
  default:
    return;
  }

  axom::sina::setDefaultCurveOrder(order);
  return;
}

extern "C" void sina_set_record_curves_order_(char *recId, int *curve_order)
{
  axom::sina::Record *sina_record = Sina_Get_Record(recId);
  if(!sina_record)
  {
    return;
  }
  axom::sina::CurveSet::CurveOrder order;
  switch(*curve_order)
  {
  case 0:
    order = axom::sina::CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST;
    break;
  case 1:
    order = axom::sina::CurveSet::CurveOrder::REGISTRATION_NEWEST_FIRST;
    break;
  case 2:
    order = axom::sina::CurveSet::CurveOrder::ALPHABETIC;
    break;
  case 3:
    order = axom::sina::CurveSet::CurveOrder::REVERSE_ALPHABETIC;
    break;
  default:
    return;
  }

  sina_record->setDefaultCurveOrder(order);
  return;
}
