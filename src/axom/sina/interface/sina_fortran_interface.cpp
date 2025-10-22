// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <string.h>

#include "axom/sina/interface/sina_fortran_interface.h"
#include "axom/sina/core/Document.hpp"
#include "axom/sina/core/Record.hpp"
#include "axom/sina/core/CurveSet.hpp"
#include <cstring>


std::vector<std::unique_ptr<axom::sina::Record>> sinaRecordsList;
axom::sina::Document *sina_document;
char default_record_type[25] = "fortran_code_output";

// Helper that returns empty string for invalid input
inline const char* validate_optional_string_safe(char *str, int str_len) {
  if (str == NULL || str_len <= 0 || str[0] == '\0') {
    return "";  // Return empty string instead of NULL
  }
  return str;
}

// Helper function to check if modifications are allowed
inline bool can_modify_records() {
  if (sina_document != nullptr) {
    std::cerr << "ERROR: Cannot modify records after document has been created. "
              << "Call sina_write_document with preserve=0 first"
              << std::endl;
    return false;
  }
  return true;
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
  // Create a record of "My Sim Code" version "1.2.3", which was run by "jdoe".
  // The run has an ID of "run1", which has to be unique to this file.
  if (!can_modify_records()) return; 
  axom::sina::ID id {recID, axom::sina::IDType::Global};
  // std::unique_ptr<axom::sina::Record> myRecord {new axom::sina::Record {id, recType}};
  const char *validated_recType = validate_optional_string_safe(recType, recType_length);
  if (validated_recType[0] == '\0') {
    validated_recType = default_record_type;
  }
  sinaRecordsList.emplace_back(std::make_unique<axom::sina::Record>(id, validated_recType));
}

extern "C" axom::sina::Record *Sina_Get_Record(char * recId=NULL)
{

    if (recId == NULL || recId[0] == '\0') {
      std::unique_ptr<axom::sina::Record> const &myRecord = sinaRecordsList.front();
      return myRecord.get();
    }
    else {
      axom::sina::ID id {recId, axom::sina::IDType::Global};
      for (const std::unique_ptr<axom::sina::Record> &myRecord : sinaRecordsList) {
        const char* current_id_str = myRecord->getId().getId().c_str();
        // Compare the input C-string (recId) with the current record's C-string ID
        if (strcmp(recId, current_id_str) == 0) { 
          return myRecord.get();
        } 
      }
      // Didn't match we will return a new record
      sina_create_record_(recId, "", strlen(recId), 0);
      return sinaRecordsList.back().get();
    }
  return nullptr;
}

extern "C" void sina_add_logical_(char *key, bool *value, char *units, char *tags, char *recId, int key_len, int units_len, int tags_len, int recId_len)
{
    if (!can_modify_records()) return; 
    const char *validated_recId = validate_optional_string_safe(recId, recId_len);
    axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
    std::string key_name = std::string(key);
      axom::sina::Datum datum {static_cast<double>(*value)};
      if(units)
      {
        std::string key_units = std::string(units);
        if(key_units != "")
        {
          datum.setUnits(key_units);
        }
      }
      if(tags)
      {
        std::vector<std::string> tagVector = {tags};
        if(tagVector.front() != "")
        {
          datum.setTags(tagVector);
        }
      }
      sina_record->add(key_name, datum);
}

extern "C" void sina_add_long_(char *key, long long int *value, char *units, char *tags, char *recId, int key_len, int units_len, int tags_len, int recId_len)
{
    if (!can_modify_records()) return; 
    const char *validated_recId = validate_optional_string_safe(recId, recId_len);
    axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
    std::string key_name = std::string(key);
      axom::sina::Datum datum {static_cast<double>(*value)};
      if(units)
      {
        std::string key_units = std::string(units);
        if(key_units != "")
        {
          datum.setUnits(key_units);
        }
      }
      if(tags)
      {
        std::vector<std::string> tagVector = {tags};
        if(tagVector.front() != "")
        {
          datum.setTags(tagVector);
        }
      }
      sina_record->add(key_name, datum);
}

extern "C" void sina_add_int_(char *key, int *value, char *units, char *tags, char *recId, int key_len, int units_len, int tags_len, int recId_len)
{
    if (!can_modify_records()) return; 
    const char *validated_recId = validate_optional_string_safe(recId, recId_len);
    axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
    std::string key_name = std::string(key);
      axom::sina::Datum datum {static_cast<double>(*value)};
      if(units)
      {
        std::string key_units = std::string(units);
        if(key_units != "")
        {
          datum.setUnits(key_units);
        }
      }
      if(tags)
      {
        std::vector<std::string> tagVector = {tags};
        if(tagVector.front() != "")
        {
          datum.setTags(tagVector);
        }
      }
      sina_record->add(key_name, datum);
}

extern "C" void sina_add_double_(char *key, double *value, char *units, char *tags, char *recId, int key_len, int units_len, int tags_len, int recId_len)
{
    if (!can_modify_records()) return; 
    const char *validated_recId = validate_optional_string_safe(recId, recId_len);
    axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
    std::string key_name = std::string(key);
      axom::sina::Datum datum {*value};
      if(units)
      {
        std::string key_units = std::string(units);
        if(key_units != "")
        {
          datum.setUnits(key_units);
        }
      }
      if(tags)
      {
        std::vector<std::string> tagVector = {tags};
        if(tagVector.front() != "")
        {
          datum.setTags(tagVector);
        }
      }
      sina_record->add(key_name, datum);
}

extern "C" void sina_add_float_(char *key, float *value, char *units, char *tags, char *recId, int key_len, int units_len, int tags_len, int recId_len)
{
    if (!can_modify_records()) return; 
    const char *validated_recId = validate_optional_string_safe(recId, recId_len);
    axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
    std::string key_name = std::string(key);
      axom::sina::Datum datum {*value};
      if(units)
      {
        std::string key_units = std::string(units);
        if(key_units != "")
        {
          datum.setUnits(key_units);
        }
      }
      if(tags)
      {
        std::vector<std::string> tagVector = {tags};
        if(tagVector.front() != "")
        {
          datum.setTags(tagVector);
        }
      }
      sina_record->add(key_name, datum);
}

extern "C" void sina_add_string_(char *key, char *value, char *units, char *tags, char *recId, int key_len, int units_len, int tags_len, int recId_len)
{
    if (!can_modify_records()) return; 
    const char *validated_recId = validate_optional_string_safe(recId, recId_len);
    axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
    std::string key_name = std::string(key);
    std::string key_value = std::string(value);
    std::string key_units = std::string(units);
      axom::sina::Datum datum {key_value};
      if(units)
      {
        std::string key_units = std::string(units);
        if(key_units != "")
        {
          datum.setUnits(key_units);
        }
      }
      if(tags)
      {
        std::vector<std::string> tagVector = {tags};
        if(tagVector.front() != "")
        {
          datum.setTags(tagVector);
        }
      }
      sina_record->add(key_name, datum);
}

extern "C" void sina_add_file_(char *filename, char *mime_type, char *recId, int file_len, int mime_len, int recId_len)
{
    if (!can_modify_records()) return; 
    std::string used_mime_type = "";
    const char *validated_recId = validate_optional_string_safe(recId, recId_len);
    axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
    if(mime_type)
    {
      used_mime_type = std::string(mime_type);
    }
    axom::sina::File my_file {filename};
    if(used_mime_type != "")
    {
      my_file.setMimeType(used_mime_type);
    }
    else
    {
      used_mime_type = Get_File_Extension(filename);
      my_file.setMimeType(used_mime_type);
    }

    if(sina_record)
    {
      sina_record->add(my_file);
    }
}

extern "C" void sina_write_document_all_args_(char *input_fn, int *protocol, int *preserve, int *mergeProtocol=0)
{

  // Lets create the document if needed
  if (sina_document == nullptr) {
    sina_document = new axom::sina::Document();
    
    // Move all records into the document
    for (auto& uniquePtr : sinaRecordsList) {
      if (uniquePtr) {
        sina_document->add(std::move(uniquePtr));
      }
    }
    sinaRecordsList.clear();  // Clear empty pointers
  }


  std::string filename(input_fn);
  axom::sina::Protocol proto = static_cast<axom::sina::Protocol>(*protocol);
  // Save everything
  axom::sina::appendDocument(*sina_document, filename.c_str(), *mergeProtocol, proto);

  // Do we want to bring it back?
  if (*preserve == 0) {
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

extern "C" void sina_write_document_protocol_preserve_nomerge_(char *input_fn, int *protocol, int *preserve)
{
    int default_merge_protocol = 0;

    sina_write_document_all_args_(input_fn, protocol, preserve, &default_merge_protocol);

}

extern "C" void sina_add_curveset_(char *name, char *recId, int name_len, int recId_len)
{
    if (!can_modify_records()) return; 
    const char *validated_recId = validate_optional_string_safe(recId, recId_len);

    axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
    if(sina_record)
    {
      axom::sina::CurveSet cs {name};
      sina_record->add(cs);
    }
}

extern "C" void sina_add_curve_long_(char *curveset_name,
                                     char *curve_name,
                                     long long int *values,
                                     int *n,
                                     bool *independent,
                                     char *recId,
                                     int curveset_len,
                                     int curvename_len,
                                     int recId_len)
{
    if (!can_modify_records()) return; 
    const char *validated_recId = validate_optional_string_safe(recId, recId_len);
    axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
    if(sina_record)
    {
      std::vector<double> y(*n);
      for(int i = 0; i < *n; i++)
      {
        y[i] = values[i];  // cast to doubles
      }
      axom::sina::Curve curve {curve_name, y};

      auto &curvesets = sina_record->getCurveSets();
      if (curvesets.find(curveset_name) == curvesets.end()) {
        // missing curveset let's create it (useful for append mode where rec was deleted)
        sina_add_curveset_(curveset_name, recId, strlen(curveset_name), strlen(recId));
        auto &curvesets = sina_record->getCurveSets();
      }

      axom::sina::CurveSet cs = curvesets.at(curveset_name);
      if(*independent)
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

extern "C" void sina_add_curve_int_(char *curveset_name,
                                    char *curve_name,
                                    int *values,
                                    int *n,
                                    bool *independent,
                                    char *recId,
                                    int curveset_len,
                                    int curvename_len,
                                    int recId_len)
{
    if (!can_modify_records()) return; 
      const char *validated_recId = validate_optional_string_safe(recId, recId_len);
      axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
      std::vector<double> y(*n);
      for(int i = 0; i < *n; i++)
      {
        y[i] = values[i];
      }
      axom::sina::Curve curve {curve_name, y};

      auto &curvesets = sina_record->getCurveSets();
      if (curvesets.find(curveset_name) == curvesets.end()) {
        // missing curveset let's create it (useful for append mode where rec was deleted)
        sina_add_curveset_(curveset_name, recId, strlen(curveset_name), strlen(recId));
        auto &curvesets = sina_record->getCurveSets();
      }

      axom::sina::CurveSet cs = curvesets.at(curveset_name);
      if(*independent)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
}

extern "C" void sina_add_curve_float_(char *curveset_name,
                                      char *curve_name,
                                      float *values,
                                      int *n,
                                      bool *independent,
                                      char *recId,
                                     int curveset_len,
                                     int curvename_len,
                                     int recId_len)
{
      if (!can_modify_records()) return; 
      const char *validated_recId = validate_optional_string_safe(recId, recId_len);
      axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
      std::vector<double> y(*n);
      for(int i = 0; i < *n; i++)
      {
        y[i] = values[i];  // cast to doubles
      }
      axom::sina::Curve curve {curve_name, y};

      auto &curvesets = sina_record->getCurveSets();
      if (curvesets.find(curveset_name) == curvesets.end()) {
        // missing curveset let's create it (useful for append mode where rec was deleted)
        sina_add_curveset_(curveset_name, recId, strlen(curveset_name), strlen(recId));
        auto &curvesets = sina_record->getCurveSets();
      }

      axom::sina::CurveSet cs = curvesets.at(curveset_name);
      if(*independent)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
}

extern "C" void sina_add_curve_double_(char *curveset_name,
                                       char *curve_name,
                                       double *values,
                                       int *n,
                                       bool *independent,
                                       char *recId,
                                       int curveset_len,
                                       int curvename_len,
                                       int recId_len)
{
      if (!can_modify_records()) return; 
      const char *validated_recId = validate_optional_string_safe(recId, recId_len);

      axom::sina::Record *sina_record = Sina_Get_Record(const_cast<char *>(validated_recId));
      axom::sina::Curve curve {curve_name, values, static_cast<size_t>(*n)};

      auto &curvesets = sina_record->getCurveSets();
      if (curvesets.find(curveset_name) == curvesets.end()) {
        // missing curveset let's create it (useful for append mode where rec was deleted)
        sina_add_curveset_(curveset_name, recId, strlen(curveset_name), strlen(recId));
        auto &curvesets = sina_record->getCurveSets();
      }

      axom::sina::CurveSet cs = curvesets.at(curveset_name);
      if(*independent)
      {
        cs.addIndependentCurve(curve);
      }
      else
      {
        cs.addDependentCurve(curve);
      }
      sina_record->add(cs);
}

//=============================================================================
// Curve Ordering Functions 
//=============================================================================


extern "C" void sina_set_curves_order_(int* curve_order)
{
        axom::sina::CurveSet::CurveOrder order;
        switch (*curve_order) {
            case 0: order = axom::sina::CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST; break;
            case 1: order = axom::sina::CurveSet::CurveOrder::REGISTRATION_NEWEST_FIRST; break;
            case 2: order = axom::sina::CurveSet::CurveOrder::ALPHABETIC; break;
            case 3: order = axom::sina::CurveSet::CurveOrder::REVERSE_ALPHABETIC; break;
            default: return;
        }
        
        axom::sina::setDefaultCurveOrder(order);
        return;
}

extern "C" void sina_set_record_curves_order_(char* recId, int* curve_order)
{
    axom::sina::Record *sina_record = Sina_Get_Record(recId);
        if (!sina_record) {
          return;
        }
        axom::sina::CurveSet::CurveOrder order;
        switch (*curve_order) {
            case 0: order = axom::sina::CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST; break;
            case 1: order = axom::sina::CurveSet::CurveOrder::REGISTRATION_NEWEST_FIRST; break;
            case 2: order = axom::sina::CurveSet::CurveOrder::ALPHABETIC; break;
            case 3: order = axom::sina::CurveSet::CurveOrder::REVERSE_ALPHABETIC; break;
            default: return;
        }
        
        sina_record->setDefaultCurveOrder(order);
        return;
}
