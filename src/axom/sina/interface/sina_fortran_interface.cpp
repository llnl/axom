// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <string.h>

#include "axom/sina/interface/sina_fortran_interface.h"

axom::sina::Document *sina_document;

extern "C" char *Get_File_Extension(char *input_fn)
{
  char *ext = strrchr(input_fn, '.');
  if(!ext)
  {
    return (new char[1] {'\0'});
  }
  return (ext + 1);
}

extern "C" void create_document_and_record_(char *recID)
{
  sina_document = new axom::sina::Document;
  // Create a record of "My Sim Code" version "1.2.3", which was run by "jdoe".
  // The run has an ID of "run1", which has to be unique to this file.
  axom::sina::ID id {recID, axom::sina::IDType::Global};
  std::unique_ptr<axom::sina::Record> myRecord {new axom::sina::Record {id, "my_type"}};
  sina_document->add(std::move(myRecord));
}

extern "C" axom::sina::Record *Sina_Get_Record()
{
  if(sina_document)
  {
    axom::sina::Document::RecordList const &allRecords = sina_document->getRecords();
    if(allRecords.size())
    {
      std::unique_ptr<axom::sina::Record> const &myRecord = allRecords.front();
      return myRecord.get();
    }
  }
  return nullptr;
}

extern "C" void sina_add_logical_(char *key, bool *value, char *units, char *tags)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    std::string key_name = std::string(key);
    if(sina_record)
    {
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
  }
}

extern "C" void sina_add_long_(char *key, long long int *value, char *units, char *tags)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    std::string key_name = std::string(key);
    if(sina_record)
    {
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
  }
}

extern "C" void sina_add_int_(char *key, int *value, char *units, char *tags)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    std::string key_name = std::string(key);
    if(sina_record)
    {
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
  }
}

extern "C" void sina_add_double_(char *key, double *value, char *units, char *tags)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    std::string key_name = std::string(key);
    if(sina_record)
    {
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
  }
}

extern "C" void sina_add_float_(char *key, float *value, char *units, char *tags)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    std::string key_name = std::string(key);
    if(sina_record)
    {
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
  }
}

extern "C" void sina_add_string_(char *key, char *value, char *units, char *tags)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    std::string key_name = std::string(key);
    std::string key_value = std::string(value);
    std::string key_units = std::string(units);
    if(sina_record)
    {
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
  }
}

extern "C" void sina_add_file_(char *filename, char *mime_type)
{
  std::string used_mime_type = "";
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
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
}

extern "C" void write_sina_document_protocol_(char *input_fn, int *protocol)
{
  std::string filename(input_fn);
  axom::sina::Protocol proto = static_cast<axom::sina::Protocol>(*protocol);
  // Save everything
  if(sina_document)
  {
    axom::sina::saveDocument(*sina_document, filename.c_str(), proto);
  }
}

extern "C" void write_sina_document_noprotocol_(char *input_fn)
{
  std::string filename(input_fn);
  // Save everything
  if(sina_document)
  {
    axom::sina::saveDocument(*sina_document, filename.c_str());
  }
}

extern "C" void sina_add_curveset_(char *name)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    if(sina_record)
    {
      axom::sina::CurveSet cs {name};
      sina_record->add(cs);
    }
  }
}

extern "C" void sina_add_curve_long_(char *curveset_name,
                                     char *curve_name,
                                     long long int *values,
                                     int *n,
                                     bool *independent)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    if(sina_record)
    {
      double y[*n];
      for(int i = 0; i < *n; i++)
      {
        y[i] = values[i];
      }
      axom::sina::Curve curve {curve_name, y, static_cast<size_t>(*n)};

      auto &curvesets = sina_record->getCurveSets();
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
}

extern "C" void sina_add_curve_int_(char *curveset_name,
                                    char *curve_name,
                                    int *values,
                                    int *n,
                                    bool *independent)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    if(sina_record)
    {
      double y[*n];
      for(int i = 0; i < *n; i++)
      {
        y[i] = values[i];
      }
      axom::sina::Curve curve {curve_name, y, static_cast<size_t>(*n)};

      auto &curvesets = sina_record->getCurveSets();
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
}

extern "C" void sina_add_curve_float_(char *curveset_name,
                                      char *curve_name,
                                      float *values,
                                      int *n,
                                      bool *independent)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    if(sina_record)
    {
      double y[*n];
      for(int i = 0; i < *n; i++)
      {
        y[i] = values[i];
      }
      axom::sina::Curve curve {curve_name, y, static_cast<size_t>(*n)};

      auto &curvesets = sina_record->getCurveSets();
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
}

extern "C" void sina_add_curve_double_(char *curveset_name,
                                       char *curve_name,
                                       double *values,
                                       int *n,
                                       bool *independent)
{
  if(sina_document)
  {
    axom::sina::Record *sina_record = Sina_Get_Record();
    if(sina_record)
    {
      axom::sina::Curve curve {curve_name, values, static_cast<size_t>(*n)};

      auto &curvesets = sina_record->getCurveSets();
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
}
