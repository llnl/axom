// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file JSONSchemaWriter.hpp
 *
 * \brief This file contains the class definition of the JSONSchemaWriter.
 *******************************************************************************
 */

#ifndef INLET_JSONSCHEMAWRITER_HPP
#define INLET_JSONSCHEMAWRITER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>

#include "axom/sidre.hpp"
#include "axom/inlet/Writer.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class JSONSchemaWriter
 *
 * \brief A Writer that is able to write the input file schema in JSON Schema
 * format for a given input file.
 *
 * \see Writer
 * \see https://json-schema.org/
 *******************************************************************************
 */
class JSONSchemaWriter : public Writer
{
public:
  /*!
  *******************************************************************************
  * \brief A constructor for JSONSchemaWriter.
  * 
  * \param [in] fileName The name of the file the schema should be written to.
  *******************************************************************************
  */
  JSONSchemaWriter(const std::string& fileName);

  void documentContainer(const Container& container) override;

  void finalize() override;

  virtual ~JSONSchemaWriter() = default;

private:
  conduit::Node m_schemaRoot;
  std::string m_fileName;
  std::vector<std::string> m_ArrayPaths;
  std::vector<std::string> m_DictionaryPaths;
};

}  // namespace inlet
}  // namespace axom

#endif
