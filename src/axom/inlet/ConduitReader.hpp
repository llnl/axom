// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file ConduitReader.hpp
 *
 * \brief This file contains the class definition of the ConduitReader.
 *******************************************************************************
 */

#ifndef INLET_CONDUITREADER_HPP
#define INLET_CONDUITREADER_HPP

#include "axom/inlet/Reader.hpp"

#include "conduit.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class ConduitReader
 *
 * \brief A Reader that is able to read variables from a file compatible with
 * one of Conduit's protocols.
 *
 * \see Reader
 * 
 * \note This class should not be used directly, instead, use the YAMLReader
 * and JSONReader derived classes
 *******************************************************************************
 */
class ConduitReader : public Reader
{
public:
  /*!
   *****************************************************************************
   * \brief Constructs a reader with the selected Conduit protocol
   * \param [in] protocol The Conduit protocol, currently only "yaml" and "json"
   * are supported
   *****************************************************************************
   */
  ConduitReader(const std::string& protocol);

  bool parseFile(const std::string& filePath) override;

  bool parseString(const std::string& stringToRead) override;

  ReaderResult getBool(const std::string& id, bool& value) override;

  ReaderResult getDouble(const std::string& id, double& value) override;

  ReaderResult getInt(const std::string& id, int& value) override;

  ReaderResult getString(const std::string& id, std::string& value) override;

  ReaderResult getIntMap(const std::string& id, std::unordered_map<int, int>& values) override;
  ReaderResult getIntMap(const std::string& id, std::unordered_map<VariantKey, int>& values) override;

  ReaderResult getDoubleMap(const std::string& id, std::unordered_map<int, double>& values) override;
  ReaderResult getDoubleMap(const std::string& id,
                            std::unordered_map<VariantKey, double>& values) override;

  ReaderResult getBoolMap(const std::string& id, std::unordered_map<int, bool>& values) override;
  ReaderResult getBoolMap(const std::string& id,
                          std::unordered_map<VariantKey, bool>& values) override;

  ReaderResult getStringMap(const std::string& id,
                            std::unordered_map<int, std::string>& values) override;
  ReaderResult getStringMap(const std::string& id,
                            std::unordered_map<VariantKey, std::string>& values) override;

  ReaderResult getIndices(const std::string& id, std::vector<int>& indices) override;
  ReaderResult getIndices(const std::string& id, std::vector<VariantKey>& indices) override;

  /*!
   *****************************************************************************
   * \brief Get a function from the input file
   *
   * \param [in]  id    The identifier to the function that will be retrieved
   * \param [in]  ret_type    The return type of the function
   * \param [in]  arg_types    The argument types of the function
   *
   * \return The function, compares false if not found
   * \note Conduit does not support functions - calling this function will result
   * in a SLIC_ERROR
   *****************************************************************************
   */
  FunctionVariant getFunction(const std::string& id,
                              const FunctionTag ret_type,
                              const std::vector<FunctionTag>& arg_types) override;

  std::vector<std::string> getAllNames() override;

  /*!
   *****************************************************************************
   * \brief The base index for Conduit arrays
   *****************************************************************************
   */
  static const int baseIndex = 0;

private:
  ReaderResult getValue(const conduit::Node* node, int& value);
  ReaderResult getValue(const conduit::Node* node, std::string& value);
  ReaderResult getValue(const conduit::Node* node, double& value);
  ReaderResult getValue(const conduit::Node* node, bool& value);

  template <typename T>
  ReaderResult getDictionary(const std::string& id, std::unordered_map<VariantKey, T>& values);

  template <typename T>
  ReaderResult getArray(const std::string& id, std::unordered_map<int, T>& values);
  conduit::Node m_root;
  const std::string m_protocol;
};

}  // end namespace inlet
}  // end namespace axom

#endif
