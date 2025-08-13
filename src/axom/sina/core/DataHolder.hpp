// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_DATAHOLDER_HPP
#define SINA_DATAHOLDER_HPP

/*!
 ******************************************************************************
 *
 * \file DataHolder.hpp
 *
 * \brief   Header file for Sina DataHolder class
 *
 ******************************************************************************
 */

#include <string>
#include <memory>
#include <unordered_map>

#include "conduit.hpp"

#include "axom/sina/core/Datum.hpp"
#include "axom/sina/core/CurveSet.hpp"

namespace axom
{
namespace sina
{

/**
 * \brief A DataHolder is a basic container for certain types of information.
 *
 * DataHolders contain curves, libraries, and data (Datum), and represent
 * all the information a library can have associated with it. Records expand
 * on DataHolders to contain additional info.
 *
 * \sa Record
 * \sa LibraryDataMap
 */
class DataHolder
{
public:
  /**
     * An unordered map of Datum objects.
     */
  using DatumMap = std::unordered_map<std::string, Datum>;

  /**
     * An unordered map of CurveSet objects.
     */
  using CurveSetMap = std::unordered_map<std::string, CurveSet>;

  /**
     * An unordered map of shared pointers to DataHolder objects.
     */
  using LibraryDataMap = std::unordered_map<std::string, std::shared_ptr<DataHolder>>;

  /**
     * Construct an empty DataHolder.
     */
  DataHolder() = default;

  /**
     * Virtual destructor to automatically clean up resources held by an instance of the DataHolder class.
     */
  virtual ~DataHolder() = default;

  /**
     * Copy constructor that disallows this constructor type.
     */
  DataHolder(DataHolder const &) = delete;

  /**
     * Disable copy assignment.
     */
  DataHolder &operator=(DataHolder const &) = delete;

  /**
     * \brief Construct a DataHolder from its conduit Node representation.
     *
     * \param asNode the DataHolder as a Node
     */
  explicit DataHolder(conduit::Node const &asNode);

  /**
     * \brief Get the DataHolder's data.
     *
     * \return the DataHolder's data
     */
  DatumMap const &getData() const noexcept { return data; }

  /**
     * \brief Add a Datum to this DataHolder.
     *
     * \param name the key for the Datum to add
     * \param datum the Datum to add
     */
  void add(std::string name, Datum datum);

  /**
     * \brief Add a CurveSet to this DataHolder.
     *
     * \param curveSet the CurveSet to add
     */
  void add(CurveSet curveSet);

  /**
     * \brief Get the curve sets associated with this DataHolder.
     *
     * \return the dataholder's curve sets
     */
  CurveSetMap const &getCurveSets() const noexcept { return curveSets; }

  /**
     * \brief Add a new library to this DataHolder.
     *
     * If you try to add a library with a name that already exists, the old
     * library will be replaced.
     *
     * \return a pointer to a new DataHolder for a library
     *         of the given name.
     */
  std::shared_ptr<DataHolder> addLibraryData(std::string const &name);

  /**
     * \brief Add a new library to this DataHolder with existing library data.
     *
     * \return a pointer to a new DataHolder for a library of the given name.
     */
  std::shared_ptr<DataHolder> addLibraryData(std::string const &name,
                                             conduit::Node existingLibraryData);

  /**
     * \brief Get all library data associated with this DataHolder.
     *
     * \return the dataholder's library data
     */
  LibraryDataMap const &getLibraryData() const noexcept { return libraryData; }

  /**
     * \brief Get a specific library associated with this DataHolder.
     *
     * \return the dataholder's library data
     */
  std::shared_ptr<DataHolder> getLibraryData(std::string const &libraryName)
  {
    return libraryData.at(libraryName);
  }

  /**
     * \brief Get a specific library associated with this DataHolder.
     *
     * \return the dataholder's library data
     */
  std::shared_ptr<DataHolder> const getLibraryData(std::string const &libraryName) const
  {
    return libraryData.at(libraryName);
  }

  /**
     * \brief Get the user-defined content of the object.
     *
     * \return the user-defined content
     */
  conduit::Node const &getUserDefinedContent() const noexcept { return userDefined; }

  /**
     * \brief Get the user-defined content of the object.
     *
     * \return the user-defined content
     */
  conduit::Node &getUserDefinedContent() noexcept { return userDefined; }

  /**
     * \brief Set the user-defined content of the object.
     *
     * \param userDefined the user-defined content. Must be an object (key/value pairs)
     */
  void setUserDefinedContent(conduit::Node userDefined);

  /**
     * \brief Convert this DataHolder to its conduit Node representation.
     *
     * \param curveOrder the order in which curves within curvesets should be written.
     *                   A default is set in Record.
     *
     * \return the Node representation of this DataHolder.
     */
  virtual conduit::Node toNode(CurveSet::CurveOrder curveOrder) const;

private:
  CurveSetMap curveSets;
  DatumMap data;
  LibraryDataMap libraryData;
  conduit::Node userDefined;
};

}  // namespace sina
}  // namespace axom

#endif  //SINA_DATAHOLDER_HPP
