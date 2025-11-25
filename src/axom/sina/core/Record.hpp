// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_RECORD_HPP
#define SINA_RECORD_HPP

/*!
 ******************************************************************************
 *
 * \file Record.hpp
 *
 * \brief   Header file for Sina Record class
 *
 * \sa DataHolder.hpp
 *
 ******************************************************************************
 */

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "conduit.hpp"

#include "axom/sina/core/ID.hpp"
#include "axom/sina/core/DataHolder.hpp"
#include "axom/sina/core/CurveSet.hpp"
#include "axom/sina/core/Datum.hpp"
#include "axom/sina/core/File.hpp"

namespace axom
{
namespace sina
{

/**
 * FileEqualByURI is used to store files in a Record.
 */
struct FileEqualByURI
{
  bool operator()(const File &file1, const File &file2) const
  {
    return file1.getUri() == file2.getUri();
  }
};

/**
 * FileHashByURI is used to store files in a Record. Files are stored according
 * to the hash of their URI.
 */
struct FileHashByURI
{
  size_t operator()(const File &file) const { return std::hash<std::string>()(file.getUri()); }
};

/**
 * \brief An object representing an entry in a Document's Record list
 *
 * The Record class represents an entry in a Document's Record list. Records represent the data to be stored
 * (as opposed to the relationships between data)--natural scopes for Records include things like a single run
 * of an application, an msub job, a cluster of runs that has some metadata attached to the cluster (this Record
 * might have a "contains" Relationship for all the runs within it), etc.
 *
 * Each Record must have a type and an id. Each Record can also have a list of
 * File objects and a map of Datum objects.
 *
 * \code
 *   axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
 *   std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};
 *   std::vector<std::string> myTags{"input"};
 *   axom::sina::Datum myDatum{12, myTags};
 *   myRecord->add("my_scalar",std::move(myDatum));
 *   std::cout << myRecord->toNode().to_json() << std::endl;
 * \endcode
 *
 * The output would be:
 * \code{.json}
 *   {"local_id":"my_record","type":"my_type","data":{"my_scalar":{"tags":["input"],"value":12.0}}}
 * \endcode
 */
class Record : public DataHolder
{
public:
  /**
     * An unordered set of File objects.
     */
  using FileSet = std::unordered_set<File, FileHashByURI, FileEqualByURI>;

  /**
    * An enum member representing the default ordering of Curves within CurveSets.
    */
  static CurveSet::CurveOrder defaultCurveOrder;

  /**    
   * \brief Set the default curve order for all Records.
   * 
   * \param order the new default curve order
  */
  void setDefaultCurveOrder(CurveSet::CurveOrder order);

  /**
   * \brief Get the effective default curve order for this Record.
   * 
   * \return the instance default if set, otherwise the global default
   */
  CurveSet::CurveOrder getDefaultCurveOrder() const;

  /**
     * \brief Construct a new Record.
     *
     * \param id the ID of the record
     * \param type the type of the record
     */
  Record(ID id, std::string type);

  /**
     * \brief Construct a Record from its conduit Node representation.
     *
     * For the purpose of the order in which curves are written, Nodes are
     * assumed to be ordered OLDEST FIRST. This should be consistent within codes.
     *
     * \param asNode the Record as a Node
     */
  explicit Record(conduit::Node const &asNode);

  /**
     * Disable the copy constructor.
     */
  Record(Record const &) = delete;

  /**
     * Disable copy assignment.
     */
  Record &operator=(Record const &) = delete;

  /**
     * \brief Get the Record's ID.
     *
     * \return the ID
     */
  ID const &getId() const noexcept { return id.getID(); }

  /**
     * \brief Get the Record's type.
     *
     * \return the Record's type
     */
  std::string const &getType() const noexcept { return type; }

  /**
     * \brief Remove a File from this record.
     *
     * \param file the File to remove
     */
  void remove(File const &file);

  using DataHolder::add;
  /**
     * \brief Add a File to this record.
     *
     * \param file the File to add
     */
  void add(File file);

  /**
     * \brief Get the files associated with this record.
     *
     * \return the record's files
     */
  FileSet const &getFiles() const noexcept { return files; }

  /**
     * \brief Convert this record to its conduit Node representation.
     *
     * \param curveOrder Optionally, specify an order that CurveSets should be written in. Options
     *                   are enumerated in CurveSets; REGISTRATION_OLDEST_FIRST is "ULTRA style".
     *
     * \return the Node representation of this record.
     */
  conduit::Node toNode(CurveSet::CurveOrder curveOrder) const override;

  /**
   * \brief Convert using this record's default curve order.
   */
  conduit::Node toNode() const;

  /**
    * \brief Add another record to this one as library data.
    *
    * Useful for libraries that can run in standalone mode; the host
    * simply calls this method on the record the library produces.
    * Merges file lists.
    *
    * \param name The host code's name for the library
    */
  void addRecordAsLibraryData(Record const &childRecord, std::string const &name);

private:
  internal::IDField id;
  std::string type;
  FileSet files;
  bool hasInstanceDefault = false;
  CurveSet::CurveOrder instanceDefaultCurveOrder = sinaDefaultCurveOrder;  // Per-record default
};

/**
 * \brief An object to convert Conduit Nodes into Records
 *
 * A RecordLoader is used to convert conduit::Node instances which represent
 * Sina Records into instances of their corresponding axom::sina::Record
 * subclasses. For convenience, a RecordLoader capable of handling Records of all known
 * types can be created using createRecordLoaderWithAllKnownTypes:
 *
 * \code
 *   axom::sina::Document myDocument = axom::sina::Document(jObj, axom::sina::createRecordLoaderWithAllKnownTypes());
 * \endcode
 */
class RecordLoader
{
public:
  /**
     * A TypeLoader is a function which converts records of a specific type
     * to their corresponding sub classes.
     */
  using TypeLoader = std::function<std::unique_ptr<Record>(conduit::Node const &)>;

  /**
     * \brief Add a function for loading records of the specified type.
     *
     * \param type the type of records this function can load
     * \param loader the function which can load the records
     */
  void addTypeLoader(std::string const &type, TypeLoader loader);

  /**
     * \brief Load a Record from its conduit Node representation.
     *
     * \param recordAsNode the Record as a Node
     * \return the Record
     */
  std::unique_ptr<Record> load(conduit::Node const &recordAsNode) const;

  /**
     * \brief Check whether this loader can load records of the given type.
     *
     * \param type the type of the records to check
     * \return whether records of the given type can be loaded
     */
  bool canLoad(std::string const &type) const;

private:
  std::unordered_map<std::string, TypeLoader> typeLoaders;
};

/**
 * \brief Create a RecordLoader which can load records of all known types.
 *
 * \return the newly-created loader
 */
RecordLoader createRecordLoaderWithAllKnownTypes();

}  // namespace sina
}  // namespace axom

#endif  //SINA_RECORD_HPP
