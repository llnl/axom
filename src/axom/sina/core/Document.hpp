// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_DOCUMENT_HPP
#define SINA_DOCUMENT_HPP

/*!
 ******************************************************************************
 *
 * \file Document.hpp
 *
 * \brief   Header file for Sina Document class
 *
 ******************************************************************************
 */

#include "axom/config.hpp"
#include "axom/sina/core/Record.hpp"
#include "axom/sina/core/Relationship.hpp"

#include "conduit.hpp"
#include "conduit_relay.hpp"
#include "conduit_relay_io.hpp"

#include <memory>
#include <vector>

#define SINA_FILE_FORMAT_VERSION_MAJOR 1
#define SINA_FILE_FORMAT_VERSION_MINOR 0

namespace axom
{
namespace sina
{

/**
 * @brief File format enum for explicit format specification
 * 
 * This enum provides type-safe format specification for saveDocument and
 * appendDocument functions. It works alongside the existing supported_types
 * for backward compatibility.
 * 
 * The enum values map to supported_types indices:
 * - AUTO_DETECT: Special value for auto-detection from file extension
 * - JSON: maps to supported_types[0] = "JSON"
 * - HDF5: maps to supported_types[1] = "HDF5" (if AXOM_USE_HDF5 is defined)
 */

enum class Protocol
{
  AUTO_DETECT = -1,  ///< Automatically detect from file extension (.json, .h5, .hdf5)
  JSON = 0,          ///< Force JSON format (supported_types[0])
  HDF5 = 1           ///< Force HDF5 format (supported_types[1], requires AXOM_USE_HDF5)
};

const std::vector<std::string> supported_types = {"JSON",
#ifdef AXOM_USE_HDF5
                                                  "HDF5"
#endif
};

/**
 * \brief The string used to replace '/' in parent node names when saving to HDF5.
 */
const std::string slashSubstitute = "__SINA_SLASHREPLACE__";

/**
 * \brief An object representing the top-level object of a Sina file
 *
 * A Document represents the top-level object of a file conforming to the
 * Sina schema. When serialized, these documents can be ingested into a
 * Sina database and used with the Sina tool. Sina files are defaulted to
 * JSON but optionally support HDF5.
 *
 * Documents contain at most two objects: a list of Records and a list of Relationships. A simple, empty document:
 * \code{.json}
 *   {
 *      "records": [],
 *      "relationships": []
 *   }
 * \endcode
 *
 * The "records" list can contain Record objects and their inheriting types, such as Run (for a full list, please see
 * the inheritance diagram in the Record documentation). The "relationships" list can contain Relationship objects.
 *
 * Documents can be assembled programatically and/or generated from existing JSON. An example of an assembled
 * Document is provided on the main page. To load a Document from an existing JSON file:
 * \code
 *   axom::sina::Document myDocument = axom::sina::loadDocument("path/to/infile.json");
 * \endcode
 *
 * To generate a Document from a JSON string and vice versa:
 * \code
 *   std::string my_json = "{\"records\":[{\"type\":\"run\",\"id\":\"test\"}],\"relationships\":[]}";
 *   axom::sina::Document myDocument = axom::sina::Document(my_json, axom::sina::createRecordLoaderWithAllKnownTypes());
 *   std::cout << myDocument.toJson() << std::endl;);
 * \endcode
 *
 * You can add further entries to the Document using add():
 * \code
 *   std::unique_ptr<sina::Record> myRun{new axom::sina::Run{someID, "My Sim Code", "1.2.3", "jdoe"}};
 *   axom::sina::Relationship myRelationship{someID, "comes before", someOtherID};
 *   myDocument.add(myRun);
 *   myDocument.add(myRelationship);
 * \endcode
 *
 * You can also export your Document to file:
 * \code
 *   axom::sina::saveDocument(myDocument, "path/to/outfile.json")
 * \endcode
 * 
 *  Loading and Saving documents will default to the JSON file type, but if an optional file type is
 *  loaded the Protocol parameter will control your file type. For example with HDF5:
 * \code
 *   axom::sina::Document myDocument = axom::sina::loadDocument("path/to/infile.hdf5, Protocol::HDF5");
 *   axom::sina::saveDocument(myDocument, "path/to/outfile.hdf5", Protocol::HDF5)
 * \endcode
 *
 * Check the Sina file format version with:
 * \code
 *   axom::sina::getSinaFileFormatVersion();
 * \endcode
 *
 */
class Document
{
public:
  /**
     * A vector of pointers to Record objects.
     */
  using RecordList = std::vector<std::unique_ptr<Record>>;

  /**
     * A vector of Relationship objects.
     */
  using RelationshipList = std::vector<Relationship>;

  /**
     * Construct an empty Document.
     */
  Document() = default;

  /**
     * Disable copying Document objects. We must do this since we hold
     * pointers to polymorphic objects.
     */
  Document(Document const &) = delete;

  /**
     * Disabling copy assignment.
     */
  Document &operator=(Document const &) = delete;

  /**
     * Move constructor which should be handled by the compiler.
     */
  Document(Document &&) = default;

  /**
     * Move assignment which should be handled by the compiler.
     */
  Document &operator=(Document &&) = default;

  /**
     * \brief Create a Document from its Conduit Node representation
     *
     * \param asNode the Document as a Node
     * \param recordLoader an RecordLoader to use to load the different
     *                     types of records which may be in the document
     */
  Document(conduit::Node const &asNode, RecordLoader const &recordLoader);

  /**
     * \brief Create a Document from a JSON string representation
     *
     * \param asJson the Document as a JSON string
     * \param recordLoader an RecordLoader to use to load the different
     *                     types of records which may be in the document
     */
  Document(std::string const &asJson, RecordLoader const &recordLoader);

  /**
     * \brief Add the given record to this document.
     *
     * \param record the record to add
     */
  void add(std::unique_ptr<Record> record);

  /**
     * \brief Get the list of records currently in this document.
     *
     * \return the list of records
     */
  RecordList const &getRecords() const noexcept { return records; }

  /**
     * \brief Add a relationship to this document
     *
     * \param relationship the relationship to add
     */
  void add(Relationship relationship);

  /**
     * \brief Get the list of relationships in this document.
     *
     * \return the list of relationships
     */
  RelationshipList const &getRelationships() const noexcept { return relationships; }

  /**
     * \brief Convert this document to a conduit Node.
     *
     * \return the contents of the document as a Node
     */
  conduit::Node toNode() const;

#ifdef AXOM_USE_HDF5

  /**
   *  \brief Convert this document to the HDF5 format (escaped slashes)
   * 
   *  \param writeTo an empty conduit node to write the new node into
   * 
   *  \return This node with slashes escaped for HDF5.
   */
  conduit::Node &toHDF5Node(conduit::Node &writeTo) const;

  /**
   *  \brief Dump this document as an HDF5 File
   * 
   *  \param filename the location of which to save the file
   */
  void toHDF5(const std::string &filename) const;
#endif

  /**
     * \brief Convert this document to a JSON string.
     *
     * \return the contents of the document as a JSON string
     */
  std::string toJson(conduit::index_t indent = 0,
                     conduit::index_t depth = 0,
                     const std::string &pad = "",
                     const std::string &eoe = "") const;

  /**
    * \brief Get the list of file types currently supported by the implementation.
    * 
    * \return a string of supported file types
    */
  std::string get_supported_file_types();

private:
  /**
     * Constructor helper method, extracts info from a conduit Node.
     */
  void createFromNode(conduit::Node const &asNode, RecordLoader const &recordLoader);
  RecordList records;
  RelationshipList relationships;
};

/**
 * \brief Save the given Document to the specified location. If the given file exists,
 *        it will be overwritten.
 * 
 * \param document the Document to save
 * \param fileName the location of which to save the file
 * \param protocol the file type requested to save as contained in supported_types, default = JSON
 * \throws std::ios::failure if there are any IO errors
 *         std::invalid_argument if the protocol given is an undefined, optional protocol
 */
void saveDocument(Document const &document,
                  std::string const &fileName,
                  Protocol protocol = Protocol::AUTO_DETECT);

/**
 * \brief Get the current file format version.
 *
 * \return A string representing the file format version.
 */
inline std::string getSinaFileFormatVersion()
{
  return std::to_string(SINA_FILE_FORMAT_VERSION_MAJOR) + "." +
    std::to_string(SINA_FILE_FORMAT_VERSION_MINOR);
}

void restoreSlashes(const conduit::Node &modifiedNode, conduit::Node &restoredNode);

/**
 * \brief Load a document from the given path. Only records which this library
 *        knows about will be able to be loaded.
 *
 * \param path the file system path from which to load the document
 * \param protocol the type of file being loaded, default = JSON
 * \return the loaded Document
 */
Document loadDocument(std::string const &path, Protocol protocol = Protocol::JSON);

/**
 * \brief Load a document from the given path.
 *
 * \param path the file system path from which to load the document
 * \param recordLoader the RecordLoader to use to load the different types
 *                     of records
 * \param protocol the type of file being loaded, default = JSON
 * \throws std::invalid_argument if the protocol given is an undefined, optional protocol
 * \return the loaded Document
 */
Document loadDocument(std::string const &path,
                      RecordLoader const &recordLoader,
                      Protocol protocol = Protocol::JSON);

/**
 * \brief Append the new records or, per-record, new data, user defined content, curves/curve sets,
 *        and library data of a Sina Document to an existing JSON file. For a full explanation of behavior,
 *        see appendDocumentToHDF5. Note that a JSON file must be entirely re-written for an append to function;
 *        it's more efficient for shorter files, but prefer HDF5 for large files (many timeseries/long timeseries)
 *
 * \param jsonFilePath the path to the JSON file
 * \param newData a Sina Document containing the new data to append
 * \param mergeProtocol protocol for handling duplicate data/files/user_defined/etc
 *                            1 = take the new value. 2 = keep the old value. 3/Other = Cancel the append.
 *                     Protocol 1 is the conduit default behavior. Ignored entirely when skipping validation.
 * \param skipValidation whether to skip the validation step entirely. Most useful for well-controlled cases,
 *                       ex: a code is appending values to every timeseries every N cycles.
 * \return a conduit Node containing a list of any errors encountered in appending. If empty, success.
 */
conduit::Node appendDocumentToJson(const std::string &jsonFilePath,
                                   const Document &newData,
                                   const int mergeProtocol = 1,
                                   const bool skipValidation = false);

/**
 * \brief Append the new records or, per-record, new data, user defined content, curves/curve sets,
 *        and library data of a Sina Document to an existing HDF5 file. Appending is meant to allow a
 *        Sina document to grow as a simulation continues, i.e. by appending new snapshots into the same
 *        file, or by writing additional points of a timeseries to disk as they become available. HDF5 is
 *        the preferred filetype for heavy appending into, as it does not require being entirely re-written
 *        to file each time, unlike JSON.
 * 
 *        When appending a new Document to one already written on disk, it will first find any records that
 *        don't exist in the file on disk and add them. Then, for any file it finds that does exist, it will
 *        go through its various fields and compare them. New data/curve sets/etc. will simply be added in.
 *        For fields that already exist:
 *         - For curves in curve sets, new values will be appended onto the end. This is useful for ex: codes that
 *           don't store curves in memory, appending the new set of curve values each cycle will allow the file on
 *           disk to grow as the simulation continues.
 *        - For everything else (data, files, user_defined...), behavior depends on mergeProtocol, with the default
 *        behavior being to overwrite the old (on-disk) field with the new value.
 *
 * \param hdf5FilePath the path to the HDF5 file
 * \param new_data a vector of the new data to append
 * \param mergeProtocol protocol for handling duplicate data/files/user_defined/etc
 *                            1 = take the new value. 2 = keep the old value. 3 = Cancel the append.
 *                     Protocol 1 is the conduit default behavior. Ignored entirely when skipping validation.
 * \param skipValidation whether to skip the validation step entirely. Most useful for well-controlled cases,
 *                       ex: a code is appending values to every timeseries every N cycles.
 * 
 * \return a conduit Node containing a list of any errors encountered in appending. If empty, success!
 */
conduit::Node appendDocumentToHDF5(const std::string &hdf5FilePath,
                                   Document const &newData,
                                   const int mergeProtocol = 1,
                                   const bool skipValidation = false);

/**
 * @brief Append a Document to an existing file with automatic format detection
 * 
 * Automatically detects the file format from the extension and appends accordingly.
 * The format can optionally be overridden.
 * 
 * \param document The Document to append
 * \param filepath Path to the existing file
 * \param mergeProtocol How to handle conflicts (1=KEEP_ORIGINAL, 2=OVERWRITE, 3=ERROR)
 * \param outputProtocol Optional format override (default: AUTO_DETECT)
 * \throws std::runtime_error If the file cannot be opened or the format is unsupported
 */
void appendDocument(const Document &document,
                    const std::string &filepath,
                    int mergeProtocol = 1,
                    Protocol Protocol = Protocol::AUTO_DETECT);

/**
 * \brief Check a node against some file handle and return a Conduit node populated with any errors that
 *        would prevent appending. Primarily useful for dry run testing of what may be going wrong with an append.
 * 
 * \param appendTo Either a conduit node or a conduit relay filehandle (prefer the latter for HDF5) representing
 *                 the data on disk
 * \param appendFrom A conduit node representing the document being appended to the file on disk
 * \param endpoint The endpoint within the record at which the node should be checked. This is primarily useful for library_data
 * \param mergeProtocol protocol for handling duplicate data/files/user_defined/etc. See append methods above.
 * \param record_num The offset of the record within appendFrom, i.e. if you'd like to test the append against the first record, pass 0
 * \param original_file_path The path to the appendFrom document on disk, required for fast lookup of some HDF5 info
 * \return a conduit Node containing a list of any errors encountered in appending. If empty, success.
 */
template <typename ConduitRelayLike>
conduit::Node validateAppendDocument(ConduitRelayLike &appendTo,
                                     const conduit::Node &appendFrom,
                                     const std::string &endpoint,
                                     const int mergeProtocol,
                                     int record_num,
                                     // default'd because it might go away with a conduit update
                                     const std::string &original_file_path = "");

}  // namespace sina
}  // namespace axom

#endif  //SINA_DOCUMENT_HPP
