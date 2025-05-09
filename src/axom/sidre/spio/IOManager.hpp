// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file IOManager.hpp
 *
 * \brief   Header file containing definition of IOManager class.
 *
 ******************************************************************************
 */

#ifndef SIDRE_IOMANAGER_HPP_
#define SIDRE_IOMANAGER_HPP_

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/sidre/core/Group.hpp"

// Spio headers
#include "axom/sidre/spio/IOBaton.hpp"

#include "mpi.h"

namespace axom
{
namespace sidre
{
/*!
 * \class IOManager
 *
 * \brief IOManager manages and organizes the I/O operations.
 *
 * This class handles the bookkeeping and organizing tasks that must be done
 * before calling Group's I/O methods.  It uses IOBaton to control the
 * parallel I/O operations, such that one rank at a time interacts with any
 * particular output file.
 */
class IOManager
{
public:
  /*!
   * \brief Constructor
   *
   * \param com               MPI communicator
   * \param use_scr           Use SCR library for scalable I/O management.
   *                          If true, the calling code must have already
   *                          called SCR_Init() after MPI_Init().
   */
  IOManager(MPI_Comm com, bool use_scr = false);

  /*!
   * \brief Destructor
   */
  ~IOManager();

  /*!
   * \brief write a Group to output files
   *
   * The Group, including all of its child groups and views, is written
   * to files according to the given protocol.
   *
   * This is an MPI collective call that must be called on all ranks in the
   * communicator used in this object's constructor.
   *
   * valid protocols:
   *
   *    sidre_hdf5
   *    sidre_conduit_json
   *    sidre_json
   *
   *    conduit_hdf5
   *    conduit_bin
   *    conduit_json
   *    json
   *
   * The output will consist of a root file and a subdirectory containing
   * data output files. The name of the root file will be the file_base
   * parameter with the suffix ".root" appended, while the subdirectory
   * will have the file_base string as its name. Exception:  If the file_base
   * parameter is provided as a string ending in ".root", that string will be
   * unchanged and used as the root file's name, while the subdirectory
   * name will be the file_base string with the ".root" suffix removed from
   * its end.
   *
   * The file_base string must be identical on all ranks making this call.
   * If it is not identical, a warning will be printed and the files produced
   * here may not be useable in a subsequent call to IOManager::read().
   *
   * \note The sidre_hdf5 and conduit_hdf5 protocols are only available
   * when Axom is configured with hdf5.
   *
   * \param group         Group to write to output
   * \param num_files     number of output data files
   * \param file_base     base name for output files
   * \param protocol      identifies I/O protocol
   * \param tree_pattern  Optional tree pattern string placed in root file,
   *                      to set search path for data in the output files
   *
   */
  void write(sidre::Group* group,
             int num_files,
             const std::string& file_base,
             const std::string& protocol,
             const std::string& tree_pattern = "datagroup");

  /*!
   * \brief write additional group to existing root file
   *
   * Should be called after write().  The native layout of the group will
   * be added to the root file.
   *
   * This may be called more than once to write multiple groups to the file.
   *
   * This is not an MPI collective call.  It writes one group from one rank to
   * one root file.
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.
   *
   * \param group         Group to add to root file
   * \param file_name     name of existing root file
   */
  void writeGroupToRootFile(sidre::Group* group, const std::string& file_name);

  /*!
   * \brief write additional group to a path inside an existing root file
   *
   * This should be called on an existing root file that already contains a
   * Group, most likely from a previous call to writeGroupToRootFile.
   *
   * group_path is a path that specifies a location somewhere in the native
   * layout of a Group already in the file.  The group added in this
   * method will be stored as a child group of the group specified by the
   * path.
   *
   * This is not an MPI collective call.  It writes one group from one rank to
   * one root file.
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.
   *
   * \param group         Group to add to root file
   * \param file_name     name of existing root file
   * \param group_path    path to a location within the root file
   */
  void writeGroupToRootFileAtPath(sidre::Group* group,
                                  const std::string& file_name,
                                  const std::string& group_path);

  /*!
   * \brief write additional group to a path inside an existing root file
   *
   * This should be called on an existing root file that already contains a
   * Group, most likely from a previous call to writeGroupToRootFile.
   *
   * group_path is a path that specifies a location somewhere in the native
   * layout of a Group already in the file.  The view added in this
   * method will be stored as a child view of the group specified by the
   * path.
   *
   * This is not an MPI collective call.  It writes one view from one rank to
   * one root file.
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.
   *
   * \param view          View to add to root file
   * \param file_name     name of existing root file
   * \param group_path    path to a location within the root file
   */
  void writeViewToRootFileAtPath(sidre::View* view,
                                 const std::string& file_name,
                                 const std::string& group_path);

  /*!
   * \brief write a Conduit Blueprint index to an existing root file
   *
   * Given a domain that adheres to Conduit's Blueprint format, this method
   * generates a Blueprint index and writes it to an existing root file.
   *
   * The domain must be stored in a Group located at the path in the
   * DataStore specified by domain_path argument.
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.  This must be called after calling write().
   *
   * The parameters domain_path and mesh_path are related. domain_path is the
   * path from the root of the DataStore to the domain that is being used to
   * generate the index.  mesh_path is the path from the group that was
   * passed into the preceding write() call to the same domain.  If write()
   * was called using the root group of the DataStore, then domain_path and
   * mesh_path will be identical.  Otherwise mesh_path is a sub-path of
   * domain_path.
   *
   * For example, the DataStore may contain a hierarchy of data that looks
   * like this, and we want to generate a blueprint index based on the mesh
   * located at "/hierarchy/domain_data/domain/blueprint_mesh":
   *
   * <root>
   * |--hierarchy
   * |  |--domain_data
   * |     |--domain
   * |     |  |--blueprint_mesh
   * |     |     |--coordsets
   * |     |     |  |--...
   * |     |     |--topologies
   * |     |     |  |--...
   * |     |     |--fields
   * |     |        |--...
   * |     |--...
   * |--
   *
   * If write() is called using the Group located at "/hierarchy/domain_data",
   * then only the Groups and Views descending from that Group are written
   * to the file.  To call this method, we would choose the full path in
   * the DataStore "hierarchy/domain_data/domain/blueprint_mesh" for
   * domain_path.  For the mesh_path argument, we choose only the path that
   * exists in the file:  "domain/blueprint_mesh".
   *
   * This is not an MPI collective call.  One rank writes a blueprint index
   * to one root file.
   *
   * \param datastore     DataStore containing Groups that hold domains
   *                      that adhere to the Blueprint format
   * \param domain_path   path in the DataStore to the domain that will be
   *                      used to generate a Blueprint index
   * \param file_name     name of existing root file
   * \param mesh_path     path in the data file to the domain that will be
   *                      used to generate a Blueprint index
   */
  void writeBlueprintIndexToRootFile(DataStore* datastore,
                                     const std::string& domain_path,
                                     const std::string& file_name,
                                     const std::string& mesh_path);

  /*!
   * \brief read from input file
   *
   * This is an MPI collective call on all ranks in the communicator used
   * to construct this object.  Calling code may also need to add
   * an MPI barrier after this call if invoking subsequent operations that
   * may change the input files.
   *
   * \param group         Group to fill with input data
   * \param root_file     root file containing input data
   * \param protocol      identifies I/O protocol
   * \param preserve_contents   Preserves group's existing contents if true
   */
  void read(sidre::Group* group,
            const std::string& root_file,
            const std::string& protocol,
            bool preserve_contents = false);

  /*!
   * \brief read from a root file
   *
   * This is an MPI collective call on all ranks in the communicator used
   * to construct this object.  Calling code may also need to add
   * an MPI barrier after this call if invoking subsequent operations that
   * may change the input files.
   *
   * \param group      Group to fill with input data
   * \param root_file  root file containing input data
   * \param preserve_contents   Preserves group's existing contents if true
   */
  void read(sidre::Group* group, const std::string& root_file, bool preserve_contents = false);

  /**
   * \brief Finds conduit relay protocol corresponding to a sidre protocol
   *
   * \param sidre_protocol String representing the sidre protocol
   * \return The conduit relay protocol corresponding to \a sidre_protocol
   * Options are: "hdf5", "json" and "conduit_json"
   * \see Group::save() for a list of valid sidre protocols
   */
  static std::string correspondingRelayProtocol(const std::string& sidre_protocol);

  /*!
   * \brief load external data into a group
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.
   *
   * The call to this method must follow a call to the IOManager::read with the
   * same group and root file.
   *
   * This is intended as the third step of the three step process to load
   * external data.  The first step is the call to IOManager::read; the
   * second step is to set valid pointers on all external views in the
   * hierarchy under the passed-in group; the third step is the call to
   * this method.
   *
   * \param group         Group to fill with external data from input
   * \param root_file     root file containing input data
   */
  void loadExternalData(sidre::Group* group, const std::string& root_file);

  /*!
   * \brief piecewise load of external data into a group
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.
   *
   * This is intended as an alternative way to load external data, in a
   * piecewise manner rather than loading all external data of a Sidre
   * hierarchy at once.  load_group is a Group somewhere in the hierarchy
   * under parent_group, and only the external data for views in the
   * subtree under load_group will be loaded.  In the second step of the
   * three step process, the calling code should set valid pointers on all
   * external views in the subtree under load_group.  The third step is
   * to call this method.  parent_group must be the same group that was
   * passed to IOManager::read with the same root file.
   *
   * This method may be called multiple times with different instances of
   * load_group representing different subtrees of the hierachy under
   * parent_group.
   *
   * \param parent_group  Group that was passed to IOManager::read
   * \param load_group    Group holding views to be filled with external data
   * \param root_file     root file containing input data
   */
  void loadExternalData(sidre::Group* parent_group,
                        sidre::Group* load_group,
                        const std::string& root_file);

  /*!
   * \brief gets the number of files in the dataset from the specified root file
   */
  int getNumFilesFromRoot(const std::string& root_file);

  /*!
   * \brief gets the number of groups in the dataset from the specified root file
   *
   * Usually this is the number of MPI ranks that wrote data to this set of files.
   */
  int getNumGroupsFromRoot(const std::string& root_file);

private:
  DISABLE_COPY_AND_ASSIGNMENT(IOManager);

  /*!
   * \brief create the root file for new output
   *
   * The root file will contain information about the organization of the
   * output subdirectory and files containing sidre data. It will be
   * named as the root_base string plus the suffix ".root", unless
   * root_base already ends with ".root", in which case no additional
   * suffix will be added.
   *
   * The returned string will be identical to root_base, except when
   * root_base is provided with the suffix ".root" already there, in
   * which case the returned string will be the root_base string with
   * the ".root" suffix removed from its end.
   */
  std::string createRootFile(const std::string& root_base,
                             int num_files,
                             const std::string& protocol,
                             const std::string& tree_pattern);

  std::string getProtocol(const std::string& root_name);

  /*!
   * Collective operation to get the file pattern from the root file.
   * The string is read on rank 0 and broadcast to the other ranks.
   * \note Works for all sidre protocols.
   */
  std::string getFilePatternFromRoot(const std::string& root_name, const std::string& protocol);

#ifdef AXOM_USE_HDF5
  std::string getHDF5FilePattern(const std::string& root_name);

  void readSidreHDF5(sidre::Group* group, const std::string& root_file, bool preserve_contents = false);
#endif /* AXOM_USE_HDF5 */

  std::string getFileNameForRank(const std::string& file_pattern,
                                 const std::string& root_name,
                                 int rankgroup_id);

  /*!
   * /brief Get a map of ranks to file ID numbers.
   *
   * This fills the given output View with a map identifying an integer ID
   * for the file that each rank will interact with during I/O.
   * The data in the map is an array indexed by rank with the values being
   * the file IDs
   */
  void getRankToFileMap(View* rank_to_file_map, int num_files);

  /*!
   * \brief If needed, get a file path created by SCR.
   *
   * When using this class with the SCR library, SCR must create a file path
   * to a storage location that it controls, where the actual I/O operations
   * will occur.  This private method invokes this SCR operation and returns
   * a string containing the SCR-controlled path.
   *
   * When SCR is not being used, either because the IOManager instance was
   * constructed with the SCR flag set to false, or because Sidre was not
   * built with SCR, the returned string is identical to the argument string.
   *
   * \param  path  The file path in the parallel file system known by the
   *               calling code.
   * \return       The path created by SCR for I/O, or a copy of the argument
   *               string when SCR is not being used.
   */
  std::string getSCRPath(const std::string& path);

  int m_comm_size;  // num procs in the MPI communicator
  int m_my_rank;    // rank of this proc

  IOBaton* m_baton;

  MPI_Comm m_mpi_comm;

  bool m_use_scr;
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_IOMANAGER_HPP_ */
