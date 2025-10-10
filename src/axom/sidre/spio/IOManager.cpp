// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Associated header file
#include "IOManager.hpp"

// Other axom headers
#include "axom/core/Macros.hpp"
#include "axom/core/Path.hpp"
#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

// SiDRe project headers
#include "axom/sidre/core/Buffer.hpp"
#include "axom/sidre/core/Group.hpp"
#include "axom/sidre/core/DataStore.hpp"
#include "axom/sidre/core/SidreTypes.hpp"
#include "axom/fmt.hpp"

// Conduit headers
#include "conduit_relay.hpp"
#include "conduit_relay_mpi.hpp"

#ifdef AXOM_USE_HDF5
  #include "conduit_relay_io_hdf5.hpp"
  #include "hdf5.h"
#endif

#ifdef AXOM_USE_SCR
  #include "scr.h"
#endif

namespace
{
/*!
 *  Utility function to broadcast a string from rank 0 to all other ranks
 */
std::string broadcastString(const std::string& str, MPI_Comm comm, int rank)
{
  // Find the string size and send to all ranks
  int buf_size = 0;
  if(rank == 0)
  {
    buf_size = str.size() + 1;
  }
  MPI_Bcast(&buf_size, 1, MPI_INT, 0, comm);

  // Allocate a buffer, copy and broadcast
  char* buf = new char[buf_size];
  if(rank == 0)
  {
    strcpy(buf, str.c_str());
  }
  MPI_Bcast(buf, buf_size, MPI_CHAR, 0, comm);

  // Create a string, delete the buffer and return
  std::string res(buf);
  delete[] buf;
  buf = nullptr;

  return res;
}

}  // end anonymous namespace

namespace axom
{
namespace sidre
{
/*
 *************************************************************************
 *
 * Create manager
 *
 *************************************************************************
 */
IOManager::IOManager(MPI_Comm comm, bool use_scr)
  : m_comm_size(1)
  , m_my_rank(0)
  , m_baton(nullptr)
  , m_mpi_comm(comm)
  , m_use_scr(use_scr)
{
  MPI_Comm_size(comm, &m_comm_size);
  MPI_Comm_rank(comm, &m_my_rank);
#ifndef AXOM_USE_SCR
  if(m_use_scr)
  {
    SLIC_WARNING(
      "IOManager constructor called with use_scr = true, but Axom was not "
      "compiled with SCR. IOManager will operate without SCR.");
  }
  m_use_scr = false;
#endif
}

/*
 *************************************************************************
 *
 * Destroy all contents
 *
 *************************************************************************
 */
IOManager::~IOManager()
{
  if(m_baton)
  {
    delete m_baton;
  }
}

/*
 ******************************************************************************
 *
 * Returns the conduit relay protocol corresponding to the given sidre protocol
 *
 ******************************************************************************
 */
std::string IOManager::correspondingRelayProtocol(const std::string& sidre_protocol)
{
#ifdef AXOM_USE_HDF5
  const std::string DEFAULT_PROTOCOL = "hdf5";
#else
  const std::string DEFAULT_PROTOCOL = "json";
#endif

  if(sidre_protocol == "sidre_hdf5" || sidre_protocol == "conduit_hdf5")
  {
    return "hdf5";
  }
  else if(sidre_protocol == "sidre_json" || sidre_protocol == "conduit_bin" ||
          sidre_protocol == "json")
  {
    return "json";
  }
  else if(sidre_protocol == "sidre_conduit_json" || sidre_protocol == "conduit_json")
  {
    return "conduit_json";
  }

  SLIC_WARNING("'" << sidre_protocol << "' is not a valid sidre protocol.");
  return DEFAULT_PROTOCOL;
}

/*
 *************************************************************************
 *
 * Write to file.
 *
 *************************************************************************
 */
void IOManager::write(sidre::Group* datagroup,
                      int num_files,
                      const std::string& file_base,
                      const std::string& protocol,
                      const std::string& tree_pattern)
{
  if(m_baton)
  {
    if(m_baton->getNumFiles() != num_files)
    {
      delete m_baton;
      m_baton = nullptr;
    }
  }

  if(!m_baton)
  {
    m_baton = new IOBaton(m_mpi_comm, num_files, m_comm_size);
  }

  SLIC_ERROR_IF(m_use_scr && num_files != m_comm_size, "SCR requires a file per process");

  std::string test_file_base(broadcastString(file_base, m_mpi_comm, m_my_rank));

  SLIC_WARNING_IF(test_file_base != file_base,
                  "IOManager::write() file_base argument is not identical "
                    << "on all ranks. This may cause the output files to be "
                    << "incompatible with a call to IOManager::read()."
                    << axom::fmt::format("\n\tfile_base: '{}'", file_base)
                    << axom::fmt::format("\n\ttest_file_base: '{}'", test_file_base));

  std::string output_base = createRootFile(file_base, num_files, protocol, tree_pattern);
  MPI_Barrier(m_mpi_comm);

  std::string root_name = output_base + ".root";

  if(protocol == "sidre_hdf5")
  {
#ifdef AXOM_USE_HDF5
    std::string file_pattern = getHDF5FilePattern(root_name);

    int set_id = m_baton->wait();

    std::string hdf5_name = getFileNameForRank(file_pattern, root_name, set_id);

    hdf5_name = getSCRPath(hdf5_name);

    hid_t h5_file_id, h5_group_id;
    if(m_baton->isFirstInGroup())
    {
      // no need to create directories in SCR
      if(!m_use_scr)
      {
        std::string dir_name;
        utilities::filesystem::getDirName(dir_name, hdf5_name);
        if(!dir_name.empty())
        {
          utilities::filesystem::makeDirsForPath(dir_name);
        }
      }
      h5_file_id = conduit::relay::io::hdf5_create_file(hdf5_name);
    }
    else
    {
      h5_file_id = conduit::relay::io::hdf5_open_file_for_read_write(hdf5_name);
    }
    SLIC_ASSERT(h5_file_id >= 0);

    std::string group_name = "datagroup";
    if(m_comm_size != num_files)
    {
      group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
    }
    h5_group_id = H5Gcreate(h5_file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    SLIC_ASSERT(h5_group_id >= 0);

    datagroup->save(h5_group_id);
    herr_t status;
    AXOM_UNUSED_VAR(status);

    status = H5Gclose(h5_group_id);
    SLIC_ASSERT(status >= 0);
    status = H5Fflush(h5_file_id, H5F_SCOPE_LOCAL);
    SLIC_ASSERT(status >= 0);
    status = H5Fclose(h5_file_id);
    SLIC_ASSERT(status >= 0);
#else
    SLIC_WARNING("'sidre_hdf5' protocol only available " << "when axom is configured with hdf5");
#endif /* AXOM_USE_HDF5 */
  }
  else
  {
    int set_id = m_baton->wait();
    std::string file_name = fmt::sprintf("%s_%07d", file_base, set_id);

    std::string obase = file_name + "." + protocol;
    datagroup->save(obase, protocol);
  }
  (void)m_baton->pass();

  MPI_Barrier(m_mpi_comm);
}

/*
 *************************************************************************
 *
 * Read from files
 *
 *************************************************************************
 */
void IOManager::read(sidre::Group* datagroup,
                     const std::string& root_file,
                     const std::string& protocol,
                     bool preserve_contents)
{
  MPI_Barrier(m_mpi_comm);

  if(protocol == "sidre_hdf5")
  {
#ifdef AXOM_USE_HDF5
    readSidreHDF5(datagroup, root_file, preserve_contents);
#else
    SLIC_WARNING("'sidre_hdf5' protocol only available " << "when axom is configured with hdf5");
#endif /* AXOM_USE_HDF5 */
  }
  else
  {
    if(m_baton)
    {
      if(m_baton->getNumFiles() != 1)
      {
        delete m_baton;
        m_baton = nullptr;
      }
    }

    if(!m_baton)
    {
      m_baton = new IOBaton(m_mpi_comm, m_comm_size, m_comm_size);
    }

    std::string file_pattern = getFilePatternFromRoot(root_file, protocol);

    int set_id = m_baton->wait();

    std::string file_name = getFileNameForRank(file_pattern, root_file, set_id);

    datagroup->load(file_name, protocol, preserve_contents);

    (void)m_baton->pass();
  }
}

/*
 *************************************************************************
 *
 * Read based on HDF5 root file.
 *
 *************************************************************************
 */
void IOManager::read(sidre::Group* datagroup, const std::string& root_file, bool preserve_contents)
{
  MPI_Barrier(m_mpi_comm);
  std::string protocol = getProtocol(root_file);
  read(datagroup, root_file, protocol, preserve_contents);
}

std::string IOManager::getSCRPath(const std::string& path)
{
#ifdef AXOM_USE_SCR
  if(m_use_scr)
  {
    char scr_name[SCR_MAX_FILENAME];
    std::string scr_path;
    if(SCR_Route_file(path.c_str(), scr_name) == SCR_SUCCESS)
    {
      // SCR gave us a path to use to open this file
      scr_path = scr_name;
    }
    else
    {
      SLIC_WARNING("SCR routing of path "
                   << path << " was unsuccessful. Attempting to proceed without routing.");
      scr_path = path;
    }
    return scr_path;
  }
  else
  {
    // SCR not enabled, return the path unchanged
    return path;
  }
#else
  // SCR support not compiled in, return path unchanged
  return path;
#endif
}

void IOManager::loadExternalData(sidre::Group* datagroup, const std::string& root_file)
{
  int num_files = getNumFilesFromRoot(root_file);
  int num_groups = getNumGroupsFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);
  SLIC_ASSERT(num_groups > 0);

  if(m_baton)
  {
    if(m_baton->getNumFiles() != num_files)
    {
      delete m_baton;
      m_baton = nullptr;
    }
  }

  if(!m_baton)
  {
    m_baton = new IOBaton(m_mpi_comm, num_files, num_groups);
  }

#ifdef AXOM_USE_HDF5
  std::string file_pattern = getHDF5FilePattern(root_file);

  int set_id = m_baton->wait();

  if(num_groups <= m_comm_size)
  {
    if(m_my_rank < num_groups)
    {
      herr_t errv;
      AXOM_UNUSED_VAR(errv);

      std::string hdf5_name = getFileNameForRank(file_pattern, root_file, set_id);

      hdf5_name = getSCRPath(hdf5_name);

      hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
      SLIC_ASSERT(h5_file_id >= 0);

      std::string group_name = "datagroup";
      if(H5Lexists(h5_file_id, group_name.c_str(), 0) <= 0)
      {
        group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
      }

      hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
      SLIC_ASSERT(h5_group_id >= 0);

      datagroup->loadExternalData(h5_group_id);

      errv = H5Gclose(h5_group_id);
      SLIC_ASSERT(errv >= 0);

      errv = H5Fclose(h5_file_id);
      SLIC_ASSERT(errv >= 0);
    }
  }
  else
  {
    for(int input_rank = m_my_rank; input_rank < num_groups; input_rank += m_comm_size)
    {
      herr_t errv;
      AXOM_UNUSED_VAR(errv);

      std::string hdf5_name = getFileNameForRank(file_pattern, root_file, input_rank);

      hdf5_name = getSCRPath(hdf5_name);

      hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
      SLIC_ASSERT(h5_file_id >= 0);

      std::string group_name = "datagroup";
      if(H5Lexists(h5_file_id, group_name.c_str(), 0) <= 0)
      {
        group_name = fmt::sprintf("datagroup_%07d", input_rank);
      }

      hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
      SLIC_ASSERT(h5_group_id >= 0);

      std::string input_name = fmt::sprintf("rank_%07d/sidre_input", input_rank);
      Group* one_rank_input = datagroup->getGroup(input_name);

      one_rank_input->loadExternalData(h5_group_id);

      errv = H5Gclose(h5_group_id);
      SLIC_ASSERT(errv >= 0);

      errv = H5Fclose(h5_file_id);
      SLIC_ASSERT(errv >= 0);
    }
  }

#else
  AXOM_UNUSED_VAR(datagroup);
  SLIC_WARNING("Loading external data only only available " << "when Axom is configured with hdf5");
#endif /* AXOM_USE_HDF5 */

  (void)m_baton->pass();
}

void IOManager::loadExternalData(sidre::Group* parent_group,
                                 sidre::Group* load_group,
                                 const std::string& root_file)
{
  int num_files = getNumFilesFromRoot(root_file);
  int num_groups = getNumGroupsFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);
  SLIC_ASSERT(num_groups > 0);

  if(m_baton)
  {
    if(m_baton->getNumFiles() != num_files)
    {
      delete m_baton;
      m_baton = nullptr;
    }
  }

  if(!m_baton)
  {
    m_baton = new IOBaton(m_mpi_comm, num_files, num_groups);
  }

#ifdef AXOM_USE_HDF5
  std::string file_pattern = getHDF5FilePattern(root_file);

  std::string parent_name(parent_group->getPathName());
  std::string load_name(load_group->getPathName());
  axom::Path parent_path(parent_name);
  axom::Path load_path(load_name);
  std::vector<std::string> parent_parts(parent_path.parts());
  std::vector<std::string> load_parts(load_path.parts());

  size_t parent_size = parent_parts.size();
  bool can_load = true;
  if(load_parts.size() < parent_size)
  {
    can_load = false;
  }

  for(size_t i = 0; i < parent_size; ++i)
  {
    if(parent_parts[i] != load_parts[i])
    {
      can_load = false;
    }
  }

  std::string subpath = load_name.substr(parent_name.size());
  char delimiter = parent_group->getPathDelimiter();
  if(!subpath.empty() && subpath[0] == delimiter)
  {
    subpath.erase(0, 1);
  }

  if(!subpath.empty() && !parent_group->hasGroup(subpath))
  {
    can_load = false;
  }

  int set_id = m_baton->wait();

  if(can_load)
  {
    if(num_groups <= m_comm_size)
    {
      if(m_my_rank < num_groups)
      {
        herr_t errv;
        AXOM_UNUSED_VAR(errv);

        std::string hdf5_name = getFileNameForRank(file_pattern, root_file, set_id);

        hdf5_name = getSCRPath(hdf5_name);

        hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
        SLIC_ASSERT(h5_file_id >= 0);

        std::string group_name = "datagroup";
        if(H5Lexists(h5_file_id, group_name.c_str(), 0) <= 0)
        {
          group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
        }

        hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
        SLIC_ASSERT(h5_group_id >= 0);

        load_group->loadExternalData(h5_group_id, subpath);

        errv = H5Gclose(h5_group_id);
        SLIC_ASSERT(errv >= 0);

        errv = H5Fclose(h5_file_id);
        SLIC_ASSERT(errv >= 0);
      }
    }
    else
    {
      for(int input_rank = m_my_rank; input_rank < num_groups; input_rank += m_comm_size)
      {
        herr_t errv;
        AXOM_UNUSED_VAR(errv);

        std::string hdf5_name = getFileNameForRank(file_pattern, root_file, input_rank);

        hdf5_name = getSCRPath(hdf5_name);

        hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
        SLIC_ASSERT(h5_file_id >= 0);

        std::string group_name = "datagroup";
        if(H5Lexists(h5_file_id, group_name.c_str(), 0) <= 0)
        {
          group_name = fmt::sprintf("datagroup_%07d", input_rank);
        }

        hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
        SLIC_ASSERT(h5_group_id >= 0);

        std::string input_name = fmt::sprintf("rank_%07d", input_rank);

        input_name = input_name + delimiter + "sidre_input" + delimiter + subpath;

        Group* one_rank_input = parent_group->getGroup(input_name);

        one_rank_input->loadExternalData(h5_group_id, subpath);

        errv = H5Gclose(h5_group_id);
        SLIC_ASSERT(errv >= 0);

        errv = H5Fclose(h5_file_id);
        SLIC_ASSERT(errv >= 0);
      }
    }
  }
  else
  {
    SLIC_WARNING("Path from parent group " << parent_group->getPathName() << " to group "
                                           << load_group->getPathName()
                                           << " was not found. No external data will be loaded.");
  }

#else
  AXOM_UNUSED_VAR(parent_group);
  AXOM_UNUSED_VAR(load_group);
  SLIC_WARNING("Loading external data only only available " << "when Axom is configured with hdf5");
#endif /* AXOM_USE_HDF5 */

  (void)m_baton->pass();
}

/*
 *************************************************************************
 *
 * Create a root file.
 *
 *************************************************************************
 */
std::string IOManager::createRootFile(const std::string& root_base,
                                      int num_files,
                                      const std::string& protocol,
                                      const std::string& tree_pattern)
{
  conduit::Node n;

  std::string file_base(utilities::string::removeSuffix(root_base, ".root"));

  if(m_my_rank == 0)
  {
    std::string root_file_name;
    std::string local_file_base;

    std::string relay_protocol = correspondingRelayProtocol(protocol);

    if(protocol == "sidre_hdf5" || protocol == "conduit_hdf5")
    {
#ifdef AXOM_USE_HDF5
      n["number_of_files"] = num_files;
      if(protocol == "sidre_hdf5")
      {
        std::string next;
        std::string delimiter = "/";
        axom::Path axom_file_path(file_base);
        local_file_base = axom_file_path.baseName();

        n["file_pattern"] = local_file_base + delimiter + local_file_base + "_" + "%07d.hdf5";
      }
      else
      {
        n["file_pattern"] = file_base + "_" + "%07d.conduit_hdf5";
      }
      n["number_of_trees"] = m_comm_size;

      n["tree_pattern"] = tree_pattern;
      n["protocol/name"] = protocol;
      n["protocol/version"] = "0.0";

      root_file_name = file_base + ".root";
#else
      SLIC_WARNING("IOManager::createRootFile() -- '" << protocol << "' protocol only available "
                                                      << "when Axom is configured with hdf5");
#endif /* AXOM_USE_HDF5 */
    }
    else
    {
      n["number_of_files"] = num_files;
      n["file_pattern"] = file_base + "_" + "%07d." + protocol;
      n["number_of_trees"] = m_comm_size;

      n["tree_pattern"] = tree_pattern;
      n["protocol/name"] = protocol;
      n["protocol/version"] = "0.0";

      root_file_name = file_base + ".root";
    }

    root_file_name = getSCRPath(root_file_name);

    conduit::relay::io::save(n, root_file_name, relay_protocol);
  }

  return file_base;
}

/*
 *************************************************************************
 *
 * Get protocol from root file
 *
 *************************************************************************
 */
std::string IOManager::getProtocol(const std::string& root_orig)
{
  std::string protocol;

  if(m_my_rank == 0)
  {
    std::string root_name = getSCRPath(root_orig);

    std::string extension;
    std::string base;
    std::string dot = ".";

    // Separate the first extension from the root file name
    // It should always be "root"
    conduit::utils::rsplit_string(root_name, dot, extension, base);

    SLIC_CHECK_MSG(
      extension == "root",
      "The root file name should always end in 'root'." << " File name was '" << root_name << "'");

    std::string relay_protocol = "json";
#ifdef AXOM_USE_HDF5
    // Attempt to open the root file using HDF5.  If it succeeds, set
    // relay_protocol to "hdf5", otherwise we assume a json protocol.
    //
    // Suppress error output for H5Fopen, since failure is acceptable here.

    H5E_auto2_t herr_func;
    void* old_client_data;
    H5Eget_auto(H5E_DEFAULT, &herr_func, &old_client_data);
    H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

    hid_t file_id = H5Fopen(root_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if(file_id > 0)
    {
      relay_protocol = "hdf5";
      herr_t errv = H5Fclose(file_id);
      AXOM_UNUSED_VAR(errv);
      SLIC_ASSERT(errv >= 0);
    }

    // Restore error output
    H5Eset_auto(H5E_DEFAULT, herr_func, old_client_data);

#endif

    conduit::Node n;
    conduit::relay::io::load(root_name, relay_protocol, n);

    conduit::Node& pnode = n["protocol/name"];
    if(pnode.dtype().is_string())
    {
      protocol = n["protocol/name"].as_string();
    }
    else if(relay_protocol == "json")
    {
      //If relay_protocol "json" didn't work, try "conduit_json"
      n.reset();
      conduit::relay::io::load(root_name, "conduit_json", n);
      protocol = n["protocol/name"].as_string();
    }

    if(protocol.empty())
    {
      // Did not find protocol name, issue warning and assign a default guess.

      SLIC_WARNING("'" << root_name << "/protocol/name' does not contain a valid Sidre protocol name.  "
                       << "Will attempt to use a default protocol.");

      if(relay_protocol == "hdf5")
      {
        protocol = "sidre_hdf5";
      }
      else
      {
        protocol = "sidre_json";
      }
    }
  }

  protocol = broadcastString(protocol, m_mpi_comm, m_my_rank);
  return protocol;
}

std::string IOManager::getFilePatternFromRoot(const std::string& root_name,
                                              const std::string& protocol)
{
  std::string file_pattern;
  if(m_my_rank == 0)
  {
    std::string root_path = getSCRPath(root_name);

    conduit::Node n;
    std::string relay_protocol = correspondingRelayProtocol(protocol);
    conduit::relay::io::load(root_path, relay_protocol, n);
    file_pattern = n["file_pattern"].as_string();
  }

  file_pattern = broadcastString(file_pattern, m_mpi_comm, m_my_rank);
  return file_pattern;
}

// Several private functions are only relevant when HDF5 is enabled
#ifdef AXOM_USE_HDF5

/*
 *************************************************************************
 *
 * Get namescheme pattern for a file holding Group data.
 *
 *************************************************************************
 */
std::string IOManager::getHDF5FilePattern(const std::string& root_name)
{
  std::string file_pattern;
  if(m_my_rank == 0)
  {
    std::string root_path = getSCRPath(root_name);

    conduit::Node n;
    conduit::relay::io::load(root_path + ":file_pattern", "hdf5", n);

    file_pattern = n.as_string();
  }
  file_pattern = broadcastString(file_pattern, m_mpi_comm, m_my_rank);

  return file_pattern;
}

/*
 *************************************************************************
 *
 * Read based on HDF5 root file.
 *
 *************************************************************************
 */
void IOManager::readSidreHDF5(sidre::Group* datagroup,
                              const std::string& root_file,
                              bool preserve_contents)
{
  int num_files = getNumFilesFromRoot(root_file);
  int num_groups = getNumGroupsFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);
  SLIC_ERROR_IF(num_groups > m_comm_size && num_groups != num_files,
                "IOManager attempted to read using a smaller number of processors "
                  << "than were used to produce the I/O files.  This only can work if "
                  << "those files were created in file-per-processor mode.");

  if(m_baton)
  {
    if(m_baton->getNumFiles() != num_files)
    {
      delete m_baton;
      m_baton = nullptr;
    }
  }

  if(!m_baton)
  {
    m_baton = new IOBaton(m_mpi_comm, num_files, num_groups);
  }

  std::string file_pattern = getHDF5FilePattern(root_file);

  int set_id = m_baton->wait();
  if(num_groups <= m_comm_size)
  {
    if(m_my_rank < num_groups)
    {
      herr_t errv;
      AXOM_UNUSED_VAR(errv);

      std::string hdf5_name = getFileNameForRank(file_pattern, root_file, set_id);

      hdf5_name = getSCRPath(hdf5_name);

      hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
      SLIC_ASSERT(h5_file_id >= 0);

      std::string group_name = "datagroup";
      if(H5Lexists(h5_file_id, group_name.c_str(), 0) <= 0)
      {
        group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
      }
      hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
      SLIC_ASSERT(h5_group_id >= 0);

      datagroup->load(h5_group_id, "sidre_hdf5", preserve_contents);

      errv = H5Gclose(h5_group_id);
      SLIC_ASSERT(errv >= 0);

      errv = H5Fclose(h5_file_id);
      SLIC_ASSERT(errv >= 0);
    }
  }
  else
  {
    datagroup->createViewScalar("reduced_input_ranks", num_groups);

    for(int input_rank = m_my_rank; input_rank < num_groups; input_rank += m_comm_size)
    {
      herr_t errv;
      AXOM_UNUSED_VAR(errv);

      std::string hdf5_name = getFileNameForRank(file_pattern, root_file, input_rank);

      hdf5_name = getSCRPath(hdf5_name);

      hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
      SLIC_ASSERT(h5_file_id >= 0);

      std::string group_name = "datagroup";
      if(H5Lexists(h5_file_id, group_name.c_str(), 0) <= 0)
      {
        group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
      }
      hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
      SLIC_ASSERT(h5_group_id >= 0);

      std::string input_name = fmt::sprintf("rank_%07d/sidre_input", input_rank);
      Group* one_rank_input = datagroup->createGroup(input_name);

      one_rank_input->load(h5_group_id, "sidre_hdf5", preserve_contents);

      errv = H5Gclose(h5_group_id);
      SLIC_ASSERT(errv >= 0);

      errv = H5Fclose(h5_file_id);
      SLIC_ASSERT(errv >= 0);
    }
  }
  (void)m_baton->pass();
}
#endif /* AXOM_USE_HDF5 */

/*
 *************************************************************************
 *
 * Get the file name based on a pattern and a group id
 *
 *************************************************************************
 */
std::string IOManager::getFileNameForRank(const std::string& file_pattern,
                                          const std::string& root_name,
                                          int rankgroup_id)
{
  std::string file_name = fmt::sprintf(file_pattern.c_str(), rankgroup_id);

  //If the root file was given as a path,
  //find the directory and add it to file_name
  char delimiter = '/';
  axom::Path root_path(root_name);
  std::string root_dir = root_path.dirName();

  file_name = axom::utilities::string::appendPrefix(root_dir, file_name, delimiter);

  return file_name;
}

void IOManager::getRankToFileMap(View* rank_to_file_map, int num_files)
{
  SLIC_ASSERT(rank_to_file_map != nullptr);

  if(m_baton)
  {
    if(m_baton->getNumFiles() != num_files)
    {
      delete m_baton;
      m_baton = nullptr;
    }
  }

  if(!m_baton)
  {
    m_baton = new IOBaton(m_mpi_comm, num_files, m_comm_size);
  }

  std::vector<int64_t> map_vec(m_comm_size, 0);

  int set_id = m_baton->wait();
  map_vec[m_my_rank] = static_cast<int64_t>(set_id);
  (void)m_baton->pass();

  conduit::Node map_local;
  map_local.set_external(&map_vec[0], map_vec.size());

  conduit::Node map_global;
  conduit::relay::mpi::max_all_reduce(map_local, map_global, m_mpi_comm);
  rank_to_file_map->importArrayNode(map_global);
}

/*
 *************************************************************************
 *
 * Query the rootfile for the number of input files.
 *
 *************************************************************************
 */
int IOManager::getNumFilesFromRoot(const std::string& root_file)
{
  /*
   * Read num_files from rootfile on rank 0.
   */
  int read_num_files = 0;
  if(m_my_rank == 0)
  {
    std::string root_path = getSCRPath(root_file);

    conduit::Node n;
    conduit::relay::io::load(root_path + ":number_of_files", "hdf5", n);
    read_num_files = n.to_int();
    SLIC_ASSERT(read_num_files > 0);
  }

  /*
   * Reduction sets num_files on all ranks.
   */
  int num_files;
  MPI_Allreduce(&read_num_files, &num_files, 1, MPI_INT, MPI_SUM, m_mpi_comm);
  SLIC_ASSERT(num_files > 0);

  return num_files;
}

int IOManager::getNumGroupsFromRoot(const std::string& root_file)
{
  /*
   * Read number_of_trees from rootfile on rank 0.
   */
  int read_num_trees = 0;
  if(m_my_rank == 0)
  {
    std::string root_path = getSCRPath(root_file);

    conduit::Node n;
    conduit::relay::io::load(root_path + ":number_of_trees", "hdf5", n);
    read_num_trees = n.to_int();
    SLIC_ASSERT(read_num_trees > 0);
  }

  /*
   * Reduction sets num_groups on all ranks.
   */
  int num_groups;
  MPI_Allreduce(&read_num_trees, &num_groups, 1, MPI_INT, MPI_SUM, m_mpi_comm);
  SLIC_ASSERT(num_groups > 0);

  return num_groups;
}

/*
 *************************************************************************
 *
 * Write a group to an existing root file
 *
 *************************************************************************
 */

void IOManager::writeGroupToRootFile(sidre::Group* group, const std::string& file_name)
{
#ifdef AXOM_USE_HDF5
  std::string tmp_name = getSCRPath(file_name);

  hid_t root_file_id = conduit::relay::io::hdf5_open_file_for_read_write(tmp_name);

  SLIC_ASSERT(root_file_id >= 0);

  hid_t group_id =
    H5Gcreate2(root_file_id, group->getName().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  SLIC_ASSERT(group_id >= 0);

  conduit::Node data_holder;
  group->createNativeLayout(data_holder);

  conduit::relay::io::hdf5_write(data_holder, group_id);

  herr_t errv;
  AXOM_UNUSED_VAR(errv);

  errv = H5Gclose(group_id);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fflush(root_file_id, H5F_SCOPE_LOCAL);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fclose(root_file_id);
  SLIC_ASSERT(errv >= 0);
#else
  AXOM_UNUSED_VAR(group);
  AXOM_UNUSED_VAR(file_name);

  SLIC_WARNING("Axom configured without hdf5. "
               << "IOManager::writeGroupToRootFile() only currently implemented "
               << "for 'sidre_hdf5' protocol. ");
#endif /* AXOM_USE_HDF5 */
}

/*
 *************************************************************************
 *
 * Write a group to an existing path inside a root file
 *
 *************************************************************************
 */

void IOManager::writeGroupToRootFileAtPath(sidre::Group* group,
                                           const std::string& file_name,
                                           const std::string& group_path)
{
#ifdef AXOM_USE_HDF5
  std::string tmp_name = getSCRPath(file_name);

  hid_t root_file_id = conduit::relay::io::hdf5_open_file_for_read_write(tmp_name);

  SLIC_ASSERT(root_file_id >= 0);

  hid_t path_id = H5Gopen(root_file_id, group_path.c_str(), 0);

  SLIC_ASSERT(path_id >= 0);

  hid_t group_id =
    H5Gcreate2(path_id, group->getName().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  SLIC_ASSERT(group_id >= 0);

  conduit::Node data_holder;
  group->createNativeLayout(data_holder);

  conduit::relay::io::hdf5_write(data_holder, group_id);

  herr_t errv;
  AXOM_UNUSED_VAR(errv);

  errv = H5Gclose(group_id);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fflush(root_file_id, H5F_SCOPE_LOCAL);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fclose(root_file_id);
  SLIC_ASSERT(errv >= 0);
#else
  AXOM_UNUSED_VAR(group);
  AXOM_UNUSED_VAR(file_name);
  AXOM_UNUSED_VAR(group_path);

  SLIC_WARNING("Axom configured without hdf5. "
               << "IOManager::writeGroupToRootFileAtPath() only currently "
               << "implemented for 'sidre_hdf5' protocol. ");
#endif /* AXOM_USE_HDF5 */
}

/*
 *************************************************************************
 *
 * Write a view to an existing path inside a root file
 *
 *************************************************************************
 */

void IOManager::writeViewToRootFileAtPath(sidre::View* view,
                                          const std::string& file_name,
                                          const std::string& group_path)
{
#ifdef AXOM_USE_HDF5
  std::string tmp_name = getSCRPath(file_name);

  hid_t root_file_id = conduit::relay::io::hdf5_open_file_for_read_write(tmp_name);

  SLIC_ASSERT(root_file_id >= 0);

  hid_t path_id = H5Gopen(root_file_id, group_path.c_str(), 0);

  SLIC_ASSERT(path_id >= 0);

  conduit::Node data_holder;
  view->createNativeLayout(data_holder[view->getName()]);

  conduit::relay::io::hdf5_write(data_holder, path_id);

  herr_t errv;
  AXOM_UNUSED_VAR(errv);

  errv = H5Fflush(root_file_id, H5F_SCOPE_LOCAL);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fclose(root_file_id);
  SLIC_ASSERT(errv >= 0);
#else
  AXOM_UNUSED_VAR(view);
  AXOM_UNUSED_VAR(file_name);
  AXOM_UNUSED_VAR(group_path);

  SLIC_WARNING("Axom configured without hdf5. "
               << "writeViewToRootFileAtPath() only currently implemented "
               << "for 'sidre_hdf5' protocol. ");
#endif /* AXOM_USE_HDF5 */
}

void IOManager::writeBlueprintIndexToRootFile(DataStore* datastore,
                                              const std::string& domain_path,
                                              const std::string& file_name,
                                              const std::string& mesh_path)
{
#ifdef AXOM_USE_HDF5
  std::string tmp_name = getSCRPath(file_name);

  hid_t root_file_id = conduit::relay::io::hdf5_open_file_for_read_write(tmp_name);

  AXOM_UNUSED_VAR(root_file_id);
  SLIC_ASSERT(root_file_id >= 0);

  std::string delimiter(1, datastore->getRoot()->getPathDelimiter());

  //The final name in mesh_path will be used as the name of the
  //blueprint index.
  axom::Path axom_mesh_path(mesh_path, delimiter[0]);
  std::string blueprint_name = axom_mesh_path.baseName();

  std::string bp_index = "blueprint_index" + delimiter + blueprint_name;

  bool multi_domain = false;
  if(m_comm_size > 1)
  {
    multi_domain = true;
  }
  else
  {
    Group* domain;
    if(domain_path == "/")
    {
      domain = datastore->getRoot();
    }
    else if(datastore->getRoot()->hasGroup(domain_path))
    {
      domain = datastore->getRoot()->getGroup(domain_path);
    }
    else
    {
      domain = nullptr;
    }

    if(domain && !domain->hasChildGroup("coordsets"))
    {
      multi_domain = true;
    }
  }

  bool success = false;

  if(multi_domain)
  {
    success = datastore->generateBlueprintIndex(MPI_COMM_WORLD, domain_path, mesh_path, bp_index);
  }
  else
  {
    success = datastore->generateBlueprintIndex(domain_path, mesh_path, bp_index, m_comm_size);
  }

  if(success)
  {
    Group* state_group = datastore->getRoot()->getGroup(bp_index)->getGroup("state");

    View* rank_to_file = state_group->createView("rank_to_file_map");
    getRankToFileMap(rank_to_file, m_baton->getNumFiles());

    if(state_group->hasGroup("partition_map"))
    {
      View* domain_to_rank = state_group->getView("partition_map/datagroup");

      View* domain_to_file = state_group->createViewAndAllocate("partition_map/file",
                                                                domain_to_rank->getTypeID(),
                                                                domain_to_rank->getNumElements());

      View* domain_ids = state_group->createViewAndAllocate("partition_map/domain",
                                                            domain_to_rank->getTypeID(),
                                                            domain_to_rank->getNumElements());

      int64_t* rank_file_map = rank_to_file->getArray();
      int64_t* domain_rank_map = domain_to_rank->getArray();
      int64_t* domain_file_map = domain_to_file->getArray();
      int64_t* domain_ids_array = domain_ids->getArray();

      for(IndexType i = 0; i < domain_to_file->getNumElements(); ++i)
      {
        domain_ids_array[i] = i;

        if(domain_rank_map[i] >= 0)
        {
          domain_file_map[i] = rank_file_map[domain_rank_map[i]];
        }
      }

      std::string file_pattern = getFilePatternFromRoot(file_name, "sidre_hdf5");

      std::string domain_pattern;
      if(datastore->getRoot()->hasView("domain_pattern") &&
         datastore->getRoot()->getView("domain_pattern")->isString())
      {
        domain_pattern = std::string(datastore->getRoot()->getView("domain_pattern")->getString());
      }

      std::string partition_pattern = file_pattern + "/datagroup";
      if(m_comm_size != m_baton->getNumFiles())
      {
        partition_pattern = partition_pattern + "_{datagroup:07d}";
      }
      if(!domain_pattern.empty())
      {
        partition_pattern = partition_pattern + "/" + domain_pattern;
      }

      state_group->createViewString("partition_pattern", partition_pattern);
    }

    if(m_my_rank == 0)
    {
      Group* ind_group = datastore->getRoot()->getGroup("blueprint_index");
      writeGroupToRootFile(ind_group, file_name);
    }
  }
  else
  {
    SLIC_WARNING("DataStore failed to generate Blueprint Index " << "based on group at path "
                                                                 << domain_path);
  }
#else
  AXOM_UNUSED_VAR(datastore);
  AXOM_UNUSED_VAR(domain_path);
  AXOM_UNUSED_VAR(file_name);
  AXOM_UNUSED_VAR(mesh_path);

  SLIC_WARNING("Axom configured without hdf5. "
               << "writeBlueprintIndexToRootFile() only currently implemented "
               << "for 'sidre_hdf5' protocol. ");
#endif /* AXOM_USE_HDF5 */
}

} /* end namespace sidre */
} /* end namespace axom */
