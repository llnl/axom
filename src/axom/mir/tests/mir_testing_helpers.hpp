// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_TESTING_HELPERS_HPP_
#define AXOM_MIR_TESTING_HELPERS_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include <conduit/conduit.hpp>
#include <conduit/conduit_relay_io.hpp>
#include <conduit/conduit_relay_io_blueprint.hpp>
#include <string>
#include <vector>

//------------------------------------------------------------------------------
// clang-format off
using seq_exec = axom::SEQ_EXEC;

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  using omp_exec = axom::OMP_EXEC;
#else
  using omp_exec = seq_exec;
#endif

#if defined (AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_CUDA)
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
  #else
    using cuda_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_HIP)
    constexpr int HIP_BLOCK_SIZE = 64;
    using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
  #else
    using hip_exec = seq_exec;
  #endif
#endif
// clang-format on

//------------------------------------------------------------------------------
// Define better names for the execution spaces. This header needs to be included
// after the ExecSpace types are defined.
template <typename ExecSpace>
struct execution_name
{
  static std::string name() { return "seq"; }
};

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_OPENMP)
template <>
struct execution_name<omp_exec>
{
  static std::string name() { return "omp"; }
};
  #endif
  #if defined(AXOM_USE_CUDA)
template <>
struct execution_name<cuda_exec>
{
  static std::string name() { return "cuda"; }
};
  #endif
  #if defined(AXOM_USE_HIP)
template <>
struct execution_name<hip_exec>
{
  static std::string name() { return "hip"; }
};
  #endif
#endif

//------------------------------------------------------------------------------
std::string pjoin(const std::string &str) { return str; }

std::string pjoin(const char *str) { return std::string(str); }

template <typename... Args>
std::string pjoin(const std::string &str, Args... args)
{
  return axom::utilities::filesystem::joinPath(str, pjoin(args...));
}

template <typename... Args>
std::string pjoin(const char *str, Args... args)
{
  return axom::utilities::filesystem::joinPath(std::string(str), pjoin(args...));
}

void psplit(const std::string &filepath, std::string &path, std::string &filename)
{
  axom::Path p(filepath);
  path = p.dirName();
  filename = p.baseName();
}

std::string dataDirectory() { return AXOM_DATA_DIR; }

std::string testData(const std::string &filename) { return pjoin(dataDirectory(), filename); }

std::string baselineDirectory();

std::string yamlRoot(const std::string &filepath)
{
  std::string retval, path, filename;
  psplit(filepath, path, filename);
  auto idx = filename.rfind(".");
  if(idx != std::string::npos)
  {
    retval = filename.substr(0, idx);
  }
  else
  {
    retval = filename;
  }
  return retval;
}

void printNode(const conduit::Node &n)
{
  conduit::Node options;
  options["num_children_threshold"] = 10000;
  options["num_elements_threshold"] = 10000;
  n.to_summary_string_stream(std::cout, options);
}

bool compareConduit(const conduit::Node &n1,
                    const conduit::Node &n2,
                    double tolerance,
                    conduit::Node &info)
{
  bool same = true;
  if(n1.dtype().id() == n2.dtype().id() && n1.dtype().is_floating_point())
  {
    const auto a1 = n1.as_double_accessor();
    const auto a2 = n2.as_double_accessor();
    double maxdiff = 0.;
    for(int i = 0; i < a1.number_of_elements() && same; i++)
    {
      double diff = fabs(a1[i] - a2[i]);
      maxdiff = std::max(diff, maxdiff);
      same &= diff <= tolerance;
      if(!same)
      {
        info.append().set(axom::fmt::format("\"{}\" fields differ at index {}.", n1.name(), i));
      }
    }
    info["maxdiff"][n1.name()] = maxdiff;
  }
  else
  {
    for(int i = 0; i < n1.number_of_children() && same; i++)
    {
      const auto &n1c = n1.child(i);
      const auto &n2c = n2.fetch_existing(n1c.name());
      same &= compareConduit(n1c, n2c, tolerance, info);
    }
  }
  return same;
}

void saveBaseline(const std::string &filename, const conduit::Node &n)
{
  std::string file_with_ext(filename + ".yaml");
  try
  {
    SLIC_INFO(axom::fmt::format("Save baseline {}", file_with_ext));
    conduit::relay::io::save(n, file_with_ext, "yaml");

#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
    SLIC_INFO(axom::fmt::format("Save visualization files..."));
    conduit::relay::io::blueprint::save_mesh(n, filename + "_hdf5", "hdf5");
      //axom::mir::utilities::blueprint::save_vtk(n, filename + "_vtk.vtk");
#endif
  }
  catch(...)
  {
    SLIC_INFO(axom::fmt::format("Could not save baseline to {}!", file_with_ext));

    printNode(n);

    // Check the data for errors.
    conduit::Node info;
    if(!conduit::blueprint::mesh::verify(n, info))
    {
      printNode(info);
    }
  }
}

void saveBaseline(const std::vector<std::string> &baselinePaths,
                  const std::string &baselineName,
                  const conduit::Node &n)
{
  for(const auto &path : baselinePaths)
  {
    axom::utilities::filesystem::makeDirsForPath(path);
    std::string filename(pjoin(path, baselineName));
    saveBaseline(filename, n);
  }
}

bool loadBaseline(const std::string &filename, conduit::Node &n)
{
  bool loaded = false;
  std::string file_with_ext(filename + ".yaml");
  //SLIC_INFO(axom::fmt::format("Load baseline {}", file_with_ext));
  if(axom::utilities::filesystem::pathExists(file_with_ext))
  {
    conduit::relay::io::load(file_with_ext, "yaml", n);
    loaded = true;
  }
  return loaded;
}

template <typename ExecSpace>
std::vector<std::string> baselinePaths()
{
  std::vector<std::string> paths;
  paths.push_back(pjoin(baselineDirectory(), execution_name<ExecSpace>::name()));
  paths.push_back(baselineDirectory());
  return paths;
}

bool compareBaseline(const std::vector<std::string> &baselinePaths,
                     const std::string &baselineName,
                     const conduit::Node &current,
                     double tolerance = 1.5e-6)
{
  bool success = false;
  int count = 0;
  for(const auto &path : baselinePaths)
  {
    try
    {
      // Load the baseline file.
      conduit::Node info, baselineNode;
      std::string filename(pjoin(path, baselineName));
      if(loadBaseline(filename, baselineNode))
      {
        // Compare the baseline to the current DC.
        SLIC_INFO(axom::fmt::format("Comparing to baseline {}", filename));
        success = compareConduit(baselineNode, current, tolerance, info);
        count++;
        if(!success)
        {
          info.print();

          std::string errFile(filename + "_err");
#if defined(AXOM_USE_HDF5)
          conduit::relay::io::blueprint::save_mesh(current, errFile, "hdf5");
#endif
          conduit::relay::io::blueprint::save_mesh(current, errFile + "_yaml", "yaml");
        }
        // We found a baseline so we can exit
        break;
      }
    }
    catch(...)
    {
      SLIC_INFO(axom::fmt::format("Could not load {} from {}!", baselineName, path));
    }
  }
  if(!success && count == 0)
  {
    SLIC_INFO(axom::fmt::format("No baselines found for {}", baselineName));
  }
  return success;
}

//------------------------------------------------------------------------------
template <typename Container1, typename Container2>
bool compare_views(const Container1 &a, const Container2 &b)
{
  bool eq = a.size() == b.size();
  for(axom::IndexType i = 0; i < a.size() && eq; i++)
  {
    eq &= a[i] == b[i];
  }
  if(!eq)
  {
    axom::fmt::format("a={{{}}}\nb={{{}}}", axom::fmt::join(a, ","), axom::fmt::join(b, ","));
  }
  return eq;
}
#endif
