// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BLUEPRINT_TESTING_HELPERS_HPP_
#define AXOM_BLUEPRINT_TESTING_HELPERS_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include <conduit/conduit.hpp>
#include <conduit/conduit_relay_io.hpp>
#include <conduit/conduit_relay_io_blueprint.hpp>
#include <algorithm>
#include <string>
#include <vector>

//------------------------------------------------------------------------------
// clang-format off
using seq_exec = axom::SEQ_EXEC;

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  using omp_exec = axom::OMP_EXEC;
#else
  using omp_exec = seq_exec;
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  constexpr int CUDA_BLOCK_SIZE = 256;
  using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
#else
  using cuda_exec = seq_exec;
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  constexpr int HIP_BLOCK_SIZE = 64;
  using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
#else
  using hip_exec = seq_exec;
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

template <typename T>
struct compareValue
{
  static inline bool compare(T v1, T v2, T AXOM_UNUSED_PARAM(tolerance)) { return v1 == v2; }
};

template <>
struct compareValue<float>
{
  static inline bool compare(float v1, float v2, float tolerance)
  {
    float diff = axom::utilities::max(v1, v2) - axom::utilities::min(v1, v2);
    return diff <= tolerance;
  }
};

template <>
struct compareValue<double>
{
  static inline bool compare(double v1, double v2, double tolerance)
  {
    double diff = axom::utilities::abs(v1 - v2);
    return diff <= tolerance;
  }
};

template <typename T>
bool compareArray(const conduit::Node &n1,
                  const conduit::Node &AXOM_UNUSED_PARAM(n2),
                  const conduit::DataAccessor<T> &a1,
                  const conduit::DataAccessor<T> &a2,
                  conduit::Node &info,
                  T tolerance = T {0})
{
  bool same = true;
  if(a1.number_of_elements() != a2.number_of_elements())
  {
    info[n1.path()]["errors"] = axom::fmt::format("Different lengths. {} != {}",
                                                  a1.number_of_elements(),
                                                  a2.number_of_elements());
    same = false;
  }
  else
  {
    T maxdiff {0};
    conduit::Node errors;
    constexpr int errorLimit = 10;
    int errorCount = 0;
    for(int i = 0; i < a1.number_of_elements(); i++)
    {
      const T diff = axom::utilities::max(a1[i], a2[i]) - axom::utilities::min(a1[i], a2[i]);
      maxdiff = std::max(diff, maxdiff);
      if(!compareValue<T>::compare(a1[i], a2[i], tolerance))
      {
        if(errorCount < errorLimit)
        {
          errors.append().set(
            axom::fmt::format("Difference at index {}. ({} != {})", i, a1[i], a2[i]));
        }
        errorCount++;
        if(errorCount == errorLimit + 1)
        {
          errors.append().set("...");
        }
        same = false;
      }
    }
    if(errorCount > 0)
    {
      info[n1.path()]["maxdiff"] = maxdiff;
      info[n1.path()]["errors"].set(errors);
    }
  }
  return same;
}

template <typename T>
bool compareScalar(const conduit::Node &n1,
                   const conduit::Node &AXOM_UNUSED_PARAM(n2),
                   const T &v1,
                   const T &v2,
                   conduit::Node &info,
                   T tolerance = T {})
{
  bool same = compareValue<T>::compare(v1, v2, tolerance);
  if(!same)
  {
    info[n1.path()]["errors"].set(axom::fmt::format("{} != {}", v1, v2));
  }
  return same;
}

bool compareNode(const conduit::Node &n1, const conduit::Node &n2, double tolerance, conduit::Node &info)
{
  bool same = false;
  // String
  if(n1.dtype().is_string() && n2.dtype().is_string())
  {
    same = n1.as_string() == n2.as_string();
    if(!same)
    {
      info[n1.path()]["errors"].set(
        axom::fmt::format("\"{}\" != \"{}\"", n1.as_string(), n2.as_string()));
    }
  }
  else if(n1.dtype().number_of_elements() > 1 || n2.dtype().number_of_elements() > 1)
  {
    // Array comparison.
    if(n1.dtype().id() == n1.dtype().id())
    {
      // Types are equal
      if(n1.dtype().is_int8())
      {
        same =
          compareArray<conduit::int8>(n1, n2, n1.as_int8_accessor(), n2.as_int8_accessor(), info);
      }
      else if(n1.dtype().is_int16())
      {
        same =
          compareArray<conduit::int16>(n1, n2, n1.as_int16_accessor(), n2.as_int16_accessor(), info);
      }
      else if(n1.dtype().is_int32())
      {
        same =
          compareArray<conduit::int32>(n1, n2, n1.as_int32_accessor(), n2.as_int32_accessor(), info);
      }
      else if(n1.dtype().is_int64())
      {
        same =
          compareArray<conduit::int64>(n1, n2, n1.as_int64_accessor(), n2.as_int64_accessor(), info);
      }
      else if(n1.dtype().is_uint8())
      {
        same =
          compareArray<conduit::uint8>(n1, n2, n1.as_uint8_accessor(), n2.as_uint8_accessor(), info);
      }
      else if(n1.dtype().is_uint16())
      {
        same =
          compareArray<conduit::uint16>(n1, n2, n1.as_uint16_accessor(), n2.as_uint16_accessor(), info);
      }
      else if(n1.dtype().is_uint32())
      {
        same =
          compareArray<conduit::uint32>(n1, n2, n1.as_uint32_accessor(), n2.as_uint32_accessor(), info);
      }
      else if(n1.dtype().is_uint64())
      {
        same =
          compareArray<conduit::uint64>(n1, n2, n1.as_uint64_accessor(), n2.as_uint64_accessor(), info);
      }
      else if(n1.dtype().is_float32())
      {
        same = compareArray<conduit::float32>(n1,
                                              n2,
                                              n1.as_float32_accessor(),
                                              n2.as_float32_accessor(),
                                              info,
                                              static_cast<conduit::float32>(tolerance));
      }
      else if(n1.dtype().is_float64())
      {
        same = compareArray<conduit::float64>(n1,
                                              n2,
                                              n1.as_float64_accessor(),
                                              n2.as_float64_accessor(),
                                              info,
                                              tolerance);
      }
      else
      {
        info[n1.path()]["errors"].set(
          axom::fmt::format("Unsupported array type {}.", n1.dtype().name()));
      }
    }
    // Array comparison - types differ
    else if(n1.dtype().is_floating_point() && n2.dtype().is_floating_point())
    {
      same = compareArray<conduit::float64>(n1,
                                            n2,
                                            n1.as_double_accessor(),
                                            n2.as_double_accessor(),
                                            info,
                                            tolerance);
    }
    else
    {
      same =
        compareArray<conduit::index_t>(n1, n2, n1.as_index_t_accessor(), n2.as_index_t_accessor(), info);
    }
  }
  else
  {
    // Scalars.
    if(n1.dtype().is_int8())
    {
      same = compareScalar<conduit::int8>(n1, n2, n1.to_int8(), n2.to_int8(), info);
    }
    else if(n1.dtype().is_int16())
    {
      same = compareScalar<conduit::int16>(n1, n2, n1.to_int16(), n2.to_int16(), info);
    }
    else if(n1.dtype().is_int32())
    {
      same = compareScalar<conduit::int32>(n1, n2, n1.to_int32(), n2.to_int32(), info);
    }
    else if(n1.dtype().is_int64())
    {
      same = compareScalar<conduit::int64>(n1, n2, n1.to_int64(), n2.to_int64(), info);
    }
    else if(n1.dtype().is_uint8())
    {
      same = compareScalar<conduit::uint8>(n1, n2, n1.to_uint8(), n2.to_uint8(), info);
    }
    else if(n1.dtype().is_uint16())
    {
      same = compareScalar<conduit::uint16>(n1, n2, n1.to_uint16(), n2.to_uint16(), info);
    }
    else if(n1.dtype().is_uint32())
    {
      same = compareScalar<conduit::uint32>(n1, n2, n1.to_uint32(), n2.to_uint32(), info);
    }
    else if(n1.dtype().is_uint64())
    {
      same = compareScalar<conduit::uint64>(n1, n2, n1.to_uint64(), n2.to_uint64(), info);
    }
    else if(n1.dtype().is_float32())
    {
      same = compareScalar<conduit::float32>(n1,
                                             n2,
                                             n1.to_float32(),
                                             n2.to_float32(),
                                             info,
                                             static_cast<conduit::float32>(tolerance));
    }
    else if(n1.dtype().is_float64())
    {
      same =
        compareScalar<conduit::float64>(n1, n2, n1.to_float64(), n2.to_float64(), info, tolerance);
    }
    else
    {
      info[n1.path()]["errors"].set(
        axom::fmt::format("Error comparing \"{}\" and \"{}\"", n1.dtype().name(), n2.dtype().name()));
      same = false;
    }
  }
  return same;
}

bool compareConduit(const conduit::Node &n1,
                    const conduit::Node &n2,
                    double tolerance,
                    conduit::Node &info)
{
  bool same = true;
  // See if n1, n2 are objects - but not both.
  if((n1.dtype().is_object() && !n2.dtype().is_object()) ||
     (!n1.dtype().is_object() && n2.dtype().is_object()))
  {
    info[n1.path()]["errors"] =
      axom::fmt::format("Object types differ. \"{}\" is a {}. \"{}\" is a {}.",
                        n1.path(),
                        n1.dtype().name(),
                        n1.path(),
                        n1.dtype().name());
    same = false;
  }
  else if(n1.dtype().is_object() && n2.dtype().is_object())
  {
    // Both are objects. Recurse.
    for(int i = 0; i < n1.number_of_children(); i++)
    {
      const auto &n1c = n1.child(i);
      const auto &n2c = n2.fetch_existing(n1c.name());
      same &= compareConduit(n1c, n2c, tolerance, info);
    }
  }
  // Arrays
  else
  {
    same = compareNode(n1, n2, tolerance, info);
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
                     conduit::Node &info,
                     double tolerance = 1.5e-6)
{
  bool success = false;
  int count = 0;
  for(const auto &path : baselinePaths)
  {
    try
    {
      // Load the baseline file.
      conduit::Node baselineNode;
      std::string filename(pjoin(path, baselineName));
      if(loadBaseline(filename, baselineNode))
      {
        // Compare the baseline to the current DC.
        SLIC_INFO(axom::fmt::format("Comparing to baseline {}", filename));
        success = compareConduit(baselineNode, current, tolerance, info);
        count++;
        if(!success)
        {
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
    catch(conduit::Error &e)
    {
      SLIC_INFO(axom::fmt::format("Could not load {} from {}! {}", baselineName, path, e.message()));
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
    SLIC_INFO(
      axom::fmt::format("a={{{}}}\nb={{{}}}", axom::fmt::join(a, ","), axom::fmt::join(b, ",")));
  }
  return eq;
}

//------------------------------------------------------------------------------
namespace axom
{
namespace blueprint
{
namespace testing
{

/*!
 * \brief Base class for MIR test applications. It provides support for operations
 *        such as comparing against Conduit baseline files.
 */
class TestApplication
{
public:
  /*!
   * \brief Constructor.
   */
  TestApplication()
    : m_app()
    , m_annotationMode("none")
    , m_handler(false)
    , m_visualize(false)
    , m_verbose(true)
    , m_rebaseline_raw()
    , m_rebaseline()
  {
    m_rebaseline.push_back("none");
  }

  /*!
   * \brief Destructor
   */
  ~TestApplication() { }

  /*!
   * \brief Parse the command line and run the tests.
   */
  int execute(int argc, char *argv[])
  {
    int result = 0;

    try
    {
      // Define command line options.
#if defined(AXOM_USE_CALIPER)
      m_app.add_option("--caliper", m_annotationMode)
        ->description(
          "caliper annotation mode. Valid options include 'none' and 'report'. "
          "Use 'help' to see full list.")
        ->capture_default_str()
        ->check(axom::utilities::ValidCaliperMode);
#endif
      m_app.add_flag("--handler", m_handler, "Enable Conduit handler.");
      m_app.add_flag("--visualize", m_visualize, "Save visualization files.");
      m_app
        .add_option("--rebaseline",
                    m_rebaseline_raw,
                    "List of comma-separated test "
                    "names, or \"all\" if you want to rebaseline all tests.")
        ->expected(1);

      // Parse command line options.
      m_app.parse(argc, argv);

      // More initialization.
      initialize();
#if defined(AXOM_USE_CALIPER)
      axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(m_annotationMode);
#endif
      axom::slic::SimpleLogger logger;
      if(m_handler)
      {
        conduit::utils::set_error_handler(conduit_debug_err_handler);
      }

      // Run all the tests.
      result = RUN_ALL_TESTS();
    }
    catch(axom::CLI::CallForHelp &e)
    {
      std::cout << m_app.help() << std::endl;
      result = 0;
    }
    catch(axom::CLI::ParseError &e)
    {
      // Handle other parsing errors
      std::cerr << e.what() << std::endl;
      result = m_app.exit(e);
    }
    return result;
  }

  /*!
   * \brief Save a mesh on the host to a file for visualization.
   *
   * \param name The root filename to use.
   * \param hostMesh A Blueprint mesh.
   */
  void saveVisualization(const std::string &name, const conduit::Node &hostMesh)
  {
    if(m_visualize)
    {
      AXOM_ANNOTATE_SCOPE("saveVisualization");
#if defined(AXOM_USE_HDF5)
      conduit::relay::io::blueprint::save_mesh(hostMesh, name, "hdf5");
#endif
      conduit::relay::io::save(hostMesh, name + ".yaml", "yaml");
    }
  }

  /*!
   * \brief Perform a test or save a baseline, depending on the mode.
   *
   * \param name The name of the test (used for filenames)
   * \param currentMesh The current mesh to be compared to a baseline (or saved, if rebaselining).
   * \param tolerance The tolerance to use when comparing the Blueprint results.
   *
   * \return true on success; false if the test did not pass.
   */
  template <typename ExecSpace = axom::SEQ_EXEC>
  bool test(const std::string &name, const conduit::Node &currentMesh, double tolerance = 2.6e-6)
  {
    AXOM_ANNOTATE_SCOPE("test");
    bool retval = true;
    std::string baselineName(yamlRoot(name));
    const auto paths = baselinePaths<ExecSpace>();
    if(rebaseline(name))
    {
      SLIC_INFO(
        axom::fmt::format("Saving new baseline for {} {}.", name, execution_name<ExecSpace>::name()));
      saveBaseline(paths, baselineName, currentMesh);
    }
    else
    {
      conduit::Node info;
      retval = compareBaseline(paths, baselineName, currentMesh, info, tolerance);
      if(!retval && m_verbose)
      {
        info.print();
      }
    }
    return retval;
  }

  void setVerbose(bool value) { m_verbose = value; }
  bool getVerbose() const { return m_verbose; }

protected:
  /*!
   * \brief Initializes the class.
   */
  virtual void initialize()
  {
    // If an argument is provided, split it into a vector
    if(!m_rebaseline_raw.empty())
    {
      const auto testNames = axom::utilities::string::split(m_rebaseline_raw, ',');
      if(testNames.size() == 1 && testNames[0] == "all")
      {
        m_rebaseline.clear();
      }
      else if(!testNames.empty())
      {
        m_rebaseline = testNames;
      }
    }
  }

  /*!
   * \brief Return whether we should rebaseline a test.
   *
   * \param name The name of the test.
   *
   * \return true if the test should be rebaselined; false otherwise.
   */
  bool rebaseline(const std::string &name)
  {
    bool retval = false;
    if(m_rebaseline.size() == 1 && m_rebaseline[0] == "none")
    {
      // It does not look like the --rebaseline argument was given.
      retval = false;
    }
    else if(m_rebaseline.empty())
    {
      // --rebaseline all was given. Rebaseline all.
      retval = true;
    }
    else
    {
      // Rebaseline if name is in m_rebaseline.
      retval = std::find(m_rebaseline.begin(), m_rebaseline.end(), name) != m_rebaseline.end();
    }
    return retval;
  }

  /// Conduit error handler that blocks (helpful for getting a stack in a debugger)
  static void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
  {
    std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
    // This is on purpose.
    while(1);
  }

  axom::CLI::App m_app;
  std::string m_annotationMode;
  bool m_handler;
  bool m_visualize;
  bool m_verbose;
  std::string m_rebaseline_raw;
  std::vector<std::string> m_rebaseline;
};

}  // end namespace testing
}  // end namespace blueprint
}  // end namespace axom
#endif
