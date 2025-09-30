// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file core_array_perf.cpp
 *  \brief This example illustrates performance of array classes.
 */

// Axom includes
#include "axom/core/StackArray.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/MDMapping.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/nested_for_exec.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/utilities/RAII.hpp"
#include "axom/core/utilities/CommandLineUtilities.hpp"
#include "axom/core/utilities/System.hpp"
#include "axom/core/utilities/Timer.hpp"
#include "axom/core/AnnotationMacros.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

// C/C++ includes
#include <iostream>
#include <vector>

///////////////////////////////////////////////////////////////
/// Struct to parse and store the input parameters
struct InputParams
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

  // Array shape.
  std::vector<axom::IndexType> shape;
  axom::IndexType ghostWidth = 1;
  std::vector<axom::IndexType> paddedShape;
  std::vector<axom::IndexType> idxBegin;
  std::vector<axom::IndexType> idxEnd;
  axom::IndexType realSize;
  axom::IndexType paddedSize;

  // Array stride order (same length as shape)
  std::vector<std::uint32_t> dataSlowestDirections;
  axom::ArrayStrideOrder dataOrder = axom::ArrayStrideOrder::ARBITRARY;

  RuntimePolicy runtimePolicy = RuntimePolicy::seq;

  axom::IndexType repCount = 10;

  std::string annotationMode {"none"};

private:
  bool _verboseOutput {false};
  const std::map<std::string, int> strideValidator {{"row", int(axom::ArrayStrideOrder::ROW)},
                                                    {"col", int(axom::ArrayStrideOrder::COLUMN)}};

public:
  bool isVerbose() const { return _verboseOutput; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-p, --policy", runtimePolicy)
      ->description("Set runtime policy for test")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("-s, --shape", shape)->description("Array shape")->expected(1, 4);

    app.add_option("-g, --ghost", ghostWidth)->description("Ghost width");

    app.add_option("-r, --repCount", repCount)->description("Number of repetitions to run");

    auto* dataOrderOption = app.add_option("--dataOrder", dataOrder)
                              ->description("Stride order of array data.")
                              ->transform(axom::CLI::CheckedTransformer(strideValidator))
                              ->expected(1, 4);

    app.add_option("--dataSlowestDirections", dataSlowestDirections)
      ->description(
        "Array data stride directions, from slowest to fastest."
        "  Must be same length as shape.")
      ->excludes(dataOrderOption);

#ifdef AXOM_USE_CALIPER
    app.add_option("--caliper", annotationMode)
      ->description(
        "caliper annotation mode. Valid options include 'none' and 'report'. "
        "Use 'help' to see full list.")
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    app.get_formatter()->column_width(60);

    app.parse(argc, argv);

    if(shape.empty())
    {
      std::cerr << "You must specify shape (1-4 integers)." << std::endl;
      std::abort();
    }

    /*
      If dataSlowestDirections is specified, it must match length of shape.
      If neither is specified, default to row-major ordering.
    */
    if(!dataSlowestDirections.empty())
    {
      if(dataSlowestDirections.size() != shape.size())
      {
        std::cerr << axom::fmt::format("slowestDimension size ({}) must match shape size ({}).",
                                       dataSlowestDirections.size(),
                                       shape.size())
                  << std::endl;
        std::abort();
      }
    }
    if(dataSlowestDirections.empty() && dataOrder == axom::ArrayStrideOrder::ARBITRARY)
    {
      dataOrder = axom::ArrayStrideOrder::ROW;
    }

    //
    // Dependent data
    //

    paddedShape.resize(shape.size());
    idxBegin.resize(shape.size());
    idxEnd.resize(shape.size());
    for(size_t i = 0; i < shape.size(); ++i)
    {
      paddedShape[i] = shape[i] + 2 * ghostWidth;
      idxBegin[i] = ghostWidth;
      idxEnd[i] = idxBegin[i] + shape[i];
    }

    realSize = shape[0];
    paddedSize = paddedShape[0];
    for(size_t i = 1; i < shape.size(); ++i)
    {
      realSize *= shape[i];
      paddedSize *= paddedShape[i];
    }
  }
};

InputParams params;

//!@brief Return allocator id suitable for the given runtime policy.
int allocatorIdFromPolicy(axom::runtime_policy::Policy policy)
{
  AXOM_UNUSED_VAR(policy);
#if defined(AXOM_USE_UMPIRE)
  int allocatorID = policy == axom::runtime_policy::Policy::seq
    ? axom::detail::getAllocatorID<axom::MemorySpace::Host>()
    :
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    policy == axom::runtime_policy::Policy::omp
    ? axom::detail::getAllocatorID<axom::MemorySpace::Host>()
    :
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    policy == axom::runtime_policy::Policy::cuda
    ? axom::detail::getAllocatorID<axom::MemorySpace::Device>()
    :
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    policy == axom::runtime_policy::Policy::hip
    ? axom::detail::getAllocatorID<axom::MemorySpace::Device>()
    :
  #endif
    axom::INVALID_ALLOCATOR_ID;
#else
  int allocatorID = axom::getDefaultAllocatorID();
#endif
  return allocatorID;
}

template <int DIM, typename ExecSpace>
class MDMappingPerfTester
{
public:
  using Element_t = std::uint64_t;
  const Element_t m_baseFactor = 1000000;
  Element_t m_testAccumulation = 0;
  const Element_t m_flatTestAdd = 1;
  const Element_t m_rowTestAdd = 10;
  const Element_t m_columnTestAdd = 100;
  const Element_t m_dynamicTestAdd = 1000;

  int m_allocatorId;

  MDMappingPerfTester()
  {
    m_allocatorId = allocatorIdFromPolicy(params.runtimePolicy);
#ifdef AXOM_USE_UMPIRE
    umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
    umpire::Allocator allocator = rm.getAllocator(m_allocatorId);
    std::cout << axom::fmt::format("Allocator id: {}, Umpire memory space {}",
                                   m_allocatorId,
                                   allocator.getName())
              << std::endl;
#else
    std::cout << axom::fmt::format("Allocator id: {}, default memory space", m_allocatorId)
              << std::endl;
#endif
  }

  /*!
    @brief Time the pointer access of every element of an array.

    This is the fastest we expect to visit every element.
  */
  void runTest_pointerAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("pointerAccess");
    auto testAdd = m_flatTestAdd;
    m_testAccumulation += testAdd;
    auto count = array.size();
    auto* ptr = array.data();
    axom::for_all<ExecSpace>(0, count, AXOM_LAMBDA(axom::IndexType i) { ptr[i] += testAdd; });
  }

  /*!
    @brief Time the flat-index access of every element of an array.

    Compared to runtest_flatAccess, this includes the flatIndex overhead.
  */
  void runTest_flatAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("flatAccess");
    auto testAdd = m_flatTestAdd;
    m_testAccumulation += testAdd;
    auto count = array.size();
    axom::for_all<ExecSpace>(
      0,
      count,
      AXOM_LAMBDA(axom::IndexType i) { array.flatIndex(i) += testAdd; });
  }

  /*!
  Methods to time the access of every element of an array.

  Multidimensional loops are capable of skipping ghost
  layers, but the flat loop used for the baseline
  performance doesn't have logic to skip them.
*/

  //
  // Row-major access tests
  //

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 1>::type runTest_rowMajorAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("rowMajorAccess-1D");
    auto testAdd = m_rowTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    axom::for_all<ExecSpace>(
      idxBegin[0],
      idxEnd[0],
      AXOM_LAMBDA(axom::IndexType i) { array[i] += testAdd; });
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type runTest_rowMajorAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("rowMajorAccess-2D");
    auto testAdd = m_rowTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    axom::StackArray<axom::IndexType, 2> jRange {{idxBegin[1], idxEnd[1]}};
    axom::StackArray<axom::IndexType, 2> iRange {{idxBegin[0], idxEnd[0]}};
    axom::for_all<ExecSpace>(
      jRange,
      iRange,
      AXOM_LAMBDA(axom::IndexType j, axom::IndexType i) { array(i, j) += testAdd; });
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type runTest_rowMajorAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("rowMajorAccess-3D");
    auto testAdd = m_rowTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    axom::StackArray<axom::IndexType, 2> kRange {{idxBegin[2], idxEnd[2]}};
    axom::StackArray<axom::IndexType, 2> jRange {{idxBegin[1], idxEnd[1]}};
    axom::StackArray<axom::IndexType, 2> iRange {{idxBegin[0], idxEnd[0]}};
    axom::for_all<ExecSpace>(
      kRange,
      jRange,
      iRange,
      AXOM_LAMBDA(axom::IndexType k, axom::IndexType j, axom::IndexType i) {
        array(i, j, k) += testAdd;
      });
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 4>::type runTest_rowMajorAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("rowMajorAccess-4D");
    auto testAdd = m_rowTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
    AXOM_UNUSED_VAR(array);
    AXOM_UNUSED_VAR(idxBegin);
    AXOM_UNUSED_VAR(idxEnd);
    std::cerr << "Cannot run higher than 3D with RAJA." << std::endl;
    std::abort();
#else
    for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
    {
      for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
      {
        for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
        {
          for(axom::IndexType l = idxBegin[3]; l < idxEnd[3]; ++l)
          {
            array(i, j, k, l) += testAdd;
          }
        }
      }
    }
#endif
  }

  //
  // Colunn-major access tests
  //

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 1>::type runTest_columnMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("columnMajorAccess-1D");
    auto testAdd = m_columnTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    axom::for_all<ExecSpace>(
      idxBegin[0],
      idxEnd[0],
      AXOM_LAMBDA(axom::IndexType i) { array[i] += testAdd; });
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type runTest_columnMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("columnMajorAccess-2D");
    auto testAdd = m_columnTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    axom::StackArray<axom::IndexType, 2> jRange {{idxBegin[1], idxEnd[1]}};
    axom::StackArray<axom::IndexType, 2> iRange {{idxBegin[0], idxEnd[0]}};
    axom::for_all<ExecSpace>(
      iRange,
      jRange,
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j) { array(i, j) += testAdd; });
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type runTest_columnMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("columnMajorAccess-3D");
    auto testAdd = m_columnTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    axom::StackArray<axom::IndexType, 2> kRange {{idxBegin[2], idxEnd[2]}};
    axom::StackArray<axom::IndexType, 2> jRange {{idxBegin[1], idxEnd[1]}};
    axom::StackArray<axom::IndexType, 2> iRange {{idxBegin[0], idxEnd[0]}};
    axom::for_all<ExecSpace>(
      iRange,
      jRange,
      kRange,
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k) {
        array(i, j, k) += testAdd;
      });
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 4>::type runTest_columnMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("columnMajorAccess-4D");
    auto testAdd = m_columnTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
    AXOM_UNUSED_VAR(array);
    AXOM_UNUSED_VAR(idxBegin);
    AXOM_UNUSED_VAR(idxEnd);
    std::cerr << "Cannot run higher than 3D with RAJA." << std::endl;
    std::abort();
#else
    for(axom::IndexType l = idxBegin[3]; l < idxEnd[3]; ++l)
    {
      for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
      {
        for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
        {
          for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
          {
            array(i, j, k, l) += testAdd;
          }
        }
      }
    }
#endif
  }

  //
  // Dynamic ordering access tests

  // The dynamic order should match the optimal order,
  // Any performance difference is due to overhead of dynamic
  // nesting of the loops.
  //

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 1>::type runTest_dynamicAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("dynamicAccess-1D");
    auto testAdd = m_dynamicTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    axom::for_all<ExecSpace>(
      idxBegin[0],
      idxEnd[0],
      AXOM_LAMBDA(axom::IndexType i) { array[i] += testAdd; });
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type runTest_dynamicAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("dynamicAccess-2D");
    auto testAdd = m_dynamicTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    const auto& mapping = array.mapping();
    const auto& slowestDirs = mapping.slowestDirs();
    axom::StackArray<axom::IndexType, DIM> invSlowestDirs;
    for(int i = 0; i < DIM; ++i)
    {
      invSlowestDirs[slowestDirs[i]] = i;
    }
    const axom::StackArray<axom::IndexType, DIM> begins {idxBegin[slowestDirs[0]],
                                                         idxBegin[slowestDirs[1]]};
    const axom::StackArray<axom::IndexType, DIM> ends {idxEnd[slowestDirs[0]],
                                                       idxEnd[slowestDirs[1]]};
    axom::StackArray<axom::IndexType, 2> nRange {{begins[1], ends[1]}};
    axom::StackArray<axom::IndexType, 2> mRange {{begins[0], ends[0]}};
    axom::for_all<ExecSpace>(
      nRange,
      mRange,
      AXOM_LAMBDA(axom::IndexType n, axom::IndexType m) {
        axom::StackArray<axom::IndexType, DIM> idx {m, n};
        auto i = idx[invSlowestDirs[0]];
        auto j = idx[invSlowestDirs[1]];
        array(i, j) += testAdd;
      });
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type runTest_dynamicAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("dynamicAccess-3D");
    auto testAdd = m_dynamicTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    const auto& mapping = array.mapping();
    const auto& slowestDirs = mapping.slowestDirs();
    axom::StackArray<axom::IndexType, DIM> invSlowestDirs;
    for(int i = 0; i < DIM; ++i)
    {
      invSlowestDirs[slowestDirs[i]] = i;
    }
    const axom::StackArray<axom::IndexType, DIM> begins {idxBegin[slowestDirs[0]],
                                                         idxBegin[slowestDirs[1]],
                                                         idxBegin[slowestDirs[2]]};
    const axom::StackArray<axom::IndexType, DIM> ends {idxEnd[slowestDirs[0]],
                                                       idxEnd[slowestDirs[1]],
                                                       idxEnd[slowestDirs[2]]};
    axom::StackArray<axom::IndexType, 2> oRange {{begins[2], ends[2]}};
    axom::StackArray<axom::IndexType, 2> nRange {{begins[1], ends[1]}};
    axom::StackArray<axom::IndexType, 2> mRange {{begins[0], ends[0]}};

    axom::for_all<ExecSpace>(
      oRange,
      nRange,
      mRange,
      AXOM_LAMBDA(axom::IndexType o, axom::IndexType n, axom::IndexType m) {
        axom::StackArray<axom::IndexType, DIM> idx {m, n, o};
        auto i = idx[invSlowestDirs[0]];
        auto j = idx[invSlowestDirs[1]];
        auto k = idx[invSlowestDirs[2]];
        array(i, j, k) += testAdd;
      });
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 4>::type runTest_dynamicAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_ANNOTATE_SCOPE("dynamicAccess-4D");
    auto testAdd = m_dynamicTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    const auto& mapping = array.mapping();
    const auto& slowestDirs = mapping.slowestDirs();
    const axom::StackArray<axom::IndexType, DIM> begins {idxBegin[slowestDirs[0]],
                                                         idxBegin[slowestDirs[1]],
                                                         idxBegin[slowestDirs[2]],
                                                         idxBegin[slowestDirs[3]]};
    const axom::StackArray<axom::IndexType, DIM> ends {idxEnd[slowestDirs[0]],
                                                       idxEnd[slowestDirs[1]],
                                                       idxEnd[slowestDirs[2]],
                                                       idxEnd[slowestDirs[3]]};
#ifdef AXOM_USE_RAJA
    AXOM_UNUSED_VAR(array);
    AXOM_UNUSED_VAR(begins);
    AXOM_UNUSED_VAR(ends);
    std::cerr << "Cannot run higher than 3D with RAJA." << std::endl;
    std::abort();
#else
    axom::StackArray<axom::IndexType, DIM> idx;
    axom::IndexType& m = idx[slowestDirs[0]];
    axom::IndexType& n = idx[slowestDirs[1]];
    axom::IndexType& o = idx[slowestDirs[2]];
    axom::IndexType& p = idx[slowestDirs[3]];
    for(m = begins[0]; m < ends[0]; ++m)
    {
      for(n = begins[1]; n < ends[1]; ++n)
      {
        for(o = begins[2]; o < ends[2]; ++o)
        {
          for(p = begins[3]; p < ends[3]; ++p)
          {
            array[idx] += testAdd;
          }
        }
      }
    }
#endif
  }

  //
  // Make test Array objects.
  //

  /*!
    @brief Make MDMapping based on command line parameters.
  */
  void getPaddedShapeAndIndexer(axom::StackArray<axom::IndexType, DIM>& paddedShape,
                                axom::MDMapping<DIM>& mapping)
  {
    for(int d = 0; d < DIM; ++d)
    {
      paddedShape[d] = params.paddedShape[d];
    }

    if(!params.dataSlowestDirections.empty())
    {
      axom::StackArray<std::uint16_t, DIM> dataSlowestDirections;
      for(int d = 0; d < DIM; ++d)
      {
        dataSlowestDirections[d] = params.dataSlowestDirections[d];
      }
      mapping.initializeShape(paddedShape, dataSlowestDirections);
    }
    else
    {
      mapping.initializeShape(paddedShape, params.dataOrder);
    }
  }

  /*!
    @brief Return an array for testing, dimension DIM,
    sized and ordered according to params values.
  */
  axom::Array<Element_t, DIM> makeArray()
  {
    assert(DIM <= params.shape.size());

    AXOM_ANNOTATE_SCOPE("makeArray");
    axom::StackArray<axom::IndexType, DIM> paddedShape;
    axom::MDMapping<DIM> mapping;
    getPaddedShapeAndIndexer(paddedShape, mapping);

    return axom::Array<Element_t, DIM>(paddedShape, mapping.slowestDirs(), m_allocatorId);
  }

  /*!
    @brief Run test with array shape of the first DIM values in
    params.shape.
  */
  void runTest_dim()
  {
    AXOM_ANNOTATE_SCOPE("MDMappingPerfTester::runTest_dim");
    axom::Array<Element_t, DIM> arrayMd = makeArray();
    axom::ArrayView<Element_t, DIM> array = arrayMd.view();

    std::cout << axom::fmt::format("Real-to-padded size: {}/{} = {:.4f}",
                                   params.realSize,
                                   params.paddedSize,
                                   double(params.realSize) / params.paddedSize)
              << std::endl;

    auto count = array.size();
    auto baseFactor = m_baseFactor;
    axom::for_all<ExecSpace>(
      count,
      AXOM_LAMBDA(axom::IndexType i) { array.flatIndex(i) = Element_t(i * baseFactor); });

    /*
      Warm-up: The first rep can be a slow outlier.  Since we are
      generally interested in the steady-state performance, get past
      this outlier before timing.
    */
    for(int r = 0; r < 1; ++r)
    {
      runTest_flatAccess(array);
      runTest_rowMajorAccess(array);
      runTest_columnMajorAccess(array);
      runTest_dynamicAccess(array);
    }

    axom::utilities::Timer flatTimer(false);
    flatTimer.start();
    for(int r = 0; r < params.repCount; ++r)
    {
      runTest_flatAccess(array);
    }
    flatTimer.stop();
    std::cout << axom::fmt::format("Avg flat-index time   {:.8f} seconds, base",
                                   flatTimer.elapsedTimeInSec() / params.repCount)
              << std::endl;

    auto baseTime = flatTimer.elapsedTimeInSec();

    axom::utilities::Timer pointerTimer(false);
    pointerTimer.start();
    for(int r = 0; r < params.repCount; ++r)
    {
      runTest_pointerAccess(array);
    }
    pointerTimer.stop();
    std::cout << axom::fmt::format("Avg pointer time      {:.8f} seconds, {:.2f}x",
                                   pointerTimer.elapsedTimeInSec() / params.repCount,
                                   pointerTimer.elapsedTimeInSec() / baseTime)
              << std::endl;

    axom::utilities::Timer rowMajorTimer(false);
    rowMajorTimer.start();
    for(int r = 0; r < params.repCount; ++r)
    {
      runTest_rowMajorAccess(array);
    }
    rowMajorTimer.stop();
    std::cout << axom::fmt::format("Avg row-major time    {:.8f} seconds, {:.2f}x",
                                   rowMajorTimer.elapsedTimeInSec() / params.repCount,
                                   rowMajorTimer.elapsedTimeInSec() / baseTime)
              << std::endl;

    axom::utilities::Timer columnMajorTimer(false);
    columnMajorTimer.start();
    for(int r = 0; r < params.repCount; ++r)
    {
      runTest_columnMajorAccess(array);
    }
    columnMajorTimer.stop();
    std::cout << axom::fmt::format("Avg column-major time {:.8f} seconds, {:.2f}x",
                                   columnMajorTimer.elapsedTimeInSec() / params.repCount,
                                   columnMajorTimer.elapsedTimeInSec() / baseTime)
              << std::endl;

    axom::utilities::Timer dynamicTimer(false);
    dynamicTimer.start();
    for(int r = 0; r < params.repCount; ++r)
    {
      runTest_dynamicAccess(array);
    }
    dynamicTimer.stop();
    std::cout << axom::fmt::format("Avg dynamic time      {:.8f} seconds, {:.2f}x",
                                   dynamicTimer.elapsedTimeInSec() / params.repCount,
                                   dynamicTimer.elapsedTimeInSec() / baseTime)
              << std::endl;

    {
      AXOM_ANNOTATE_SCOPE("MDMappingPerfTester::runTest_dim-verify");
      // Verify that the elements are touched the correct number of times.
      // (Bring the data to the host so this test doesn't rely on RAJA.)
      auto hostAllocatorId =
        axom::detail::getAllocatorID<axom::execution_space<axom::SEQ_EXEC>::memory_space>();
      axom::Array<Element_t, DIM> hostArray(array, hostAllocatorId);
      axom::IndexType matchCount = 0;
      for(axom::IndexType i = 0; i < count; ++i)
      {
        matchCount += (hostArray.flatIndex(i) == Element_t(i * m_baseFactor) + m_testAccumulation);
      }
      if(matchCount != params.realSize)
      {
        std::cerr << axom::fmt::format(
                       "Unexpected error in tests: counted match ({}) != expected ({})",
                       matchCount,
                       params.realSize)
                  << std::endl;
      }
    }
  }

};  // class MDMappingPerfTester

/*!
  @brief Run test based on dimension specified in params.

  On devices, we skip DIM > 3.  We don't have the RAJA set-up for
  higher dimensions.
*/
template <typename ExecSpace>
void runTest()
{
  AXOM_ANNOTATE_SCOPE("runTest");
  if(params.shape.size() == 1)
  {
    MDMappingPerfTester<1, ExecSpace> tester1D;
    tester1D.runTest_dim();
  }
  else if(params.shape.size() == 2)
  {
    MDMappingPerfTester<2, ExecSpace> tester2D;
    tester2D.runTest_dim();
  }
  else if(params.shape.size() == 3)
  {
    MDMappingPerfTester<3, ExecSpace> tester3D;
    tester3D.runTest_dim();
  }
  else if(params.shape.size() == 4 && !axom::execution_space<ExecSpace>::onDevice())
  {
    MDMappingPerfTester<4, ExecSpace> tester4D;
    tester4D.runTest_dim();
  }
}

int main(int argc, char** argv)
{
  axom::CLI::App app {"Driver for array indexing performace tests"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    auto retval = app.exit(e);
    exit(retval);
  }

  axom::utilities::raii::AnnotationsWrapper annotation_raii_wrapper(params.annotationMode);
  AXOM_ANNOTATE_SCOPE("core array performance example");

  auto hostname = axom::utilities::getHostName();
  std::cout << axom::fmt::format("Host: {}", hostname) << std::endl;

  using RuntimePolicy = axom::runtime_policy::Policy;

  std::cout << axom::fmt::format("Runtime policy: {}",
                                 axom::runtime_policy::policyToName(params.runtimePolicy))
            << std::endl;
  std::cout << axom::fmt::format("Array shape: {::d}", params.shape) << std::endl;
  std::cout << axom::fmt::format("Data order: {}",
                                 (params.dataOrder == axom::ArrayStrideOrder::ROW ? "row" : "col"))
            << std::endl;
  std::cout << axom::fmt::format("Repetition count: {}", params.repCount) << std::endl;

  if(params.runtimePolicy == RuntimePolicy::seq)
  {
    runTest<axom::SEQ_EXEC>();
  }
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  else if(params.runtimePolicy == RuntimePolicy::omp)
  {
    runTest<axom::OMP_EXEC>();
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  else if(params.runtimePolicy == RuntimePolicy::cuda)
  {
    runTest<axom::CUDA_EXEC<256>>();
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(params.runtimePolicy == RuntimePolicy::hip)
  {
    runTest<axom::HIP_EXEC<256>>();
  }
#endif
  return 0;
}
