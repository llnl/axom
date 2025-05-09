// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 \file quest_mesh_view_util.cpp
 \brief Test for MeshViewUtil class.
*/

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT

  // Axom includes
  #include "axom/core.hpp"
  #include "axom/slic.hpp"
  #include "axom/primal.hpp"
  #include "axom/quest/MeshViewUtil.hpp"
  #include "axom/core/Types.hpp"

  #include "conduit_blueprint.hpp"
  #include "conduit_relay_io_blueprint.hpp"

  #include "axom/fmt.hpp"
  #include "axom/CLI11.hpp"

  // C/C++ includes
  #include <string>

namespace quest = axom::quest;
namespace slic = axom::slic;
namespace primal = axom::primal;
namespace numerics = axom::numerics;

///////////////////////////////////////////////////////////////
// converts the input string into an 80 character string
// padded on both sides with '=' symbols
std::string banner(const std::string& str) { return axom::fmt::format("{:=^80}", str); }

///////////////////////////////////////////////////////////////
/// Struct to parse and store the input parameters
struct Input
{
private:
  bool _verboseOutput {false};

public:
  bool isVerbose() const { return _verboseOutput; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.get_formatter()->column_width(60);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(_verboseOutput ? slic::message::Debug : slic::message::Info);
  }
};

Input params;

/// Utility function to initialize the logger
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  slic::LogStream* logStream;

  std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);

  slic::addStreamToAllMsgLevels(logStream);

  conduit::utils::set_error_handler(
    [](auto& msg, auto& file, int line) { slic::logErrorMessage(msg, file, line); });
  conduit::utils::set_warning_handler(
    [](auto& msg, auto& file, int line) { slic::logWarningMessage(msg, file, line); });
  conduit::utils::set_info_handler([](auto& msg, auto& file, int line) {
    slic::logMessage(slic::message::Info, msg, file, line);
  });
}

/// Utility function to finalize the logger
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

//! @brief Add a scalar to all elements of a StackArray.
template <typename T, int DIM, typename U>
axom::StackArray<T, DIM> operator+(const axom::StackArray<T, DIM>& rhs, U lhs)
{
  axom::StackArray<T, DIM> ret = rhs;
  for(int d = 0; d < DIM; ++d)
  {
    ret[d] += lhs;
  }
  return ret;
}

//! @brief Add two StackArrays of the same dimension.
template <typename T, int DIM, typename U>
axom::StackArray<T, DIM> operator+(const axom::StackArray<T, DIM>& rhs,
                                   const axom::StackArray<U, DIM>& lhs)
{
  axom::StackArray<T, DIM> ret = rhs;
  for(int d = 0; d < DIM; ++d)
  {
    ret[d] += lhs[d];
  }
  return ret;
}

template <typename T, int DIM>
T product(const axom::StackArray<T, DIM>& a)
{
  T ret = 1;
  for(int d = 0; d < DIM; ++d)
  {
    ret *= a[d];
  }
  return ret;
}

template <int DIM, typename A, typename B>
bool isEqual(const axom::StackArray<A, DIM>& a, const axom::StackArray<B, DIM>& b)
{
  for(int d = 0; d < DIM; ++d)
  {
    if(a[d] != b[d])
    {
      return false;
    }
  }
  return true;
}

using MVU2Type = axom::quest::MeshViewUtil<2, axom::MemorySpace::Dynamic>;
using MVU3Type = axom::quest::MeshViewUtil<3, axom::MemorySpace::Dynamic>;
using IndexCoords = axom::StackArray<axom::IndexType, 3>;

//---------------------------------------------------------------------------
// Create a mesh using a Conduit example.
// ---------------------------------------------------------------------------
void createConduitBlueprintDomain(const axom::StackArray<std::int64_t, 3>& domainShape,
                                  const axom::StackArray<std::int64_t, 3>& loPads,
                                  const axom::StackArray<std::int64_t, 3>& hiPads,
                                  conduit::Node& domain)
{
  // Conduit's example requires int64_t inputs.
  using Int64Coords = axom::StackArray<std::int64_t, 3>;
  Int64Coords elemOrigin = loPads;
  Int64Coords elemPaddedShape = domainShape + loPads + hiPads;
  Int64Coords vertOrigin = loPads;
  Int64Coords vertPaddedShape = elemPaddedShape + 1;
  Int64Coords vertShape = domainShape + 1;
  conduit::Node desc;
  desc["element_data/shape"].set(conduit::DataType::int64(3), elemPaddedShape);
  desc["element_data/origin"].set(conduit::DataType::int64(3), elemOrigin);
  desc["vertex_data/shape"].set(conduit::DataType::int64(3), vertPaddedShape);
  desc["vertex_data/origin"].set(conduit::DataType::int64(3), vertOrigin);

  conduit::blueprint::mesh::examples::strided_structured(desc,
                                                         vertShape[0],
                                                         vertShape[1],
                                                         vertShape[2],
                                                         domain);

  return;
}

template <typename T>
struct ModifyNumber
{
  void operator()(T& t) const { t = 1000000 + 2 * t; }
};

// Check that all elements of an ArrayView can be accessed
// and found at the expected distance from the first element.
template <typename T, typename FUNCTOR>
int testDataAccess(axom::ArrayView<T, 3>& da,
                   const std::string& daName,
                   const IndexCoords& strideOrder,
                   const FUNCTOR& functor)
{
  const auto& shape = da.shape();

  const IndexCoords zero {0, 0, 0};
  IndexCoords ijk;
  axom::IndexType& l = ijk[strideOrder[2]];
  axom::IndexType& m = ijk[strideOrder[1]];
  axom::IndexType& n = ijk[strideOrder[0]];

  const axom::IndexType lMax = shape[strideOrder[2]];
  const axom::IndexType mMax = shape[strideOrder[1]];
  const axom::IndexType nMax = shape[strideOrder[0]];

  int errCount = 0;

  // Check ArrayView access.
  for(l = 0; l < lMax; ++l)
  {
    for(m = 0; m < mMax; ++m)
    {
      for(n = 0; n < nMax; ++n)
      {
        functor(da[ijk]);
      }
    }
  }

  // Check expected monotonicy.
  bool monotone = true;
  axom::IndexType prevDist = -1;
  for(l = 0; l < lMax; ++l)
  {
    for(m = 0; m < mMax; ++m)
    {
      for(n = 0; n < nMax; ++n)
      {
        auto dist = &da[ijk] - &da[zero];
        // We loop in striding order, so dist should always increase.
        // Can't require strict sequential order because we skip ghosts.
        monotone &= (prevDist < dist);
        prevDist = dist;
      }
    }
  }
  if(!monotone)
  {
    ++errCount;
    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("Failed monotone ordering for ArrayView '{}'", daName));
  }

  return errCount;
}

//---------------------------------------------------------------------------
// Test methods that convert between shapes and strides-and-offsets.
//---------------------------------------------------------------------------
int testConversionMethods(const MVU3Type::MdIndices& realShape,
                          const MVU3Type::MdIndices& loPads,
                          const MVU3Type::MdIndices& hiPads,
                          const MVU3Type::MdIndices& strideOrder,
                          axom::IndexType minStride = 1)
{
  SLIC_INFO_IF(params.isVerbose(),
               axom::fmt::format("Conversion test for realShape={}, loPads={}, "
                                 "hiPads={}, strideOrder={}, minStride={}",
                                 realShape,
                                 loPads,
                                 hiPads,
                                 strideOrder,
                                 minStride));
  int errCount = 0;

  MVU3Type::MdIndices paddedShape = realShape + loPads + hiPads;

  MVU3Type::MdIndices offsets;
  MVU3Type::MdIndices strides;
  axom::IndexType valuesCount = 0;
  axom::quest::internal::shapesToStridesAndOffsets(realShape,
                                                   loPads,
                                                   hiPads,
                                                   strideOrder,
                                                   minStride,
                                                   offsets,
                                                   strides,
                                                   valuesCount);

  MVU3Type::MdIndices loPads1;
  MVU3Type::MdIndices hiPads1;
  MVU3Type::MdIndices paddedShape1;
  MVU3Type::MdIndices strideOrder1;
  axom::quest::internal::stridesAndOffsetsToShapes(realShape,
                                                   offsets,
                                                   strides,
                                                   valuesCount,
                                                   paddedShape1,
                                                   loPads1,
                                                   hiPads1,
                                                   strideOrder1);

  // Conversions should return original values.
  if(!isEqual(loPads, loPads1))
  {
    ++errCount;
    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("Mismatched loPads: {} vs {}.", loPads1, loPads));
  }
  if(!isEqual(hiPads, hiPads1))
  {
    ++errCount;
    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("Mismatched hiPads: {} vs {}.", hiPads1, hiPads));
  }
  if(paddedShape1 != paddedShape)
  {
    ++errCount;
    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("Mismatched valuesCount: {} vs {}.", paddedShape1, paddedShape));
  }
  if(!isEqual(strideOrder, strideOrder1))
  {
    ++errCount;
    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("Mismatched strideOrder: {} vs {}.", strideOrder1, strideOrder));
  }

  return errCount;
}

//---------------------------------------------------------------------------
// Create a mesh using a Conduit example, take its view and verify
// expectations.
// ---------------------------------------------------------------------------
int testByConduitExample(const IndexCoords& domainShape,
                         const IndexCoords& loPads,
                         const IndexCoords& hiPads,
                         const IndexCoords& strideOrder,
                         axom::IndexType minStride = 1)
{
  SLIC_INFO_IF(params.isVerbose(),
               axom::fmt::format("Check conduit mesh with domainShape={}, loPads={}, "
                                 "hiPads={}, strideOrder={}, minStride={}",
                                 domainShape,
                                 loPads,
                                 hiPads,
                                 strideOrder,
                                 minStride));
  int errCount = 0;

  conduit::Node domain;
  createConduitBlueprintDomain({domainShape[0], domainShape[1], domainShape[2]},
                               {loPads[0], loPads[1], loPads[2]},
                               {hiPads[0], hiPads[1], hiPads[2]},
                               domain);
  if(params.isVerbose())
  {
    SLIC_INFO("Testing with this domain:");
    domain.print();
  }

  axom::quest::MeshViewUtil<3, axom::MemorySpace::Dynamic> mview(domain);

  // Data parameters not dependent on striding order.

  IndexCoords elemShape = mview.getCellShape();
  IndexCoords vertShape = mview.getNodeShape();

  IndexCoords elemPaddedShape = domainShape + loPads + hiPads;
  IndexCoords vertPaddedShape = elemPaddedShape + 1;

  IndexCoords vertOrigin = loPads;

  // Access calls shouldn't abort.
  mview.getTopology();
  mview.getCoordSet();

  {
    // Verify data created by the Conduit example,
    // by comparing to results from separate calls to converter.
    // This data always uses column-major order,
    // regardless of the specified strideOrder.
    IndexCoords conduitStrideOrder = {0, 1, 2};

    IndexCoords elemOffsets;
    IndexCoords elemStrides;
    axom::IndexType elemValuesCount;
    axom::quest::internal::shapesToStridesAndOffsets(elemShape,
                                                     loPads,
                                                     hiPads,
                                                     conduitStrideOrder,
                                                     minStride,
                                                     elemOffsets,
                                                     elemStrides,
                                                     elemValuesCount);

    IndexCoords vertOffsets;
    IndexCoords vertStrides;
    axom::IndexType vertValuesCount;
    axom::quest::internal::shapesToStridesAndOffsets(vertShape,
                                                     loPads,
                                                     hiPads,
                                                     conduitStrideOrder,
                                                     minStride,
                                                     vertOffsets,
                                                     vertStrides,
                                                     vertValuesCount);

    if(!isEqual(mview.getCellShape(), elemShape))
    {
      ++errCount;
      SLIC_INFO_IF(
        params.isVerbose(),
        axom::fmt::format("Mismatched domain shape: {} vs {}", mview.getCellShape(), elemShape));
    }

    if(!isEqual(mview.getRealShape("element"), elemShape))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched real element shape: {} vs {}",
                                     mview.getRealShape("element"),
                                     elemShape));
    }

    if(!isEqual(mview.getRealShape("vertex"), vertShape))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched real vertex shape: {} vs {}",
                                     mview.getRealShape("vertex"),
                                     vertShape));
    }

    if(!isEqual(mview.getCoordsOffsets(), vertOrigin))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched coords offsets: {} vs {}",
                                     mview.getCoordsOffsets(),
                                     vertOrigin));
    }

    if(!isEqual(mview.getCoordsStrides(), vertStrides))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched coords strides: {} vs {}",
                                     mview.getCoordsStrides(),
                                     vertStrides));
    }

    if(mview.getCoordsCountWithGhosts() != vertValuesCount)
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched coords count with ghosts: {} vs {}",
                                     mview.getCoordsCountWithGhosts(),
                                     vertValuesCount));
    }

    // Check field views sizes and strides.
    // These may change if Conduit example changes.
    auto vertField = mview.getFieldView<double>("vert_vals", false);
    auto elemField = mview.getFieldView<double>("ele_vals", false);
    auto vertFieldWithGhosts = mview.getFieldView<double>("vert_vals", true);
    auto elemFieldWithGhosts = mview.getFieldView<double>("ele_vals", true);

    if(!isEqual(elemField.shape(), elemShape))
    {
      ++errCount;
      SLIC_INFO_IF(
        params.isVerbose(),
        axom::fmt::format("Mismatched elemField shape: {} vs {}", elemField.shape(), elemShape));
    }

    if(!isEqual(vertField.shape(), vertShape))
    {
      ++errCount;
      SLIC_INFO_IF(
        params.isVerbose(),
        axom::fmt::format("Mismatched vertField shape: {} vs {}", vertField.shape(), vertShape));
    }

    if(!isEqual(elemFieldWithGhosts.shape(), elemPaddedShape))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched elemField padded shape: {} vs {}",
                                     elemField.shape(),
                                     elemPaddedShape));
    }

    if(!isEqual(elemFieldWithGhosts.strides(), elemStrides))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched elemField padded stride: {} vs {}",
                                     elemFieldWithGhosts.strides(),
                                     elemStrides));
    }

    if(!isEqual(vertFieldWithGhosts.shape(), vertPaddedShape))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched vertField padded shape: {} vs {}",
                                     vertField.shape(),
                                     vertPaddedShape));
    }

    if(!isEqual(vertFieldWithGhosts.strides(), vertStrides))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched vertField padded stride: {} vs {}",
                                     vertField.strides(),
                                     vertStrides));
    }

    errCount += testDataAccess(elemField, "elemField", conduitStrideOrder, ModifyNumber<double>());
    errCount += testDataAccess(vertField, "vertField", conduitStrideOrder, ModifyNumber<double>());
    errCount += testDataAccess(elemFieldWithGhosts,
                               "elemFieldWithGhosts",
                               conduitStrideOrder,
                               ModifyNumber<double>());
    errCount += testDataAccess(vertFieldWithGhosts,
                               "vertFieldWithGhosts",
                               conduitStrideOrder,
                               ModifyNumber<double>());

    auto coordsViews = mview.getCoordsViews();
    for(int d = 0; d < 3; ++d)
    {
      errCount += testDataAccess(coordsViews[d],
                                 axom::fmt::format("coordsViews[{}]", d),
                                 conduitStrideOrder,
                                 ModifyNumber<double>());
    }
  }

  {
    // Verify data created by by us,
    // by comparing to results from separate calls to converter.
    // Unlike the data created by the conduit example, this
    // data always uses the specified striding order strideOrder.

    IndexCoords elemOffsets;
    IndexCoords elemStrides;
    axom::IndexType elemValuesCount;
    axom::quest::internal::shapesToStridesAndOffsets(elemShape,
                                                     loPads,
                                                     hiPads,
                                                     strideOrder,
                                                     minStride,
                                                     elemOffsets,
                                                     elemStrides,
                                                     elemValuesCount);

    IndexCoords vertOffsets;
    IndexCoords vertStrides;
    axom::IndexType vertValuesCount;
    axom::quest::internal::shapesToStridesAndOffsets(vertShape,
                                                     loPads,
                                                     hiPads,
                                                     strideOrder,
                                                     minStride,
                                                     vertOffsets,
                                                     vertStrides,
                                                     vertValuesCount);

    //---------------------------------------------------------------------------
    // Create a cell-centered field and verify shape and strides.
    //---------------------------------------------------------------------------
    mview.createField("inCells",
                      "element",
                      conduit::DataType::int32(elemValuesCount),
                      elemStrides,
                      elemOffsets);
    auto inCellsView = mview.getFieldView<std::int32_t>("inCells", false);
    auto inCellsPaddedView = mview.getFieldView<std::int32_t>("inCells", true);

    if(!isEqual(inCellsView.shape(), elemShape))
    {
      ++errCount;
      SLIC_INFO_IF(
        params.isVerbose(),
        axom::fmt::format("Mismatched vertField shape: {} vs {}", inCellsView.shape(), elemShape));
    }

    if(!isEqual(inCellsView.strides(), elemStrides))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched vertField strides: {} vs {}",
                                     inCellsView.strides(),
                                     elemStrides));
    }

    if(!isEqual(inCellsPaddedView.shape(), elemPaddedShape))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched vertField padded shape: {} vs {}",
                                     inCellsPaddedView.shape(),
                                     elemPaddedShape));
    }

    if(!isEqual(inCellsPaddedView.strides(), elemStrides))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched vertField padded strides: {} vs {}",
                                     inCellsPaddedView.strides(),
                                     elemStrides));
    }

    //---------------------------------------------------------------------------
    // Create a node-centered field and verify shape and strides.
    //---------------------------------------------------------------------------
    mview.createField("onNodes",
                      "vertex",
                      conduit::DataType::int32(vertValuesCount),
                      vertStrides,
                      vertOffsets);
    auto onNodesView = mview.getFieldView<std::int32_t>("onNodes", false);
    auto onNodesPaddedView = mview.getFieldView<std::int32_t>("onNodes", true);

    if(!isEqual(onNodesView.shape(), vertShape))
    {
      ++errCount;
      SLIC_INFO_IF(
        params.isVerbose(),
        axom::fmt::format("Mismatched vertField shape: {} vs {}", onNodesView.shape(), vertShape));
    }

    if(!isEqual(onNodesView.strides(), vertStrides))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched vertField strides: {} vs {}",
                                     onNodesView.strides(),
                                     vertStrides));
    }

    if(!isEqual(onNodesPaddedView.shape(), vertPaddedShape))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched vertField padded shape: {} vs {}",
                                     onNodesPaddedView.shape(),
                                     vertPaddedShape));
    }

    if(!isEqual(onNodesPaddedView.strides(), vertStrides))
    {
      ++errCount;
      SLIC_INFO_IF(params.isVerbose(),
                   axom::fmt::format("Mismatched vertField padded strides: {} vs {}",
                                     onNodesPaddedView.strides(),
                                     vertStrides));
    }

    errCount +=
      testDataAccess(inCellsPaddedView, "inCellsPaddedView", strideOrder, ModifyNumber<std::int32_t>());
    errCount += testDataAccess(inCellsView, "inCellsView", strideOrder, ModifyNumber<std::int32_t>());
    errCount +=
      testDataAccess(onNodesPaddedView, "onNodesPaddedView", strideOrder, ModifyNumber<std::int32_t>());
    errCount += testDataAccess(onNodesView, "onNodesView", strideOrder, ModifyNumber<std::int32_t>());
  }

  return errCount;
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  int errCount = 0;

  initializeLogger();

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  axom::CLI::App app {"Test code for MeshViewUtil class"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = app.exit(e);
    exit(retval);
  }

  {
    // Test conversion methods.
    errCount += testConversionMethods({2, 3, 5}, {0, 0, 0}, {0, 0, 0}, {0, 1, 2});
    errCount += testConversionMethods({2, 3, 5}, {0, 0, 0}, {0, 0, 0}, {2, 1, 0});
    errCount += testConversionMethods({2, 3, 5}, {1, 1, 1}, {1, 1, 1}, {0, 1, 2});
    errCount += testConversionMethods({2, 3, 5}, {1, 1, 1}, {1, 1, 1}, {2, 1, 0});
    errCount += testConversionMethods({2, 3, 5}, {2, 2, 2}, {1, 1, 1}, {0, 1, 2});
    errCount += testConversionMethods({2, 3, 5}, {2, 2, 2}, {1, 1, 1}, {2, 1, 0});
    errCount += testConversionMethods({2, 3, 5}, {1, 1, 1}, {2, 2, 2}, {0, 1, 2});
    errCount += testConversionMethods({2, 3, 5}, {1, 1, 1}, {2, 2, 2}, {2, 1, 0});
    errCount += testConversionMethods({2, 3, 5}, {2, 2, 2}, {1, 1, 1}, {2, 0, 1});
    errCount += testConversionMethods({2, 3, 5}, {2, 2, 2}, {1, 1, 1}, {1, 0, 2});
    errCount += testConversionMethods({2, 3, 5}, {2, 2, 2}, {1, 1, 1}, {1, 2, 0});
  }

  {
    IndexCoords domainShape {5, 3, 2};
    IndexCoords loPads {2, 2, 2};
    IndexCoords hiPads {1, 1, 1};
    IndexCoords strideOrder {0, 1, 2};
    errCount += testByConduitExample(domainShape, loPads, hiPads, strideOrder);
  }

  {
    IndexCoords domainShape {5, 3, 2};
    IndexCoords loPads {2, 2, 2};
    IndexCoords hiPads {1, 1, 1};
    IndexCoords strideOrder {2, 1, 0};
    errCount += testByConduitExample(domainShape, loPads, hiPads, strideOrder);
  }

  {
    IndexCoords domainShape {5, 3, 2};
    IndexCoords loPads {2, 1, 0};
    IndexCoords hiPads {1, 0, 2};
    IndexCoords strideOrder {2, 0, 1};
    errCount += testByConduitExample(domainShape, loPads, hiPads, strideOrder);
  }

  SLIC_INFO(axom::fmt::format("Test found {} errors.", errCount));

  finalizeLogger();

  return (errCount != 0);
}

#endif  // AXOM_USE_CONDUIT
