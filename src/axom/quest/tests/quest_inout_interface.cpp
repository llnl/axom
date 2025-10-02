// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include "axom/quest/interface/inout.hpp"
// for test that reads in a mesh without the interface
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include <string>

/// Helper class to wrap a template for the dimension.
/// Appears to be required to use non-type template parameters with TYPED_TEST
template <int NDIMS>
struct DimWrapper
{
  static const int DIM = NDIMS;
};
template <int NDIMS>
const int DimWrapper<NDIMS>::DIM;

/// Test fixture for quest::inout_query interface
template <typename DimWrapper>
class InOutInterfaceTest : public ::testing::Test
{
public:
  static const int DIM = DimWrapper::DIM;

#ifdef AXOM_USE_C2C
  static_assert(DIM == 2 || DIM == 3, "Dimension must be either 2 or 3");
#else
  static_assert(DIM == 3, "Dimension must be 3");
#endif

  using CoordType = double;
  using InOutPoint = axom::primal::Point<CoordType, DIM>;
  using InOutBBox = axom::primal::BoundingBox<CoordType, DIM>;

protected:
  virtual void SetUp()
  {
    // Setup meshfile
#ifdef AXOM_DATA_DIR
    if(DIM == 2)
    {
      namespace fs = axom::utilities::filesystem;
      const std::string DATA_DIR = fs::joinPath(AXOM_DATA_DIR, "contours");
      const std::string fileName = "unit_circle.contour";
      meshfile = fs::joinPath(DATA_DIR, fileName);
    }
    else  // DIM == 3
    {
      namespace fs = axom::utilities::filesystem;
      const std::string DATA_DIR = fs::joinPath(AXOM_DATA_DIR, "quest");
      const std::string fileName = "sphere.stl";
      meshfile = fs::joinPath(DATA_DIR, fileName);
    }
#else
    FAIL() << "quest_inout_interface test requires AXOM_DATA_DIR";
#endif
  }

  std::string meshfile;
};

// ----------------------------------------------------------------------------
// Set up the list of types we want to test.
// ----------------------------------------------------------------------------
#ifdef AXOM_USE_C2C
using MyDims = ::testing::Types<DimWrapper<2>, DimWrapper<3>>;
#else
using MyDims = ::testing::Types<DimWrapper<3>>;
#endif
TYPED_TEST_SUITE(InOutInterfaceTest, MyDims);

// ----------------------------------------------------------------------------

TYPED_TEST(InOutInterfaceTest, initialize_and_finalize)
{
  const int DIM = TestFixture::DIM;

  EXPECT_TRUE(axom::utilities::filesystem::pathExists(this->meshfile));

  // InOut begins uninitialized
  EXPECT_FALSE(axom::quest::inout_initialized());

  // Initialize the InOut query
  EXPECT_EQ(0, axom::quest::inout_set_dimension(DIM));
  EXPECT_EQ(0, axom::quest::inout_init(this->meshfile));

  // InOut should now be initialized
  EXPECT_TRUE(axom::quest::inout_initialized());

  // Finalize the InOut query
  EXPECT_EQ(0, axom::quest::inout_finalize());

  // InOut should no longer  be initialized
  EXPECT_FALSE(axom::quest::inout_initialized());
}

TYPED_TEST(InOutInterfaceTest, logger_inited)
{
  const int DIM = TestFixture::DIM;
  const bool origSlicInited = axom::slic::isInitialized();

  auto origLogLevel = axom::slic::getLoggingMsgLevel();

  // Initialize the InOut query
  EXPECT_EQ(0, axom::quest::inout_set_dimension(DIM));
  EXPECT_EQ(0, axom::quest::inout_init(this->meshfile));

  // slic is initialized now
  EXPECT_TRUE(axom::slic::isInitialized());

  // Finalize the InOut query
  EXPECT_EQ(0, axom::quest::inout_finalize());

  // After finalizing, slic should be in the same status that it was originally
  EXPECT_EQ(origSlicInited, axom::slic::isInitialized());

  if(origSlicInited)
  {
    EXPECT_EQ(origLogLevel, axom::slic::getLoggingMsgLevel());
  }
}

TYPED_TEST(InOutInterfaceTest, initialize_from_mesh)
{
  const int DIM = TestFixture::DIM;
  const int failCode = axom::quest::QUEST_INOUT_FAILED;

  EXPECT_TRUE(axom::utilities::filesystem::pathExists(this->meshfile));

  axom::mint::Mesh* tmpMeshPtr {nullptr};

  int rc = failCode;

  if(DIM == 2)
  {
#ifdef AXOM_USE_C2C
    int segmentsPerKnotSpan = 10;
    double weldThreshold = 1E-9;
    double revolvedVolume = 0.;
    const bool uniform = true;
    const double percentError = 0.;
    auto identity = axom::numerics::Matrix<double>::identity(4);
    rc = axom::quest::internal::read_c2c_mesh(this->meshfile,
                                              uniform,
                                              identity,
                                              segmentsPerKnotSpan,
                                              weldThreshold,
                                              percentError,
                                              tmpMeshPtr,
                                              revolvedVolume);
#endif  // AXOM_USE_C2C
  }
  else  // DIM == 3
  {
    rc = axom::quest::internal::read_stl_mesh(this->meshfile, tmpMeshPtr);
  }

  std::shared_ptr<axom::mint::Mesh> mesh {tmpMeshPtr};
  tmpMeshPtr = nullptr;

  EXPECT_EQ(0, rc);

  ASSERT_NE(nullptr, mesh);
  EXPECT_GT(mesh->getNumberOfNodes(), 0);
  EXPECT_GT(mesh->getNumberOfCells(), 0);

  // Initialize the InOut query
  EXPECT_EQ(0, axom::quest::inout_set_dimension(DIM));
  EXPECT_EQ(0, axom::quest::inout_init(mesh));

  // InOut should now be initialized
  EXPECT_TRUE(axom::quest::inout_initialized());

  // Finalize the InOut query
  EXPECT_EQ(0, axom::quest::inout_finalize());

  // InOut should no longer  be initialized
  EXPECT_FALSE(axom::quest::inout_initialized());
}

TYPED_TEST(InOutInterfaceTest, query_properties)
{
  const int failCode = axom::quest::QUEST_INOUT_FAILED;
  const int successCode = axom::quest::QUEST_INOUT_SUCCESS;

  using PointType = typename TestFixture::InOutPoint;
  using BBoxType = typename TestFixture::InOutBBox;
  PointType lo, hi, cm;

  // first, test before initializing
  {
    EXPECT_FALSE(axom::quest::inout_initialized());

    // The following should return an error since initialized is false
    SLIC_INFO("--[==[ \n"
              << "\t The following three calls might emit warning messages."
              << " This is expected.\n");

    EXPECT_EQ(failCode, axom::quest::inout_mesh_min_bounds(lo.data()));
    EXPECT_EQ(failCode, axom::quest::inout_mesh_max_bounds(hi.data()));
    EXPECT_EQ(failCode, axom::quest::inout_mesh_center_of_mass(cm.data()));
    SLIC_INFO("--]==]");
  }

  // next, test after initialization
  {
    const int DIM = TestFixture::DIM;
    EXPECT_EQ(0, axom::quest::inout_set_dimension(DIM));
    axom::quest::inout_init(this->meshfile);

    EXPECT_EQ(successCode, axom::quest::inout_mesh_min_bounds(lo.data()));
    EXPECT_EQ(successCode, axom::quest::inout_mesh_max_bounds(hi.data()));
    EXPECT_EQ(successCode, axom::quest::inout_mesh_center_of_mass(cm.data()));

    const int expectedSpatialDim = DIM;
    int spatialDim = axom::quest::inout_get_dimension();
    EXPECT_EQ(expectedSpatialDim, spatialDim);

    SLIC_INFO("Mesh bounding box is " << BBoxType(lo, hi));
    SLIC_INFO("Mesh center of mass is " << cm);
    SLIC_INFO("Spatial dimension of query is " << spatialDim);

    axom::quest::inout_finalize();
  }
}

TYPED_TEST(InOutInterfaceTest, set_params)
{
  const int DIM = TestFixture::DIM;
  const int failCode = axom::quest::QUEST_INOUT_FAILED;
  const int successCode = axom::quest::QUEST_INOUT_SUCCESS;
  const double EPS = 1E-6;

  // Set parameters at the appropriate time
  {
    EXPECT_FALSE(axom::quest::inout_initialized());

    EXPECT_EQ(successCode, axom::quest::inout_set_verbose(true));
    EXPECT_EQ(successCode, axom::quest::inout_set_vertex_weld_threshold(EPS));
    EXPECT_EQ(successCode, axom::quest::inout_set_dimension(DIM));
    // The following is not used in 3D, but we can still invoke it
    EXPECT_EQ(successCode, axom::quest::inout_set_segments_per_knot_span(10));
  }

  // Initialize the query
  {
    SLIC_INFO("*** About to initialize with verbose output.");
    axom::quest::inout_init(this->meshfile);
    SLIC_INFO("*** End verbose output.\n");
  }

  // Set parameters at an inappropriate time
  {
    SLIC_INFO("*** The following calls might emit warning messages.");
    SLIC_INFO("--[==[");

    EXPECT_EQ(failCode, axom::quest::inout_set_verbose(true));
    EXPECT_EQ(failCode, axom::quest::inout_set_vertex_weld_threshold(EPS));
    EXPECT_EQ(failCode, axom::quest::inout_set_dimension(DIM));
    // The following is not used in 3D, but we can still invoke it, and get a warning
    EXPECT_EQ(failCode, axom::quest::inout_set_segments_per_knot_span(10));

    SLIC_INFO("--]==]");
  }

  axom::quest::inout_finalize();
}

TYPED_TEST(InOutInterfaceTest, query)
{
  using PointType = typename TestFixture::InOutPoint;
  PointType query;

  const int DIM = TestFixture::DIM;
  EXPECT_EQ(0, axom::quest::inout_set_verbose(true));
  EXPECT_EQ(0, axom::quest::inout_set_dimension(DIM));
  axom::quest::inout_init(this->meshfile);

  // test an inside point
  EXPECT_TRUE(axom::quest::inout_evaluate(0, 0, 0));

  // test an outside point
  EXPECT_FALSE(axom::quest::inout_evaluate(10, 10, 10));

  axom::quest::inout_finalize();
}

int main(int argc, char** argv)
{
#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
#endif
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  int result = RUN_ALL_TESTS();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return result;
}
