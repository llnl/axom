// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IntersectionShaper.hpp
 *
 * \brief Helper class for intersection-based shaping queries
 */

#ifndef AXOM_QUEST_INTERSECTION_SHAPER__HPP_
#define AXOM_QUEST_INTERSECTION_SHAPER__HPP_

#include "axom/config.hpp"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/sidre/core/Group.hpp"
#include "axom/sidre/core/View.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/mint/utils/vtk_utils.hpp"
#include "axom/klee.hpp"
#include "axom/spin/BVH.hpp"

#include "axom/quest/detail/shaping/shaping_helpers.hpp"
#include "axom/quest/Shaper.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include "axom/fmt.hpp"

#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

#include "axom/fmt.hpp"

#if defined(AXOM_USE_CONDUIT)
  #include "conduit_node.hpp"
  #include "conduit_blueprint_mesh.hpp"
  #include "conduit_blueprint_mcarray.hpp"
#endif

namespace axom
{
namespace quest
{

#if defined(AXOM_USE_64BIT_INDEXTYPE) && !defined(AXOM_NO_INT64_T)
  #if defined(AXOM_USE_CONDUIT)
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType = conduit::DataType::INT64_ID;
  #endif
#else
  #if defined(AXOM_USE_CONDUIT)
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType = conduit::DataType::INT32_ID;
  #endif
#endif

/*!
 * \class TempArrayView
 *
 * \brief Provides a view over an MFEM grid function. MFEM grid functions are
 *        assumed to live in host memory. This class performs data movement
 *        needed to access the grid function data within a GPU device lambda. This
 *        view is limited in scope, though could be expanded in the future.
 *
 * \tparam ExecSpace The execution space where the grid function data will
 *                   be accessed.
 */
template <typename ExecSpace>
class TempArrayView
{
public:
  /*!
   * \brief Host constructor that accepts data in an Array.
   *
   * \param gf The data that will be accessed/modified by this view.
   * \param _needResult Whether the data needs to be brought back to the host
   *                    from the device.
   */
  AXOM_HOST TempArrayView(axom::Array<double>& gf, bool _needResult = true)
  {
    initialize(gf.data(), gf.size(), _needResult);
  }
  /*!
   * \brief Host constructor that accepts data in an ArrayView.
   *
   * \param gf The data that will be accessed/modified by this view.
   * \param _needResult Whether the data needs to be brought back to the host
   *                    from the device.
   */
  AXOM_HOST TempArrayView(axom::ArrayView<double>& gf, bool _needResult = true)
  {
    initialize(gf.data(), gf.size(), _needResult);
  }

  /*!
   * \brief Copy constructor, which is called to make a copy of the host
   *        object so it is accessible inside a RAJA kernel. Any data movement
   *        happened in the host constructor. This version sets hostData to
   *        nullptr so we know not to clean up in the destructor.
   */
  AXOM_HOST_DEVICE TempArrayView(const TempArrayView& obj)
    : m_hostData(nullptr)
    , m_deviceData(obj.m_deviceData)
    , m_numElements(obj.m_numElements)
    , m_needResult(obj.m_needResult)
  { }

  /*!
   * \brief Destructor. On the host, this method may move data from the 
            device and deallocate device storage.
   */
  AXOM_HOST_DEVICE ~TempArrayView() { finalize(); }

  /*!
   * \brief Indexing operator for accessing the data.
   *
   * \param i The index at which to access the data.
   *
   * \return A reference to the data at index i.
   */
  AXOM_HOST_DEVICE double& operator[](int i) { return m_deviceData[i]; }
  // non-const return on purpose.
  AXOM_HOST_DEVICE double& operator[](int i) const { return m_deviceData[i]; }

private:
  /*!
   * \brief Initializes members using data from the grid function. This method
   *        is called on the host.
   *
   * \param hostPtr The grid function data pointer on the host.
   * \param nElem   The grid function size.
   * \param _needResult Whether any data are copied from device.
   */
  AXOM_HOST void initialize(double* hostPtr, int nElem, bool _needResult)
  {
    m_hostData = m_deviceData = hostPtr;
    m_numElements = nElem;
    m_needResult = _needResult;
  }

  /*!
   * \brief Helps during destruction.
   */
  AXOM_HOST_DEVICE void finalize() { m_deviceData = nullptr; }

#if (defined(AXOM_RUNTIME_POLICY_USE_CUDA) || defined(AXOM_RUNTIME_POLICY_USE_HIP))
  /*!
   * \brief Initializes members using data from the grid function. This method
   *        is called on the host and it copies data to the device.
   *
   * \param hostPtr The grid function data pointer on the host.
   * \param nElem   The grid function size.
   * \param _needResult Whether any data are copied from device.
   */
  AXOM_HOST void initializeDevice(double* hostPtr, int nElem, bool _needResult)
  {
    m_hostData = hostPtr;
    m_numElements = nElem;
    m_needResult = _needResult;
    int execSpaceAllocatorID = axom::execution_space<ExecSpace>::allocatorID();

    auto dataSize = sizeof(double) * m_numElements;
    m_deviceData = axom::allocate<double>(dataSize, execSpaceAllocatorID);
    axom::copy(m_deviceData, m_hostData, dataSize);
  }

  /*!
   * \brief Helps during destruction. On the host, it copies device data back
   *        into the grid function on the host.
   */
  AXOM_HOST_DEVICE void finalizeDevice()
  {
  #ifndef AXOM_DEVICE_CODE
    // Only the host will do this work.
    if(m_hostData != nullptr)
    {
      if(m_needResult)
      {
        auto dataSize = sizeof(double) * m_numElements;
        axom::copy(m_hostData, m_deviceData, dataSize);
      }
      axom::deallocate(m_deviceData);
      m_deviceData = nullptr;
    }
  #endif
  }
#endif

private:
  double* m_hostData {nullptr};
  double* m_deviceData {nullptr};
  int m_numElements {0};
  bool m_needResult {false};
};

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
/*!
 * \brief CUDA specialization that calls initializeDevice to copy data
 *        from the host to the device.
 *
 * \param hostPtr The grid function data pointer on the host.
 * \param nElem   The grid function size.
 * \param _needResult Whether any data are copied from device.
 */
template <>
AXOM_HOST inline void TempArrayView<cuda_exec>::initialize(double* hostPtr, int nElem, bool _needResult)
{
  initializeDevice(hostPtr, nElem, _needResult);
}

/*!
 * \brief CUDA specialization that may copy data back from the device
 *        and deallocate any associated device data.
 */
template <>
AXOM_HOST_DEVICE inline void TempArrayView<cuda_exec>::finalize()
{
  finalizeDevice();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
/*!
 * \brief HIP specialization that calls initializeDevice to copy data
 *        from the host to the device.
 *
 * \param hostPtr The grid function data pointer on the host.
 * \param nElem   The grid function size.
 * \param _needResult Whether any data are copied from device.
 */
template <>
AXOM_HOST inline void TempArrayView<hip_exec>::initialize(double* hostPtr, int nElem, bool _needResult)
{
  initializeDevice(hostPtr, nElem, _needResult);
}

/*!
 * \brief HIP specialization that may copy data back from the device
 *        and deallocate any associated device data.
 */
template <>
AXOM_HOST_DEVICE inline void TempArrayView<hip_exec>::finalize()
{
  finalizeDevice();
}
#endif

//---------------------------------------------------------------------------
/**
 * \class
 * \brief Intersects a solid with an input mesh to produce volume fractions.
 * The solid may be any supported by the \c klee::Geometry class.
 *
 * The IntersectionShaper generates material volume fractions:
 *
 * - (3D) For c2c, an input set of 2D contours and replacement rules. Each contour
 *   covers an area from the curve down to the axis of revolution about which
 *   the area is revolved to produce a volume. Contours are refined into smaller
 *   linear spans that are revolved to produce a set of truncated cones, which
 *   are divided into a set set of progressively refined octahedra that can be
 *   intersected with the mesh. The octahedra are intersected with the input
 *   mesh to produce volume fractions.
 *
 * - (3D) For tetrahedral mesh (including Pro/E), an input mesh of 3D tetrahedra is loaded in.
 *   Each tetrahedron has its own respective volume. The tetrahedra are
 *   intersected with the input mesh to produce volume fractions.
 *
 * - (3D) For analytical geometries, the shapes are discretized into a tetrahedral
 *   mesh first.  Sphere and surfaces-of-revolution discretization uses
 *   the refinement level specified in the \c Geometry.
 *
 * - (2D) For c2c, support is a work-in-progress. Initial support for a
 *   single contour that covers an area from the curve down to the x-axis (z-axis).
 *   The contour cannot overlap, and is expected to be entirely above the x-axis.
 *   The contour is refined into smaller linear segments that form triangles with the
 *   x-axis. The trianges are intersected with the input mesh to product volume fractions.
 *
 * - (2D) For triangle mesh, an input STL mesh is loaded in (z-coordinates must be 0).
 *   Each triangle has its own respective area/volume. The triangles are
 *   intersected with the input mesh to produce volume fractions.
 *
 *
 * The input mesh can be an MFEM mesh stored as a \c
 * sidre::MFEMSidreDataCollection or be a Blueprint mesh stored as a
 * \c conduit::Node or a \c sidre::Group.
 *
 * Support for replacement rules exists for MFEM input meshes.
 * Replacement rules for Blueprint meshes is not yet supported.
 * The following comments apply to replacement rules.
 *
 * Volume fractions are represented in the input mesh as a GridFunction with a special prefix,
 * currently "vol_frac_", followed by a material name. Volume fractions
 * can be present in the input data collection prior to shaping and the
 * IntersectionShaper will augment them when changes are needed such as when
 * a material overwrites them. If a new material is not yet represented in
 * the mesh, one will be added.
 *
 * In addition to user-specified materials, the IntersectionShaper creates
 * a "free" material that is used to account for volume fractions that are
 * not assigned to any other material. The free material mainly is used to
 * account for materials when using replacement rules. The free material
 * starts out as all 1's indicating that it contains 100% of all possible
 * material in a zone. Volume fractions for other materials are then
 * subtracted from the free material so no zone exceeds 100% of material.
 */
class IntersectionShaper : public Shaper
{
public:
  using BoundingBox2D = primal::BoundingBox<double, 2>;
  using BoundingBox3D = primal::BoundingBox<double, 3>;
  using HexahedronType = primal::Hexahedron<double, 3>;
  using OctahedronType = primal::Octahedron<double, 3>;
  using PolyhedronType = primal::Polyhedron<double, 3>;
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using TetrahedronType = primal::Tetrahedron<double, 3>;
  using SegmentMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
  using TetMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

  // Use default for MAX_VERTS (8).
  // Assume intersection between triangles and quads,
  // max vertices for overlap is 7.
  using PolygonStaticType = primal::Polygon<double, 2, axom::primal::PolygonArray::Static>;

  using RuntimePolicy = axom::runtime_policy::Policy;

  static constexpr int DEFAULT_CIRCLE_REFINEMENT_LEVEL {7};
  static constexpr double DEFAULT_REVOLVED_VOLUME {0.};

public:
#if defined(AXOM_USE_MFEM)
  /*!
    \brief Construct Shaper to operate on an MFEM mesh.
  */
  IntersectionShaper(RuntimePolicy runtimePolicy,
                     int allocatorId,
                     const klee::ShapeSet& shapeSet,
                     sidre::MFEMSidreDataCollection* dc)
    : Shaper(runtimePolicy, allocatorId, shapeSet, dc)
  {
    m_free_mat_name = "free";
  }
#endif

#if defined(AXOM_USE_CONDUIT)
  /*!
    \brief Construct Shaper to operate on a blueprint-formatted mesh
    stored in a sidre Group.
    \param [in] runtimePolicy A value from RuntimePolicy.
                The simplest policy is RuntimePolicy::seq, which specifies
                running sequentially on the CPU.
    \param [in] allocatorID Data allocator ID.  Choose something compatible
                with \c runtimePolicy.  See \c execution_space.
  */
  IntersectionShaper(RuntimePolicy runtimePolicy,
                     int allocatorId,
                     const klee::ShapeSet& shapeSet,
                     sidre::Group* bpGrp,
                     const std::string& topo = "")
    : Shaper(runtimePolicy, allocatorId, shapeSet, bpGrp, topo)
    , m_free_mat_name("free")
  { }

  /*!
    \brief Construct Shaper to operate on a blueprint-formatted mesh
    stored in a Conduit Node.
  */
  IntersectionShaper(RuntimePolicy runtimePolicy,
                     int allocatorId,
                     const klee::ShapeSet& shapeSet,
                     conduit::Node& bpNode,
                     const std::string& topo = "")
    : Shaper(runtimePolicy, allocatorId, shapeSet, bpNode, topo)
    , m_free_mat_name("free")
  { }
#endif

  //!@brief Set data that depends on mesh (but not on shapes).
  template <typename ShapeType>
  void setMeshDependentData()
  {
    AXOM_ANNOTATE_SCOPE("IntersectionShaper::setMeshDependentData");

    // Setup 2D mesh
    if(std::is_same<ShapeType, PolygonStaticType>::value)
    {
      switch(m_execPolicy)
      {
      case RuntimePolicy::seq:
        setMeshDependentDataImpl2D<seq_exec>();
        break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case RuntimePolicy::omp:
        setMeshDependentDataImpl2D<omp_exec>();
        break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case RuntimePolicy::cuda:
        setMeshDependentDataImpl2D<cuda_exec>();
        break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case RuntimePolicy::hip:
        setMeshDependentDataImpl2D<hip_exec>();
        break;
#endif
      default:
        SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
        break;
      }
    }
    // Setup 3D mesh
    else
    {
      switch(m_execPolicy)
      {
      case RuntimePolicy::seq:
        setMeshDependentDataImpl3D<seq_exec>();
        break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case RuntimePolicy::omp:
        setMeshDependentDataImpl3D<omp_exec>();
        break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case RuntimePolicy::cuda:
        setMeshDependentDataImpl3D<cuda_exec>();
        break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case RuntimePolicy::hip:
        setMeshDependentDataImpl3D<hip_exec>();
        break;
#endif
      default:
        SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
        break;
      }
    }
  }

  /*!
    @brief Set 2D mesh-dependent data, using the given ExecSpace for execution.

    This method has proven to be a potential bottleneck on devices.
    The performance annotations will be removed once it is robustly fixed.
  */
  template <typename ExecSpace>
  void setMeshDependentDataImpl2D()
  {
    populateQuadsFromMesh<ExecSpace>();
    auto quads_device_view = m_quads.view();

    AXOM_ANNOTATE_BEGIN("allocate m_cell_volumes");
    m_cell_volumes = axom::Array<double>(m_cellCount, m_cellCount, m_allocatorId);
    AXOM_ANNOTATE_END("allocate m_cell_volumes");
    m_cell_volumes.fill(0.0);

    SLIC_INFO(axom::fmt::format("{:-^80}", " Calculating quadrilateral element volume "));
    auto cell_volumes_device_view = m_cell_volumes.view();
    AXOM_ANNOTATE_BEGIN("cell_volume");
    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) { cell_volumes_device_view[i] = quads_device_view[i].area(); });
    AXOM_ANNOTATE_END("cell_volume");

    AXOM_ANNOTATE_BEGIN("populate m_quad_bbs");
    m_quad_bbs = axom::Array<BoundingBox2D>(m_cellCount, m_cellCount, m_allocatorId);
    axom::ArrayView<BoundingBox2D> quad_bbs_device_view = m_quad_bbs.view();

    // Get bounding boxes for quadrilateral elements
    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        BoundingBox2D res;

        int num_verts = quads_device_view[i].numVertices();
        for(int j = 0; j < num_verts; ++j)
        {
          res.addPoint(quads_device_view[i][j]);
        }

        quad_bbs_device_view[i] = res;
      });
    AXOM_ANNOTATE_END("populate m_quad_bbs");

    AXOM_ANNOTATE_BEGIN("allocate m_overlap_volumes");
    m_overlap_volumes = axom::Array<double>(m_cellCount, m_cellCount, m_allocatorId);
    AXOM_ANNOTATE_END("allocate m_overlap_volumes");
  }

  /*!
    @brief Set 3D mesh-dependent data, using the given ExecSpace for execution.

    This method has proven to be a potential bottleneck on devices.
    The performance annotations will be removed once it is robustly fixed.
  */
  template <typename ExecSpace>
  void setMeshDependentDataImpl3D()
  {
    constexpr int NUM_TETS_PER_HEX = 24;

    AXOM_ANNOTATE_BEGIN("allocate m_tets_from_hexes_device");
    m_tets_from_hexes_device = axom::Array<TetrahedronType>(ArrayOptions::Uninitialized(),
                                                            m_cellCount * NUM_TETS_PER_HEX,
                                                            m_cellCount * NUM_TETS_PER_HEX,
                                                            m_allocatorId);
    AXOM_ANNOTATE_END("allocate m_tets_from_hexes_device");

    populateHexesFromMesh<ExecSpace>();
    auto hexesView = m_hexes.view();

    AXOM_ANNOTATE_BEGIN("allocate m_cell_volumes");
    m_cell_volumes = axom::Array<double>(m_cellCount, m_cellCount, m_allocatorId);
    AXOM_ANNOTATE_END("allocate m_cell_volumes");
    m_cell_volumes.fill(0.0);

    SLIC_INFO(axom::fmt::format("{:-^80}", " Calculating hexahedron element volume "));
    auto cellVolumesView = m_cell_volumes.view();
    AXOM_ANNOTATE_BEGIN("cell_volume");
    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) { cellVolumesView[i] = hexesView[i].volume(); });
    AXOM_ANNOTATE_END("cell_volume");

    SLIC_INFO(
      axom::fmt::format("{:-^80}", " Decomposing each hexahedron element into 24 tetrahedrons "));

    AXOM_ANNOTATE_BEGIN("populate m_hex_bbs");
    m_hex_bbs = axom::Array<BoundingBox3D>(m_cellCount, m_cellCount, m_allocatorId);

    // Get bounding boxes for hexahedral elements
    axom::ArrayView<BoundingBox3D> hexBbsView = m_hex_bbs.view();
    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        hexBbsView[i] = primal::compute_bounding_box<double, 3>(hexesView[i]);
      });  // end of loop to initialize hexahedral elements and bounding boxes
    AXOM_ANNOTATE_END("populate m_hex_bbs");

    using TetHexArray = axom::StackArray<TetrahedronType, NUM_TETS_PER_HEX>;

    auto tetsFromHexesView = m_tets_from_hexes_device.view();
    AXOM_ANNOTATE_BEGIN("init_tets");
    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        TetHexArray cur_tets;
        hexesView[i].triangulate(cur_tets);

        for(int j = 0; j < NUM_TETS_PER_HEX; j++)
        {
          tetsFromHexesView[i * NUM_TETS_PER_HEX + j] = cur_tets[j];
        }
      });
    AXOM_ANNOTATE_END("init_tets");

    AXOM_ANNOTATE_BEGIN("allocate m_overlap_volumes");
    m_overlap_volumes = axom::Array<double>(m_cellCount, m_cellCount, m_allocatorId);
    AXOM_ANNOTATE_END("allocate m_overlap_volumes");
  }

  ///@{
  //!  @name Functions to get and set shaping parameters related to intersection; supplements parameters in base class

  void setLevel(int level) { m_level = level; }

  /*!
   * \brief Set the name of the material used to account for free volume fractions.
   * \param name The new name of the material. This name cannot contain 
   *             underscores and it cannot be set once shaping has started.
   * \note This should not be called once any shaping has occurred.
   */
  void setFreeMaterialName(const std::string& name)
  {
    if(name.find("_") != std::string::npos)
    {
      SLIC_ERROR("The free material name cannot contain underscores.");
    }
    if(m_cellCount > 0)
    {
      SLIC_ERROR("The free material name cannot be set once shaping has occurred.");
    }
    m_free_mat_name = name;
  }
  ///@}

  /*!
   * \brief Return the revolved volume that was computed during dynamic refinement.
   * \return The revolved volume (or zero).
   */
  double getRevolvedVolume() const { return m_revolvedVolume; }

  /*!
   * \brief Return the revolved volume for the m_surfaceMesh at m_level circle refinement.
   * \note loadShape should have been called before this method.
   * \return The revolved volume (or zero).
   */
  double getApproximateRevolvedVolume() const { return volume(m_surfaceMesh.get(), m_level); }

  virtual void loadShape(const klee::Shape& shape) override
  {
    // Make sure we can store the revolved volume in member m_revolvedVolume.
    loadShapeInternal(shape, m_percentError, m_revolvedVolume);

    // Filter the mesh, store in m_surfaceMesh.
    if(this->shapeFormat(shape) == "c2c")
    {
      SegmentMesh* newm = filterMesh(dynamic_cast<const SegmentMesh*>(m_surfaceMesh.get()));
      m_surfaceMesh.reset(newm);
    }
  }

// The following private methods are made public due to CUDA compilers
// requirements for methods that call device functions.
#if defined(__CUDACC__)
public:
#else
private:
#endif

  ///@{
  //!  @name Private functions related to the stages for a given shape

  template <typename ExecSpace>
  void prepareTriCells()
  {
    const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    const int device_allocator = axom::execution_space<ExecSpace>::allocatorID();

    // Number of triangles in mesh
    m_tricount = m_surfaceMesh->getNumberOfCells();

    axom::Array<PolygonStaticType> tris_host(m_tricount, m_tricount, host_allocator);

    // Initialize 2D triangles from mesh (ignore z coordinate)
    axom::Array<IndexType> nodeIds(3);

    // Buffer is 3D for stl mesh
    axom::Array<Point3D> pts(3);

    for(int i = 0; i < m_tricount; i++)
    {
      m_surfaceMesh->getCellNodeIDs(i, nodeIds.data());

      m_surfaceMesh->getNode(nodeIds[0], pts[0].data());
      m_surfaceMesh->getNode(nodeIds[1], pts[1].data());
      m_surfaceMesh->getNode(nodeIds[2], pts[2].data());

      // Verify that z-coordinates are unused by the mesh.
      // (0 for stl mesh, undefined by in-memory triangle mesh)
      if(pts[0][2] != 0 || pts[1][2] != 0 || pts[2][2] != 0)
      {
        SLIC_ERROR(axom::fmt::format("2D triangles must have undefined or value 0 z-coordinates"));
      }

      Point2D p1({pts[0][0], pts[0][1]});
      Point2D p2({pts[1][0], pts[1][1]});
      Point2D p3({pts[2][0], pts[2][1]});

      tris_host[i] = PolygonStaticType({p1, p2, p3});
    }

    // Copy triangles to device
    m_tris = axom::Array<PolygonStaticType>(tris_host, device_allocator);

    if(this->isVerbose())
    {
      // Print out the bounding box containing all the triangles
      BoundingBox2D all_tris_bb;
      for(int i = 0; i < m_tricount; i++)
      {
        // Use non-static Polygon to match template
        axom::primal::Polygon<double, 2> tempPoly({tris_host[i][0], tris_host[i][1], tris_host[i][2]});
        all_tris_bb.addBox(primal::compute_bounding_box(tempPoly));
      }
      SLIC_INFO(
        axom::fmt::format("VERBOSE: Bounding box containing all generated triangles "
                          "has dimensions:\n\t{}",
                          all_tris_bb));

      auto tri_device_view = m_tris.view();

      // Print out the total volume of all the triangles
      axom::ReduceSum<ExecSpace, double> total_tri_area(0.0);
      axom::for_all<ExecSpace>(
        m_tricount,
        AXOM_LAMBDA(axom::IndexType i) { total_tri_area += tri_device_view[i].area(); });

      SLIC_INFO(axom::fmt::format("VERBOSE: Total area of all generated triangles is {}",
                                  total_tri_area.get()));

      // Check if any Triangles are degenerate with zero area
      axom::ReduceSum<ExecSpace, int> num_degenerate(0);
      axom::for_all<ExecSpace>(
        m_tricount,
        AXOM_LAMBDA(axom::IndexType i) {
          if(axom::utilities::isNearlyEqual(tri_device_view[i].area(), 0.0))
          {
            num_degenerate += 1;
          }
        });

      SLIC_INFO(axom::fmt::format("VERBOSE: Degenerate {} triangles found with zero area",
                                  num_degenerate.get()));

    }  // end of verbose output for triangles

  }  // end of prepareTriCells()

  // Prepares the tet mesh cells for the spatial index
  template <typename ExecSpace>
  void prepareTetCells()
  {
    const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    const int device_allocator = m_allocatorId;

    // Number of tets in mesh
    m_tetcount = m_surfaceMesh->getNumberOfCells();

    axom::Array<TetrahedronType> tets_host(m_tetcount, m_tetcount, host_allocator);

    // Initialize tetrahedra
    axom::Array<IndexType> nodeIds(4);
    axom::Array<Point3D> pts(4);

    for(int i = 0; i < m_tetcount; i++)
    {
      m_surfaceMesh->getCellNodeIDs(i, nodeIds.data());

      m_surfaceMesh->getNode(nodeIds[0], pts[0].data());
      m_surfaceMesh->getNode(nodeIds[1], pts[1].data());
      m_surfaceMesh->getNode(nodeIds[2], pts[2].data());
      m_surfaceMesh->getNode(nodeIds[3], pts[3].data());

      tets_host[i] = TetrahedronType({pts[0], pts[1], pts[2], pts[3]});
    }

    // Copy tets to device
    m_tets = axom::Array<TetrahedronType>(tets_host, device_allocator);

    if(this->isVerbose())
    {
      // Print out the bounding box containing all the tetrahedra
      BoundingBox3D all_tet_bb;
      for(int i = 0; i < m_tetcount; i++)
      {
        all_tet_bb.addBox(primal::compute_bounding_box(tets_host[i]));
      }
      SLIC_INFO(
        axom::fmt::format("VERBOSE: Bounding box containing all generated tetrahedra "
                          "has dimensions:\n\t{}",
                          all_tet_bb));

      auto tets_device_view = m_tets.view();

      // Print out the total volume of all the tetrahedra
      axom::ReduceSum<ExecSpace, double> total_tet_vol(0.0);
      axom::for_all<ExecSpace>(
        m_tetcount,
        AXOM_LAMBDA(axom::IndexType i) { total_tet_vol += tets_device_view[i].volume(); });

      SLIC_INFO(axom::fmt::format("VERBOSE: Total volume of all generated tetrahedra is {}",
                                  total_tet_vol.get()));

      // Check if any Tetrahedron are degenerate with zero volume
      axom::ReduceSum<ExecSpace, int> num_degenerate(0);
      axom::for_all<ExecSpace>(
        m_tetcount,
        AXOM_LAMBDA(axom::IndexType i) {
          if(tets_device_view[i].degenerate())
          {
            num_degenerate += 1;
          }
        });

      SLIC_INFO(axom::fmt::format("VERBOSE: Degenerate {} tetrahedra found with zero volume",
                                  num_degenerate.get()));

      // Dump tet mesh as a vtk mesh
      axom::mint::write_vtk(m_surfaceMesh.get(), "proe_tet.vtk");

    }  // end of verbose output for Pro/E
  }

  // Prepares the C2C mesh cells for the spatial index.
  // Produced octahedra for the revolved contour in the 3D case,
  // or triangles for the contour bounded by the x-axis in the 2D case.
  template <typename ExecSpace>
  void prepareC2CCells()
  {
    const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

    // Number of points in polyline
    int pointcount = getSurfaceMesh()->getNumberOfNodes();

    axom::Array<Point2D> polyline(pointcount, pointcount);

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format(axom::utilities::locale(),
                        " Checking contour with {:L} points for degenerate segments",
                        pointcount)));

    // The mesh points are filtered like we want. We need only copy
    // them into the polyline array.
    for(int i = 0; i < pointcount; ++i)
    {
      m_surfaceMesh->getNode(i, polyline[i].data());
    }
    int polyline_size = pointcount;

    int compMeshDim = getCompMeshDim();

    // Produce octahedra for revolved contour (3D)
    if(compMeshDim == 3)
    {
      SLIC_INFO(
        axom::fmt::format("{:-^80}", axom::fmt::format(" Refinement level set to {} ", m_level)));

      // Generate the Octahedra
      // (octahedra m_octs will be on device)
      const bool disc_status =
        axom::quest::discretize<ExecSpace>(polyline, polyline_size, m_level, m_octs, m_octcount);

      axom::ArrayView<OctahedronType> octs_device_view = m_octs.view();

      AXOM_UNUSED_VAR(disc_status);  // silence warnings in release configs
      SLIC_ASSERT_MSG(disc_status,
                      "Discretization of contour has failed. Check that contour is valid");

      SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                  "Contour has been discretized into {:L} octahedra ",
                                  m_octcount));

      if(this->isVerbose())
      {
        // Print out the bounding box containing all the octahedra
        BoundingBox3D all_oct_bb;
        axom::Array<OctahedronType> octs_host = axom::Array<OctahedronType>(m_octs, host_allocator);
        auto octs_host_view = octs_host.view();

        for(int i = 0; i < m_octcount; i++)
        {
          all_oct_bb.addBox(primal::compute_bounding_box(octs_host[i]));
        }
        SLIC_INFO(
          axom::fmt::format("VERBOSE: Bounding box containing all generated octahedra "
                            "has dimensions:\n\t{}",
                            all_oct_bb));

        // Print out the total volume of all the octahedra
        axom::ReduceSum<ExecSpace, double> total_oct_vol(0.0);
        axom::for_all<ExecSpace>(
          m_octcount,
          AXOM_LAMBDA(axom::IndexType i) {
            // Convert Octahedron into Polyhedron
            PolyhedronType octPoly = PolyhedronType::from_primitive(octs_device_view[i]);

            total_oct_vol += octPoly.volume();
          });

        SLIC_INFO(axom::fmt::format("VERBOSE: Total volume of all generated octahedra is {}",
                                    total_oct_vol.get()));

        // Check if any Octahedron are degenerate with all points {0,0,0}
        axom::ReduceSum<ExecSpace, int> num_degenerate(0);

        const int device_allocator = m_allocatorId;
        axom::Array<OctahedronType> degenerate_oct_host(1, 1, host_allocator);
        degenerate_oct_host[0] = OctahedronType();
        axom::Array<OctahedronType> degenerate_oct_device =
          axom::Array<OctahedronType>(degenerate_oct_host, device_allocator);
        auto degenerate_oct_device_view = degenerate_oct_device.view();

        axom::for_all<ExecSpace>(
          m_octcount,
          AXOM_LAMBDA(axom::IndexType i) {
            if(octs_device_view[i].equals(degenerate_oct_device_view[0]))
            {
              num_degenerate += 1;
            }
          });

        SLIC_INFO(axom::fmt::format("VERBOSE: {} Octahedron found with all points (0,0,0)",
                                    num_degenerate.get()));

        // Dump discretized octs as a tet mesh
        axom::mint::Mesh* tetmesh;
        axom::quest::mesh_from_discretized_polyline(octs_host_view,
                                                    m_octcount,
                                                    polyline_size - 1,
                                                    tetmesh);
        axom::mint::write_vtk(tetmesh, "discretized_surface_of_revolution.vtk");
        delete tetmesh;

      }  // end of verbose output for contour
    }  // end of 3D case

    // Produce triangles for contour bounded by x-axis (2D)
    else
    {
      SLIC_INFO(
        axom::fmt::format("{:-^80}",
                          axom::fmt::format(axom::utilities::locale(),
                                            "{:L} linear segments to generate per NURBS knot span",
                                            m_samplesPerKnotSpan)));

      const int device_allocator = axom::execution_space<ExecSpace>::allocatorID();

      // Number of triangles in mesh (2 triangles per segment/quad cell)
      m_tricount = m_surfaceMesh->getNumberOfCells() * 2;

      axom::Array<PolygonStaticType> tris_host(m_tricount, m_tricount, host_allocator);

      // Initialize 2D triangles from segment mesh (3rd point is on the x-axis)
      axom::Array<IndexType> nodeIds(2);

      // Buffer to store 2D points
      axom::Array<Point2D> pts(2);

      for(int i = 0; i < m_tricount / 2; i++)
      {
        m_surfaceMesh->getCellNodeIDs(i, nodeIds.data());

        m_surfaceMesh->getNode(nodeIds[0], pts[0].data());
        m_surfaceMesh->getNode(nodeIds[1], pts[1].data());

        constexpr double EPS = 1e-10;
        if(axom::primal::detail::isLt(pts[0][1], 0.0, EPS) ||
           axom::primal::detail::isLt(pts[1][1], 0.0, EPS))
        {
          SLIC_ERROR("2D: Contour with non-negative r values expected");
        }

        // Right ear of quad
        tris_host[i * 2] = PolygonStaticType({pts[1], pts[0], Point2D {pts[0][0], 0}});

        // Left ear of quad
        tris_host[(i * 2) + 1] = PolygonStaticType({pts[1], pts[0], Point2D {pts[1][0], 0}});
      }

      // Copy triangles to device
      m_tris = axom::Array<PolygonStaticType>(tris_host, device_allocator);

      SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                  "Contour has been discretized into {:L} triangles ",
                                  m_tricount));

    }  // end of 2D case
  }

  /// Initializes the spatial index for shaping
  template <typename ExecSpace>
  void prepareShapeQueryImpl(klee::Dimensions shapeDimension, const klee::Shape& shape)
  {
    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format("Running intersection-based shaper in execution Space: {}",
                        axom::execution_space<ExecSpace>::name())));

    const auto& shapeName = shape.getName();
    AXOM_UNUSED_VAR(shapeDimension);
    AXOM_UNUSED_VAR(shapeName);

    std::string shapeFormat = shape.getGeometry().getFormat();

    // C2C mesh is not discretized into tets, but all others are.
    if(shapeFormat == "c2c")
    {
      prepareC2CCells<ExecSpace>();
    }
    else if(surfaceMeshIsTet())
    {
      prepareTetCells<ExecSpace>();
    }
    // 2D STL Triangle mesh
    else if(shapeFormat == "stl" || surfaceMeshIsTri())
    {
      prepareTriCells<ExecSpace>();
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("The shape format {} is unsupported", shapeFormat));
    }
  }

  template <typename ExecSpace>
  void runShapeQuery2DImpl(const klee::Shape& shape,
                           axom::Array<PolygonStaticType>& shapes,
                           int shape_count)
  {
    AXOM_ANNOTATE_SCOPE("IntersectionShaper::runShapeQuery2DImpl");
    AXOM_UNUSED_VAR(shape);

    if(m_quads.empty())
    {
      setMeshDependentData<PolygonStaticType>();
    }

    const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    const int device_allocator = axom::execution_space<ExecSpace>::allocatorID();

    SLIC_INFO(axom::fmt::format("{:-^80}", " Inserting shapes' bounding boxes into BVH "));

    // Generate the BVH tree over the shapes
    // Access-aligned bounding boxes
    m_aabbs_2d = axom::Array<BoundingBox2D>(shape_count, shape_count, device_allocator);

    axom::ArrayView<PolygonStaticType> shapes_device_view = shapes.view();

    axom::ArrayView<BoundingBox2D> aabbs_device_view = m_aabbs_2d.view();

    // Get the bounding boxes for the shapes
    axom::for_all<ExecSpace>(
      shape_count,
      AXOM_LAMBDA(axom::IndexType i) {
        BoundingBox2D res;

        int num_verts = shapes_device_view[i].numVertices();
        for(int j = 0; j < num_verts; ++j)
        {
          res.addPoint(shapes_device_view[i][j]);
        }

        aabbs_device_view[i] = res;
      });

    // Insert shapes' Bounding Boxes into BVH.
    //bvh.setAllocatorID(poolID);
    spin::BVH<2, ExecSpace, double> bvh;
    bvh.initialize(aabbs_device_view, shape_count);

    SLIC_INFO(axom::fmt::format("{:-^80}", " Querying the BVH tree "));

    auto quad_bbs_device_view = m_quad_bbs.view();

    // Find which shape bounding boxes intersect quadrilateral bounding boxes
    SLIC_INFO(axom::fmt::format("{:-^80}", " Finding shape candidates for each quad element "));

    axom::Array<IndexType> offsets(m_cellCount, m_cellCount, device_allocator);
    axom::Array<IndexType> counts(m_cellCount, m_cellCount, device_allocator);
    axom::Array<IndexType> candidates;
    bvh.findBoundingBoxes(offsets, counts, candidates, m_cellCount, quad_bbs_device_view);

    // Get the total number of candidates
    const auto counts_device_view = counts.view();
    AXOM_ANNOTATE_BEGIN("populate totalCandidates");
    axom::ReduceSum<ExecSpace, int> totalCandidates(0);
    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) { totalCandidates += counts_device_view[i]; });
    AXOM_ANNOTATE_END("populate totalCandidates");

    AXOM_ANNOTATE_BEGIN("allocate scratch space");
    // Initialize quadrilateral indices and shape candidates
    AXOM_ANNOTATE_BEGIN("allocate quad_indices_device");
    axom::Array<IndexType> quad_indices_device(totalCandidates.get(),
                                               totalCandidates.get(),
                                               device_allocator);
    AXOM_ANNOTATE_END("allocate quad_indices_device");
    auto quad_indices_device_view = quad_indices_device.view();

    // Quad elements
    auto quads_device_view = m_quads.view();

    AXOM_ANNOTATE_BEGIN("allocate shape_candidates_device");
    axom::Array<IndexType> shape_candidates_device(totalCandidates.get(),
                                                   totalCandidates.get(),
                                                   device_allocator);
    AXOM_ANNOTATE_END("allocate shape_candidates_device");
    auto shape_candidates_device_view = shape_candidates_device.view();
    AXOM_ANNOTATE_END("allocate scratch space");

    // New total number of candidates after omitting degenerate shapes
    AXOM_ANNOTATE_BEGIN("newTotalCandidates memory");
    axom::Array<IndexType> newTotalCandidates_host(1, 1, host_allocator);
    newTotalCandidates_host[0] = 0;
    axom::Array<IndexType> newTotalCandidates_device =
      axom::Array<IndexType>(newTotalCandidates_host, device_allocator);
    auto newTotalCandidates_device_view = newTotalCandidates_device.view();
    AXOM_ANNOTATE_END("newTotalCandidates memory");

    SLIC_INFO(axom::fmt::format("{:-^80}", " Creating an array of candidate pairs for shaping "));

    const auto offsets_device_view = offsets.view();
    const auto candidates_device_view = candidates.view();
    {
      AXOM_ANNOTATE_SCOPE("init_candidates");
      axom::for_all<ExecSpace>(
        m_cellCount,
        AXOM_LAMBDA(axom::IndexType i) {
          for(int j = 0; j < counts_device_view[i]; j++)
          {
            int shapeIdx = candidates_device_view[offsets_device_view[i] + j];

            IndexType idx =
              axom::atomicAdd<ExecSpace>(&newTotalCandidates_device_view[0], IndexType {1});
            quad_indices_device_view[idx] = i;
            shape_candidates_device_view[idx] = shapeIdx;
          }
        });
    }

    // Overlap volume is the area of clip(tri,quad) for stl mesh
    m_overlap_volumes.fill(0.0);
    axom::ArrayView<double> overlap_volumes_device_view = m_overlap_volumes.view();

    SLIC_INFO(axom::fmt::format("{:-^80}",
                                " Calculating element overlap volume from each quad-shape pair "));

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = true;

    {
      AXOM_ANNOTATE_SCOPE("clipLoop");
      // Copy calculated total back to host
      axom::Array<IndexType> newTotalCandidates_calc_host =
        axom::Array<IndexType>(newTotalCandidates_device, host_allocator);

      axom::for_all<ExecSpace>(
        newTotalCandidates_calc_host[0],
        AXOM_LAMBDA(axom::IndexType i) {
          const int index = quad_indices_device_view[i];
          const int shapeIndex = shape_candidates_device_view[i];

          const PolygonStaticType poly = primal::clip(shapes_device_view[shapeIndex],
                                                      quads_device_view[index],
                                                      EPS,
                                                      tryFixOrientation);

          // Polygon is valid
          if(poly.numVertices() >= 3)
          {
            // Workaround - intermediate volume variable needed for
            // CUDA Pro/E test case correctness
            double area = poly.area();
            axom::atomicAdd<ExecSpace>(overlap_volumes_device_view.data() + index, area);
          }
        });
    }

    axom::ReduceSum<ExecSpace, double> totalOverlap(0);
    axom::ReduceSum<ExecSpace, double> totalQuad(0);

    auto cell_volumes_device_view = m_cell_volumes.view();

    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        totalOverlap += overlap_volumes_device_view[i];
        totalQuad += cell_volumes_device_view[i];
      });

    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Total overlap volume with shape is {:.3Lf}",
                                this->allReduceSum(totalOverlap.get())));
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Total mesh volume is {:.3Lf}",
                                this->allReduceSum(totalQuad.get())));

  }  // end of runShapeQuery2DImpl() function

  /*!
    \tparam ShapeType either TetrahedeonType or OctahedronType.
            depending on whether shape is tet or c2c.
    \param shape the input shape to be superimposed on the mesh.
    \param shapes either the array m_tets (if surface mesh is tet)
           or m_octs (if surface mesh is c2c).
    \param shape_count the count for the shapes array, maintained
           separately from the array.
  */
  template <typename ExecSpace, typename ShapeType>
  void runShapeQuery3DImpl(const klee::Shape& shape, axom::Array<ShapeType>& shapes, int shape_count)

  {
    AXOM_ANNOTATE_SCOPE("IntersectionShaper::runShapeQueryImpl");

    // No need for parameter shape, because it has been converted into
    // m_tets or m_octs, which is what the parameter shapes is.
    AXOM_UNUSED_VAR(shape);

    if(m_hexes.empty())
    {
      setMeshDependentData<ShapeType>();
    }

    const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    const int device_allocator = m_allocatorId;

    constexpr int NUM_TETS_PER_HEX = 24;

    SLIC_INFO(axom::fmt::format("{:-^80}", " Inserting shapes' bounding boxes into BVH "));

    // Generate the BVH tree over the shapes
    // Axis-aligned bounding boxes
    axom::Array<BoundingBox3D> aabbs(shape_count, shape_count, device_allocator);

    axom::ArrayView<ShapeType> shapes_device_view = shapes.view();

    axom::ArrayView<BoundingBox3D> aabbs_device_view = aabbs.view();

    // Get the bounding boxes for the shapes
    axom::for_all<ExecSpace>(
      shape_count,
      AXOM_LAMBDA(axom::IndexType i) {
        aabbs_device_view[i] = primal::compute_bounding_box<double, 3>(shapes_device_view[i]);
      });

    // Insert shapes' Bounding Boxes into BVH.
    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(aabbs_device_view, shape_count);

    SLIC_INFO(axom::fmt::format("{:-^80}", " Querying the BVH tree "));

    axom::ArrayView<const BoundingBox3D> hex_bbs_device_view = m_hex_bbs.view();

    // Find which shape bounding boxes intersect hexahedron bounding boxes
    SLIC_INFO(
      axom::fmt::format("{:-^80}", " Finding shape candidates for each hexahedral element "));

    axom::Array<IndexType> offsets(m_cellCount, m_cellCount, device_allocator);
    axom::Array<IndexType> counts(m_cellCount, m_cellCount, device_allocator);
    axom::Array<IndexType> candidates;
    AXOM_ANNOTATE_BEGIN("bvh.findBoundingBoxes");
    bvh.findBoundingBoxes(offsets, counts, candidates, m_cellCount, hex_bbs_device_view);
    AXOM_ANNOTATE_END("bvh.findBoundingBoxes");

    // Get the total number of candidates
    const auto counts_device_view = counts.view();
    AXOM_ANNOTATE_BEGIN("populate totalCandidates");
    axom::ReduceSum<ExecSpace, int> totalCandidates(0);
    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) { totalCandidates += counts_device_view[i]; });
    AXOM_ANNOTATE_END("populate totalCandidates");

    AXOM_ANNOTATE_BEGIN("allocate scratch space");
    // Initialize hexahedron indices and shape candidates
    AXOM_ANNOTATE_BEGIN("allocate hex_indices");
    axom::Array<IndexType> hex_indices_device(totalCandidates.get() * NUM_TETS_PER_HEX,
                                              totalCandidates.get() * NUM_TETS_PER_HEX,
                                              device_allocator);
    AXOM_ANNOTATE_END("allocate hex_indices");
    auto hex_indices_device_view = hex_indices_device.view();

    AXOM_ANNOTATE_BEGIN("allocate shape_candidates");
    axom::Array<IndexType> shape_candidates_device(totalCandidates.get() * NUM_TETS_PER_HEX,
                                                   totalCandidates.get() * NUM_TETS_PER_HEX,
                                                   device_allocator);
    AXOM_ANNOTATE_END("allocate shape_candidates");
    auto shape_candidates_device_view = shape_candidates_device.view();

    // Tetrahedrons from hexes (24 for each hex)
    axom::ArrayView<TetrahedronType> tets_from_hexes_device_view = m_tets_from_hexes_device.view();

    // Index into 'tets'
    AXOM_ANNOTATE_BEGIN("allocate tet_indices_device");
    axom::Array<IndexType> tet_indices_device(totalCandidates.get() * NUM_TETS_PER_HEX,
                                              totalCandidates.get() * NUM_TETS_PER_HEX,
                                              device_allocator);
    AXOM_ANNOTATE_END("allocate tet_indices_device");
    auto tet_indices_device_view = tet_indices_device.view();
    AXOM_ANNOTATE_END("allocate scratch space");

    // New total number of candidates after omitting degenerate shapes
    AXOM_ANNOTATE_BEGIN("newTotalCandidates memory");
    axom::Array<IndexType> newTotalCandidates_host(1, 1, host_allocator);
    newTotalCandidates_host[0] = 0;
    axom::Array<IndexType> newTotalCandidates_device =
      axom::Array<IndexType>(newTotalCandidates_host, device_allocator);
    auto newTotalCandidates_device_view = newTotalCandidates_device.view();
    AXOM_ANNOTATE_END("newTotalCandidates memory");

    SLIC_INFO(axom::fmt::format("{:-^80}", " Creating an array of candidate pairs for shaping "));

    const auto offsets_device_view = offsets.view();
    const auto candidates_device_view = candidates.view();
    {
      AXOM_ANNOTATE_SCOPE("init_candidates");
      axom::for_all<ExecSpace>(
        m_cellCount,
        AXOM_LAMBDA(axom::IndexType i) {
          for(int j = 0; j < counts_device_view[i]; j++)
          {
            int shapeIdx = candidates_device_view[offsets_device_view[i] + j];

            for(int k = 0; k < NUM_TETS_PER_HEX; k++)
            {
              IndexType idx =
                axom::atomicAdd<ExecSpace>(&newTotalCandidates_device_view[0], IndexType {1});
              hex_indices_device_view[idx] = i;
              shape_candidates_device_view[idx] = shapeIdx;
              tet_indices_device_view[idx] = i * NUM_TETS_PER_HEX + k;
            }
          }
        });
    }

    // Overlap volume is the volume of clip(oct,tet) for c2c
    // or clip(tet,tet) for Pro/E meshes
    m_overlap_volumes.fill(0.0);

    axom::ArrayView<double> overlap_volumes_device_view = m_overlap_volumes.view();

    SLIC_INFO(axom::fmt::format("{:-^80}",
                                " Calculating element overlap volume from each tet-shape pair "));

    constexpr double EPS = 1e-10;
    constexpr bool tryFixOrientation = true;

    {
      AXOM_ANNOTATE_SCOPE("clipLoop");
      // Copy calculated total back to host
      axom::Array<IndexType> newTotalCandidates_calc_host =
        axom::Array<IndexType>(newTotalCandidates_device, host_allocator);

      axom::for_all<ExecSpace>(
        newTotalCandidates_calc_host[0],  // Number of candidates found.
        AXOM_LAMBDA(axom::IndexType i) {
          const int index = hex_indices_device_view[i];
          const int shapeIndex = shape_candidates_device_view[i];
          const int tetIndex = tet_indices_device_view[i];

          const PolyhedronType poly = primal::clip(shapes_device_view[shapeIndex],
                                                   tets_from_hexes_device_view[tetIndex],
                                                   EPS,
                                                   tryFixOrientation);

          // Poly is valid
          if(poly.numVertices() >= 4)
          {
            // Workaround - intermediate volume variable needed for
            // CUDA Pro/E test case correctness
            double volume = poly.volume();
            axom::atomicAdd<ExecSpace>(overlap_volumes_device_view.data() + index, volume);
          }
        });
    }
  }  // end of runShapeQueryImpl() function

  // These methods are private in support of replacement rules.
private:
  /*!
   * \brief Turn a material name into a grid function name.
   *
   * \param materialName The name of the material.
   *
   * \return The name of the material's grid function.
   */
  std::string materialNameToFieldName(const std::string& materialName) const
  {
    return axom::fmt::format("vol_frac_{}", materialName);
  }

  /*!
   * \brief Turn a grid function name into a material name.
   *
   * \param fieldName The name of the grid function.
   *
   * \return The name of the material material.
   */
  std::string fieldNameToMaterialName(const std::string& fieldName) const
  {
    const std::string vol_frac_("vol_frac_");
    std::string name;
    if(fieldName.find(vol_frac_) == 0)
    {
      name = fieldName.substr(vol_frac_.size());
    }
    return name;
  }

  /*!
   * \brief Scans the grid functions in the data collection and creates
   *        a material entry for any that do not already exist. We maintain
   *        our own vectors because we assume that the order of materials
   *        does not change. The mfem::DataCollection uses a map internally
   *        so if we add materials, it could change the traversal order.
   *
   */
  void populateMaterials()
  {
    std::vector<std::string> materialNames = getMaterialNames();
    // Add any of these existing fields to this class' bookkeeping.
    for(const auto& materialName : materialNames)
    {
      (void)getMaterial(materialName);
    }
  }

  // Switch back to public. This is done here because the CUDA compiler
  // does not like the following template functions to be private.

public:
  /*!
   * \brief Make a new grid function that contains all of the free space not
   *        occupied by existing materials.
   *
   * \note We currently leave the free material in the data collection,
   *       even after the shaper has executed.
   *
   * \tparam ExecSpace The execution space where the data are computed.
   *
   * \return The grid function that represents the amount of completely
   *         free space in each zone.
   */
  template <typename ExecSpace>
  axom::ArrayView<double> getCompletelyFree()
  {
    // Add the material prefix so the mesh will automatically
    // consider the free material something it needs to write as a matset.
    const std::string fieldName(materialNameToFieldName(m_free_mat_name));

    bool makeNewData = !hasData(fieldName);
    axom::ArrayView<double> cfgf = getScalarCellData(fieldName);
    SLIC_ASSERT(!cfgf.empty());

    if(makeNewData)
    {
      AXOM_ANNOTATE_SCOPE("compute_free");

      int dataSize = cfgf.size();
      TempArrayView<ExecSpace> cfView(cfgf, true);

      axom::for_all<ExecSpace>(dataSize, AXOM_LAMBDA(axom::IndexType i) { cfView[i] = 1.; });

      // Iterate over all materials and subtract off their VFs from cfgf.
      for(axom::ArrayView<double>& gf : m_vf_grid_functions)
      {
        TempArrayView<ExecSpace> matVFView(gf, false);
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) {
            cfView[i] -= matVFView[i];
            cfView[i] = (cfView[i] < 0.) ? 0. : cfView[i];
          });
      }
    }

    return cfgf;
  }

  /*!
   * \brief Set the volume fractions for the current shape into the grid
   *        function for the material and adjust any other material volume
   *        fraction grid functions.
   *
   * \tparam ExecSpace The execution space where the data are computed.
   * \param shape The shape whose volume fractions are being stored.
   *
   * The replacement rules operate on the current shape's material as well as the
   * rest of the materials since they need to be updated to sum to 1. When adding
   * a volume fraction to a zone, the code adds the allowed amount to the zone.
   * In most cases, this is just the computed material volume, though it can be
   * restricted by amounts of other materials if replacing materials. Once a
   * volume fraction has been added, a list of update materials is iterated and
   * a corresponding amount is subtracted from them. This update list consists of
   * "completely free", all non-shape  materials, and finally the shape material
   * itself. The shape material is added at the end in case the initial volume
   * fraction addition exceeded 1.
   *
   * The "completely free" material (CF) represents any volume fraction in the
   * zones that has not been assigned to any material. We subtract from CF first
   * as part of the update process so we do not have to subtract from real materials
   * in the event that a zone is not full.
   *
   * The list of update materials depends on the shape's material replacement rule.
   *
   * Example:
   *                        update mats
   *                        |---->
   *
   *                  mat1  |  CF     mat0
   *                  ---------------------
   * add 0.4 mat1  -> | 0.4 | 0.2  | 0.4  |
   *                  ---------------------
   *
   *                  mat1     CF     mat0     subtract
   *                  ---------------------    -------
   *                  | 0.8 | 0.2  | 0.4  |    | 0.4 |  (added 0.4 to mat1)
   *                  ---------------------    -------
   *
   *                  mat1     CF     mat0     subtract
   *                  ---------------------    -------
   *                  | 0.8 | 0.0  | 0.4  |    | 0.2 |  (subtracted 0.2 from CF)
   *                  ---------------------    -------
   *
   *                  mat1     CF     mat0     subtract
   *                  ---------------------    -------
   *                  | 0.8 | 0.0  | 0.2  |    | 0.0 |  (subtracted 0.2 from mat0)
   *                  ---------------------    -------
   */
  template <typename ExecSpace>
  void applyReplacementRulesImpl(const klee::Shape& shape)
  {
    // Make sure the material lists are up to date.
    populateMaterials();

    // Get the free material so it is created first.
    // mfem::GridFunction* freeMat = getCompletelyFree<ExecSpace>();
    axom::ArrayView<double> freeMat = getCompletelyFree<ExecSpace>();

    // Get this shape's material, creating the GridFunction if needed.
    auto matVF = getMaterial(shape.getMaterial());
    // int dataSize = matVF.first->Size();
    int dataSize = matVF.first.size();

    // Get this shape's array.
    auto shapeVolFracName = axom::fmt::format("shape_vol_frac_{}", shape.getName());
    // auto* shapeVolFrac = this->getDC()->GetField(shapeVolFracName);
    auto shapeVolFrac = getScalarCellData(shapeVolFracName);
    SLIC_ERROR_IF(shapeVolFrac.empty(),
                  "Field '" + shapeVolFracName +
                    "' must be pre-allocated"
                    " in the Conduit Node computational mesh before using"
                    " IntersectionShaper::applyReplacementRules.");

    // Allocate some memory for the replacement rule data arrays.
    Array<double> vf_subtract_array(dataSize, dataSize, m_allocatorId);
    Array<double> vf_writable_array(dataSize, dataSize, m_allocatorId);
    ArrayView<double> vf_subtract(vf_subtract_array);
    ArrayView<double> vf_writable(vf_writable_array);

    // Determine which grid functions need to be considered for VF updates.
    // std::vector<std::pair<mfem::GridFunction*, int>> gf_order_by_matnumber;
    // std::vector<mfem::GridFunction*> updateVFs, excludeVFs;
    std::vector<std::pair<axom::ArrayView<double>, int>> gf_order_by_matnumber;
    std::vector<axom::ArrayView<double>> updateVFs, excludeVFs;
    if(!shape.getMaterialsReplaced().empty())
    {
      // Include materials replaced in updateVFs.
      for(const auto& name : shape.getMaterialsReplaced())
      {
        gf_order_by_matnumber.emplace_back(getMaterial(name));
      }
    }
    else
    {
      // Include all materials except those in "does_not_replace".
      // We'll also sort them by material number since the field map
      // sorts them by name rather than order added.
      std::vector<std::string> materialNames = getMaterialNames();
      for(auto name : materialNames)
      {
        // Check whether the field name is not the
        // "free" field, which we handle specially)
        if(name != m_free_mat_name)
        {
          // See if the field is in the exclusion list. For the normal
          // case, the list is empty so we'd add the material.
          auto it2 = std::find(shape.getMaterialsNotReplaced().cbegin(),
                               shape.getMaterialsNotReplaced().cend(),
                               name);
          // The field is not in the exclusion list so add it to vfs.
          if(it2 == shape.getMaterialsNotReplaced().cend())
          {
            // Do not add the current shape material since it should
            // not end up in updateVFs.
            if(name != shape.getMaterial())
            {
              gf_order_by_matnumber.emplace_back(getMaterial(name));
            }
          }
          else
          {
            // The material was in the exclusion list. This means that
            // cannot write to materials that have volume fraction in
            // that zone.
            excludeVFs.emplace_back(getMaterial(name).first);
          }
        }
      }
    }
    // Sort eligible update materials by material number.
    std::sort(
      gf_order_by_matnumber.begin(),
      gf_order_by_matnumber.end(),
      [&](const std::pair<axom::ArrayView<double>, int>& lhs,
          const std::pair<axom::ArrayView<double>, int>& rhs) { return lhs.second < rhs.second; });

    // Append the completely free grid function to the materials we update
    // Add it first so it is the highest priority material. This helps us
    // account for how much we need to deduct from the real materials during
    // subtraction.
    updateVFs.push_back(freeMat);

    // Append the grid functions in mat number order.
    for(const auto& mat : gf_order_by_matnumber)
    {
      updateVFs.push_back(mat.first);
    }

    // We put the current shape's material at the end of the update list,
    // in the case that adding to its VF initially exceeds 1.
    updateVFs.push_back(getMaterial(shape.getMaterial()).first);

    // Figure out how much VF in each zone can be written.
    if(!shape.getMaterialsReplaced().empty())
    {
      // Replaces - We'll sum up the VFs that we can replace in a zone.
      AXOM_ANNOTATE_SCOPE("compute_vf_writable");
      axom::for_all<ExecSpace>(dataSize, AXOM_LAMBDA(axom::IndexType i) { vf_writable[i] = 0.; });

      for(const auto& name : shape.getMaterialsReplaced())
      {
        auto mat = getMaterial(name);
        TempArrayView<ExecSpace> matVFView(mat.first, false);
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) {
            vf_writable[i] += matVFView[i];
            vf_writable[i] = (vf_writable[i] > 1.) ? 1. : vf_writable[i];
          });
      }
    }
    else
    {
      // Does not replace. We can replace all except for listed mats.
      AXOM_ANNOTATE_SCOPE("compute_vf_writable");
      axom::for_all<ExecSpace>(dataSize, AXOM_LAMBDA(axom::IndexType i) { vf_writable[i] = 1.; });

      for(auto& gf : excludeVFs)
      {
        TempArrayView<ExecSpace> matVFView(gf, false);
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) {
            vf_writable[i] -= matVFView[i];
            vf_writable[i] = (vf_writable[i] < 0.) ? 0. : vf_writable[i];
          });
      }
    }

    // Compute the volume fractions for the current shape's material.
    {
      AXOM_ANNOTATE_SCOPE("compute_vf");

      TempArrayView<ExecSpace> matVFView(matVF.first, true);
      TempArrayView<ExecSpace> shapeVFView(shapeVolFrac, true);

      axom::ArrayView<double> overlap_volumes_view = m_overlap_volumes.view();
      axom::ArrayView<double> cell_volumes_view = m_cell_volumes.view();

      axom::for_all<ExecSpace>(
        dataSize,
        AXOM_LAMBDA(axom::IndexType i) {
          // Update this material's VF and vf_subtract, which is the
          // amount to subtract from the gf's in updateVF.
          double vf = (overlap_volumes_view[i] / cell_volumes_view[i]);

          // Write at most the writable amount.
          double vf_actual = (vf <= vf_writable[i]) ? vf : vf_writable[i];

          // NOTE: if matVFView[i] temporarily exceeds 1, it will be corrected
          //       during the subtraction stage.
          matVFView[i] += vf_actual;
          vf_subtract[i] = vf_actual;

          // Store the max shape VF.
          shapeVFView[i] = vf;
        });
    }

    // Iterate over updateVFs to subtract off VFs we allocated to the current shape's material.
    {
      AXOM_ANNOTATE_SCOPE("update_vf");
      for(auto& gf : updateVFs)
      {
        TempArrayView<ExecSpace> matVFView(gf, true);
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) {
            constexpr double INSIGNIFICANT_VOLFRAC = 1.e-14;
            double s = (matVFView[i] <= vf_subtract[i]) ? matVFView[i] : vf_subtract[i];
            matVFView[i] -= s;
            // Turn any slight negatives or positive insignificant volume fractions to zero.
            matVFView[i] = (matVFView[i] < INSIGNIFICANT_VOLFRAC) ? 0. : matVFView[i];
            vf_subtract[i] -= s;
          });
      }
    }
  }

  /*!
   * \brief Apply material replacement rules for the current shape, using
   *        the appropriate execution policy.
   */
  void applyReplacementRules(const klee::Shape& shape) override
  {
    AXOM_ANNOTATE_SCOPE("applyReplacementRules");

    switch(m_execPolicy)
    {
    case RuntimePolicy::seq:
      applyReplacementRulesImpl<seq_exec>(shape);
      break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    case RuntimePolicy::omp:
      applyReplacementRulesImpl<omp_exec>(shape);
      break;
#endif  // AXOM_USE_OPENMP
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    case RuntimePolicy::cuda:
      applyReplacementRulesImpl<cuda_exec>(shape);
      break;
#endif  // AXOM_USE_CUDA
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    case RuntimePolicy::hip:
      applyReplacementRulesImpl<hip_exec>(shape);
      break;
#endif  // AXOM_USE_HIP
    }
    AXOM_UNUSED_VAR(shape);
  }

  void finalizeShapeQuery() override
  {
    AXOM_ANNOTATE_SCOPE("finalizeShapeQuery");

    // Implementation here -- destroy BVH tree and other shape-based data structures
    m_surfaceMesh.reset();
  }

  ///@}

public:
  /*!
    \brief Prepares the shaping query, based on the policy member set
    (default is sequential)

    \internal This method populates m_tets or m_octs from the given \c
    shape.  These arrays are used in runShapeQuery.
  */
  void prepareShapeQuery(klee::Dimensions shapeDimension, const klee::Shape& shape) override
  {
    AXOM_ANNOTATE_SCOPE("prepareShapeQuery");
    const std::string shapeFormat = shape.getGeometry().getFormat();

    // Save m_percentError and m_level in case refineShape needs to change them
    // to meet the overall desired error tolerance for the volume.
    const double saved_percentError = m_percentError;
    const double saved_level = m_level;
    if(shapeFormat == "c2c")
    {
      // Refine the shape, potentially reloading it more refined.
      refineShape(shape);
    }

    // Now that the mesh is refined, dispatch to device implementations.
    switch(m_execPolicy)
    {
    case RuntimePolicy::seq:
      prepareShapeQueryImpl<seq_exec>(shapeDimension, shape);
      break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    case RuntimePolicy::omp:
      prepareShapeQueryImpl<omp_exec>(shapeDimension, shape);
      break;
#endif  // AXOM_USE_OPENMP
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    case RuntimePolicy::cuda:
      prepareShapeQueryImpl<cuda_exec>(shapeDimension, shape);
      break;
#endif  // AXOM_USE_CUDA
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    case RuntimePolicy::hip:
      prepareShapeQueryImpl<hip_exec>(shapeDimension, shape);
      break;
#endif  // AXOM_USE_HIP
    }
    AXOM_UNUSED_VAR(shapeDimension);
    AXOM_UNUSED_VAR(shape);

    // Restore m_percentError, m_level in case refineShape changed them.
    m_percentError = saved_percentError;
    m_level = saved_level;
  }

  // Runs the shaping query, based on the policy member and shape format set
  // (default is sequential)
  // Fills m_overlap_volumes and m_cell_volumes, whose data
  // will be in the default memory space for m_execPolicy.
  // The data will be used in applyReplacementRules and can
  // also be accessed by getOverlapVolumes() and getCellVolumes().
  void runShapeQuery(const klee::Shape& shape) override
  {
    AXOM_ANNOTATE_SCOPE("runShapeQuery");
    const std::string shapeFormat = shape.getGeometry().getFormat();

    // C2C mesh is not discretized into tets, but all others are.
    if(surfaceMeshIsTet())
    {
      switch(m_execPolicy)
      {
      case RuntimePolicy::seq:
        runShapeQuery3DImpl<seq_exec, TetrahedronType>(shape, m_tets, m_tetcount);
        break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case RuntimePolicy::omp:
        runShapeQuery3DImpl<omp_exec, TetrahedronType>(shape, m_tets, m_tetcount);
        break;
#endif  // AXOM_USE_OPENMP
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case RuntimePolicy::cuda:
        runShapeQuery3DImpl<cuda_exec, TetrahedronType>(shape, m_tets, m_tetcount);
        break;
#endif  // AXOM_USE_CUDA
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case RuntimePolicy::hip:
        runShapeQuery3DImpl<hip_exec, TetrahedronType>(shape, m_tets, m_tetcount);
        break;
#endif  // AXOM_USE_HIP
      }
    }
    else if(shapeFormat == "c2c" && getCompMeshDim() == 3)
    {
      switch(m_execPolicy)
      {
      case RuntimePolicy::seq:
        runShapeQuery3DImpl<seq_exec, OctahedronType>(shape, m_octs, m_octcount);
        break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case RuntimePolicy::omp:
        runShapeQuery3DImpl<omp_exec, OctahedronType>(shape, m_octs, m_octcount);
        break;
#endif  // AXOM_USE_OPENMP
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case RuntimePolicy::cuda:
        runShapeQuery3DImpl<cuda_exec, OctahedronType>(shape, m_octs, m_octcount);
        break;
#endif  // AXOM_USE_CUDA
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case RuntimePolicy::hip:
        runShapeQuery3DImpl<hip_exec, OctahedronType>(shape, m_octs, m_octcount);
        break;
#endif  // AXOM_USE_HIP
      }
    }
    else if(shapeFormat == "c2c" && getCompMeshDim() == 2)
    {
      switch(m_execPolicy)
      {
      case RuntimePolicy::seq:
        runShapeQuery2DImpl<seq_exec>(shape, m_tris, m_tricount);
        break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case RuntimePolicy::omp:
        runShapeQuery2DImpl<omp_exec>(shape, m_tris, m_tricount);
        break;
#endif  // AXOM_USE_OPENMP
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case RuntimePolicy::cuda:
        runShapeQuery2DImpl<cuda_exec>(shape, m_tris, m_tricount);
        break;
#endif  // AXOM_USE_CUDA
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case RuntimePolicy::hip:
        runShapeQuery2DImpl<hip_exec>(shape, m_tris, m_tricount);
        break;
#endif  // AXOM_USE_HIP
      }
    }
    else if(shapeFormat == "stl" || surfaceMeshIsTri())
    {
      switch(m_execPolicy)
      {
      case RuntimePolicy::seq:
        runShapeQuery2DImpl<seq_exec>(shape, m_tris, m_tricount);
        break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case RuntimePolicy::omp:
        runShapeQuery2DImpl<omp_exec>(shape, m_tris, m_tricount);
        break;
#endif  // AXOM_USE_OPENMP
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case RuntimePolicy::cuda:
        runShapeQuery2DImpl<cuda_exec>(shape, m_tris, m_tricount);
        break;
#endif  // AXOM_USE_CUDA
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case RuntimePolicy::hip:
        runShapeQuery2DImpl<hip_exec>(shape, m_tris, m_tricount);
        break;
#endif  // AXOM_USE_HIP
      }
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("The shape format {} is unsupported", shapeFormat));
    }
  }

  axom::ArrayView<const double> getOverlapVolumes() const { return m_overlap_volumes.view(); }

  axom::ArrayView<const double> getCellVolumes() const { return m_cell_volumes.view(); }

  double sumOverlapVolumes(bool global = true) const
  {
    double overlapVol = 0.0;
    switch(m_execPolicy)
    {
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    case RuntimePolicy::omp:
      overlapVol = sumArray<omp_exec>(m_overlap_volumes.data(), m_overlap_volumes.size());
      break;
#endif  // AXOM_USE_OPENMP
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    case RuntimePolicy::cuda:
      overlapVol = sumArray<cuda_exec>(m_overlap_volumes.data(), m_overlap_volumes.size());
      break;
#endif  // AXOM_USE_CUDA
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    case RuntimePolicy::hip:
      overlapVol = sumArray<hip_exec>(m_overlap_volumes.data(), m_overlap_volumes.size());
      break;
#endif  // AXOM_USE_HIP
    case RuntimePolicy::seq:
    default:
      overlapVol = sumArray<seq_exec>(m_overlap_volumes.data(), m_overlap_volumes.size());
      break;
    }

    if(global)
    {
      overlapVol = this->allReduceSum(overlapVol);
    }
    return overlapVol;
  }

  template <typename ExecSpace, typename Summable>
  Summable sumArray(const Summable* a, axom::IndexType count) const
  {
    axom::ReduceSum<ExecSpace, Summable> vsum {0};
    axom::for_all<ExecSpace>(count, AXOM_LAMBDA(axom::IndexType i) { vsum += a[i]; });
    Summable sum = static_cast<Summable>(vsum.get());
    return sum;
  }

  void adjustVolumeFractions() override
  {
    // Implementation here -- not sure if this will require anything for intersection-based shaping
  }

  std::vector<std::string> getMaterialNames() const
  {
    std::vector<std::string> materialNames;
#if defined(AXOM_USE_MFEM)
    if(m_dc)
    {
      for(auto it : this->getDC()->GetFieldMap())
      {
        std::string materialName = fieldNameToMaterialName(it.first);
        if(!materialName.empty())
        {
          materialNames.emplace_back(materialName);
        }
      }
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bpGrp)
    {
      auto fieldsGrp = m_bpGrp->getGroup("fields");
      SLIC_ERROR_IF(fieldsGrp == nullptr, "Input blueprint mesh lacks the 'fields' Group/Node.");
      for(auto& group : fieldsGrp->groups())
      {
        std::string materialName = fieldNameToMaterialName(group.getName());
        if(!materialName.empty())
        {
          materialNames.emplace_back(materialName);
        }
      }
    }
#endif
    return materialNames;
  }

  /*!
   * \brief Gets the grid function and material number for a material name.
   *
   * \param materialName The name of the material.
   *
   * \return A pair containing the associated grid function (an element
   *         in the m_vf_grid_functions array) and material
   *         number (its index in the array).
   */
  std::pair<axom::ArrayView<double>, int> getMaterial(const std::string& materialName)
  {
    // If we already know about the material, return it.
    for(size_t i = 0; i < m_vf_material_names.size(); i++)
    {
      if(m_vf_material_names[i] == materialName)
      {
        return std::make_pair(m_vf_grid_functions[i], static_cast<int>(i));
      }
    }

    // Get or create the volume fraction field for this shape's material
    auto materialVolFracName = materialNameToFieldName(materialName);

    bool makeNewData = !hasData(materialVolFracName);
    auto matVolFrac = getScalarCellData(materialVolFracName);
    SLIC_ASSERT(!matVolFrac.empty());
    if(makeNewData)
    {
      // Zero out the volume fractions (on host).
#ifdef AXOM_USE_UMPIRE
      auto allocId = matVolFrac.getAllocatorID();
      const axom::MemorySpace memorySpace = axom::detail::getAllocatorSpace(allocId);
      const bool onDevice =
        memorySpace == axom::MemorySpace::Device || memorySpace == axom::MemorySpace::Unified;
#else
      const bool onDevice = false;
#endif
      if(onDevice)
      {
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
        if(m_execPolicy == RuntimePolicy::cuda)
        {
          axom::for_all<axom::CUDA_EXEC<256>>(
            matVolFrac.size(),
            AXOM_LAMBDA(axom::IndexType i) { matVolFrac[i] = 0.0; });
        }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
        if(m_execPolicy == RuntimePolicy::hip)
        {
          axom::for_all<axom::HIP_EXEC<256>>(
            matVolFrac.size(),
            AXOM_LAMBDA(axom::IndexType i) { matVolFrac[i] = 0.0; });
        }
#endif
      }
      else
      {
        memset(matVolFrac.data(), 0, matVolFrac.size() * sizeof(double));
      }
    }

    // Add the material to our vectors.
    int idx = static_cast<int>(m_vf_grid_functions.size());
    m_vf_grid_functions.push_back(matVolFrac);
    m_vf_material_names.push_back(materialName);

    return std::make_pair(matVolFrac, idx);
  }

private:
  /*!
   * \brief Filter the input mesh so it does not have any degenerate segments
   *        (consecutive points that are the same) and the points are increasing
   *        order on the Z (really X axis).
   *
   * \param m The input mesh.
   * \return A new mint::Mesh that has been cleaned up.
   */
  SegmentMesh* filterMesh(const SegmentMesh* m) const
  {
    // Number of points in polyline
    int pointcount = m->getNumberOfNodes();

    // We'll be filtering the points in the mesh. Make a new mesh to contain
    // those points.
    SegmentMesh* newm = new SegmentMesh(m->getDimension(), mint::SEGMENT);
    newm->reserveNodes(pointcount);

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format(" Checking contour with {} points for degenerate segments ", pointcount)));

    constexpr int R = 1;
    constexpr int Z = 0;

    // Add contour points
    int polyline_size = 0;
    const double EPS = m_vertexWeldThreshold;
    const double EPS_SQ = EPS * EPS;
    Point2D cur_point, prev_point;
    for(int i = 0; i < pointcount; ++i)
    {
      // Get the current point.
      m->getNode(i, cur_point.data());

      // Check for degenerate segments
      if(polyline_size > 0)
      {
        using axom::utilities::isNearlyEqual;

        // both r are (nearly) equal to 0
        if(isNearlyEqual(cur_point[R], 0.0, EPS) && isNearlyEqual(prev_point[R], 0.0, EPS))
        {
          continue;
        }

        // both Z are (nearly) the same
        if(isNearlyEqual(cur_point[Z], prev_point[Z], EPS))
        {
          continue;
        }

        // both points are (nearly) the same
        if(squared_distance(cur_point, prev_point) < EPS_SQ)
        {
          continue;
        }
      }

      // Add the current point to the new mesh.
      newm->appendNode(cur_point[Z], cur_point[R]);
      prev_point = cur_point;
      polyline_size += 1;
    }

    // Get the ZR coordinates.
    double* z = newm->getCoordinateArray(mint::X_COORDINATE);
    double* r = newm->getCoordinateArray(mint::Y_COORDINATE);

    // Check if we need to flip the points order.
    // discretize() is only valid if x increases as index increases
    bool flip = false;
    if(polyline_size > 1)
    {
      if(z[0] > z[1])
      {
        flip = true;
        SLIC_INFO("Order of contour points has been reversed!"
                  << " Discretization algorithm expects Z values of contour "
                     "points to be increasing");
      }
    }

    SLIC_INFO(axom::fmt::format("{:-^80}",
                                axom::fmt::format(axom::utilities::locale(),
                                                  " Discretizing contour with {:L} points ",
                                                  polyline_size)));

    // Flip point order
    if(flip)
    {
      int i = polyline_size - 1;
      int j = 0;
      while(i > j)
      {
        axom::utilities::swap(z[i], z[j]);
        axom::utilities::swap(r[i], r[j]);
        i -= 1;
        j += 1;
      }
    }

    // Now make polyline_size - 1 segments.
    int numNewSegments = polyline_size - 1;
    newm->reserveCells(numNewSegments);
    for(int i = 0; i < numNewSegments; ++i)
    {
      IndexType seg[2] = {i, i + 1};
      newm->appendCell(seg, mint::SEGMENT);
    }

    return newm;
  }

  /*!
   * \brief Compute the area of the polygon given by \a pts.
   * \param pts The list of points that make up the polygon.
   * \return The polygon area.
   */
  double area(const std::vector<Point2D>& pts) const
  {
    auto npts = static_cast<int>(pts.size());
    // Compute the areas
    double A = 0.;
    for(int i = 0; i < npts; i++)
    {
      int nexti = (i + 1) % npts;
      A += pts[i][0] * pts[nexti][1] - pts[i][1] * pts[nexti][0];
    }
    // Take absolute value just in case the polygon had reverse orientation.
    return fabs(A * 0.5);
  }

  /*!
   * \brief Compute the circle area of a circle of radius at \a level.
   *        The area is approximated by a set of line segments around
   *        the perimeter of the circle where the number of line segments
   *        is determined by \a level.
   *
   * \param radius The radius of the circle.
   * \return The circle area
   */
  double calcCircleArea(double radius, int level) const
  {
    int npts = 3 * pow(2, level);
    std::vector<Point2D> pts;
    pts.reserve(npts);
    for(int i = 0; i < npts; i++)
    {
      double angle = 2. * M_PI * static_cast<double>(i) / static_cast<double>(npts);
      pts.push_back(Point2D {radius * cos(angle), radius * sin(angle)});
    }
    return area(pts);
  }

  /*!
   * \brief Compute the area of a circle of radius at \a level.
   * \param radius The radius of the circle.
   * \param level The refinement level for the circle.
   * \return The circle area
   * \note If the requested circle area is not in the table, it will be
   *       computed but that gets SLOW.
   */
  double circleArea(double radius, int level) const
  {
    // The lut covers most values we'd use.
    static const double lut[] = {
      /*level 0*/ 1.29903810568,   // diff=1.84255454791
      /*level 1*/ 2.59807621135,   // diff=0.543516442236
      /*level 2*/ 3,               // diff=0.14159265359
      /*level 3*/ 3.10582854123,   // diff=0.0357641123595
      /*level 4*/ 3.13262861328,   // diff=0.00896404030856
      /*level 5*/ 3.13935020305,   // diff=0.00224245054292
      /*level 6*/ 3.14103195089,   // diff=0.000560702699288
      /*level 7*/ 3.14145247229,   // diff=0.000140181304342
      /*level 8*/ 3.14155760791,   // diff=3.50456779015e-05
      /*level 9*/ 3.14158389215,   // diff=8.76144145456e-06
      /*level 10*/ 3.14159046323,  // diff=2.19036162807e-06
      /*level 11*/ 3.141592106,    // diff=5.47590411681e-07
      /*level 12*/ 3.14159251669,  // diff=1.36898179903e-07
      /*level 13*/ 3.14159261936,  // diff=3.4225411838e-08
      /*level 14*/ 3.14159264503,  // diff=8.55645243547e-09
    };
    constexpr int MAX_LEVELS = sizeof(lut) / sizeof(double);
    return (level < MAX_LEVELS) ? (lut[level] * radius * radius) : calcCircleArea(radius, level);
  }

  /*!
   * \brief Iterate over line segments in the surface mesh and compute
   *        truncated cone values with circles at given refinement and
   *        add them all up. This will be the approximate revolved
   *        volume that is possible with the current curve to line segment
   *        refinement that was done.
   *
   * \param m The mint mesh that contains the line segments.
   * \param level The circle refinement level.
   *
   * \note This function assumes that the line segments have been
   *       processed and increase in x.
   *
   * \return The approximate revolved volume for the linearized curve
   *         and circle refinement level.
   */
  double volume(const mint::Mesh* m, int level) const
  {
    const int numSurfaceVertices = m->getNumberOfNodes();
    const int nSegments = numSurfaceVertices - 1;
    const double* z = m->getCoordinateArray(mint::X_COORDINATE);
    const double* r = m->getCoordinateArray(mint::Y_COORDINATE);
    double vol_approx = 0.;
    for(int seg = 0; seg < nSegments; seg++)
    {
      double r0 = r[seg];
      double r1 = r[seg + 1];
      double h = fabs(z[seg + 1] - z[seg]);
      if(axom::utilities::isNearlyEqual(r0, r1, m_vertexWeldThreshold))
      {
        // cylinder
        double A = circleArea(r0, level);

        // Add the cylinder to the volume total.
        vol_approx += A * h;
      }
      else
      {
        // truncated cone
        double h2 = (r0 * h / (r0 - r1)) - h;

        // Approximate cone volume
        double A2 = circleArea(r0, level);
        double A3 = circleArea(r1, level);
        double approx_cone_vol = (1. / 3.) * (A2 * h + (A2 - A3) * h2);

        // Add the truncated cone to the volume total.
        vol_approx += approx_cone_vol;
      }
    }

    return vol_approx;
  }

  /*!
   * \brief Refine the shape to get under a certain errorPercent in the
   *        linearized revolved volume. To do this, we may end up reloading
   *        the shape.
   *
   * \param shape The shape to refine.
   *
   * \note The m_surfaceMesh will be set to a mesh that should be fine
   *       enough that its linearized representation is close enough to
   *       the specified error percent. m_level will be set to the circle
   *       refinement level that let us meet the percent error target.
   */
  void refineShape(const klee::Shape& shape)
  {
    // If we are not refining dynamically, return.
    if(m_percentError <= MINIMUM_PERCENT_ERROR || m_refinementType != DiscreteShape::RefinementDynamic)
    {
      return;
    }

    // If the prior loadShape call was unable to create a revolved volume for
    // the shape then we can't do any better than the current mesh.
    if(m_revolvedVolume <= DEFAULT_REVOLVED_VOLUME)
    {
      return;
    }

    /*!
     * \brief Examines the history values to determine if the deltas between
     *        iterations are sufficiently small that the iterations should
     *        terminate, even though it might not have reached the desired
     *        target error.
     *
     * \param iteration The solve iteration.
     * \param history A circular buffer of the last several history values.
     * \param nhistory The number of history values.
     * \param percentError The percent error to achieve.
     *
     * \note This function assumes that history values increase.
     */
    auto diminishing_returns =
      [](int iteration, const double* history, int nhistory, double percentError) -> bool {
      bool dr = false;
      // We have enough history to decide if there are diminishing returns.
      if(iteration >= nhistory)
      {
        // Compute a series of percents between pairs of history values.
        double sum = 0.;
        for(int i = 0; i < nhistory - 1; i++)
        {
          int cur = (iteration - i) % nhistory;
          int prev = (iteration - i - 1) % nhistory;
          double pct = 100. * (1. - history[prev] / history[cur]);
          sum += fabs(pct);
        }
        double avg_pct = sum / static_cast<double>(nhistory - 1);
        // Check whether the pct is less than the percentError. It may be
        // good enough.
        dr = avg_pct < percentError;
        if(dr)
        {
          SLIC_INFO(axom::fmt::format("Dimishing returns triggered: {} < {}.", avg_pct, percentError));
        }
      }
      return dr;
    };

    // We have gone through loadShape once at this point. For C2C contours,
    // we should have a reasonable revolved volume of the shape. We can
    // compute the volume of the shape using the discretized mesh and see
    // how close we are to m_percentError. If we can't get there by increasing
    // the circle refinement, we can loadShape again using a smaller error
    // tolerance.

    // Initial values for the refinement knobs.
    double curvePercentError = m_percentError;
    int circleLevel = m_level;

    // Save the revolvedVolume
    double revolvedVolume = m_revolvedVolume;

    // Try refining the curve different ways to see if we get to a refinement
    // strategy that meets the error tolerance.
    bool refine = true;
    // Limit level refinement for now since it makes too many octahedra.
    int MAX_LEVELS = m_level + 1;
    constexpr int MAX_ITERATIONS = 20;
    constexpr int MAX_HISTORY = 4;
    constexpr double MINIMUM_CURVE_PERCENT_ERROR = 1.e-10;
    constexpr double CURVE_PERCENT_SCALING = 0.5;
    double history[MAX_HISTORY] = {0., 0., 0., 0.};
    for(int iteration = 0; iteration < MAX_ITERATIONS && refine; iteration++)
    {
      // Increase the circle refinement level to see if we can get to an acceptable
      // error percentage using the line segment mesh's revolved volume.
      double pct = 0., currentVol = 0.;
      for(int level = m_level; level < MAX_LEVELS; level++)
      {
        // Compute the revolved volume of the surface mesh at level.
        currentVol = volume(m_surfaceMesh.get(), level);
        pct = 100. * (1. - currentVol / revolvedVolume);

        SLIC_INFO(
          axom::fmt::format("Refining... "
                            "revolvedVolume = {}"
                            ", currentVol = {}"
                            ", pct = {}"
                            ", level = {}"
                            ", curvePercentError = {}",
                            revolvedVolume,
                            currentVol,
                            pct,
                            level,
                            curvePercentError));

        if(pct <= m_percentError)
        {
          SLIC_INFO(
            axom::fmt::format("Contour refinement complete. "
                              "revolvedVolume = {}"
                              ", currentVol = {}"
                              ", pct = {}"
                              ", level = {}"
                              ", curvePercentError = {}",
                              revolvedVolume,
                              currentVol,
                              pct,
                              level,
                              curvePercentError));

          circleLevel = level;
          refine = false;
          break;
        }
      }

      // Check whether there are diminishing returns. In other words, no level
      // of refinement can quite get to the target error percent.
      history[iteration % MAX_HISTORY] = currentVol;
      if(diminishing_returns(iteration - 1, history, MAX_HISTORY, m_percentError))
      {
        circleLevel = MAX_LEVELS - 1;
        refine = false;

        SLIC_INFO(
          axom::fmt::format("Stop refining due to diminishing returns. "
                            "revolvedVolume = {}"
                            ", currentVol = {}"
                            ", pct = {}"
                            ", level = {}"
                            ", curvePercentError = {}",
                            revolvedVolume,
                            currentVol,
                            pct,
                            circleLevel,
                            curvePercentError));

        // NOTE: Trying to increase circleLevel at this point does not help.
      }

      if(refine)
      {
        // We could not get to an acceptable error percentage and we should refine.

        // Check whether to permit the curve to further refine.
        double ce = curvePercentError * CURVE_PERCENT_SCALING;
        if(ce > MINIMUM_CURVE_PERCENT_ERROR)
        {
          // Set the new curve error.
          curvePercentError = ce;

          // Free the previous surface mesh.
          m_surfaceMesh.reset();

          // Reload the shape using new curvePercentError. This will cause
          // a new m_surfaceMesh to be created.
          double rv = 0.;
          SLIC_INFO(axom::fmt::format("Reloading shape {} with curvePercentError = {}.",
                                      shape.getName(),
                                      curvePercentError));
          loadShapeInternal(shape, curvePercentError, rv);

          // Filter the mesh, store in m_surfaceMesh.
          SegmentMesh* newm = filterMesh(dynamic_cast<const SegmentMesh*>(m_surfaceMesh.get()));
          m_surfaceMesh.reset(newm);
        }
        else
        {
          SLIC_INFO(
            axom::fmt::format("Stopping refinement due to curvePercentError {} being too small.", ce));
          refine = false;
        }
      }
    }

    // We arrived at a mesh refinement that satisfied m_percentError.
    m_level = circleLevel;
  }

  /// Whether the given field exists in the mesh.
  bool hasData(const std::string& fieldName)
  {
    bool has = false;
#if defined(AXOM_USE_MFEM)
    if(m_dc != nullptr)
    {
      has = m_dc->HasField(fieldName);
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bpGrp != nullptr)
    {
      std::string fieldPath = axom::fmt::format("fields/{}", fieldName);
      has = m_bpGrp->hasGroup(fieldPath);
    }
#endif
    return has;
  }

  /*!  \brief Get a scalar double-type field data from the mesh,
    "fields/fieldName/values", creating it if it doesn't exist.

    Also add the corresponding entry in the blueprint field
    "matsets/fieldName/volume_fractions".

    If mesh is in an external Conduit Node, and the field
    doesn't exist, emit a warning and return an empty ArrayView.
    We don't add fields to the external Conduit Node.
    \see Shaper::Shaper().
  */
  axom::ArrayView<double> getScalarCellData(const std::string& fieldName, bool volumeDependent = false)
  {
    axom::ArrayView<double> rval;

#if defined(AXOM_USE_MFEM)
    if(m_dc != nullptr)
    {
      mfem::GridFunction* gridFunc = nullptr;
      if(m_dc->HasField(fieldName))
      {
        gridFunc = m_dc->GetField(fieldName);
      }
      else
      {
        gridFunc = newVolFracGridFunction();
        m_dc->RegisterField(fieldName, gridFunc);
      }
      rval = axom::ArrayView<double>(gridFunc->GetData(), gridFunc->Size());
    }
#endif

#if defined(AXOM_USE_CONDUIT)
    if(m_bpGrp != nullptr)
    {
      std::string fieldPath = "fields/" + fieldName;
      auto dtype = conduit::DataType::float64(m_cellCount);
      axom::sidre::View* valuesView = nullptr;
      if(m_bpGrp->hasGroup(fieldPath))
      {
        auto* fieldGrp = m_bpGrp->getGroup(fieldPath);
        valuesView = fieldGrp->getView("values");
        SLIC_ASSERT(fieldGrp->getView("association")->getString() == std::string("element"));
        SLIC_ASSERT(fieldGrp->getView("topology")->getString() == m_bpTopo);
        SLIC_ASSERT(valuesView->getNumElements() == m_cellCount);
        SLIC_ASSERT(valuesView->getNode().dtype().id() == dtype.id());
      }
      else
      {
        if(m_bpNodeExt != nullptr)
        {
          /*
            If the computational mesh is an external conduit::Node, it
            must have all necessary fields.  We will only generate
            fields for meshes in sidre::Group, where the user can set
            the allocator id for only array data.  conduit::Node doesn't
            have this capability.
          */
          SLIC_WARNING_IF(m_bpNodeExt != nullptr,
                          "For a computational mesh in a conduit::Node, all"
                          " output fields must be preallocated before shaping."
                          "  IntersectionShaper will NOT contravene the user's"
                          " memory management.  The cell-centered field '" +
                            fieldPath +
                            "' is missing.  Please pre-allocate"
                            " this output memory, or to have IntersectionShaper"
                            " allocate it, construct the IntersectionShaper"
                            " with the mesh as a sidre::Group  with your"
                            " specific allocator id.");
        }
        else
        {
          constexpr axom::IndexType componentCount = 1;
          axom::IndexType shape[2] = {m_cellCount, componentCount};
          auto* fieldGrp = m_bpGrp->createGroup(fieldPath);
          // valuesView = fieldGrp->createView("values");
          valuesView =
            fieldGrp->createViewWithShape("values", axom::sidre::DataTypeId::FLOAT64_ID, 2, shape);
          fieldGrp->createView("association")->setString("element");
          fieldGrp->createView("topology")->setString(m_bpTopo);
          fieldGrp->createView("volume_dependent")
            ->setString(std::string(volumeDependent ? "true" : "false"));
          valuesView->allocate();
        }
      }

      rval = axom::ArrayView<double>(static_cast<double*>(valuesView->getVoidPtr()), m_cellCount);
    }
#endif
    return rval;
  }

#if defined(__CUDACC__)
public:
    // These methods should be private, but NVCC complains unless they're public.
#endif

  template <typename ExecSpace>
  void populateQuadsFromMesh()
  {
    AXOM_ANNOTATE_SCOPE("populateQuadsFromMesh");
    constexpr int NUM_VERTS_PER_QUAD = 4;
    constexpr int NUM_COMPS_PER_VERT = 2;
    const int allocId = m_allocatorId;

    axom::Array<double> vertCoords(m_cellCount * NUM_VERTS_PER_QUAD * NUM_COMPS_PER_VERT,
                                   m_cellCount * NUM_VERTS_PER_QUAD * NUM_COMPS_PER_VERT,
                                   allocId);

#if defined(AXOM_USE_MFEM)
    if(m_dc != nullptr)
    {
      populateVertCoordsFromMFEMMesh<ExecSpace>(vertCoords, 2);
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bpGrp != nullptr)
    {
      populateVertCoordsFromBlueprintMesh2D<ExecSpace>(vertCoords);
    }
#endif

    auto vertCoords_device_view = vertCoords.view();

    // Initialize quad elements
    m_quads = axom::Array<PolygonStaticType>(m_cellCount, m_cellCount, m_allocatorId);
    axom::ArrayView<PolygonStaticType> quads_device_view = m_quads.view();

    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        // Set each quad element vertices
        quads_device_view[i] = PolygonStaticType();
        for(int j = 0; j < NUM_VERTS_PER_QUAD; ++j)
        {
          int vertIndex = (i * NUM_VERTS_PER_QUAD * NUM_COMPS_PER_VERT) + j * NUM_COMPS_PER_VERT;
          quads_device_view[i].addVertex(
            Point2D({vertCoords_device_view[vertIndex], vertCoords_device_view[vertIndex + 1]}));
        }
      });
  }  // end of populateQuadsFromMesh()

  template <typename ExecSpace>
  void populateHexesFromMesh()
  {
    AXOM_ANNOTATE_SCOPE("populateHexesFromMesh");

    constexpr int NUM_VERTS_PER_HEX = 8;
    constexpr int NUM_COMPS_PER_VERT = 3;
    const int allocId = m_allocatorId;

    axom::Array<double> vertCoords(m_cellCount * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT,
                                   m_cellCount * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT,
                                   allocId);

#if defined(AXOM_USE_MFEM)
    if(m_dc != nullptr)
    {
      populateVertCoordsFromMFEMMesh<ExecSpace>(vertCoords, 3);
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bpGrp != nullptr)
    {
      populateVertCoordsFromBlueprintMesh3D<ExecSpace>(vertCoords);
    }
#endif

    auto vertCoords_device_view = vertCoords.view();

    m_hexes = axom::Array<HexahedronType>(m_cellCount, m_cellCount, allocId);
    axom::ArrayView<HexahedronType> hexes_device_view = m_hexes.view();
    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        // Set each hexahedral element vertices
        hexes_device_view[i] = HexahedronType();
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          int vertIndex = (i * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT) + j * NUM_COMPS_PER_VERT;
          hexes_device_view[i][j] = Point3D({vertCoords_device_view[vertIndex],
                                             vertCoords_device_view[vertIndex + 1],
                                             vertCoords_device_view[vertIndex + 2]});
        }
      });  // end of loop to initialize hexahedral elements and bounding boxes
  }

#if defined(AXOM_USE_CONDUIT)
  template <typename ExecSpace>
  void populateVertCoordsFromBlueprintMesh2D(axom::Array<double>& vertCoords)
  {
    using XS = axom::execution_space<ExecSpace>;

    const int allocId = m_allocatorId;

    // Initialize vertices from blueprint mesh and
    // set each shape volume fraction to 1
    // Allocation size is:
    // # of elements * # of vertices per quad * # of components per vertex
    constexpr int NUM_VERTS_PER_QUAD = 4;
    constexpr int NUM_COMPS_PER_VERT = 2;

    // Put mesh in Node so we can use conduit::blueprint utilities.
    // conduit::Node meshNode;
    // m_bpGrp->createNativeLayout(m_bpNodeInt);

    const conduit::Node& topoNode = m_bpNodeInt.fetch_existing("topologies").fetch_existing(m_bpTopo);
    const std::string coordsetName = topoNode.fetch_existing("coordset").as_string();

    // Assume unstructured and hexahedral
    SLIC_ERROR_IF(topoNode["type"].as_string() != "unstructured",
                  "topology type must be 'unstructured'");
    SLIC_ERROR_IF(topoNode["elements/shape"].as_string() != "quad", "element shape must be 'quad'");

    const auto& connNode = topoNode["elements/connectivity"];
    SLIC_ERROR_IF(
      !XS::usesAllocId(axom::getAllocatorIDFromPointer(connNode.data_ptr())),
      std::string(XS::name()) +
        axom::fmt::format(" execution space cannot use the connectivity allocator id {}",
                          axom::getAllocatorIDFromPointer(connNode.data_ptr())));
    SLIC_ERROR_IF(connNode.dtype().id() != conduitDataIdOfAxomIndexType,
                  "IntersectionShaper error: connectivity data type must be "
                  "axom::IndexType.");
    const auto* connPtr = static_cast<const axom::IndexType*>(connNode.data_ptr());
    axom::ArrayView<const axom::IndexType, 2> conn(connPtr, m_cellCount, NUM_VERTS_PER_QUAD);

    const conduit::Node& coordNode = m_bpNodeInt["coordsets"][coordsetName];
    const conduit::Node& coordValues = coordNode.fetch_existing("values");
    axom::IndexType vertexCount = coordValues["x"].dtype().number_of_elements();
    bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
    int stride = isInterleaved ? NUM_COMPS_PER_VERT : 1;

    axom::StackArray<axom::ArrayView<const double>, 2> coordArrays {
      axom::ArrayView<const double>(coordValues["x"].as_double_ptr(), {vertexCount}, stride),
      axom::ArrayView<const double>(coordValues["y"].as_double_ptr(), {vertexCount}, stride)};

    vertCoords = axom::Array<double>(m_cellCount * NUM_VERTS_PER_QUAD * NUM_COMPS_PER_VERT,
                                     m_cellCount * NUM_VERTS_PER_QUAD * NUM_COMPS_PER_VERT,
                                     allocId);
    auto vertCoordsView = vertCoords.view();

    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        // Get the indices of this element's vertices
        auto quadVerts = conn[i];

        // Get the coordinates for the vertices
        for(int j = 0; j < NUM_VERTS_PER_QUAD; ++j)
        {
          auto vertId = quadVerts[j];
          for(int k = 0; k < NUM_COMPS_PER_VERT; k++)
          {
            vertCoordsView[(i * NUM_VERTS_PER_QUAD * NUM_COMPS_PER_VERT) + (j * NUM_COMPS_PER_VERT) + k] =
              coordArrays[k][vertId];
          }
        }
      });
  }

  template <typename ExecSpace>
  void populateVertCoordsFromBlueprintMesh3D(axom::Array<double>& vertCoords)
  {
    using XS = axom::execution_space<ExecSpace>;

    const int allocId = m_allocatorId;

    // Initialize vertices from blueprint mesh and
    // set each shape volume fraction to 1
    // Allocation size is:
    // # of elements * # of vertices per hex * # of components per vertex
    constexpr int NUM_VERTS_PER_HEX = 8;
    constexpr int NUM_COMPS_PER_VERT = 3;

    // Put mesh in Node so we can use conduit::blueprint utilities.
    // conduit::Node meshNode;
    // m_bpGrp->createNativeLayout(m_bpNodeInt);

    const conduit::Node& topoNode = m_bpNodeInt.fetch_existing("topologies").fetch_existing(m_bpTopo);
    const conduit::Node& topoCoordsetNode = topoNode.fetch_existing("coordset");
    const std::string coordsetName = topoCoordsetNode.as_string();

    // Assume unstructured and hexahedral
    SLIC_ERROR_IF(topoNode["type"].as_string() != "unstructured",
                  "topology type must be 'unstructured'");
    SLIC_ERROR_IF(topoNode["elements/shape"].as_string() != "hex", "element shape must be 'hex'");

    const auto& connNode = topoNode["elements/connectivity"];
    SLIC_ERROR_IF(
      !XS::usesAllocId(axom::getAllocatorIDFromPointer(connNode.data_ptr())),
      std::string(XS::name()) +
        axom::fmt::format(" execution space cannot use the connectivity allocator id {}",
                          axom::getAllocatorIDFromPointer(connNode.data_ptr())));
    SLIC_ERROR_IF(connNode.dtype().id() != conduitDataIdOfAxomIndexType,
                  "IntersectionShaper error: connectivity data type must be "
                  "axom::IndexType.");
    const auto* connPtr = static_cast<const axom::IndexType*>(connNode.data_ptr());
    axom::ArrayView<const axom::IndexType, 2> conn(connPtr, m_cellCount, NUM_VERTS_PER_HEX);

    const conduit::Node& coordNode = m_bpNodeInt["coordsets"][coordsetName];
    const conduit::Node& coordValues = coordNode.fetch_existing("values");
    axom::IndexType vertexCount = coordValues["x"].dtype().number_of_elements();
    bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
    int stride = isInterleaved ? NUM_COMPS_PER_VERT : 1;

    axom::StackArray<axom::ArrayView<const double>, 3> coordArrays {
      axom::ArrayView<const double>(coordValues["x"].as_double_ptr(), {vertexCount}, stride),
      axom::ArrayView<const double>(coordValues["y"].as_double_ptr(), {vertexCount}, stride),
      axom::ArrayView<const double>(coordValues["z"].as_double_ptr(), {vertexCount}, stride)};

    vertCoords = axom::Array<double>(m_cellCount * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT,
                                     m_cellCount * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT,
                                     allocId);
    auto vertCoordsView = vertCoords.view();

    axom::for_all<ExecSpace>(
      m_cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        // Get the indices of this element's vertices
        auto hexVerts = conn[i];

        // Get the coordinates for the vertices
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          auto vertId = hexVerts[j];
          for(int k = 0; k < NUM_COMPS_PER_VERT; k++)
          {
            vertCoordsView[(i * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT) + (j * NUM_COMPS_PER_VERT) + k] =
              coordArrays[k][vertId];
          }
        }
      });
  }
#endif  // AXOM_USE_CONDUIT

#if defined(AXOM_USE_MFEM)

  template <typename ExecSpace>
  void populateVertCoordsFromMFEMMesh(axom::Array<double>& vertCoords, int dim)
  {
    mfem::Mesh* mesh = getDC()->GetMesh();
    // Intersection algorithm only works on linear elements
    SLIC_ASSERT(mesh != nullptr);

    // Can only construct from 2D or 3D mesh
    SLIC_ASSERT(dim == 2 || dim == 3);

    if(m_cellCount > 0)
    {
      SLIC_ASSERT(mesh->GetNodes() == nullptr || mesh->GetNodes()->FESpace()->GetOrder(0));
    }

    // Allocation size is:
    // # of elements * # of vertices per cell * # of components per vertex
    int num_verts_per_cell;
    int num_comps_per_vert;

    // Quads
    if(dim == 2)
    {
      num_verts_per_cell = 4;
      num_comps_per_vert = 2;
    }
    // Hexes
    else
    {
      num_verts_per_cell = 8;
      num_comps_per_vert = 3;
    }

    // The MFEM mesh interface works only on host.
    // If on device, fill temporary host array then copy to device.
    axom::Array<double> tmpVertCoords;

    axom::Array<double>& fillVertCoords =
      axom::execution_space<ExecSpace>::onDevice() ? tmpVertCoords : vertCoords;
    fillVertCoords = axom::Array<double>(m_cellCount * num_verts_per_cell * num_comps_per_vert,
                                         m_cellCount * num_verts_per_cell * num_comps_per_vert);

    // Initialize vertices from mfem mesh and
    // set each shape volume fraction to 1

    auto fillVertCoordsView = fillVertCoords.view();
    for(int i = 0; i < m_cellCount; i++)
    {
      // Get the indices of this element's vertices
      mfem::Array<int> verts;
      mesh->GetElementVertices(i, verts);
      SLIC_ASSERT(verts.Size() == num_verts_per_cell);

      // Get the coordinates for the vertices
      for(int j = 0; j < num_verts_per_cell; ++j)
      {
        for(int k = 0; k < num_comps_per_vert; k++)
        {
          fillVertCoordsView[(i * num_verts_per_cell * num_comps_per_vert) +
                             (j * num_comps_per_vert) + k] = (mesh->GetVertex(verts[j]))[k];
        }
      }
    }
    if(vertCoords.data() != fillVertCoords.data())
    {
      axom::copy(vertCoords.data(), fillVertCoords.data(), sizeof(double) * vertCoords.size());
    }
  }

private:
  /// Create and return a new volume fraction grid function for the current mesh
  mfem::GridFunction* newVolFracGridFunction()
  {
    mfem::Mesh* mesh = getDC()->GetMesh();
    SLIC_ASSERT(mesh != nullptr);

    const int vfOrder = 0;
    const int dim = mesh->Dimension();
    mfem::L2_FECollection* coll = new mfem::L2_FECollection(vfOrder, dim, mfem::BasisType::Positive);
    mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, coll);
    mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
    volFrac->MakeOwner(coll);

    return volFrac;
  }
#endif  // AXOM_USE_MFEM

  // Check that surface mesh is composed of 3D Tetrahedra
  bool surfaceMeshIsTet() const
  {
    bool isTet = m_surfaceMesh != nullptr && m_surfaceMesh->getDimension() == 3 &&
      !m_surfaceMesh->hasMixedCellTypes() && m_surfaceMesh->getCellType() == mint::TET;
    return isTet;
  }

  // Check that surface mesh is composed of 2D Triangles
  bool surfaceMeshIsTri() const
  {
    bool isTri = m_surfaceMesh != nullptr && m_surfaceMesh->getDimension() == 2 &&
      !m_surfaceMesh->hasMixedCellTypes() && m_surfaceMesh->getCellType() == mint::TRIANGLE;
    return isTri;
  }

  /*!
   * \brief Returns the dimension of the computational mesh
   * \return dim The dimensions of the mesh (expected is 2 or 3, -1 in case of failure)
   */
  int getCompMeshDim()
  {
    int dim = -1;
#if defined(AXOM_USE_MFEM)
    if(m_dc != nullptr)
    {
      dim = this->getDC()->GetMesh()->SpaceDimension();
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bpGrp != nullptr)
    {
      std::string mesh_type = m_bpGrp->getView("topologies/mesh/elements/shape")->getString();
      if(mesh_type == "hex")
      {
        dim = 3;
      }
      else if(mesh_type == "quad")
      {
        dim = 2;
      }
    }
#endif

    SLIC_ERROR_IF(!(dim == 2 || dim == 3), "Invalid computational mesh dimension");
    return dim;
  }

private:
  int m_level {DEFAULT_CIRCLE_REFINEMENT_LEVEL};
  double m_revolvedVolume {DEFAULT_REVOLVED_VOLUME};
  std::string m_free_mat_name;

  //! \brief Volumes of cells in the computational mesh.
  axom::Array<double> m_cell_volumes;

  //! \brief Overlap volumes of cells in the computational mesh for the last shape.
  axom::Array<double> m_overlap_volumes;

  double m_vertexWeldThreshold {1.e-10};
  // Guard these to prevent warnings.
  int m_octcount {0};
  int m_tetcount {0};
  int m_tricount {0};

  axom::Array<OctahedronType> m_octs;
  axom::Array<TetrahedronType> m_tets;
  axom::Array<PolygonStaticType> m_tris;

  axom::Array<BoundingBox2D> m_aabbs_2d;
  axom::Array<BoundingBox3D> m_aabbs_3d;

  axom::Array<PolygonStaticType> m_quads;
  axom::Array<BoundingBox2D> m_quad_bbs;

  axom::Array<HexahedronType> m_hexes;
  axom::Array<BoundingBox3D> m_hex_bbs;
  axom::Array<TetrahedronType> m_tets_from_hexes_device;

  // Views of volume-fraction data owned by grid.
  std::vector<axom::ArrayView<double>> m_vf_grid_functions;
  std::vector<std::string> m_vf_material_names;
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_INTERSECTION_SHAPER__HPP_
