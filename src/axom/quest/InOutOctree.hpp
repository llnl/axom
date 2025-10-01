// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InOutOctree.hpp
 *
 * \brief Defines an InOutOctree for containment queries on a surface.
 */

#ifndef AXOM_QUEST_INOUT_OCTREE__HPP_
#define AXOM_QUEST_INOUT_OCTREE__HPP_

#include "axom/core.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"

#include "detail/inout/BlockData.hpp"
#include "detail/inout/MeshWrapper.hpp"
#include "detail/inout/InOutOctreeValidator.hpp"
#include "detail/inout/InOutOctreeStats.hpp"

#include "axom/fmt.hpp"

#include <vector>
#include <iterator>
#include <sstream>

#ifndef DUMP_VTK_MESH
//  #define DUMP_VTK_MESH
#endif

#ifndef DUMP_OCTREE_INFO
//  #define DUMP_OCTREE_INFO 1
#endif

#ifndef DEBUG_OCTREE_ACTIVE
//  #define DEBUG_OCTREE_ACTIVE
#endif

#if defined(DEBUG_OCTREE_ACTIVE) && defined(AXOM_DEBUG)
  #define QUEST_OCTREE_DEBUG_LOG_IF(_cond, _msg) \
    if(_cond) SLIC_DEBUG(_msg)
#else
  #define QUEST_OCTREE_DEBUG_LOG_IF(_cond, _msg) ((void)0)
#endif

// The following four variables are used in the QUEST_OCTREE_DEBUG_LOG_IF macro
// and are only active when DEBUG_OCTREE_ACTIVE is defined
#define DEBUG_VERT_IDX -2  // 1160
#define DEBUG_TRI_IDX -2   // 187820

#define DEBUG_BLOCK_1 BlockIndex::invalid_index()
//                     BlockIndex( {1346,1972,1691}, 12)
#define DEBUG_BLOCK_2 BlockIndex::invalid_index()
//                     BlockIndex( {336,493,423}, 10)

namespace axom
{
namespace quest
{
namespace detail
{
template <int DIM>
class InOutOctreeMeshDumper;

template <int DIM, typename Derived>
class InOutOctreeMeshDumperBase;

}  // namespace detail

/**
 * \class
 * \brief Handles generation of a point containment spatial index over a surface mesh
 *
 * The point containment queries determine whether a given arbitrary point in
 * space lies inside or outside of the surface.  This class depends on a
 * watertight surface mesh.  In order to repair common mesh defects, this
 * class modifies the Mesh passed to it.  Please discard all other copies
 * of the Mesh pointer.
 */
template <int DIM>
class InOutOctree : public spin::SpatialOctree<DIM, InOutBlockData>
{
private:
  friend class detail::InOutOctreeStats<DIM>;
  friend class detail::InOutOctreeValidator<DIM>;
  friend class detail::InOutOctreeMeshDumper<DIM>;
  friend class detail::InOutOctreeMeshDumperBase<DIM, detail::InOutOctreeMeshDumper<DIM>>;

public:
  using OctreeBaseType = spin::OctreeBase<DIM, InOutBlockData>;
  using SpatialOctreeType = spin::SpatialOctree<DIM, InOutBlockData>;

  using GeometricBoundingBox = typename SpatialOctreeType::GeometricBoundingBox;
  using SpacePt = typename SpatialOctreeType::SpacePt;
  using SpaceVector = typename SpatialOctreeType::SpaceVector;
  using BlockIndex = typename SpatialOctreeType::BlockIndex;
  using GridPt = typename OctreeBaseType::GridPt;
  using SpaceRay = primal::Ray<double, DIM>;

private:
  enum GenerationState
  {
    INOUTOCTREE_UNINITIALIZED,
    INOUTOCTREE_VERTICES_INSERTED,
    INOUTOCTREE_MESH_REORDERED,
    INOUTOCTREE_ELEMENTS_INSERTED,
    INOUTOCTREE_LEAVES_COLORED
  };

  static double DEFAULT_VERTEX_WELD_THRESHOLD;
  static double DEFAULT_BOUNDING_BOX_SCALE_FACTOR;

public:
  using SurfaceMesh = typename MeshWrapper<DIM>::SurfaceMesh;

  using VertexIndex = axom::IndexType;
  using CellIndex = axom::IndexType;
  using IndexRegistry = slam::FieldRegistry<slam::Set<VertexIndex>, VertexIndex>;

  using SpaceCell = typename MeshWrapper<DIM>::SpaceCell;

  using MeshVertexSet = typename MeshWrapper<DIM>::MeshVertexSet;
  using MeshElementSet = typename MeshWrapper<DIM>::MeshElementSet;
  using CellVertIndices = typename MeshWrapper<DIM>::CellVertIndices;
  using VertexIndexMap = typename MeshWrapper<DIM>::VertexIndexMap;

  // Type aliases for the Relations from Gray leaf blocks to mesh entities
  static const int MAX_VERTS_PER_BLOCK = 1;
  using VertexBlockMap = slam::Map<BlockIndex>;
  using STLIndirection = slam::policies::STLVectorIndirection<VertexIndex, VertexIndex>;

  using GrayLeafSet = slam::PositionSet<>;
  using BVStride = slam::policies::CompileTimeStride<VertexIndex, MAX_VERTS_PER_BLOCK>;
  using ConstantCardinality = slam::policies::ConstantCardinality<VertexIndex, BVStride>;
  using GrayLeafVertexRelation =
    slam::StaticRelation<VertexIndex, VertexIndex, ConstantCardinality, STLIndirection, GrayLeafSet, MeshVertexSet>;

  using VariableCardinality = slam::policies::VariableCardinality<VertexIndex, STLIndirection>;
  using GrayLeafElementRelation =
    slam::StaticRelation<VertexIndex, VertexIndex, VariableCardinality, STLIndirection, GrayLeafSet, MeshElementSet>;
  using CellIndexSet = typename GrayLeafElementRelation::RelationSubset;

  using GrayLeafsLevelMap = slam::Map<GrayLeafSet>;
  using GrayLeafVertexRelationLevelMap = slam::Map<GrayLeafVertexRelation>;
  using GrayLeafElementRelationLevelMap = slam::Map<GrayLeafElementRelation>;

public:
  /**
   * \brief Construct an InOutOctree to handle containment queries on a surface mesh
   *
   * \param [in] bb The spatial extent covered by the octree
   * \note We slightly scale the bounding box so all mesh elements are
   * guaranteed to be enclosed by the octree and to remedy problems we've
   * encountered related to meshes that are aligned with the octree grid
   * \note The InOutOctree modifies its mesh in an effort to repair common
   * problems. Please make sure to discard all old copies of the meshPtr.
   */
  InOutOctree(const GeometricBoundingBox& bb, std::shared_ptr<SurfaceMesh>& meshPtr)
    : SpatialOctreeType(GeometricBoundingBox(bb).scale(DEFAULT_BOUNDING_BOX_SCALE_FACTOR))
    , m_meshWrapper(meshPtr)
    , m_vertexToBlockMap(&m_meshWrapper.vertexSet())
    //
    , m_grayLeafsMap(&this->m_levels)
    , m_grayLeafToVertexRelationLevelMap(&this->m_levels)
    , m_grayLeafToElementRelationLevelMap(&this->m_levels)
    //
    , m_generationState(INOUTOCTREE_UNINITIALIZED)
  {
    setVertexWeldThreshold(DEFAULT_VERTEX_WELD_THRESHOLD);
  }

  /// \brief Generate the spatial index over the surface mesh
  void generateIndex();

  /**
   * \brief The point containment query.
   *
   * \param pt The point at which we are checking for containment
   * \return True if the point is within (or on) the surface, false otherwise
   * \note Points outside the octree bounding box are considered outside
   */
  bool within(const SpacePt& pt) const;

  /**
   * \brief Sets the threshold for welding vertices during octree construction
   *
   * \param [in] thresh The cutoff distance at which we consider two vertices
   * to be identical during the octree construction
   *
   * \pre thresh >= 0
   * \pre This function cannot be called after the octree has been constructed
   *
   * \note The InOutOctree requires the input surface to be watertight so this
   * parameter should be set with care. A welding threshold that is too
   * high could unnecessarily merge vertices and create topological defects,
   * while a value that is too low risks leaving gaps in meshes with tolerances
   * between vertices. The default value tends to work well in practice.
   *
   * \note The code actually uses the square of the threshold for comparisons
   */
  void setVertexWeldThreshold(double thresh)
  {
    SLIC_WARNING_IF(thresh < 0., "Distance threshold for vertices cannot be negative.");

    SLIC_WARNING_IF(
      m_generationState > INOUTOCTREE_UNINITIALIZED,
      "Can only set the vertex welding threshold " << "before initializing the InOutOctree");

    m_vertexWeldThresholdSquared = thresh * thresh;
  }

private:
  /**
   * \brief Helper function to insert a vertex into the octree
   *
   * \param idx The index of the vertex that we are inserting
   * \param startingLevel (optional, default = 0) The octree level at which
   * to begin the search for the leaf node covering this vertex
   */
  void insertVertex(VertexIndex idx, int startingLevel = 0);

  /// \brief Insert all mesh cells into the octree, generating a PM octree
  void insertMeshCells();

  /**
   * \brief Set a color for each leaf block of the octree.
   *
   * Black blocks are entirely within the surface, white blocks are entirely
   * outside the surface and Gray blocks intersect the surface.
   */
  void colorOctreeLeaves();

  /**
   * Use octree index over mesh vertices to convert the input 'soup of cells'
   * into an indexed mesh representation.
   * In particular, all vertices in the mesh that are nearly coincident will be
   * merged and degenerate cells (where the three vertices do not have unique
   * indices) will be removed.
   */
  void updateSurfaceMeshVertices();

private:
  /**
   * \brief Checks if all indexed cells in the block share a common vertex
   *
   * \param leafBlock [in] The current octree block
   * \param leafData [inout] The data associated with this block
   * \note A side effect of this function is that we set the leafData's vertex
   * to the common vertex if one is found
   * \return True, if all cells indexed by this leaf share a common vertex, false otherwise.
   */
  bool allCellsIncidentInCommonVertex(const BlockIndex& leafBlock,
                                      DynamicGrayBlockData& leafData) const;

  /**
   * \brief Finds a color for the given block \a blk and propagates to neighbors
   *
   * \note Propagates color to same-level and coarser level neighbors
   * \param blk The block to color
   * \param blkData The data associated with this block
   * \return True if we were able to find a color for \a blk, false otherwise
   */
  bool colorLeafAndNeighbors(const BlockIndex& blk, InOutBlockData& blkData);

  /**
   * \brief Predicate to determine if the vertex is indexed by the blk
   *
   * \pre This function assumes the vertices have been inserted and the mesh has been reordered
   * \param vIdx The index of the vertex to check
   * \param blk The block that we are checking against
   * \return true if \a vIdx is indexed by \a blk, false otherwise
   */
  bool blockIndexesVertex(VertexIndex vIdx, const BlockIndex& blk) const
  {
    SLIC_ASSERT(m_generationState >= INOUTOCTREE_MESH_REORDERED);

    // Needs to account for non-leaf ancestors of the block
    return vIdx >= 0 && m_vertexToBlockMap[vIdx].isDescendantOf(blk);
  }

  /**
   * \brief Predicate to determine if any of the elements vertices are indexed by the given BlockIndex
   *
   * \pre This function assumes the vertices have been inserted and the mesh has been reordered
   * \param tIdx The index of the cell to check
   * \param blk The block that we are checking against
   * \return true if one of the cell's vertices are indexed by \a blk, false otherwise
   */
  bool blockIndexesElementVertex(CellIndex tIdx, const BlockIndex& blk) const
  {
    SLIC_ASSERT(m_generationState >= INOUTOCTREE_MESH_REORDERED);

    CellVertIndices tVerts = m_meshWrapper.cellVertexIndices(tIdx);
    for(int i = 0; i < tVerts.size(); ++i)
    {
      // Using the vertex-to-block cache to avoid numerical degeneracies
      if(blockIndexesVertex(tVerts[i], blk))
      {
        return true;
      }
    }
    return false;
  }

  /**
   * \brief Determines whether the specified 3D point is within the gray leaf
   *
   * \param queryPt The point we are querying
   * \param leafBlk The block of the gray leaf
   * \param data The data associated with the leaf block
   * \return True, if the point is inside the local surface associated with this
   * block, false otherwise
   */
  template <int TDIM>
  typename std::enable_if<TDIM == 3, bool>::type withinGrayBlock(const SpacePt& queryPt,
                                                                 const BlockIndex& leafBlk,
                                                                 const InOutBlockData& data) const;

  /**
   * \brief Determines whether the specified 2D point is within the gray leaf
   *
   * \param queryPt The point we are querying
   * \param leafBlk The block of the gray leaf
   * \param data The data associated with the leaf block
   * \return True, if the point is inside the local surface associated with this
   * block, false otherwise
   */
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2, bool>::type withinGrayBlock(const SpacePt& queryPt,
                                                                 const BlockIndex& leafBlk,
                                                                 const InOutBlockData& data) const;

  /**
   * \brief Returns the index of the mesh vertex associated with the given leaf block
   *
   * \pre leafBlk is a leaf block of the octree
   * \param leafBlk The BlockIndex of a leaf block in the octree
   * \param leafData The data associated with this leaf block
   * \return The index of the mesh vertex associated with this leaf block
   */
  VertexIndex leafVertex(const BlockIndex& leafBlk, const InOutBlockData& leafData) const;

  /**
   * \brief Returns the set of mesh cell indices associated with the given leaf block
   *
   * \pre leafBlk is a leaf block of the octree
   * \param leafBlk The BlockIndex of a leaf block in the octree
   * \param leafData The data associated with this leaf block
   * \return The set of mesh cell indices associated with this leaf block
   */
  CellIndexSet leafCells(const BlockIndex& leafBlk, const InOutBlockData& leafData) const;

private:
  DISABLE_COPY_AND_ASSIGNMENT(InOutOctree);
  DISABLE_MOVE_AND_ASSIGNMENT(InOutOctree);

  /// \brief Checks internal consistency of the octree representation
  void checkValid() const;

  /// \brief Helper function to verify that all leaves at the given level have a color
  void checkAllLeavesColoredAtLevel(int AXOM_DEBUG_PARAM(level)) const;

  void dumpOctreeMeshVTK(const std::string& name) const;
  void dumpSurfaceMeshVTK(const std::string& name) const;

  /**
   * \brief Utility function to dump any Inside blocks whose neighbors are
   * outside (and vice-versa)
   *
   * \note There should not be any such blocks in a valid InOutOctree
   */
  void dumpDifferentColoredNeighborsMeshVTK(const std::string& name) const;

  /// \brief Utility function to print some statistics about the InOutOctree instance
  void printOctreeStats() const;

protected:
  MeshWrapper<DIM> m_meshWrapper;

  VertexBlockMap m_vertexToBlockMap;

  GrayLeafsLevelMap m_grayLeafsMap;
  GrayLeafVertexRelationLevelMap m_grayLeafToVertexRelationLevelMap;
  GrayLeafElementRelationLevelMap m_grayLeafToElementRelationLevelMap;

  GenerationState m_generationState;

  IndexRegistry m_indexRegistry;

  double m_vertexWeldThresholdSquared;

  /// Bounding box scaling factor for dealing with grazing triangles
  double m_boundingBoxScaleFactor {DEFAULT_BOUNDING_BOX_SCALE_FACTOR};
};

template <int DIM>
double InOutOctree<DIM>::DEFAULT_VERTEX_WELD_THRESHOLD = 1E-9;

template <int DIM>
double InOutOctree<DIM>::DEFAULT_BOUNDING_BOX_SCALE_FACTOR = 1.000123;

namespace
{
#ifdef AXOM_DEBUG
/// \brief Utility function to print the vertex indices of a cell
inline std::ostream& operator<<(std::ostream& os, const InOutOctree<3>::CellVertIndices& tvInd)
{
  os << axom::fmt::format("[{}]", fmt::join(tvInd, ","));
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const InOutOctree<3>::CellIndexSet& tSet)
{
  os << axom::fmt::format("[{}]", fmt::join(tSet, ","));
  return os;
}

#endif
}  // namespace

template <int DIM>
void InOutOctree<DIM>::generateIndex()
{
  using Timer = axom::utilities::Timer;

  // Loop through mesh vertices
  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "  Generating InOutOctree over surface mesh with "
                                   "{:L} vertices and {:L} elements.",
                                   m_meshWrapper.numMeshVertices(),
                                   m_meshWrapper.numMeshCells()));

  Timer timer;
  AXOM_ANNOTATE_SCOPE("InOutOctree::generateIndex");

  // STEP 1 -- Add mesh vertices to octree
  {
    AXOM_ANNOTATE_SCOPE("insert vertices");
    timer.reset();
    timer.start();
    int numMeshVerts = m_meshWrapper.numMeshVertices();
    for(int idx = 0; idx < numMeshVerts; ++idx)
    {
      insertVertex(idx);
    }
    timer.stop();
    m_generationState = INOUTOCTREE_VERTICES_INSERTED;
  }
  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "\t--Inserting vertices took {:.3Lf} seconds.",
                                   timer.elapsed()));

  // STEP 1(b) -- Update the mesh vertices and cells with after vertex welding from octree
  {
    AXOM_ANNOTATE_SCOPE("update surface mesh vertices");
    timer.reset();
    timer.start();
    updateSurfaceMeshVertices();
    timer.stop();
    m_generationState = INOUTOCTREE_MESH_REORDERED;
  }
  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "\t--Updating mesh took {:.3Lf} seconds.",
                                   timer.elapsed()));
  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "  After inserting vertices, reindexed mesh has {:L} "
                                   "vertices and {:L} cells.",
                                   m_meshWrapper.numMeshVertices(),
                                   m_meshWrapper.numMeshCells()));

#ifdef DUMP_OCTREE_INFO
  // -- Print some stats about the octree
  SLIC_INFO_ROOT("** Octree stats after inserting vertices");
  {
    AXOM_ANNOTATE_SCOPE("dump stats after inserting vertices");
    dumpSurfaceMeshVTK("surfaceMesh");
    dumpOctreeMeshVTK("prOctree");
    printOctreeStats();
  }
#endif
  {
    AXOM_ANNOTATE_SCOPE("validate after inserting vertices");
    checkValid();
  }

  // STEP 2 -- Add mesh cells (segments in 2D; triangles in 3D) to octree
  {
    AXOM_ANNOTATE_SCOPE("insert mesh cells");
    timer.reset();
    timer.start();
    insertMeshCells();
    timer.stop();
    m_generationState = INOUTOCTREE_ELEMENTS_INSERTED;
  }
  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "\t--Inserting cells took {:.3Lf} seconds.",
                                   timer.elapsed()));

  // STEP 3 -- Color the blocks of the octree
  // -- Black (in), White(out), Gray(Intersects surface)
  {
    AXOM_ANNOTATE_SCOPE("color octree leaves");
    timer.reset();
    timer.start();
    colorOctreeLeaves();

    timer.stop();
    m_generationState = INOUTOCTREE_LEAVES_COLORED;
  }
  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "\t--Coloring octree leaves took {:.3Lf} seconds.",
                                   timer.elapsed()));

#ifdef DUMP_OCTREE_INFO
  // -- Print some stats about the octree
  SLIC_INFO_ROOT("** Octree stats after inserting cells");
  {
    AXOM_ANNOTATE_SCOPE("dump stats after inserting cells");
    dumpOctreeMeshVTK("pmOctree");
    dumpDifferentColoredNeighborsMeshVTK("differentNeighbors");
    printOctreeStats();
  }
#endif
  {
    AXOM_ANNOTATE_SCOPE("validate after inserting cells");
    checkValid();
  }

  // CLEANUP -- Finally, fix up the surface mesh after octree operations
  {
    AXOM_ANNOTATE_SCOPE("regenerate surface mesh");
    timer.reset();
    timer.start();
    m_meshWrapper.regenerateSurfaceMesh();
    timer.stop();
  }
  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "\t--Regenerating the mesh took {:.3Lf} seconds.",
                                   timer.elapsed()));

  SLIC_INFO_ROOT("  Finished generating the InOutOctree.");
}

template <int DIM>
void InOutOctree<DIM>::insertVertex(VertexIndex idx, int startingLevel)
{
  const SpacePt pt = m_meshWrapper.getMeshVertexPosition(idx);

  BlockIndex block = this->findLeafBlock(pt, startingLevel);
  InOutBlockData& blkData = (*this)[block];

  QUEST_OCTREE_DEBUG_LOG_IF(
    idx == DEBUG_VERT_IDX,
    axom::fmt::format("\t -- inserting pt {} with index {}. "
                      "Looking at block {} w/ blockBB {} indexing leaf vertex {}",
                      pt,
                      idx,
                      axom::fmt::streamed(block),
                      this->blockBoundingBox(block),
                      blkData.dataIndex()));

  if(!blkData.hasData())
  {
    blkData.setData(idx);

    // Update the vertex-to-block map for this vertex
    if(m_generationState >= INOUTOCTREE_MESH_REORDERED)
    {
      m_vertexToBlockMap[idx] = block;
    }
  }
  else
  {
    // check if we should merge the vertices
    VertexIndex origVertInd = blkData.dataIndex();
    if(squared_distance(pt, m_meshWrapper.getMeshVertexPosition(origVertInd)) >=
       m_vertexWeldThresholdSquared)
    {
      blkData.clear();
      this->refineLeaf(block);

      insertVertex(origVertInd, block.childLevel());
      insertVertex(idx, block.childLevel());
    }
  }

  QUEST_OCTREE_DEBUG_LOG_IF(
    blkData.dataIndex() == DEBUG_VERT_IDX,
    axom::fmt::format("-- vertex {} is indexed in block {}. Leaf vertex is {}",
                      idx,
                      axom::fmt::streamed(block),
                      blkData.dataIndex()));
}

template <int DIM>
void InOutOctree<DIM>::insertMeshCells()
{
  using Timer = axom::utilities::Timer;
  using LeavesLevelMap = typename OctreeBaseType::OctreeLevelType;

  SLIC_ASSERT(m_meshWrapper.meshWasReindexed());

  // Temporary arrays of DyamicGrayBlockData for current and next level
  using DynamicLevelData = std::vector<DynamicGrayBlockData>;
  const int NUM_INIT_DATA_ENTRIES = 1 << 10;
  DynamicLevelData currentLevelData;
  DynamicLevelData nextLevelData;
  currentLevelData.reserve(NUM_INIT_DATA_ENTRIES);
  nextLevelData.reserve(NUM_INIT_DATA_ENTRIES);

  /// --- Initialize root level data
  BlockIndex rootBlock = this->root();
  InOutBlockData& rootData = (*this)[rootBlock];

  currentLevelData.push_back(DynamicGrayBlockData());
  DynamicGrayBlockData& dynamicRootData = currentLevelData[0];
  if(rootData.hasData())
  {
    dynamicRootData.setVertex(rootData.dataIndex());
  }
  dynamicRootData.setLeafFlag(rootData.isLeaf());
  rootData.setData(0);

  // Add all cell references to the root
  {
    int const numCells = m_meshWrapper.numMeshCells();
    dynamicRootData.cells().reserve(numCells);
    for(int idx = 0; idx < numCells; ++idx)
    {
      dynamicRootData.addCell(idx);
    }
  }

  // Iterate through octree levels
  // and insert cells into the blocks that they intersect
  for(int lev = 0; lev < this->m_levels.size(); ++lev)
  {
    Timer levelTimer(true);

    auto& gvRelData = m_indexRegistry.addNamelessBuffer();
    auto& geIndRelData = m_indexRegistry.addNamelessBuffer();
    auto& geSizeRelData = m_indexRegistry.addNamelessBuffer();
    geSizeRelData.push_back(0);

    int nextLevelDataBlockCounter = 0;

    auto& levelLeafMap = this->getOctreeLevel(lev);
    auto itEnd = levelLeafMap.end();
    for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
    {
      InOutBlockData& blkData = *it;

      if(!blkData.hasData())
      {
        continue;
      }

      BlockIndex blk(it.pt(), lev);
      DynamicGrayBlockData& dynamicLeafData = currentLevelData[blkData.dataIndex()];

      bool isInternal = !dynamicLeafData.isLeaf();
      bool isLeafThatMustRefine =
        !isInternal && !allCellsIncidentInCommonVertex(blk, dynamicLeafData);

      QUEST_OCTREE_DEBUG_LOG_IF(
        DEBUG_BLOCK_1 == blk || DEBUG_BLOCK_2 == blk,
        axom::fmt::format("Attempting to insert cells from block {}."
                          "\n\tDynamic data: {}"
                          "\n\tBlock data: {}"
                          "\n\tAbout to finalize? {}",
                          axom::fmt::streamed(blk),
                          dynamicLeafData,
                          blkData,
                          (!isInternal && !isLeafThatMustRefine ? " yes" : "no")));

      // Leaf blocks that don't refine are 'finalized'
      // -- add  them to the current level's relations
      if(!isInternal && !isLeafThatMustRefine)
      {
        if(dynamicLeafData.hasCells())
        {
          // Set the leaf data in the octree
          blkData.setData(static_cast<int>(gvRelData.size()));

          // Add the vertex index to the gray blocks vertex relation
          gvRelData.push_back(dynamicLeafData.vertexIndex());

          // Add the cells to the gray block's element relations
          std::copy(dynamicLeafData.cells().begin(),
                    dynamicLeafData.cells().end(),
                    std::back_inserter(geIndRelData));
          geSizeRelData.push_back(static_cast<int>(geIndRelData.size()));

          QUEST_OCTREE_DEBUG_LOG_IF(DEBUG_BLOCK_1 == blk || DEBUG_BLOCK_2 == blk,
                                    axom::fmt::format("[Added block {} into tree as a gray leaf]."
                                                      "\n\tDynamic data: {}"
                                                      "\n\tBlock data: {}",
                                                      axom::fmt::streamed(blk),
                                                      dynamicLeafData,
                                                      blkData));
        }
      }
      else
      {
        /// Otherwise, we must distribute the block data among the children

        // Refine the leaf if necessary
        if(isLeafThatMustRefine)
        {
          const VertexIndex vIdx = dynamicLeafData.vertexIndex();

          this->refineLeaf(blk);
          dynamicLeafData.setLeafFlag(false);

          // Reinsert the vertex into the tree, if vIdx was indexed by blk
          if(blockIndexesVertex(vIdx, blk))
          {
            insertVertex(vIdx, blk.childLevel());
          }
        }
        else if(isInternal)
        {
          // Need to mark the leaf as internal since we were using its data
          // as an index into the DynamicGrayBlockData array
          blkData.setInternal();
        }

        SLIC_ASSERT_MSG(
          this->isInternal(blk),
          axom::fmt::format("Block {} was refined, so it should be marked as internal.",
                            fmt::streamed(blk)));

        /// Setup caches for data associated with children
        BlockIndex childBlk[BlockIndex::NUM_CHILDREN];
        GeometricBoundingBox childBB[BlockIndex::NUM_CHILDREN];
        DynamicGrayBlockData childData[BlockIndex::NUM_CHILDREN];
        DynamicGrayBlockData* childDataPtr[BlockIndex::NUM_CHILDREN];

        const typename LeavesLevelMap::BroodData& broodData =
          this->getOctreeLevel(lev + 1).getBroodData(blk.pt());

        for(int j = 0; j < BlockIndex::NUM_CHILDREN; ++j)
        {
          childBlk[j] = blk.child(j);
          childBB[j] = this->blockBoundingBox(childBlk[j]);

          // expand bounding box slightly to deal with grazing cells
          childBB[j].scale(m_boundingBoxScaleFactor);

          const InOutBlockData& childBlockData = broodData[j];
          if(!childBlockData.hasData())
          {
            childData[j] = DynamicGrayBlockData();
            childData[j].setLeafFlag(childBlockData.isLeaf());
          }
          else
          {
            childData[j] = DynamicGrayBlockData(childBlockData.dataIndex(), childBlockData.isLeaf());
          }

          childDataPtr[j] = &childData[j];
        }

        // Check that the vector has enough capacity for all children
        // This ensures that our child data pointers will not be invalidated
        if(nextLevelData.capacity() < (nextLevelData.size() + BlockIndex::NUM_CHILDREN))
        {
          nextLevelData.reserve(nextLevelData.size() * 4);
        }

        // Add all cells to intersecting children blocks
        DynamicGrayBlockData::CellList& parentCells = dynamicLeafData.cells();
        int numCells = static_cast<int>(parentCells.size());
        for(int i = 0; i < numCells; ++i)
        {
          CellIndex tIdx = parentCells[i];
          SpaceCell spaceTri = m_meshWrapper.cellPositions(tIdx);
          GeometricBoundingBox tBB = m_meshWrapper.cellBoundingBox(tIdx);

          for(int j = 0; j < BlockIndex::numChildren(); ++j)
          {
            bool shouldAddCell = blockIndexesElementVertex(tIdx, childBlk[j]) ||
              (childDataPtr[j]->isLeaf() ? intersect(spaceTri, childBB[j])
                                         : intersect(tBB, childBB[j]));

            QUEST_OCTREE_DEBUG_LOG_IF(DEBUG_BLOCK_1 == childBlk[j] || DEBUG_BLOCK_2 == childBlk[j],
                                      //&& tIdx == DEBUG_TRI_IDX
                                      axom::fmt::format("Attempting to insert cell {} @ {} w/ BB {}"
                                                        "\n\t into block {} w/ BB {} and data {} "
                                                        "\n\tShould add? {}",
                                                        tIdx,
                                                        spaceTri,
                                                        tBB,
                                                        axom::fmt::streamed(childBlk[j]),
                                                        childBB[j],
                                                        *childDataPtr[j],
                                                        (shouldAddCell ? " yes" : "no")));

            if(shouldAddCell)
            {
              // Place the DynamicGrayBlockData in the array before adding its data
              if(!childDataPtr[j]->hasCells())
              {
                // Copy the DynamicGrayBlockData into the array
                nextLevelData.push_back(childData[j]);

                // Update the child data pointer
                childDataPtr[j] = &nextLevelData[nextLevelDataBlockCounter];

                // Set the data in the octree to this index and update the index
                (*this)[childBlk[j]].setData(nextLevelDataBlockCounter++);
              }

              childDataPtr[j]->addCell(tIdx);

              QUEST_OCTREE_DEBUG_LOG_IF(
                DEBUG_BLOCK_1 == childBlk[j] || DEBUG_BLOCK_2 == childBlk[j],
                //&& tIdx == DEBUG_TRI_IDX
                axom::fmt::format("Added cell {} @ {} with verts [{}]"
                                  "\n\tinto block {} with data {}.",
                                  tIdx,
                                  spaceTri,
                                  fmt::join(m_meshWrapper.cellVertexIndices(tIdx), ", "),
                                  axom::fmt::streamed(childBlk[j]),
                                  *(childDataPtr[j])));
            }
          }
        }
      }
    }

    if(!levelLeafMap.empty())
    {
      // Create the relations from gray leaves to mesh vertices and elements
      m_grayLeafsMap[lev] = GrayLeafSet(static_cast<int>(gvRelData.size()));

      m_grayLeafToVertexRelationLevelMap[lev] =
        GrayLeafVertexRelation(&m_grayLeafsMap[lev], &m_meshWrapper.vertexSet());
      m_grayLeafToVertexRelationLevelMap[lev].bindIndices(static_cast<int>(gvRelData.size()),
                                                          &gvRelData);

      m_grayLeafToElementRelationLevelMap[lev] =
        GrayLeafElementRelation(&m_grayLeafsMap[lev], &m_meshWrapper.elementSet());
      m_grayLeafToElementRelationLevelMap[lev].bindBeginOffsets(m_grayLeafsMap[lev].size(),
                                                                &geSizeRelData);
      m_grayLeafToElementRelationLevelMap[lev].bindIndices(static_cast<int>(geIndRelData.size()),
                                                           &geIndRelData);
    }

    currentLevelData.clear();
    nextLevelData.swap(currentLevelData);

    if(!levelLeafMap.empty())
    {
      SLIC_DEBUG("\tInserting cells into level " << lev << " took " << levelTimer.elapsed()
                                                 << " seconds.");
    }
  }
}

template <int DIM>
void InOutOctree<DIM>::colorOctreeLeaves()
{
  // Note (KW): Existence of leaf implies that either
  // * it is gray
  // * one of its siblings is gray
  // * one of its siblings has a gray descendant

  using Timer = axom::utilities::Timer;
  using GridPtVec = std::vector<GridPt>;
  GridPtVec uncoloredBlocks;

  // Bottom-up traversal of octree
  for(int lev = this->maxLeafLevel() - 1; lev >= 0; --lev)
  {
    uncoloredBlocks.clear();
    Timer levelTimer(true);

    auto& levelLeafMap = this->getOctreeLevel(lev);
    auto itEnd = levelLeafMap.end();
    for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
    {
      if(!it->isLeaf())
      {
        continue;
      }

      BlockIndex leafBlk(it.pt(), lev);
      InOutBlockData& blockData = *it;
      if(!colorLeafAndNeighbors(leafBlk, blockData))
      {
        uncoloredBlocks.push_back(leafBlk.pt());
      }
    }

    // Iterate through the uncolored blocks until all have a color
    // This terminates since we know that one of its siblings
    // (or their descendants) is gray
    while(!uncoloredBlocks.empty())
    {
      int prevCount = static_cast<int>(uncoloredBlocks.size());
      AXOM_UNUSED_VAR(prevCount);

      GridPtVec prevVec;
      prevVec.swap(uncoloredBlocks);
      auto end = prevVec.end();
      for(auto it = prevVec.begin(); it < end; ++it)
      {
        BlockIndex leafBlk(*it, lev);
        if(!colorLeafAndNeighbors(leafBlk, (*this)[leafBlk]))
        {
          uncoloredBlocks.push_back(*it);
        }
      }

      SLIC_ASSERT_MSG(static_cast<int>(uncoloredBlocks.size()) < prevCount,
                      axom::fmt::format("Problem coloring leaf blocks at level {}. "
                                        "There are {} blocks that are still not colored. "
                                        "First problem block is: {}",
                                        lev,
                                        uncoloredBlocks.size(),
                                        fmt::streamed(BlockIndex(uncoloredBlocks[0], lev))));
    }

    if(!levelLeafMap.empty())
    {
      checkAllLeavesColoredAtLevel(lev);
      SLIC_DEBUG(axom::fmt::format(axom::utilities::locale(),
                                   "\tColoring level {} took {:.3Lf} seconds.",
                                   lev,
                                   levelTimer.elapsed()));
    }
  }
}

template <int DIM>
bool InOutOctree<DIM>::colorLeafAndNeighbors(const BlockIndex& leafBlk, InOutBlockData& leafData)
{
  bool isColored = leafData.isColored();

  QUEST_OCTREE_DEBUG_LOG_IF(
    leafBlk == DEBUG_BLOCK_1 || leafBlk == DEBUG_BLOCK_2,
    axom::fmt::format("Trying to color {} with data: {}", axom::fmt::streamed(leafBlk), leafData));

  if(!isColored)
  {
    // Leaf does not yet have a color... try to find its color from same-level face neighbors
    for(int i = 0; !isColored && i < leafBlk.numFaceNeighbors(); ++i)
    {
      BlockIndex neighborBlk = leafBlk.faceNeighbor(i);
      if(this->isLeaf(neighborBlk))
      {
        const InOutBlockData& neighborData = (*this)[neighborBlk];

        QUEST_OCTREE_DEBUG_LOG_IF(
          DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk ||
            DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
          axom::fmt::format("Spreading color to block {} with data {}, "
                            "bounding box {} w/ midpoint {}"
                            "\n\t\t from block {} with data {}, "
                            "bounding box {} w/ midpoint {}.",
                            axom::fmt::streamed(leafBlk),
                            leafData,
                            this->blockBoundingBox(leafBlk),
                            this->blockBoundingBox(leafBlk).getCentroid(),
                            axom::fmt::streamed(neighborBlk),
                            neighborData,
                            this->blockBoundingBox(neighborBlk),
                            this->blockBoundingBox(neighborBlk).getCentroid()));

        switch(neighborData.color())
        {
        case InOutBlockData::Black:
          leafData.setBlack();
          break;
        case InOutBlockData::White:
          leafData.setWhite();
          break;
        case InOutBlockData::Gray:
        {
          SpacePt faceCenter = SpacePt::midpoint(this->blockBoundingBox(leafBlk).getCentroid(),
                                                 this->blockBoundingBox(neighborBlk).getCentroid());
          if(withinGrayBlock<DIM>(faceCenter, neighborBlk, neighborData))
          {
            leafData.setBlack();
          }
          else
          {
            leafData.setWhite();
          }
        }
        break;
        case InOutBlockData::Undetermined:
          break;
        }

        isColored = leafData.isColored();

        QUEST_OCTREE_DEBUG_LOG_IF(
          isColored && (DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk),
          axom::fmt::format("Leaf block was colored -- {} now has data {}",
                            axom::fmt::streamed(leafBlk),
                            leafData));
      }
    }
  }

  // If the block has a color, try to color its face neighbors at the same or coarser resolution
  if(isColored)
  {
    for(int i = 0; i < leafBlk.numFaceNeighbors(); ++i)
    {
      BlockIndex neighborBlk = this->coveringLeafBlock(leafBlk.faceNeighbor(i));
      if(neighborBlk != BlockIndex::invalid_index())
      {
        InOutBlockData& neighborData = (*this)[neighborBlk];
        if(!neighborData.isColored())
        {
          QUEST_OCTREE_DEBUG_LOG_IF(
            DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk ||
              DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
            axom::fmt::format("Spreading color from block {} with data {}, "
                              "bounding box {} w/ midpoint {}"
                              "\n\t\t to block {} with data {}, "
                              "bounding box {} w/ midpoint {}.",
                              axom::fmt::streamed(leafBlk),
                              leafData,
                              this->blockBoundingBox(leafBlk),
                              this->blockBoundingBox(leafBlk).getCentroid(),
                              axom::fmt::streamed(neighborBlk),
                              neighborData,
                              this->blockBoundingBox(neighborBlk),
                              this->blockBoundingBox(neighborBlk).getCentroid()));

          switch(leafData.color())
          {
          case InOutBlockData::Black:
            neighborData.setBlack();
            break;
          case InOutBlockData::White:
            neighborData.setWhite();
            break;
          case InOutBlockData::Gray:
          {
            // Check the center of the shared face between the same-level neighbor of the gray block
            SpacePt faceCenter =
              SpacePt::midpoint(this->blockBoundingBox(leafBlk).getCentroid(),
                                this->blockBoundingBox(leafBlk.faceNeighbor(i)).getCentroid());

            if(withinGrayBlock<DIM>(faceCenter, leafBlk, leafData))
            {
              neighborData.setBlack();
            }
            else
            {
              neighborData.setWhite();
            }
          }
          break;
          case InOutBlockData::Undetermined:
            break;
          }

          QUEST_OCTREE_DEBUG_LOG_IF(
            neighborData.isColored() && (DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk),
            axom::fmt::format("Neighbor block was colored -- {} now has data {}",
                              axom::fmt::streamed(neighborBlk),
                              neighborData));
        }
      }
    }
  }

  return isColored;
}

template <int DIM>
typename InOutOctree<DIM>::VertexIndex InOutOctree<DIM>::leafVertex(const BlockIndex& leafBlk,
                                                                    const InOutBlockData& leafData) const
{
  if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
  {
    SLIC_ASSERT(leafData.hasData());
    return m_grayLeafToVertexRelationLevelMap[leafBlk.level()][leafData.dataIndex()][0];
  }
  else
  {
    return leafData.dataIndex();
  }
}

template <int DIM>
typename InOutOctree<DIM>::CellIndexSet InOutOctree<DIM>::leafCells(const BlockIndex& leafBlk,
                                                                    const InOutBlockData& leafData) const
{
  SLIC_ASSERT(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED && leafData.hasData());

  return m_grayLeafToElementRelationLevelMap[leafBlk.level()][leafData.dataIndex()];
}

template <int DIM>
template <int TDIM>
typename std::enable_if<TDIM == 3, bool>::type InOutOctree<DIM>::withinGrayBlock(
  const SpacePt& queryPt,
  const BlockIndex& leafBlk,
  const InOutBlockData& leafData) const
{
  /// Finds a ray from queryPt to a point of a triangle within leafBlk.
  /// Then find the first triangle along this ray. The orientation of the ray
  /// against this triangle's normal indicates queryPt's containment.
  /// It is inside when the dot product is positive.

  SLIC_ASSERT(leafData.color() == InOutBlockData::Gray);
  SLIC_ASSERT(leafData.hasData());

  GeometricBoundingBox blockBB = this->blockBoundingBox(leafBlk);

  SpacePt triPt;

  CellIndexSet triSet = leafCells(leafBlk, leafData);
  const int numTris = triSet.size();
  for(int i = 0; i < numTris; ++i)
  {
    /// Get the triangle
    CellIndex idx = triSet[i];
    SpaceCell tri = m_meshWrapper.cellPositions(idx);

    /// Find a point from this triangle within the bounding box of the mesh
    primal::Polygon<double, DIM> poly = primal::clip(tri, blockBB);
    if(poly.numVertices() == 0)
    {
      // Account for cases where the triangle only grazes the bounding box.
      // Here, intersect(tri,blockBB) is true, but the clipping algorithm
      // produces an empty polygon.  To resolve this, clip against a
      // slightly expanded bounding box
      GeometricBoundingBox expandedBB = blockBB;
      expandedBB.scale(10 * m_boundingBoxScaleFactor);

      poly = primal::clip(tri, expandedBB);

      // If that still doesn't work, move on to the next triangle
      if(poly.numVertices() == 0)
      {
        continue;
      }
    }

    triPt = poly.vertexMean();

    /// Use a ray from the query point to the triangle point to find an
    /// intersection. Note: We have to check all triangles to ensure that
    /// there is not a closer triangle than tri along this direction.
    CellIndex tIdx = MeshWrapper<DIM>::NO_CELL;
    double minRayParam = axom::numeric_limits<double>::infinity();
    SpaceRay ray(queryPt, SpaceVector(queryPt, triPt));

    QUEST_OCTREE_DEBUG_LOG_IF(
      DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
      axom::fmt::format("Checking if pt {} is within block {} with data {}, "
                        "ray is {}, triangle point is {} on triangle with index {}.",
                        queryPt,
                        axom::fmt::streamed(leafBlk),
                        leafData,
                        ray,
                        triPt,
                        idx));

    double rayParam = 0;
    if(primal::intersect(tri, ray, rayParam))
    {
      minRayParam = rayParam;
      tIdx = idx;

      QUEST_OCTREE_DEBUG_LOG_IF(
        DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
        axom::fmt::format("... intersection for triangle w/ index {} at ray "
                          "parameter {} at point {}",
                          tIdx,
                          minRayParam,
                          ray.at(minRayParam)));
    }

    for(int j = 0; j < numTris; ++j)
    {
      CellIndex localIdx = triSet[j];
      if(localIdx == idx)
      {
        continue;
      }

      if(primal::intersect(m_meshWrapper.cellPositions(localIdx), ray, rayParam))
      {
        if(rayParam < minRayParam)
        {
          minRayParam = rayParam;
          tIdx = localIdx;

          QUEST_OCTREE_DEBUG_LOG_IF(
            DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
            axom::fmt::format("... intersection for triangle w/ index {} at ray "
                              "parameter {} at point {}",
                              tIdx,
                              minRayParam,
                              ray.at(minRayParam)));
        }
      }
    }

    if(tIdx == MeshWrapper<DIM>::NO_CELL)
    {
      continue;
    }

    // Inside when the dot product of the normal with this triangle is positive
    SpaceVector normal = (tIdx == idx) ? tri.normal() : m_meshWrapper.cellPositions(tIdx).normal();

    return normal.dot(ray.direction()) > 0.;
  }

  // SLIC_DEBUG("Could not determine inside/outside for point "
  //            << queryPt << " on block " << leafBlk);

  return false;  // query points on boundary might get here -- revisit this.
}

template <int DIM>
template <int TDIM>
typename std::enable_if<TDIM == 2, bool>::type InOutOctree<DIM>::withinGrayBlock(
  const SpacePt& queryPt,
  const BlockIndex& leafBlk,
  const InOutBlockData& leafData) const
{
  /// Finds a ray from queryPt to a point of a segment within leafBlk.
  /// Then finds the first segment along this ray. The orientation of the ray
  /// against this segment's normal indicates queryPt's containment.
  /// It is inside when the dot product is positive.

  SLIC_ASSERT(leafData.color() == InOutBlockData::Gray);
  SLIC_ASSERT(leafData.hasData());

  GeometricBoundingBox blockBB = this->blockBoundingBox(leafBlk);
  GeometricBoundingBox expandedBB = blockBB;
  expandedBB.scale(m_boundingBoxScaleFactor);

  SpacePt segmentPt;

  CellIndexSet segmentSet = leafCells(leafBlk, leafData);
  const int numSegments = segmentSet.size();
  for(int i = 0; i < numSegments; ++i)
  {
    /// Get the segment
    CellIndex idx = segmentSet[i];
    SpaceCell seg = m_meshWrapper.cellPositions(idx);

    /// Find a point from this segment within the expanded bounding box of the mesh
    // We'll use the midpoint of the segment after clipping it against the bounding box
    double pMin, pMax;
    const bool intersects = primal::intersect(seg, expandedBB, pMin, pMax);
    if(!intersects)
    {
      continue;
    }
    segmentPt = seg.at(0.5 * (pMin + pMax));

    // Using a ray from query pt to point on this segment
    // Find closest intersection to surface within cell inside this bounding box
    CellIndex tIdx = MeshWrapper<DIM>::NO_CELL;
    double minRayParam = axom::numeric_limits<double>::infinity();
    double minSegParam = axom::numeric_limits<double>::infinity();
    SpaceRay ray(queryPt, SpaceVector(queryPt, segmentPt));

    QUEST_OCTREE_DEBUG_LOG_IF(
      DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
      axom::fmt::format("Checking if pt {} is within block {} with data {}, "
                        "ray is {}, segment point is {} on segment with index {}.",
                        queryPt,
                        axom::fmt::streamed(leafBlk),
                        leafData,
                        ray,
                        segmentPt,
                        idx));

    double rayParam = 0;
    double segParam = 0;
    if(primal::intersect(ray, seg, rayParam, segParam))
    {
      minRayParam = rayParam;
      minSegParam = segParam;
      tIdx = idx;

      QUEST_OCTREE_DEBUG_LOG_IF(DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
                                axom::fmt::format("... intersection for segment w/ index {} "
                                                  "at ray parameter {} and segment parameter {} "
                                                  "is at point {}",
                                                  tIdx,
                                                  minRayParam,
                                                  minSegParam,
                                                  ray.at(minRayParam)));
    }

    for(int j = 0; j < numSegments; ++j)
    {
      CellIndex localIdx = segmentSet[j];
      if(localIdx == idx)
      {
        continue;
      }

      if(primal::intersect(ray, m_meshWrapper.cellPositions(localIdx), rayParam, segParam))
      {
        if(rayParam < minRayParam)
        {
          minRayParam = rayParam;
          minSegParam = segParam;
          tIdx = localIdx;

          QUEST_OCTREE_DEBUG_LOG_IF(
            DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
            axom::fmt::format("... intersection for segment w/ index {} "
                              "at ray parameter {} and segment parameter {} "
                              "is at point {}",
                              tIdx,
                              minRayParam,
                              minSegParam,
                              ray.at(minRayParam)));
        }
      }
    }
    if(tIdx == MeshWrapper<DIM>::NO_CELL)
    {
      continue;
    }

    // Get the surface normal at the intersection point
    // If the latter is a vertex, the normal is the average of its two incident segments
    SpaceVector normal = m_meshWrapper.surfaceNormal(tIdx, minSegParam, segmentSet);

    QUEST_OCTREE_DEBUG_LOG_IF(
      DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
      axom::fmt::format("... normal {}, ray {}, dot {}", normal, ray, normal.dot(ray.direction())));

    // Query point is inside when the dot product of the normal with ray is positive
    return normal.dot(ray.direction()) > 0.;
  }

  SLIC_DEBUG("Could not determine inside/outside for point " << queryPt << " on block " << leafBlk);

  return false;  // query points on boundary might get here -- revisit this.
}

template <int DIM>
void InOutOctree<DIM>::updateSurfaceMeshVertices()
{
  // Create a map from old vertex indices to new vertex indices
  MeshVertexSet origVerts(m_meshWrapper.numMeshVertices());
  VertexIndexMap vertexIndexMap(&origVerts, MeshWrapper<DIM>::NO_VERTEX);

  // Generate unique indices for new mesh vertices
  int uniqueVertexCounter = 0;
  for(int i = 0; i < origVerts.size(); ++i)
  {
    // Find the block and its indexed vertex in the octree
    BlockIndex leafBlock = this->findLeafBlock(m_meshWrapper.getMeshVertexPosition(i));
    SLIC_ASSERT((*this)[leafBlock].hasData());
    VertexIndex vInd = (*this)[leafBlock].dataIndex();

    // If the indexed vertex doesn't have a new id, give it one
    if(vertexIndexMap[vInd] == MeshWrapper<DIM>::NO_VERTEX)
    {
      vertexIndexMap[vInd] = uniqueVertexCounter++;
    }

    // If this is not the indexed vertex of the block, set the new index
    if(vInd != i)
    {
      vertexIndexMap[i] = vertexIndexMap[vInd];
    }
  }

  // Use the index map to reindex the mesh verts and elements
  m_meshWrapper.reindexMesh(uniqueVertexCounter, vertexIndexMap);

  // Update the octree leaf vertex IDs to the new mesh IDs
  // and create the map from the new vertices to their octree blocks
  m_vertexToBlockMap = VertexBlockMap(&m_meshWrapper.vertexSet());
  for(int i = 0; i < m_meshWrapper.numMeshVertices(); ++i)
  {
    const SpacePt& pos = m_meshWrapper.vertexPosition(i);
    BlockIndex leafBlock = this->findLeafBlock(pos);
    SLIC_ASSERT(this->isLeaf(leafBlock) && (*this)[leafBlock].hasData());

    (*this)[leafBlock].setData(i);
    m_vertexToBlockMap[i] = leafBlock;
  }
}

template <int DIM>
bool InOutOctree<DIM>::allCellsIncidentInCommonVertex(const BlockIndex& leafBlock,
                                                      DynamicGrayBlockData& leafData) const
{
  bool shareCommonVert = false;

  VertexIndex commonVert = leafData.vertexIndex();

  const int numCells = leafData.numCells();
  const auto& cells = leafData.cells();

  if(blockIndexesVertex(commonVert, leafBlock))
  {
    // This is a leaf node containing the indexed vertex
    // Loop through the triangles and check that all are incident with this
    // vertex
    for(int i = 0; i < numCells; ++i)
    {
      if(!m_meshWrapper.incidentInVertex(m_meshWrapper.cellVertexIndices(cells[i]), commonVert))
      {
        return false;
      }
    }
    shareCommonVert = true;
  }
  else
  {
    SLIC_ASSERT(numCells > 0);
    switch(numCells)
    {
    case 1:
      /// Choose an arbitrary vertex from this cell
      commonVert = m_meshWrapper.cellVertexIndices(cells[0])[0];
      shareCommonVert = true;
      break;
    case 2:
      /// Find a vertex that both triangles share
      shareCommonVert = m_meshWrapper.haveSharedVertex(cells[0], cells[1], commonVert);
      break;
    default:  // numCells > 3
      if(DIM == 3)
      {
        /// Find a vertex that the first three triangles share
        shareCommonVert = m_meshWrapper.haveSharedVertex(cells[0], cells[1], cells[2], commonVert);

        /// Check that all other triangles have this vertex
        for(int i = 3; shareCommonVert && i < numCells; ++i)
        {
          if(!m_meshWrapper.incidentInVertex(m_meshWrapper.cellVertexIndices(cells[i]), commonVert))
          {
            shareCommonVert = false;
          }
        }
      }
      break;
    }

    if(shareCommonVert)
    {
      leafData.setVertex(commonVert);
    }
  }

  return shareCommonVert;
}

template <int DIM>
bool InOutOctree<DIM>::within(const SpacePt& pt) const
{
  if(this->boundingBox().contains(pt))
  {
    const BlockIndex block = this->findLeafBlock(pt);
    const InOutBlockData& data = (*this)[block];

    switch(data.color())
    {
    case InOutBlockData::Black:
      return true;
    case InOutBlockData::White:
      return false;
    case InOutBlockData::Gray:
      return withinGrayBlock<DIM>(pt, block, data);
    case InOutBlockData::Undetermined:
      SLIC_ASSERT_MSG(false,
                      axom::fmt::format("Error -- All leaf blocks must have a color. The color of "
                                        "leafBlock {} was 'Undetermined' when querying point {}",
                                        fmt::streamed(block),
                                        pt));
      break;
    }
  }

  return false;
}

template <int DIM>
void InOutOctree<DIM>::printOctreeStats() const
{
  detail::InOutOctreeStats<DIM> octreeStats(*this);
  SLIC_INFO_ROOT(octreeStats.summaryStats());

#ifdef DUMP_VTK_MESH
  // Print out some debug meshes for vertex, triangle and/or blocks defined in
  // DEBUG_XXX macros
  if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
  {
    detail::InOutOctreeMeshDumper<DIM> meshDumper(*this);

    if(DEBUG_VERT_IDX >= 0 && DEBUG_VERT_IDX < m_meshWrapper.numMeshVertices())
    {
      meshDumper.dumpLocalOctreeMeshesForCell("debug_", DEBUG_VERT_IDX);
    }
    if(DEBUG_TRI_IDX >= 0 && DEBUG_TRI_IDX < m_meshWrapper.numMeshCells())
    {
      meshDumper.dumpLocalOctreeMeshesForCell("debug_", DEBUG_TRI_IDX);
    }

    if(DEBUG_BLOCK_1 != BlockIndex::invalid_index() && this->hasBlock(DEBUG_BLOCK_1))
    {
      meshDumper.dumpLocalOctreeMeshesForBlock("debug_", DEBUG_BLOCK_1);
    }

    if(DEBUG_BLOCK_2 != BlockIndex::invalid_index() && this->hasBlock(DEBUG_BLOCK_2))
    {
      meshDumper.dumpLocalOctreeMeshesForBlock("debug_", DEBUG_BLOCK_2);
    }
  }
#endif
}

template <int DIM>
void InOutOctree<DIM>::checkAllLeavesColoredAtLevel(int AXOM_DEBUG_PARAM(level)) const
{
#ifdef AXOM_DEBUG
  detail::InOutOctreeValidator<DIM> validator(*this);
  validator.checkAllLeavesColoredAtLevel(level);
#endif
}

template <int DIM>
void InOutOctree<DIM>::checkValid() const
{
#ifdef AXOM_DEBUG
  SLIC_DEBUG(axom::fmt::format("Inside InOutOctree::checkValid() to verify state of {} octree",
                               (m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED ? "PM" : "PR")));

  detail::InOutOctreeValidator<DIM> validator(*this);
  validator.checkValid();

  SLIC_DEBUG("done.");
#endif
}

template <int DIM>
void InOutOctree<DIM>::dumpSurfaceMeshVTK(const std::string& name) const
{
#ifdef DUMP_VTK_MESH

  detail::InOutOctreeMeshDumper<DIM> meshDumper(*this);
  meshDumper.dumpSurfaceMeshVTK(name);

#else
  AXOM_UNUSED_VAR(name);
#endif
}

template <int DIM>
void InOutOctree<DIM>::dumpOctreeMeshVTK(const std::string& name) const
{
#ifdef DUMP_VTK_MESH

  detail::InOutOctreeMeshDumper<DIM> meshDumper(*this);
  meshDumper.dumpOctreeMeshVTK(name);

#else
  AXOM_UNUSED_VAR(name);
#endif
}

template <int DIM>
void InOutOctree<DIM>::dumpDifferentColoredNeighborsMeshVTK(const std::string& name) const
{
#ifdef DUMP_VTK_MESH

  detail::InOutOctreeMeshDumper<DIM> meshDumper(*this);
  meshDumper.dumpDifferentColoredNeighborsMeshVTK(name);

#else
  AXOM_UNUSED_VAR(name);
#endif
}

}  // end namespace quest
}  // end namespace axom

// Note: The following needs to be included after InOutOctree is defined
#include "detail/inout/InOutOctreeMeshDumper.hpp"

#endif  // AXOM_QUEST_INOUT_OCTREE__HPP_
