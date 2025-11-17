// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_IMPLICIT_GRID__HPP_
#define AXOM_SPIN_IMPLICIT_GRID__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/spin/RectangularLattice.hpp"

#include <type_traits>
#include <vector>

namespace axom
{
namespace spin
{
/*!
 * \class ImplicitGrid
 *
 * \brief An implicit grid is an occupancy-based spatial index over an indexed
 * set of objects in space.
 *
 * An ImplicitGrid divides a given rectilinear slab of space (defined by an
 * axis aligned bounding box) into uniformly sized cells
 * of a specified resolution.
 * The GridCells of the ImplicitGrid index a subset of the items from an indexed
 * set whose cardinality is specified during the ImplicitGrid's initialization.
 * Users can insert items from the indexed set into an ImplicitGrid by providing
 * the item's bounding box and index.
 *
 * In contrast to a spin::UniformGrid, which encodes an array of indices
 * for each cell in the underlying multidimensional grid,
 * the ImplicitGrid encodes a single array of bins per dimension, each of which
 * has a bitset over the index space.  Thus, the storage overhead of an
 * ImplicitGrid is fixed at initialization time to
 *   \f$ numElts * sum_i { res[i] } \f$ bits.
 * Queries are implemented in terms of unions and intersections of bitsets.
 *
 * One might prefer an ImplicitGrid over a UniformGrid when one expects
 * a relatively dense index relative to the grid resolution (i.e. that
 * there will be many items indexed per bucket).  The ImplicitGrid
 * is designed for quick indexing and searching over a static (and relatively
 * small index space) in a relatively coarse grid.
 */
template <int NDIMS, typename ExecSpace = axom::SEQ_EXEC, typename TheIndexType = int>
class ImplicitGrid
{
public:
  using IndexType = TheIndexType;
  using GridCell = primal::Point<IndexType, NDIMS>;
  using SpacePoint = primal::Point<double, NDIMS>;
  using SpaceVec = primal::Vector<double, NDIMS>;

  using SpatialBoundingBox = primal::BoundingBox<double, NDIMS>;
  using LatticeType = RectangularLattice<NDIMS, double, IndexType>;

  using SizePolicy = slam::policies::RuntimeSize<IndexType>;
  using ElementSet = slam::OrderedSet<IndexType, IndexType, SizePolicy>;
  using BinSet = slam::OrderedSet<IndexType, IndexType, SizePolicy>;

  using BitsetType = slam::BitSet;
  using BinBitMap = slam::Map<BitsetType,
                              slam::Set<IndexType, IndexType>,
                              slam::policies::ArrayIndirection<IndexType, BitsetType>,
                              slam::policies::StrideOne<IndexType>>;

  struct QueryObject;

  /*!
   * \brief Default constructor for an ImplicitGrid
   *
   * \note Users must call initialize() to initialize the ImplicitGrid
   *       after constructing with the default constructor
   */
  ImplicitGrid() : m_initialized(false) { }

  /*!
   * \brief Constructor for an implicit grid from a bounding box,
   * a grid resolution a number of elements.
   *
   * \param [in] boundingBox Bounding box of domain to index
   * \param [in] gridRes Pointer to resolution for lattice covering bounding box
   * \param [in] numElts The number of elements to be indexed
   *
   * \pre \a gridRes is either NULL or has \a NDIMS coordinates
   * \sa initialize() for details on setting grid resolution
   * when \a gridRes is NULL
   */
  ImplicitGrid(const SpatialBoundingBox& boundingBox,
               const GridCell* gridRes,
               int numElts,
               int allocatorID = axom::execution_space<ExecSpace>::allocatorID())
    : m_bb(boundingBox)
    , m_initialized(false)
  {
    initialize(m_bb, gridRes, numElts, allocatorID);
  }

  /*!
   * \brief Constructor for an implicit grid from arrays of primitive types
   *
   * \param [in] bbMin Lower bounds of mesh bounding box
   * \param [in] bbMax Upper bounds of mesh bounding box
   * \param [in] gridRes Resolution for lattice covering mesh bounding box
   * \param [in] numElts The number of elements in the index space
   *
   * \pre \a bbMin and \a bbMax are not NULL and have \a NDIMS coordinates
   * \pre \a gridRes is either NULL or has \a NDIMS coordinates
   * \sa initialize() for details on setting grid resolution
   * when \a gridRes is NULL
   */
  ImplicitGrid(const double* bbMin,
               const double* bbMax,
               const int* gridRes,
               int numElts,
               int allocatorID = axom::execution_space<ExecSpace>::allocatorID())
    : m_initialized(false)
  {
    SLIC_ASSERT(bbMin != nullptr);
    SLIC_ASSERT(bbMax != nullptr);

    // Set up the grid resolution from the gridRes array
    //   if NULL, GridCell parameter to initialize should also be NULL

    GridCell res;
    if(gridRes != nullptr)
    {
      res = GridCell(gridRes, NDIMS);
    }

    initialize(SpatialBoundingBox(SpacePoint(bbMin), SpacePoint(bbMax)),
               (gridRes != nullptr) ? &res : nullptr,
               numElts,
               allocatorID);
  }

  /*! Predicate to check if the ImplicitGrid has been initialized */
  bool isInitialized() const { return m_initialized; }

  QueryObject getQueryObject() const;

  /*!
   * \brief Initializes an implicit grid or resolution gridRes over an axis
   * aligned domain covered by boundingBox. The implicit grid indexes a set
   * with numElts elements.
   *
   * \param [in] boundingBox Bounding box of domain to index
   * \param [in] gridRes Resolution for lattice covering bounding box
   * \param [in] numElts The number of elements to be indexed
   * \pre The ImplicitGrid has not already been initialized
   *
   * \note When \a gridRes is NULL, we use a heuristic to set the grid
   * resolution to the N^th root of \a numElts. We also ensure that
   * the resolution along any dimension is at least one.
   *
   */
  void initialize(const SpatialBoundingBox& boundingBox,
                  const GridCell* gridRes,
                  int numElts,
                  int allocatorID = axom::execution_space<ExecSpace>::allocatorID())
  {
    SLIC_ASSERT(!m_initialized);

    m_allocatorId = allocatorID;

    // Setup Grid Resolution, dealing with possible null pointer
    if(gridRes == nullptr)
    {
      // Heuristic: use the n^th root of the number of elements
      // add 0.5 to round to nearest integer
      double nthRoot = std::pow(static_cast<double>(numElts), 1. / NDIMS);
      int dimRes = static_cast<int>(nthRoot + 0.5);
      m_gridRes = GridCell(dimRes);
    }
    else
    {
      m_gridRes = GridCell(*gridRes);
    }

    // ensure that resolution in each dimension is at least one
    for(int i = 0; i < NDIMS; ++i)
    {
      m_gridRes[i] = axom::utilities::max(m_gridRes[i], IndexType(1));
    }

    // Setup lattice
    m_bb = boundingBox;
    m_lattice = spin::rectangular_lattice_from_bounding_box(boundingBox, m_gridRes.array());
    m_elementSet = ElementSet(numElts);

    for(int i = 0; i < NDIMS; ++i)
    {
      m_bins[i] = BinSet(m_gridRes[i]);
      m_binData[i] = BinBitMap(&m_bins[i], BitsetType(numElts, allocatorID), 1, allocatorID);

      axom::IndexType gridResDim = m_gridRes[i];
      m_minBlockBin[i] = axom::Array<IndexType>(gridResDim, gridResDim, allocatorID);
      m_maxBlockBin[i] = axom::Array<IndexType>(gridResDim, gridResDim, allocatorID);

      // We set initial min/max word indices to dummy values. These will be
      // set correctly on the first call to ImplicitGrid::insert().
      m_minBlockBin[i].fill(numElts);
      m_maxBlockBin[i].fill(0);
    }

    // Set the expansion factor for each element to a small fraction of the
    // grid's bounding boxes diameter
    // TODO: Add a constructor that allows users to set the expansion factor
    constexpr double EPS = 1e-8;
    m_expansionFactor = m_bb.range().norm() * EPS;
    SLIC_ASSERT(m_expansionFactor > 0.);

    m_initialized = true;
  }

  /*! Accessor for ImplicitGrid's resolution */
  const GridCell& gridResolution() const { return m_gridRes; }

  /*! Returns the number of elements in the ImplicitGrid's index set */
  int numIndexElements() const { return m_elementSet.size(); }

  /*!
   * \brief Inserts an element with index \a idx and bounding box \a bbox
   * into the implicit grid
   *
   * \param [in] bbox The bounding box of the element
   * \param [in] idx  The index of the element
   */
  void insert(const SpatialBoundingBox& bbox, IndexType idx)
  {
    if(axom::execution_space<ExecSpace>::onDevice())
    {
      int deviceAllocId = axom::execution_space<ExecSpace>::allocatorID();
      // Copy host box to device array
      ArrayView<const SpatialBoundingBox> bbox_host(&bbox, 1);
      Array<SpatialBoundingBox> bbox_device(bbox_host, deviceAllocId);
      insert(1, bbox_device.data(), idx);
    }
    else
    {
      insert(1, &bbox, idx);
    }
  }

  /*!
   * \brief Inserts a set of elements represented by bounding boxes \a bboxes
   *  and beginning at index \a startIdx into the implicit grid
   *
   * \param [in] nelems the number of elements to insert
   * \param [in] bboxes an array of bounding boxes for each element
   * \param [in] startIdx the index of the first bounding box in bboxes
   */
  void insert(IndexType nelems, const SpatialBoundingBox* bboxes, IndexType startIdx = 0)
  {
    SLIC_ASSERT(m_initialized);
    const double expansionFactor = m_expansionFactor;
    LatticeType lattice = m_lattice;

    BitsetType* binData[NDIMS];
    IndexType* minBlkBins[NDIMS];
    IndexType* maxBlkBins[NDIMS];
    IndexType highestBins[NDIMS];
    for(int i = 0; i < NDIMS; i++)
    {
      binData[i] = m_binData[i].data().data();
      highestBins[i] = m_binData[i].set()->size() - 1;
      minBlkBins[i] = m_minBlockBin[i].data();
      maxBlkBins[i] = m_maxBlockBin[i].data();
    }

    for_all<ExecSpace>(
      nelems,
      AXOM_LAMBDA(axom::IndexType ibox) {
        IndexType elemIdx = startIdx + ibox;

        SpatialBoundingBox scaledBox = bboxes[ibox];
        // Note: We slightly inflate the bbox of the objects.
        //       This effectively ensures that objects on grid boundaries are added
        //       in all nearby grid cells.
        scaledBox.expand(expansionFactor);

        const GridCell lowerCell = lattice.gridCell(scaledBox.getMin());
        const GridCell upperCell = lattice.gridCell(scaledBox.getMax());

        for(int idim = 0; idim < NDIMS; idim++)
        {
          const IndexType lower = axom::utilities::clampLower(lowerCell[idim], IndexType());
          const IndexType upper = axom::utilities::clampUpper(upperCell[idim], highestBins[idim]);

          const IndexType word = elemIdx / BitsetType::BitsPerWord;

          for(int j = lower; j <= upper; ++j)
          {
            binData[idim][j].atomicSet(elemIdx);
            axom::atomicMin<ExecSpace>(&minBlkBins[idim][j], word);
            axom::atomicMax<ExecSpace>(&maxBlkBins[idim][j], word);
          }
        }
      });
  }

  /*!
   * Finds the candidate elements in the vicinity of query point \a pt
   *
   * \param [in] pt The query point
   * \return A bitset \a bSet whose bits correspond to
   * the elements of the IndexSet.
   * The bits of \a bSet are set if their corresponding element bounding boxes
   * overlap the grid cell containing \a pt.
   */
  BitsetType getCandidates(const SpacePoint& pt) const
  {
    if(!m_initialized || !m_bb.contains(pt))
    {
      return BitsetType(0);
    }

    const GridCell gridCell = m_lattice.gridCell(pt);

    // Note: Need to clamp the upper range of the gridCell
    //       to handle points on the upper boundaries of the bbox
    //       This is valid since we've already ensured that pt is in the bbox.
    IndexType idx = axom::utilities::clampUpper(gridCell[0], highestBin(0));
    BitsetType res = m_binData[0][idx];

    for(int i = 1; i < NDIMS; ++i)
    {
      idx = axom::utilities::clampUpper(gridCell[i], highestBin(i));
      res &= m_binData[i][idx];
    }

    return res;
  }

  /*!
   * Finds the candidate elements in the given \a gridCell of the grid
   *
   * \param [in] gridCell The cell of the grid
   * \return A bitset whose bits correspond to the elements of the IndexSet.
   * The bits are set if their corresponding element bounding boxes overlap \a gridCell
   */
  BitsetType getCandidates(const GridCell& gridCell) const
  {
    // Perform some validity checks
    if(!m_initialized)
    {
      return BitsetType(0);
    }
    for(int i = 0; i < NDIMS; ++i)
    {
      if(gridCell[i] < 0 || gridCell[i] > highestBin(i))
      {
        return BitsetType(0);
      }
    }

    // Note: Due to above checks, gridCell[i] is always valid
    BitsetType res = m_binData[0][gridCell[0]];
    for(int i = 1; i < NDIMS; ++i)
    {
      res &= m_binData[i][gridCell[i]];
    }

    return res;
  }

  /*!
   * Finds the candidate elements in the vicinity of query box \a box
   *
   * \param [in] box The query box
   * \return A bitset \a bSet whose bits correspond to
   * the elements of the IndexSet.
   * The bits of \a bSet are set if their corresponding element bounding boxes
   * overlap the grid cell containing \a box
   */
  BitsetType getCandidates(const SpatialBoundingBox& box) const
  {
    if(!m_initialized || !m_bb.intersectsWith(box))
    {
      return BitsetType(0);
    }

    const GridCell lowerCell = m_lattice.gridCell(box.getMin());
    const GridCell upperCell = m_lattice.gridCell(box.getMax());

    BitsetType bits = getBitsInRange(0, lowerCell[0], upperCell[0]);

    for(int dim = 1; dim < NDIMS; ++dim)
    {
      bits &= getBitsInRange(dim, lowerCell[dim], upperCell[dim]);
    }

    return bits;
  }

  /*!
   * Returns an explicit list of candidates in the vicinity of a query object
   *
   * \tparam QueryGeom The type of the query object (e.g. point or box)
   * \param [in] query The query object
   * \return A list of indexes from the IndexSet whose corresponding
   * bounding boxes overlap the grid cell containing \a query
   *
   * \pre This function is implemented in terms of
   * ImplicitGrid::getCandidates(const QueryGeom& ). An overload for the actual
   * \a QueryGeom type (e.g. \a SpacePoint or \a SpatialBoundingBox) must exist.
   *
   * \note This function returns the same information as \a getCandidates(),
   * but in a different format. While the latter returns a bitset of the
   * candidates, this function returns an explicit list of indices.
   *
   * \sa getCandidates()
   */
  template <typename QueryGeom>
  std::vector<IndexType> getCandidatesAsArray(const QueryGeom& query) const
  {
    std::vector<IndexType> candidatesVec;

    BitsetType candidateBits = getCandidates(query);
    candidatesVec.reserve(candidateBits.count());
    for(IndexType eltIdx = candidateBits.find_first(); eltIdx != BitsetType::npos;
        eltIdx = candidateBits.find_next(eltIdx))
    {
      candidatesVec.push_back(eltIdx);
    }

    return candidatesVec;
  }

  /*!
   * \brief Returns a list of candidates in the vicinity of a set of query
   *  objects.
   *
   * \tparam QueryGeom The type of the query object (e.g. point or box)
   * \param [in] qsize The number of objects to query against the implicit grid
   * \param [in] queryObjs The array of query objects, of length qsize
   * \param [out] outOffsets Offsets into the candidates array for each query
   *  object
   * \param [out] outCounts The number of candidates for each query object
   * \param [out] outCandidates The candidate IDs for each query object
   *
   * \note outOffsets and outCounts should point to arrays allocated in a
   *  memory space accessible from the given execution space, and be of
   *  length qsize.
   *
   * \note Upon completion, the ith query point has:
   *  * counts[ i ] candidates
   *  * Stored in the candidates array in the following range:
   *    [ offsets[ i ], offsets[ i ]+counts[ i ] ]
   */
  template <typename QueryGeom>
  void getCandidatesAsArray(axom::IndexType qsize,
                            const QueryGeom* queryObjs,
                            axom::ArrayView<IndexType> outOffsets,
                            axom::ArrayView<IndexType> outCounts,
                            axom::Array<IndexType>& outCandidates) const;

  /// \overload
  void getCandidatesAsArray(axom::ArrayView<const SpacePoint> queryObjs,
                            axom::ArrayView<IndexType> outOffsets,
                            axom::ArrayView<IndexType> outCounts,
                            axom::Array<IndexType>& outCandidates) const
  {
    getCandidatesAsArray(queryObjs.size(), queryObjs.data(), outOffsets, outCounts, outCandidates);
  }

  /// \overload
  void getCandidatesAsArray(axom::ArrayView<const SpatialBoundingBox> queryObjs,
                            axom::ArrayView<IndexType> outOffsets,
                            axom::ArrayView<IndexType> outCounts,
                            axom::Array<IndexType>& outCandidates) const
  {
    getCandidatesAsArray(queryObjs.size(), queryObjs.data(), outOffsets, outCounts, outCandidates);
  }

  /*!
   * Tests whether grid cell gridPt indexes the element with index idx
   *
   * \param [in] gridCell The cell within the ImplicitGrid that we are testing
   * \param [in] idx An element index from the IndexSet to test
   *
   * \pre \a idx should be in the IndexSet.  I.e. 0 <= idx < numIndexElements()
   * \return True if the bounding box of element \a idx overlaps
   * with GridCell \a gridCell.
   */
  bool contains(const GridCell& gridCell, IndexType idx) const
  {
    bool ret = true;

    if(!m_elementSet.isValidIndex(idx))
    {
      ret = false;
    }

    for(int i = 0; i < NDIMS; ++i)
    {
      ret = ret && m_bins[i].isValidIndex(gridCell[i]) && m_binData[i][gridCell[i]].test(idx);
    }

    return ret;
  }

private:
  /*!
   * \brief Returns the bin index in the given dimension dim
   *
   * \pre 0 <= dim < NDIMS
   */
  IndexType highestBin(int dim) const
  {
    SLIC_ASSERT(0 <= dim && dim < NDIMS);
    return m_bins[dim].size() - 1;
  }

  /*!
   * \brief Queries the bits that are set for dimension \a dim
   * within the range of boxes \a lower to \a upper
   *
   * \param dim The dimension to check
   * \param lower The index of the lower bin in the range (inclusive)
   * \param upper The index of the upper bin in the range (inclusive)
   *
   * \return A bitset whose bits are set if they are set in
   * any of the boxes between \a lower and \a upper for
   * dimension \a dim
   *
   * \note We perform range checking to ensure that \a lower
   * is at least 0 and \a upper is at most \a highestBin(dim)
   *
   * \sa highestBin()
   */
  BitsetType getBitsInRange(int dim, int lower, int upper) const
  {
    // Note: Need to clamp the gridCell ranges since the input box boundaries
    //       are not restricted to the implicit grid's bounding box
    lower = axom::utilities::clampLower(lower, IndexType());
    upper = axom::utilities::clampUpper(upper, highestBin(dim));

    BitsetType bits = m_binData[dim][lower];
    for(int i = lower + 1; i <= upper; ++i)
    {
      bits |= m_binData[dim][i];
    }

    return bits;
  }

private:
  //! The bounding box of the ImplicitGrid
  SpatialBoundingBox m_bb;

  //! A lattice to help in converting from points in space to GridCells
  LatticeType m_lattice;

  //! The amount by which to expand bounding boxes
  double m_expansionFactor;

  //! Resolution of the ImplicitGrid
  GridCell m_gridRes;

  //! The index set of the elements
  ElementSet m_elementSet;

  //! A set of bins, per dimension
  BinSet m_bins[NDIMS];

  //! The data associated with each bin
  BinBitMap m_binData[NDIMS];

  //! The lowest word index in each bin with at least one bit set
  axom::Array<IndexType> m_minBlockBin[NDIMS];

  //! The highest word index in each bin with at least one bit set
  axom::Array<IndexType> m_maxBlockBin[NDIMS];

  //! The allocator ID to use
  int m_allocatorId;

  //! Tracks whether the ImplicitGrid has been initialized
  bool m_initialized;
};

}  // end namespace spin
}  // end namespace axom

//------------------------------------------------------------------------------
//  ImplicitGrid Implementation
//------------------------------------------------------------------------------

namespace axom
{
namespace spin
{
/// Device-copyable query object for running implicit grid queries on the GPU
template <int NDIMS, typename ExecSpace, typename IndexType>
struct ImplicitGrid<NDIMS, ExecSpace, IndexType>::QueryObject
{
public:
  using SpacePoint = primal::Point<double, NDIMS>;
  using SpatialBoundingBox = primal::BoundingBox<double, NDIMS>;

  using LatticeType = RectangularLattice<NDIMS, double, IndexType>;

  using BitsetType = slam::BitSet;
  using BinBitMap = slam::Map<BitsetType,
                              slam::Set<IndexType, IndexType>,
                              slam::policies::ArrayIndirection<IndexType, BitsetType>,
                              slam::policies::StrideOne<IndexType>>;

  QueryObject(const SpatialBoundingBox& spaceBb,
              const LatticeType& lattice,
              const BinBitMap (&binData)[NDIMS],
              const axom::Array<IndexType> (&minBlkBins)[NDIMS],
              const axom::Array<IndexType> (&maxBlkBins)[NDIMS])
    : m_bb(spaceBb)
    , m_lattice(lattice)
  {
    for(int idim = 0; idim < NDIMS; idim++)
    {
      m_highestBins[idim] = binData[idim].set()->size() - 1;
      m_binData[idim] = binData[idim].data().view();
      m_minBlkBin[idim] = minBlkBins[idim].view();
      m_maxBlkBin[idim] = maxBlkBins[idim].view();
    }
  }

  /*!
   * \brief Counts the number of elements in the implicit grid which may
   *  intersect with the given point.
   *
   * \param [in] pt the point to query the implicit grid against.
   *
   * \return numCandidates the number of candidates for the given point.
   */
  AXOM_HOST_DEVICE IndexType countCandidates(const SpacePoint& pt) const;

  /*!
   * \brief Counts the number of elements in the implicit grid which may
   *  intersect with the given bounding box.
   *
   * \param [in] bbox the bounding box to query the implicit grid against.
   *
   * \return numCandidates the number of candidates for the given bounding box.
   */
  AXOM_HOST_DEVICE IndexType countCandidates(const SpatialBoundingBox& bbox) const;

  /*!
   * \brief Iterates through the implicit grid, calling a given function for
   *  each candidate element which potentially intersects the given point.
   *
   * \param [in] pt the point to query the implicit grid against.
   * \param [in] candidateFunc the function object to be called for each
   *  intersection candidate
   *
   * \note The supplied functor `candidateFunc` is expected to take one argument,
   *  the index of the candidate element.
   *  The functor may optionally return a boolean, where a value of `true`
   *  terminates the candidate search early.
   */
  template <typename FuncType>
  AXOM_HOST_DEVICE void visitCandidates(const SpacePoint& pt, FuncType&& candidateFunc) const;

  /*!
   * \brief Iterates through the implicit grid, calling a given function for
   *  each candidate element which potentially intersects the given bounding box.
   *
   * \param [in] bbox the bounding box to query the implicit grid against.
   * \param [in] candidateFunc the function object to be called for each
   *  intersection candidate
   *
   * \note The supplied functor `candidateFunc` is expected to take one argument,
   *  the index of the candidate element.
   *  The functor may optionally return a boolean, where a value of `true`
   *  terminates the candidate search early.
   */
  template <typename FuncType>
  AXOM_HOST_DEVICE void visitCandidates(const SpatialBoundingBox& bbox,
                                        FuncType&& candidateFunc) const;

private:
  template <typename FuncType, typename ReturnType>
  struct VisitDispatch;

  template <typename FuncType>
  struct VisitDispatch<FuncType, void>
  {
    AXOM_HOST_DEVICE static bool getResult(FuncType&& type, int arg)
    {
      type(arg);
      return false;
    }
  };

  template <typename FuncType>
  struct VisitDispatch<FuncType, bool>
  {
    AXOM_HOST_DEVICE static bool getResult(FuncType&& type, int arg) { return type(arg); }
  };

  template <typename FuncType>
  AXOM_HOST_DEVICE bool getVisitResult(FuncType&& type, int arg) const
  {
    using ReturnType = std::invoke_result_t<FuncType, int>;
    return VisitDispatch<FuncType, ReturnType>::getResult(type, arg);
  }

  /*!
   * \brief Gets the expected range of word indices where bits may be set for
   *  a given bin coordinate.
   *
   * \param [in] cellIdx the bin indices
   * \param [out] minWord the lowest-indexed word where a bit may be set
   * \param [out] maxWord the highest-indexed word where a bit may be set
   */
  AXOM_HOST_DEVICE void getWordBounds(const GridCell& cellIdx,
                                      IndexType& minWord,
                                      IndexType& maxWord) const
  {
    minWord = m_minBlkBin[0][cellIdx[0]];
    maxWord = m_maxBlkBin[0][cellIdx[0]];
    for(int idim = 1; idim < NDIMS; idim++)
    {
      // Intersect with word ranges of other dimensions
      minWord = axom::utilities::max(m_minBlkBin[idim][cellIdx[idim]], minWord);
      maxWord = axom::utilities::min(m_maxBlkBin[idim][cellIdx[idim]], maxWord);
    }
  }

  /*!
   * \brief Gets the expected range of word indices where bits may be set for
   *  a given range of bin coordinate.
   *
   * \param [in] lowerRange the lower bound of bin coordinates
   * \param [in] upperRange the upper bound of bin coordinates
   * \param [out] minWord the lowest-indexed word where a bit may be set
   * \param [out] maxWord the highest-indexed word where a bit may be set
   */
  AXOM_HOST_DEVICE void getWordBounds(const GridCell& lowerRange,
                                      const GridCell& upperRange,
                                      IndexType& minWord,
                                      IndexType& maxWord) const
  {
    for(int idim = 0; idim < NDIMS; idim++)
    {
      IndexType minWordDim = m_minBlkBin[idim][lowerRange[idim]],
                maxWordDim = m_maxBlkBin[idim][upperRange[idim]];
      for(int ibin = lowerRange[idim] + 1; ibin <= upperRange[idim]; ibin++)
      {
        // Take the union of candidate bins within a dimension
        minWordDim = axom::utilities::min(m_minBlkBin[idim][ibin], minWordDim);
        maxWordDim = axom::utilities::max(m_maxBlkBin[idim][ibin], maxWordDim);
      }

      if(idim == 0)
      {
        minWord = minWordDim;
        maxWord = maxWordDim;
      }
      else
      {
        // Intersect the resulting ranges between dimensions
        minWord = axom::utilities::max(minWordDim, minWord);
        maxWord = axom::utilities::min(maxWordDim, maxWord);
      }
    }
  }

  //! The bounding box of the ImplicitGrid
  SpatialBoundingBox m_bb;

  //! A lattice to help in converting from points in space to GridCells
  LatticeType m_lattice;

  //! The highest bin index in each dimension
  IndexType m_highestBins[NDIMS];

  //! The data associated with each bin
  axom::ArrayView<const BitsetType> m_binData[NDIMS];

  //! The lowest word index in each bin with at least one bit set
  axom::ArrayView<const IndexType> m_minBlkBin[NDIMS];

  //! The highest word index in each bin with at least one bit set
  axom::ArrayView<const IndexType> m_maxBlkBin[NDIMS];
};

template <int NDIMS, typename ExecSpace, typename IndexType>
typename ImplicitGrid<NDIMS, ExecSpace, IndexType>::QueryObject
ImplicitGrid<NDIMS, ExecSpace, IndexType>::getQueryObject() const
{
  static_assert(std::is_copy_constructible<ImplicitGrid::QueryObject>::value,
                "ImplicitGrid::QueryObject must be copy-constructible.");

  SLIC_ASSERT(m_initialized);
  return QueryObject {m_bb, m_lattice, m_binData, m_minBlockBin, m_maxBlockBin};
}

template <int NDIMS, typename ExecSpace, typename IndexType>
template <typename QueryGeom>
void ImplicitGrid<NDIMS, ExecSpace, IndexType>::getCandidatesAsArray(
  axom::IndexType qsize,
  const QueryGeom* queryObjs,
  axom::ArrayView<IndexType> outOffsets,
  axom::ArrayView<IndexType> outCounts,
  axom::Array<IndexType>& outCandidates) const
{
  SLIC_ERROR_IF(outOffsets.size() < qsize, "outOffsets must have at least qsize elements");
  SLIC_ERROR_IF(outCounts.size() < qsize, "outCounts must have at least qsize elements");
  auto gridQuery = getQueryObject();

#ifdef AXOM_USE_RAJA
  axom::ReduceSum<ExecSpace, IndexType> totalCountReduce(0);
  // Step 1: count number of candidate intersections for each point
  for_all<ExecSpace>(
    qsize,
    AXOM_LAMBDA(IndexType i) {
      outCounts[i] = gridQuery.countCandidates(queryObjs[i]);
      totalCountReduce += outCounts[i];
    });

  // Step 2: exclusive scan for offsets in candidate array
  axom::exclusive_scan<ExecSpace>(outCounts, outOffsets);

  axom::IndexType totalCount = totalCountReduce.get();

  // Step 3: allocate memory for all candidates
  outCandidates = axom::Array<IndexType>(totalCount, totalCount, m_allocatorId);
  const auto candidates_v = outCandidates.view();

  // Step 4: fill candidate array for each query box
  for_all<ExecSpace>(
    qsize,
    AXOM_LAMBDA(IndexType i) {
      int startIdx = outOffsets[i];
      int currCount = 0;
      auto onCandidate = [&](int candidateIdx) -> bool {
        candidates_v[startIdx] = candidateIdx;
        currCount++;
        startIdx++;
        return currCount >= outCounts[i];
      };
      gridQuery.visitCandidates(queryObjs[i], onCandidate);
    });
#else
  outOffsets[0] = 0;
  for(int i = 0; i < qsize; i++)
  {
    outCounts[i] = 0;
    gridQuery.visitCandidates(queryObjs[i], [&](int candidateIdx) {
      outCounts[i]++;
      outCandidates.push_back(candidateIdx);
    });
    if(i + 1 < qsize)
    {
      outOffsets[i + 1] = outOffsets[i] + outCounts[i];
    }
  }
#endif
}

template <int NDIMS, typename ExecSpace, typename IndexType>
AXOM_HOST_DEVICE IndexType
ImplicitGrid<NDIMS, ExecSpace, IndexType>::QueryObject::countCandidates(const SpacePoint& pt) const
{
  if(!m_bb.contains(pt))
  {
    return 0;
  }

  GridCell gridCell = m_lattice.gridCell(pt);

  IndexType ncandidates {0};

  // Note: Need to clamp the upper range of the gridCell
  //       to handle points on the upper boundaries of the bbox
  //       This is valid since we've already ensured that pt is in the bbox.

  for(int idim = 0; idim < NDIMS; idim++)
  {
    gridCell[idim] = axom::utilities::clampUpper(gridCell[idim], m_highestBins[idim]);
  }

  const GridCell cellIdx = gridCell;

  // HACK: we use the underlying word data in the bitsets
  // is it possible to lazy-evaluate whole-bitset operations?
  IndexType minWord, maxWord;
  getWordBounds(cellIdx, minWord, maxWord);
  for(int iword = minWord; iword <= maxWord; iword++)
  {
    BitsetType::Word currWord = ~(BitsetType::Word {0});
    for(int idim = 0; idim < NDIMS; idim++)
    {
      currWord &= m_binData[idim][cellIdx[idim]].data()[iword];
    }
    if(currWord == BitsetType::Word {0})
    {
      continue;
    }
    // currWord now contains the resulting candidacy information for our given point
    ncandidates += axom::utilities::popcount(currWord);
  }
  return ncandidates;
}

template <int NDIMS, typename ExecSpace, typename IndexType>
AXOM_HOST_DEVICE IndexType ImplicitGrid<NDIMS, ExecSpace, IndexType>::QueryObject::countCandidates(
  const SpatialBoundingBox& bbox) const
{
  if(!m_bb.intersectsWith(bbox))
  {
    return 0;
  }

  GridCell lowerCell = m_lattice.gridCell(bbox.getMin());
  GridCell upperCell = m_lattice.gridCell(bbox.getMax());

  for(int idim = 0; idim < NDIMS; idim++)
  {
    // Note: Need to clamp the gridCell ranges since the input box boundaries
    //       are not restricted to the implicit grid's bounding box
    lowerCell[idim] = axom::utilities::clampLower(lowerCell[idim], IndexType {0});
    upperCell[idim] = axom::utilities::clampUpper(upperCell[idim], m_highestBins[idim]);
  }

  const GridCell lowerRange = lowerCell;
  const GridCell upperRange = upperCell;

  IndexType ncandidates {0};

  // HACK: we use the underlying word data in the bitsets
  // is it possible to lazy-evaluate whole-bitset operations?
  IndexType minWord, maxWord;
  getWordBounds(lowerRange, upperRange, minWord, maxWord);
  for(int iword = minWord; iword <= maxWord; iword++)
  {
    BitsetType::Word currWord = ~(BitsetType::Word {0});
    for(int idim = 0; idim < NDIMS; idim++)
    {
      // Compute candidates across all bins for current word
      BitsetType::Word dimWord {0};
      for(int ibin = lowerRange[idim]; ibin <= upperRange[idim]; ibin++)
      {
        dimWord |= m_binData[idim][ibin].data()[iword];
      }
      // Intersect with candidate sets from other dimensions
      currWord &= dimWord;
    }
    if(currWord == BitsetType::Word {0})
    {
      continue;
    }
    // currWord now contains the resulting candidacy information for our given point
    ncandidates += axom::utilities::popcount(currWord);
  }
  return ncandidates;
}

template <int NDIMS, typename ExecSpace, typename IndexType>
template <typename FuncType>
AXOM_HOST_DEVICE void ImplicitGrid<NDIMS, ExecSpace, IndexType>::QueryObject::visitCandidates(
  const SpacePoint& pt,
  FuncType&& candidatePredicate) const
{
  if(!m_bb.contains(pt))
  {
    return;
  }

  GridCell gridCell = m_lattice.gridCell(pt);

  constexpr int bitsPerWord = BitsetType::BitsPerWord;

  // Note: Need to clamp the upper range of the gridCell
  //       to handle points on the upper boundaries of the bbox
  //       This is valid since we've already ensured that pt is in the bbox.

  for(int idim = 0; idim < NDIMS; idim++)
  {
    gridCell[idim] = axom::utilities::clampUpper(gridCell[idim], m_highestBins[idim]);
  }

  const GridCell cellIdx = gridCell;

  // HACK: we use the underlying word data in the bitsets
  // is it possible to lazy-evaluate whole-bitset operations?
  const int nbits = m_binData[0][0].size();
  IndexType minWord, maxWord;
  getWordBounds(cellIdx, minWord, maxWord);
  for(int iword = minWord; iword <= maxWord; iword++)
  {
    BitsetType::Word currWord = ~(BitsetType::Word {0});
    for(int idim = 0; idim < NDIMS; idim++)
    {
      currWord &= m_binData[idim][cellIdx[idim]].data()[iword];
    }
    if(currWord == BitsetType::Word {0})
    {
      continue;
    }
    // currWord now contains the resulting candidacy information for our given point
    const int numBits = axom::utilities::min(bitsPerWord, nbits - (iword * bitsPerWord));
    int currBit = axom::utilities::countr_zero(currWord);
    while(currBit < numBits)
    {
      bool found = getVisitResult(candidatePredicate, iword * bitsPerWord + currBit);
      currBit++;
      currBit += (currBit < numBits) ? axom::utilities::countr_zero(currWord >> currBit) : 0;

      if(found)
      {
        return;
      }
    }
  }
}

template <int NDIMS, typename ExecSpace, typename IndexType>
template <typename FuncType>
AXOM_HOST_DEVICE void ImplicitGrid<NDIMS, ExecSpace, IndexType>::QueryObject::visitCandidates(
  const SpatialBoundingBox& bbox,
  FuncType&& candidatePredicate) const
{
  if(!m_bb.intersectsWith(bbox))
  {
    return;
  }

  GridCell lowerCell = m_lattice.gridCell(bbox.getMin());
  GridCell upperCell = m_lattice.gridCell(bbox.getMax());

  for(int idim = 0; idim < NDIMS; idim++)
  {
    // Note: Need to clamp the gridCell ranges since the input box boundaries
    //       are not restricted to the implicit grid's bounding box
    lowerCell[idim] = axom::utilities::clampLower(lowerCell[idim], IndexType {0});
    upperCell[idim] = axom::utilities::clampUpper(upperCell[idim], m_highestBins[idim]);
  }

  const GridCell lowerRange = lowerCell;
  const GridCell upperRange = upperCell;

  constexpr int bitsPerWord = BitsetType::BitsPerWord;

  // HACK: we use the underlying word data in the bitsets
  // is it possible to lazy-evaluate whole-bitset operations?
  const int nbits = m_binData[0][0].size();
  IndexType minWord, maxWord;
  getWordBounds(lowerRange, upperRange, minWord, maxWord);
  for(int iword = minWord; iword <= maxWord; iword++)
  {
    BitsetType::Word currWord = ~(BitsetType::Word {0});
    for(int idim = 0; idim < NDIMS; idim++)
    {
      // Compute candidates across all bins for current word
      BitsetType::Word dimWord {0};
      for(int ibin = lowerRange[idim]; ibin <= upperRange[idim]; ibin++)
      {
        dimWord |= m_binData[idim][ibin].data()[iword];
      }
      // Intersect with candidate sets from other dimensions
      currWord &= dimWord;
    }
    if(currWord == BitsetType::Word {0})
    {
      continue;
    }
    // currWord now contains the resulting candidacy information for our given point
    const int numBits = axom::utilities::min(bitsPerWord, nbits - (iword * bitsPerWord));
    int currBit = axom::utilities::countr_zero(currWord);
    while(currBit < numBits)
    {
      bool found = getVisitResult(candidatePredicate, iword * bitsPerWord + currBit);
      currBit++;
      currBit += (currBit < numBits) ? axom::utilities::countr_zero(currWord >> currBit) : 0;

      if(found)
      {
        return;
      }
    }
  }
}

}  // namespace spin
}  // namespace axom

#endif  // AXOM_SPIN_IMPLICIT_GRID__HPP_
