// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_H_
#define AXOM_SPIN_BVH_H_

// axom core includes
#include "axom/config.hpp"                               // for Axom compile-time definitions
#include "axom/core/Macros.hpp"                          // for Axom macros
#include "axom/core/Types.hpp"                           // for axom::IndexType
#include "axom/core/numerics/floating_point_limits.hpp"  // floating_point_limits

#include "axom/core/execution/execution_space.hpp"  // for execution spaces

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"

#include "axom/primal/operators/intersect.hpp"  // for detail::intersect_ray()

#include "axom/spin/policy/LinearBVH.hpp"

// slic includes
#include "axom/slic/interface/slic.hpp"  // for SLIC macros

// C/C++ includes
#include <type_traits>  // for std::is_floating_point(), std::is_same()
#include <memory>

namespace axom
{
namespace spin
{
/*!
 * \brief Enumerates the list of return codes for various BVH operations.
 */
enum BVHReturnCodes
{
  BVH_BUILD_FAILED = -1,  //!< indicates that generation of the BVH failed
  BVH_BUILD_OK,           //!< indicates that the BVH was generated successfully
};

enum class BVHType
{
  LinearBVH
};

template <typename FloatType, int NDIMS, typename ExecType, BVHType Policy>
struct BVHPolicy;

template <typename FloatType, int NDIMS, typename ExecType>
struct BVHPolicy<FloatType, NDIMS, ExecType, BVHType::LinearBVH>
{
  using ImplType = policy::LinearBVH<FloatType, NDIMS, ExecType>;
};

/*!
 * \class BVH
 *
 * \brief Defines a Bounding Volume Hierarchy (BVH) spatial acceleration
 *  data structure over a set of geometric entities.
 *
 * The BVH class provides functionality for generating a hierarchical spatial
 * partitioning over a set of geometric entities. Each entity in the BVH is
 * represented by a bounding volume, in this case an axis-aligned bounding box.
 * Once the BVH structure is generated, it is used to accelerate various spatial
 * queries, such as, collision detection, ray tracing, etc., by reducing the
 * search space for a given operation to an abbreviated list of candidate
 * geometric entities to check for a particular query.
 *
 * \tparam NDIMS the number of dimensions, e.g., 2 or 3.
 * \tparam ExecSpace the execution space to use, e.g. SEQ_EXEC, CUDA_EXEC, etc.
 * \tparam FloatType floating precision, e.g., `double` or `float`. Optional.
 *
 * \note The last template parameter is optional. Defaults to double precision
 *  if not specified.
 *
 * \pre The spin::BVH class requires RAJA and Umpire with CUDA_EXEC.
 *
 * \note The Execution Space, supplied as the 2nd template argument, specifies
 *
 *  1. Where and how the BVH is generated and stored
 *  2. Where and how subsequent queries are performed
 *  3. The default memory space, bound to the corresponding execution space
 *
 * \see axom::execution_space for more details.
 *
 *  A simple example illustrating how to use the BVH class is given below:
 *  \code
 *
 *     namespace spin = axom::spin;
 *     constexpr int DIMENSION = 3;
 *
 *     // get a list of axis-aligned bounding boxes in a flat array
 *     const primal::BoundingBox<float, DIMENSION>* aabbs = ...
 *
 *     // create a 3D BVH instance in parallel on the CPU using OpenMP
 *     spin::BVH< DIMENSION, axom::OMP_EXEC > bvh;
 *     bvh.initialize( aabbs, numItems );
 *
 *     // query points supplied in arrays, qx, qy, qz
 *     const axom::IndexType numPoints = ...
 *     const double* qx = ...
 *     const double* qy = ...
 *     const double* qz = ...
 *     //use the ZipPoint class to tie them together
 *     ZipPoint qpts {{qx,qy,qz}};
 *
 *     // output array buffers
 *     axom::IndexType* offsets    = axom::allocate< IndexType >( numPoints );
 *     axom::IndexType* counts     = axom::allocate< IndexType >( numPoints );
 *     axom::IndexType* candidates = nullptr;
 *
 *     // find candidates in parallel, allocates and populates the supplied
 *     // candidates array
 *     bvh.findPoints( offsets, counts, candidates, numPoints, qpts );
 *     SLIC_ASSERT( candidates != nullptr );
 *
 *     ...
 *
 *     // caller is responsible for properly de-allocating the candidates array
 *     axom::deallocate( candidates );
 *
 *  \endcode
 *  \accelerated
 */
template <int NDIMS, typename ExecSpace = axom::SEQ_EXEC, typename FloatType = double, BVHType BVHImpl = BVHType::LinearBVH>
class BVH
{
private:
  // compile time checks
  AXOM_STATIC_ASSERT_MSG(((NDIMS == 2) || (NDIMS == 3)),
                         "The BVH class may be used only in 2D or 3D.");
  AXOM_STATIC_ASSERT_MSG(std::is_floating_point<FloatType>::value,
                         "A valid FloatingType must be used for the BVH.");
  AXOM_STATIC_ASSERT_MSG(axom::execution_space<ExecSpace>::valid(),
                         "A valid execution space must be supplied to the BVH.");

  using ImplType = typename BVHPolicy<FloatType, NDIMS, ExecSpace, BVHImpl>::ImplType;

  template <typename It>
  class IteratorTraits
  {
  private:
    template <typename U, typename Ret = decltype(std::declval<U&>()[0])>
    static Ret array_operator_type(U AXOM_UNUSED_PARAM(obj))
    { }

    static std::false_type array_operator_type(...) { return std::false_type {}; }

  public:
    // The base object of the return type of a call to iter[], or std::false_type
    // if operator[] does not exist.
    using BaseType = typename std::decay<decltype(array_operator_type(std::declval<It>()))>::type;

    // The iterator must be an array-like type.
    static_assert(!std::is_same<BaseType, std::false_type>::value,
                  "Iterator type must be accessible with the array operator[].");

    // The iterator object is copy-captured in a lambda and thus needs to be
    // copy-constructible.
    AXOM_STATIC_ASSERT_MSG(std::is_copy_constructible<It>::value,
                           "Iterator type must be copy constructible.");

    // We need to be able to copy-construct a primitve BaseType within the
    // traversal kernel.
    AXOM_STATIC_ASSERT_MSG(std::is_copy_constructible<BaseType>::value,
                           "Iterator's value type must be copy constructible.");
  };

public:
  using BoxType = typename primal::BoundingBox<FloatType, NDIMS>;
  using PointType = typename primal::Point<FloatType, NDIMS>;
  using RayType = typename primal::Ray<FloatType, NDIMS>;

  using TraverserType = typename ImplType::TraverserType;
  using ExecSpaceType = ExecSpace;

public:
  /*!
   * \brief Default constructor.
   */
  BVH() : m_AllocatorID(axom::execution_space<ExecSpace>::allocatorID()) { }

  /*!
   * \brief Constructs a BVH instance of specified dimension over a given set
   *  of geometric entities, each represented by its corresponding axis-aligned
   *  bounding box.
   *
   * \param [in] boxes buffer consisting of bounding boxes for each entity.
   * \param [in] numItems the total number of items to store in the BVH.
   * \param [in] allocatorID the ID of the allocator to use in BVH construction.
   * \param [in] scaleFactor the scale factor for the given bounding boxes.
   * \param [in] tolerance the tolerance to use when querying the BVH.
   *
   * \warning The supplied boxes array must point to a buffer in a memory space
   *  that is compatible with the execution space. For example, when using
   *  CUDA_EXEC, boxes must be in unified memory or GPU memory. The supplied
   *  allocator ID must also correspond to an umpire allocator and be compatible
   *  with the execution space. The code currently does not check for these
   *  conditions.
   *
   * \pre boxes != nullptr
   * \pre numItems > 0
   * \pre allocatorID must correspond to an Umpire allocator and be compatible
   *      with the execution space
   */
  template <typename BoxIndexable>
  BVH(const BoxIndexable boxes,
      IndexType numItems,
      int allocatorID = axom::execution_space<ExecSpace>::allocatorID(),
      FloatType tolerance = DEFAULT_TOLERANCE,
      FloatType scaleFactor = DEFAULT_SCALE_FACTOR)
    : m_AllocatorID {allocatorID}
    , m_tolerance {tolerance}
    , m_scaleFactor {scaleFactor}
  {
    initialize(boxes, numItems);
  }

  /*!
   * \brief Initializes a BVH instance of specified dimension over a given set
   *  of geometric entities, each represented by its corresponding axis-aligned
   *  bounding box.
   *
   * \param [in] boxes buffer consisting of bounding boxes for each entity.
   * \param [in] numItems the total number of items to store in the BVH.
   *
   * \return status set to BVH_BUILD_OK on success.
   *
   * \note If an allocatorID has not been set in a call to setAllocatorID(),
   *  the code will use the default allocator ID for the execution space
   *  specified via axom::execution_space<ExecSpace>::allocatorID() when the
   *  BVH object is instantiated.
   *
   * \warning The supplied boxes array must point to a buffer in a memory space
   *  that is compatible with the execution space. For example, when using
   *  CUDA_EXEC, boxes must be in unified memory or GPU memory. The code
   *  currently does not check for that.
   *
   * \pre boxes != nullptr
   * \pre numItems > 0
   */
  template <typename BoxIndexable>
  int initialize(const BoxIndexable boxes, IndexType numItems);

  bool isInitialized() const { return m_bvh != nullptr; }

  /*!
   * \brief Sets the ID of the allocator used by the BVH.
   * \param [in] allocatorID the ID of the allocator to use in BVH construction
   */
  void setAllocatorID(int allocatorID) { m_AllocatorID = allocatorID; };

  /*!
   * \brief Get the ID of the allocator used by the BVH.
   * \return allocatorID the ID of the allocator used by the BVH.
   */
  int getAllocatorID() const { return m_AllocatorID; };

  /*!
   * \brief Sets the scale factor for scaling the supplied bounding boxes.
   * \param [in] scale_factor the scale factor
   *
   * \note The default scale factor is set to 1.000123
   */
  void setScaleFactor(FloatType scale_factor) { m_scaleFactor = scale_factor; };

  /*!
   * \brief Returns the scale factor used when constructing the BVH.
   * \return scale_factor the scale factor
   */
  FloatType getScaleFactor() const { return m_scaleFactor; };

  /*!
   * \brief Sets the tolerance used for querying the BVH.
   * \param [in] TOL the tolerance to use.
   *
   * \note Default tolerance set to floating_point_limits<FloatType>::epsilon()
   */
  void setTolerance(FloatType EPS) { m_tolerance = EPS; };

  /*!
   * \brief Returns the tolerance value used for BVH queries.
   * \return TOL the tolerance
   */
  FloatType getTolerance() const { return m_tolerance; };

  /*!
   * \brief Returns the bounds of the BVH, given by the the root bounding box.
   *
   * \return box the bounding box of the constructed BVH, or an invalid
   *  bounding box if the BVH has not been initialized.
   */
  BoxType getBounds() const
  {
    if(m_bvh)
    {
      return m_bvh->getBoundsImpl();
    }
    else
    {
      // return invalid box
      return BoxType {};
    }
  }

  /*!
   * \brief Returns a device-copyable object that can be used to traverse the
   *  BVH from inside a device kernel.
   *
   * \return it the traverser object for the current BVH.
   *
   * \node The traverser object may only be used in the same execution space as
   *  the one that the BVH class was instantiated with.
   */
  TraverserType getTraverser() const { return m_bvh->getTraverserImpl(); }

  /*!
   * \brief Finds the candidate bins that contain each of the query points.
   *
   * \param [out] offsets offset to the candidates array for each query point
   * \param [out] counts stores the number of candidates per query point
   * \param [out] candidates array of the candidate IDs for each query point
   * \param [in]  numPts the total number of query points supplied
   * \param [in]  points array of points to query against the BVH
   *
   * \note Upon completion, the ith query point has:
   *  * counts[ i ] candidates
   *  * Stored in the candidates array in the following range:
   *    [ offsets[ i ], offsets[ i ]+counts[ i ] ]
   *  * The sum of all counts is the size of the candidates array,
   *    candidates.size()
   *
   * \pre offsets.size() == numPts
   * \pre counts.size()  == numPts
   * \pre points != nullptr
   */
  template <typename PointIndexable>
  void findPoints(axom::ArrayView<IndexType> offsets,
                  axom::ArrayView<IndexType> counts,
                  axom::Array<IndexType>& candidates,
                  IndexType numPts,
                  PointIndexable points) const;

  /*!
   * \brief Finds the candidate bins that intersect the given rays.
   *
   * \param [out] offsets offset to the candidates array for each ray
   * \param [out] counts stores the number of candidates for each ray
   * \param [out] candidates array of candidate IDs for each ray
   * \param [in] numRays the total number of rays
   * \param [in] rays array of the rays to query against the BVH
   *
   * \note After the call to findRays(), the ith ray has:
   *  * counts[ i ] candidates
   *  * candidates stored in [ offsets[ i ], offsets[i]+counts[i] ]
   *  * The sum of all counts is the size of the candidates array,
   *    candidates.size()
   *
   * \pre offsets.size() == numRays
   * \pre counts.size()  == numRays
   * \pre rays != nullptr
   */
  template <typename RayIndexable>
  void findRays(axom::ArrayView<IndexType> offsets,
                axom::ArrayView<IndexType> counts,
                axom::Array<IndexType>& candidates,
                IndexType numRays,
                RayIndexable rays) const;

  /*!
   * \brief Finds the candidate bins that intersect the given bounding boxes.
   *
   * \param [out] offsets offset to the candidates array for each bounding box
   * \param [out] counts stores the number of candidates for each bounding box
   * \param [out] candidates array of candidate IDs for each bounding box
   * \param [in]  numBoxes the total number of bounding boxes
   * \param [in]  boxes array of boxes to query against the BVH
   *
   * \note After the call to findBoundingBoxes(), the ith bounding box has:
   *  * counts[ i ] candidates
   *  * candidates stored in [ offsets[ i ], offsets[i]+counts[i] ]
   *  * The sum of all counts is the size of the candidates array,
   *    candidates.size()
   *
   * \pre offsets.size() == numBoxes
   * \pre counts.size()  == numBoxes
   * \pre boxes != nullptr
   */
  template <typename BoxIndexable>
  void findBoundingBoxes(axom::ArrayView<IndexType> offsets,
                         axom::ArrayView<IndexType> counts,
                         axom::Array<IndexType>& candidates,
                         IndexType numBoxes,
                         BoxIndexable boxes) const;

  /*!
   * \brief Writes the BVH to the specified VTK file for visualization.
   * \param [in] fileName the name of VTK file.
   * \note Primarily used for debugging.
   */
  void writeVtkFile(const std::string& fileName) const;

private:
  /// \name Private Members
  /// @{
  static constexpr FloatType DEFAULT_SCALE_FACTOR = 1.000123;
  static constexpr FloatType DEFAULT_TOLERANCE =
    axom::numerics::floating_point_limits<FloatType>::epsilon();

  int m_AllocatorID;
  FloatType m_tolerance {DEFAULT_TOLERANCE};
  FloatType m_scaleFactor {DEFAULT_SCALE_FACTOR};
  std::unique_ptr<ImplType> m_bvh {};
  /// @}
};

//------------------------------------------------------------------------------
//  BVH implementation
//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
template <typename BoxIndexable>
int BVH<NDIMS, ExecSpace, FloatType, Impl>::initialize(const BoxIndexable boxes, IndexType numBoxes)
{
  AXOM_ANNOTATE_SCOPE("BVH::initialize");

  using IterBase = typename IteratorTraits<BoxIndexable>::BaseType;

  // Ensure that the iterator returns objects convertible to primal::BoundingBox.
  static_assert(std::is_convertible<IterBase, BoxType>::value,
                "Iterator must return objects convertible to primal::BoundingBox.");

  // STEP 1: Allocate a BVH, potentially deleting the existing BVH if it exists
  m_bvh.reset(new ImplType);

  // STEP 1: Handle case when user supplied 0 or 1 bounding boxes.
  BoxType* boxesptr = nullptr;
  if(numBoxes <= 1)
  {
    const bool copyFirst = numBoxes == 1;
    numBoxes = 2;
    boxesptr = axom::allocate<BoxType>(numBoxes, m_AllocatorID);

    // copy first box and add a fake 2nd box
    for_all<ExecSpace>(
      2,
      AXOM_LAMBDA(IndexType i) {
        if(copyFirst && i == 0)
        {
          boxesptr[i] = boxes[i];
        }
        else
        {
          BoxType empty_box;
          // Make the box invalid.
          empty_box.clear();
          boxesptr[i] = empty_box;
        }
      });
    m_bvh->buildImpl(boxesptr, numBoxes, m_scaleFactor, m_AllocatorID);
  }
  else
  {
    m_bvh->buildImpl(boxes, numBoxes, m_scaleFactor, m_AllocatorID);
  }

  // STEP 5: deallocate boxesptr if user supplied a single box
  if(boxesptr)
  {
    SLIC_ASSERT(boxesptr != nullptr);
    axom::deallocate(boxesptr);
  }
  return BVH_BUILD_OK;
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
template <typename PointIndexable>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::findPoints(axom::ArrayView<IndexType> offsets,
                                                        axom::ArrayView<IndexType> counts,
                                                        axom::Array<IndexType>& candidates,
                                                        IndexType numPts,
                                                        PointIndexable pts) const
{
  AXOM_ANNOTATE_SCOPE("BVH::findPoints");

  using IterBase = typename IteratorTraits<PointIndexable>::BaseType;

  // Ensure that the iterator returns objects convertible to primal::Point.
  static_assert(std::is_convertible<IterBase, PointType>::value,
                "Iterator must return objects convertible to primal::Point.");

  SLIC_ASSERT(m_bvh != nullptr);

  // Define traversal predicates
  auto predicate = [=] AXOM_HOST_DEVICE(const PointType& p, const BoxType& bb) -> bool {
    return bb.contains(p);
  };

  candidates =
    m_bvh->template findCandidatesImpl<PointType>(predicate, offsets, counts, numPts, pts, m_AllocatorID);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
template <typename RayIndexable>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::findRays(axom::ArrayView<IndexType> offsets,
                                                      axom::ArrayView<IndexType> counts,
                                                      axom::Array<IndexType>& candidates,
                                                      IndexType numRays,
                                                      RayIndexable rays) const
{
  AXOM_ANNOTATE_SCOPE("BVH::findRays");

  using IterBase = typename IteratorTraits<RayIndexable>::BaseType;

  // Ensure that the iterator returns objects convertible to primal::Ray.
  static_assert(std::is_convertible<IterBase, RayType>::value,
                "Iterator must return objects convertible to primal::Ray.");

  SLIC_ASSERT(m_bvh != nullptr);

  const FloatType TOL = m_tolerance;

  // Define traversal predicates
  auto predicate = [=] AXOM_HOST_DEVICE(const RayType& r, const BoxType& bb) -> bool {
    primal::Point<FloatType, NDIMS> tmp;
    return primal::detail::intersect_ray(r, bb, tmp, TOL);
  };

  candidates =
    m_bvh->template findCandidatesImpl<RayType>(predicate, offsets, counts, numRays, rays, m_AllocatorID);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
template <typename BoxIndexable>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::findBoundingBoxes(axom::ArrayView<IndexType> offsets,
                                                               axom::ArrayView<IndexType> counts,
                                                               axom::Array<IndexType>& candidates,
                                                               IndexType numBoxes,
                                                               BoxIndexable boxes) const
{
  AXOM_ANNOTATE_SCOPE("BVH::findBoundingBoxes");

  using IterBase = typename IteratorTraits<BoxIndexable>::BaseType;

  // Ensure that the iterator returns objects convertible to primal::BoundingBox.
  static_assert(std::is_convertible<IterBase, BoxType>::value,
                "Iterator must return objects convertible to primal::BoundingBox.");

  SLIC_ASSERT(m_bvh != nullptr);

  // STEP 2: define traversal predicates
  auto predicate = [=] AXOM_HOST_DEVICE(const BoxType& bb1, const BoxType& bb2) -> bool {
    return bb1.intersectsWith(bb2);
  };

  candidates = m_bvh->template findCandidatesImpl<BoxType>(predicate,
                                                           offsets,
                                                           counts,
                                                           numBoxes,
                                                           boxes,
                                                           m_AllocatorID);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType, BVHType Impl>
void BVH<NDIMS, ExecSpace, FloatType, Impl>::writeVtkFile(const std::string& fileName) const
{
  SLIC_ASSERT(m_bvh != nullptr);

  m_bvh->writeVtkFileImpl(fileName);
}

}  // namespace spin
}  // namespace axom

#endif  // AXOM_SPIN_BVH_H_
