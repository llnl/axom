// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_GEOMETRYCLIPPER_HPP
#define AXOM_QUEST_GEOMETRYCLIPPER_HPP

#include "axom/config.hpp"

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/quest/ShapeeMesh.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief Abstract base class encapsulating the clipping task in
  intersection shaping, to be specialized for geometry-specific
  implementations.

  Each object operates is associated with a single klee::Geometry object.

  This class defines the interfaces for specializing computations
  that can run fast on specific geometries.  This allows for
  geometry-specific optimizations.  If no such specialization
  is provided, default methods are used.
*/
class GeometryClipper
{
public:
  //!@brief Whether an element as in, out or on shape boundary.
  using LabelType = GeometryClipperStrategy::LabelType;

  /*!
    @brief Construct a shape clipper

    @param [in/out] bpMesh Single-domain Blueprint mesh
      to shape into.
    @param [in] strategy Strategy where external
      shape-dependent operations are implemented.

    @c bpMesh must be an unstructured hex mesh.
    That is the only type currently supported.
  */
  GeometryClipper(quest::ShapeeMesh& shapeeMesh,
                  const std::shared_ptr<GeometryClipperStrategy>& strategy);

  //!@brief The mesh.
  ShapeeMesh& getShapeeMesh() { return m_shapeeMesh; }

  //!@brief Allocator id to be used for all array data.
  int getAllocatorID() const { return m_shapeeMesh.getAllocatorID(); }

  void setVerbose(bool verbose) { m_verbose = verbose; }

  /*!
    @brief Clip

    @param ovlap [out] Shape overlap volume of each cell
      in the shapee mesh.
  */
  void clip(axom::Array<double>& ovlap);

  /*!
    @brief Clip

    @param ovlap [out] Shape overlap volume of each cell
      in the shapee mesh.
  */
  void clip(axom::ArrayView<double> ovlap);

  //!@brief Dimension of the shape (2 or 3)
  int dimension() const { return m_shapeeMesh.dimension(); }

  /*!
    @brief Single interface for some methods delegated out of
    GeometryClipper

    Delegated methods are those requiring messy instantiation of
    execution spaces and their runtime selection.

    The implementations are in class detail::GeometryClipperDelegateExec,
    which is templated on execution space.
  */
  struct Delegate
  {
    Delegate(GeometryClipper& delegator) : m_delegator(delegator) { }
    virtual ~Delegate() = default;

    /*!
      @brief Initialize overlap volumes to full for cells completely
      inside the shape and zero for cells outside or on shape boundary.
    */
    virtual void initVolumeOverlaps(
      const axom::ArrayView<GeometryClipperStrategy::LabelType>& labels,
      axom::ArrayView<double> ovlap) = 0;

    //!@brief Collect unlabeled cells indices into an index list.
    virtual void collectUnlabeledCellIndices(const axom::ArrayView<LabelType>& labels,
                                             axom::Array<axom::IndexType>& unlabeledCells) = 0;

    //!@brief Compute clip volumes for every cell.
    virtual void computeClipVolumes3D(axom::ArrayView<double> ovlap) = 0;

    //!@brief Compute clip volumes for cell in an index list.
    virtual void computeClipVolumes3D(const axom::ArrayView<axom::IndexType>& cellIndices,
                                      axom::ArrayView<double> ovlap) = 0;

    //!@brief Delegate for getLabelCounts.
    virtual void getLabelCounts(axom::ArrayView<const LabelType> labels,
                                axom::IndexType& inCount,
                                axom::IndexType& onCount,
                                axom::IndexType& outCount) = 0;

    ShapeeMesh& getShapeeMesh() { return m_delegator.m_shapeeMesh; }

    GeometryClipperStrategy& getStrategy() { return *m_delegator.m_strategy; }

    GeometryClipper& getDelegator() { return m_delegator; }

  private:
    GeometryClipper& m_delegator;
  };

private:
  friend Delegate;

  quest::ShapeeMesh& m_shapeeMesh;

  //! @brief Shape-specific operations in clipping.
  std::shared_ptr<quest::GeometryClipperStrategy> m_strategy;

  //! @brief Delegate object handling execution space templates.
  std::unique_ptr<Delegate> m_delegate;

  /* NOTE: GeometryClipperStrategy is for shape-specific functions,
     implemented externally.  Delegate implements internal algorithms
     for multiple execution spaces.
     TODO: Change delegate term to policy.
  */

  bool m_verbose;

#if defined(__CUDACC__)
public:
#endif
  //!@brief Allocate a delegate for m_shapeeMesh's runtime policy.
  std::unique_ptr<Delegate> newDelegate();

  //@{
  //!@name Convenience methods
  void getLabelCounts(const axom::Array<LabelType>& labels,
                      axom::IndexType& inCount,
                      axom::IndexType& onCount,
                      axom::IndexType& outCount)
  {
    m_delegate->getLabelCounts(labels, inCount, onCount, outCount);
  }
  //@}
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_GEOMETRYCLIPPER_HPP
