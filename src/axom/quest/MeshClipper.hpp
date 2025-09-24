// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_MESHCLIPPER_HPP
#define AXOM_QUEST_MESHCLIPPER_HPP

#include "axom/config.hpp"

#include "axom/klee/Geometry.hpp"
#include "axom/quest/MeshClipperStrategy.hpp"
#include "axom/quest/ShapeeMesh.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

/*!
 * @brief Abstract base class encapsulating the clipping task in
 * intersection shaping, to be specialized for geometry-specific
 * implementations.
 *
 * Each object operates is associated with a single klee::Geometry object.*
 *
 * This class defines the interfaces for specializing computations
 * that can run fast on specific geometries.  This allows for
 * geometry-specific optimizations.  If no such specialization
 * is provided, default methods are used.
 */
class MeshClipper
{
public:
  //!@brief Whether an element is in, out or on shape boundary.
  using LabelType = MeshClipperStrategy::LabelType;

  static constexpr axom::IndexType TETS_PER_HEXAHEDRON = MeshClipperStrategy::TETS_PER_HEXAHEDRON;

  /*!
   * @brief Construct a shape clipper
   *
   * @param [in/out] bpMesh Single-domain Blueprint mesh
   *   to shape into.
   * @param [in] strategy Strategy where external
   *   shape-dependent operations are implemented.
   *
   * @c bpMesh must be an unstructured hex mesh.
   * That is the only type currently supported.
   */
  MeshClipper(quest::experimental::ShapeeMesh& shapeeMesh, const std::shared_ptr<MeshClipperStrategy>& strategy);

  //!@brief The mesh.
  ShapeeMesh& getShapeeMesh() { return m_shapeeMesh; }

  //!@brief Allocator id to be used for all array data.
  int getAllocatorID() const { return m_shapeeMesh.getAllocatorID(); }

  void setVerbose(bool verbose) { m_verbose = verbose; }

  /*!
   * @brief Clip
   *
   * @param ovlap [out] Shape overlap volume of each cell
   *   in the shapee mesh.
   */
  void clip(axom::Array<double>& ovlap);

  /*!
   * @brief Clip
   *
   * @param ovlap [out] Shape overlap volume of each cell
   *   in the shapee mesh.
   */
  void clip(axom::ArrayView<double> ovlap);

  //!@brief Dimension of the shape (2 or 3)
  int dimension() const { return m_shapeeMesh.dimension(); }

  /*!
   * @brief Single interface for methods implemented with
   * execution space templates.
   *
   * These methods require messy instantiation of
   * execution spaces and their runtime selection.
   *
   * The implementations are in class detail::MeshClipperImpl,
   * which is templated on execution space.
   */
  struct Impl
  {
    Impl(MeshClipper& impl) : m_myClipper(impl) { }
    virtual ~Impl() = default;

    static constexpr axom::IndexType TETS_PER_HEXAHEDRON = MeshClipperStrategy::TETS_PER_HEXAHEDRON;

    /*!
     * @brief Initialize overlap volumes to full for cells completely
     * inside the shape and zero for cells outside or on shape boundary.
     */
    virtual void initVolumeOverlaps(const axom::ArrayView<MeshClipperStrategy::LabelType>& labels,
                                    axom::ArrayView<double> ovlap) = 0;

    //! @brief Initialize overlap volumes to zero.
    virtual void initVolumeOverlaps(axom::ArrayView<double> ovlap) = 0;

    //!@brief Collect unlabeled LABEL_ON indices into an index list.
    virtual void collectOnIndices(const axom::ArrayView<LabelType>& labels,
                                  axom::Array<axom::IndexType>& onIndices) = 0;

    //!@brief Change tet indices from sparse-set values to full-set values.
    virtual void remapTetIndices(axom::ArrayView<axom::IndexType> tetsOnBdry,
                                 axom::ArrayView<const axom::IndexType> cellsOnBdry) = 0;

    //!@brief Add volumes of tets inside the geometry to the volume data.
    virtual void addVolumesOfInteriorTets(axom::ArrayView<const axom::IndexType> cellsOnBdry,
                                          axom::ArrayView<const LabelType> tetLabels,
                                          axom::ArrayView<double> ovlap) = 0;

    //!@brief Compute clip volumes for every cell.
    virtual void computeClipVolumes3D(axom::ArrayView<double> ovlap) = 0;

    //!@brief Compute clip volumes for cell in an index list.
    virtual void computeClipVolumes3D(const axom::ArrayView<axom::IndexType>& cellIndices,
                                      axom::ArrayView<double> ovlap) = 0;

    /*!
     * @brief Compute clip volumes for cell tets in an index list.
     *
     * The tets are the results from decomposing each cell hex into
     * TETS_PER_HEXAHEDRON tets and stored consecutively.
     */
    virtual void computeClipVolumes3DTets(const axom::ArrayView<axom::IndexType>& tetIndices,
                                          axom::ArrayView<double> ovlap) = 0;

    //!@brief Count the number of labels of each type.
    virtual void getLabelCounts(axom::ArrayView<const LabelType> labels,
                                axom::IndexType& inCount,
                                axom::IndexType& onCount,
                                axom::IndexType& outCount) = 0;

    ShapeeMesh& getShapeeMesh() { return m_myClipper.m_shapeeMesh; }

    MeshClipperStrategy& getStrategy() { return *m_myClipper.m_strategy; }

  private:
    //!@brief The MeshClipper that owns this Impl.
    MeshClipper& m_myClipper;
  };

  //! @brief For assessments, not general use.
  void getClippingStats(axom::IndexType& localCellInCount,
                        axom::IndexType& globalCellInCount,
                        axom::IndexType& maxLocalCellInCount) const;

private:
  friend Impl;

  quest::experimental::ShapeeMesh& m_shapeeMesh;

  //! @brief Shape-specific operations in clipping.
  std::shared_ptr<quest::experimental::MeshClipperStrategy> m_strategy;

  //! @brief Object where execution space code is instantiated.
  std::unique_ptr<Impl> m_impl;

  /* NOTE: MeshClipperStrategy is for shape-specific functions,
   * implemented externally.  Impl implements internal algorithms
   * for multiple execution spaces.
   */

  ///@{
  //! @name Statistics
  axom::IndexType m_localCellInCount {0};
  ///@}

  bool m_verbose;

#if defined(__CUDACC__)
public:
#endif
  //!@brief Allocate a delegate for m_shapeeMesh's runtime policy.
  std::unique_ptr<Impl> newImpl();

  //@{
  //!@name Convenience methods
  //!@brief Count the number of labels of each type.
  void getLabelCounts(const axom::Array<LabelType>& labels,
                      axom::IndexType& inCount,
                      axom::IndexType& onCount,
                      axom::IndexType& outCount)
  {
    m_impl->getLabelCounts(labels, inCount, onCount, outCount);
  }

  void logLabelStats(axom::ArrayView<const LabelType> labels, const std::string& labelType);
  //@}
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_MESHCLIPPER_HPP
