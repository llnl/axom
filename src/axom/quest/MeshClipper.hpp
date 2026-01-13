// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_MESHCLIPPER_HPP
#define AXOM_QUEST_MESHCLIPPER_HPP

#include "axom/config.hpp"

#include "axom/klee/Geometry.hpp"
#include "axom/quest/MeshClipperStrategy.hpp"
#include "axom/quest/ShapeMesh.hpp"
#include "conduit/conduit_node.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

/*!
 * @brief Entry point for computing clipping a computational mesh
 * by overlaying a geometry and computing the intersection volume
 * the geometry makes with each mesh cell.
 *
 * To construct:
 * - Wrap the computational mesh in a ShapeMesh object.
 * - Allocate a MeshClipperStrategy implementation to provide
 *   geometry-specific operations.  Axom natively provides
 *   implementations for some common geometries.
 *
 * To compute the intersection volumes, use one of the clip() methods.
 */
class MeshClipper
{
public:
  //!@brief Whether an element is in, out or on shape boundary.
  using LabelType = MeshClipperStrategy::LabelType;

  static constexpr axom::IndexType NUM_TETS_PER_HEX = ShapeMesh::NUM_TETS_PER_HEX;

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
  MeshClipper(quest::experimental::ShapeMesh& shapeMesh,
              const std::shared_ptr<MeshClipperStrategy>& strategy);

  //!@brief The mesh.
  ShapeMesh& getShapeMesh() { return m_shapeMesh; }

  //!@brief Allocator id to be used for all array data.
  int getAllocatorID() const { return m_shapeMesh.getAllocatorID(); }

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
  int dimension() const { return m_shapeMesh.dimension(); }

  //@{

  /*!
   * @brief Log clipping statistics.
   * Intended for developer use.
   *
   * @param local [in] Log the MPI local statistics.
   * @param sum [in] Log the statistics summed over MPI ranks.
   * @param max [in] Log the max statistics over MPI ranks.
   *
   * This is a collective method if MPI-parallel.
   */
  void logClippingStats(bool local = false, bool sum = true, bool max = false) const;

  /*!
   * @brief Get local assorted clipping statistics,
   * intended for developer use.
   */
  const conduit::Node& getClippingStats() const { return m_counterStats; }

  /*!
   * @brief Get global assorted clipping statistics,
   * intended for developer use.
   *
   * This is a collective method if MPI-parallel.
   */
  conduit::Node getGlobalClippingStats() const;

  /*!
   * @brief Set the level of screening,
   * intended for developer use.
   * @param screenLevel [in] the value to control how much screening
   *   is done to avoid using the slower clipping functions.
   *
   * Screenlevel values:
   * - 0 allows no screening
   * - 1 allows screening cells.  @see MeshClipperStrategy::labelCellsInOut.
   * - 2 also allows screening tets.  @see MeshClipperStrategy::labelTetsInOut.
   * - 3 also allows a check on the tet from mesh and a piece of the discretized
   *     geometry to (conservatively) determine whether the piece is completely
   *     inside or completely outside of the tst.
   */
  void setScreenLevel(int screenLevel) { m_screenLevel = screenLevel; }

  /*!
   * @brief Get the level of screening,
   * intended for developer use.
   */
  int getScreenLevel() const { return m_screenLevel; }

  /*!
   * @brief Add new stats to current stats,
   * intended for developer use.
   */
  static void accumulateClippingStats(conduit::Node& curStats, const conduit::Node& newStats);

  //@}

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

    static constexpr axom::IndexType NUM_TETS_PER_HEX = MeshClipperStrategy::NUM_TETS_PER_HEX;

    /*!
     * @brief Initialize overlap volumes to full for cells completely
     * inside the shape and zero for cells outside or on shape boundary.
     */
    virtual void initVolumeOverlaps(const axom::ArrayView<MeshClipperStrategy::LabelType>& labels,
                                    axom::ArrayView<double> ovlap) = 0;

    //! @brief Initialize overlap volumes to zero.
    virtual void zeroVolumeOverlaps(axom::ArrayView<double> ovlap) = 0;

    //!@brief Collect unlabeled LABEL_ON indices into an index list.
    virtual void collectOnIndices(const axom::ArrayView<LabelType>& labels,
                                  axom::Array<axom::IndexType>& onIndices) = 0;

    /*!
     * @brief Remap tet indices by de-referencing the cell indices they refer to.
     *
     * @param cellIndices [in] a list of cell indices.
     * @param tetIndices [in,out] a list of tet indices.
     *
     * On input, tetIndices have values in [0, cellIndices.size()*NUM_TETS_PER_HEX),
     * as though the cells have indices in [0, cellIndices.size()).
     *
     * On output, tetIndices values are remapped to match real cell indices.
     * For example, tet index values in
     * [i*NUM_TETS_PER_HEX, (i+1)*NUM_TETS_PER_HEX) are remapped to
     * [j*NUM_TETS_PER_HEX, (j+1)*NUM_TETS_PER_HEX), where j = cellIndices[i],
     * the real cell index.
     */
    virtual void remapTetIndices(axom::ArrayView<const axom::IndexType> cellIndices,
                                 axom::ArrayView<axom::IndexType> tetIndices) = 0;

    //!@brief Add volumes of tets inside the geometry to the volume data.
    virtual void addVolumesOfInteriorTets(axom::ArrayView<const axom::IndexType> cellsOnBdry,
                                          axom::ArrayView<const LabelType> tetLabels,
                                          axom::ArrayView<double> ovlap) = 0;

    //!@brief Compute clip volumes for every cell.
    virtual void computeClipVolumes3D(axom::ArrayView<double> ovlap, conduit::Node& statistics) = 0;

    //!@brief Compute clip volumes for cell in an index list.
    virtual void computeClipVolumes3D(const axom::ArrayView<axom::IndexType>& cellIndices,
                                      axom::ArrayView<double> ovlap,
                                      conduit::Node& statistics) = 0;

    /*!
     * @brief Compute clip volumes for cell tets in an index list.
     *
     * The tets are the results from decomposing each cell hex into
     * NUM_TETS_PER_HEX tets and stored consecutively.
     */
    virtual void computeClipVolumes3DTets(const axom::ArrayView<axom::IndexType>& tetIndices,
                                          axom::ArrayView<double> ovlap,
                                          conduit::Node& statistics) = 0;

    //!@brief Count the number of labels of each type.
    virtual void getLabelCounts(axom::ArrayView<const LabelType> labels,
                                std::int64_t& inCount,
                                std::int64_t& onCount,
                                std::int64_t& outCount) = 0;

    ShapeMesh& getShapeMesh() { return m_myClipper.m_shapeMesh; }

    MeshClipperStrategy& getStrategy() { return *m_myClipper.m_strategy; }

    MeshClipper& myClipper() { return m_myClipper; }

  private:
    //!@brief The MeshClipper that owns this Impl.
    MeshClipper& m_myClipper;
  };

private:
  friend Impl;

  quest::experimental::ShapeMesh& m_shapeMesh;

  //! @brief Shape-specific operations in clipping.
  std::shared_ptr<quest::experimental::MeshClipperStrategy> m_strategy;

  //! @brief Object where execution space code is instantiated.
  std::unique_ptr<Impl> m_impl;

  /* NOTE: MeshClipperStrategy is for shape-specific functions,
   * implemented externally.  Impl implements internal algorithms
   * for multiple execution spaces.
   */

  //! @brief Statistics
  conduit::Node m_counterStats;

  bool m_verbose;

  int m_screenLevel;

#if defined(__CUDACC__)
public:
#endif
  //!@brief Allocate a delegate for m_shapeMesh's runtime policy.
  std::unique_ptr<Impl> newImpl();

  //@{
  //!@name Convenience methods
  //!@brief Count the number of labels of each type.
  void getLabelCounts(const axom::Array<LabelType>& labels,
                      std::int64_t& inCount,
                      std::int64_t& onCount,
                      std::int64_t& outCount)
  {
    m_impl->getLabelCounts(labels, inCount, onCount, outCount);
  }
  //@}
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_MESHCLIPPER_HPP
