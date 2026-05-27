// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Shaper.hpp
 *
 * \brief Helper class for shaping queries
 */

#ifndef AXOM_QUEST_SHAPER__HPP_
#define AXOM_QUEST_SHAPER__HPP_

#include "axom/config.hpp"
#ifndef AXOM_USE_KLEE
  #error Shaping functionality requires Axom to be configured with the Klee component
#endif

#if !defined(AXOM_USE_MFEM) && !defined(AXOM_USE_CONDUIT)
  #error Shaping functionality requires Axom to be configured with Conduit or MFEM
#endif

#include "axom/sidre.hpp"
#include "axom/klee.hpp"
#include "axom/mint.hpp"
#include "axom/quest/DiscreteShape.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"
#include "axom/core/execution/runtime_policy.hpp"

#if defined(AXOM_USE_MFEM)
  #include "mfem.hpp"
#endif
#if defined(AXOM_USE_CONDUIT)
  #include "conduit_node.hpp"
#endif

#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"

namespace axom
{
namespace quest
{
/**
 * Abstract base class for shaping material volume fractions
 *
 * Shaper requires Axom to be configured with Conduit or MFEM or both.
 */
class Shaper
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

#if defined(AXOM_USE_MFEM)
  /// @brief Construct Shaper to operate on an MFEM mesh.
  Shaper(RuntimePolicy execPolicy,
         int allocatorId,
         const klee::ShapeSet& shapeSet,
         sidre::MFEMSidreDataCollection* dc);
#endif

#if defined(AXOM_USE_CONDUIT)
  /*!
   * @brief Construct Shaper to operate on a blueprint-formatted mesh
   * stored in a sidre Group.
   */
  Shaper(RuntimePolicy execPolicy,
         int allocatorId,
         const klee::ShapeSet& shapeSet,
         sidre::Group* bpMesh,
         const std::string& topo = "");

  /*!
   * @brief Construct Shaper to operate on a blueprint-formatted mesh
   * stored in a conduit Node.
   * 
   * Because \c conduit::Node doesn't support application-specified
   * allocator id for (only) arrays, the incoming \c bpNode must have
   * all arrays pre-allocated in a space accessible by the runtime
   * policy.  Any needed-but-missing space would lead to an exception.
  */
  Shaper(RuntimePolicy execPolicy,
         int allocatorId,
         const klee::ShapeSet& shapeSet,
         conduit::Node& bpNode,
         const std::string& topo = "");
#endif

  virtual ~Shaper();

public:
  // Some default values.
  static constexpr int DEFAULT_SAMPLES_PER_KNOT_SPAN {25};
  static constexpr double MINIMUM_PERCENT_ERROR {0.};
  static constexpr double MAXIMUM_PERCENT_ERROR {100.};
  static constexpr double DEFAULT_VERTEX_WELD_THRESHOLD {1e-9};

  /// Refinement type.
  using RefinementType = DiscreteShape::RefinementType;

  //! @brief Verify the input mesh is okay for this backend to work with.
  bool verifyInputMesh(std::string& whyBad) const;

  ///@{
  //!  @name Functions to get and set shaping parameters

  void setSamplesPerKnotSpan(int nSamples);
  void setVertexWeldThreshold(double threshold);
  void setVerbosity(bool isVerbose) { m_verboseOutput = isVerbose; }
  void setPercentError(double percent);
  void setRefinementType(RefinementType t);

  ///@}

  /*!
   * @brief Set path of shape input file.
   *
   * The path is used to resolve relative paths that may have been
   * specified in the file.
  */
  void setFilePath(const std::string& filePath);

  mint::Mesh* getSurfaceMesh() const { return m_surfaceMesh.get(); }

  bool isVerbose() const { return m_verboseOutput; }

#ifdef AXOM_USE_MFEM
  sidre::MFEMSidreDataCollection* getDC()
  {
    return m_mfem_state != nullptr ? m_mfem_state->m_dc : nullptr;
  }
  const sidre::MFEMSidreDataCollection* getDC() const
  {
    return m_mfem_state != nullptr ? m_mfem_state->m_dc : nullptr;
  }
#endif

#if defined(AXOM_USE_CONDUIT)
  conduit::Node* getBlueprintMeshNode()
  {
    return m_bp_state != nullptr ? &m_bp_state->m_internal_node : nullptr;
  }
  const conduit::Node* getBlueprintMeshNode() const
  {
    return m_bp_state != nullptr ? &m_bp_state->m_internal_node : nullptr;
  }
#endif

  /*!
   * \brief Predicate to determine if the specified format is valid
   *
   * \param format A string listing the format to check
   */
  virtual bool isValidFormat(const std::string& format) const;

  /// \brief Returns the format type of the supplied \a shape
  std::string shapeFormat(const klee::Shape& shape) const
  {
    return shape.getGeometry().getFormat();
  }

  /// \brief Returns the execution policy used by the Shaper
  RuntimePolicy getExecutionPolicy() const { return m_execPolicy; }

public:
  ///@{
  ///  @name Functions related to the stages for a given shape

  /// Loads the shape from file into m_surfaceMesh
  virtual void loadShape(const klee::Shape& shape);

  virtual void prepareShapeQuery(klee::Dimensions shapeDimension, const klee::Shape& shape) = 0;

  virtual void runShapeQuery(const klee::Shape& shape) = 0;

  virtual void applyReplacementRules(const klee::Shape& shape) = 0;

  virtual void finalizeShapeQuery() = 0;

  ///@}

public:
  ///@{
  ///  @name Functions to generate/adjust volume fractions after all shapes have been applied

  virtual void adjustVolumeFractions() = 0;

  ///@}

  /*!
   * \brief Save the shaping results to disk.
   *
   * \param extra Save extra data when available.
   */
  virtual void saveResults(bool extra);

  /*!
   * \brief Helper to apply a parallel sum reduction to a quantity
   *
   * \note This is the identity function when running without MPI 
   */
  double allReduceSum(double val) const;

  /*!
   * \brief Helper to apply a parallel min reduction to a quantity
   *
   * \note This is the identity function when running without MPI 
   */
  double allReduceMin(double val) const;

  /*!
   * \brief Helper to apply a parallel max reduction to a quantity
   *
   * \note This is the identity function when running without MPI 
   */
  double allReduceMax(double val) const;

protected:
  /*!
   * \brief Loads the shape into m_surfaceMesh.
   * \param shape The shape.
   * \param percentError A percent error to use when refining the shape. If it
   *                     positive then Axom will try to refine dynamically
   *                     according to this error. Otherwise, it will use the
   *                     segmentsPerKnotSpan value.
   * \param[out] revolvedvolume A revolved volume for the shape, if the shape
   *             is from a C2C contour.
   */
  void loadShapeInternal(const klee::Shape& shape, double percentError, double& revolvedVolume);

  /*!
   * \brief Computes transforms for the shape and applies them to the surface mesh.
   * \param shape The shape.
   */
  void applyTransforms(const klee::Shape& shape);

  /*!
   * \brief Computes transforms for the shape and applies them to the surface mesh.
   * \param shape The shape.
   * \param transform A 4x4 matrix containing the transformation to apply.
   */
  void applyTransforms(const numerics::Matrix<double>& transform);

  /*!
   * \brief Get a matrix that contains the shape's concatenated transforms.
   *
   * \param shape The shape whose transforms are being concatenated.
   *
   * \return A 4x4 matrix that represents the transforms.
   */
  numerics::Matrix<double> getTransforms(const klee::Shape& shape) const;

  /*!
   * \brief Helper function to get the rank associated with the current process
   *
   * \note This function can be called even in non-mpi configurations
   */
  int getRank() const;

  /*!
   * \brief Backend-specific input-mesh validation hook.
   *
   * Derived shapers may support different Blueprint mesh representations, so
   * the validation policy lives with the concrete backend.
   */
  virtual bool verifyInputMeshImpl(std::string& whyBad) const = 0;

#if defined(AXOM_USE_CONDUIT)
  /*!
   * \brief Selects the Blueprint topology name to use and verifies it exists.
   */
  std::string resolveBlueprintTopologyName(const sidre::Group* bpMesh, const std::string& topo) const;

  /*!
   * \brief Selects the Blueprint topology name to use and verifies it exists.
   */
  std::string resolveBlueprintTopologyName(const conduit::Node& bpMesh, const std::string& topo) const;

  /*!
   * \brief Rebuilds the internal Conduit view and cached cell count from the
   *        current Sidre-owned Blueprint mesh.
   */
  void refreshBlueprintMeshState();

  /*!
   * \brief Returns the active Blueprint topology node.
   */
  const conduit::Node& getBlueprintTopologyNode() const;

  /*!
   * \brief Returns the active Blueprint coordset node.
   */
  const conduit::Node& getBlueprintCoordsetNode() const;

  /*!
   * \brief Returns the active Blueprint cell shape name.
   */
  std::string getBlueprintCellShape() const;

  /*!
   * \brief Returns the active Blueprint mesh dimension for supported quad/hex
   *        meshes.
   */
  int getBlueprintMeshDimension() const;

  /*!
   * \brief Helper for Blueprint meshes supported directly by sampling or by
   *        lazy conversion in the intersection backend.
   *
   * This helper verifies the internal Blueprint mesh uses a structured or
   * unstructured quad/hex topology over an explicit coordset.
   */
  bool verifyBlueprintMeshIsStructuredOrUnstructuredQuadHex(std::string& whyBad) const;

  /*!
   * \brief Converts a structured explicit Blueprint quad/hex mesh to an
   *        unstructured working representation if needed.
   */
  void ensureBlueprintMeshIsUnstructured();
#endif

  /*!
   * \brief Get the protocol to use for Blueprint output.
   *
   * \return "hdf5" when possible, otherwise "yaml".
   */
  std::string outputProtocol() const;

#if defined(AXOM_USE_MFEM)
  //! \brief MFEM meshes currently have no additional validation here.
  bool verifyMFEMInputMesh(std::string& whyBad) const;
#endif

#if defined(AXOM_USE_MFEM)
  virtual std::unique_ptr<shaping::MFEMState> createMFEMState()
  {
    return std::make_unique<shaping::MFEMState>();
  }
#endif
#if defined(AXOM_USE_CONDUIT)
  virtual std::unique_ptr<shaping::BlueprintState> createBlueprintState()
  {
    return std::make_unique<shaping::BlueprintState>();
  }
#endif

protected:
  RuntimePolicy m_execPolicy;
  int m_allocatorId;

  // For any mesh represented in Conduit or sidre
  sidre::DataStore m_dataStore;

  const klee::ShapeSet& m_shapeSet;

  //! \brief Prefix path for shape file names with relative path.
  std::string m_prefixPath;

#if defined(AXOM_USE_MFEM)
  std::unique_ptr<shaping::MFEMState> m_mfem_state;
#endif
#if defined(AXOM_USE_CONDUIT)
  std::unique_ptr<shaping::BlueprintState> m_bp_state;
#endif

  //! @brief Number of cells in the computational mesh.
  axom::IndexType m_cellCount;

  std::shared_ptr<mint::Mesh> m_surfaceMesh;

  int m_samplesPerKnotSpan {DEFAULT_SAMPLES_PER_KNOT_SPAN};
  double m_percentError {MINIMUM_PERCENT_ERROR};
  RefinementType m_refinementType {DiscreteShape::RefinementUniformSegments};
  double m_vertexWeldThreshold {DEFAULT_VERTEX_WELD_THRESHOLD};
  bool m_verboseOutput {false};

  MPI_Comm m_comm {MPI_COMM_SELF};
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SHAPER__HPP_
