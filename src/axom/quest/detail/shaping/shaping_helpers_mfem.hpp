// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SHAPING_HELPERS_MFEM__HPP_
#define AXOM_QUEST_SHAPING_HELPERS_MFEM__HPP_

#include "shaping_helpers.hpp"

#if defined(AXOM_USE_MFEM)

  #include "axom/fmt.hpp"

  #include "mfem.hpp"
  #include "mfem/linalg/dtensor.hpp"

  #include <set>
  #include <string>
  #include <vector>

namespace axom
{
namespace quest
{
namespace shaping
{

int to_mfem_quadrature_type(axom::numerics::QuadratureType quadratureType);

using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;
using MFEMArrayCollection = mfem::NamedFieldsMap<mfem::Array<int>>;

/// MFEM state shared by Quest shapers and MFEM-backed sampling helpers.
struct MFEMState
{
  ~MFEMState()
  {
    m_inoutShapeQFuncs.DeleteData(true);
    m_inoutShapeQFuncs.clear();

    m_inoutMaterialQFuncs.DeleteData(true);
    m_inoutMaterialQFuncs.clear();

    m_inoutTensors.DeleteData(true);
    m_inoutTensors.clear();

    m_inoutArrays.DeleteData(true);
    m_inoutArrays.clear();
  }

  int meshDimension() const { return m_dc->GetMesh()->Dimension(); }

  sidre::MFEMSidreDataCollection* m_dc {nullptr};

  mfem::QuadratureFunction* getShapeFunction(const std::string& name)
  {
    return m_inoutShapeQFuncs.Get(name);
  }

  const mfem::QuadratureFunction* getShapeFunction(const std::string& name) const
  {
    return m_inoutShapeQFuncs.Get(name);
  }

  void deleteShapeFunction(const std::string& AXOM_UNUSED_PARAM(name))
  {
    // TODO: remove the function from m_inoutShapeQFuncs if it exists.
  }

  mfem::QuadratureFunction* getMaterialFunction(const std::string& name)
  {
    return m_inoutMaterialQFuncs.Get(name);
  }

  const mfem::QuadratureFunction* getMaterialFunction(const std::string& name) const
  {
    return m_inoutMaterialQFuncs.Get(name);
  }

  mfem::QuadratureFunction* createMaterialFunction(const std::string& name)
  {
    auto* positions = m_inoutShapeQFuncs.Get("positions");
    SLIC_ERROR_IF(positions == nullptr,
                  std::string("Cannot create material function '") + name + "' without positions.");

    auto* qfunc = new mfem::QuadratureFunction(positions->GetSpace(), 1);
    qfunc->HostWrite();
    *qfunc = 0.;
    m_inoutMaterialQFuncs.Register(name, qfunc, true);
    return qfunc;
  }

  QFunctionCollection& shapeQFuncs() { return m_inoutShapeQFuncs; }
  const QFunctionCollection& shapeQFuncs() const { return m_inoutShapeQFuncs; }

  QFunctionCollection& materialQFuncs() { return m_inoutMaterialQFuncs; }
  const QFunctionCollection& materialQFuncs() const { return m_inoutMaterialQFuncs; }

  DenseTensorCollection& tensors() { return m_inoutTensors; }
  const DenseTensorCollection& tensors() const { return m_inoutTensors; }

  MFEMArrayCollection& arrays() { return m_inoutArrays; }
  const MFEMArrayCollection& arrays() const { return m_inoutArrays; }

  QFunctionCollection m_inoutShapeQFuncs;
  QFunctionCollection m_inoutMaterialQFuncs;
  DenseTensorCollection m_inoutTensors;
  MFEMArrayCollection m_inoutArrays;
};

/*!
 * \brief Print the registered field names in the \a mfemState.
 *
 * \param mfemState The MFEM state.
 * \param knownMaterials A set of known material names.
 * \param vfSampling The type of volume fraction sampling being performed.
 * \param initialMessage A string to prepend to the printed message.
 */
void printRegisteredFieldNames(const MFEMState& mfemState,
                               const std::set<std::string>& knownMaterials,
                               VolFracSampling vfSampling,
                               const std::string& initialMessage);

/*!
 * \brief Utility function to either return a grid function from the DataCollection \a dc, 
 * or to allocate the grud function through the dc, ensuring the memory doesn't leak
 * 
 * \return A pointer to the (allocated) grid function. nullptr if it cannot be allocated
 */
mfem::GridFunction* getOrAllocateL2GridFunction(mfem::DataCollection* dc,
                                                const std::string& gf_name,
                                                int order,
                                                int dim,
                                                const int basis);

/*!
 * Utility function to zero out inout quadrature points for a material replaced by a shape
 *
 * Each location in space can only be covered by one material.
 * When \a shouldReplace is true, we clear all values in \a materialQFunc 
 * that are set in \a shapeQFunc. When it is false, we do the opposite.
 *
 * \param shapeQFunc The inout quadrature function for the shape samples
 * \param materialQFunc The inout quadrature function for the material samples
 * \param shouldReplace Flag for whether the shape replaces the material 
 *   or whether the material remains and we should zero out the shape sample (when false)
 */
void replaceMaterial(mfem::QuadratureFunction* shapeQFunc,
                     mfem::QuadratureFunction* materialQFunc,
                     bool shouldReplace);

/*!
 * \brief Utility function to copy inout quadrature point values from \a shapeQFunc to \a materialQFunc
 *
 * \param shapeQFunc The inout samples for the current shape
 * \param materialQFunc The inout samples for the material we're writing into
 * \param reuseExisting When a value is not set in \a shapeQFunc, should we retain existing values 
 * from \a materialQFunc or overwrite them based on \a shapeQFunc. The default is to retain values
 */
void copyShapeIntoMaterial(const mfem::QuadratureFunction* shapeQFunc,
                           mfem::QuadratureFunction* materialQFunc,
                           bool reuseExisting = true);

/*!
 * \brief Create a copy of the supplied function.
 *
 * \param qfunc A pointer to the function to clone.
 *
 * \return A pointer to a new copy of the supplied function.
 */
mfem::QuadratureFunction* cloneInOutFunction(const mfem::QuadratureFunction* qfunc);

/*!
 * \brief Generates sampling positions within each zone based on element quadrature.
 *
 * \param mesh The MFEM mesh.
 * \param inoutQFuncs A function collection in which to place the position function.
 * \param sampleResolution The number of samples in each dimension.
 * \param quadratureType The quadrature type that determines the sample locations.
 */
void generatePositionsQFunction(mfem::Mesh* mesh,
                                QFunctionCollection& inoutQFuncs,
                                axom::ArrayView<int> sampleResolution,
                                axom::numerics::QuadratureType quadratureType);

/*!
 * \brief Generates sampling positions within each zone based on element quadrature.
 *
 * \param mfemState The MFEM state.
 * \param sampleResolution The number of samples in each dimension.
 * \param quadratureType The quadrature type that determines the sample locations.
 *
 * \note The sample points are stored as a function corresponding to the mesh positions
 */
void generateSamplingPositions(MFEMState& mfemState,
                               axom::ArrayView<int> sampleResolution,
                               axom::numerics::QuadratureType quadratureType);

/*!
 * \brief Import initial volume fractions from the map into the quadrature
 *        "mat_inout_" fields in \a mfemState.
 *
 * \param mfemState The MFEM state.
 * \param initialVolumeFractions A map of initial volume fraction fields used to
 *                               initialize mat_inout fields over the quadrature
 *                               points.
 * \param anisotropic Whether the quadrature points are anisotropic.
 */
void importInitialVolumeFractions(MFEMState& mfemState,
                                  const std::map<std::string, mfem::GridFunction*>& initialVolumeFractions,
                                  bool anisotropic);

/*!
 * \brief Create volume fractions for a material using the existing material field
 *        (mat_inout_{matField}) to make the new field (vol_fract_{matField}).
 *
 * \param mfemState The MFEM state that contains the mesh and functions.
 * \param matField The name of the material field.
 * \param volfracOrder The order of the volume fraction function to create.
 * \param sampleResolution The number of samples in each mesh dimension.
 * \param quadratureType The quadrature type that determines the sample point locations.
 */
void computeVolumeFractionsForMaterial(MFEMState& mfemState,
                                       const std::string& matField,
                                       int volfracOrder,
                                       axom::ArrayView<int> sampleResolution,
                                       axom::numerics::QuadratureType quadratureType,
                                       axom::runtime_policy::Policy execPolicy);

void computeVolumeFractionsIdentity(mfem::DataCollection* dc,
                                    mfem::QuadratureFunction* inout,
                                    const std::string& name);

/*!
 * \brief Examines the sample resolution and quadrature rule to decide whether the
 *        requested quadrature is a custom-tensor / anisotropic. Some algorithms for
 *        MFEM must take special paths in this case.
 *
 * \param mesh The MFEM mesh.
 * \param sampleResolution The number of samples in each mesh dimension.
 * \param quadratureType The quadrature type that determines the sample point locations.
 *
 * \return True if the quadrature is custom / anisotropic; false otherwise.
 */
bool usesAnisotropicCustomTensorQuadrature(const mfem::Mesh& mesh,
                                           axom::ArrayView<int> sampleResolution,
                                           axom::numerics::QuadratureType quadratureType);

/*!
  * \brief Samples the inout field over the indexed geometry, possibly using a
  * callback function to project the input points (from the computational mesh)
  * to query points on the spatial index
  *
  * \tparam FromDim The dimension of points from the input mesh
  * \tparam ToDim The dimension of points on the indexed shape
  * \tparam InsideFunc A function that takes a point and returns a bool indicating whether the
  *                    point is inside or outside of relevant shapes.
  *
  * \param [in] shapeName The name of the shape used in making data array names.
  * \param [in] mfemState The data collection containing the mesh, associated query points
  *                       and a collection of quadrature functions for the shape and material
  *                       inout samples.
  * \param [in] checkInside The function that determines whether a point is inside.
  * \param [in] projector A callback function to apply to points from the input mesh
  * before querying them on the spatial index
  *
  * \note A projector callback must be supplied when \a FromDim is not equal
  *       to \a ToDim.
  */
template <int FromDim, int ToDim, typename InsideFunc>
void sampleInOutField(const std::string shapeName,
                      shaping::MFEMState& mfemState,
                      InsideFunc&& checkInside,
                      PointProjector<FromDim, ToDim> projector = {})
{
  using FromPoint = primal::Point<double, FromDim>;
  using ToPoint = primal::Point<double, ToDim>;
  AXOM_ANNOTATE_SCOPE("sampleInOutField");

  SLIC_ERROR_IF(FromDim != ToDim && !projector,
                "A projector callback function is required when FromDim != ToDim");

  auto* mesh = mfemState.m_dc->GetMesh();
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  auto& inoutQFuncs = mfemState.m_inoutShapeQFuncs;
  SLIC_ASSERT(inoutQFuncs.Has("positions"));

  mfem::QuadratureFunction* pos_coef = inoutQFuncs.Get("positions");
  auto* sp = pos_coef->GetSpace();
  const int nq = sp->GetIntRule(0).GetNPoints();
  const int numQueryPoints = sp->GetSize();
  SLIC_ASSERT(numQueryPoints == NE * nq);

  const auto pos = mfem::Reshape(pos_coef->HostRead(), dim, nq, NE);

  const std::string inoutName = shaping::shapeInOutFieldName(shapeName);
  auto* inout = new mfem::QuadratureFunction(sp, 1);
  inoutQFuncs.Register(inoutName, inout, true);
  auto inout_vals = mfem::Reshape(inout->HostWrite(), nq, NE);

  axom::utilities::Timer timer(true);
  if(projector)
  {
    for(int i = 0; i < NE; ++i)
    {
      for(int p = 0; p < nq; ++p)
      {
        const ToPoint pt = projector(FromPoint(&pos(0, p, i), dim));
        inout_vals(p, i) = checkInside(pt) ? 1. : 0.;
      }
    }
  }
  else
  {
    for(int i = 0; i < NE; ++i)
    {
      for(int p = 0; p < nq; ++p)
      {
        const ToPoint pt(&pos(0, p, i), dim);
        inout_vals(p, i) = checkInside(pt) ? 1. : 0.;
      }
    }
  }
  timer.stop();

  SLIC_INFO_ROOT(axom::fmt::format(
    axom::utilities::locale(),
    "\t Sampling inout field '{}' took {:.3Lf} seconds (@ {:L} queries per second)",
    inoutName,
    timer.elapsed(),
    static_cast<int>(numQueryPoints / timer.elapsed())));
}

/*!
 * \brief Called when sampling shapes at dofs.
 *
 * \tparam FromDim The dimension of points from the input mesh
 * \tparam ToDim The dimension of points on the indexed shape
 * \tparam InsideFunc A function that takes a point and returns a bool indicating whether the
 *                    point is inside or outside of relevant shapes.
 *
 * \param [in] shapeName The name of the shape used in making data array names.
 * \param [in] mfemState The data collection containing the mesh, associated query points
 *                       and a collection of quadrature functions for the shape and material
 *                       inout samples.
 * \param [in] outputOrder The order of the volume fraction function.
 * \param [in] checkInside The function that determines whether a point is inside.
 * \param [in] projector A callback function to apply to points from the input mesh
 * before querying them on the spatial index
 *
 * \note A projector callback must be supplied when \a FromDim is not equal
 *       to \a ToDim.
 */
template <int FromDim, int ToDim, typename InsideFunc>
void computeVolumeFractionsBaseline(const std::string& shapeName,
                                    shaping::MFEMState& mfemState,
                                    int outputOrder,
                                    InsideFunc&& checkInside,
                                    PointProjector<FromDim, ToDim> projector = {})
{
  using FromPoint = primal::Point<double, FromDim>;
  using ToPoint = primal::Point<double, ToDim>;
  AXOM_ANNOTATE_SCOPE("computeVolumeFractionsBaseline");

  mfem::DataCollection* dc = mfemState.m_dc;
  mfem::Mesh* mesh = dc->GetMesh();
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return;
  }

  const auto volFracName = shaping::volumeFractionFieldName(shapeName);
  mfem::GridFunction* volFrac =
    shaping::getOrAllocateL2GridFunction(dc, volFracName, outputOrder, dim, mfem::BasisType::Positive);
  const mfem::FiniteElementSpace* fes = volFrac->FESpace();

  auto* fe = fes->GetFE(0);
  auto& ir = fe->GetNodes();

  const int nq = ir.GetNPoints();
  const auto* geomFactors = mesh->GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);

  mfem::DenseTensor pos_coef(dim, nq, NE);
  for(int i = 0; i < NE; ++i)
  {
    for(int j = 0; j < dim; ++j)
    {
      for(int k = 0; k < nq; ++k)
      {
        pos_coef(j, k, i) = geomFactors->X((i * nq * dim) + (j * nq) + k);
      }
    }
  }

  mfem::Vector res(nq);
  mfem::Array<int> dofs;
  if(projector)
  {
    for(int i = 0; i < NE; ++i)
    {
      const mfem::DenseMatrix& m = pos_coef(i);
      for(int p = 0; p < nq; ++p)
      {
        const ToPoint pt = projector(FromPoint(m.GetColumn(p), dim));
        res(p) = checkInside(pt) ? 1. : 0.;
      }

      fes->GetElementDofs(i, dofs);
      volFrac->SetSubVector(dofs, res);
    }
  }
  else
  {
    for(int i = 0; i < NE; ++i)
    {
      const mfem::DenseMatrix& m = pos_coef(i);
      for(int p = 0; p < nq; ++p)
      {
        const ToPoint pt(m.GetColumn(p), dim);
        res(p) = checkInside(pt) ? 1. : 0.;
      }

      fes->GetElementDofs(i, dofs);
      volFrac->SetSubVector(dofs, res);
    }
  }
}

/** 
 * Implements flux-corrected transport (FCT) to correct the solution obtained
 * when converting from inout samples (ones and zeros) to a grid function 
 * on the degrees of freedom such that the volume fractions are doubles
 * between 0 and 1 ( \a y_min and \a y_max )
 */
void FCT_correct(const double* M,
                 const int s,
                 const double* m,
                 const double y_min,  // 0
                 const double y_max,  // 1
                 double* xy,
                 double* fct_mat);  // scratch buffer

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#endif  // defined(AXOM_USE_MFEM)

#endif  // AXOM_QUEST_SHAPING_HELPERS_MFEM__HPP_
