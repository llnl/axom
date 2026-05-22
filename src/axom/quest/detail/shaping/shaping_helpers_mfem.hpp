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

struct MFEMState
{
  virtual ~MFEMState() = default;

  int meshDimension() const { return m_dc->GetMesh()->Dimension(); }

  sidre::MFEMSidreDataCollection* m_dc {nullptr};
};

struct SamplingMFEMState : public MFEMState
{
  ~SamplingMFEMState() override
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

  QFunctionCollection m_inoutShapeQFuncs;
  QFunctionCollection m_inoutMaterialQFuncs;
  DenseTensorCollection m_inoutTensors;
  MFEMArrayCollection m_inoutArrays;
};

void printRegisteredFieldNames(const SamplingMFEMState& mfemState,
                               const std::set<std::string>& knownMaterials,
                               VolFracSampling vfSampling,
                               const std::string& initialMessage);

mfem::GridFunction* getOrAllocateL2GridFunction(mfem::DataCollection* dc,
                                                const std::string& gf_name,
                                                int order,
                                                int dim,
                                                const int basis);

void replaceMaterial(mfem::QuadratureFunction* shapeQFunc,
                     mfem::QuadratureFunction* materialQFunc,
                     bool shouldReplace);

void copyShapeIntoMaterial(const mfem::QuadratureFunction* shapeQFunc,
                           mfem::QuadratureFunction* materialQFunc,
                           bool reuseExisting = true);

mfem::QuadratureFunction* cloneInOutFunction(const mfem::QuadratureFunction* qfunc);

void generatePositionsQFunction(mfem::Mesh* mesh,
                                QFunctionCollection& inoutQFuncs,
                                axom::ArrayView<int> sampleResolution,
                                axom::numerics::QuadratureType quadratureType);

void generateSamplingPositions(SamplingMFEMState& mfemState,
                               axom::ArrayView<int> sampleResolution,
                               axom::numerics::QuadratureType quadratureType);

void computeVolumeFractionsForMaterial(SamplingMFEMState& mfemState,
                                       const std::string& matField,
                                       int volfracOrder,
                                       axom::ArrayView<int> sampleResolution,
                                       axom::numerics::QuadratureType quadratureType);

void computeVolumeFractionsIdentity(mfem::DataCollection* dc,
                                    mfem::QuadratureFunction* inout,
                                    const std::string& name);

bool usesAnisotropicCustomTensorQuadrature(const mfem::Mesh& mesh,
                                           axom::ArrayView<int> sampleResolution,
                                           axom::numerics::QuadratureType quadratureType);

template <int FromDim, int ToDim, typename InsideFunc>
void sampleInOutField(const std::string shapeName,
                      shaping::SamplingMFEMState& mfemState,
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

  const std::string inoutName = axom::fmt::format("inout_{}", shapeName);
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

template <int FromDim, int ToDim, typename InsideFunc>
void computeVolumeFractionsBaseline(const std::string& shapeName,
                                    shaping::SamplingMFEMState& mfemState,
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

  const auto volFracName = axom::fmt::format("vol_frac_{}", shapeName);
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

void FCT_correct(const double* M,
                 const int s,
                 const double* m,
                 const double y_min,
                 const double y_max,
                 double* xy,
                 double* fct_mat);

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#endif  // defined(AXOM_USE_MFEM)

#endif  // AXOM_QUEST_SHAPING_HELPERS_MFEM__HPP_
