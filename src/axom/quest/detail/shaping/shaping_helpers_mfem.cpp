// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "shaping_helpers_mfem.hpp"

#if defined(AXOM_USE_MFEM)

  #include <memory>
  #include <vector>

namespace axom
{
namespace quest
{
namespace shaping
{

namespace
{

class OwnedQuadratureSpace : public mfem::QuadratureSpace
{
public:
  OwnedQuadratureSpace(mfem::Mesh& mesh, std::unique_ptr<mfem::IntegrationRule> ir)
    : mfem::QuadratureSpace(mesh, *ir)
    , m_ir(std::move(ir))
  { }

private:
  std::unique_ptr<mfem::IntegrationRule> m_ir;
};

}  // namespace

bool usesAnisotropicCustomTensorQuadrature(const mfem::Mesh& mesh,
                                           axom::ArrayView<int> sampleResolution,
                                           axom::numerics::QuadratureType quadratureType)
{
  if(quadratureType == axom::numerics::QuadratureType::Invalid)
  {
    return false;
  }

  const auto dim = mesh.Dimension();
  SLIC_ERROR_IF(sampleResolution.size() != static_cast<axom::IndexType>(dim),
                "Sample resolution dimension does not match mesh dimension");

  if(mesh.GetNE() > 0)
  {
    switch(mesh.GetTypicalElementGeometry())
    {
    case mfem::Geometry::SQUARE:
      return sampleResolution[0] != sampleResolution[1];
    case mfem::Geometry::CUBE:
      return sampleResolution[0] != sampleResolution[1] || sampleResolution[0] != sampleResolution[2];
    default:
      return false;
    }
  }
  return mesh.Dimension();
}

int to_mfem_quadrature_type(axom::numerics::QuadratureType quadratureType)
{
  switch(quadratureType)
  {
  case axom::numerics::QuadratureType::Invalid:
    return mfem::Quadrature1D::Invalid;
  case axom::numerics::QuadratureType::GaussLegendre:
    return mfem::Quadrature1D::GaussLegendre;
  case axom::numerics::QuadratureType::GaussLobatto:
    return mfem::Quadrature1D::GaussLobatto;
  case axom::numerics::QuadratureType::OpenUniform:
    return mfem::Quadrature1D::OpenUniform;
  case axom::numerics::QuadratureType::ClosedUniform:
    return mfem::Quadrature1D::ClosedUniform;
  case axom::numerics::QuadratureType::OpenHalfUniform:
    return mfem::Quadrature1D::OpenHalfUniform;
  case axom::numerics::QuadratureType::ClosedGL:
    return mfem::Quadrature1D::ClosedGL;
  }

  SLIC_ERROR(axom::fmt::format("Invalid quadrature type {}.", static_cast<int>(quadratureType)));
  return mfem::Quadrature1D::Invalid;
}

mfem::GridFunction* getOrAllocateL2GridFunction(mfem::DataCollection* dc,
                                                const std::string& gf_name,
                                                int order,
                                                int dim,
                                                const int basis)
{
  if(dc == nullptr)
  {
    SLIC_WARNING("Cannot allocate grid function into null data collection");
    return nullptr;
  }

  mfem::GridFunction* gf = nullptr;

  if(dc->HasField(gf_name))
  {
    gf = dc->GetField(gf_name);
  }
  else
  {
    auto* fec = new mfem::L2_FECollection(order, dim, basis);
    auto* mesh = dc->GetMesh();
    mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, fec);

    auto* sidreDC = dynamic_cast<sidre::MFEMSidreDataCollection*>(dc);
    if(sidreDC)
    {
      const int sz = fes->GetVSize();
      auto* vw = sidreDC->AllocNamedBuffer(gf_name, sz);
      gf = new mfem::GridFunction();
      gf->MakeRef(fes, vw->getData());
    }
    else
    {
      gf = new mfem::GridFunction(fes);
    }

    gf->MakeOwner(fec);
    gf->HostReadWrite();
    *gf = 0.;

    dc->RegisterField(gf_name, gf);
  }

  return gf;
}

void replaceMaterial(mfem::QuadratureFunction* shapeQFunc,
                     mfem::QuadratureFunction* materialQFunc,
                     bool shapeReplacesMaterial)
{
  SLIC_ASSERT(shapeQFunc != nullptr);
  SLIC_ASSERT(materialQFunc != nullptr);
  SLIC_ASSERT(materialQFunc->Size() == shapeQFunc->Size());

  const int SZ = materialQFunc->Size();
  double* mData = materialQFunc->HostReadWrite();
  double* sData = shapeQFunc->HostReadWrite();

  if(shapeReplacesMaterial)
  {
    for(int j = 0; j < SZ; ++j)
    {
      mData[j] = sData[j] > 0 ? 0 : mData[j];
    }
  }
  else
  {
    for(int j = 0; j < SZ; ++j)
    {
      sData[j] = mData[j] > 0 ? 0 : sData[j];
    }
  }
}

void copyShapeIntoMaterial(const mfem::QuadratureFunction* shapeQFunc,
                           mfem::QuadratureFunction* materialQFunc,
                           bool reuseExisting)
{
  SLIC_ASSERT(shapeQFunc != nullptr);
  SLIC_ASSERT(materialQFunc != nullptr);
  SLIC_ASSERT(materialQFunc->Size() == shapeQFunc->Size());

  const int SZ = materialQFunc->Size();
  double* mData = materialQFunc->HostReadWrite();
  const double* sData = shapeQFunc->HostRead();

  if(reuseExisting)
  {
    for(int j = 0; j < SZ; ++j)
    {
      mData[j] = sData[j] > 0 ? 1 : mData[j];
    }
  }
  else
  {
    for(int j = 0; j < SZ; ++j)
    {
      mData[j] = sData[j];
    }
  }
}

mfem::QuadratureFunction* cloneInOutFunction(const mfem::QuadratureFunction* qfunc)
{
  SLIC_ASSERT(qfunc != nullptr);
  return new mfem::QuadratureFunction(*qfunc);
}

void printRegisteredFieldNames(const SamplingMFEMState& mfemState,
                               const std::set<std::string>& knownMaterials,
                               VolFracSampling vfSampling,
                               const std::string& initialMessage)
{
  SLIC_ASSERT(mfemState.m_dc != nullptr);

  auto extractKeys = [](const auto& map) {
    std::vector<std::string> keys;
    for(const auto& kv : map)
    {
      keys.push_back(kv.first);
    }
    return keys;
  };

  axom::fmt::memory_buffer out;
  axom::fmt::format_to(std::back_inserter(out),
                       "List of registered fields in the SamplingShaper {}"
                       "\n\t* Data collection grid funcs: {}"
                       "\n\t* Data collection qfuncs: {}"
                       "\n\t* Known materials: {}",
                       initialMessage,
                       axom::fmt::join(extractKeys(mfemState.m_dc->GetFieldMap()), ", "),
                       axom::fmt::join(extractKeys(mfemState.m_dc->GetQFieldMap()), ", "),
                       axom::fmt::join(knownMaterials, ", "));

  if(vfSampling == VolFracSampling::SAMPLE_AT_QPTS)
  {
    axom::fmt::format_to(std::back_inserter(out),
                         "\n\t* Shape qfuncs: {}"
                         "\n\t* Mat qfuncs: {}",
                         axom::fmt::join(extractKeys(mfemState.m_inoutShapeQFuncs), ", "),
                         axom::fmt::join(extractKeys(mfemState.m_inoutMaterialQFuncs), ", "));
  }
  else if(vfSampling == VolFracSampling::SAMPLE_AT_DOFS)
  {
    axom::fmt::format_to(std::back_inserter(out),
                         "\n\t* Shaping tensors: {}",
                         axom::fmt::join(extractKeys(mfemState.m_inoutTensors), ", "));
  }

  SLIC_INFO_ROOT(axom::fmt::to_string(out));
}

namespace
{

mfem::QuadratureSpace* makeDefaultQuadratureSpace(mfem::Mesh* mesh, int sampleRes)
{
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return nullptr;
  }

  const int sampleOrder = 2 * sampleRes - 1;
  return new mfem::QuadratureSpace(mesh, sampleOrder);
}

mfem::QuadratureSpace* makeCustomQuadratureSpace(mfem::Mesh* mesh,
                                                 axom::ArrayView<int> sampleRes,
                                                 axom::numerics::QuadratureType quadratureType)
{
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  SLIC_ERROR_IF(sampleRes.size() != static_cast<axom::IndexType>(dim),
                "Sample resolution dimension does not match mesh dimension");

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return nullptr;
  }

  mfem::IntegrationRule ird[3];
  for(int d = 0; d < dim; d++)
  {
    SLIC_ERROR_IF(sampleRes[d] < 1,
                  axom::fmt::format("Invalid sample value {} for dimension {}.", sampleRes[d], d));
    switch(quadratureType)
    {
    case axom::numerics::QuadratureType::GaussLegendre:
      mfem::QuadratureFunctions1D::GaussLegendre(sampleRes[d], &ird[d]);
      break;
    case axom::numerics::QuadratureType::GaussLobatto:
      mfem::QuadratureFunctions1D::GaussLobatto(sampleRes[d], &ird[d]);
      break;
    case axom::numerics::QuadratureType::OpenUniform:
      mfem::QuadratureFunctions1D::OpenUniform(sampleRes[d], &ird[d]);
      break;
    case axom::numerics::QuadratureType::ClosedUniform:
      mfem::QuadratureFunctions1D::ClosedUniform(sampleRes[d], &ird[d]);
      break;
    case axom::numerics::QuadratureType::OpenHalfUniform:
      mfem::QuadratureFunctions1D::OpenHalfUniform(sampleRes[d], &ird[d]);
      break;
    case axom::numerics::QuadratureType::ClosedGL:
      mfem::QuadratureFunctions1D::ClosedGL(sampleRes[d], &ird[d]);
      break;
    case axom::numerics::QuadratureType::Invalid:
    default:
      SLIC_ERROR(axom::fmt::format("Invalid quadrature type {}.", static_cast<int>(quadratureType)));
      break;
    }
  }
  std::unique_ptr<mfem::IntegrationRule> ir;
  if(dim == 1)
  {
    ir = std::make_unique<mfem::IntegrationRule>(ird[0]);
  }
  else if(dim == 2)
  {
    ir = std::make_unique<mfem::IntegrationRule>(ird[0], ird[1]);
  }
  else if(dim == 3)
  {
    ir = std::make_unique<mfem::IntegrationRule>(ird[0], ird[1], ird[2]);
  }

  return new OwnedQuadratureSpace(*mesh, std::move(ir));
}

void assembleVolumeFractionRHS(const mfem::FiniteElementSpace& fes,
                               mfem::QuadratureFunction& inout,
                               const mfem::IntegrationRule& sampleIR,
                               bool useAnisotropicAssembly,
                               mfem::Vector& b)
{
  mfem::QuadratureFunctionCoefficient qfc(inout);
  mfem::DomainLFIntegrator rhs(qfc, &sampleIR);

  if(useAnisotropicAssembly)
  {
    mfem::Vector elemVec;
    mfem::Array<int> elemVDofs;
    const int NE = fes.GetNE();
    for(int elem = 0; elem < NE; ++elem)
    {
      rhs.AssembleRHSElementVect(*fes.GetFE(elem), *fes.GetElementTransformation(elem), elemVec);
      fes.GetElementVDofs(elem, elemVDofs);
      b.AddElementVector(elemVDofs, elemVec);
    }
  }
  else
  {
    mfem::Array<int> elem_marker(fes.GetNE());
    elem_marker.HostWrite();
    elem_marker = 1;
    elem_marker.ReadWrite();
    rhs.AssembleDevice(fes, elem_marker, b);
  }
}

}  // namespace

void generatePositionsQFunction(mfem::Mesh* mesh,
                                QFunctionCollection& inoutQFuncs,
                                axom::ArrayView<int> sampleResolution,
                                axom::numerics::QuadratureType quadratureType)
{
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return;
  }

  mfem::QuadratureSpace* sp = nullptr;
  if(quadratureType == axom::numerics::QuadratureType::Invalid)
  {
    SLIC_ERROR_IF(sampleResolution.empty(), "Invalid sampleResolution.");
    sp = makeDefaultQuadratureSpace(mesh, sampleResolution[0]);
  }
  else
  {
    sp = makeCustomQuadratureSpace(mesh, sampleResolution, quadratureType);
  }
  SLIC_ERROR_IF(sp == nullptr, "Null QuadratureSpace.");

  const auto& ir = sp->GetElementIntRule(0);
  const int nq = ir.GetNPoints();

  auto* pos_coef = new mfem::QuadratureFunction(sp, dim);
  pos_coef->SetOwnsSpace(true);
  auto pos = mfem::Reshape(pos_coef->HostWrite(), dim, nq, NE);

  if(!usesAnisotropicCustomTensorQuadrature(*mesh, sampleResolution, quadratureType))
  {
    const auto* geomFactors = mesh->GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);
    geomFactors->X.HostRead();

    for(int i = 0; i < NE; ++i)
    {
      const int gf_elStartIdx = i * nq * dim;
      for(int j = 0; j < dim; ++j)
      {
        for(int k = 0; k < nq; ++k)
        {
          pos(j, k, i) = geomFactors->X(gf_elStartIdx + (j * nq) + k);
        }
      }
    }

    mesh->DeleteGeometricFactors();
  }
  else
  {
    mfem::DenseMatrix pointMat(dim, nq);
    for(int i = 0; i < NE; ++i)
    {
      auto* transform = sp->GetTransformation(i);
      transform->Transform(ir, pointMat);

      for(int j = 0; j < dim; ++j)
      {
        for(int k = 0; k < nq; ++k)
        {
          pos(j, k, i) = pointMat(j, k);
        }
      }
    }
  }

  inoutQFuncs.Register("positions", pos_coef, true);
}

void generateSamplingPositions(SamplingMFEMState& mfemState,
                               axom::ArrayView<int> sampleResolution,
                               axom::numerics::QuadratureType quadratureType)
{
  AXOM_ANNOTATE_SCOPE("generateSamplingPositions");

  checkSampleResolution(mfemState, sampleResolution, quadratureType);

  if(mfemState.m_inoutShapeQFuncs.Has("positions"))
  {
    return;
  }

  generatePositionsQFunction(mfemState.m_dc->GetMesh(),
                             mfemState.m_inoutShapeQFuncs,
                             sampleResolution,
                             quadratureType);
}

void importInitialVolumeFractions(SamplingMFEMState& mfemState,
                                  const std::map<std::string, mfem::GridFunction*>& initialGridFunctions,
                                  bool anisotropic)
{
  auto* positionsQSpace = mfemState.shapeQFuncs().Get("positions")->GetSpace();
  auto* mesh = mfemState.m_dc->GetMesh();

  // Interpolate grid functions at quadrature points & register material quad functions
  // assume all elements have same integration rule
  for(auto& entry : initialGridFunctions)
  {
    const auto& name = entry.first;
    auto* gf = entry.second;

    SLIC_INFO_ROOT(axom::fmt::format("Importing volume fraction field for '{}' material", name));

    if(gf == nullptr)
    {
      SLIC_WARNING(
        axom::fmt::format("Skipping missing volume fraction field for material '{}'", name));
      continue;
    }

    auto* matQFunc = new mfem::QuadratureFunction(*positionsQSpace);
    const auto& ir = matQFunc->GetSpace()->GetIntRule(0);

    if(anisotropic)
    {
      // Avoid MFEM's tensor quadrature interpolation path only for
      // anisotropic custom quad/hex rules. MFEM infers a single q1d from
      // ir.GetNPoints(), which cannot represent per-direction sample counts
      // such as 3 x 5 or 3 x 5 x 2.
      mfem::Vector elemValues;
      mfem::Vector qfuncValues;
      for(int elem = 0; elem < mesh->GetNE(); ++elem)
      {
        gf->GetValues(elem, ir, elemValues);
        matQFunc->GetValues(elem, qfuncValues);
        qfuncValues = elemValues;
      }
    }
    else
    {
      const auto* interp = gf->FESpace()->GetQuadratureInterpolator(ir);
      SLIC_ERROR_IF(interp == nullptr,
                    axom::fmt::format("Could not create a quadrature interpolator while "
                                      "importing volume fractions for '{}'.",
                                      name));
      interp->Values(*gf, *matQFunc);
    }

    const auto matName = axom::fmt::format("mat_inout_{}", name);
    mfemState.materialQFuncs().Register(matName, matQFunc, true);
  }
}

void computeVolumeFractionsForMaterial(SamplingMFEMState& mfemState,
                                       const std::string& matField,
                                       int volfracOrder,
                                       axom::ArrayView<int> sampleResolution,
                                       axom::numerics::QuadratureType quadratureType)
{
  AXOM_ANNOTATE_SCOPE("computeVolumeFractionsForMaterial");

  SLIC_ASSERT(axom::utilities::string::startsWith(matField, "mat_inout_"));
  auto* inout = mfemState.getMaterialFunction(matField);
  SLIC_ASSERT(inout != nullptr);

  auto* dc = mfemState.m_dc;
  SLIC_ASSERT(dc != nullptr);

  const auto& sampleIR = inout->GetSpace()->GetIntRule(0);
  const int sampleOrder = sampleIR.GetOrder();
  const int sampleNQ = sampleIR.GetNPoints();
  const int sampleSZ = inout->GetSpace()->GetSize();

  mfem::Mesh* mesh = dc->GetMesh();
  const int dim = mesh->Dimension();
  const int NE = mesh->GetNE();

  auto samples_per_dim = [=](auto sampleRes, int dimValue) -> std::string {
    switch(dimValue)
    {
    case 2:
      return axom::fmt::format(" ({} * {})", sampleRes[0], sampleRes[1]);
    case 3:
      return axom::fmt::format(" ({} * {} * {})", sampleRes[0], sampleRes[1], sampleRes[2]);
    default:
      return std::string();
    }
  };

  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "In computeVolumeFractions(): num samples per element {}{} | "
                                   "sample polynomial order {} | total samples {:L}",
                                   sampleNQ,
                                   samples_per_dim(sampleResolution, dim),
                                   sampleOrder,
                                   sampleSZ));

  SLIC_INFO_ROOT(
    axom::fmt::format(axom::utilities::locale(), "Mesh has dim {} and {:L} elements", dim, NE));

  const auto vf_name = axom::fmt::format("vol_frac_{}", matField.substr(10));
  mfem::GridFunction* vf =
    getOrAllocateL2GridFunction(dc, vf_name, volfracOrder, dim, mfem::BasisType::Positive);
  const mfem::FiniteElementSpace* fes = vf->FESpace();
  const int dofs = fes->GetTypicalFE()->GetDof();

  mfem::DenseTensor* mass_mat {nullptr};
  const std::string mass_matrix_name = "shaping_mass_matrix";
  if(mfemState.m_inoutTensors.Has(mass_matrix_name))
  {
    mass_mat = mfemState.m_inoutTensors.Get(mass_matrix_name);
  }
  else
  {
    AXOM_ANNOTATE_SCOPE("mass integrator assemble");

    mass_mat = new mfem::DenseTensor(dofs, dofs, NE);
    mass_mat->HostWrite();
    (*mass_mat) = 0.;
    mass_mat->ReadWrite();

    mfem::ConstantCoefficient one_coef(1.0);
    mfem::MassIntegrator mass_integrator(one_coef, &sampleIR);

    if(usesAnisotropicCustomTensorQuadrature(*fes->GetMesh(), sampleResolution, quadratureType))
    {
      mfem::DenseMatrix elemMat;
      mass_mat->HostWrite();
      for(int elem = 0; elem < NE; ++elem)
      {
        mass_integrator.AssembleElementMatrix(*fes->GetFE(elem),
                                              *fes->GetElementTransformation(elem),
                                              elemMat);
        for(int j = 0; j < dofs; ++j)
        {
          for(int i = 0; i < dofs; ++i)
          {
            (*mass_mat)(i, j, elem) = elemMat(i, j);
          }
        }
      }
    }
    else
    {
      const int sz = mass_mat->TotalSize();
      mfem::Vector mass_vec;
      mfem::Swap(mass_mat->GetMemory(), mass_vec.GetMemory());
      mass_vec.SetSize(sz);
      mass_integrator.AssembleEA(*fes, mass_vec, false);
      mfem::Swap(mass_mat->GetMemory(), mass_vec.GetMemory());
    }

    mfemState.m_inoutTensors.Register(mass_matrix_name, mass_mat, true);
  }

  mfem::DenseTensor* mass_mat_inv {nullptr};
  mfem::Array<int>* mass_mat_pivots {nullptr};
  const std::string minv_name = "shaping_mass_matrix_inv";
  const std::string pivots_name = "shaping_mass_matrix_pivots";
  if(mfemState.m_inoutTensors.Has(minv_name) && mfemState.m_inoutArrays.Has(pivots_name))
  {
    mass_mat_inv = mfemState.m_inoutTensors.Get(minv_name);
    mass_mat_pivots = mfemState.m_inoutArrays.Get(pivots_name);
  }
  else
  {
    AXOM_ANNOTATE_SCOPE("batch lu factor");

    mass_mat->ReadWrite();
    mass_mat_inv = new mfem::DenseTensor(*mass_mat);
    mass_mat_pivots = new mfem::Array<int>(dofs * NE);

    mass_mat_inv->ReadWrite();
    mass_mat_pivots->Write();
    mfem::BatchLUFactor(*mass_mat_inv, *mass_mat_pivots);

    mfemState.m_inoutTensors.Register(minv_name, mass_mat_inv, true);
    mfemState.m_inoutArrays.Register(pivots_name, mass_mat_pivots, true);
  }

  mfem::DenseTensor* shaping_scratch_buffer {nullptr};
  const std::string scratch_buffer_name = "shaping_scratch_buffer";
  if(mfemState.m_inoutTensors.Has(scratch_buffer_name))
  {
    shaping_scratch_buffer = mfemState.m_inoutTensors.Get(scratch_buffer_name);
  }
  else
  {
    shaping_scratch_buffer = new mfem::DenseTensor(dofs, dofs, NE);
    shaping_scratch_buffer->HostWrite();
    (*shaping_scratch_buffer) = 0.;
    mfemState.m_inoutTensors.Register(scratch_buffer_name, shaping_scratch_buffer, true);
  }

  axom::utilities::Timer timer(true);
  {
    mfem::Vector b(fes->GetVSize());
    SLIC_ASSERT(b.Size() == dofs * NE);
    {
      AXOM_ANNOTATE_SCOPE("domain lf integrator assemble");

      inout->ReadWrite();
      b.HostWrite();
      b = 0.;
      b.ReadWrite();

      assembleVolumeFractionRHS(
        *fes,
        *inout,
        sampleIR,
        usesAnisotropicCustomTensorQuadrature(*fes->GetMesh(), sampleResolution, quadratureType),
        b);
    }
    inout->HostReadWrite();

    {
      AXOM_ANNOTATE_SCOPE("batch lu solve");

      mass_mat_inv->Read();
      mass_mat_pivots->Read();

      vf->HostReadWrite();
      (*vf) = b;
      vf->ReadWrite();
      mfem::BatchLUSolve(*mass_mat_inv, *mass_mat_pivots, *vf);
    }
    mass_mat_inv->HostReadWrite();
    mass_mat_pivots->HostReadWrite();

    constexpr double minY = 0.;
    constexpr double maxY = 1.;

    auto m_d = mfem::Reshape(mass_mat->HostReadWrite(), dofs, dofs, NE);
    auto b_d = mfem::Reshape(b.HostReadWrite(), dofs, NE);
    auto vf_d = mfem::Reshape(vf->HostReadWrite(), dofs, NE);
    auto fct_mat_d = mfem::Reshape(shaping_scratch_buffer->HostReadWrite(), dofs, dofs, NE);

    AXOM_ANNOTATE_BEGIN("fct project");
    axom::for_all<axom::SEQ_EXEC>(0, NE, [=](int i) {
      FCT_correct(&m_d(0, 0, i), dofs, &b_d(0, i), minY, maxY, &vf_d(0, i), &fct_mat_d(0, 0, i));
    });
    AXOM_ANNOTATE_END("fct project");
  }
  timer.stop();

  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "\t Generating volume fractions '{}' took {:.3f} seconds (@ "
                                   "{:L} dofs processed per second)",
                                   vf_name,
                                   timer.elapsed(),
                                   static_cast<int>(fes->GetNDofs() / timer.elapsed())));

  vf->HostReadWrite();
}

void FCT_correct(const double* M,
                 const int s,
                 const double* m,
                 const double y_min,
                 const double y_max,
                 double* xy,
                 double* fct_mat)
{
  constexpr int STACK_CAPACITY = 64;
  using StackArray = axom::StackArray<double, STACK_CAPACITY>;

  if(s == 1)
  {
    return;
  }

  StackArray ML_stack;
  StackArray z_stack;
  StackArray beta_stack;
  axom::Array<double> ML_heap;
  axom::Array<double> z_heap;
  axom::Array<double> beta_heap;

  double* ML = nullptr;
  double* z = nullptr;
  double* beta = nullptr;

  if(s <= STACK_CAPACITY)
  {
    ML = ML_stack.data();
    z = z_stack.data();
    beta = beta_stack.data();
  }
  else
  {
    ML_heap.resize(s);
    z_heap.resize(s);
    beta_heap.resize(s);

    ML = ML_heap.data();
    z = z_heap.data();
    beta = beta_heap.data();
  }

  for(int r = 0; r < s; ++r)
  {
    double dot = 0.;
    for(int c = 0; c < s; ++c)
    {
      dot += M[r + c * s];
    }
    ML[r] = dot;
  }

  double sum_ML = 0.;
  double sum_m = 0.;
  for(int i = 0; i < s; ++i)
  {
    sum_ML += ML[i];
    sum_m += m[i];
  }

  const double y_avg = sum_m / sum_ML;

  #ifdef AXOM_DEBUG
  constexpr double EPS = 1e-12;
  SLIC_WARNING_IF(
    !(y_min < y_avg + EPS && y_avg < y_max + EPS),
    axom::fmt::format("Average ({}) is out of bounds [{},{}]: ", y_avg, y_min - EPS, y_max + EPS));
  #endif

  double sum_beta = 0.;
  for(int i = 0; i < s; ++i)
  {
    beta[i] = ML[i];
    z[i] = m[i] - ML[i] * y_avg;
    sum_beta += beta[i];
  }

  for(int i = 0; i < s; ++i)
  {
    beta[i] /= sum_beta;
  }

  for(int i = 1; i < s; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      const int idx = i + j * s;
      fct_mat[idx] = M[idx] * (xy[i] - xy[j]) + (beta[j] * z[i] - beta[i] * z[j]);
    }
  }

  auto* gp = z;
  auto* gm = beta;
  for(int t = 0; t < s; ++t)
  {
    gp[t] = 0.0;
    gm[t] = 0.0;
  }

  for(int i = 1; i < s; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      const int idx = i + j * s;
      const double fij = fct_mat[idx];
      if(fij >= 0.0)
      {
        gp[i] += fij;
        gm[j] -= fij;
      }
      else
      {
        gm[i] += fij;
        gp[j] -= fij;
      }
    }
  }

  for(int i = 0; i < s; ++i)
  {
    xy[i] = y_avg;
  }

  for(int i = 0; i < s; ++i)
  {
    const double mi = ML[i];
    const double xyLi = xy[i];
    const double rp = axom::utilities::max(mi * (y_max - xyLi), 0.0);
    const double rm = axom::utilities::min(mi * (y_min - xyLi), 0.0);
    const double sp = gp[i];
    const double sm = gm[i];

    gp[i] = (rp < sp) ? rp / sp : 1.0;
    gm[i] = (rm > sm) ? rm / sm : 1.0;
  }

  for(int i = 1; i < s; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      double fij = fct_mat[i + j * s];

      const double aij =
        fij >= 0.0 ? axom::utilities::min(gp[i], gm[j]) : axom::utilities::min(gm[i], gp[j]);
      fij *= aij;
      xy[i] += fij / ML[i];
      xy[j] -= fij / ML[j];
    }
  }

  #ifdef AXOM_DEBUG
  for(int i = 0; i < s; ++i)
  {
    SLIC_WARNING_IF(!(y_min < xy[i] + EPS && xy[i] < y_max + EPS),
                    axom::fmt::format("Volume fraction {} w/ value {} is out of bounds [{},{}]: ",
                                      i,
                                      xy[i],
                                      y_min - EPS,
                                      y_max + EPS));
  }
  #endif
}

void computeVolumeFractionsIdentity(mfem::DataCollection* dc,
                                    mfem::QuadratureFunction* inout,
                                    const std::string& name)
{
  const int order = inout->GetSpace()->GetIntRule(0).GetOrder();

  mfem::Mesh* mesh = dc->GetMesh();
  const int dim = mesh->Dimension();
  const int NE = mesh->GetNE();

  std::cout << axom::fmt::format("Mesh has dim {} and {} elements", dim, NE) << std::endl;

  auto* fec = new mfem::L2_FECollection(order, dim, mfem::BasisType::Positive);
  auto* fes = new mfem::FiniteElementSpace(mesh, fec);
  auto* volFrac = new mfem::GridFunction(fes);
  volFrac->MakeOwner(fec);
  volFrac->HostReadWrite();
  dc->RegisterField(name, volFrac);

  (*volFrac) = (*inout);
}

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#endif  // defined(AXOM_USE_MFEM)
