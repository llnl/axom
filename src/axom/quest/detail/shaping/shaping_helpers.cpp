// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "shaping_helpers.hpp"
#include "GenerateQuadratureMesh.hpp"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/numerics/quadrature.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#include "axom/fmt.hpp"

#include <memory>

#if defined(AXOM_USE_CONDUIT)
  #include "axom/bump/views/dispatch_coordset.hpp"
  #include "axom/bump/views/dispatch_unstructured_topology.hpp"
#endif

#if defined(AXOM_USE_MFEM)
  #include "mfem/linalg/dtensor.hpp"
#endif

namespace axom
{
namespace quest
{
namespace shaping
{
#if defined(AXOM_USE_CONDUIT)
namespace
{

constexpr const char* QUADRATURE_COORDSET_NAME = "quadrature_points";
constexpr const char* QUADRATURE_TOPOLOGY_NAME = "quadrature_points";
constexpr const char* ORIGINAL_ELEMENTS_FIELD_NAME = "originalElements";
constexpr const char* QUADRATURE_WEIGHTS_FIELD_NAME = "quadratureWeights";

numerics::QuadratureRule getBlueprintQuadratureRule(axom::numerics::QuadratureType quadratureType,
                                                    int npts,
                                                    int allocatorID)
{
  SLIC_ERROR_IF(npts < 1, axom::fmt::format("Invalid sample resolution {}.", npts));
  SLIC_ERROR_IF(!axom::numerics::is_supported_quadrature_type(quadratureType),
                axom::fmt::format(
                  "Quadrature type {} is not yet supported for Blueprint quadrature meshes.",
                  static_cast<int>(quadratureType)));

  return numerics::get_quadrature_rule(quadratureType, npts, allocatorID);
}

template <typename ExecSpace, typename CoordsetView>
void buildBlueprintQuadratureMesh(const conduit::Node& topoNode,
                                  const conduit::Node& coordsetNode,
                                  const CoordsetView& coordsetView,
                                  int allocatorID,
                                  const numerics::QuadratureRule& ruleX,
                                  const numerics::QuadratureRule& ruleY,
                                  const numerics::QuadratureRule& ruleZ,
                                  conduit::Node& meshNode)
{
  namespace views = axom::bump::views;
  constexpr int SupportedShapes = views::select_shapes(views::Quad_ShapeID, views::Hex_ShapeID);

  views::dispatch_unstructured_topology<SupportedShapes>(topoNode, [&](const auto&, auto topoView) {
    GenerateQuadratureMesh<ExecSpace, decltype(topoView), CoordsetView> generator(topoView,
                                                                                   coordsetView);
    generator.setAllocatorID(allocatorID);
    generator.execute(topoNode,
                      coordsetNode,
                      QUADRATURE_TOPOLOGY_NAME,
                      QUADRATURE_COORDSET_NAME,
                      ORIGINAL_ELEMENTS_FIELD_NAME,
                      QUADRATURE_WEIGHTS_FIELD_NAME,
                      ruleX,
                      ruleY,
                      ruleZ,
                      meshNode);
  });
}

}  // namespace
#endif

#if defined(AXOM_USE_MFEM)

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

bool usesAnisotropicCustomTensorQuadrature(const mfem::Mesh& mesh,
                                           const int sampleResolution[3],
                                           axom::numerics::QuadratureType quadratureType)
{
  if(quadratureType == axom::numerics::QuadratureType::Invalid)
  {
    return false;
  }

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

}  // namespace

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

// Utility function to either return a gf from the dc, or to allocate it through the dc
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

    // allocate data through sidre and tell the grid function to use it
    // the grid function will manage memory for the fec and fes
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
    // If shapeReplacesMaterial, clear material samples that are inside current shape
    for(int j = 0; j < SZ; ++j)
    {
      mData[j] = sData[j] > 0 ? 0 : mData[j];
    }
  }
  else
  {
    // Otherwise, clear current shape samples that are in the material
    for(int j = 0; j < SZ; ++j)
    {
      sData[j] = mData[j] > 0 ? 0 : sData[j];
    }
  }
}

/// Utility function to copy in_out quadrature samples from one QFunc to another
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

  // When reuseExisting, don't reset material values; otherwise, just copy values over
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

mfem::QuadratureSpace* makeDefaultQuadratureSpace(mfem::Mesh* mesh, int sampleRes)
{
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return nullptr;
  }

  // convert requested samples into a compatible polynomial order
  // that will use that many samples: 2n-1 and 2n-2 will work
  // NOTE: Might be different for simplices
  const int sampleOrder = 2 * sampleRes - 1;
  return new mfem::QuadratureSpace(mesh, sampleOrder);
}

mfem::QuadratureSpace* makeCustomQuadratureSpace(mfem::Mesh* mesh,
                                                 int sampleRes[3],
                                                 axom::numerics::QuadratureType quadratureType)
{
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return nullptr;
  }

  // Make custom integration rule
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

    for(int elem = 0; elem < fes.GetNE(); ++elem)
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

/// Generates a quadrature function corresponding to the mesh "positions" field
void generatePositionsQFunction(mfem::Mesh* mesh,
                                QFunctionCollection& inoutQFuncs,
                                int sampleResolution[3],
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

  // Make a quadrature space to determine the point locations in each element.
  mfem::QuadratureSpace* sp = nullptr;
  if(quadratureType == axom::numerics::QuadratureType::Invalid)
  {
    sp = makeDefaultQuadratureSpace(mesh, sampleResolution[0]);
  }
  else
  {
    sp = makeCustomQuadratureSpace(mesh, sampleResolution, quadratureType);
  }
  SLIC_ERROR_IF(sp == nullptr, "Null QuadratureSpace.");

  // Assume all elements have the same integration rule
  const auto& ir = sp->GetElementIntRule(0);
  const int nq = ir.GetNPoints();

  mfem::QuadratureFunction* pos_coef = new mfem::QuadratureFunction(sp, dim);
  pos_coef->SetOwnsSpace(true);
  auto pos = mfem::Reshape(pos_coef->HostWrite(), dim, nq, NE);

  if(!usesAnisotropicCustomTensorQuadrature(*mesh, sampleResolution, quadratureType))
  {
    const auto* geomFactors = mesh->GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);
    geomFactors->X.HostRead();

    // Rearrange positions into quadrature function
    for(int i = 0; i < NE; ++i)
    {
      const int gf_elStartIdx = i * nq * dim;
      for(int j = 0; j < dim; ++j)
      {
        for(int k = 0; k < nq; ++k)
        {
          // X has dims nqpts x sdim x ne
          pos(j, k, i) = geomFactors->X(gf_elStartIdx + (j * nq) + k);
        }
      }
    }

    // Delete the geometric factors associated w/ our quadrature rule
    mesh->DeleteGeometricFactors();
  }
  else
  {
    // MFEM's tensor quadrature interpolation assumes the same number of
    // points in each logical dimension. For anisotropic custom tensor-product
    // rules, map the integration points explicitly through each element.
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

  // register positions with the QFunction collection, which will handle its deletion
  inoutQFuncs.Register("positions", pos_coef, true);
}

void generateSamplingPositions(SamplingMFEMState& mfemState,
                               int sampleResolution[3],
                               axom::numerics::QuadratureType quadratureType)
{
  if(mfemState.m_inoutShapeQFuncs.Has("positions"))
  {
    return;
  }

  generatePositionsQFunction(mfemState.m_dc->GetMesh(),
                             mfemState.m_inoutShapeQFuncs,
                             sampleResolution,
                             quadratureType);
}

void computeVolumeFractionsForMaterial(SamplingMFEMState& mfemState,
                                       const std::string& matField,
                                       int volfracOrder,
                                       int sampleResolution[3],
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
  const auto geom = mesh->GetTypicalElementGeometry();

  auto samples_per_dim = [=](int sampleRes[3], mfem::Geometry::Type geomType) -> std::string {
    switch(geomType)
    {
    case mfem::Geometry::SQUARE:
      return axom::fmt::format(" ({} * {})", sampleRes[0], sampleRes[1]);
    case mfem::Geometry::CUBE:
      return axom::fmt::format(" ({} * {} * {})", sampleRes[0], sampleRes[1], sampleRes[2]);
    default:
      return std::string();
    }
  };

  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "In computeVolumeFractions(): num samples per element {}{} | "
                                   "sample polynomial order {} | total samples {:L}",
                                   sampleNQ,
                                   samples_per_dim(sampleResolution, geom),
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

      assembleVolumeFractionRHS(*fes,
                                *inout,
                                sampleIR,
                                usesAnisotropicCustomTensorQuadrature(*fes->GetMesh(),
                                                                      sampleResolution,
                                                                      quadratureType),
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
      FCT_correct(&m_d(0, 0, i),
                  dofs,
                  &b_d(0, i),
                  minY,
                  maxY,
                  &vf_d(0, i),
                  &fct_mat_d(0, 0, i));
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

#if defined(AXOM_USE_CONDUIT)
void generateQuadraturePointMesh(conduit::Node& bpMeshNode,
                                 const std::string& topologyName,
                                 int allocatorID,
                                 int sampleResolution[3],
                                 axom::numerics::QuadratureType quadratureType)
{
  if(bpMeshNode.has_path(axom::fmt::format("topologies/{}", QUADRATURE_TOPOLOGY_NAME)))
  {
    return;
  }

  const conduit::Node& topoNode =
    bpMeshNode.fetch_existing("topologies").fetch_existing(topologyName);
  const std::string topoType = topoNode.fetch_existing("type").as_string();
  SLIC_ERROR_IF(topoType != "unstructured",
                axom::fmt::format("Unsupported Blueprint topology type '{}' for quadrature mesh generation.",
                                  topoType));

  const std::string shape = topoNode.fetch_existing("elements/shape").as_string();
  SLIC_ERROR_IF(shape != "quad" && shape != "hex",
                axom::fmt::format("Unsupported Blueprint element shape '{}' for quadrature mesh generation.",
                                  shape));

  const std::string coordsetName = topoNode.fetch_existing("coordset").as_string();
  const conduit::Node& coordsetNode = bpMeshNode.fetch_existing("coordsets").fetch_existing(coordsetName);
  const std::string coordsetType = coordsetNode.fetch_existing("type").as_string();
  SLIC_ERROR_IF(coordsetType != "explicit",
                axom::fmt::format("Unsupported Blueprint coordset type '{}' for quadrature mesh generation.",
                                  coordsetType));

  int selectedAllocatorID = allocatorID;
  if(!axom::execution_space<seq_exec>::usesAllocId(selectedAllocatorID) &&
     !axom::execution_space<omp_exec>::usesAllocId(selectedAllocatorID)
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
     && !axom::execution_space<cuda_exec>::usesAllocId(selectedAllocatorID)
#endif
#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
     && !axom::execution_space<hip_exec>::usesAllocId(selectedAllocatorID)
#endif
  )
  {
    selectedAllocatorID = axom::execution_space<seq_exec>::allocatorID();
  }

  auto ruleX = getBlueprintQuadratureRule(quadratureType, sampleResolution[0], selectedAllocatorID);
  auto ruleY = getBlueprintQuadratureRule(quadratureType, sampleResolution[1], selectedAllocatorID);
  auto ruleZ = getBlueprintQuadratureRule(quadratureType, sampleResolution[2], selectedAllocatorID);

  axom::bump::views::dispatch_explicit_coordset(coordsetNode, [&](auto coordsetView) {
#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    if(axom::execution_space<hip_exec>::usesAllocId(selectedAllocatorID))
    {
      buildBlueprintQuadratureMesh<hip_exec>(topoNode,
                                             coordsetNode,
                                             coordsetView,
                                             selectedAllocatorID,
                                             ruleX,
                                             ruleY,
                                             ruleZ,
                                             bpMeshNode);
      return;
    }
#endif
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    if(axom::execution_space<cuda_exec>::usesAllocId(selectedAllocatorID))
    {
      buildBlueprintQuadratureMesh<cuda_exec>(topoNode,
                                              coordsetNode,
                                              coordsetView,
                                              selectedAllocatorID,
                                              ruleX,
                                              ruleY,
                                              ruleZ,
                                              bpMeshNode);
      return;
    }
#endif
    if(axom::execution_space<omp_exec>::usesAllocId(selectedAllocatorID))
    {
      buildBlueprintQuadratureMesh<omp_exec>(topoNode,
                                             coordsetNode,
                                             coordsetView,
                                             selectedAllocatorID,
                                             ruleX,
                                             ruleY,
                                             ruleZ,
                                             bpMeshNode);
      return;
    }

    buildBlueprintQuadratureMesh<seq_exec>(topoNode,
                                           coordsetNode,
                                           coordsetView,
                                           selectedAllocatorID,
                                           ruleX,
                                           ruleY,
                                           ruleZ,
                                           bpMeshNode);
  });
}

void generateSamplingPositions(BlueprintState& bpState,
                               int sampleResolution[3],
                               axom::numerics::QuadratureType quadratureType)
{
  if(bpState.m_internal_node.has_path(
       axom::fmt::format("topologies/{}", QUADRATURE_TOPOLOGY_NAME)))
  {
    return;
  }

  generateQuadraturePointMesh(bpState.m_internal_node,
                              bpState.m_topology_name,
                              bpState.m_allocator_id,
                              sampleResolution,
                              quadratureType);
}

void computeVolumeFractionsForMaterial(BlueprintState& bpState, const std::string& matField)
{
  AXOM_ANNOTATE_SCOPE("computeVolumeFractionsForMaterial");

  SLIC_ASSERT(axom::utilities::string::startsWith(matField, "mat_inout_"));

  conduit::Node* inout = bpState.getMaterialFunction(matField);
  SLIC_ERROR_IF(inout == nullptr,
                axom::fmt::format("Missing Blueprint material field '{}' for volume fraction projection.",
                                  matField));

  conduit::Node& bpMeshNode = bpState.m_internal_node;
  SLIC_ERROR_IF(!bpMeshNode.has_path("fields/originalElements/values"),
                "Missing Blueprint originalElements field for volume fraction projection.");
  SLIC_ERROR_IF(!bpMeshNode.has_path("fields/quadratureWeights/values"),
                "Missing Blueprint quadratureWeights field for volume fraction projection.");

  const conduit::Node& topoNode =
    bpMeshNode.fetch_existing("topologies").fetch_existing(bpState.m_topology_name);

  axom::IndexType numZones = 0;
  axom::bump::views::dispatch_unstructured_topology(
    topoNode, [&](const auto&, auto topoView) { numZones = topoView.numberOfZones(); });

  namespace utils = axom::bump::utilities;
  const auto originalElements =
    utils::make_array_view<conduit::index_t>(bpMeshNode["fields/originalElements/values"]);
  const auto quadratureWeights =
    utils::make_array_view<double>(bpMeshNode["fields/quadratureWeights/values"]);
  const auto inoutValues = utils::make_array_view<double>(inout->fetch_existing("values"));

  SLIC_ASSERT(originalElements.size() == quadratureWeights.size());
  SLIC_ASSERT(originalElements.size() == inoutValues.size());

  const std::string vfName = axom::fmt::format("vol_frac_{}", matField.substr(10));
  conduit::Node& vfNode = bpMeshNode["fields/" + vfName];
  vfNode.reset();
  vfNode["association"] = "element";
  vfNode["topology"] = bpState.m_topology_name;

  const auto conduitAllocatorId =
    axom::sidre::ConduitMemory::axomAllocIdToConduit(bpState.m_allocator_id);
  conduit::Node& valuesNode = vfNode["values"];
  valuesNode.set_allocator(conduitAllocatorId);
  valuesNode.set(conduit::DataType::float64(numZones));
  auto vfValues = utils::make_array_view<double>(valuesNode);

  for(axom::IndexType zoneIdx = 0; zoneIdx < vfValues.size(); ++zoneIdx)
  {
    vfValues[zoneIdx] = 0.;
  }

  for(axom::IndexType pointIdx = 0; pointIdx < inoutValues.size(); ++pointIdx)
  {
    const conduit::index_t zoneIdx = originalElements[pointIdx];
    SLIC_ASSERT(zoneIdx >= 0);
    SLIC_ASSERT(zoneIdx < vfValues.size());
    vfValues[zoneIdx] += inoutValues[pointIdx] * quadratureWeights[pointIdx];
  }
}

void replaceMaterial(conduit::Node* shapeNode,
                     conduit::Node* materialNode,
                     bool shapeReplacesMaterial)
{
  SLIC_ASSERT(shapeNode != nullptr);
  SLIC_ASSERT(materialNode != nullptr);

  namespace utils = axom::bump::utilities;
  auto shapeValues = utils::make_array_view<double>(shapeNode->fetch_existing("values"));
  auto materialValues = utils::make_array_view<double>(materialNode->fetch_existing("values"));

  SLIC_ASSERT(shapeValues.size() == materialValues.size());

  for(axom::IndexType i = 0; i < materialValues.size(); ++i)
  {
    if(shapeReplacesMaterial)
    {
      materialValues[i] = shapeValues[i] > 0. ? 0. : materialValues[i];
    }
    else
    {
      shapeValues[i] = materialValues[i] > 0. ? 0. : shapeValues[i];
    }
  }
}

void copyShapeIntoMaterial(const conduit::Node* shapeNode,
                           conduit::Node* materialNode,
                           bool reuseExisting)
{
  SLIC_ASSERT(shapeNode != nullptr);
  SLIC_ASSERT(materialNode != nullptr);

  namespace utils = axom::bump::utilities;
  const auto shapeValues = utils::make_array_view<double>(shapeNode->fetch_existing("values"));
  auto materialValues = utils::make_array_view<double>(materialNode->fetch_existing("values"));

  SLIC_ASSERT(shapeValues.size() == materialValues.size());

  if(reuseExisting)
  {
    for(axom::IndexType i = 0; i < materialValues.size(); ++i)
    {
      materialValues[i] = shapeValues[i] > 0. ? 1. : materialValues[i];
    }
  }
  else
  {
    for(axom::IndexType i = 0; i < materialValues.size(); ++i)
    {
      materialValues[i] = shapeValues[i];
    }
  }
}

conduit::Node* cloneInOutFunction(const conduit::Node* node)
{
  SLIC_ASSERT(node != nullptr);
  return new conduit::Node(*node);
}
#endif

void FCT_correct(const double* M,     // Mass matrix
                 const int s,         // num dofs
                 const double* m,     // rhs (incorporating the inout samples)
                 const double y_min,  // lower bound for FCT
                 const double y_max,  // upper bound for FCt
                 double* xy,          // uncorrected volume fraction dofs
                 double* fct_mat)     // use as scratch buffer
{
  // [IN]  - M, s, m, y_min, y_max
  // [INOUT] - xy

  constexpr int STACK_CAPACITY = 64;
  using StackArray = axom::StackArray<double, STACK_CAPACITY>;

  // Q0 solutions can't be adjusted conservatively. It is what it is.
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

  // Compute the lumped mass matrix in ML:  M.GetRowSums(ML);
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
    // Some different options for beta:
    //beta[i] = 1.0;
    beta[i] = ML[i];
    //beta[i] = ML[i]*(1. + 1e-14);

    // The low order flux correction
    z[i] = m[i] - ML[i] * y_avg;
    sum_beta += beta[i];
  }

  // Make beta_i sum to 1
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

  // NOTE: `z' and `beta' are no longer used.
  // Zero them out and reuse their memory under different aliases: gp and gm
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
  // check that volume fractions are in bounds
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

// Note: This function is not currently being used, but might be in the near future
void computeVolumeFractionsIdentity(mfem::DataCollection* dc,
                                    mfem::QuadratureFunction* inout,
                                    const std::string& name)
{
  const int order = inout->GetSpace()->GetIntRule(0).GetOrder();

  mfem::Mesh* mesh = dc->GetMesh();
  const int dim = mesh->Dimension();
  const int NE = mesh->GetNE();

  std::cout << axom::fmt::format("Mesh has dim {} and {} elements", dim, NE) << std::endl;

  mfem::L2_FECollection* fec = new mfem::L2_FECollection(order, dim, mfem::BasisType::Positive);
  mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, fec);
  mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
  volFrac->MakeOwner(fec);
  volFrac->HostReadWrite();
  dc->RegisterField(name, volFrac);

  (*volFrac) = (*inout);
}

#endif  // defined(AXOM_USE_MFEM)

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom
