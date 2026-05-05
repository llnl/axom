// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "shaping_helpers.hpp"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#include "axom/fmt.hpp"

#include <memory>

#if defined(AXOM_USE_MFEM)
  #include "mfem/linalg/dtensor.hpp"
#endif

namespace axom
{
namespace quest
{
namespace shaping
{
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

}  // namespace

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

mfem::QuadratureSpace* makeCustomQuadratureSpace(mfem::Mesh* mesh, int sampleRes[3], int quadratureType)
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
    case mfem::Quadrature1D::GaussLegendre:
      mfem::QuadratureFunctions1D::GaussLegendre(sampleRes[d], &ird[d]);
      break;
    case mfem::Quadrature1D::GaussLobatto:
      mfem::QuadratureFunctions1D::GaussLobatto(sampleRes[d], &ird[d]);
      break;
    case mfem::Quadrature1D::OpenUniform:
      mfem::QuadratureFunctions1D::OpenUniform(sampleRes[d], &ird[d]);
      break;
    case mfem::Quadrature1D::ClosedUniform:
      mfem::QuadratureFunctions1D::ClosedUniform(sampleRes[d], &ird[d]);
      break;
    case mfem::Quadrature1D::OpenHalfUniform:
      mfem::QuadratureFunctions1D::OpenHalfUniform(sampleRes[d], &ird[d]);
      break;
    case mfem::Quadrature1D::ClosedGL:
      mfem::QuadratureFunctions1D::ClosedGL(sampleRes[d], &ird[d]);
      break;
    default:
      SLIC_ERROR(axom::fmt::format("Invalid quadrature type {}.", quadratureType));
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

/// Generates a quadrature function corresponding to the mesh "positions" field
void generatePositionsQFunction(mfem::Mesh* mesh,
                                QFunctionCollection& inoutQFuncs,
                                int sampleResolution[3],
                                int quadratureType)
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
  if(quadratureType == static_cast<int>(mfem::Quadrature1D::Invalid))
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

  if(quadratureType == static_cast<int>(mfem::Quadrature1D::Invalid))
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
    // points in each logical dimension. For custom anisotropic tensor-product
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
