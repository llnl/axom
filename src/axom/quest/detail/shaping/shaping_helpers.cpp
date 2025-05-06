// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "shaping_helpers.hpp"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#include "axom/fmt.hpp"

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

    // allocate a (par) finite element space on the mesh
    mfem::FiniteElementSpace* fes = nullptr;
    {
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
      auto* pmesh = dynamic_cast<mfem::ParMesh*>(mesh);
      if(pmesh)
      {
        fes = new mfem::ParFiniteElementSpace(pmesh, fec);
      }
      else
      {
        fes = new mfem::FiniteElementSpace(mesh, fec);
      }
  #else
      fes = new mfem::FiniteElementSpace(mesh, fec);
  #endif
    }

    // allocate data through sidre and tell the grid function to use it
    // the grid function will manage memory for the fec and fes
    auto* sidreDC = dynamic_cast<sidre::MFEMSidreDataCollection*>(dc);
    if(sidreDC)
    {
      const int sz = fes->GetVSize();
      auto* vw = sidreDC->AllocNamedBuffer(gf_name, sz);
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
      gf = new mfem::ParGridFunction();
  #else
      gf = new mfem::GridFunction();
  #endif
      gf->MakeRef(fes, vw->getData());
    }
    else
    {
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
      auto* pfes = dynamic_cast<mfem::ParFiniteElementSpace*>(fes);
      if(pfes)
      {
        gf = new mfem::ParGridFunction(pfes);
      }
      else
      {
        gf = new mfem::GridFunction(fes);
      }
  #else
      gf = new mfem::GridFunction(fes);
  #endif
    }

    gf->MakeOwner(fec);
    gf->HostReadWrite();

    dc->RegisterField(gf_name, gf);
    *gf = 0.;
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

/// Generates a quadrature function corresponding to the mesh "positions" field
void generatePositionsQFunction(mfem::Mesh* mesh, QFunctionCollection& inoutQFuncs, int sampleRes)
{
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return;
  }

  // convert requested samples into a compatible polynomial order
  // that will use that many samples: 2n-1 and 2n-2 will work
  // NOTE: Might be different for simplices
  const int sampleOrder = 2 * sampleRes - 1;
  mfem::QuadratureSpace* sp = new mfem::QuadratureSpace(mesh, sampleOrder);

  // TODO: Should the samples be along a uniform grid
  //       instead of Guassian quadrature?
  //       This would need quadrature weights for the uniform
  //       samples -- Newton-Cotes ?
  //       With uniform points, we could do HO polynomial fitting
  //       Using 0s and 1s is non-oscillatory in Bernstein basis

  // Assume all elements have the same integration rule
  const auto& ir = sp->GetElementIntRule(0);
  const int nq = ir.GetNPoints();
  const auto* geomFactors = mesh->GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);
  geomFactors->X.HostRead();

  mfem::QuadratureFunction* pos_coef = new mfem::QuadratureFunction(sp, dim);
  pos_coef->SetOwnsSpace(true);
  pos_coef->HostReadWrite();

  // Rearrange positions into quadrature function
  {
    for(int i = 0; i < NE; ++i)
    {
      const int elStartIdx = i * nq * dim;
      for(int j = 0; j < dim; ++j)
      {
        for(int k = 0; k < nq; ++k)
        {
          //X has dims nqpts x sdim x ne
          (*pos_coef)(elStartIdx + (k * dim) + j) = geomFactors->X(elStartIdx + (j * nq) + k);
        }
      }
    }
  }

  // register positions with the QFunction collection, which wil handle its deletion
  inoutQFuncs.Register("positions", pos_coef, true);
}

/** 
 * Implements flux-corrected transport (FCT) to convert the inout samples (ones and zeros)
 * to a grid function on the degrees of freedom such that the volume fractions are doubles
 * between 0 and 1 ( \a y_min and \a y_max )
 */
void FCT_project(const double* M,
                 const double* M_inv,
                 const int* P,
                 const int s,
                 const double* m,
                 const double* x,  // indicators
                 const double y_min,
                 const double y_max,
                 double* xy,
                 double* fct_mat)  // use as scratch buffer
{
  // [IN]  - M, M_inv, P, s, m, x, y_min, y_max
  // [OUT] - xy

  constexpr int ND = 64;
  using StackArray = axom::StackArray<double, ND>;
  SLIC_ASSERT(s <= ND);

  // Compute the lumped mass matrix in ML:  GetRowSums(M, s, s, ML);
  StackArray ML;
  for(int r = 0; r < s; ++r)
  {
    double dot {0.};
    for(int c = 0; c < s; ++c)
    {
      dot += M[r + c * s];
    }
    ML[r] = dot;
  }

  // Compute the high-order projection in xy: M_inv.Mult(m, xy);
  {
    for(int t = 0; t < s; ++t)
    {
      xy[t] = m[t];
    }

    // xy <- P xy
    for(int i = 0; i < s; ++i)
    {
      axom::utilities::swap(xy[i], xy[P[i]]);
    }

    // xy <- L^{-1} xy
    for(int j = 0; j < s; ++j)
    {
      const double x_j = xy[j];
      for(int i = j + 1; i < s; ++i)
      {
        xy[i] -= M_inv[i + j * s] * x_j;
      }
    }

    // xy <- U^{-1} xy
    for(int j = s - 1; j >= 0; --j)
    {
      const double x_j = (xy[j] /= M_inv[j + j * s]);
      for(int i = 0; i < j; ++i)
      {
        xy[i] -= M_inv[i + j * s] * x_j;
      }
    }
  }

  // Q0 solutions can't be adjusted conservatively. It's what it is.
  if(s == 1)
  {
    return;
  }

  double dMLX = 0.;
  double mSum = 0.;
  for(int i = 0; i < s; ++i)
  {
    dMLX += ML[i] * x[i];
    mSum += m[i];
  }

  const double y_avg = mSum / dMLX;

  #ifdef AXOM_DEBUG
  constexpr double EPS = 1e-12;
  if(!(y_min < y_avg + EPS && y_avg < y_max + EPS))
  {
    SLIC_WARNING(
      axom::fmt::format("Average ({}) is out of bounds [{},{}]: ", y_avg, y_min - EPS, y_max + EPS));
  }
  #endif

  StackArray z;
  StackArray beta;
  double betaSum = 0.;
  for(int i = 0; i < s; ++i)
  {
    // Some different options for beta:
    //beta[i] = 1.0;
    beta[i] = ML[i] * x[i];
    //beta[i] = ML[i]*(x[i] + 1e-14);
    //beta[i] = ML[i];

    // The low order flux correction
    z[i] = m[i] - ML[i] * x[i] * y_avg;
    betaSum += beta[i];
  }

  // Make beta_i sum to 1
  for(int i = 0; i < s; ++i)
  {
    beta[i] /= betaSum;
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
  auto& gp = z;
  auto& gm = beta;
  for(int t = 0; t < s; ++t)
  {
    gp[t] = 0.0;
    gm[t] = 0.0;
  };

  for(int i = 1; i < s; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      const int idx = i + j * s;
      double fij = fct_mat[idx];
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
    xy[i] = x[i] * y_avg;
  }

  for(int i = 0; i < s; ++i)
  {
    const double mi = ML[i];
    const double xyLi = xy[i];
    const double rp = axom::utilities::max(mi * (x[i] * y_max - xyLi), 0.0);
    const double rm = axom::utilities::min(mi * (x[i] * y_min - xyLi), 0.0);
    const double sp = gp[i];
    const double sm = gm[i];

    gp[i] = (rp < sp) ? rp / sp : 1.0;
    gm[i] = (rm > sm) ? rm / sm : 1.0;
  }

  for(int i = 1; i < s; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      const int idx = i + j * s;
      double fij = fct_mat[idx];

      const double aij =
        fij >= 0.0 ? axom::utilities::min(gp[j], gm[j]) : axom::utilities::min(gm[i], gp[j]);
      fij *= aij;
      xy[i] += fij / ML[i];
      xy[j] -= fij / ML[j];
    }
  }
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
