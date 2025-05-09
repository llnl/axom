// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "shaping_helpers.hpp"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
namespace shaping
{
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

/// Generate a volume fraction from a quadrature field for \a matField using FCT
void computeVolumeFractions(const std::string& matField,
                            mfem::DataCollection* dc,
                            QFunctionCollection& inoutQFuncs,
                            int outputOrder)
{
  SLIC_ASSERT(axom::utilities::string::startsWith(matField, "mat_inout_"));

  const auto volFracName = axom::fmt::format("vol_frac_{}", matField.substr(10));

  // Grab a pointer to the inout samples QFunc
  mfem::QuadratureFunction* inout = inoutQFuncs.Get(matField);

  const int sampleOrder = inout->GetSpace()->GetIntRule(0).GetOrder();
  const int sampleNQ = inout->GetSpace()->GetIntRule(0).GetNPoints();
  const int sampleSZ = inout->GetSpace()->GetSize();
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "In computeVolumeFractions(): sample order {} | "
                              "sample num qpts {} |  total samples {:L}",
                              sampleOrder,
                              sampleNQ,
                              sampleSZ));

  mfem::Mesh* mesh = dc->GetMesh();
  const int dim = mesh->Dimension();
  const int NE = mesh->GetNE();

  SLIC_INFO(
    axom::fmt::format(axom::utilities::locale(), "Mesh has dim {} and {:L} elements", dim, NE));

  // Access or create a registered volume fraction grid function from the data collection
  mfem::FiniteElementSpace* fes = nullptr;
  mfem::GridFunction* volFrac = nullptr;
  if(dc->HasField(volFracName))
  {
    volFrac = dc->GetField(volFracName);
    fes = volFrac->FESpace();
  }
  else
  {
    const auto basis = mfem::BasisType::Positive;
    auto* fec = new mfem::L2_FECollection(outputOrder, dim, basis);
    fes = new mfem::FiniteElementSpace(mesh, fec);
    volFrac = new mfem::GridFunction(fes);
    volFrac->MakeOwner(fec);
    volFrac->HostReadWrite();

    dc->RegisterField(volFracName, volFrac);
  }

  // Project QField onto volume fractions field using flux corrected transport (FCT)
  // to keep the range of values between 0 and 1
  axom::utilities::Timer timer(true);
  {
    mfem::MassIntegrator mass_integrator;
    mfem::QuadratureFunctionCoefficient qfc(*inout);
    mfem::DomainLFIntegrator rhs(qfc);

    const auto& ir = inout->GetSpace()->GetIntRule(0);  // assume all elements are the same
    rhs.SetIntRule(&ir);

    mfem::DenseMatrix m;
    mfem::DenseMatrixInverse mInv;
    mfem::Vector b, x;
    mfem::Array<int> dofs;
    mfem::Vector one;
    const double minY = 0.;
    const double maxY = 1.;

    for(int i = 0; i < NE; ++i)
    {
      auto* T = mesh->GetElementTransformation(i);
      auto* el = fes->GetFE(i);

      mass_integrator.AssembleElementMatrix(*el, *T, m);
      rhs.AssembleRHSElementVect(*el, *T, b);
      mInv.Factor(m);

      if(one.Size() != b.Size())
      {
        one.SetSize(b.Size());
        one = 1.0;
      }
      FCT_project(m, mInv, b, one, minY, maxY, x);

      fes->GetElementDofs(i, dofs);
      volFrac->SetSubVector(dofs, x);
    }
  }
  timer.stop();
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "\t Generating volume fractions '{}' took {:.3f} seconds (@ "
                              "{:L} dofs processed per second)",
                              volFracName,
                              timer.elapsed(),
                              static_cast<int>(fes->GetNDofs() / timer.elapsed())));
}

/** 
 * Implements flux-corrected transport (FCT) to convert the inout samples (ones and zeros)
 * to a grid function on the degrees of freedom such that the volume fractions are doubles
 * between 0 and 1 ( \a y_min and \a y_max )
 */
void FCT_project(mfem::DenseMatrix& M,
                 mfem::DenseMatrixInverse& M_inv,
                 mfem::Vector& m,
                 mfem::Vector& x,  // indicators
                 double y_min,
                 double y_max,
                 mfem::Vector& xy)  // indicators * rho
{
  // [IN]  - M, M_inv, m, x, y_min, y_max
  // [OUT] - xy

  using namespace mfem;

  const int s = M.Size();

  xy.SetSize(s);

  // Compute the high-order projection in xy
  M_inv.Mult(m, xy);

  // Q0 solutions can't be adjusted conservatively. It's what it is.
  if(xy.Size() == 1)
  {
    return;
  }

  // Compute the lumped mass matrix in ML
  Vector ML(s);
  M.GetRowSums(ML);

  //Ensure dot product is done on the CPU
  double dMLX(0);
  for(int i = 0; i < x.Size(); ++i)
  {
    dMLX += ML(i) * x(i);
  }

  const double y_avg = m.Sum() / dMLX;

#ifdef AXOM_DEBUG
  SLIC_WARNING_IF(
    !(y_min < y_avg + 1e-12 && y_avg < y_max + 1e-12),
    axom::fmt::format("Average ({}) is out of bounds [{},{}]: ", y_avg, y_min - 1e-12, y_max + 1e-12));
#endif

  Vector z(s);
  Vector beta(s);
  Vector Mxy(s);
  M.Mult(xy, Mxy);
  for(int i = 0; i < s; i++)
  {
    // Some different options for beta:
    //beta(i) = 1.0;
    beta(i) = ML(i) * x(i);
    //beta(i) = ML(i)*(x(i) + 1e-14);
    //beta(i) = ML(i);
    //beta(i) = Mxy(i);

    // The low order flux correction
    z(i) = m(i) - ML(i) * x(i) * y_avg;
  }

  // Make beta_i sum to 1
  beta /= beta.Sum();

  DenseMatrix F(s);
  // Note: indexing F(i,j) where  0 <= j < i < s
  for(int i = 1; i < s; i++)
  {
    for(int j = 0; j < i; j++)
    {
      F(i, j) = M(i, j) * (xy(i) - xy(j)) + (beta(j) * z(i) - beta(i) * z(j));
    }
  }

  Vector gp(s), gm(s);
  gp = 0.0;
  gm = 0.0;
  for(int i = 1; i < s; i++)
  {
    for(int j = 0; j < i; j++)
    {
      double fij = F(i, j);
      if(fij >= 0.0)
      {
        gp(i) += fij;
        gm(j) -= fij;
      }
      else
      {
        gm(i) += fij;
        gp(j) -= fij;
      }
    }
  }

  for(int i = 0; i < s; i++)
  {
    xy(i) = x(i) * y_avg;
  }

  for(int i = 0; i < s; i++)
  {
    double mi = ML(i), xyLi = xy(i);
    double rp = std::max(mi * (x(i) * y_max - xyLi), 0.0);
    double rm = std::min(mi * (x(i) * y_min - xyLi), 0.0);
    double sp = gp(i), sm = gm(i);

    gp(i) = (rp < sp) ? rp / sp : 1.0;
    gm(i) = (rm > sm) ? rm / sm : 1.0;
  }

  for(int i = 1; i < s; i++)
  {
    for(int j = 0; j < i; j++)
    {
      double fij = F(i, j), aij;

      if(fij >= 0.0)
      {
        aij = std::min(gp(i), gm(j));
      }
      else
      {
        aij = std::min(gm(i), gp(j));
      }

      fij *= aij;
      xy(i) += fij / ML(i);
      xy(j) -= fij / ML(j);
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

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom
