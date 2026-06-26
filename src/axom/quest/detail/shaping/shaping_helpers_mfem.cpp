// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "shaping_helpers_mfem.hpp"

#if defined(AXOM_USE_MFEM)

  #include <algorithm>
  #include <cstdint>
  #include <limits>
  #include <memory>
  #include <vector>

  #include "mfem/linalg/kernels.hpp"

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

struct VolumeFractionMassConfig
{
  int dofs {};
  int numElements {};
  bool usesAnisotropicQuadrature {false};
  bool useChunkedMassProcessing {false};
  std::int64_t elemTensorEntries {};
  std::int64_t cachedMassBytes {};
};

axom::runtime_policy::Policy selectFCTExecutionPolicy(axom::runtime_policy::Policy execPolicy,
                                                      const std::string& vfName)
{
  using RuntimePolicy = axom::runtime_policy::Policy;

  switch(execPolicy)
  {
  case RuntimePolicy::seq:
    return RuntimePolicy::seq;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    return RuntimePolicy::omp;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA) || defined(AXOM_RUNTIME_POLICY_USE_HIP)
    SLIC_WARNING_ROOT(axom::fmt::format(
      "FCT projection for '{}' uses MFEM host data and currently falls back to sequential execution for device runtime policies.",
      vfName));
    return RuntimePolicy::seq;
#endif
  default:
    SLIC_WARNING_ROOT(axom::fmt::format(
      "FCT projection for '{}' falls back to sequential execution because the requested runtime policy is not available in this build.",
      vfName));
    return RuntimePolicy::seq;
  }
}

void assembleVolumeFractionRHS(const mfem::FiniteElementSpace& fes,
                               mfem::QuadratureFunction& inout,
                               const mfem::IntegrationRule& sampleIR,
                               bool usesAnisotropicQuadrature,
                               mfem::Vector& rhs);

std::string formatSamplesPerDimension(axom::ArrayView<int> sampleResolution, int dim)
{
  switch(dim)
  {
  case 2:
    return axom::fmt::format(" ({} * {})", sampleResolution[0], sampleResolution[1]);
  case 3:
    return axom::fmt::format(" ({} * {} * {})",
                             sampleResolution[0],
                             sampleResolution[1],
                             sampleResolution[2]);
  default:
    return std::string();
  }
}

void logVolumeFractionInputs(int sampleNQ,
                             int sampleOrder,
                             int sampleSZ,
                             int dim,
                             int numElements,
                             axom::ArrayView<int> sampleResolution)
{
  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "In computeVolumeFractions(): num samples per element {}{} | "
                                   "sample polynomial order {} | total samples {:L}",
                                   sampleNQ,
                                   formatSamplesPerDimension(sampleResolution, dim),
                                   sampleOrder,
                                   sampleSZ));

  SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                   "Mesh has dim {} and {:L} elements",
                                   dim,
                                   numElements));
}

VolumeFractionMassConfig makeVolumeFractionMassConfig(const mfem::FiniteElementSpace& fes,
                                                      axom::ArrayView<int> sampleResolution,
                                                      axom::numerics::QuadratureType quadratureType)
{
  constexpr std::int64_t MAX_CACHED_MASS_BYTES = 1LL << 30;

  VolumeFractionMassConfig config;
  config.dofs = fes.GetTypicalFE()->GetDof();
  config.numElements = fes.GetMesh()->GetNE();
  config.usesAnisotropicQuadrature =
    usesAnisotropicCustomTensorQuadrature(*fes.GetMesh(), sampleResolution, quadratureType);
  config.elemTensorEntries = static_cast<std::int64_t>(config.dofs) * config.dofs;

  const std::int64_t totalTensorEntries = config.elemTensorEntries * config.numElements;
  config.cachedMassBytes = totalTensorEntries * 3 * sizeof(double);
  config.useChunkedMassProcessing = totalTensorEntries > std::numeric_limits<int>::max() ||
    config.cachedMassBytes > MAX_CACHED_MASS_BYTES;

  return config;
}

void logChunkedMassProcessing(const std::string& vfName, const VolumeFractionMassConfig& config)
{
  if(!config.useChunkedMassProcessing)
  {
    return;
  }

  SLIC_INFO_ROOT(
    axom::fmt::format(axom::utilities::locale(),
                      "Using chunked local mass assembly for '{}' (dofs={}, elements={:L}, "
                      "dense cache would require {:.2f} GiB)",
                      vfName,
                      config.dofs,
                      config.numElements,
                      static_cast<double>(config.cachedMassBytes) / (1024.0 * 1024.0 * 1024.0)));
}

void assembleVolumeFractionRHSVector(const mfem::FiniteElementSpace& fes,
                                     mfem::QuadratureFunction& inout,
                                     const mfem::IntegrationRule& sampleIR,
                                     bool usesAnisotropicQuadrature,
                                     mfem::Vector& rhs)
{
  AXOM_ANNOTATE_SCOPE("domain lf integrator assemble");

  inout.ReadWrite();
  rhs.HostWrite();
  rhs = 0.;
  rhs.ReadWrite();

  assembleVolumeFractionRHS(fes, inout, sampleIR, usesAnisotropicQuadrature, rhs);
  inout.HostReadWrite();
}

mfem::DenseTensor* getOrAssembleMassMatrix(MFEMState& mfemState,
                                           const mfem::FiniteElementSpace& fes,
                                           const mfem::IntegrationRule& sampleIR,
                                           const VolumeFractionMassConfig& config)
{
  const std::string massMatrixName = "shaping_mass_matrix";
  if(mfemState.m_inoutTensors.Has(massMatrixName))
  {
    return mfemState.m_inoutTensors.Get(massMatrixName);
  }

  AXOM_ANNOTATE_SCOPE("mass integrator assemble");

  auto* massMat = new mfem::DenseTensor(config.dofs, config.dofs, config.numElements);
  massMat->HostWrite();
  (*massMat) = 0.;
  massMat->ReadWrite();

  mfem::ConstantCoefficient one_coef(1.0);
  mfem::MassIntegrator mass_integrator(one_coef, &sampleIR);

  if(config.usesAnisotropicQuadrature)
  {
    mfem::DenseMatrix elemMat;
    massMat->HostWrite();
    for(int elem = 0; elem < config.numElements; ++elem)
    {
      mass_integrator.AssembleElementMatrix(*fes.GetFE(elem),
                                            *fes.GetElementTransformation(elem),
                                            elemMat);
      for(int j = 0; j < config.dofs; ++j)
      {
        for(int i = 0; i < config.dofs; ++i)
        {
          (*massMat)(i, j, elem) = elemMat(i, j);
        }
      }
    }
  }
  else
  {
    const int sz = massMat->TotalSize();
    mfem::Vector mass_vec;
    mfem::Swap(massMat->GetMemory(), mass_vec.GetMemory());
    mass_vec.SetSize(sz);
    mass_integrator.AssembleEA(fes, mass_vec, false);
    mfem::Swap(massMat->GetMemory(), mass_vec.GetMemory());
  }

  mfemState.m_inoutTensors.Register(massMatrixName, massMat, true);
  return massMat;
}

std::pair<mfem::DenseTensor*, mfem::Array<int>*> getOrFactorMassMatrix(
  MFEMState& mfemState,
  mfem::DenseTensor& massMat,
  const VolumeFractionMassConfig& config)
{
  const std::string minvName = "shaping_mass_matrix_inv";
  const std::string pivotsName = "shaping_mass_matrix_pivots";
  if(mfemState.m_inoutTensors.Has(minvName) && mfemState.m_inoutArrays.Has(pivotsName))
  {
    return {mfemState.m_inoutTensors.Get(minvName), mfemState.m_inoutArrays.Get(pivotsName)};
  }

  AXOM_ANNOTATE_SCOPE("batch lu factor");

  massMat.ReadWrite();
  auto* massMatInv = new mfem::DenseTensor(massMat);
  auto* massMatPivots = new mfem::Array<int>(config.dofs * config.numElements);

  massMatInv->ReadWrite();
  massMatPivots->Write();
  mfem::BatchLUFactor(*massMatInv, *massMatPivots);

  mfemState.m_inoutTensors.Register(minvName, massMatInv, true);
  mfemState.m_inoutArrays.Register(pivotsName, massMatPivots, true);

  return {massMatInv, massMatPivots};
}

mfem::DenseTensor* getOrAllocateScratchBuffer(MFEMState& mfemState,
                                              const VolumeFractionMassConfig& config)
{
  const std::string scratchBufferName = "shaping_scratch_buffer";
  if(mfemState.m_inoutTensors.Has(scratchBufferName))
  {
    return mfemState.m_inoutTensors.Get(scratchBufferName);
  }

  auto* scratchBuffer = new mfem::DenseTensor(config.dofs, config.dofs, config.numElements);
  scratchBuffer->HostWrite();
  (*scratchBuffer) = 0.;
  mfemState.m_inoutTensors.Register(scratchBufferName, scratchBuffer, true);

  return scratchBuffer;
}

void applyFCTProjection(mfem::DenseTensor& massMat,
                        axom::runtime_policy::Policy execPolicy,
                        mfem::Vector& rhs,
                        const int dofs,
                        const int numElements,
                        mfem::Vector& vfData,
                        mfem::DenseTensor& scratchBuffer);

template <typename ExecSpace>
void applyFCTProjectionImpl(mfem::DenseTensor& massMat,
                            mfem::Vector& rhs,
                            const int dofs,
                            const int numElements,
                            mfem::Vector& vfData,
                            mfem::DenseTensor& scratchBuffer)
{
  constexpr double minY = 0.;
  constexpr double maxY = 1.;

  AXOM_ANNOTATE_SCOPE("fct project");

  auto m_d = mfem::Reshape(massMat.HostReadWrite(), dofs, dofs, numElements);
  auto b_d = mfem::Reshape(rhs.HostReadWrite(), dofs, numElements);
  auto vf_d = mfem::Reshape(vfData.HostReadWrite(), dofs, numElements);
  auto fct_mat_d = mfem::Reshape(scratchBuffer.HostReadWrite(), dofs, dofs, numElements);

  axom::for_all<ExecSpace>(0, numElements, [=](int elem) {
    FCT_correct(&m_d(0, 0, elem), dofs, &b_d(0, elem), minY, maxY, &vf_d(0, elem), &fct_mat_d(0, 0, elem));
  });
}

void applyFCTProjection(mfem::DenseTensor& massMat,
                        axom::runtime_policy::Policy execPolicy,
                        mfem::Vector& rhs,
                        const int dofs,
                        const int numElements,
                        mfem::Vector& vfData,
                        mfem::DenseTensor& scratchBuffer)
{
  using RuntimePolicy = axom::runtime_policy::Policy;

  switch(execPolicy)
  {
  case RuntimePolicy::seq:
    applyFCTProjectionImpl<seq_exec>(massMat, rhs, dofs, numElements, vfData, scratchBuffer);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    applyFCTProjectionImpl<omp_exec>(massMat, rhs, dofs, numElements, vfData, scratchBuffer);
    break;
#endif
  default:
    applyFCTProjectionImpl<seq_exec>(massMat, rhs, dofs, numElements, vfData, scratchBuffer);
    break;
  }
}

void solveVolumeFractionsCached(MFEMState& mfemState,
                                const mfem::FiniteElementSpace& fes,
                                const mfem::IntegrationRule& sampleIR,
                                const VolumeFractionMassConfig& config,
                                axom::runtime_policy::Policy execPolicy,
                                mfem::Vector& rhs,
                                mfem::GridFunction& vf)
{
  auto* massMat = getOrAssembleMassMatrix(mfemState, fes, sampleIR, config);
  auto [massMatInv, massMatPivots] = getOrFactorMassMatrix(mfemState, *massMat, config);
  auto* scratchBuffer = getOrAllocateScratchBuffer(mfemState, config);

  {
    AXOM_ANNOTATE_SCOPE("batch lu solve");

    massMatInv->Read();
    massMatPivots->Read();

    vf.HostReadWrite();
    vf = rhs;
    vf.ReadWrite();
    mfem::BatchLUSolve(*massMatInv, *massMatPivots, vf);
  }
  massMatInv->HostReadWrite();
  massMatPivots->HostReadWrite();

  applyFCTProjection(*massMat, execPolicy, rhs, config.dofs, config.numElements, vf, *scratchBuffer);
}

int computeChunkSize(const VolumeFractionMassConfig& config)
{
  constexpr std::int64_t TARGET_CHUNK_BYTES = 256LL * 1024 * 1024;
  const std::int64_t bytesPerElement = (3 * config.elemTensorEntries * sizeof(double)) +
    (static_cast<std::int64_t>(config.dofs) * sizeof(int));
  const std::int64_t elemsPerChunk =
    std::max<std::int64_t>(1, TARGET_CHUNK_BYTES / std::max<std::int64_t>(1, bytesPerElement));

  return static_cast<int>(std::min<std::int64_t>(config.numElements, elemsPerChunk));
}

/// Serial version that does not have to create MFEM objects each iteration.
void assembleChunkMassMatricesSequential(const mfem::FiniteElementSpace& fes,
                                         mfem::MassIntegrator& massIntegrator,
                                         const VolumeFractionMassConfig& config,
                                         int elemBegin,
                                         int chunkNE,
                                         mfem::DenseTensor& massMat,
                                         mfem::DenseTensor& massMatInv)
{
  mfem::DenseMatrix elemMat;
  massMat.HostWrite();
  massMatInv.HostWrite();

  for(int elem = 0; elem < chunkNE; ++elem)
  {
    const int globalElem = elemBegin + elem;
    massIntegrator.AssembleElementMatrix(*fes.GetFE(globalElem),
                                         *fes.GetElementTransformation(globalElem),
                                         elemMat);
    for(int j = 0; j < config.dofs; ++j)
    {
      for(int i = 0; i < config.dofs; ++i)
      {
        const double value = elemMat(i, j);
        massMat(i, j, elem) = value;
        massMatInv(i, j, elem) = value;
      }
    }
  }
}

/// OpenMP parallelizable version that builds MFEM objects in the for_all.
template <typename ExecSpace>
void assembleChunkMassMatricesImpl(const mfem::FiniteElementSpace& fes,
                                   const mfem::IntegrationRule& sampleIR,
                                   const VolumeFractionMassConfig& config,
                                   int elemBegin,
                                   int chunkNE,
                                   mfem::DenseTensor& massMat,
                                   mfem::DenseTensor& massMatInv)
{
  auto massMat_d = mfem::Reshape(massMat.HostWrite(), config.dofs, config.dofs, chunkNE);
  auto massMatInv_d = mfem::Reshape(massMatInv.HostWrite(), config.dofs, config.dofs, chunkNE);
  auto* mesh = fes.GetMesh();

  axom::for_all<ExecSpace>(0, chunkNE, [=, &fes, &sampleIR](int elem) {
    mfem::ConstantCoefficient oneCoef(1.0);
    mfem::MassIntegrator massIntegrator(oneCoef, &sampleIR);
    mfem::DenseMatrix elemMat;
    mfem::IsoparametricTransformation tr;

    const int globalElem = elemBegin + elem;
    mesh->GetElementTransformation(globalElem, &tr);
    massIntegrator.AssembleElementMatrix(*fes.GetFE(globalElem), tr, elemMat);

    for(int j = 0; j < config.dofs; ++j)
    {
      for(int i = 0; i < config.dofs; ++i)
      {
        const double value = elemMat(i, j);
        massMat_d(i, j, elem) = value;
        massMatInv_d(i, j, elem) = value;
      }
    }
  });
}

void assembleChunkMassMatrices(const mfem::FiniteElementSpace& fes,
                               const mfem::IntegrationRule& sampleIR,
                               mfem::MassIntegrator& massIntegrator,
                               const VolumeFractionMassConfig& config,
                               axom::runtime_policy::Policy execPolicy,
                               int elemBegin,
                               int chunkNE,
                               mfem::DenseTensor& massMat,
                               mfem::DenseTensor& massMatInv)
{
  using RuntimePolicy = axom::runtime_policy::Policy;

  switch(execPolicy)
  {
  case RuntimePolicy::seq:
    assembleChunkMassMatricesSequential(
      fes, massIntegrator, config, elemBegin, chunkNE, massMat, massMatInv);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    assembleChunkMassMatricesImpl<omp_exec>(
      fes, sampleIR, config, elemBegin, chunkNE, massMat, massMatInv);
    break;
#endif
  default:
    assembleChunkMassMatricesSequential(
      fes, massIntegrator, config, elemBegin, chunkNE, massMat, massMatInv);
    break;
  }
}

void copyChunkRHS(mfem::Vector& rhs,
                  const VolumeFractionMassConfig& config,
                  int elemBegin,
                  int chunkNE,
                  mfem::Vector& rhsChunk,
                  mfem::Vector& vfChunk)
{
  auto rhs_d = mfem::Reshape(rhs.HostReadWrite(), config.dofs, config.numElements);
  auto rhs_chunk_d = mfem::Reshape(rhsChunk.HostWrite(), config.dofs, chunkNE);
  auto vf_chunk_d = mfem::Reshape(vfChunk.HostWrite(), config.dofs, chunkNE);

  for(int elem = 0; elem < chunkNE; ++elem)
  {
    const int globalElem = elemBegin + elem;
    for(int i = 0; i < config.dofs; ++i)
    {
      const double value = rhs_d(i, globalElem);
      rhs_chunk_d(i, elem) = value;
      vf_chunk_d(i, elem) = value;
    }
  }
}

template <typename ExecSpace>
void factorChunkMassMatricesImpl(const VolumeFractionMassConfig& config,
                                 int chunkNE,
                                 mfem::DenseTensor& massMatInv,
                                 mfem::Array<int>& massMatPivots)
{
  auto massMatInv_d = mfem::Reshape(massMatInv.HostReadWrite(), config.dofs, config.dofs, chunkNE);
  auto massMatPivots_d = mfem::Reshape(massMatPivots.Write(), config.dofs, chunkNE);

  axom::for_all<ExecSpace>(0, chunkNE, [=](int elem) {
    mfem::kernels::LUFactor(&massMatInv_d(0, 0, elem), config.dofs, &massMatPivots_d(0, elem));
  });
}

void factorChunkMassMatrices(const VolumeFractionMassConfig& config,
                             axom::runtime_policy::Policy execPolicy,
                             int chunkNE,
                             mfem::DenseTensor& massMatInv,
                             mfem::Array<int>& massMatPivots)
{
  using RuntimePolicy = axom::runtime_policy::Policy;

  switch(execPolicy)
  {
  case RuntimePolicy::seq:
    massMatInv.ReadWrite();
    massMatPivots.Write();
    mfem::BatchLUFactor(massMatInv, massMatPivots);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    factorChunkMassMatricesImpl<omp_exec>(config, chunkNE, massMatInv, massMatPivots);
    break;
#endif
  default:
    massMatInv.ReadWrite();
    massMatPivots.Write();
    mfem::BatchLUFactor(massMatInv, massMatPivots);
    break;
  }
}

void copyChunkResultToGridFunction(const VolumeFractionMassConfig& config,
                                   int elemBegin,
                                   int chunkNE,
                                   mfem::Vector& vfChunk,
                                   mfem::GridFunction& vf)
{
  auto vf_chunk_d = mfem::Reshape(vfChunk.HostReadWrite(), config.dofs, chunkNE);
  auto vf_d = mfem::Reshape(vf.HostReadWrite(), config.dofs, config.numElements);

  for(int elem = 0; elem < chunkNE; ++elem)
  {
    const int globalElem = elemBegin + elem;
    for(int i = 0; i < config.dofs; ++i)
    {
      vf_d(i, globalElem) = vf_chunk_d(i, elem);
    }
  }
}

void solveVolumeFractionsChunked(const mfem::FiniteElementSpace& fes,
                                 const mfem::IntegrationRule& sampleIR,
                                 const VolumeFractionMassConfig& config,
                                 axom::runtime_policy::Policy execPolicy,
                                 mfem::Vector& rhs,
                                 mfem::GridFunction& vf)
{
  AXOM_ANNOTATE_SCOPE("chunked mass solve");

  mfem::ConstantCoefficient one_coef(1.0);
  mfem::MassIntegrator massIntegrator(one_coef, &sampleIR);
  const int chunkSize = computeChunkSize(config);

  for(int elemBegin = 0; elemBegin < config.numElements; elemBegin += chunkSize)
  {
    const int chunkNE = std::min(chunkSize, config.numElements - elemBegin);

    mfem::DenseTensor massMat(config.dofs, config.dofs, chunkNE);
    mfem::DenseTensor massMatInv(config.dofs, config.dofs, chunkNE);
    mfem::DenseTensor scratchBuffer(config.dofs, config.dofs, chunkNE);
    mfem::Array<int> massMatPivots(config.dofs * chunkNE);
    mfem::Vector rhsChunk(config.dofs * chunkNE);
    mfem::Vector vfChunk(config.dofs * chunkNE);

    {
      AXOM_ANNOTATE_SCOPE("chunk mass assembly");
      assembleChunkMassMatrices(fes,
                                sampleIR,
                                massIntegrator,
                                config,
                                execPolicy,
                                elemBegin,
                                chunkNE,
                                massMat,
                                massMatInv);
    }

    {
      AXOM_ANNOTATE_SCOPE("chunk rhs copy");
      copyChunkRHS(rhs, config, elemBegin, chunkNE, rhsChunk, vfChunk);
    }

    {
      AXOM_ANNOTATE_SCOPE("chunk batch lu factor");
      factorChunkMassMatrices(config, execPolicy, chunkNE, massMatInv, massMatPivots);
    }

    {
      AXOM_ANNOTATE_SCOPE("chunk batch lu solve");
      massMatInv.Read();
      massMatPivots.Read();
      vfChunk.ReadWrite();
      mfem::BatchLUSolve(massMatInv, massMatPivots, vfChunk);
    }

    applyFCTProjection(massMat, execPolicy, rhsChunk, config.dofs, chunkNE, vfChunk, scratchBuffer);
    {
      AXOM_ANNOTATE_SCOPE("chunk result copy");
      copyChunkResultToGridFunction(config, elemBegin, chunkNE, vfChunk, vf);
    }
  }
}

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

void printRegisteredFieldNames(const MFEMState& mfemState,
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

void generateSamplingPositions(MFEMState& mfemState,
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

void importInitialVolumeFractions(MFEMState& mfemState,
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

    const auto matName = shaping::materialInOutFieldName(name);
    mfemState.materialQFuncs().Register(matName, matQFunc, true);
  }
}

void computeVolumeFractionsForMaterial(MFEMState& mfemState,
                                       const std::string& matField,
                                       int volfracOrder,
                                       axom::ArrayView<int> sampleResolution,
                                       axom::numerics::QuadratureType quadratureType,
                                       axom::runtime_policy::Policy execPolicy)
{
  AXOM_ANNOTATE_SCOPE("computeVolumeFractionsForMaterial");

  const std::string materialName = shaping::materialNameFromMaterialInOutFieldName(matField);
  SLIC_ASSERT(!materialName.empty());
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
  logVolumeFractionInputs(sampleNQ, sampleOrder, sampleSZ, dim, NE, sampleResolution);

  const auto vf_name = shaping::volumeFractionFieldName(materialName);
  mfem::GridFunction* vf =
    getOrAllocateL2GridFunction(dc, vf_name, volfracOrder, dim, mfem::BasisType::Positive);
  const mfem::FiniteElementSpace* fes = vf->FESpace();
  const VolumeFractionMassConfig config =
    makeVolumeFractionMassConfig(*fes, sampleResolution, quadratureType);
  logChunkedMassProcessing(vf_name, config);
  const auto fctExecPolicy = selectFCTExecutionPolicy(execPolicy, vf_name);

  axom::utilities::Timer timer(true);
  {
    mfem::Vector b(fes->GetVSize());
    SLIC_ASSERT(b.Size() == config.dofs * config.numElements);
    assembleVolumeFractionRHSVector(*fes, *inout, sampleIR, config.usesAnisotropicQuadrature, b);

    if(config.useChunkedMassProcessing)
    {
      solveVolumeFractionsChunked(*fes, sampleIR, config, fctExecPolicy, b, *vf);
    }
    else
    {
      solveVolumeFractionsCached(mfemState, *fes, sampleIR, config, fctExecPolicy, b, *vf);
    }
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
