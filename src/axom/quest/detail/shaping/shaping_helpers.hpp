// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file shaping_helpers.hpp
 *
 * \brief Common shaping helper utilities and backend-specific helper facade
 */

#ifndef AXOM_QUEST_SHAPING_HELPERS__HPP_
#define AXOM_QUEST_SHAPING_HELPERS__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/numerics/quadrature.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/slic.hpp"

#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>

namespace axom
{

template <typename Signature, size_t MaxSize = 16>
class function;

/**
 * \brief Basic implementation of a host/device compatible analogue to std::function
 *
 * \tparam R The return type of the callable object
 * \tparam Args The parameter types of the callable object
 * \tparam MaxSize The maximum size of the callable (including its captured variables)
 *
 * \note We will extend this and move it to the core component
 */
template <typename R, typename... Args, size_t MaxSize>
class function<R(Args...), MaxSize>
{
private:
  using Storage = typename std::aligned_storage<MaxSize>::type;

public:
  AXOM_HOST_DEVICE function() : invoke(nullptr) { }

  template <typename Callable>
  AXOM_HOST_DEVICE function(Callable callable)
  {
    static_assert(sizeof(Callable) <= MaxSize, "Callable object too large!");
    static_assert(std::is_trivially_copyable<Callable>::value,
                  "Callable must be trivially copyable!");

    invoke = [](const void* storage, Args... args) -> R {
      return (*reinterpret_cast<const Callable*>(storage))(std::forward<Args>(args)...);
    };
    new(&storage) Callable(std::move(callable));
  }

  AXOM_HOST_DEVICE R operator()(Args... args) const
  {
    if(!invoke)
    {
      return R();
    }
    return invoke(&storage, std::forward<Args>(args)...);
  }

  AXOM_HOST_DEVICE explicit operator bool() const { return invoke != nullptr; }

private:
  Storage storage;

  R (*invoke)(const void*, Args... args) = nullptr;
};

template <typename Lambda>
auto make_host_device_function(Lambda&& lambda)
{
  using Signature = decltype(&Lambda::operator());
  return function<Signature>(std::forward<Lambda>(lambda));
}

namespace quest
{

// clang-format off
using seq_exec = axom::SEQ_EXEC;

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  using omp_exec = axom::OMP_EXEC;
#else
  using omp_exec = seq_exec;
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  constexpr int CUDA_BLOCK_SIZE = 256;
  using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
#else
  using cuda_exec = seq_exec;
#endif

#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  constexpr int HIP_BLOCK_SIZE = 64;
  using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
#else
  using hip_exec = seq_exec;
#endif
// clang-format on

namespace shaping
{

template <int FromDim, int ToDim>
using PointProjector =
  axom::function<primal::Point<double, ToDim>(const primal::Point<double, FromDim>&)>;

enum class VolFracSampling : int
{
  SAMPLE_AT_DOFS,
  SAMPLE_AT_QPTS
};

std::string shapeInOutFieldName(const std::string& shapeName);
std::string materialInOutFieldName(const std::string& materialName);
std::string volumeFractionFieldName(const std::string& materialName);
std::string shapeVolumeFractionFieldName(const std::string& shapeName);

std::string materialNameFromMaterialInOutFieldName(const std::string& fieldName);
std::string materialNameFromVolumeFractionFieldName(const std::string& fieldName);

template <typename MeshState>
void checkSampleResolution(const MeshState& meshState,
                           axom::ArrayView<int> sampleResolution,
                           axom::numerics::QuadratureType quadratureType)
{
  SLIC_ERROR_IF(quadratureType != axom::numerics::QuadratureType::Invalid &&
                  sampleResolution.size() != meshState.meshDimension(),
                "Inconsistent mesh dimension and sample resolutions.");
}

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#if defined(AXOM_USE_MFEM)
  #include "shaping_helpers_mfem.hpp"
#endif

#if defined(AXOM_USE_CONDUIT)
  #include "shaping_helpers_blueprint.hpp"
#endif

#endif  // AXOM_QUEST_SHAPING_HELPERS__HPP_
