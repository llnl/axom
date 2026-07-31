// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * \file shaping_helpers.hpp
 *
 * \brief Common shaping helper utilities and backend-specific helper facade
 */

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

/*!
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

  /*!
   * \brief Constructs a function object from a callable object
   *
   * \tparam Callable The type of the callable object
   * \param callable The callable object to store and invoke
   *
   * This constructor stores the callable object in the internal storage
   * and sets up the invoke function pointer to call the stored object.
   * The callable object must be trivially copyable and its size must not
   * exceed the maximum storage size.
   */
  template <typename Callable>
  AXOM_HOST_DEVICE function(Callable callable)
  {
    static_assert(sizeof(Callable) <= MaxSize, "Callable object too large!");
    static_assert(std::is_trivially_copyable<Callable>::value,
                  "Callable must be trivially copyable!");
    //SLIC_WARNING("sizeof(Callable): " << sizeof(Callable));

    invoke = [](const void* storage, Args... args) -> R {
      return (*reinterpret_cast<const Callable*>(storage))(std::forward<Args>(args)...);
    };
    new(&storage) Callable(std::move(callable));
  }

  /*!
   * \brief invoke the stored callable object
   *
   * \param args The arguments to be forwarded to the callable object
   *
   * \return The result of invoking the callable object with the provided arguments.
   *         If the callable object is not set (i.e., `invoke` is null), a default-constructed
   *         value of type R is returned.
   */
  AXOM_HOST_DEVICE R operator()(Args... args) const
  {
    if(!invoke)
    {
      return R();
    }
    return invoke(&storage, std::forward<Args>(args)...);
  }

  /*!
   * \brief Explicit conversion operator to check the validity of the object
   *
   * \return True if `invoke` is not null, false otherwise
   */
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

/*!
 * \brief Return the registered shape in/out field name for a shape.
 *
 * \param shapeName The shape name.
 *
 * \return The corresponding shape in/out field name.
 */
std::string shapeInOutFieldName(const std::string& shapeName);

/*!
 * \brief Return the registered material in/out field name for a material.
 *
 * \param materialName The material name.
 *
 * \return The corresponding material in/out field name.
 */
std::string materialInOutFieldName(const std::string& materialName);

/*!
 * \brief Return the registered volume-fraction field name for a material.
 *
 * \param materialName The material name.
 *
 * \return The corresponding volume-fraction field name.
 */
std::string volumeFractionFieldName(const std::string& materialName);

/*!
 * \brief Return the per-shape volume-fraction field name for a shape.
 *
 * \param shapeName The shape name.
 *
 * \return The corresponding per-shape volume-fraction field name.
 */
std::string shapeVolumeFractionFieldName(const std::string& shapeName);

/*!
 * \brief Return whether a field name is a registered shape in/out field name.
 *
 * \param fieldName The field name to inspect.
 *
 * \return True if the field name belongs to the shape in/out family.
 */
bool isShapeInOutFieldName(const std::string& fieldName);

/*!
 * \brief Return whether a field name is a registered material in/out field name.
 *
 * \param fieldName The field name to inspect.
 *
 * \return True if the field name belongs to the material in/out family.
 */
bool isMaterialInOutFieldName(const std::string& fieldName);

/*!
 * \brief Return whether a field name is a registered material volume-fraction field name.
 *
 * \param fieldName The field name to inspect.
 *
 * \return True if the field name belongs to the material volume-fraction family.
 */
bool isVolumeFractionFieldName(const std::string& fieldName);

/*!
 * \brief Return whether a field name is a registered per-shape volume-fraction field name.
 *
 * \param fieldName The field name to inspect.
 *
 * \return True if the field name belongs to the per-shape volume-fraction family.
 */
bool isShapeVolumeFractionFieldName(const std::string& fieldName);

/*!
 * \brief Extract the material name from a material in/out field name.
 *
 * \param fieldName The field name to inspect.
 *
 * \return The material name, or an empty string if the field name does not match.
 */
std::string materialNameFromMaterialInOutFieldName(const std::string& fieldName);

/*!
 * \brief Extract the material name from a volume-fraction field name.
 *
 * \param fieldName The field name to inspect.
 *
 * \return The material name, or an empty string if the field name does not match.
 */
std::string materialNameFromVolumeFractionFieldName(const std::string& fieldName);

/*!
 * \brief Checks that input sample resolution array view is appropriate for the
 *        mesh dimension and quadrature type.
 *
 * \tparam MeshState A type that contains mesh state.
 *
 * \param meshState The mesh state
 * \param sampleResolution An array view that contains the sample resolutions in each dimension.
 * \param quadratureType The quadrature type.
 */
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
