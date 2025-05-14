// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file shaping_helpers.hpp
 *
 * \brief Free-standing helper functions in support of shaping query
 */

#ifndef AXOM_QUEST_SHAPING_HELPERS__HPP_
#define AXOM_QUEST_SHAPING_HELPERS__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"

#if defined(AXOM_USE_MFEM)
  #include "mfem.hpp"
#endif

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

  /**
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

  /**
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

  /**
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

/// Alias to function pointer that projects a \a FromDim dimensional input point to
/// a \a ToDim dimensional query point when sampling the InOut field
template <int FromDim, int ToDim>
using PointProjector =
  axom::function<primal::Point<double, ToDim>(const primal::Point<double, FromDim>&)>;

#if defined(AXOM_USE_MFEM)

using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;
using MFEMArrayCollection = mfem::NamedFieldsMap<mfem::Array<int>>;

enum class VolFracSampling : int
{
  SAMPLE_AT_DOFS,
  SAMPLE_AT_QPTS
};

/**
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

/**
 * Utility function to zero out inout quadrature points for a material replaced by a shape
 *
 * Each location in space can only be covered by one material.
 * When \a shouldReplace is true, we clear all values in \a materialQFunc 
 * that are set in \a shapeQFunc. When it is false, we do the opposite.
 *
 * \param shapeQFunc The inout quadrature function for the shape samples
 * \param materialQFunc The inout quadrature function for the material samples
 * \param shapeReplacesMaterial Flag for whether the shape replaces the material 
 *   or whether the material remains and we should zero out the shape sample (when false)
 */
void replaceMaterial(mfem::QuadratureFunction* shapeQFunc,
                     mfem::QuadratureFunction* materialQFunc,
                     bool shouldReplace);

/**
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

/// Generates a quadrature function corresponding to the mesh positions
void generatePositionsQFunction(mfem::Mesh* mesh, QFunctionCollection& inoutQFuncs, int sampleRes);

/** 
 * Implements flux-corrected transport (FCT) to convert the inout samples (ones and zeros)
 * to a grid function on the degrees of freedom such that the volume fractions are doubles
 * between 0 and 1 ( \a y_min and \a y_max )
 */
void FCT_project(const double* M,
                 const int s,
                 const double* m,
                 const double y_min,  // 0
                 const double y_max,  // 1
                 double* xy,
                 double* fct_mat);  // scratch buffer

/**
 * \brief Identity transform for volume fractions from inout samples
 *
 * Copies \a inout samples from the quadrature function directly into volume fraction DOFs.
 * \param dc The data collection to which we will add the volume fractions
 * \param inout The inout samples
 * \param name The name of the generated volume fraction function
 * \note Assumes that the inout samples are co-located with the grid function DOFs.
 */
void computeVolumeFractionsIdentity(mfem::DataCollection* dc,
                                    mfem::QuadratureFunction* inout,
                                    const std::string& name);

#endif  // defined(AXOM_USE_MFEM)

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SHAPING_HELPERS__HPP_
