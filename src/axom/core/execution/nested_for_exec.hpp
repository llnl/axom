// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_NESTED_FOR_EXEC_HPP_
#define AXOM_CORE_NESTED_FOR_EXEC_HPP_

#include "axom/core/execution/execution_space.hpp"

// RAJA includes
#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"

  // NOTE: add RAJA alias for older versions of RAJA prior to RAJA-v0.11.0
  #if (RAJA_VERSION_MAJOR == 0) && (RAJA_VERSION_MINOR < 11)

/* clang-format off */

    namespace RAJA
    {
      template < int SIZE >
      using tile_fixed = ::RAJA::statement::tile_fixed< SIZE >;
    } // namespace RAJA

    /* clang-format on */

  #endif

#endif

namespace axom
{
namespace internal
{
template <typename ExecSpace>
struct nested_for_exec
{
  using loop2d_policy = void;
  using loop3d_policy = void;
};

//--------------------------------------------------------| SEQ_EXEC |----------
template <>
struct nested_for_exec<SEQ_EXEC>
{
#ifdef AXOM_USE_RAJA
  /* clang-format off */
  using loop2d_policy = RAJA::KernelPolicy<
#if RAJA_VERSION_MAJOR > 2022
    RAJA::statement::For< 1, RAJA::seq_exec,   // j
      RAJA::statement::For< 0, RAJA::seq_exec, // i
#else
    RAJA::statement::For< 1, RAJA::loop_exec,   // j
      RAJA::statement::For< 0, RAJA::loop_exec, // i
#endif
        RAJA::statement::Lambda< 0 >
      > // END i
    > // END j
  >; // END kernel

  using loop3d_policy = RAJA::KernelPolicy<
#if RAJA_VERSION_MAJOR > 2022
    RAJA::statement::For< 2, RAJA::seq_exec,     // k
      RAJA::statement::For< 1, RAJA::seq_exec,   // j
        RAJA::statement::For< 0, RAJA::seq_exec, // i
#else
    RAJA::statement::For< 2, RAJA::loop_exec,     // k
      RAJA::statement::For< 1, RAJA::loop_exec,   // j
        RAJA::statement::For< 0, RAJA::loop_exec, // i
#endif
          RAJA::statement::Lambda< 0 >
        > // END i
      > // END j
    > // END k
  >; // END kernel
    /* clang-format on */

#else
  using loop2d_policy = void;
  using loop3d_policy = void;
#endif
};

//--------------------------------------------------------| OMP_EXEC |----------
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
template <>
struct nested_for_exec<OMP_EXEC>
{
  /* clang-format off */

  using loop2d_policy = RAJA::KernelPolicy<
    RAJA::statement::For< 1, RAJA::omp_parallel_for_exec, // j
#if RAJA_VERSION_MAJOR > 2022
      RAJA::statement::For< 0, RAJA::seq_exec,            // i
#else
      RAJA::statement::For< 0, RAJA::loop_exec,           // i
#endif
        RAJA::statement::Lambda< 0 >
     > // END i
    > // END j
  >; // END kernel

  using loop3d_policy = RAJA::KernelPolicy<
    RAJA::statement::For< 2, RAJA::omp_parallel_for_exec, // k
#if RAJA_VERSION_MAJOR > 2022
      RAJA::statement::For< 1, RAJA::seq_exec,            // j
        RAJA::statement::For< 0, RAJA::seq_exec,          // i
#else
      RAJA::statement::For< 1, RAJA::loop_exec,           // j
        RAJA::statement::For< 0, RAJA::loop_exec,         // i
#endif
          RAJA::statement::Lambda< 0 >
        > // END i
      > // END j
    > // END k
  >; // END kernel

  /* clang-format on */
};
#endif

// Device Kernel settings:
//
// CudaKernel/HipKernel launches 256 threads total
// - In 2D, the launch configuration is 16 x 16
// - In 3D, the launch configuration is 8 x 8 x 4
constexpr int DEVICE_KERNEL_FIXED_SIZE = 256;
constexpr int TILE_SIZE_2D = 16;
constexpr int TILE_SIZE_X = 8;
constexpr int TILE_SIZE_Y = 8;
constexpr int TILE_SIZE_Z = 4;

//--------------------------------------------------------| CUDA_EXEC |---------
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__CUDACC__)

template <int BLOCK_SIZE>
struct nested_for_exec<CUDA_EXEC<BLOCK_SIZE, SYNCHRONOUS>>
{
  /* clang-format off */

  using loop2d_policy = RAJA::KernelPolicy<
    RAJA::statement::CudaKernelFixed< DEVICE_KERNEL_FIXED_SIZE,
      RAJA::statement::Tile<1, RAJA::tile_fixed< TILE_SIZE_2D >, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<0, RAJA::tile_fixed< TILE_SIZE_2D >, RAJA::cuda_block_x_loop,
          RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
            RAJA::statement::For<0, RAJA::cuda_thread_x_direct,
              RAJA::statement::Lambda<0>
            >
          >
        >
      >
    >
  >;

  using loop3d_policy = RAJA::KernelPolicy<
    RAJA::statement::CudaKernelFixed< DEVICE_KERNEL_FIXED_SIZE,
      RAJA::statement::Tile<2, RAJA::tile_fixed< TILE_SIZE_Z >, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<1, RAJA::tile_fixed< TILE_SIZE_Y >, RAJA::cuda_block_y_loop,
          RAJA::statement::Tile<0, RAJA::tile_fixed< TILE_SIZE_X >, RAJA::cuda_block_x_loop,
            RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
              RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
                RAJA::statement::For<0, RAJA::cuda_thread_x_direct,
                  RAJA::statement::Lambda<0>
                >
              >
            >
          >
        >
      >
    >
  >;

  /* clang-format on */
};

template <int BLOCK_SIZE>
struct nested_for_exec<CUDA_EXEC<BLOCK_SIZE, ASYNC>>
{
  /* clang-format off */

  using loop2d_policy = RAJA::KernelPolicy <
    RAJA::statement::CudaKernelFixedAsync< DEVICE_KERNEL_FIXED_SIZE,
      RAJA::statement::Tile<1, RAJA::tile_fixed< TILE_SIZE_2D >, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<0, RAJA::tile_fixed< TILE_SIZE_2D >, RAJA::cuda_block_x_loop,
          RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
            RAJA::statement::For<0, RAJA::cuda_thread_x_direct,
              RAJA::statement::Lambda<0>
            >
          >
        >
      >
    >
  >;

  using loop3d_policy = RAJA::KernelPolicy<
    RAJA::statement::CudaKernelFixedAsync< DEVICE_KERNEL_FIXED_SIZE,
      RAJA::statement::Tile<2, RAJA::tile_fixed< TILE_SIZE_Z >, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<1, RAJA::tile_fixed< TILE_SIZE_Y >, RAJA::cuda_block_y_loop,
          RAJA::statement::Tile<0, RAJA::tile_fixed< TILE_SIZE_X >, RAJA::cuda_block_x_loop,
            RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
              RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
                RAJA::statement::For<0, RAJA::cuda_thread_x_direct,
                  RAJA::statement::Lambda<0>
                >
              >
            >
          >
        >
      >
    >
  >;

  /* clang-format on */
};

#endif

//--------------------------------------------------------| HIP_EXEC |---------
#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__HIPCC__)

template <int BLOCK_SIZE>
struct nested_for_exec<HIP_EXEC<BLOCK_SIZE, SYNCHRONOUS>>
{
  /* clang-format off */

  using loop2d_policy = RAJA::KernelPolicy<
    RAJA::statement::HipKernelFixed< DEVICE_KERNEL_FIXED_SIZE,
      RAJA::statement::Tile<1, RAJA::tile_fixed< TILE_SIZE_2D >, RAJA::hip_block_y_loop,
        RAJA::statement::Tile<0, RAJA::tile_fixed< TILE_SIZE_2D >, RAJA::hip_block_x_loop,
          RAJA::statement::For<1, RAJA::hip_thread_y_direct,
            RAJA::statement::For<0, RAJA::hip_thread_x_direct,
              RAJA::statement::Lambda<0>
            >
          >
        >
      >
    >
  >;

  using loop3d_policy = RAJA::KernelPolicy<
    RAJA::statement::HipKernelFixed< DEVICE_KERNEL_FIXED_SIZE,
      RAJA::statement::Tile<2, RAJA::tile_fixed< TILE_SIZE_Z >, RAJA::hip_block_z_loop,
        RAJA::statement::Tile<1, RAJA::tile_fixed< TILE_SIZE_Y >, RAJA::hip_block_y_loop,
          RAJA::statement::Tile<0, RAJA::tile_fixed< TILE_SIZE_X >, RAJA::hip_block_x_loop,
            RAJA::statement::For<2, RAJA::hip_thread_z_direct,
              RAJA::statement::For<1, RAJA::hip_thread_y_direct,
                RAJA::statement::For<0, RAJA::hip_thread_x_direct,
                  RAJA::statement::Lambda<0>
                >
              >
            >
          >
        >
      >
    >
  >;

  /* clang-format on */
};

template <int BLOCK_SIZE>
struct nested_for_exec<HIP_EXEC<BLOCK_SIZE, ASYNC>>
{
  /* clang-format off */

  using loop2d_policy = RAJA::KernelPolicy <
    RAJA::statement::HipKernelFixedAsync< DEVICE_KERNEL_FIXED_SIZE,
      RAJA::statement::Tile<1, RAJA::tile_fixed< TILE_SIZE_2D >, RAJA::hip_block_y_loop,
        RAJA::statement::Tile<0, RAJA::tile_fixed< TILE_SIZE_2D >, RAJA::hip_block_x_loop,
          RAJA::statement::For<1, RAJA::hip_thread_y_direct,
            RAJA::statement::For<0, RAJA::hip_thread_x_direct,
              RAJA::statement::Lambda<0>
            >
          >
        >
      >
    >
  >;

  using loop3d_policy = RAJA::KernelPolicy<
    RAJA::statement::HipKernelFixedAsync< DEVICE_KERNEL_FIXED_SIZE,
      RAJA::statement::Tile<2, RAJA::tile_fixed< TILE_SIZE_Z >, RAJA::hip_block_z_loop,
        RAJA::statement::Tile<1, RAJA::tile_fixed< TILE_SIZE_Y >, RAJA::hip_block_y_loop,
          RAJA::statement::Tile<0, RAJA::tile_fixed< TILE_SIZE_X >, RAJA::hip_block_x_loop,
            RAJA::statement::For<2, RAJA::hip_thread_z_direct,
              RAJA::statement::For<1, RAJA::hip_thread_y_direct,
                RAJA::statement::For<0, RAJA::hip_thread_x_direct,
                  RAJA::statement::Lambda<0>
                >
              >
            >
          >
        >
      >
    >
  >;

  /* clang-format on */
};

#endif

} /* namespace internal */

} /* namespace axom */

#endif /* AXOM_CORE_NESTED_FOR_EXEC_HPP_ */
