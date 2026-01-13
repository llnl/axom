// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file AxomMacros.hpp
 *
 * \brief Contains several useful macros for the axom project
 */

#ifndef AXOM_MACROS_HPP_
#define AXOM_MACROS_HPP_

#include "axom/config.hpp"
#include <cassert>  // for assert()

// Header for assert() in HIP device kernels
#if defined(AXOM_USE_HIP)
  #include <hip/hip_runtime.h>
#endif

// _guarding_macros_start
/*!
 * \def AXOM_USE_GPU
 *
 * \brief Convenience macro used for GPU-enabled checks
 *
 * \note AXOM_USE_CUDA is defined if Axom is built with CUDA.
 *       AXOM_USE_HIP is defined if Axom is built with HIP.
 */
#if defined(AXOM_USE_CUDA) || defined(AXOM_USE_HIP)
  #define AXOM_USE_GPU
#endif

/*!
 * \def AXOM_DEVICE_CODE
 *
 * \brief Convenience macro used for kernel code
 */
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  #define AXOM_DEVICE_CODE
#endif

/*!
 * \def AXOM_GPUCC
 *
 * \brief Convenience macro for compiling CUDA/HIP source files
 */
#if defined(__CUDACC__) || defined(__HIPCC__)
  #define AXOM_GPUCC
#endif
// _guarding_macros_end

// _decorating_macros_start
/*!
 * \def AXOM_DEVICE
 * \def AXOM_HOST_DEVICE
 *
 * \brief CUDA or HIP host/device macros for decorating functions/lambdas
 *
 * \note These will expand to the corresponding CUDA/HIP decorations when
 *  compiled with -DAXOM_USE_CUDA or -DAXOM_USE_HIP
 */
#if defined(__CUDACC__) || defined(__HIPCC__)
  #define AXOM_DEVICE __device__
  #define AXOM_HOST_DEVICE __host__ __device__
  #define AXOM_HOST __host__
#else
  #define AXOM_DEVICE
  #define AXOM_HOST_DEVICE
  #define AXOM_HOST
#endif

/*!
 * \def AXOM_LAMBDA
 *
 * \brief Convenience macro used for lambda capture by value.
 * \note When CUDA or HIP is used, the macro always expands to a host/device lambda.
 */
#if defined(AXOM_USE_CUDA) || defined(AXOM_USE_HIP)
  #define AXOM_LAMBDA [=] AXOM_HOST_DEVICE
  #define AXOM_DEVICE_LAMBDA [=] AXOM_DEVICE
  #define AXOM_HOST_LAMBDA [=] AXOM_HOST
#else
  #define AXOM_LAMBDA [=]
  #define AXOM_DEVICE_LAMBDA [=]
  #define AXOM_HOST_LAMBDA [=]
#endif
// _decorating_macros_end

/*
 * \def AXOM_STRINGIFY
 *
 * \brief Helper Macro to specify strings inside pragmas.
 */
#define AXOM_STRINGIFY(x) AXOM_DO_STRINGIFY(x)
#define AXOM_DO_STRINGIFY(x) #x

/*
 * \def AXOM_PRAGMA
 *
 * \brief Macro used to specify pragma directive
 */
#ifdef _WIN32
  #define AXOM_PRAGMA(x) __pragma(x)
#else
  #define AXOM_PRAGMA(x) _Pragma(AXOM_STRINGIFY(x))
#endif

/*
 * \def AXOM_SUPPRESS_HD_WARN
 *
 * \brief Macro used to silence __host__ __device__ compiler warnings
 *  when calling __host__ function from __host__ __device__ function.
 */
#if defined(__CUDACC__)
  #define AXOM_SUPPRESS_HD_WARN AXOM_PRAGMA(nv_exec_check_disable)
#else
  #define AXOM_SUPPRESS_HD_WARN
#endif

/*!
 * \def AXOM_CUDA_TEST
 *
 * \brief Convenience macro used for a gtest that uses cuda.
 */
#if defined(AXOM_USE_CUDA)
  #define AXOM_CUDA_TEST(X, Y)         \
    static void cuda_test_##X##Y();    \
    TEST(X, Y) { cuda_test_##X##Y(); } \
    static void cuda_test_##X##Y()
#else
  #define AXOM_CUDA_TEST(X, Y) TEST(X, Y)
#endif

/*
 * \def AXOM_SUPPRESS_UBSAN
 *
 * \brief Macro used to silence UndefinedBehaviorSanitizer errors
 *  when compiling and linking with -fsanitize=undefined
 */
#if defined(__has_attribute)
  #if __has_attribute(no_sanitize)
    #define AXOM_SUPPRESS_UBSAN __attribute__((no_sanitize("undefined")))
  #else
    #define AXOM_SUPPRESS_UBSAN
  #endif
#else
  #define AXOM_SUPPRESS_UBSAN
#endif

/*!
 *
 * \def AXOM_UNUSED_PARAM(x)
 * \brief Macro used to silence compiler warnings in methods with unused arguments.
 * \note The intent is to use this macro in the function signature. For example:
 * \code
 *
 *  void my_function(int x, int AXOM_UNUSED_PARAM(y))
 *  {
 *    // my implementation
 *  }
 *
 * \endcode
 */
#define AXOM_UNUSED_PARAM(x)

/*!
 *
 * \def AXOM_STATIC_ASSERT(cond)
 * \def AXOM_STATIC_ASSERT_MSG(cond, MSG)
 *
 * \brief This macro wraps the compile time static_assert functionality
 *  so you don't have to provide a message.
 */
#define AXOM_STATIC_ASSERT(cond) static_assert(cond, #cond)
#define AXOM_STATIC_ASSERT_MSG(cond, MSG) static_assert(cond, MSG)

/*!
 *
 * \def AXOM_UNUSED_VAR(x)
 * \brief Macro used to silence compiler warnings about variables
 *        that are defined but not used.
 * \note The intent is to use this macro for variables that are only used
 *       within compiler defines (e.g. in debug assertions). For example:
 * \code
 *
 *  double myVar = ...
 *  AXOM_UNUSED_VAR(myVar);
 *
 *  // code emits the following warning in release builds
 *  // if extra warnings are enabled and this macro is not called
 *  // warning: unused variable 'myVar' [-Wunused-variable]
 *
 *  SLIC_ASSERT(myVar > 0)
 *
 * \endcode
 */
#define AXOM_UNUSED_VAR(_x) static_cast<void>(_x)

/*!
 * \def AXOM_DEBUG_PARAM(x)
 * \brief Macro used to silence compiler warnings about parameters
 *        that are only used when AXOM_DEBUG is defined
 * \note Default values are ok
 * \code
 *
 *  void my_function(int x, int AXOM_DEBUG_PARAM(y))
 *  {
 *    // my implementation
 *    SLIC_ASSERT(y > 0)
 *  }
 *
 * \endcode
 */
#ifdef AXOM_DEBUG
  #define AXOM_DEBUG_PARAM(_x) _x
#else
  #define AXOM_DEBUG_PARAM(_x)
#endif

/*!
 * \def AXOM_MAYBE_UNUSED
 * \brief Macro that is used to silence compiler warnings about variables that
 *        might be unused in some build configurations.
 *
 * \note This macro can be used for the same purpose as AXOM_UNUSED_VAR and AXOM_DEBUG_PARAM but it
 *       is applied before the variable type declaration.
 * \code
 *
 *  AXOM_MAYBE_UNUSED double myVar = ...
 *  SLIC_ASSERT(myVar > 0)
 *
 * \endcode
 */
#if __cplusplus >= 201703L
  // C++17 and later.
  #define AXOM_MAYBE_UNUSED [[maybe_unused]]
#else
  #define AXOM_MAYBE_UNUSED
#endif

/*!
 * \def AXOM_LIKELY
 * \brief Macro that is used to annotate for compiler optimizations a
 *  conditional that is likely to be true.
 *
 * \note This macro is placed before the conditional in an if-statement.
 * \warning Use with caution, as unwarranted usage may result in pessimistic
 *  optimizations.
 * \code
 *
 *  bool success = UnlikelyToFail();
 *  if AXOM_LIKELY(success) { ... }
 *
 * \endcode
 */
#if __cplusplus >= 202002L
  // C++20 and later
  #define AXOM_LIKELY(cond) (cond) [[likely]]
#elif defined(__GNUC__)
  // GCC/Clang compilers have __builtin_expect
  #define AXOM_LIKELY(cond) (__builtin_expect(!!(cond), 1))
#else
  #define AXOM_LIKELY(cond) (cond)
#endif

/*!
 * \def AXOM_UNLIKELY
 * \brief Macro that is used to annotate for compiler optimizations a
 *  conditional that is likely to be false.
 *
 * \note This macro is placed before the conditional in an if-statement.
 * \warning Use with caution, as unwarranted usage may result in pessimistic
 *  optimizations.
 * \code
 *
 *  bool error = UnlikelyToFail();
 *  if AXOM_UNLIKELY(error) { ... }
 *
 * \endcode
 */
#if __cplusplus >= 202002L
  // C++20 and later
  #define AXOM_UNLIKELY(cond) (cond) [[unlikely]]
#elif defined(__GNUC__)
  // GCC/Clang compilers have __builtin_expect
  #define AXOM_UNLIKELY(cond) (__builtin_expect(!!(cond), 0))
#else
  #define AXOM_UNLIKELY(cond) (cond)
#endif

/*!
 * \def DISABLE_DEFAULT_CTOR(className)
 * \brief Macro to disable default constructor for the given class.
 * \note This macro should only be used within the private section of a class,
 *  as indicated in the example below.
 *
 * \code
 *
 *   class Foo
 *   {
 *   public:
 *
 *       // Public methods here
 *
 *   private:
 *      DISABLE_DEFAULT_CTOR(Foo);
 *   };
 *
 * \endcode
 */
#define DISABLE_DEFAULT_CTOR(className) className() = delete

/*!
 * \def DISABLE_COPY_AND_ASSIGNMENT(className)
 * \brief Macro to disable copy and assignment operations for the given class.
 * \note This macro should only be used within the private section of a class,
 *  as indicated in the example below.
 *
 * \code
 *
 *   class Foo
 *   {
 *   public:
 *      Foo();
 *      ~Foo();
 *
 *       // Other methods here
 *
 *   private:
 *      DISABLE_COPY_AND_ASSIGNMENT(Foo);
 *   };
 *
 * \endcode
 */
#define DISABLE_COPY_AND_ASSIGNMENT(className) \
  className(const className&) = delete;        \
  className& operator=(const className&) = delete

/*!
 * \def DISABLE_MOVE_AND_ASSIGNMENT(className)
 * \brief Macro to disable move constructor and move assignment operations for
 * the given class.
 * \note This macro should only be used within the private section of a class,
 *  as indicated in the example below.
 *
 * \code
 *
 *   class Foo
 *   {
 *   public:
 *      Foo();
 *      ~Foo();
 *
 *       // Other methods here
 *
 *   private:
 *      DISABLE_MOVE_AND_ASSIGNMENT(Foo);
 *   };
 *
 * \endcode
 */
#define DISABLE_MOVE_AND_ASSIGNMENT(className) \
  className(className&&) = delete;             \
  className& operator=(className&&) = delete

/*!
 * \def AXOM_TYPED_TEST(CaseName, TestName)
 * \brief Minor tweak of gtest's TYPED_TEST macro to work with device code
 * \note Can be used in test files after including `gtest/gtest.h`
 */
// Specifically, we expose the `TestBody()` method as public instead of private
#define AXOM_TYPED_TEST(CaseName, TestName)                                                       \
  static_assert(sizeof(GTEST_STRINGIFY_(TestName)) > 1, "test-name must not be empty");           \
  template <typename gtest_TypeParam_>                                                            \
  class GTEST_TEST_CLASS_NAME_(CaseName, TestName) : public CaseName<gtest_TypeParam_>            \
  {                                                                                               \
  public:                                                                                         \
    typedef CaseName<gtest_TypeParam_> TestFixture;                                               \
    typedef gtest_TypeParam_ TypeParam;                                                           \
    void TestBody() override;                                                                     \
  };                                                                                              \
  GTEST_INTERNAL_ATTRIBUTE_MAYBE_UNUSED static bool gtest_##CaseName##_##TestName##_registered_ = \
    ::testing::internal::TypeParameterizedTest<                                                   \
      CaseName,                                                                                   \
      ::testing::internal::TemplateSel<GTEST_TEST_CLASS_NAME_(CaseName, TestName)>,               \
      GTEST_TYPE_PARAMS_(                                                                         \
        CaseName)>::Register("",                                                                  \
                             ::testing::internal::CodeLocation(__FILE__, __LINE__),               \
                             GTEST_STRINGIFY_(CaseName),                                          \
                             GTEST_STRINGIFY_(TestName),                                          \
                             0,                                                                   \
                             ::testing::internal::GenerateNames<GTEST_NAME_GENERATOR_(CaseName),  \
                                                                GTEST_TYPE_PARAMS_(CaseName)>()); \
  template <typename gtest_TypeParam_>                                                            \
  void GTEST_TEST_CLASS_NAME_(CaseName, TestName)<gtest_TypeParam_>::TestBody()

#endif  // AXOM_MACROS_HPP_
