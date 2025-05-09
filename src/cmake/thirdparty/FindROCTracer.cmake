# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sets up target for roctracer and roctx libraries
#------------------------------------------------------------------------------

if(NOT ENABLE_HIP)
  message(FATAL_ERROR "Please set up hip before setting up ROC-tracer")
endif()

find_path(ROCTRACER_INCLUDE_DIR
          NAMES roctracer.h
          HINTS ${ROCM_PATH}/include/roctracer)

blt_find_libraries(
    FOUND_LIBS ROCTRACER_LIBRARIES
    NAMES roctracer64 roctx64
    REQUIRED FALSE
    PATHS ${ROCM_PATH}/lib)

find_package_handle_standard_args(ROCTRACER
  DEFAULT_MSG
  ROCTRACER_LIBRARIES
  ROCTRACER_INCLUDE_DIR)

mark_as_advanced(
    ROCTRACER_LIBRARIES
    ROCTRACER_INCLUDE_DIR)

blt_import_library(
  NAME roctracer
  INCLUDES ${ROCTRACER_INCLUDE_DIR}
  LIBRARIES  ${ROCTRACER_LIBRARIES}
  TREAT_INCLUDES_AS_SYSTEM ON
  EXPORTABLE ON
)

