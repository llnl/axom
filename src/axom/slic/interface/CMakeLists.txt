# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Generate slic's C and Fortran interfaces
#------------------------------------------------------------------------------

#
# Generate source
# slic/src/interface/CMakeLists.txt
#
set(YAML_INPUT_FILE slic_shroud.yaml)

add_shroud(
    YAML_INPUT_FILE      ${YAML_INPUT_FILE}
    C_FORTRAN_OUTPUT_DIR c_fortran
)
