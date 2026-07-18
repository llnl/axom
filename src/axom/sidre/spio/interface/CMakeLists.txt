# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Generate spio's C and Fortran interface 
#------------------------------------------------------------------------------

set(YAML_INPUT_FILE spio_shroud.yaml)

add_shroud(
    YAML_INPUT_FILE      ${YAML_INPUT_FILE}
    YAML_OUTPUT_DIR      yaml
    C_FORTRAN_OUTPUT_DIR c_fortran
)
