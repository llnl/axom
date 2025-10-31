# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# Axom version information that go into the generated axom/config.hpp
# NOTE: if we are building from a Git repo these values will be autodetected,
#       otherwise, the hard-coded values will go in to the config.
#------------------------------------------------------------------------------
set(AXOM_VERSION_MAJOR 0)
set(AXOM_VERSION_MINOR 12)
set(AXOM_VERSION_PATCH 0)
string(CONCAT AXOM_VERSION_FULL
    "v${AXOM_VERSION_MAJOR}"
    ".${AXOM_VERSION_MINOR}"
    ".${AXOM_VERSION_PATCH}" )

#------------------------------------------------------------------------------
# extract_version_components( <tag>, <prefix> )
#
# Helper function to extract the version components from a Git tag. The
# extracted version components are stored in the following variables:
#   <prefix>_VERSION_MAJOR
#   <prefix>_VERSION_MINOR
#   <prefix>_VERSION_PATCH
#   <prefix>_VERSION_FULL
#
# Note: This function assumes tags of the form "vMAJOR.MINOR.PATCH"
#------------------------------------------------------------------------------
function(axom_extract_version tag_string prefix)

  string(REGEX MATCH "v([0-9]+)\\.([0-9]+)\\.([0-9]+)" tmp "${tag_string}" )

  if ( CMAKE_MATCH_0 )
    set(full ${CMAKE_MATCH_0})
    set(major ${CMAKE_MATCH_1})
    set(minor ${CMAKE_MATCH_2})
    set(patch ${CMAKE_MATCH_3})

    set(${prefix}_VERSION_MAJOR ${major} PARENT_SCOPE)
    set(${prefix}_VERSION_MINOR ${minor} PARENT_SCOPE)
    set(${prefix}_VERSION_PATCH ${patch} PARENT_SCOPE)
    set(${prefix}_VERSION_FULL ${full} PARENT_SCOPE)
  endif()

endfunction()


if (Git_FOUND)

  ## check to see if we are building from a Git repo or an exported tarball
  blt_is_git_repo( OUTPUT_STATE is_git_repo )

  if ( ${is_git_repo} )

    ## get latest tag from main
    blt_git_tag( OUTPUT_TAG axom_tag RETURN_CODE rc ON_BRANCH main )
    if ( {rc} EQUAL 0)
      axom_extract_version( "${axom_tag}" AXOM )
    endif()

    blt_git_hashcode( HASHCODE sha1 RETURN_CODE rc )
    if ( NOT ${rc} EQUAL 0 )
      message( FATAL_ERROR "blt_git_hashcode failed!" )
    endif()

    set(AXOM_GIT_SHA ${sha1})

  endif()

endif()
