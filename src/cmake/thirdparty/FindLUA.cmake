# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup Lua
#------------------------------------------------------------------------------
# find_package(LUA) can use the LUA_DIR env variable as a hint
# ------------------------------------------------------------------------------

# first Check for LUA_DIR
if (NOT EXISTS "${LUA_DIR}")
    message(FATAL_ERROR "Given LUA_DIR does not exist: ${LUA_DIR}")
endif()

if (NOT IS_DIRECTORY "${LUA_DIR}")
    message(FATAL_ERROR "Given LUA_DIR is not a directory: ${LUA_DIR}")
endif()

# Add LUA_DIR as env variable
set(ENV{LUA_DIR} ${LUA_DIR})

# Uncomment the following for more debug output in FindLUA
#set(LUA_Debug TRUE)

# HACK: Workaround for lua@5.4 and older versions of cmake (e.g. 3.16)
# which did not account for versions of lua beyond 5.3
string(FIND ${LUA_DIR} lua-5.4 _is_lua_5_4)
if(NOT ${_is_lua_5_4} EQUAL -1)
    find_package(Lua 5.4 EXACT REQUIRED) 
else()
    find_package(Lua REQUIRED) 
endif()

if(NOT LUA_FOUND)
    message(FATAL_ERROR "LUA_DIR is not a path to a valid Lua install")
endif()

message(STATUS "Lua Includes: ${LUA_INCLUDE_DIR}")
message(STATUS "Lua Libraries: ${LUA_LIBRARIES}")
