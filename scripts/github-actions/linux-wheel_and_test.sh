#!/bin/bash
##############################################################################
# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
##############################################################################

# Build the thin, pip/uv-installable Axom wheel and exercise it end to end,
# without using the run_python_with_axom.sh wrapper or updating the PYTHONPATH:
#
#   1. configure + build + install Axom (with Python bindings) from a docker t-config;
#   2. build the wheel from src/python against that install (find_package(axom));
#   3. install the wheel into a fresh uv venv;
#   4. verify the wheel installed conduit.pth for the same-build Conduit python module;
#   5. run the Sidre Python test suite with plain pytest.
#
# Intended for the gcc docker image, which is nanobind-enabled.

# Fail on the first error, including inside pipelines, and trace every command so
# a CI failure is readable from the log alone.
set -e
set -o pipefail
set -x

HOST_CONFIG="${HOST_CONFIG:-gcc@13.3.1.cmake}"
BUILD_TYPE="${BUILD_TYPE:-Debug}"
BUILD_DIR="${BUILD_DIR:-builddir_wheel}"

echo "~~~~ helpful info ~~~~"
echo "USER=$(id -u -n)"
echo "PWD=$(pwd)"
echo "HOST_CONFIG=${HOST_CONFIG}"
echo "BUILD_TYPE=${BUILD_TYPE}"
echo "~~~~~~~~~~~~~~~~~~~~~~"

NUM_BUILD_PROCS=$(python3 -c 'import os; print(max(2, os.cpu_count() * 8 // 10))')

echo "~~~~~~ CONFIGURE + BUILD + INSTALL AXOM (+python) ~~~~~~"
python3 ./config-build.py \
    -bp "${BUILD_DIR}" \
    -hc "./host-configs/docker/${HOST_CONFIG}" \
    -bt "${BUILD_TYPE}"
cmake --build "${BUILD_DIR}" -j "${NUM_BUILD_PROCS}"
cmake --install "${BUILD_DIR}"

# Resolve the Axom install prefix from the CMake cache
CACHE="${BUILD_DIR}/CMakeCache.txt"
AXOM_INSTALL=$(awk -F= '/^CMAKE_INSTALL_PREFIX:[A-Z]*=/{print $2}' "${CACHE}")
echo "AXOM_INSTALL=${AXOM_INSTALL}"

if [[ -z "${AXOM_INSTALL}" || ! -d "${AXOM_INSTALL}" ]]; then
    echo "ERROR: Axom install prefix not found (${AXOM_INSTALL})."
    exit 1
fi

cache_value() {
    awk -F= -v key="$1" '$1 ~ "^" key ":[^=]*$" {print $2; exit}' "${CACHE}"
}

cache_bool_is_on() {
    local value
    value=$(cache_value "$1")
    case "${value^^}" in
        ON|TRUE|YES|1)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

require_cmake_define_from_cache() {
    local name="$1"
    local value
    value=$(cache_value "${name}")
    if [[ -z "${value}" ]]; then
        echo "ERROR: Required CMake cache entry ${name} was not found in ${CACHE}."
        exit 1
    fi
    WHEEL_BUILD_ARGS+=("-C" "cmake.define.${name}=${value}")
}

AXOM_WHEEL_ENABLE_MPI=OFF
if cache_bool_is_on ENABLE_MPI || cache_bool_is_on AXOM_ENABLE_MPI; then
    AXOM_WHEEL_ENABLE_MPI=ON
fi
echo "AXOM_WHEEL_ENABLE_MPI=${AXOM_WHEEL_ENABLE_MPI}"

echo "~~~~~~ ENSURE uv IS AVAILABLE ~~~~~~"
if ! command -v uv >/dev/null 2>&1; then
    python3 -m pip install --user uv
    export PATH="${HOME}/.local/bin:${PATH}"
fi
uv --version

echo "~~~~~~ BUILD THE THIN WHEEL FROM src/python ~~~~~~"
# Point find_package at the install with AXOM_DIR
# Conduit resolves transitively from axom's config, which records its Conduit prefix
rm -rf dist
WHEEL_BUILD_ARGS=(
    "-C" "cmake.define.AXOM_DIR=${AXOM_INSTALL}/lib/cmake"
)
require_cmake_define_from_cache CMAKE_C_COMPILER
require_cmake_define_from_cache CMAKE_CXX_COMPILER
if [[ "${AXOM_WHEEL_ENABLE_MPI}" == "ON" ]]; then
    require_cmake_define_from_cache MPI_C_COMPILER
    require_cmake_define_from_cache MPI_CXX_COMPILER
fi
uv build --wheel "${WHEEL_BUILD_ARGS[@]}" --out-dir dist src/python
ls -l dist
AXOM_WHEEL=$(find dist -maxdepth 1 -name 'axom-*.whl' -print -quit)
if [[ -z "${AXOM_WHEEL}" ]]; then
    echo "ERROR: Axom wheel not found in dist/."
    exit 1
fi

echo "~~~~~~ FRESH VENV + INSTALL THE WHEEL ~~~~~~"
# Pin the interpreter that built the wheel, so the venv cannot pick a different one.
VENV_DIR=/tmp/axom-wheel-venv
rm -rf "${VENV_DIR}"
uv venv --python "$(command -v python3)" "${VENV_DIR}"
VENV_PY="${VENV_DIR}/bin/python"
AXOM_WHEEL_EXTRAS="test"
if [[ "${AXOM_WHEEL_ENABLE_MPI}" == "ON" ]]; then
    AXOM_WHEEL_EXTRAS="test,mpi"
fi
uv pip install --python "${VENV_PY}" "${AXOM_WHEEL}[${AXOM_WHEEL_EXTRAS}]"

echo "~~~~~~ VERIFY WHEEL-INSTALLED CONDUIT .pth ~~~~~~"
PLATLIB=$("${VENV_PY}" -c 'import sysconfig; print(sysconfig.get_paths()["platlib"])')
CONDUIT_PTH="${PLATLIB}/conduit.pth"
if [[ ! -f "${CONDUIT_PTH}" ]]; then
    echo "ERROR: Expected wheel to install ${CONDUIT_PTH}."
    echo "       The wheel should expose the same-build Conduit python module without a manual PYTHONPATH update."
    exit 1
fi
CONDUIT_PY_DIR=$(sed -n '1p' "${CONDUIT_PTH}")
if [[ -z "${CONDUIT_PY_DIR}" || ! -d "${CONDUIT_PY_DIR}" ]]; then
    echo "ERROR: ${CONDUIT_PTH} points to missing Conduit python module directory '${CONDUIT_PY_DIR}'."
    exit 1
fi
echo "verified ${CONDUIT_PTH} -> ${CONDUIT_PY_DIR}"

echo "~~~~~~ IMPORT SMOKE TEST ~~~~~~"
"${VENV_PY}" -c \
    "import axom, axom.sidre, conduit, numpy; print('axom', axom.__version__); print('axom.sidre', axom.sidre.__version__)"
if [[ "${AXOM_WHEEL_ENABLE_MPI}" == "ON" ]]; then
    "${VENV_PY}" -c "import mpi4py, axom.sidre as sidre; assert sidre.AXOM_ENABLE_MPI"
else
    "${VENV_PY}" -c "import axom.sidre as sidre; assert not sidre.AXOM_ENABLE_MPI"
fi

echo "~~~~~~ RUN THE SIDRE PYTHON SUITE VIA PLAIN pytest ~~~~~~"
# Axom's Python tests are named *_Py.py, which pytest's default python_files patterns do not match
TEST_DIR="$(pwd)/src/axom/sidre/tests"
SCRATCH="$(mktemp -d)"
pushd "${SCRATCH}" > /dev/null
"${VENV_PY}" -m pytest -s -p no:cacheprovider \
    -o python_files='*_Py.py' \
    "${TEST_DIR}"
popd > /dev/null
