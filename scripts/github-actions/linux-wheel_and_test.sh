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

set -e
set -o pipefail

function or_die () {
    "$@"
    local status=$?
    if [[ $status != 0 ]]; then
        echo "ERROR $status command: $*"
        exit $status
    fi
}

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
or_die python3 ./config-build.py \
    -bp "${BUILD_DIR}" \
    -hc "./host-configs/docker/${HOST_CONFIG}" \
    -bt "${BUILD_TYPE}"
or_die cmake --build "${BUILD_DIR}" -j "${NUM_BUILD_PROCS}"
or_die cmake --install "${BUILD_DIR}"

# Resolve the Axom install prefix from the CMake cache
CACHE="${BUILD_DIR}/CMakeCache.txt"
AXOM_INSTALL=$(awk -F= '/^CMAKE_INSTALL_PREFIX:[A-Z]*=/{print $2}' "${CACHE}")
echo "AXOM_INSTALL=${AXOM_INSTALL}"

if [[ -z "${AXOM_INSTALL}" || ! -d "${AXOM_INSTALL}" ]]; then
    echo "ERROR: Axom install prefix not found (${AXOM_INSTALL})."
    exit 1
fi

echo "~~~~~~ ENSURE uv IS AVAILABLE ~~~~~~"
if ! command -v uv >/dev/null 2>&1; then
    or_die python3 -m pip install --user uv
    export PATH="${HOME}/.local/bin:${PATH}"
fi
uv --version

echo "~~~~~~ BUILD THE THIN WHEEL FROM src/python ~~~~~~"
# Point find_package at the install with AXOM_DIR
# Conduit resolves transitively from axom's config, which records its Conduit prefix
rm -rf dist
or_die uv build --wheel \
    -C cmake.define.AXOM_DIR="${AXOM_INSTALL}/lib/cmake" \
    --out-dir dist \
    src/python
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
or_die uv venv --python "$(command -v python3)" "${VENV_DIR}"
VENV_PY="${VENV_DIR}/bin/python"
or_die uv pip install --python "${VENV_PY}" "${AXOM_WHEEL}[test]"

echo "~~~~~~ VERIFY WHEEL-INSTALLED CONDUIT .pth ~~~~~~"
PURELIB=$("${VENV_PY}" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')
CONDUIT_PTH="${PURELIB}/conduit.pth"
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
or_die "${VENV_PY}" -c \
    "import axom, axom.sidre, conduit, numpy; print('axom', axom.__version__); print('axom.sidre', axom.sidre.__version__)"

echo "~~~~~~ RUN THE SIDRE PYTHON SUITE VIA PLAIN pytest ~~~~~~"
# Axom's Python tests are named *_Py.py, which pytest's default python_files patterns do not match
TEST_DIR="$(pwd)/src/axom/sidre/tests"
SCRATCH="$(mktemp -d)"
pushd "${SCRATCH}" > /dev/null
or_die "${VENV_PY}" -m pytest -s -p no:cacheprovider \
    -o python_files='*_Py.py' \
    "${TEST_DIR}"
popd > /dev/null
