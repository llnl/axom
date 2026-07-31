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
#   1. build the wheel from src/python against the prebuilt Axom install
#      specified by AXOM_DIR or AXOM_INSTALL (find_package(axom));
#   2. install the wheel into a fresh uv venv;
#   3. verify the wheel installed conduit.pth for the same-build Conduit python module;
#   4. run the Sidre Python test suite with plain pytest.
#
# Intended for the gcc docker image, which is nanobind-enabled.

# Fail on the first error, including inside pipelines, and trace every command so
# a CI failure is readable from the log alone.
set -e
set -o pipefail
set -x

HOST_CONFIG="${HOST_CONFIG:-host-configs/docker/gcc@13.3.1.cmake}"

echo "~~~~ helpful info ~~~~"
echo "USER=$(id -u -n)"
echo "PWD=$(pwd)"
echo "HOST_CONFIG=${HOST_CONFIG}"
echo "~~~~~~~~~~~~~~~~~~~~~~"

absolute_path() {
    local path="$1"
    if [[ ! -e "${path}" ]]; then
        echo "ERROR: Path does not exist: ${path}" >&2
        return 1
    fi
    local dir
    local base
    dir=$(dirname "${path}")
    base=$(basename "${path}")
    printf "%s/%s" "$(cd "${dir}" && pwd -P)" "${base}"
}

if [[ ! -f "${HOST_CONFIG}" ]]; then
    echo "ERROR: Host-config not found: ${HOST_CONFIG}" >&2
    exit 1
fi
HOST_CONFIG_PATH=$(absolute_path "${HOST_CONFIG}")
echo "HOST_CONFIG_PATH=${HOST_CONFIG_PATH}"

# extract value from host-config line of the form `set(${name} ON CACHE BOOL "")`
# then capitalizes it and looks for true-like patterns
cmake_bool_from_file_is_on() {
    local name="$1"
    local value
    value=$(awk -v name="${name}" '
        $0 ~ "set\\(" name "[ \t\"]+" {
            line = $0
            sub("^[ \t]*set\\(" name "[ \t\"]+", "", line)
            sub("[ \t\"\\)].*$", "", line)
            print line
            exit
        }
    ' "${HOST_CONFIG_PATH}")
    value="${value^^}"
    [[ "${value}" == "ON" || "${value}" == "TRUE" || "${value}" == "YES" || "${value}" == "1" ]]
}

if [[ -n "${AXOM_INSTALL:-}" && -z "${AXOM_DIR:-}" ]]; then
    AXOM_DIR="${AXOM_INSTALL%/}/lib/cmake"
fi

if [[ -z "${AXOM_DIR:-}" || ! -f "${AXOM_DIR}/axom-config.cmake" ]]; then
    echo "ERROR: Axom CMake package not found." >&2
    echo "       Set AXOM_DIR to the directory containing axom-config.cmake," >&2
    echo "       or set AXOM_INSTALL to an Axom install prefix." >&2
    exit 1
fi
AXOM_DIR=$(absolute_path "${AXOM_DIR}")
echo "AXOM_DIR=${AXOM_DIR}"

AXOM_WHEEL_ENABLE_MPI=OFF
if cmake_bool_from_file_is_on ENABLE_MPI; then
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
uv build --wheel \
    -C cmake.args=-C \
    -C "cmake.args=${HOST_CONFIG_PATH}" \
    -C "cmake.define.AXOM_DIR=${AXOM_DIR}" \
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
