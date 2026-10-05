#!/bin/bash
##############################################################################
# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
##############################################################################

# Build a wheel against an installed Axom, install it in a fresh uv venv,
# check axom-conduit.pth, and run the Sidre Python tests with pytest.
# The test avoids PYTHONPATH changes and the run_python_with_axom.sh wrapper.
# Run in the nanobind-enabled GCC Docker image.
#
# Required environment:
#   HOST_CONFIG            host-config the Axom install was configured with; it supplies the
#                          compilers, which the installed axom-config.cmake does not record
#   AXOM_DIR/AXOM_INSTALL  Axom install prefix

# Stop on errors, including pipeline failures, and log each command.
set -e
set -o pipefail
set -x

# Require the Axom install's host-config to select matching compilers.
if [[ -z "${HOST_CONFIG:-}" ]]; then
    echo "ERROR: HOST_CONFIG is not set." >&2
    echo "       Set it to the host-config used to configure the Axom install," >&2
    echo "       e.g. HOST_CONFIG=host-configs/docker/gcc@13.3.1.cmake" >&2
    exit 1
fi

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

# Read the first set() value for a variable from a generated CMake file.
cmake_value_from_file() {
    local file="$1"
    local name="$2"
    awk -v name="${name}" '
        $0 ~ "set\\(" name "[ \t\"]+" {
            line = $0
            sub("^[ \t]*set\\(" name "[ \t\"]+", "", line)
            sub("[ \t\"\\)].*$", "", line)
            print line
            exit
        }
    ' "${file}"
}

# Check a generated CMake value for ON, TRUE, YES, or 1.
cmake_bool_from_file_is_on() {
    local value
    value=$(cmake_value_from_file "$1" "$2")
    value="${value^^}"
    [[ "${value}" == "ON" || "${value}" == "TRUE" || "${value}" == "YES" || "${value}" == "1" ]]
}

AXOM_DIR="${AXOM_DIR:-${AXOM_INSTALL:-}}"
if [[ -z "${AXOM_DIR}" || ! -f "${AXOM_DIR%/}/lib/cmake/axom-config.cmake" ]]; then
    echo "ERROR: Axom install not found." >&2
    echo "       Set AXOM_DIR (or AXOM_INSTALL) to an Axom install prefix," >&2
    echo "       i.e. the directory whose lib/cmake holds axom-config.cmake." >&2
    exit 1
fi
AXOM_DIR=$(absolute_path "${AXOM_DIR}")
echo "AXOM_DIR=${AXOM_DIR}"
AXOM_CONFIG="${AXOM_DIR}/lib/cmake/axom-config.cmake"

# Conduit includes a CPython extension, so the wheel must use the same Python
# minor version as the Axom/Conduit installation.
AXOM_PYTHON_EXECUTABLE=$(cmake_value_from_file "${AXOM_CONFIG}" AXOM_PYTHON_EXECUTABLE)
if [[ -z "${AXOM_PYTHON_EXECUTABLE}" || ! -x "${AXOM_PYTHON_EXECUTABLE}" ]]; then
    echo "ERROR: ${AXOM_CONFIG} records an invalid AXOM_PYTHON_EXECUTABLE:" >&2
    echo "       '${AXOM_PYTHON_EXECUTABLE}'" >&2
    exit 1
fi
echo "AXOM_PYTHON_EXECUTABLE=${AXOM_PYTHON_EXECUTABLE}"
"${AXOM_PYTHON_EXECUTABLE}" --version

# Read MPI support from the installed Axom configuration.
AXOM_WHEEL_ENABLE_MPI=OFF
if cmake_bool_from_file_is_on "${AXOM_CONFIG}" AXOM_USE_MPI; then
    AXOM_WHEEL_ENABLE_MPI=ON
fi
echo "AXOM_WHEEL_ENABLE_MPI=${AXOM_WHEEL_ENABLE_MPI}"

echo "~~~~~~ ENSURE uv IS AVAILABLE ~~~~~~"
# Use uv from PATH or install the pinned version in a separate directory.
# --target avoids modifying a system-managed Python environment.
AXOM_UV_VERSION="${AXOM_UV_VERSION:-0.12.21}"
if ! command -v uv >/dev/null 2>&1; then
    UV_BOOTSTRAP_DIR="${RUNNER_TEMP:-${TMPDIR:-/tmp}}/axom-uv-${AXOM_UV_VERSION}"
    python3 -m pip install --disable-pip-version-check --target "${UV_BOOTSTRAP_DIR}" "uv==${AXOM_UV_VERSION}"
    export PATH="${UV_BOOTSTRAP_DIR}/bin:${PATH}"
fi
uv --version

echo "~~~~~~ BUILD THE THIN WHEEL FROM src/python ~~~~~~"
# axom-config.cmake records the Conduit install to use.
rm -rf dist
uv build --wheel \
    --python "${AXOM_PYTHON_EXECUTABLE}" \
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
# Use Axom's Python for the test venv so Conduit's compiled module has the same ABI.
VENV_DIR=/tmp/axom-wheel-venv
rm -rf "${VENV_DIR}"
uv venv --python "${AXOM_PYTHON_EXECUTABLE}" "${VENV_DIR}"
VENV_PY="${VENV_DIR}/bin/python"
AXOM_WHEEL_EXTRAS="test"
if [[ "${AXOM_WHEEL_ENABLE_MPI}" == "ON" ]]; then
    AXOM_WHEEL_EXTRAS="test,mpi"
fi
uv pip install --python "${VENV_PY}" "${AXOM_WHEEL}[${AXOM_WHEEL_EXTRAS}]"

echo "~~~~~~ VERIFY WHEEL-INSTALLED CONDUIT .pth ~~~~~~"
PLATLIB=$("${VENV_PY}" -c 'import sysconfig; print(sysconfig.get_paths()["platlib"])')
CONDUIT_PTH="${PLATLIB}/axom-conduit.pth"
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
