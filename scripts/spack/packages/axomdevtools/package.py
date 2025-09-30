# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.build_systems.bundle import BundlePackage

class Axomdevtools(BundlePackage):
    """This is a set of tools necessary for the developers of Axom"""

    version('fakeversion')

    maintainers = ['white238']

    depends_on("c")
    depends_on("cxx")
    depends_on("fortran")
    depends_on("mpi")

    # 3.13 is bugfix/stable release
    depends_on("python@3.13")
    depends_on("doxygen")
    depends_on("cppcheck+rules")
    depends_on("graphviz")
    depends_on("py-sphinx")
    depends_on("py-shroud")
    depends_on("py-sphinxcontrib-jquery")

    # 4.18 builds py-rpds-py, which then needs rust...
    depends_on("py-jsonschema@4.17")
    depends_on("py-nanobind@2.7.0")
    depends_on("py-pytest")
    depends_on("py-numpy")

    depends_on("llvm+clang@19")
