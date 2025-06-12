# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *

class Axomdevtools(BundlePackage):
    """This is a set of tools necessary for the developers of Axom"""

    version('fakeversion')

    maintainers = ['white238']

    depends_on("python")
    depends_on("doxygen")
    depends_on("cppcheck+rules+htmlreport")
    depends_on("graphviz")
    depends_on("py-sphinx")
    depends_on("py-shroud")
    depends_on("py-sphinxcontrib-jquery")
    depends_on("py-jsonschema")
    depends_on("llvm+clang@10.0.0")
