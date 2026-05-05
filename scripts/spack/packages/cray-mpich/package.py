# Copyright 2013-2025 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
from spack.package import *
from spack_repo.builtin.packages.cray_mpich.package import CrayMpich as BuiltinCrayMpich

class CrayMpich(BuiltinCrayMpich):

    version('8.1.29')

    variant("slurm", default=True, description="Added to get MPIEXEC_NUMPROC_FLAG right")
