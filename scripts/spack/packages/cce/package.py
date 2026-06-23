# Copyright 2013-2025 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
from spack.package import *
from spack_repo.builtin.packages.cce.package import Cce as BuiltinCce

class Cce(BuiltinCce):

    # Add C++20 standard
    def _standard_flag(self, *, language, standard):
        flags = {
            "cxx": {"11": "-std=c++11", "14": "-std=c++14", "17": "-std=c++17", "20": "-std=c++20"},
            "c": {"99": "-std=c99", "11": "-std=c11"},
        }
        return flags[language][standard]
