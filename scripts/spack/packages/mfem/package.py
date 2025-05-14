import os

from spack.package import *
from spack.pkg.builtin.mfem import Mfem as BuiltinMfem


class Mfem(BuiltinMfem):

    version(
        "4.8.0",
        sha256="49bd2a076b0d87863092cb55f8524b5292d9afb2e48c19f80222ada367819016",
        url="https://bit.ly/mfem-4-8",
        extension="tar.gz",
    )

