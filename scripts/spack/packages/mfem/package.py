import os

from spack.package import *
from spack_repo.builtin.packages.mfem.package import Mfem as BuiltinMfem

class Mfem(BuiltinMfem):

    version("4.9.0", 
            url="https://bit.ly/mfem-4-9",
            sha256="6904974c8d5a6bcd127419c7b7adff873170d397ed2f0bccdf438e940e713af2",
            extension="tar.gz")

    ## mfem fails to ld hypre otherwise
    ### BEGIN AXOM PATCH
    depends_on("hypre~shared", when="+cuda~shared")
    ### END AXOM PATCH
