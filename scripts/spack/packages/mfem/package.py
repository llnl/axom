import os

from spack.package import *
from spack_repo.builtin.packages.mfem.package import Mfem as BuiltinMfem

class Mfem(BuiltinMfem):

    ## mfem fails to ld hypre otherwise
    ### BEGIN AXOM PATCH
    depends_on("hypre~shared", when="+cuda~shared")
    ### END AXOM PATCH
