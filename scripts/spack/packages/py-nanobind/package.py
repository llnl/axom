import os

from spack.package import *
from spack_repo.builtin.packages.py_nanobind.package import PyNanobind as BuiltinPyNanobind

class PyNanobind(BuiltinPyNanobind):

    # Workaround for "gmake: *** internal error: invalid --jobserver-auth string":
    # https://github.com/spack/spack-packages/issues/5106#issuecomment-4679593276
    depends_on("gmake", type="build")
