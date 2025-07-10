import os

from spack.package import *
from spack_repo.builtin.packages.redset.package import Redset as BuiltinRedset


class Redset(BuiltinRedset):

  depends_on("cxx", type="build")
