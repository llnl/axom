import os

from spack.package import *
from spack.pkg.builtin.py_nanobind import PyNanobind as BuiltinNanobind

class PyNanobind(BuiltinNanobind):

    version(
        "2.7.0", tag="v2.7.0", commit="44ad9a9e5729abda24ef8dc9d76233d801e651e9", submodules=True
    )
    version(
        "2.6.1", tag="v2.6.1", commit="9b3afa9dbdc23641daf26fadef7743e7127ff92f", submodules=True
    )