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

    ## Hack (KW,14May2024): We needed to add the `hipblas` library to the hip_libs
    ## I've copied the entire `get_make_config_options` function
    ## to inject the appropriate call below.
    ## We should remove this entire function once we can update 
    ## to a spack package that contains this change.
    ## It looks like spack@develop already has a fix for this
    def get_make_config_options(self, spec, prefix):
        def yes_no(varstr):
            return "YES" if varstr in self.spec else "NO"

        xcompiler = "" if "~cuda" in spec else "-Xcompiler="

        # We need to add rpaths explicitly to allow proper export of link flags
        # from within MFEM. We use the following two functions to do that.
        ld_flags_from_library_list = self.ld_flags_from_library_list
        ld_flags_from_dirs = self.ld_flags_from_dirs

        def find_optional_library(name, prefix):
            for shared in [True, False]:
                for path in ["lib64", "lib"]:
                    lib = find_libraries(
                        name, join_path(prefix, path), shared=shared, recursive=False
                    )
                    if lib:
                        return lib
            return LibraryList([])

        # Determine how to run MPI tests, e.g. when using '--test=root', when
        # Spack is run inside a batch system job.
        mfem_mpiexec = "mpirun"
        mfem_mpiexec_np = "-np"
        if "SLURM_JOBID" in os.environ:
            mfem_mpiexec = "srun"
            mfem_mpiexec_np = "-n"
        elif "LSB_JOBID" in os.environ:
            if "LLNL_COMPUTE_NODES" in os.environ:
                mfem_mpiexec = "lrun"
                mfem_mpiexec_np = "-n"
            else:
                mfem_mpiexec = "jsrun"
                mfem_mpiexec_np = "-p"
        elif "FLUX_EXEC_PATH" in os.environ:
            mfem_mpiexec = "flux run"
            mfem_mpiexec_np = "-n"
        elif "PBS_JOBID" in os.environ:
            mfem_mpiexec = "mpiexec"
            mfem_mpiexec_np = "-n"

        metis5_str = "NO"
        if ("+metis" in spec) and spec["metis"].satisfies("@5:"):
            metis5_str = "YES"

        zlib_var = "MFEM_USE_ZLIB" if (spec.satisfies("@4.1.0:")) else "MFEM_USE_GZSTREAM"

        options = [
            "PREFIX=%s" % prefix,
            "MFEM_USE_MEMALLOC=YES",
            "MFEM_DEBUG=%s" % yes_no("+debug"),
            # NOTE: env["CXX"] is the spack c++ compiler wrapper. The real
            # compiler is defined by env["SPACK_CXX"].
            "CXX=%s" % env["CXX"],
            "MFEM_USE_LIBUNWIND=%s" % yes_no("+libunwind"),
            "%s=%s" % (zlib_var, yes_no("+zlib")),
            "MFEM_USE_METIS=%s" % yes_no("+metis"),
            "MFEM_USE_METIS_5=%s" % metis5_str,
            "MFEM_THREAD_SAFE=%s" % yes_no("+threadsafe"),
            "MFEM_USE_MPI=%s" % yes_no("+mpi"),
            "MFEM_USE_LAPACK=%s" % yes_no("+lapack"),
            "MFEM_USE_SUPERLU=%s" % yes_no("+superlu-dist"),
            "MFEM_USE_STRUMPACK=%s" % yes_no("+strumpack"),
            "MFEM_USE_SUITESPARSE=%s" % yes_no("+suite-sparse"),
            "MFEM_USE_SUNDIALS=%s" % yes_no("+sundials"),
            "MFEM_USE_PETSC=%s" % yes_no("+petsc"),
            "MFEM_USE_SLEPC=%s" % yes_no("+slepc"),
            "MFEM_USE_PUMI=%s" % yes_no("+pumi"),
            "MFEM_USE_GSLIB=%s" % yes_no("+gslib"),
            "MFEM_USE_NETCDF=%s" % yes_no("+netcdf"),
            "MFEM_USE_MPFR=%s" % yes_no("+mpfr"),
            "MFEM_USE_GNUTLS=%s" % yes_no("+gnutls"),
            "MFEM_USE_OPENMP=%s" % yes_no("+openmp"),
            "MFEM_USE_CONDUIT=%s" % yes_no("+conduit"),
            "MFEM_USE_CUDA=%s" % yes_no("+cuda"),
            "MFEM_USE_HIP=%s" % yes_no("+rocm"),
            "MFEM_USE_OCCA=%s" % yes_no("+occa"),
            "MFEM_USE_RAJA=%s" % yes_no("+raja"),
            "MFEM_USE_AMGX=%s" % yes_no("+amgx"),
            "MFEM_USE_CEED=%s" % yes_no("+libceed"),
            "MFEM_USE_UMPIRE=%s" % yes_no("+umpire"),
            "MFEM_USE_FMS=%s" % yes_no("+fms"),
            "MFEM_USE_GINKGO=%s" % yes_no("+ginkgo"),
            "MFEM_USE_HIOP=%s" % yes_no("+hiop"),
            "MFEM_MPIEXEC=%s" % mfem_mpiexec,
            "MFEM_MPIEXEC_NP=%s" % mfem_mpiexec_np,
            "MFEM_USE_EXCEPTIONS=%s" % yes_no("+exceptions"),
            "MFEM_USE_MUMPS=%s" % yes_no("+mumps"),
        ]
        if spec.satisfies("@4.7.0:"):
            options += ["MFEM_PRECISION=%s" % spec.variants["precision"].value]

        # Determine C++ standard to use:
        cxxstd = None
        if self.spec.satisfies("@4.0.0:"):
            cxxstd = "11"
        if self.spec.satisfies("^raja@2022.03.0:"):
            cxxstd = "14"
        if self.spec.satisfies("^umpire@2022.03.0:"):
            cxxstd = "14"
        if self.spec.satisfies("^sundials@6.4.0:"):
            cxxstd = "14"
        if self.spec.satisfies("^ginkgo"):
            cxxstd = "14"
        # When rocPRIM is used (e.g. by PETSc + ROCm) we need C++14:
        if self.spec.satisfies("^rocprim@5.5.0:"):
            cxxstd = "14"
        cxxstd_req = spec.variants["cxxstd"].value
        if cxxstd_req != "auto":
            # Constraints for valid standard level should be imposed during
            # concretization based on 'conflicts' or other directives.
            cxxstd = cxxstd_req
        cxxstd_flag = None
        if cxxstd:
            if "+cuda" in spec:
                cxxstd_flag = "-std=c++" + cxxstd
            else:
                cxxstd_flag = getattr(self.compiler, "cxx" + cxxstd + "_flag")

        cuda_arch = None if "~cuda" in spec else spec.variants["cuda_arch"].value

        cxxflags = spec.compiler_flags["cxxflags"].copy()

        if cxxflags:
            # Add opt/debug flags if they are not present in global cxx flags
            opt_flag_found = any(f in self.compiler.opt_flags for f in cxxflags)
            debug_flag_found = any(f in self.compiler.debug_flags for f in cxxflags)

            if "+debug" in spec:
                if not debug_flag_found:
                    cxxflags.append("-g")
                if not opt_flag_found:
                    cxxflags.append("-O0")
            else:
                if not opt_flag_found:
                    cxxflags.append("-O2")

            cxxflags = [(xcompiler + flag) for flag in cxxflags]
            if "+cuda" in spec:
                cxxflags += [
                    "-x=cu --expt-extended-lambda -arch=sm_%s" % cuda_arch,
                    "-ccbin %s" % (spec["mpi"].mpicxx if "+mpi" in spec else env["CXX"]),
                ]
            if cxxstd_flag:
                cxxflags.append(cxxstd_flag)
            # The cxxflags are set by the spack c++ compiler wrapper. We also
            # set CXXFLAGS explicitly, for clarity, and to properly export the
            # cxxflags in the variable MFEM_CXXFLAGS in config.mk.
            options += ["CXXFLAGS=%s" % " ".join(cxxflags)]

        elif cxxstd_flag:
            options += ["BASE_FLAGS=%s" % cxxstd_flag]

        # Treat any 'CXXFLAGS' in the environment as extra c++ flags which are
        # handled through the 'CPPFLAGS' makefile variable in MFEM. Also, unset
        # 'CXXFLAGS' from the environment to prevent it from overriding the
        # defaults.
        if "CXXFLAGS" in env:
            options += ["CPPFLAGS=%s" % env["CXXFLAGS"]]
            del env["CXXFLAGS"]

        if "~static" in spec:
            options += ["STATIC=NO"]
        if "+shared" in spec:
            options += ["SHARED=YES", "PICFLAG=%s" % (xcompiler + self.compiler.cxx_pic_flag)]

        if "+mpi" in spec:
            options += ["MPICXX=%s" % spec["mpi"].mpicxx]
            hypre = spec["hypre"]
            all_hypre_libs = hypre.libs
            if "+lapack" in hypre:
                all_hypre_libs += hypre["lapack"].libs + hypre["blas"].libs

            hypre_gpu_libs = ""
            if "+cuda" in hypre:
                hypre_gpu_libs = " -lcusparse -lcurand -lcublas"
            elif "+rocm" in hypre:
                hypre_rocm_libs = LibraryList([])
                if "^rocsparse" in hypre:
                    hypre_rocm_libs += hypre["rocsparse"].libs
                if "^rocrand" in hypre:
                    hypre_rocm_libs += hypre["rocrand"].libs
                hypre_gpu_libs = " " + ld_flags_from_library_list(hypre_rocm_libs)
            options += [
                "HYPRE_OPT=-I%s" % hypre.prefix.include,
                "HYPRE_LIB=%s%s" % (ld_flags_from_library_list(all_hypre_libs), hypre_gpu_libs),
            ]

        if "+metis" in spec:
            options += [
                "METIS_OPT=-I%s" % spec["metis"].prefix.include,
                "METIS_LIB=%s" % ld_flags_from_library_list(spec["metis"].libs),
            ]

        if "+lapack" in spec:
            lapack_blas = spec["lapack"].libs + spec["blas"].libs
            options += [
                # LAPACK_OPT is not used
                "LAPACK_LIB=%s"
                % ld_flags_from_library_list(lapack_blas)
            ]

        if "+superlu-dist" in spec:
            lapack_blas = spec["lapack"].libs + spec["blas"].libs
            options += [
                "SUPERLU_OPT=-I%s -I%s"
                % (spec["superlu-dist"].prefix.include, spec["parmetis"].prefix.include),
                "SUPERLU_LIB=%s %s"
                % (
                    ld_flags_from_dirs(
                        [spec["superlu-dist"].prefix.lib, spec["parmetis"].prefix.lib],
                        ["superlu_dist", "parmetis"],
                    ),
                    ld_flags_from_library_list(lapack_blas),
                ),
            ]

        if "+strumpack" in spec:
            strumpack = spec["strumpack"]
            sp_opt = ["-I%s" % strumpack.prefix.include]
            sp_lib = [ld_flags_from_library_list(strumpack.libs)]
            # Parts of STRUMPACK use fortran, so we need to link with the
            # fortran library and also the MPI fortran library:
            if "~shared" in strumpack:
                if os.path.basename(env["FC"]) == "gfortran":
                    gfortran = Executable(env["FC"])
                    libext = "dylib" if sys.platform == "darwin" else "so"
                    libfile = os.path.abspath(
                        gfortran("-print-file-name=libgfortran.%s" % libext, output=str).strip()
                    )
                    gfortran_lib = LibraryList(libfile)
                    sp_lib += [ld_flags_from_library_list(gfortran_lib)]
                if "+mpi" in strumpack:
                    mpi = strumpack["mpi"]
                    if ("^mpich" in strumpack) or ("^mvapich2" in strumpack):
                        sp_lib += [ld_flags_from_dirs([mpi.prefix.lib], ["mpifort"])]
                    elif "^openmpi" in strumpack:
                        sp_lib += [ld_flags_from_dirs([mpi.prefix.lib], ["mpi_mpifh"])]
                    elif "^spectrum-mpi" in strumpack:
                        sp_lib += [ld_flags_from_dirs([mpi.prefix.lib], ["mpi_ibm_mpifh"])]
            if "+openmp" in strumpack:
                # The "+openmp" in the spec means strumpack will TRY to find
                # OpenMP; if not found, we should not add any flags -- how do
                # we figure out if strumpack found OpenMP?
                if not self.spec.satisfies("%apple-clang"):
                    sp_opt += [xcompiler + self.compiler.openmp_flag]
            if "^parmetis" in strumpack:
                parmetis = strumpack["parmetis"]
                sp_opt += [parmetis.headers.cpp_flags]
                sp_lib += [ld_flags_from_library_list(parmetis.libs)]
            if "^netlib-scalapack" in strumpack:
                scalapack = strumpack["scalapack"]
                sp_opt += ["-I%s" % scalapack.prefix.include]
                sp_lib += [ld_flags_from_dirs([scalapack.prefix.lib], ["scalapack"])]
            elif "^scalapack" in strumpack:
                scalapack = strumpack["scalapack"]
                sp_opt += [scalapack.headers.cpp_flags]
                sp_lib += [ld_flags_from_library_list(scalapack.libs)]
            if "+butterflypack" in strumpack:
                bp = strumpack["butterflypack"]
                sp_opt += ["-I%s" % bp.prefix.include]
                bp_libs = find_libraries(
                    ["libdbutterflypack", "libzbutterflypack"],
                    bp.prefix,
                    shared=("+shared" in bp),
                    recursive=True,
                )
                sp_lib += [ld_flags_from_library_list(bp_libs)]
            if "+zfp" in strumpack:
                zfp = strumpack["zfp"]
                sp_opt += ["-I%s" % zfp.prefix.include]
                zfp_lib = find_libraries(
                    "libzfp", zfp.prefix, shared=("+shared" in zfp), recursive=True
                )
                sp_lib += [ld_flags_from_library_list(zfp_lib)]
            if "+cuda" in strumpack:
                # assuming also ("+cuda" in spec)
                sp_lib += ["-lcusolver", "-lcublas"]
            options += [
                "STRUMPACK_OPT=%s" % " ".join(sp_opt),
                "STRUMPACK_LIB=%s" % " ".join(sp_lib),
            ]

        if "+suite-sparse" in spec:
            ss_spec = "suite-sparse:" + self.suitesparse_components
            options += [
                "SUITESPARSE_OPT=-I%s" % spec[ss_spec].prefix.include,
                "SUITESPARSE_LIB=%s" % ld_flags_from_library_list(spec[ss_spec].libs),
            ]

        if "+sundials" in spec:
            sun_spec = "sundials:" + self.sundials_components
            options += [
                "SUNDIALS_OPT=%s" % spec[sun_spec].headers.cpp_flags,
                "SUNDIALS_LIB=%s" % ld_flags_from_library_list(spec[sun_spec].libs),
            ]

        if "+petsc" in spec:
            petsc = spec["petsc"]
            if "+shared" in petsc:
                options += [
                    "PETSC_OPT=%s" % petsc.headers.cpp_flags,
                    "PETSC_LIB=%s" % ld_flags_from_library_list(petsc.libs),
                ]
            else:
                options += ["PETSC_DIR=%s" % petsc.prefix]

        if "+slepc" in spec:
            slepc = spec["slepc"]
            options += [
                "SLEPC_OPT=%s" % slepc.headers.cpp_flags,
                "SLEPC_LIB=%s" % ld_flags_from_library_list(slepc.libs),
            ]

        if "+pumi" in spec:
            pumi_libs = [
                "pumi",
                "crv",
                "ma",
                "mds",
                "apf",
                "pcu",
                "gmi",
                "parma",
                "lion",
                "mth",
                "apf_zoltan",
                "spr",
            ]
            pumi_dep_zoltan = ""
            pumi_dep_parmetis = ""
            if "+zoltan" in spec["pumi"]:
                pumi_dep_zoltan = ld_flags_from_dirs([spec["zoltan"].prefix.lib], ["zoltan"])
                if "+parmetis" in spec["zoltan"]:
                    pumi_dep_parmetis = ld_flags_from_dirs(
                        [spec["parmetis"].prefix.lib], ["parmetis"]
                    )
            options += [
                "PUMI_OPT=-I%s" % spec["pumi"].prefix.include,
                "PUMI_LIB=%s %s %s"
                % (
                    ld_flags_from_dirs([spec["pumi"].prefix.lib], pumi_libs),
                    pumi_dep_zoltan,
                    pumi_dep_parmetis,
                ),
            ]

        if "+gslib" in spec:
            options += [
                "GSLIB_OPT=-I%s" % spec["gslib"].prefix.include,
                "GSLIB_LIB=%s" % ld_flags_from_dirs([spec["gslib"].prefix.lib], ["gs"]),
            ]

        if "+netcdf" in spec:
            lib_flags = ld_flags_from_dirs([spec["netcdf-c"].prefix.lib], ["netcdf"])
            hdf5 = spec["hdf5:hl"]
            if hdf5.satisfies("~shared"):
                hdf5_libs = hdf5.libs
                hdf5_libs += LibraryList(find_system_libraries("libdl"))
                lib_flags += " " + ld_flags_from_library_list(hdf5_libs)
            options += [
                "NETCDF_OPT=-I%s" % spec["netcdf-c"].prefix.include,
                "NETCDF_LIB=%s" % lib_flags,
            ]

        if "+zlib" in spec:
            if "@:3.3.2" in spec:
                options += ["ZLIB_DIR=%s" % spec["zlib-api"].prefix]
            else:
                options += [
                    "ZLIB_OPT=-I%s" % spec["zlib-api"].prefix.include,
                    "ZLIB_LIB=%s" % ld_flags_from_library_list(spec["zlib-api"].libs),
                ]

        if "+mpfr" in spec:
            options += [
                "MPFR_OPT=-I%s" % spec["mpfr"].prefix.include,
                "MPFR_LIB=%s" % ld_flags_from_dirs([spec["mpfr"].prefix.lib], ["mpfr"]),
            ]

        if "+gnutls" in spec:
            options += [
                "GNUTLS_OPT=-I%s" % spec["gnutls"].prefix.include,
                "GNUTLS_LIB=%s" % ld_flags_from_dirs([spec["gnutls"].prefix.lib], ["gnutls"]),
            ]

        if "+libunwind" in spec:
            libunwind = spec["unwind"]
            headers = find_headers("libunwind", libunwind.prefix.include)
            headers.add_macro("-g")
            libs = find_optional_library("libunwind", libunwind.prefix)
            # When mfem uses libunwind, it also needs "libdl".
            libs += LibraryList(find_system_libraries("libdl"))
            options += [
                "LIBUNWIND_OPT=%s" % headers.cpp_flags,
                "LIBUNWIND_LIB=%s" % ld_flags_from_library_list(libs),
            ]

        if "+openmp" in spec:
            options += ["OPENMP_OPT=%s" % (xcompiler + self.compiler.openmp_flag)]

        if "+cuda" in spec:
            options += [
                "CUDA_CXX=%s" % join_path(spec["cuda"].prefix, "bin", "nvcc"),
                "CUDA_ARCH=sm_%s" % cuda_arch,
            ]
            # Check if we are using a CUDA installation where the math libs are
            # in a separate directory:
            culibs = ["libcusparse"]
            cuda_libs = find_optional_library(culibs, spec["cuda"].prefix)
            if not cuda_libs:
                p0 = os.path.realpath(join_path(spec["cuda"].prefix, "bin", "nvcc"))
                p0 = os.path.dirname(p0)
                p1 = os.path.dirname(p0)
                while p1 != p0:
                    cuda_libs = find_optional_library(culibs, join_path(p1, "math_libs"))
                    if cuda_libs:
                        break
                    p0, p1 = p1, os.path.dirname(p1)
                if not cuda_libs:
                    raise InstallError("Required CUDA libraries not found: %s" % culibs)
                options += ["CUDA_LIB=%s" % ld_flags_from_library_list(cuda_libs)]

        if "+rocm" in spec:
            amdgpu_target = ",".join(spec.variants["amdgpu_target"].value)
            options += ["HIP_CXX=%s" % spec["hip"].hipcc, "HIP_ARCH=%s" % amdgpu_target]
            hip_headers = HeaderList([])
            hip_libs = LibraryList([])
            # To use a C++ compiler that supports -xhip flag one can use
            # something like this:
            #   options += [
            #       "HIP_CXX=%s" % (spec["mpi"].mpicxx if "+mpi" in spec else spack_cxx),
            #       "HIP_FLAGS=-xhip --offload-arch=%s" % amdgpu_target,
            #   ]
            #   hip_libs += find_libraries("libamdhip64", spec["hip"].prefix.lib)
            if "^hipsparse" in spec:  # hipsparse is needed @4.4.0:+rocm
                hipsparse = spec["hipsparse"]
                hip_headers += hipsparse.headers
                hip_libs += hipsparse.libs
                # Note: MFEM's defaults.mk wants to find librocsparse.* in
                # $(HIP_DIR)/lib, so we set HIP_DIR to be $ROCM_PATH when using
                # external HIP, or the prefix of rocsparse (which is a
                # dependency of hipsparse) when using Spack-built HIP.
                if spec["hip"].external:
                    options += ["HIP_DIR=%s" % env["ROCM_PATH"]]
                else:
                    options += ["HIP_DIR=%s" % hipsparse["rocsparse"].prefix]
            if "^rocthrust" in spec and not spec["hip"].external:
                # petsc+rocm needs the rocthrust header path
                hip_headers += spec["rocthrust"].headers
            if "^rocprim" in spec and not spec["hip"].external:
                # rocthrust [via petsc+rocm] has a dependency on rocprim
                hip_headers += spec["rocprim"].headers
            if "^hipblas" in spec and not spec["hip"].external:
                # superlu-dist+rocm needs the hipblas header path
                hip_headers += spec["hipblas"].headers

            ### BEGIN AXOM PATCH                
            if "^hipblas" in spec:
                hip_libs += spec["hipblas"].libs
            ### END AXOM PATCH

            if "%cce" in spec:
                # We assume the proper Cray CCE module (cce) is loaded:
                proc = str(spec.target.family)
                craylibs_var = "CRAYLIBS_" + proc.upper()
                craylibs_path = env.get(craylibs_var, None)
                if not craylibs_path:
                    raise InstallError(
                        f"The environment variable {craylibs_var} is not defined.\n"
                        "\tMake sure the 'cce' module is in the compiler spec."
                    )
                craylibs = [
                    "libmodules",
                    "libfi",
                    "libcraymath",
                    "libf",
                    "libu",
                    "libcsup",
                    "libpgas-shmem",
                ]
                hip_libs += find_libraries(craylibs, craylibs_path)
                craylibs_path2 = join_path(craylibs_path, "../../../cce-clang", proc, "lib")
                hip_libs += find_libraries("libunwind", craylibs_path2)

            if hip_headers:
                options += ["HIP_OPT=%s" % hip_headers.cpp_flags]
            if hip_libs:
                options += ["HIP_LIB=%s" % ld_flags_from_library_list(hip_libs)]

        if "+occa" in spec:
            options += [
                "OCCA_OPT=-I%s" % spec["occa"].prefix.include,
                "OCCA_LIB=%s" % ld_flags_from_dirs([spec["occa"].prefix.lib], ["occa"]),
            ]

        if "+raja" in spec:
            raja = spec["raja"]
            raja_opt = "-I%s" % raja.prefix.include
            raja_lib = find_libraries(
                "libRAJA", raja.prefix, shared=("+shared" in raja), recursive=True
            )
            if raja.satisfies("^camp"):
                camp = raja["camp"]
                raja_opt += " -I%s" % camp.prefix.include
                raja_lib += find_optional_library("libcamp", camp.prefix)
            options += [
                "RAJA_OPT=%s" % raja_opt,
                "RAJA_LIB=%s" % ld_flags_from_library_list(raja_lib),
            ]

        if "+amgx" in spec:
            amgx = spec["amgx"]
            if "+shared" in amgx:
                options += [
                    "AMGX_OPT=-I%s" % amgx.prefix.include,
                    "AMGX_LIB=%s" % ld_flags_from_library_list(amgx.libs),
                ]
            else:
                options += ["AMGX_DIR=%s" % amgx.prefix]

        if "+libceed" in spec:
            options += [
                "CEED_OPT=-I%s" % spec["libceed"].prefix.include,
                "CEED_LIB=%s" % ld_flags_from_dirs([spec["libceed"].prefix.lib], ["ceed"]),
            ]

        if "+umpire" in spec:
            umpire = spec["umpire"]
            umpire_opts = umpire.headers
            umpire_libs = umpire.libs
            if "^camp" in umpire:
                umpire_opts += umpire["camp"].headers
            if "^fmt" in umpire:
                umpire_opts += umpire["fmt"].headers
                umpire_libs += umpire["fmt"].libs
            options += [
                "UMPIRE_OPT=%s" % umpire_opts.cpp_flags,
                "UMPIRE_LIB=%s" % ld_flags_from_library_list(umpire_libs),
            ]

        timer_ids = {"std": "0", "posix": "2", "mac": "4", "mpi": "6"}
        timer = spec.variants["timer"].value
        if timer != "auto":
            options += ["MFEM_TIMER_TYPE=%s" % timer_ids[timer]]

        if "+conduit" in spec:
            conduit = spec["conduit"]
            headers = HeaderList(find(conduit.prefix.include, "conduit.hpp", recursive=True))
            conduit_libs = ["libconduit", "libconduit_relay", "libconduit_blueprint"]
            libs = find_libraries(conduit_libs, conduit.prefix.lib, shared=("+shared" in conduit))
            libs += LibraryList(find_system_libraries("libdl"))
            if "+hdf5" in conduit:
                hdf5 = conduit["hdf5"]
                headers += find_headers("hdf5", hdf5.prefix.include)
                libs += hdf5.libs

            ##################
            # cyrush note:
            ##################
            # spack's HeaderList is applying too much magic, undermining us:
            #
            #  It applies a regex to strip back to the last "include" dir
            #  in the path. In our case we need to pass the following
            #  as part of the CONDUIT_OPT flags:
            #
            #    -I<install_path>/include/conduit
            #
            #  I tried several ways to present this path to the HeaderList,
            #  but the regex always kills the trailing conduit dir
            #  breaking build.
            #
            #  To resolve the issue, we simply join our own string with
            #  the headers results (which are important b/c they handle
            #  hdf5 paths when enabled).
            ##################

            # construct proper include path
            conduit_include_path = conduit.prefix.include.conduit
            # add this path to the found flags
            conduit_opt_flags = "-I{0} {1}".format(conduit_include_path, headers.cpp_flags)

            options += [
                "CONDUIT_OPT=%s" % conduit_opt_flags,
                "CONDUIT_LIB=%s" % ld_flags_from_library_list(libs),
            ]

        if "+fms" in spec:
            libfms = spec["libfms"]
            options += [
                "FMS_OPT=%s" % libfms.headers.cpp_flags,
                "FMS_LIB=%s" % ld_flags_from_library_list(libfms.libs),
            ]

        if "+ginkgo" in spec:
            ginkgo = spec["ginkgo"]
            options += [
                "GINKGO_DIR=%s" % ginkgo.prefix,
                "GINKGO_BUILD_TYPE=%s" % ginkgo.variants["build_type"].value,
            ]

        if "+hiop" in spec:
            hiop = spec["hiop"]
            hiop_hdrs = hiop.headers
            hiop_libs = hiop.libs
            hiop_hdrs += spec["lapack"].headers + spec["blas"].headers
            hiop_libs += spec["lapack"].libs + spec["blas"].libs
            hiop_opt_libs = ["magma", "umpire", "hipblas", "hiprand"]
            for opt_lib in hiop_opt_libs:
                if "^" + opt_lib in hiop:
                    hiop_hdrs += hiop[opt_lib].headers
                    hiop_libs += hiop[opt_lib].libs
            # raja's libs property does not work
            if "^raja" in hiop:
                raja = hiop["raja"]
                hiop_hdrs += raja.headers
                hiop_libs += find_libraries(
                    "libRAJA", raja.prefix, shared=("+shared" in raja), recursive=True
                )
                if raja.satisfies("^camp"):
                    camp = raja["camp"]
                    hiop_hdrs += camp.headers
                    hiop_libs += find_optional_library("libcamp", camp.prefix)
            if hiop.satisfies("@0.6:+cuda"):
                hiop_libs += LibraryList(["cublas", "curand"])
            options += [
                "HIOP_OPT=%s" % hiop_hdrs.cpp_flags,
                "HIOP_LIB=%s" % ld_flags_from_library_list(hiop_libs),
            ]

        if "+mumps" in spec:
            mumps = spec["mumps"]
            mumps_opt = ["-I%s" % mumps.prefix.include]
            if "+openmp" in mumps:
                if not self.spec.satisfies("%apple-clang"):
                    mumps_opt += [xcompiler + self.compiler.openmp_flag]
            options += [
                "MUMPS_OPT=%s" % " ".join(mumps_opt),
                "MUMPS_LIB=%s" % ld_flags_from_library_list(mumps.libs),
            ]

        return options