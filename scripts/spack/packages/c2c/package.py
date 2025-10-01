import os
import socket

from spack.package import *
from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack_repo.builtin.build_systems.cmake import CMakeBuilder
from spack_repo.builtin.build_systems.generic import Package

class C2c(CMakePackage):
    """Contour Parser Library"""

    homepage = 'https://rzlc.llnl.gov/c2c'
    url = 'file:///collab/usr/gapps/c2c/dist/c2c-0.11.0.tgz'

    version('develop', git='ssh://git@rzgitlab.llnl.gov:7999/c2c/c2c.git',
            submodules=True, branch='develop')

    version('1.8.0', sha256='6421a5314572f9e241514b36b1738289f2b556393db47ca3cac66f3428165cf4')
    version('1.3.0', sha256='8dd3ba73651f61d2f6d0dea43ab9568a566fd73c591c8bcae8792b44cb654c59')
    version('0.12.0', '250aa61251589215dd7f90f7b31dd443')
    version('0.11.0', '28cb091711af15bcb27e01b5ed2a4921')
    version('0.10.0', '99d797383f32ece809c12682fd5ce5cc')
    version('0.9.1', 'a2de3857e9581f63eb1f411e165d7e3d')
    version('0.9.0', 'a7ad2f3c8121c585d2bcd7edc22542a5')
    version('0.8.2', '934244b225e85ed335455c9dd22df961')
    version('0.8.0', '5d40e39def358115232d8b7fd6b0032a')

    variant('dev', default=False,
            description='For developers. Needed for modifying parser')
    variant('docs', default=False,
            description='Allow generating documentation')
    variant('tools', default=False,
            description='Build command-line tools in addition to library')

    depends_on("c", type="build")
    depends_on("cxx", type="build")

    depends_on('cmake@3.8.0:', type='build')
    depends_on('doxygen', type='build', when='+docs')
    depends_on('bison@3.6.0:', type='build', when='+dev')
    depends_on('flex@2.6.0:', type='build', when='+dev')
    depends_on('nlohmann-json', when='+tools')

    phases = ('hostconfig',) + CMakeBuilder.phases

    def configure_args(self):
        spec = self.spec if self.spec is not None else ""
        return [
            '-DBUILD_TOOLS={0}'.format('YES' if '+tools' in spec else 'NO'),
            '-DGENERATE_PARSER={0}'.format('YES' if '+dev' in spec else 'NO'),
            '-DRUN_TESTS={0}'.format('YES' if '+dev' in spec else 'NO'),
            '-DBUILD_DOCS={0}'.format('YES' if '+docs' in spec else 'NO'),
        ]


    def hostconfig(self, spec, prefix, py_site_pgs_dir=None):
        """
        Create a 'host-config' file with all the options used to configure
        this package.
        """

        host_config_path = self._get_host_config_path(spec)
        with open(host_config_path, 'w') as cfg:
            self._write_header(spec, cfg)
            self._write_compiler_settings(spec, cfg)
            self._write_dev_tools(spec, prefix, cfg)

    def _write_header(self, spec, cfg):
            cfg.write("####################################################################\n")
            cfg.write("# Generated host-config - Edit at own risk!\n")
            cfg.write("####################################################################\n")

            cmake_exe = spec['cmake'].command.path
            cmake_exe = os.path.realpath(cmake_exe)
            cfg.write("#---------------------------------------\n")
            cfg.write("# SYS_TYPE: {0}\n".format(self._get_sys_type(spec)))
            cfg.write("# Compiler Spec: {0}\n".format(spec.compiler))
            cfg.write("# CMake executable path: %s\n" % cmake_exe)
            cfg.write("#---------------------------------------\n\n")

    def _write_compiler_settings(self, spec, cfg):
            cfg.write("#---------------------------------------\n")
            cfg.write("# Compilers\n")
            cfg.write("#---------------------------------------\n")
            cfg.write(cmake_cache_path("CMAKE_C_COMPILER", env['SPACK_CC']))
            cfg.write(cmake_cache_path("CMAKE_CXX_COMPILER", env['SPACK_CXX']))

            # use global spack compiler flags
            cflags = ' '.join(spec.compiler_flags['cflags'])
            if cflags:
                cfg.write(cmake_cache_entry("CMAKE_C_FLAGS", cflags))
            cxxflags = ' '.join(spec.compiler_flags['cxxflags'])
            if cxxflags:
                cfg.write(cmake_cache_entry("CMAKE_CXX_FLAGS", cxxflags))

    def _write_dev_tools(self, spec, prefix, cfg):
            cfg.write("#---------------------------------------\n")
            cfg.write("# Library Dependencies\n")
            cfg.write("#---------------------------------------\n")

            path_replacements = {}

            # Try to find the common prefix of the TPL directory, including the compiler
            # If found, we will use this in the TPL paths
            compiler_str = str(spec.compiler).replace('@','-')
            prefix_paths = prefix.split(compiler_str)
            tpl_root = ""
            if len(prefix_paths) == 2:
                tpl_root = os.path.join( prefix_paths[0], compiler_str )
                path_replacements[tpl_root] = "${TPL_ROOT}"
                cfg.write(cmake_cache_path("TPL_ROOT", tpl_root))

            cfg.write("#------------------{}\n".format("-"*60))
            cfg.write("# Devtools\n")
            cfg.write("#------------------{}\n\n".format("-"*60))

            # Add common prefix to path replacement list
            if '+dev' in spec:
                # Grab common devtools root and strip the trailing slash
                path1 = os.path.realpath(spec['flex'].prefix)
                path2 = os.path.realpath(spec['bison'].prefix)
                devtools_root = os.path.commonprefix([path1, path2])
                if len(devtools_root) > 1:
                    devtools_root = devtools_root[:-1]
                    path_replacements[devtools_root] = "${DEVTOOLS_ROOT}"
                    cfg.write("# Root directory for generated developer tools\n")
                    cfg.write(cmake_cache_path("DEVTOOLS_ROOT", devtools_root))

                cfg.write(cmake_cache_option('GENERATE_PARSER', True))
                cfg.write(cmake_cache_option('RUN_TESTS', True))

                flex_bin_dir = get_spec_path(spec, 'flex', path_replacements, use_bin=True)
                cfg.write(cmake_cache_path('FLEX_EXECUTABLE', os.path.join(flex_bin_dir, 'flex')))

                bison_bin_dir = get_spec_path(spec, 'bison', path_replacements, use_bin=True)
                cfg.write(cmake_cache_path('BISON_EXECUTABLE', os.path.join(bison_bin_dir, 'bison')))

            if '+tools' in spec:
                cfg.write(cmake_cache_option('BUILD_TOOLS', True))

                # If we just set NLOHMANN_JSON_DIR, cmake can't find the
                # library. We can set this directly, but sometimes need it
                # to point under the lib directory, and others uners lib64.
                # Appending this to CMAKE_PREFIX_PATH works.
                nlohmann_dir = get_spec_path(spec, 'nlohmann-json', path_replacements)
                cfg.write(cmake_cache_entry('CMAKE_PREFIX_PATH', '${CMAKE_PREFIX_PATH};' + nlohmann_dir))

            if '+docs' in spec:
                cfg.write(cmake_cache_option('BUILD_DOCS', True))

                doxygen_bin_dir = get_spec_path(spec, 'doxygen', path_replacements, use_bin=True)
                cfg.write(cmake_cache_path('DOXYGEN_EXECUTABLE', os.path.join(doxygen_bin_dir, 'doxygen')))

    def _get_host_config_path(self, spec):
        var=''
        host_config_path = "%s-%s-%s%s.cmake" % (socket.gethostname().rstrip('1234567890'),
                                                 self._get_sys_type(spec),
                                                 spec.compiler.name, var)
        dest_dir = self.stage.source_path
        host_config_path = os.path.abspath(os.path.join(dest_dir, host_config_path))
        return host_config_path

    def _get_sys_type(self, spec):
        # if on llnl systems, we can use the SYS_TYPE
        return os.environ.get('SYS_TYPE', spec.architecture)


def cmake_cache_path(name, value, comment=""):
    """Generate a string for a cmake cache variable"""
    return 'set(%s "%s" CACHE PATH "%s")\n\n' % (name, value, comment)


def cmake_cache_entry(name, value, comment=""):
    """Generate a string for a cmake cache variable"""
    return 'set(%s "%s" CACHE STRING "%s")\n\n' % (name, value, comment)


def cmake_cache_option(name, boolean_value, comment=""):
    """Generate a string for a cmake configuration option"""
    value = "ON" if boolean_value else "OFF"
    return 'set(%s %s CACHE BOOL "%s")\n\n' % (name, value, comment)


def get_spec_path(spec, package_name, path_replacements = {}, use_bin = False) :
    """Extracts the prefix path for the given spack package
       path_replacements is a dictionary with string replacements for the path.
    """

    if not use_bin:
        path = spec[package_name].prefix
    else:
        path = spec[package_name].prefix.bin

    path = os.path.realpath(path)

    for key in path_replacements:
        path = path.replace(key, path_replacements[key])

    return path

