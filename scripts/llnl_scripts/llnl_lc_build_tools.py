# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: llnl_lc_uberenv_install_tools.py

 description:
  helpers for installing axom tpls on llnl lc systems.

"""

import os
import socket
import sys
import subprocess
import datetime
import glob
import json
import getpass
import shutil
import time

from os.path import join as pjoin

def sexe(cmd,
         ret_output=False,
         output_file = None,
         echo = False,
         error_prefix = "ERROR:"):
    """ Helper for executing shell commands. """
    if echo:
        print("[exe: %s]" % cmd)
    if ret_output:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        res =p.communicate()[0]
        if isinstance(res, bytes):
            res = res.decode()
        return p.returncode,res
    elif output_file != None:
        ofile = open(output_file,"w")
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout= ofile,
                             stderr=subprocess.STDOUT)
        res =p.communicate()[0]
        return p.returncode
    else:
        rcode = subprocess.call(cmd,shell=True)
        if rcode != 0:
            print("[{0} [return code: {1}] from command: {2}]".format(error_prefix, rcode,cmd))
        return rcode


def get_timestamp(t=None,sep="_"):
    """ Creates a timestamp that can easily be included in a filename. """
    if t is None:
        t = datetime.datetime.now()
    sargs = (t.year,t.month,t.day,t.hour,t.minute,t.second)
    sbase = "".join(["%04d",sep,"%02d",sep,"%02d",sep,"%02d",sep,"%02d",sep,"%02d"])
    return  sbase % sargs


def build_info():
    res = {}
    res["built_by"] = os.environ["USER"]
    res["built_from_branch"] = "unknown"
    res["built_from_sha1"]   = "unknown"
    res["platform"] = get_platform()
    res["hostname"] = "unknown"
    rc, out = sexe('git branch -a | grep \"*\"',ret_output=True,error_prefix="WARNING:")
    out = out.strip()
    if rc == 0 and out != "":
        res["built_from_branch"]  = out.split()[1]
    rc,out = sexe('git rev-parse --verify HEAD',ret_output=True,error_prefix="WARNING:")
    out = out.strip()
    if rc == 0 and out != "":
        res["built_from_sha1"] = out
    rc, out = sexe('hostname',ret_output=True,error_prefix="WARNING:")
    out = out.strip()
    if rc == 0 and out != "":
        res["hostname"] = out
    return res


def write_build_info(ofile):
    print("[build info]")
    binfo_str = json.dumps(build_info(),indent=2)
    print(binfo_str)
    open(ofile,"w").write(binfo_str)


def log_success(prefix, msg, timestamp=""):
    """
    Called at the end of the process to signal success.
    """
    print(msg)
    info = {}
    info["prefix"] = prefix
    info["platform"] = get_platform()
    info["status"] = "success"
    info["message"] = msg
    if timestamp == "":
        info["timestamp"] = get_timestamp()
    else:
        info["timestamp"] = timestamp
    json.dump(info,open(pjoin(prefix,"success.json"),"w"),indent=2)


def log_failure(prefix, msg, timestamp=""):
    """
    Called when the process failed.
    """
    print(msg)
    info = {}
    info["prefix"] = prefix
    info["platform"] = get_platform()
    info["status"] = "failed"
    info["message"] = msg
    if timestamp == "":
        info["timestamp"] = get_timestamp()
    else:
        info["timestamp"] = timestamp
    json.dump(info,open(pjoin(prefix,"failed.json"),"w"),indent=2)


def assertUberenvExists():
    if not os.path.exists(get_uberenv_path()):
        print("[ERROR: {0} does not exist".format(get_uberenv_path()))
        print("  run 'git submodule update --init'")
        print("]")
        sys.exit(1)


def uberenv_create_mirror(prefix, spec, project_file, mirror_path):
    """
    Calls uberenv to create a spack mirror.
    """
    assertUberenvExists()
    cmd  = "{0} {1} --create-mirror -k ".format(sys.executable, get_uberenv_path())
    cmd += "--prefix=\"{0}\" --mirror=\"{1}\" ".format(prefix, mirror_path)
    cmd += "--spec=\"{0}\" ".format(spec)
    if project_file:
        cmd += "--project-json=\"{0}\" ".format(project_file)

    print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")
    print("[ It is expected for 'spack --create-mirror' to throw warnings.                ]")
    print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")
    res = sexe(cmd, echo=True, error_prefix="WARNING:")
    print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")
    print("[ End of expected warnings from 'spack --create-mirror'                        ]")
    print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")
    set_group_and_perms(mirror_path)
    return res


def uberenv_build(prefix, spec, project_file, mirror_path):
    """
    Calls uberenv to install tpls for a given spec to given prefix.
    """
    assertUberenvExists()
    cmd  = "{0} {1} -k ".format(sys.executable, get_uberenv_path())
    cmd += "--prefix=\"{0}\" --spec=\"{1}\" ".format(prefix, spec)
    cmd += "--mirror=\"{0}\" ".format(mirror_path)
    if project_file:
        cmd += "--project-json=\"{0}\" ".format(project_file)

    # File names have a 255 character limit - take the first 200 characters if spec is too long
    truncated_spec = spec[:200]

    spack_tpl_build_log = pjoin(prefix,"output.log.spack.tpl.build.%s.txt" % truncated_spec.replace(" ", "_"))

    print("[starting tpl install of spec %s]" % spec)
    print("[log file: %s]" % spack_tpl_build_log)
    res = sexe(cmd,
               echo=True,
               output_file = spack_tpl_build_log)

    # Move files generated by spack in source directory to TPL install directory
    print("[Moving spack generated files to TPL build directory]")
    repo_dir = get_repo_dir()
    for file in ["spack-build-env.txt", "spack-build-out.txt", "spack-configure-args.txt"]:
        src = pjoin(repo_dir, file)
        dst = pjoin(prefix, "{0}-{1}".format(truncated_spec.replace(" ", "_"),file))
        if os.path.exists(src) and not os.path.exists(dst):
            shutil.move(src, dst)

    if res != 0:
        log_failure(prefix,"[ERROR: uberenv/spack build of spec: %s failed]" % spec)
    return res


############################################################
# helpers for testing a set of host configs
############################################################

def build_and_test_host_config(test_root, host_config, 
                               report_to_stdout = False,
                               extra_cmake_options = "",
                               build_type = "Debug",
                               test_serial = False):
    host_config_root = get_host_config_root(host_config)
    # setup build and install dirs
    build_dir   = pjoin(test_root,"build-%s"   % host_config_root)
    install_dir = pjoin(test_root,"install-%s" % host_config_root)
    print("[Testing build, test, and install of host config file: %s]" % host_config)
    print("[ build dir: %s]"   % build_dir)
    print("[ install dir: %s]" % install_dir)

    # configure
    cfg_output_file = pjoin(test_root,"output.log.%s.configure.txt" % host_config_root)
    print("[starting configure of %s]" % host_config)
    print("[log file: %s]" % cfg_output_file)
    res = sexe("%s config-build.py -bp %s -ip %s -bt %s -hc %s %s" % (sys.executable, build_dir, install_dir, build_type, host_config, extra_cmake_options),
               output_file = cfg_output_file,
               echo=True)

    if report_to_stdout:
        with open(cfg_output_file, 'r', encoding='utf8') as build_out:
            print(build_out.read())

    if res != 0:
        print("[ERROR: Configure for host-config: %s failed]\n" % host_config)
        return res

    ####
    # build, test, and install
    ####

    # build the code
    bld_output_file =  pjoin(build_dir,"output.log.make.txt")
    print("[starting build]")
    print("[log file: %s]" % bld_output_file)
    res = sexe("cd %s && make -j 16 VERBOSE=1 " % build_dir,
                output_file = bld_output_file,
                echo=True)

    if report_to_stdout:
        with open(bld_output_file, 'r', encoding='utf8') as build_out:
            print(build_out.read())

    if res != 0:
        print("[ERROR: Build for host-config: %s failed]\n" % host_config)
        return res

    # test the code
    tst_output_file = pjoin(build_dir,"output.log.make.test.txt")
    print("[starting unit tests]")
    print("[log file: %s]" % tst_output_file)

    parallel_test = "" if test_serial else "-j16"
    tst_cmd = "cd %s && make CTEST_OUTPUT_ON_FAILURE=1 test ARGS=\"--no-compress-output -T Test -VV %s\"" % (build_dir, parallel_test)

    res = sexe(tst_cmd,
               output_file = tst_output_file,
               echo=True)

    if report_to_stdout:
        with open(tst_output_file, 'r', encoding='utf8') as test_out:
            print(test_out.read())

    # Convert CTest output to JUnit, do not overwrite previous res
    print("[Checking to see if xsltproc exists...]")
    test_xsltproc_res = sexe("xsltproc --version", echo=True)
    if test_xsltproc_res != 0:
        print("[WARNING: xsltproc does not exist skipping JUnit conversion]")
    else:
        junit_file = pjoin(build_dir, "junit.xml")
        xsl_file = pjoin(get_blt_dir(), "tests/ctest-to-junit.xsl")
        ctest_file = pjoin(build_dir, "Testing/*/Test.xml")

        print("[Converting CTest XML to JUnit XML]")
        convert_cmd  = "xsltproc -o {0} {1} {2}".format(junit_file, xsl_file, ctest_file)
        convert_res = sexe(convert_cmd, echo=True)
        if convert_res != 0:
            print("[WARNING: Converting to JUnit failed.]")

    if res != 0:
        print("[ERROR: Tests for host-config: %s failed]\n" % host_config)
        return res

    # build the docs
    docs_output_file = pjoin(build_dir,"output.log.make.docs.txt")
    print("[starting docs generation]")
    print("[log file: %s]" % docs_output_file)

    res = sexe("cd %s && make -j16 docs " % build_dir,
               output_file = docs_output_file,
               echo=True)

    if report_to_stdout:
        with open(docs_output_file, 'r', encoding='utf8') as docs_out:
            print(docs_out.read())

    if res != 0:
        print("[ERROR: Docs generation for host-config: %s failed]\n\n" % host_config)
        return res

    # install the code
    inst_output_file = pjoin(build_dir,"output.log.make.install.txt")
    print("[starting install]")
    print("[log file: %s]" % inst_output_file)

    res = sexe("cd %s && make -j16 install " % build_dir,
               output_file = inst_output_file,
               echo=True)

    if res != 0:
        print("[ERROR: Install for host-config: %s failed]\n\n" % host_config)
        return res

    # simple sanity check for make install
    print("[checking install dir %s]" % install_dir)
    sexe("ls %s/include" % install_dir, echo=True, error_prefix="WARNING:")
    sexe("ls %s/lib" %     install_dir, echo=True, error_prefix="WARNING:")
    sexe("ls %s/bin" %     install_dir, echo=True, error_prefix="WARNING:")

    # test the installation using installed cmake examples
    # TODO: enable tests for installed makefile-based example
    should_test_installed_cmake_example = True
    should_test_installed_blt_example = True
    should_test_installed_make_example = False
    should_test_installed_tutorial = True

    if should_test_installed_cmake_example:
        install_example_dir = pjoin(install_dir, "examples", "axom", "using-with-cmake")
        install_example_output_file = pjoin(build_dir,"output.log.install_example.cmake.txt")
        print("[testing installed 'using-with-cmake' example]")
        print("[log file: %s]" % install_example_output_file)

        example_commands = [
            "cd {0}".format(install_example_dir),
            "rm -rf build",
            "mkdir build",
            "cd build",
            """echo "[Configuring '{}' example]" """.format("using-with-cmake"),
            "cmake -C ../host-config.cmake ..",
            """echo "[Building '{}' example]" """.format("using-with-cmake"),
            "make ",
            """echo "[Running '{}' example]" """.format("using-with-cmake"),
            "./example",
            """echo "[Done]" """
        ]

        res = sexe(" && ".join(example_commands),
                output_file = install_example_output_file,
                echo=True)

        if res != 0:
            print("[ERROR: Installed 'using-with-cmake' example for host-config: %s failed]\n\n" % host_config)
            return res


    if should_test_installed_blt_example:
        install_example_dir = pjoin(install_dir, "examples", "axom", "using-with-blt")
        install_example_output_file = pjoin(build_dir,"output.log.install_example.blt.txt")
        print("[testing installed 'using-with-blt' example]")
        print("[log file: %s]" % install_example_output_file)

        example_commands = [
            "cd {0}".format(install_example_dir),
            "rm -rf build",
            "mkdir build",
            "cd build",
            """echo "[Configuring '{}' example]" """.format("using-with-blt"),
            "cmake -C ../host-config.cmake ..",
            """echo "[Building '{}' example]" """.format("using-with-blt"),
            "make ",
            """echo "[Running '{}' example]" """.format("using-with-blt"),
            "./bin/example",
            """echo "[Done]" """
        ]

        res = sexe(" && ".join(example_commands),
                output_file = install_example_output_file,
                echo=True)

        if res != 0:
            print("[ERROR: Installed 'using-with-blt' example for host-config: %s failed]\n\n" % host_config)
            return res

    if should_test_installed_tutorial:
        # Notes: Might require loading module for more recent cmake than the system default for some platforms
        install_example_dir = pjoin(install_dir, "examples", "axom", "radiuss_tutorial")
        install_example_output_file = pjoin(build_dir,"output.log.install_example.radiuss_tutorial.txt")
        print("[testing installed 'radiuss_tutorial' example]")
        print("[log file: %s]" % install_example_output_file)

        example_commands = [
            "cd {0}".format(install_example_dir),
            "rm -rf build",
            "mkdir build",
            "cd build",
            """echo "[Configuring 'radiuss_tutorial']" """,
            "cmake -C ../host-config.cmake -DCMAKE_BUILD_TYPE={0} ..".format(build_type),
            """echo "[Building 'radiuss_tutorial']" """,
            "make -j16",
            """echo "[Running lessons for radiuss_tutorial]" """,
            "ctest -j16",
            """echo "[Done]" """
        ]

        res = sexe(" && ".join(example_commands),
                output_file = install_example_output_file,
                echo=True)

        if res != 0:
            print("[ERROR: Installed 'radiuss_tutorial' for host-config: %s failed]\n\n" % host_config)
            return res

    print("[SUCCESS: Build, test, and install for host-config: %s complete]\n" % host_config)

    set_group_and_perms(build_dir)
    set_group_and_perms(install_dir)

    return 0


def build_and_test_host_configs(prefix,
                                timestamp,
                                use_generated_host_configs,
                                report_to_stdout = False,
                                extra_cmake_options = "",
                                build_type = "Debug",
                                test_serial = False):
    host_configs = get_host_configs_for_current_machine(prefix, use_generated_host_configs)
    if len(host_configs) == 0:
        log_failure(prefix,"[ERROR: No host configs found at %s]" % prefix)
        return 1
    print("Found Host-configs:")
    for host_config in host_configs:
        print("    " + host_config)
    print("\n")

    test_root =  get_build_and_test_root(prefix, timestamp)
    os.mkdir(test_root)
    write_build_info(pjoin(test_root,"info.json")) 
    ok  = []
    bad = []
    for host_config in host_configs:
        build_dir = get_build_dir(test_root, host_config)

        start_time = time.time()
        if build_and_test_host_config(test_root, host_config,
                                      report_to_stdout = report_to_stdout,
                                      extra_cmake_options=extra_cmake_options,
                                      build_type = build_type,
                                      test_serial = test_serial) == 0:
            ok.append(host_config)
            log_success(build_dir, "[Success: Built host-config: {0}]".format(host_config), timestamp)
        else:
            bad.append(host_config)
            log_failure(build_dir, "[Error: Failed to build host-config: {0}]".format(host_config), timestamp)
        end_time = time.time()
        print("[build time: {0}]\n".format(convertSecondsToReadableTime(end_time - start_time)))

    # Log overall job success/failure
    if len(bad) != 0:
        log_failure(test_root, "[Error: Failed to build host-configs: {0}]".format(bad), timestamp)
    else:
        log_success(test_root,"[Success: Built all host-configs: {0}]".format(ok), timestamp)

    # Output summary of failure/succesful builds
    if len(ok) > 0:
        print("Succeeded:")
        for host_config in ok:
            print("    " + host_config)

    if len(bad) > 0:
        print("Failed:")
        for host_config in bad:
            print("    " + host_config)
        print("\n")
        return 1

    print("\n")

    return 0


def set_group_and_perms(directory):
    """
    Sets the proper group and access permissions of given input directory.
    """

    skip = True
    shared_dirs = [get_shared_base_dir(), get_shared_collab_dir()]
    for shared_dir in shared_dirs:
        if directory.startswith(shared_dir):
            skip = False
            break

    if skip:
        print("[Skipping update of group and access permissions. Provided directory was not a known shared location: {0}]".format(directory))
    else:
        print("[changing group and access perms of: %s]" % directory)
        print("[changing group to axom]")
        sexe("chgrp -f -R axom %s" % (directory),echo=True,error_prefix="WARNING:")
        print("[changing perms for axom members to 'rwX' and all to 'rX']")
        sexe("chmod -f -R g+rwX,a+rX %s" % (directory),echo=True,error_prefix="WARNING:")
        print("[done setting perms for: %s]" % directory)
    return 0


def full_build_and_test_of_tpls(builds_dir, timestamp, spec, report_to_stdout = False, mirror_location = ''):
    if spec:
        specs = [spec]
    else:
        specs = get_specs_for_current_machine()
    print("[Building and testing tpls for specs: ")
    for spec in specs:
        print("{0}".format(spec))
    print("]\n")

    # Use shared network mirror location otherwise create local one
    if mirror_location:
        mirror_dir = mirror_location
    else:
        mirror_dir = get_shared_mirror_dir()
        if not os.path.exists(mirror_dir):
            mirror_dir = pjoin(builds_dir,"mirror")
        print("[using mirror location: %s]" % mirror_dir)

    # unique install location
    prefix = pjoin(builds_dir, get_system_type())
    if not os.path.exists(prefix):
        os.mkdir(prefix)
    prefix = pjoin(prefix, timestamp)
    if not os.path.exists(prefix):
        os.mkdir(prefix)

    # create a mirror
    uberenv_create_mirror(prefix, spec, "", mirror_dir)
    # write info about this build
    write_build_info(pjoin(prefix, "info.json"))

    repo_dir = get_repo_dir()
    # Clean previously generated host-configs into TPL install directory
    print("[Cleaning previously generated host-configs if they exist]")
    host_configs = get_host_configs_for_current_machine(repo_dir, True)
    for host_config in host_configs:
        os.remove(host_config)

    # use uberenv to install for all specs
    tpl_build_failed = False
    for spec in specs:
        start_time = time.time()
        fullspec = "{0}".format(spec)
        res = uberenv_build(prefix, fullspec, "", mirror_dir)
        end_time = time.time()
        print("[build time: {0}]".format(convertSecondsToReadableTime(end_time - start_time)))
        if res != 0:
            print("[ERROR: Failed build of tpls for spec %s]\n" % spec)
            tpl_build_failed = True
            break
        else:
            print("[SUCCESS: Finished build tpls for spec %s]\n" % spec)

    # Copy generated host-configs into TPL install directory
    print("[Copying spack generated host-configs to TPL build directory]")
    host_configs = get_host_configs_for_current_machine(repo_dir, True)
    for host_config in host_configs:
        dst = pjoin(prefix, os.path.basename(host_config))
        if os.path.exists(host_config) and not os.path.exists(dst):
            shutil.copy2(host_config, dst)

    src_build_failed = False
    if not tpl_build_failed:
        # build the src against the new tpls
        res = build_and_test_host_configs(prefix, timestamp, True, report_to_stdout)
        if res != 0:
            print("[ERROR: Build and test of src vs tpls test failed.]\n")
            src_build_failed = True
        else:
            print("[SUCCESS: Build and test of src vs tpls test passed.]\n")

    # set proper perms for installed tpls
    set_group_and_perms(prefix)

    if tpl_build_failed:
        print("[ERROR: Failed to build all specs of third party libraries]")
    if src_build_failed:
        print("[ERROR: Failed to build all specs of source code against new host-configs]")
    return res


def build_devtools(builds_dir, timestamp):
    sys_type = get_system_type()
    project_file = "scripts/spack/devtools.json"

    if "toss_4" in sys_type:
        compiler_spec = "%gcc@10.3.1"
    elif "blueos" in sys_type:
        compiler_spec = "%gcc@8.3.1"

    print("[Building devtools using compiler spec: {0}]".format(compiler_spec))

    # unique install location
    prefix = pjoin(builds_dir, sys_type)
    if not os.path.exists(prefix):
        os.mkdir(prefix)
    prefix = pjoin(prefix, timestamp)
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    # Use shared mirror
    mirror_dir = get_shared_mirror_dir()
    print("[Using mirror location: {0}]".format(mirror_dir))
    uberenv_create_mirror(prefix, compiler_spec, project_file, mirror_dir)

    # write info about this build
    write_build_info(pjoin(prefix,"info.json"))

    # use uberenv to install devtools
    start_time = time.time()
    res = uberenv_build(prefix, compiler_spec, project_file, mirror_dir)
    end_time = time.time()

    print("[Build time: {0}]".format(convertSecondsToReadableTime(end_time - start_time)))
    if res != 0:
        print("[ERROR: Failed build of devtools for spec %s]\n" % compiler_spec)
    else:
        # Only update the latest symlink if successful
        link_path = pjoin(builds_dir, sys_type, "latest")
        view_dir = pjoin(prefix, "view")
        print("[Creating symlink to latest devtools view:\n{0}\n->\n{1}]".format(link_path, view_dir))
        if os.path.exists(link_path) or os.path.islink(link_path):
            if not os.path.islink(link_path):
                print("[ERROR: Latest devtools link path exists and is not a link: {0}".format(link_path))
                return 1
            os.unlink(link_path)
        os.symlink(view_dir, link_path)

        print("[SUCCESS: Finished build devtools for spec %s]\n" % compiler_spec)

    # set proper perms for installed devtools
    set_group_and_perms(prefix)

    return res


def get_specs_for_current_machine():
    repo_dir = get_repo_dir()
    specs_json_path = pjoin(repo_dir, "scripts/spack/specs.json")

    with open(specs_json_path, 'r') as f:
        specs_json = json.load(f)

    sys_type = get_system_type()
    machine_name = get_machine_name()

    specs = []
    if machine_name in specs_json.keys():
        specs = specs_json[machine_name]
    else:
        specs = specs_json[sys_type]

    specs = ['%' + spec for spec in specs]

    return specs


def get_host_configs_for_current_machine(src_dir, use_generated_host_configs):
    host_configs = []

    # Generated host-configs will be at the base of the source repository
    host_configs_dir = src_dir
    if not use_generated_host_configs:
        host_configs_dir = pjoin(src_dir, "host-configs")

    hostname_base = get_machine_name()
    host_configs = glob.glob(pjoin(host_configs_dir, hostname_base + "*.cmake"))

    return host_configs


def get_host_config_root(host_config):
    return os.path.splitext(os.path.basename(host_config))[0]


def get_blt_dir():
    return "src/cmake/blt"


def get_build_dir(prefix, host_config):
    host_config_root = get_host_config_root(host_config)
    return pjoin(prefix, "build-" + host_config_root)


def get_repo_dir():
    script_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.abspath(pjoin(script_dir, "../.."))


def get_build_and_test_root(prefix, timestamp):
    dirname = "_axom_build_and_test_{0}".format(timestamp)
    return pjoin(prefix, dirname)


def get_machine_name():
    return socket.gethostname().rstrip('1234567890')


def get_system_type():
    return os.environ["SYS_TYPE"]


def get_platform():
    return get_system_type() if "SYS_TYPE" in os.environ else get_machine_name()


def get_supported_sys_types():
    return ["blueos_3_ppc64le_ib_p9", "darwin-x86_64", "toss_4_x86_64_ib", "toss_4_x86_64_ib_cray"]

def get_username():
    return getpass.getuser()


def get_shared_base_dir():
    return "/usr/WS1/axom"


def get_shared_collab_dir():
    return "/collab/usr/gapps/axom"


def get_shared_mirror_dir():
    return pjoin(get_shared_base_dir(), "mirror")


def get_shared_libs_dir():
    return pjoin(get_shared_base_dir(), "libs")


def get_shared_devtool_dir():
    return pjoin(get_shared_collab_dir(), "devtools")


def get_uberenv_path():
    return pjoin(get_script_dir(), "../uberenv/uberenv.py")


def on_rz():
    machine_name = get_machine_name()
    if machine_name.startswith("rz"):
        return True
    return False


def get_script_dir():
    return os.path.dirname(os.path.abspath(__file__))


def get_compiler_from_spec(spec):
    compiler = spec
    for c in ['~', '+']:
        index = compiler.find(c)
        if index != -1:
            compiler = compiler[:index]
    return compiler


def convertSecondsToReadableTime(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)
