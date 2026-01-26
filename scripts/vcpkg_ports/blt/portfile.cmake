vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/blt
    REF ccd6d31c7072ff984c7c8b723fc51f69e677ca72
    SHA512 dc3ab1b01293772a048d1531f6953e82cb9b38a162153c2002f89f1dfd07c0cf087e40153024ceeee2f5040ae2611ed8936292ce0bec568b068313ec2dfd853c
    HEAD_REF develop
)

file(MAKE_DIRECTORY
    ${CURRENT_PACKAGES_DIR}/include/
    ${CURRENT_PACKAGES_DIR}/share/blt
)

# Add an empty file to include to make vcpkg happy
file(WRITE ${CURRENT_PACKAGES_DIR}/include/empty.txt "")

# Copy files over to blt's 'share' directory
set(_src_dir ${SOURCE_PATH})
set(_dest_dir ${CURRENT_PACKAGES_DIR}/share/blt)

foreach(_f cmake LICENSE RELEASE-NOTES.md SetupBLT.cmake thirdparty_builtin tests )
  file(COPY ${_src_dir}/${_f} DESTINATION ${_dest_dir} )
endforeach()

# Rename the LICENSE file to 'copyright'
file(RENAME ${_dest_dir}/LICENSE ${_dest_dir}/copyright)
