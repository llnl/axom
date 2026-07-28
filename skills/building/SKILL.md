---
name: building
description: Instructions for building Axom
---

## Compute-node requirement (Codex)

Only run compilation and test commands on a compute node. Determine this by running:

```bash
./skills/building/scripts/is_compute_node
```

If it prints `login`, do **not** run `./config-build.py`, `cmake --build ...`, or `ctest ...`. Stop and ask the user to switch to/allocate a compute node, then continue once `is_compute_node` prints `compute`.

## Axom build

This repository builds Axom against externally provided TPLs (e.g., via Spack/uberenv + a host-config).

1) Configure (recommended wrapper around CMake):

```bash
./config-build.py -bp build -ip install -hc "$(./skills/building/scripts/determine_host_config)" --exportcompilercommands
```

2) Build and test:

```bash
cmake --build build -j
ctest --test-dir build
```

## Common build options

Common CMake options (and their defaults) live in `src/cmake/AxomOptions.cmake` and `src/cmake/CMakeBasics.cmake`. Pass them at configure time as `-D<OPTION>=ON|OFF`, for example:

```bash
./config-build.py -hc "$(./skills/building/scripts/determine_host_config)" -DAXOM_ENABLE_ASAN=ON
```

### Build types (when requested)

Only set a build type when the user explicitly asks for one. Use `config-build.py`'s `-bt` option, which takes a CMake build type (e.g., `Debug`, `Release`, `RelWithDebInfo`, `MinSizeRel`):

```bash
./config-build.py -bp build -ip install -hc "$(./skills/building/scripts/determine_host_config)" -bt <CMAKE_BUILD_TYPE> --exportcompilercommands
```

### AddressSanitizer (`AXOM_ENABLE_ASAN`)

AddressSanitizer is available via the `AXOM_ENABLE_ASAN` CMake option (default: `OFF`). It is supported with GCC or Clang.
Use it with the Debug CMake build type unless the user requests a different build type.

```bash
./config-build.py -hc "$(./skills/building/scripts/determine_host_config)" -bt Debug -DAXOM_ENABLE_ASAN=ON
```

### Host-config selection (`-hc`)

If the user does **not** specify a host-config file, determine the best match by running:

```bash
./skills/building/scripts/determine_host_config
```

Use its output as the `-hc` argument (as shown in the examples above). If the user explicitly provides a host-config file/path, use that instead.
