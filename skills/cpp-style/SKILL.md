---
name: axom-cpp-style
description: C++ coding style and naming conventions for Axom (primarily src/axom), including formatting via the CMake `style` target.
---

# Axom C++ Style & Naming

Use this skill when adding or editing C++ in `src/axom/`, `src/tools/`, and `src/examples/`.

## Formatting (clang-format)

- Axom's CMake source root is `src/`, and the formatter config is `src/.clang-format`.
- The config is Google-derived, currently requires clang-format 19 from the build system, uses `ColumnLimit: 100`, and does not auto-sort includes.
- Prefer letting formatting be enforced by the build-system target instead of hand-formatting.

Run auto-format from an existing build directory configured from `src/`:

```bash
cmake --build <build_dir> --target style
# example:
cmake --build build-debug --target style
```

Useful narrower targets include `clangformat_style`, `clangformat_check`, and component variants such as `core_clangformat_style`.

## File structure

- **License header**: Keep the existing copyright + SPDX block at the top of C++ files.
- **Header guards**: Use `#pragma once` for headers.
- **Doxygen file header**: Use a Doxygen file prologue with `\file` and `\brief`.
- **Namespaces**: Prefer C++17 nested namespace syntax, e.g. `namespace axom::slam { ... }`.
- **Namespace closing comment**: Close with `}  // namespace axom::slam` (match the opened namespace).
- **Trailing newline**: Ensure each file ends with a newline.

## Includes

`src/.clang-format` has `SortIncludes: false`, so keep the project’s existing conventions:

- In a `.cpp`, include the corresponding header first (e.g., `#include "axom/foo.hpp"`).
- Include `axom/config.hpp` before any Axom headers when a source file needs build configuration macros.
- Group includes with blank lines in this order: same Axom component, other Axom/project headers, TPL headers, C/C++ standard library, then system headers.
- Don’t churn includes solely to “sort” them; keep diffs minimal and consistent.
- Use `// clang-format off` / `// clang-format on` only for tightly controlled formatting blocks (e.g., initializer tables).

## Naming conventions (as used in `src/axom`)

- When editing existing files, follow the naming conventions already used in that file for consistency.
- When creating new files under `src/axom/<component>/`, follow the naming style already used in that component directory.
- **Namespaces**: `axom` at the root; submodules use lowercase names such as `axom::slam`.
- **Types** (`class`, `struct`, `enum class`, `using` aliases): mixed/Pascal case (e.g., `GeometryOperator`).
- **Functions/methods**: `camelCase` (e.g., `defineAndParse`, `findMeshFilePath`).
- **Local variables / parameters**: `snake_case`/pot-hole (e.g., `input_file_path`, `restart_cycle`).
- **Non-static data members**: `m_` plus snake_case (e.g., `m_foo`, `m_use_warm_start`).
- **Static data members / file-scope statics**: `s_` plus snake_case (e.g., `s_default_allocator`).
- **Macros / compile options**: `SCREAMING_SNAKE_CASE` (e.g., `AXOM_ENABLE_MPI`, `AXOM_DEBUG_PARAM`).
- **Enumeration values**: Initial capital, with underscores between words when needed (e.g., `Num_Orange_Types`).
- **Constants**: `SCREAMING_SNAKE_CASE` for API constants; `snake_case` for local `constexpr` values.

## Comments & documentation

- Prefer Doxygen-style docstrings for public APIs:
  - `\brief` for summary
  - `\param[in]` / `\param[in,out]` for parameters
  - `\tparam` for template parameters
- Use `///` for short Doxygen comments on declarations when a full block is overkill.
