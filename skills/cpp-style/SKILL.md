---
name: axom-cpp-style
description: C++ coding style and naming conventions for Axom (primarily src/axom), including formatting via the CMake `style` target.
---

# Axom C++ Style & Naming

Use this skill when adding or editing C++ in `src/axom/`, `src/tools/`, and `src/examples/`.

## Formatting (clang-format)

- Axom uses repo-root `.clang-format` (Google-derived, `ColumnLimit: 120`, includes are not auto-sorted).
- Prefer letting formatting be enforced by the build-system target instead of hand-formatting.

Run auto-format (from an existing build directory):

```bash
cmake --build <build_dir> --target style
# example:
cmake --build build --target style
```

## File structure

- **License header**: Keep the existing copyright + SPDX block at the top of C++ files.
- **Header guards**: Use `#pragma once` for headers.
- **Doxygen file header**: Use a `/** ... */` block with `@file` and `@brief`.
- **Namespaces**: Prefer nested namespaces like `namespace axom::slam { ... }`.
- **Namespace closing comment**: Close with `}  // namespace axom::slam` (match the opened namespace).
- **Trailing newline**: Ensure each file ends with a newline.

## Includes

`.clang-format` has `SortIncludes: false`, so keep the project’s existing conventions:

- In a `.cpp`, include the corresponding header first (e.g., `#include "axom/foo.hpp"`).
- Include `axom/config.hpp` before any Axom headers when a source file needs build configuration macros.
- Group includes with blank lines (typical grouping: C++ standard library, third-party, Axom headers).
- Don’t churn includes solely to “sort” them; keep diffs minimal and consistent.
- Use `// clang-format off` / `// clang-format on` only for tightly controlled formatting blocks (e.g., initializer tables).

## Naming conventions (as used in `src/axom`)

- **Namespaces**: `axom` at the root; submodules use `axom::<module>`.
- **Types** (`class`, `struct`, `enum class`, `using` aliases): `PascalCase` (e.g., `GeometryOperator`).
- **Functions/methods**: `camelCase` (e.g., `defineAndParse`, `findMeshFilePath`).
- **Local variables / parameters**: `snake_case` (e.g., `input_file_path`, `restart_cycle`).
- **Data members**: `snake_case` preceded with `m_` (e.g., `m_foo`, `m_use_warm_start`).
- **Macros / compile options**: `SCREAMING_SNAKE_CASE` (e.g., `AXOM_ENABLE_MPI`, `AXOM_DEBUG_PARAM`).
- **Constants**: prefer `SCREAMING_SNAKE_CASE` when it’s part of a type’s API (e.g., `NUM_STATE_VARS`); for small local `constexpr` values, use a descriptive name that reads well at the callsite.

## Comments & documentation

- Prefer Doxygen-style docstrings for public APIs:
  - `@brief` for summary
  - `@param[in]` / `@param[in,out]` for parameters
  - `@tparam` for template parameters
- Use `///` for short Doxygen comments on declarations when a full block is overkill.
