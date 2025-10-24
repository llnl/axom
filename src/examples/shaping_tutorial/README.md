[comment]: # (#################################################################)
[comment]: # (Copyright 2017-2025, Lawrence Livermore National Security, LLC)
[comment]: # (and Axom Project Developers. See the top-level LICENSE file)
[comment]: # (for details.)
[comment]: #
[comment]: # (# SPDX-License-Identifier: BSD-3-Clause)
[comment]: # (#################################################################)

# <img src="/share/axom/logo/axom_logo_transparent.png?raw=true" width="250" valign="middle" alt="Axom"/>

# Axom shaping tutorial
This tutorial introduces features of **Axom**, an open source Computer Science infrastructure framework designed to facilitate the development of multi-physics applications and computational tools.

Our focus application in this tutorial is a **"Shaping"** application, which pertains to volume fraction initialization in multimaterial simulations. This use case is representative of Axom’s broad capabilities, as it spans several key components and subsystems within Axom.

## Tutorial roadmap

This tutorial showcases several features of Axom by incrementally building up an application that performs shaping on a computational mesh.

* [Lesson 00](lesson_00/README.md) provides a brief overview of Axom and develops a simple application against an installed version of Axom.
* [Lesson 01](lesson_01/README.md) uses Axom's `Sidre` component to structure Cartesian mesh metadata (bounding box and resolution) and generate a Conduit Mesh Blueprint representation.
* [Lesson 02](lesson_02/README.md) focuses on Axom's `Inlet` component to define, parse, and validate simulation input for mesh metadata.
* [Lesson 03](lesson_03/README.md) describes geometry setup for multimaterial simulation using Axom's `Klee` component.
* [Lesson 04](lesson_04/README.md) "shapes" Klee-based geometry onto a computational mesh to compute per-material volume fractions using Axom's `Quest` component.
