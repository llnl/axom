.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _portability-label:

Host/device portability tiers
==============================

Slam's types are designed to run both on the host and inside GPU kernels. 

Not every C++ facility is safe in device code across all of Slam's backends
(sequential, OpenMP, CUDA, HIP), so Slam categorizes the constructs it uses into
three portability tiers. Compile-time features generate no device code 
and are therefore unconditionally kernel-safe.

.. list-table:: Slam's portability tiers
   :widths: 8 62 30
   :header-rows: 1

   * - Tier
     - Contents
     - Allowed where
   * - A
     - All compile-time language features (concepts, ``if constexpr``,
       class template argument deduction (CTAD), non-type template parameters,
       fold expressions, type traits, ``constexpr`` evaluation); 
       Axom host-device types (``StackArray``, ``ArrayView``, ``NumericLimits``, ``utilities::*``); 
       and Slam's own host-device types, including ``slam::Optional``.
     - everywhere, including kernels
   * - B
     - ``std::optional``, ``std::string_view``, ``std::variant``,
       ``std::ranges`` views/algorithms, ``std::vector``, ``std::map``,
       exceptions, and iostreams.
     - host only: builders, registries, ``isValid(verbose)``, and I/O
   * - C
     - Virtual functions on host-device types; standard containers held inside
       view types; throwing accessors on host-device paths.
     - nowhere (existing instances are migration targets)

Why not a device ``std::optional``?
-----------------------------------

``libcu++`` and ``libhipcxx`` provide device-capable ``tuple``/``optional``/ ``variant``/``span``,
but we cannot depend on them for our non-GPU sequential and OpenMP builds.
Instead, Slam uses small internal host-device types in the spirit of ``axom::StackArray``.

:cpp:class:`axom::slam::Optional` is the one such type Slam adds for this purpose. 
It is a trivially-copyable aggregate of an engaged flag and storage,
is ``AXOM_HOST_DEVICE`` throughout, and has no throwing ``value()`` accessor.
The contract is to check ``has_value()`` (or use ``value_or``) before dereferencing. 

The host-side counterpart is unchanged: host-only registries and find APIs
(for example :cpp:class:`axom::slam::FieldRegistry`) return ``std::optional``.

