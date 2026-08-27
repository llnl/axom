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
       fold expressions, type traits, ``constexpr`` evaluation)
       as well as Axom host-device types (``StackArray``, ``ArrayView``, 
       ``NumericLimits``, ``utilities::*``). 
     - everywhere, including kernels
   * - B
     - ``std::string_view``, ``std::variant``,
       ``std::ranges`` views/algorithms, ``std::vector``, ``std::map``,
       exceptions, and iostreams.
     - host only: builders, registries, ``isValid(verbose)``, and I/O
   * - C
     - Virtual functions on host-device types; standard containers held inside
       view types; throwing accessors on host-device paths.
     - nowhere (existing instances are migration targets)

Device ``std::optional``
------------------------

Slam uses ``std::optional`` for optional-returning APIs.
The CUDA host-configs enable ``--expt-relaxed-constexpr``,
and HIP's compiler accepts these standard library calls from host-device paths in the supported builds.

The device-side contract is to check ``has_value()`` before dereferencing with ``operator*``,
or to use ``value_or``. Avoid ``value()`` in kernels because standard library implementations
can route the disengaged case through a throwing host-only helper.
