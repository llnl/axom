.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core numerics
******************************************************

The ``axom::numerics`` namespace provides numerical algorithms including:

* **Vectors and Matrices:** Convenient representation with manipulation and solver routines
* **Numerical Quadrature:** Gauss-Legendre and rational Fejer integration rules

Numerical Quadrature
==============================

Axom provides 1D quadrature rules for numerical integration.
Two families are supported:

**Gauss-Legendre Quadrature**
  Standard polynomial-based rules for smooth functions. An N-point Gauss-Legendre
  rule can exactly integrate polynomials of degree 2N-1.

**Rational Fejer Quadrature**
  Pole-adapted rules for functions with known singularities or sharp features.
  By specifying pole locations (where the function or its approximation is singular),
  the quadrature adapts its basis to achieve high accuracy with fewer points.

Both families use the same ``QuadratureRule`` interface. In the examples below,
we obtain a rule with ``get_gauss_legendre()`` or ``get_rational_fejer()``,
then either inspect its nodes and weights directly or pass it to a small helper
that evaluates the weighted sum:

.. literalinclude:: ../../examples/core_quadrature.cpp
   :start-after: _quadrature_rule_start
   :end-before: _quadrature_rule_end
   :language: C++

Both ``get_*`` entry points return a lightweight non-owning
``QuadratureRuleView`` backed by cached quadrature data.

Gauss-Legendre Quadrature
------------------------------

Gauss-Legendre rules are optimal for integrating smooth, polynomial-like
functions on [0, 1].

We can access a Gauss-Legendre rule with ``n`` nodes as follows.

.. code-block:: cpp

  auto rule = axom::numerics::get_gauss_legendre(5);

This returns a cached ``QuadratureRuleView``. On the first request for a given
order, Axom computes the nodes and weights with
``axom::numerics::compute_gauss_legendre_data()`` and stores them for reuse.

Gauss-Legendre rules have polynomial exactness. A 5-point rule integrates
polynomials up to degree 9 exactly:

.. literalinclude:: ../../examples/core_quadrature.cpp
   :start-after: _gauss_legendre_exactness_start
   :end-before: _gauss_legendre_exactness_end
   :language: C++

For smooth functions, higher-order rules converge rapidly:

.. literalinclude:: ../../examples/core_quadrature.cpp
   :start-after: _gauss_legendre_smooth_start
   :end-before: _gauss_legendre_smooth_end
   :language: C++

Example output:

.. code-block:: text

   Integrate sin(pi*x) from 0 to 1 with varying orders:
      3 points: error = 4.421e-04
      5 points: error = 3.510e-08
     10 points: error = 1.110e-16
     20 points: error = 2.220e-16



Rational Fejer Quadrature
------------------------------

Rational Fejer quadrature excels at integrating functions with known singularities
or sharp gradients at specific locations called **poles**. A pole is a complex
value where the function (or its rational approximation) becomes singular.

For example, integrating f(x) = 1/(x - 1.5) on [0, 1] with a pole at x = 1.5:

.. literalinclude:: ../../examples/core_quadrature.cpp
   :start-after: _rational_fejer_singularity_start
   :end-before: _rational_fejer_singularity_end
   :language: C++

Comparison with Gauss-Legendre
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For functions with near-singularities, rational Fejer can achieve machine
precision with far fewer points than Gauss-Legendre:

.. literalinclude:: ../../examples/core_quadrature.cpp
   :start-after: _rational_fejer_vs_gauss_start
   :end-before: _rational_fejer_vs_gauss_end
   :language: C++


Example output:

.. code-block:: text

   Comparison: Gauss-Legendre vs Rational Fejer
   Function: f(x) = 1/(x - 1.2) on [0, 1]
   Singularity at x = 1.2 (close to domain boundary)
   
   Gauss-Legendre:
      5 points: error = 4.241e-04
     10 points: error = 7.516e-08
     20 points: error = 8.882e-16
     50 points: error = 8.882e-16
   
   Rational Fejer (pole at x = 1.2):
      2 points: error = 0.000e+00
     (Achieves machine precision with only 2 points!)

Multiple Singularities
^^^^^^^^^^^^^^^^^^^^^^^

Functions with multiple singularities can specify multiple poles:

.. literalinclude:: ../../examples/core_quadrature.cpp
   :start-after: _rational_fejer_multiple_poles_start
   :end-before: _rational_fejer_multiple_poles_end
   :language: C++

Pole Selection Guidelines
^^^^^^^^^^^^^^^^^^^^^^^^^^

**When to use rational Fejer:**

* Functions with known singularities outside [0, 1]
* Sharp gradients near domain boundaries (geometric corners)
* Standard Gauss-Legendre requires very high orders
* Integrand is well-approximated by rational functions

**Pole requirements:**

* Finite real poles must lie **outside** [0, 1] (e.g., -0.2, 1.5, 2.0)
* Complex poles are auto-completed to conjugate pairs
* Infinite poles: use
  ``std::complex<double>(std::numeric_limits<double>::infinity(), 0.0)``
  for the polynomial limit
* The rule has (number of canonical poles + 1) points

**Pole placement strategies:**

1. **Known singularities:** Place poles at actual singularity locations
2. **Geometric features:** Place poles near sharp corners or discontinuities
3. **Multiple poles:** Repeat a pole to model higher-order singularities
4. **Infinite poles:** All poles at infinity gives polynomial Fejer (Chebyshev-based)
5. **Mixed finite/infinite poles:** Endpoint algebraic singularities such as
   ``sqrt(x)`` can benefit from combining nearby real poles with a few poles at
   infinity, rather than using only nearby real poles

Advanced Usage
^^^^^^^^^^^^^^

QuadratureRule vs QuadratureRuleView
"""""""""""""""""""""""""""""""""""""

Similar to Array/ArrayView, quadrature rules have owning and non-owning versions:

.. literalinclude:: ../../examples/core_quadrature.cpp
   :start-after: _quadrature_rule_copy_start
   :end-before: _quadrature_rule_copy_end
   :language: C++

Cached vs Direct Computation
"""""""""""""""""""""""""""""

Two APIs are provided for rational Fejer:

.. literalinclude:: ../../examples/core_quadrature.cpp
   :start-after: _quadrature_compute_vs_cached_start
   :end-before: _quadrature_compute_vs_cached_end
   :language: C++

The ``get_rational_fejer()`` function uses an LRU cache of 65,536 rules. When
the cache is full, the least recently used rule is evicted, invalidating views to
it. Use ``.copy()`` if you need stable storage beyond immediate use.

Infinite Poles
""""""""""""""

All poles at infinity produce the standard polynomial Fejer
(Chebyshev-based) rule:

.. literalinclude:: ../../examples/core_quadrature.cpp
   :start-after: _quadrature_infinite_poles_start
   :end-before: _quadrature_infinite_poles_end
   :language: C++

Matrices and Vectors
==============================

The following sections describe matrix and vector operations.

The following example shows some basic vector operations.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _vecops_start
   :end-before: _vecops_end
   :language: C++

When run, the example produces the following output::

  Originally, u and v are
  u = [4, 1, 0]
  v = [1, 2, 3]

  The dot product is 6 and the cross product is
  [3, -12, 7]

  Now orthogonal u and normalized v are
  u = [3.57143, 0.142857, -1.28571]
  v = [0.267261, 0.534522, 0.801784]

  Root-finding returned 0 (should be 0, success).  Found 3 roots (should be 3)
  at x = 1.5, 1, -2 (should be x = -2, 1, 1.5 in arbitrary order).


The following example code shows how to construct a matrix.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _matctor_start
   :end-before: _matctor_end
   :language: C++

We can add and multiply matrices, vectors, and scalars, find the determinant,
and extract upper and lower triangular matrices as is shown in the next 
example.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _matops_start
   :end-before: _matops_end
   :language: C++

The example generates the following output::

  Originally, the matrix A = 
  [ 0.6 2.4 1.1 ]
  [ 2.4 0.6 -0.1 ]
  [ 1.1 -0.1 0.6 ]

  A + identity matrix = 
  [ 1.6 2.4 1.1 ]
  [ 2.4 1.6 -0.1 ]
  [ 1.1 -0.1 1.6 ]

  A * 2*(identity matrix) = 
  [ 1.2 4.8 2.2 ]
  [ 4.8 1.2 -0.2 ]
  [ 2.2 -0.2 1.2 ]

  Vector x1 = [1, 2, -0.5]
  A * x1 = [4.85, 3.65, 0.6]

  Determinant of A = -4.5
  A's lower triangle = 
  [ 0.6 0 0 ]
  [ 2.4 0.6 0 ]
  [ 1.1 -0.1 0.6 ]

  A's upper triangle (with 1s in the main diagonal) = 
  [ 1 2.4 1.1 ]
  [ 0 1 -0.1 ]
  [ 0 0 1 ]

  A's column 1 is [2.4, 0.6, -0.1]


We can also extract rows and columns.  The preceding example shows how to get
a column.  Since the underlying storage layout of Matrix is column-based, 
retrieving a row is a little more involved: the call to `getRow()` retrieves 
the stride for accessing row elements `p` as well the upper bound for element 
indexes in the row. The next selection shows how to sum the entries in a row.

.. literalinclude:: ../../numerics/internal/matrix_norms.hpp
   :start-after: _rowsum_start
   :end-before: _rowsum_end
   :language: C++

We can use the power method or the Jacobi method to find the eigenvalues and
vectors of a matrix.  The power method is a stochastic algorithm, computing many
matrix-vector multiplications to produce approximations of a matrix's
eigenvalues and vectors.  The Jacobi method is also an iterative algorithm, but
it is not stochastic, and tends to converge much more quickly and stably than
other methods.  However, the Jacobi method is only applicable to symmetric
matrices.  In the following snippet, we show both the power method and the
Jacobi method to show that they get the same answer.

.. note::
   As of August 2020, the API of ``eigen_solve`` is not consistent
   with ``jacobi_eigensolve`` (``eigen_solve`` takes a ``double`` pointer as 
   input instead of a ``Matrix`` and the return codes differ).  This is an 
   issue we're fixing.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _eigs_start
   :end-before: _eigs_end
   :language: C++

Here is the output of the code example::

  Tried to find 3 eigenvectors and values for matrix 
  [ 0.6 2.4 1.1 ]
  [ 2.4 0.6 -0.1 ]
  [ 1.1 -0.1 0.6 ]

  and the result code was 1 (1 = success).

  Eigenvalue 0 = 3.2033 Eigenvector 0 = [0.711931, 0.645731, 0.276015]
  Eigenvalue 1 = -2.07901 Eigenvector 1 = [0.701812, -0.64037, -0.312067]
  Eigenvalue 2 = 0.675707 Eigenvector 2 = [0.0247596, -0.415881, 0.909082]

  Using the Jacobi method, tried to find eigenvectors and eigenvalues of matrix 
  [ 0.6 2.4 1.1 ]
  [ 2.4 0.6 -0.1 ]
  [ 1.1 -0.1 0.6 ]

  and the result code was 0 (0 = success).

  Eigenvalue 0 = -2.07901 Eigenvector 0 = [0.701812, -0.64037, -0.312067]
  Eigenvalue 1 = 0.675707 Eigenvector 1 = [0.0247596, -0.415881, 0.909082]
  Eigenvalue 2 = 3.2033 Eigenvector 2 = [0.711931, 0.645731, 0.276015]


We can also solve a linear system directly or by using LU decomposition and
back-substitution, as shown in the next example.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _solve_start
   :end-before: _solve_end
   :language: C++

The example produces the following output::

  Solved for x in the linear system Ax = b, where
  A = 
  [ 3 2.66667 4.66667 ]
  [ 2 0.666667 5.5 ]
  [ 1 -0.666667 3 ]
   and b = [3, 13, 4].

  Result code is 0 (0 = success)
  Found x = [3, 4, -2]

  Decomposed to 
  [ 3 2.66667 4.66667 ]
  [ 2 0.666667 5.5 ]
  [ 1 -0.666667 3 ]
   with pivots [1, 2, 2] with result 0 (0 is success)
  Found x = [3, 4, -2] 
