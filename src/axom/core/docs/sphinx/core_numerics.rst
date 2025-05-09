.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core numerics
******************************************************

The ``axom::numerics`` namespace was designed for convenient representation
and use of matrices and vectors, with accompanying manipulation and solver 
routines.

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

