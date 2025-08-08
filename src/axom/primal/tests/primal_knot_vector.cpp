// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_knot_vector.cpp
 * \brief This file tests primal's knot vector functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/KnotVector.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_knotvector, default_constructor)
{
  primal::KnotVector<double> kvector;

  constexpr axom::IndexType zero {0};

  EXPECT_EQ(-1, kvector.getDegree());
  EXPECT_EQ(zero, kvector.getNumControlPoints());
  EXPECT_EQ(zero, kvector.getNumKnots());
  EXPECT_EQ(zero, kvector.getNumKnotSpans());
  EXPECT_FALSE(kvector.isValid());
}

//------------------------------------------------------------------------------
TEST(primal_knotvector, degree_constructor)
{
  constexpr int degree = 2;
  constexpr int npts = 5;
  primal::KnotVector<double> kvector(npts, degree);
  EXPECT_TRUE(kvector.isValid());

  // Check size of knot vector
  EXPECT_EQ(degree, kvector.getDegree());
  EXPECT_EQ(npts, kvector.getNumControlPoints());
  EXPECT_EQ(npts + degree + 1, kvector.getNumKnots());

  // Check for monotonicity
  for(int i = 0; i < kvector.getNumKnots() - 1; ++i)
  {
    EXPECT_LE(kvector[i], kvector[i + 1]);
  }

  // Check for clamped-ness
  for(int i = 0; i <= degree; ++i)
  {
    EXPECT_DOUBLE_EQ(kvector[i], 0.0);
    EXPECT_DOUBLE_EQ(kvector[kvector.getNumKnots() - 1 - i], 1.0);
  }
}

//------------------------------------------------------------------------------
TEST(primal_knotvector, array_constructors)
{
  constexpr int degree = 2;

  constexpr int nkts = 9;
  double knots[] = {0.0, 0.0, 0.0, 0.2, 0.5, 0.8, 1.0, 1.0, 1.0};
  constexpr int npts = nkts - degree - 1;
  constexpr int nspans = 4;

  axom::Array<double> knot_arr;
  knot_arr.assign(std::begin(knots), std::end(knots));

  // lambda to check validity of KnotVector instances
  auto check_knot_vector = [=](const primal::KnotVector<double>& kv) {
    EXPECT_TRUE(kv.isValid());

    // Check sizeof knot vector
    EXPECT_EQ(degree, kv.getDegree());
    EXPECT_EQ(npts, kv.getNumControlPoints());
    EXPECT_EQ(npts + degree + 1, kv.getNumKnots());
    EXPECT_EQ(nspans, kv.getNumKnotSpans());

    // Check for array access
    for(int i = 0; i < degree; ++i)
    {
      EXPECT_DOUBLE_EQ(kv[i], knots[i]);
      EXPECT_DOUBLE_EQ(kv.getArray()[i], knots[i]);
    }
  };

  // test C-array constructor
  {
    primal::KnotVector<double> kvector(knots, nkts, degree);
    check_knot_vector(kvector);
  }

  // test ArrayView from C-array constructor
  {
    axom::ArrayView<double> knots_v(knots, nkts);
    primal::KnotVector<double> kvector(knots_v, degree);
    check_knot_vector(kvector);
  }

  // test Array constructor
  {
    primal::KnotVector<double> kvector(knot_arr, degree);
    check_knot_vector(kvector);
  }

  // test ArrayView from axom::Array
  {
    primal::KnotVector<double> kvector(knot_arr.view(), degree);
    check_knot_vector(kvector);
  }
}

//------------------------------------------------------------------------------
TEST(primal_knotvector, set_degree)
{
  primal::KnotVector<double> kvector;

  // Change the degree of the knot vector
  EXPECT_EQ(-1, kvector.getDegree());
  EXPECT_EQ(0, kvector.getNumControlPoints());
  EXPECT_EQ(0, kvector.getNumKnots());

  // Will increase the number of implied control points
  kvector.setDegree(4);
  EXPECT_EQ(4, kvector.getDegree());
  EXPECT_EQ(10, kvector.getNumKnots());

  // Will not change the number of implied control points
  kvector.setDegree(3);
  EXPECT_EQ(3, kvector.getDegree());
  EXPECT_EQ(9, kvector.getNumKnots());

  // Will not change the number of implied control points
  kvector.setDegree(4);
  EXPECT_EQ(4, kvector.getDegree());
  EXPECT_EQ(10, kvector.getNumKnots());

  // Will again increase the number of implied control points
  kvector.setDegree(5);
  EXPECT_EQ(5, kvector.getDegree());
  EXPECT_EQ(12, kvector.getNumKnots());

  // Reset everything
  kvector.clear();
  EXPECT_EQ(-1, kvector.getDegree());
  EXPECT_EQ(0, kvector.getNumControlPoints());
  EXPECT_EQ(0, kvector.getNumKnots());
}

//------------------------------------------------------------------------------
TEST(primal_knotvector, make_uniform)
{
  primal::KnotVector<double> kvector;

  const int degree = 2;
  const int npts = 10;

  // Make a uniform knot vector
  kvector.makeUniform(npts, degree);

  // Check size of knot vector
  EXPECT_EQ(degree, kvector.getDegree());
  EXPECT_EQ(npts, kvector.getNumControlPoints());
  EXPECT_EQ(npts + degree + 1, kvector.getNumKnots());

  // Check for clamped-ness
  for(int i = 0; i <= degree; ++i)
  {
    EXPECT_DOUBLE_EQ(kvector[i], 0.0);
    EXPECT_DOUBLE_EQ(kvector[kvector.getNumKnots() - 1 - i], 1.0);
  }

  // Check for uniformity
  for(int i = degree; i < npts - 1; ++i)
  {
    EXPECT_DOUBLE_EQ(kvector[i + 1] - kvector[i], kvector[i + 2] - kvector[i + 1]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_knotvector, find_knot_span)
{
  primal::KnotVector<double> kvector;

  const int degree = 3;

  // clang-format off
  axom::Array<double> knots {0.0, 0.0, 0.0, 0.0, 
                             0.2, 
                             0.3, 0.3, 
                             0.5, 
                             0.8, 0.8, 0.8, 
                             0.9, 
                             1.0, 1.0, 1.0, 1.0};
  // clang-format on

  kvector = primal::KnotVector<double>(knots, degree);

  // Check for the number of knot spans
  EXPECT_EQ(6, kvector.getNumKnotSpans());

  // Check for knot spans
  EXPECT_EQ(3, kvector.findSpan(0.0));
  EXPECT_EQ(3, kvector.findSpan(0.1));
  EXPECT_EQ(4, kvector.findSpan(0.2));
  EXPECT_EQ(6, kvector.findSpan(0.4));
  EXPECT_EQ(10, kvector.findSpan(0.85));
  EXPECT_EQ(11, kvector.findSpan(0.9));
  EXPECT_EQ(11, kvector.findSpan(1.0));

  // Check for knot spans with multiplicity
  int multiplicity;
  EXPECT_EQ(3, kvector.findSpan(0.0, multiplicity));
  EXPECT_EQ(multiplicity, 4);
  EXPECT_EQ(3, kvector.findSpan(0.1, multiplicity));
  EXPECT_EQ(multiplicity, 0);
  EXPECT_EQ(4, kvector.findSpan(0.2, multiplicity));
  EXPECT_EQ(multiplicity, 1);
  EXPECT_EQ(6, kvector.findSpan(0.3, multiplicity));
  EXPECT_EQ(multiplicity, 2);
  EXPECT_EQ(10, kvector.findSpan(0.8, multiplicity));
  EXPECT_EQ(multiplicity, 3);
  EXPECT_EQ(11, kvector.findSpan(1.0, multiplicity));
  EXPECT_EQ(multiplicity, 4);
}

//------------------------------------------------------------------------------
TEST(primal_knotvector, normalize)
{
  const int degree = 2;
  axom::Array<double> knots {-1, -1, -1, 0, 1, 2, 2, 4, 4, 4};

  primal::KnotVector<double> kvector(knots, degree);
  primal::KnotVector<double> kvector_norm(knots, degree);
  kvector_norm.normalize();

  // Check that both vectors are valid
  EXPECT_TRUE(kvector.isValid());
  EXPECT_TRUE(kvector_norm.isValid());

  // Check for the number of knot spans
  EXPECT_EQ(4, kvector.getNumKnotSpans());
  EXPECT_EQ(4, kvector_norm.getNumKnotSpans());

  // Check that the knot vectors are equivalent
  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    auto basis = kvector.calculateBasisFunctions(5.0 * t - 1.0);
    auto basis_norm = kvector_norm.calculateBasisFunctions(t);

    for(int i = 0; i < basis.size(); ++i)
    {
      EXPECT_NEAR(basis[i], basis_norm[i], 1e-10);
    }
  }

  // Check rescaling the knot vector
  kvector_norm.rescale(-1.0, 4.0);
  for(int i = 0; i < kvector.getNumKnots(); ++i)
  {
    EXPECT_DOUBLE_EQ(kvector[i], kvector_norm[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_knotvector, insert_knot)
{
  const int degree = 3;

  // clang-format off
  axom::Array<double> knots {0.0, 0.0, 0.0, 0.0, 
                             0.2, 
                             0.3, 0.3, 
                             0.5, 
                             0.8, 0.8, 0.8, 
                             0.9, 
                             1.0, 1.0, 1.0, 1.0};
  // clang-format on

  primal::KnotVector<double> kvector(knots, degree);

  // Insert knots into the vector to have specific multiplicities
  EXPECT_EQ(16, kvector.getNumKnots());
  EXPECT_EQ(6, kvector.getNumKnotSpans());

  kvector.insertKnot(0.0, 2);
  EXPECT_EQ(16, kvector.getNumKnots());
  EXPECT_EQ(6, kvector.getNumKnotSpans());

  kvector.insertKnot(0.1, 1);
  EXPECT_EQ(17, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  kvector.insertKnot(0.2, 2);
  EXPECT_EQ(18, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  kvector.insertKnot(0.3, 2);
  EXPECT_EQ(18, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  kvector.insertKnot(0.5, 3);
  EXPECT_EQ(20, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  kvector.insertKnot(0.8, 1);
  EXPECT_EQ(20, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  kvector.insertKnot(1.0, 1);
  EXPECT_EQ(20, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  // Already inserted knots shouldn't be removed
  kvector.insertKnot(0.5, 1);
  EXPECT_EQ(20, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  kvector.insertKnot(0.5, 2);
  EXPECT_EQ(20, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  kvector.insertKnot(0.5, 3);
  EXPECT_EQ(20, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  // Insert knots into the vector a given number of times
  kvector = primal::KnotVector<double>(knots, degree);

  EXPECT_EQ(16, kvector.getNumKnots());
  EXPECT_EQ(6, kvector.getNumKnotSpans());

  kvector.insertKnotBySpan(3, 0.1, 2);
  EXPECT_EQ(18, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  kvector.insertKnotBySpan(5, 0.2, 2);
  EXPECT_EQ(20, kvector.getNumKnots());
  EXPECT_EQ(7, kvector.getNumKnotSpans());

  kvector.insertKnotBySpan(8, 0.25, 3);
  EXPECT_EQ(23, kvector.getNumKnots());
  EXPECT_EQ(8, kvector.getNumKnotSpans());
}

//------------------------------------------------------------------------------
TEST(primal_knotvector, split)
{
  const int degree = 2;
  axom::Array<double> knots {0.0, 0.0, 0.0, 0.2, 0.5, 0.8, 0.8, 1.0, 1.0, 1.0};
  axom::Array<double> knots1 {0.0, 0.0, 0.0, 0.2, 0.5, 0.5, 0.5};
  axom::Array<double> knots2 {0.5, 0.5, 0.5, 0.8, 0.8, 0.8};
  axom::Array<double> knots3 {0.8, 0.8, 0.8, 1.0, 1.0, 1.0};

  primal::KnotVector<double> kvector(knots, degree);
  primal::KnotVector<double> kvector1, kvector2, kvector3;

  // Split the knot vector at an input value
  kvector.split(0.5, kvector1, kvector2);

  // Split the knot vector at a knot span
  kvector2.splitBySpan(4, kvector2, kvector3);

  // Check the knot vectors
  for(int i = 0; i < kvector1.getNumKnots(); ++i)
  {
    EXPECT_DOUBLE_EQ(kvector1[i], knots1[i]);
  }

  for(int i = 0; i < kvector2.getNumKnots(); ++i)
  {
    EXPECT_DOUBLE_EQ(kvector2[i], knots2[i]);
  }

  for(int i = 0; i < kvector3.getNumKnots(); ++i)
  {
    EXPECT_DOUBLE_EQ(kvector3[i], knots3[i]);
  }
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
