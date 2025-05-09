// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/slic.hpp"
#include "axom/core/NumericArray.hpp"

#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Point.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"

#include <array>

namespace primal = axom::primal;

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_vector_policy()
{
  const int DIM = 3;
  using VectorType = primal::Vector<double, DIM>;

  VectorType* vec = axom::allocate<VectorType>(1, axom::execution_space<ExecSpace>::allocatorID());

  VectorType vec_host;

  axom::for_all<ExecSpace>(
    1,
    AXOM_LAMBDA(int /*i*/) {
      vec[0] = VectorType(-1.0);
      vec[0].negate();
    });

  axom::copy(&vec_host, vec, sizeof(VectorType));

  EXPECT_EQ(vec_host, VectorType(1));
  axom::deallocate(vec);
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_constructors)
{
  constexpr int DIM = 5;
  using CoordType = double;
  using QArray = axom::NumericArray<CoordType, DIM>;
  using QVec = primal::Vector<CoordType, DIM>;

  QVec vec1;
  EXPECT_EQ(vec1.dimension(), DIM);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec1[i], CoordType());
  }

  CoordType val = 5.;
  QVec vec2(val);
  EXPECT_EQ(vec2.dimension(), DIM);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec2[i], val);
  }

  CoordType val3 = 5.;
  int halfPt = DIM / 2;
  QVec vec3(val3, halfPt);
  for(int i = 0; i < DIM; ++i)
  {
    CoordType expVal = i < halfPt ? val3 : CoordType();
    EXPECT_EQ(vec3[i], expVal);
  }

  // Compare vector initialized from arbitrary
  CoordType valsArr[DIM] = {12., 23., 34., 45., 56.432};
  QVec arrVec(valsArr);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(arrVec[i], valsArr[i]);
  }

  QVec arrHalfVec(valsArr, halfPt);
  for(int i = 0; i < DIM; ++i)
  {
    CoordType expVal = i < halfPt ? valsArr[i] : CoordType();
    EXPECT_EQ(arrHalfVec[i], expVal);
  }

  // Construct from numeric array
  QArray arr(valsArr);
  QVec vFromNA(arr);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vFromNA[i], valsArr[i]);
  }

  // Initializer list constructor
  primal::Vector<int, 3> fromInitializerListRightSize = {10, 20, 30};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(10 * (i + 1), fromInitializerListRightSize[i]);
  }

  primal::Vector<int, 3> fromInitializerListTooLong = {10, 20, 30, 40};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(10 * (i + 1), fromInitializerListTooLong[i]);
  }

  primal::Vector<int, 5> fromInitializerListTooShort = {10, 20};
  for(int i = 0; i < 2; ++i)
  {
    EXPECT_EQ(10 * (i + 1), fromInitializerListTooShort[i]);
  }
  for(int i = 2; i < 5; ++i)
  {
    EXPECT_EQ(0, fromInitializerListTooShort[i]);
  }

  primal::Vector<int, 3> fromInitializerNoEqualsSign {10, 20, 30};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(10 * (i + 1), fromInitializerNoEqualsSign[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_from_points_constructor)
{
  constexpr int DIM = 5;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVec = primal::Vector<CoordType, DIM>;

  const double aVal = 12;
  const double bVal = 5;
  QPoint a(aVal);
  QPoint b(bVal);

  QVec vec1(a, b);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec1[i], bVal - aVal);
  }

  QVec vec2(b, a);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec2[i], aVal - bVal);
  }

  QVec vec3(b);
  QVec vec4(QPoint::zero(), b);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec3[i], bVal);
    EXPECT_EQ(vec4[i], bVal);
  }
  EXPECT_EQ(vec3, vec4);
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_normalize)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVec = primal::Vector<CoordType, DIM>;

  EXPECT_DOUBLE_EQ(QVec().unitVector().norm(), 1.0);

  // Find the norm of an arbitrary point
  QPoint p {543.5435, 1566.423532, -432.4};
  QVec vec2(QPoint(), p);
  EXPECT_DOUBLE_EQ(vec2.unitVector().norm(), 1.0);

  // Find the norm of the zero vector
  // Zero vector should become the unit vector when normalized
  QVec vecZero;
  QVec unitVec(1, 1);
  EXPECT_DOUBLE_EQ(vecZero.unitVector().norm(), 1.0);
  EXPECT_EQ(vecZero.unitVector(), unitVec);
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_norm)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVec = primal::Vector<CoordType, DIM>;

  QPoint p1 {3, 0};
  QPoint p2 {0, 4};
  QVec vec(p1, p2);
  EXPECT_DOUBLE_EQ(vec.norm(), 5.0);
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_arithmetic)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVec = primal::Vector<CoordType, DIM>;

  QPoint p1 {3, 0, 1.2};
  QPoint p2 {0, 4, 1.2};
  QPoint pSum {3, 4, 2.4};
  QPoint pDiff {-3, 4, 0};

  CoordType scalar = 5.3;
  QPoint pScalar {scalar * 3, scalar * 4, scalar * 2.4};

  QVec v1(p1);
  QVec v2(p2);

  // testing vector addition and subtraction
  EXPECT_EQ(QVec(p1, p2), v2 - v1);
  EXPECT_EQ(QVec(pDiff), v2 - v1);
  EXPECT_EQ(QVec(pSum), v1 + v2);
  EXPECT_EQ(v1 + v2, v2 + v1);

  QVec vSum(v1);
  vSum += v2;
  EXPECT_EQ(vSum, QVec(pSum));

  QVec vDiff(v2);
  vDiff -= v1;
  EXPECT_EQ(vDiff, QVec(pDiff));

  // Testing scalar multiplication
  EXPECT_EQ(QVec(pScalar), QVec(pSum) * scalar);
  EXPECT_EQ(QVec(pScalar), scalar * QVec(pSum));
  EXPECT_EQ(QVec(pScalar), (v1 + v2) * scalar);

  // Testing unary negation
  EXPECT_EQ(-v1, QVec(p1, QPoint::zero()));
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_inner_product)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVec = primal::Vector<CoordType, DIM>;

  QPoint p1 {3, 0, 1.2};
  QPoint p2 {0, 4, 1.2};

  CoordType expDot = 1.2 * 1.2;

  QVec v1(p1);
  QVec v2(p2);
  EXPECT_EQ(expDot, v1.dot(v2));
  EXPECT_EQ(v2.dot(v1), v1.dot(v2));
  EXPECT_EQ(v1.dot(v2), QVec::dot_product(v1, v2));
  EXPECT_EQ(QVec::dot_product(v1, v2), QVec::dot_product(v2, v1));

  QVec zer;
  EXPECT_EQ(0., zer.dot(zer));
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector3_outer_product)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QVec = primal::Vector<CoordType, DIM>;

  {
    QVec v1 {3, 0};
    QVec v2 {0, 4};

    QVec vRes {0, 0, 12};

    EXPECT_EQ(QVec::cross_product(v1, v2), vRes);
    EXPECT_EQ(QVec::cross_product(v2, v1), -vRes);
    EXPECT_EQ(QVec::cross_product(v1, v2), -QVec::cross_product(v2, v1));
  }

  {
    QVec v1 {1, 1, 1};
    QVec v2 {1, 1, 1};

    QVec vRes {0, 0, 0};

    EXPECT_EQ(QVec::cross_product(v1, v2), vRes);
    EXPECT_EQ(QVec::cross_product(v2, v1), -vRes);
    EXPECT_EQ(QVec::cross_product(v1, v2), -QVec::cross_product(v2, v1));
  }

  {
    QVec v1 {1, -1, 1};
    QVec v2 {-.125, .25, -.5};

    QVec vRes {0.25, 0.375, 0.125};

    EXPECT_EQ(QVec::cross_product(v1, v2), vRes);
    EXPECT_EQ(QVec::cross_product(v2, v1), -vRes);
    EXPECT_EQ(QVec::cross_product(v1, v2), -QVec::cross_product(v2, v1));
  }
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector2_outer_product)
{
  using CoordType = double;
  using QVec2 = primal::Vector<CoordType, 2>;
  using QVec3 = primal::Vector<CoordType, 3>;

  {
    QVec2 v1 {3, 0};
    QVec2 v2 {0, 4};

    QVec3 vRes {0, 0, 12};

    EXPECT_EQ(QVec2::cross_product(v1, v2), vRes);
    EXPECT_EQ(QVec2::cross_product(v2, v1), -vRes);
    EXPECT_EQ(QVec2::cross_product(v1, v2), -QVec2::cross_product(v2, v1));
  }

  {
    QVec2 v1 {1, -1};
    QVec2 v2 {.25, 4};

    QVec3 vRes {0, 0, 4.25};

    EXPECT_EQ(QVec2::cross_product(v1, v2), vRes);
    EXPECT_EQ(QVec2::cross_product(v2, v1), -vRes);
    EXPECT_EQ(QVec2::cross_product(v1, v2), -QVec2::cross_product(v2, v1));
  }

  {
    QVec2 v1 {1, 1};
    QVec2 v2 {-1, -1};

    QVec3 vRes {0, 0, 0};

    EXPECT_EQ(QVec2::cross_product(v1, v2), vRes);
    EXPECT_EQ(QVec2::cross_product(v2, v1), -vRes);
    EXPECT_EQ(QVec2::cross_product(v1, v2), -QVec2::cross_product(v2, v1));
  }
}

//------------------------------------------------------------------------------
TEST(primal_vector, scalar_triple_product)
{
  using CoordType = double;
  using QVec = primal::Vector<CoordType, 3>;
  using QVecTriple = std::array<QVec, 3>;

  QVec o {0, 0, 0};
  QVec e_0 {1, 0, 0};
  QVec e_1 {0, 1, 0};
  QVec e_2 {0, 0, 1};

  // explicitly test some examples
  EXPECT_EQ(1., QVec::scalar_triple_product(e_0, e_1, e_2));
  EXPECT_EQ(-1., QVec::scalar_triple_product(e_0, e_2, e_1));
  EXPECT_DOUBLE_EQ(-10, QVec::scalar_triple_product(-e_0, 2.5 * e_1, 4. * e_2));

  // test properties for a few examples
  std::vector<QVecTriple> data = {QVecTriple {e_0, e_1, e_2},
                                  QVecTriple {e_0, e_2, e_1},
                                  QVecTriple {e_0, e_0, e_1},
                                  QVecTriple {o, e_0, e_1},
                                  QVecTriple {e_0, e_0 + e_1, e_0 + e_1 + e_2},
                                  QVecTriple {-3.33 * e_0, 2.25 * e_1, .1111 * e_2}};

  for(auto triplet : data)
  {
    const auto& i = triplet[0];
    const auto& j = triplet[1];
    const auto& k = triplet[2];

    // check definition
    EXPECT_EQ(QVec::dot_product(i, QVec::cross_product(j, k)), QVec::scalar_triple_product(i, j, k));

    // check rotations
    EXPECT_EQ(QVec::scalar_triple_product(i, j, k), QVec::scalar_triple_product(j, k, i));

    EXPECT_EQ(QVec::scalar_triple_product(i, j, k), QVec::scalar_triple_product(k, i, j));

    // check swap permutation
    EXPECT_EQ(-QVec::scalar_triple_product(i, j, k), QVec::scalar_triple_product(i, k, j));
  }
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_zero)
{
  constexpr int DIM = 5;
  using QVec = primal::Vector<double, DIM>;

  QVec zero {0.0};
  EXPECT_TRUE(zero.is_zero());

  for(int i = 0; i < DIM; ++i)
  {
    QVec notZero = zero;
    notZero[i] = 1e-7;
    EXPECT_FALSE(notZero.is_zero()) << "Wrong when changing index " << i;
  }
}

//------------------------------------------------------------------------------
TEST(primal_vector, add_point_vector)
{
  constexpr int DIM = 3;
  using QPt = primal::Point<double, DIM>;
  using QVec = primal::Vector<double, DIM>;

  {
    QPt one {1.};
    QVec zero {0.0};
    EXPECT_EQ(one, one + zero);
    EXPECT_EQ(one, zero + one);
    EXPECT_EQ(zero + one, one + zero);
  }

  {
    QPt pt {1.23, 4.56, 7.89};
    QVec vec {.23, .56, .89};

    EXPECT_EQ(pt + vec, vec + pt);

    QPt pv = pt + -vec;
    QPt vp = -vec + pt;
    QPt exp {1, 4, 7};
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(pv[i], vp[i]);
      EXPECT_DOUBLE_EQ(exp[i], pv[i]);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_vector, subtract_points)
{
  constexpr int DIM = 3;
  using QPt = primal::Point<double, DIM>;
  using QVec = primal::Vector<double, DIM>;

  {
    QPt zeroPt {0.0};
    QPt onePt {1.0};

    QVec oneVec {1.0};

    EXPECT_EQ(oneVec, onePt - zeroPt);
    EXPECT_EQ(-oneVec, zeroPt - onePt);

    EXPECT_EQ(onePt - oneVec, zeroPt);
  }

  {
    QPt head {1.23, 4.56, 7.89};
    QPt tail {.23, .56, .89};

    QVec diff = head - tail;
    QVec exp {1, 4, 7};
    QVec ctor(tail, head);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(exp[i], diff[i]);
      EXPECT_DOUBLE_EQ(ctor[i], diff[i]);
      EXPECT_DOUBLE_EQ((head - diff)[i], tail[i]);
    }
  }
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(primal_numeric_array, numeric_array_check_policies)
{
  using seq_exec = axom::SEQ_EXEC;
  check_vector_policy<seq_exec>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)

  using omp_exec = axom::OMP_EXEC;
  check_vector_policy<omp_exec>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  check_vector_policy<cuda_exec>();
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  check_vector_policy<hip_exec>();
#endif
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
