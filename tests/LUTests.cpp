#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/matrix/LU.h"
#include "linear/matrix/Matrix.h"

TEST(LUTests, SimpleLU) {
  Matrix<Rational> A = {{1, 2}, {3, 4}};

  auto [L, U] = linalg::get_lu(A);

  ASSERT_EQ(L, (Matrix<Rational>{{1, 0}, {3, 1}}));
  ASSERT_EQ(U, (Matrix<Rational>{{1, 2}, {0, -2}}));
}

TEST(LUTests, SimpleInPlaceLU) {
  Matrix<Rational> A = {{1, 2}, {3, 4}};

  linalg::inplace_lu(A);

  ASSERT_EQ(A, (Matrix<Rational>{{1, 2}, {3, -2}}));
}

TEST(LUTests, BigLU) {
  const size_t N = 50;

  Matrix<Rational> A(N, N);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      A[i, j] = static_cast<long long>(i + j) % 71 + 123;
    }
  }

  A = A + Matrix<Rational>::unity(N);

  auto [L, U] = linalg::get_lu(A);

  ASSERT_EQ(L * U, A);
}

TEST(LUTests, BigInPlaceLU) {
  const size_t N = 50;

  Matrix<Rational> A(N, N);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      A[i, j] = static_cast<long long>(i + j) % 71 + 123;
    }
  }

  A = A + Matrix<Rational>::unity(N);

  auto A_copy = A;
  linalg::inplace_lu(A);

  Matrix<Rational> L(N, N, 0);
  Matrix<Rational> U(N, N, 0);

  for (size_t i = 0; i < N; ++i) {
    L[i, {0, i}] = A[i, {0, i}];
    L[i, i] = 1;

    U[i, {i, N}] = A[i, {i, N}];
  }

  ASSERT_EQ(L * U, A_copy);
}
