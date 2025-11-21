#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/matrix/LU.h"
#include "linear/matrix/Matrix.h"

Matrix<Rational> big_matrix(size_t N) {
  Matrix<Rational> A(N, N);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      A[i, j] = static_cast<long long>(i + j) % 71 + 123;
    }
  }

  A = A + Matrix<Rational>::unity(N);

  return A;
}

template <typename Field>
void check_U(const Matrix<Field>& matrix) {
  auto [n, d] = matrix.shape();

  ASSERT_EQ(n, d);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      ASSERT_EQ((matrix[i, j]), 0);
    }
  }
}

template <typename Field>
void check_L(const Matrix<Field>& matrix) {
  auto [n, d] = matrix.shape();

  ASSERT_EQ(n, d);

  for (size_t i = 0; i < n; ++i) {
    ASSERT_EQ((matrix[i, i]), 1);

    for (size_t j = i + 1; j < n; ++j) {
      ASSERT_EQ((matrix[i, j]), 0);
    }
  }
}

template <typename Field>
Matrix<Field> apply_permutation(const Matrix<Field>& matrix,
                                const std::vector<size_t>& permutation) {
  auto [n, d] = matrix.shape();
  Matrix<Field> result(n, d);

  for (size_t i = 0; i < n; ++i) {
    result[permutation[i], {0, d}] = matrix[i, {0, d}];
  }

  return result;
}

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

  auto A = big_matrix(N);
  auto [L, U] = linalg::get_lu(A);

  check_L(L);
  check_U(U);

  ASSERT_EQ(L * U, A);
}

TEST(LUTests, BigInPlaceLU) {
  const size_t N = 50;

  auto A = big_matrix(N);
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

TEST(LUTests, SimpleLUP) {
  Matrix<Rational> A = {{1, 2}, {3, 4}};

  auto [L, U, P] = linalg::get_lup(A);

  auto product = apply_permutation(L * U, P);
  ASSERT_EQ(product, A);
}

TEST(LUTests, BigLUP) {
  size_t N = 50;
  auto A = big_matrix(N);

  auto [L, U, P] = linalg::get_lup(A);
  check_L(L);
  check_U(U);

  auto product = apply_permutation(L * U, P);
  ASSERT_EQ(product, A);
}

TEST(LUTests, BigInPlaceLUP) {
  const size_t N = 50;

  auto A = big_matrix(N);
  auto A_copy = A;
  auto P = linalg::inplace_lup(A);

  Matrix<Rational> L(N, N, 0);
  Matrix<Rational> U(N, N, 0);

  for (size_t i = 0; i < N; ++i) {
    L[i, {0, i}] = A[i, {0, i}];
    L[i, i] = 1;

    U[i, {i, N}] = A[i, {i, N}];
  }

  ASSERT_EQ(apply_permutation(L * U, P), A_copy);
}
