#include <gtest/gtest.h>

#include "TestMatrices.h"
#include "linear/BigInteger.h"
#include "linear/matrix/LU.h"
#include "linear/matrix/Matrix.h"

TEST(DenseLUTests, SimpleInPlaceLU) {
  Matrix<Rational> A = {{1, 2}, {3, 4}};

  linalg::inplace_lu(A);

  ASSERT_EQ(A, (Matrix<Rational>{{1, 2}, {3, -2}}));
}

TEST(DenseLUTests, BigInPlaceLU) {
  const size_t N = 50;

  auto A = big_dense_matrix(N);
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

TEST(DenseLUTests, SimpleLUP) {
  Matrix<Rational> A = {{1, 2}, {3, 4}};

  auto [L, U, P] = linalg::get_lup(A);

  auto product = linalg::apply_permutation(L * U, P);
  ASSERT_EQ(product, A);
}

TEST(DenseLUTests, BigInPlaceLUP) {
  const size_t N = 50;

  auto A = big_dense_matrix(N);
  auto A_copy = A;
  auto P = linalg::inplace_lup(A);

  Matrix<Rational> L(N, N, 0);
  Matrix<Rational> U(N, N, 0);

  for (size_t i = 0; i < N; ++i) {
    L[i, {0, i}] = A[i, {0, i}];
    L[i, i] = 1;

    U[i, {i, N}] = A[i, {i, N}];
  }

  // LU = PA
  ASSERT_EQ(L * U, linalg::apply_permutation(A_copy, P));
}

TEST(DenseLUTests, SolveL) {
  const size_t N = 50;

  auto A = big_dense_matrix(N);
  auto [L, U, P] = linalg::get_lup(A);

  auto b = Matrix<Rational>(N, 1, 1);

  auto s = linalg::solve_lower(L, b, std::true_type{});

  ASSERT_EQ(s.shape(), (std::pair{N, 1}));
  ASSERT_EQ(L * s, b);
}

TEST(DenseLUTests, SolveSmallU) {
  Matrix<Rational> U = {{1, 2, 3}, {0, 4, 5}, {0, 0, 6}};
  auto b = Matrix<Rational>(3, 1, 1);

  auto s = linalg::solve_upper(U, b, std::false_type{});
  Matrix expected = {{Rational{5} / 12}, {Rational{1} / 24}, {Rational{1} / 6}};

  ASSERT_EQ(s.shape(), (std::pair{3, 1}));
  ASSERT_EQ(s, expected);
}

TEST(DenseLUTests, SolveBigU) {
  const size_t N = 40;

  auto A = big_dense_matrix(N);
  auto [L, U, P] = linalg::get_lup(A);

  auto b = Matrix<Rational>(N, 1, 1);

  auto s = linalg::solve_upper(U, b, std::false_type{});

  ASSERT_EQ(s.shape(), (std::pair{N, 1}));
  ASSERT_EQ(U * s, b);
}

TEST(DenseLUTests, SolveUsingLUChain) {
  const size_t N = 50;

  auto A = big_dense_matrix(N);
  auto [L, U, P] = linalg::get_lup(A);

  auto b = Matrix<Rational>(N, 1);
  for (size_t i = 0; i < N; ++i) {
    b[i, 0] = i;
  }

  // solving Ax = b
  // PA = LU
  // Ax = b iff LUx = Pb

  // calculate Pb
  auto Pb = linalg::apply_permutation(b, P);

  // solve Ly = Pb
  auto y = linalg::solve_lower(L, Pb, std::true_type{});

  // solve Ux = y
  auto x = linalg::solve_upper(U, y, std::false_type{});

  ASSERT_EQ(x.shape(), (std::pair{N, 1}));
  ASSERT_EQ(A * x, b);
}

TEST(DenseLUTests, SolveTransposedUsingLUChain) {
  const size_t N = 50;

  auto A = big_dense_matrix(N);
  auto [L, U, P] = linalg::get_lup(A);

  auto b = Matrix<Rational>(N, 1);
  for (size_t i = 0; i < N; ++i) {
    b[i, 0] = i;
  }

  // solving A^T x = b
  // PA = LU
  // Ax = b iff U^T L^T P^-T x = b

  // solve U^T z = b
  auto z = linalg::solve_lower(linalg::transposed(U), b, std::false_type{});

  // solve L^T y = z
  auto y = linalg::solve_upper(linalg::transposed(L), z, std::true_type{});

  // solve P^-T x = y
  // P^-T x = y iff x = P^T y

  // transpose permutation
  std::vector<size_t> transposed_P(N);
  for (size_t i = 0; i < N; ++i) {
    transposed_P[P[i]] = i;
  }

  auto x = linalg::apply_permutation(y, transposed_P);

  ASSERT_EQ(x.shape(), (std::pair{N, 1}));
  ASSERT_EQ(linalg::transposed(A) * x, b);
}
