#include <gtest/gtest.h>

#include "TestMatrices.h"
#include "linear/BigInteger.h"
#include "linear/sparse/CSCMatrix.h"
#include "linear/sparse/LU.h"

TEST(SparseLUTests, DoesNotContainZeros) {
  auto matrix = sparse_matrix(50, 10);
  auto sparse = CSCMatrix(matrix);

  std::vector<size_t> columns(50);
  std::iota(columns.begin(), columns.end(), 0);

  auto [L, U, P] = linalg::sparse_lup(sparse, columns);

  for (size_t col = 0; col < matrix.get_width(); ++col) {
    for (const Rational& value : L.get_column(col) | std::views::values) {
      ASSERT_NE(value, 0);
    }

    for (const Rational& value : U.get_column(col) | std::views::values) {
      ASSERT_NE(value, 0);
    }
  }
}

TEST(SparseLUTests, SolvesLinearSystem) {
  for (size_t N = 10; N < 1000; N *= 2) {
    auto matrix = sparse_matrix(N, 1);
    auto sparse = CSCMatrix(matrix);

    std::vector<size_t> columns(N);
    std::iota(columns.begin(), columns.end(), 0);

    auto [L, U, P] = linalg::sparse_lup(sparse, columns);

    Matrix<Rational> b(N, 1, 123);
    auto x = linalg::solve_linear(L, U, P, b);

    ASSERT_EQ(matrix * x, b);
  }
}

TEST(SparseLUTests, SolvesTransposedLinearSystem) {
  for (size_t N = 10; N < 1000; N *= 2) {
    auto matrix = sparse_matrix(N, 1);
    auto sparse = CSCMatrix(matrix);

    std::vector<size_t> columns(N);
    std::iota(columns.begin(), columns.end(), 0);

    auto [L, U, P] = linalg::sparse_lup(sparse, columns);

    Matrix<Rational> b(N, 1, 123);
    auto x = linalg::solve_transposed_linear(L, U, P, b);

    ASSERT_EQ(linalg::transposed(matrix) * x, b);
  }
}

TEST(SparseLUTests, DoubleMatrix) {
  size_t N = 50;

  auto matrix = Matrix<double>::unity(N);
  auto sparse = CSCMatrix(matrix);

  std::vector<size_t> columns(N);
  std::iota(columns.begin(), columns.end(), 0);

  auto [L, U, P] = linalg::sparse_lup(sparse, columns);

  auto dense_L = linalg::to_dense(L);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      ASSERT_DOUBLE_EQ((dense_L[i, j]), 0);
    }
  }

  auto dense_U = linalg::to_dense(U);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      ASSERT_DOUBLE_EQ((dense_U[i, j]), i == j ? 1 : 0);
    }
  }
}

TEST(SparseLUTests, LUPATest) {
  size_t N = 10;

  auto matrix = sparse_matrix(N, 5);
  auto sparse = CSCMatrix(matrix);

  linalg::LUPA<Rational> lupa(sparse);

  std::vector<size_t> columns(matrix.get_width());
  std::iota(columns.begin(), columns.end(), 0);

  CSCMatrix<Rational> L(N);
  CSCMatrix<Rational> U(N);
  std::vector<size_t> P(N);

  for (size_t i = 0; i < 5; ++i) {
    lupa.get_lup(columns, L, U, P);

    auto dense_L = linalg::to_dense(L);
    auto dense_U = linalg::to_dense(U);

    // add ones on the diagonal
    for (size_t i = 0; i < matrix.get_height(); ++i) {
      dense_L[i, i] = 1;
    }

    auto expected = linalg::to_dense(linalg::apply_permutation(sparse, P));

    ASSERT_EQ(dense_L * dense_U, expected);
  }
}
