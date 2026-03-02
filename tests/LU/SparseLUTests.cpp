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

    // returns a point y, such that x[i] = y[P[i]]
    auto y = linalg::solve_transposed_linear(L, U, P, b);

    Matrix<Rational> x(N, 1);
    for (size_t i = 0; i < N; ++i) {
      x[i, 0] = y[P[i], 0];
    }

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

  std::vector<size_t> columns(matrix.get_width());
  std::iota(columns.begin(), columns.end(), 0);

  for (size_t i = 0; i < 5; ++i) {
    auto [L, U, P] = linalg::sparse_lup(sparse, columns);

    auto dense_L = linalg::to_dense(L);
    auto dense_U = linalg::to_dense(U);

    // add ones on the diagonal
    for (size_t i = 0; i < matrix.get_height(); ++i) {
      dense_L[i, i] = 1;
    }

    auto expected = linalg::to_dense(P.apply(sparse));

    ASSERT_EQ(dense_L * dense_U, expected);
  }
}

TEST(SparseLUTests, SmallSolveLinearTransposedTest) {
  Matrix<Rational> A = {
      {3, -7, -2, 2},
      {-3, 5, 1, 0},
      {6, -4, 0, -5},
      {-9, 5, -5, 12},
  };

  auto sparse = CSCMatrix(A);
  auto lupa = linalg::LUPA(sparse);

  lupa.set_columns(std::vector<size_t>{0, 1, 2, 3});

  Matrix<Rational> b = {{0}, {2}, {0}, {1}};
  auto solution = lupa.solve_linear_transposed(b);

  Matrix<Rational> expected = {{30}, {50}, {7}, {-2}};

  ASSERT_EQ(solution, expected);
}

TEST(SparseLUTests, SmallGetRowTest) {
  Matrix<Rational> A = {
      {3, -7, -2, 2},
      {-3, 5, 1, 0},
      {6, -4, 0, -5},
      {-9, 5, -5, 12},
  };

  auto sparse = CSCMatrix(A);
  auto lupa = linalg::LUPA(sparse);

  lupa.set_columns(std::vector<size_t>{0, 1, 2, 3});

  auto row = lupa.get_row(2);
  Matrix<Rational> expected = {{11}, {18}, {2}, {-1}};

  ASSERT_EQ(row, expected);
}

TEST(SparseLUTests, ChangeColumnTest) {
  Matrix<Rational> A = {
    {3, -7, -2, 2, 1, 1},
    {-3, 5, 1, 0, 0, 2},
    {6, -4, 0, -5, 2, 3},
    {-9, 5, -5, 12, 3, 4},
};

  auto sparse = CSCMatrix(A);
  auto lupa = linalg::LUPA(sparse);

  lupa.set_columns(std::vector<size_t>{0, 1, 2, 3});

  lupa.change_column(1, 4);
  lupa.change_column(2, 5);

  auto inverse = lupa.get_inverse();

  auto expected = linalg::hstack(A[{0, 4}, 0], A[{0, 4}, 4], A[{0, 4}, 5], A[{0, 4}, 3]);

  ASSERT_EQ(inverse * expected, Matrix<Rational>::unity(4));
}
