#include <gtest/gtest.h>

#include "TestMatrices.h"
#include "linear/BigInteger.h"
#include "linear/matrix/Elimination.h"
#include "linear/matrix/Random.h"
#include "linear/sparse/CSCMatrix.h"
#include "linear/sparse/LU.h"

TEST(SparseLUTests, SolvesLinearSystem) {
  for (size_t N = 10; N < 1000; N *= 2) {
    auto matrix = sparse_matrix(N, 1);
    auto sparse = CSCMatrix(matrix);

    std::vector<size_t> columns(N);
    std::iota(columns.begin(), columns.end(), 0);

    auto [P, Q, ls, us] =
        linalg::FullPivotingLU<Rational>(N).get(sparse, columns);

    Matrix<Rational> b(N, 1, 123);
    auto x = linalg::solve_linear(b, P, Q, ls, us);

    ASSERT_EQ(matrix * x, b);
  }
}

TEST(SparseLUTests, SolvesTransposedLinearSystem) {
  for (size_t N = 10; N < 1000; N *= 2) {
    auto matrix = sparse_matrix(N, 1);
    auto sparse = CSCMatrix(matrix);

    std::vector<size_t> columns(N);
    std::iota(columns.begin(), columns.end(), 0);

    auto [P, Q, ls, us] =
        linalg::FullPivotingLU<Rational>(N).get(sparse, columns);

    Matrix<Rational> b(N, 1, 123);

    auto x = linalg::solve_linear_transposed(b, P, Q, ls, us);

    ASSERT_EQ(linalg::transposed(matrix) * x, b);
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

  auto expected =
      linalg::hstack(A[{0, 4}, 0], A[{0, 4}, 4], A[{0, 4}, 5], A[{0, 4}, 3]);

  ASSERT_EQ(inverse * expected, Matrix<Rational>::unity(4));
}
TEST(SparseLUTests, GetInverseMatrixTest) {
  CSCMatrix<Rational> A = {
      {1, 0, 0},
      {0, 2, 1},
      {0, 1, 0},
  };

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 1, 2});

  const auto inverse = lupa.get_inverse();
  const Matrix<Rational> expected = {
      {1, 0, 0},
      {0, 0, 1},
      {0, 1, -2},
  };

  ASSERT_EQ(inverse, expected);
}

TEST(SparseLUTests, GetMatrixTest) {
  CSCMatrix<Rational> A = {
      {1, 1, 0},
      {0, 2, 1},
      {0, 3, 0},
  };

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 1, 2});

  const auto matrix = lupa.get_matrix();
  const auto expected = linalg::to_dense(A);

  ASSERT_EQ(matrix, expected);
}

TEST(SparseLUTests, GetMatrixRandomTest) {
  constexpr size_t size = 10;

  std::default_random_engine random;

  for (size_t i = 0; i < 100; ++i) {
    CSCMatrix<Rational> A(linalg::random_invertible<Rational>(size, random));

    auto lupa = linalg::LUPA(A);

    std::vector<size_t> columns(size);
    std::iota(columns.begin(), columns.end(), 0);
    lupa.set_columns(columns);

    const auto matrix = lupa.get_matrix();
    const auto expected = linalg::to_dense(A);

    ASSERT_EQ(matrix, expected);
  }
}

TEST(SparseLUTests, ChangeColumnsRandomTest) {
  constexpr size_t size = 10;

  std::default_random_engine random;

  for (size_t i = 0; i < 100; ++i) {
    SCOPED_TRACE(std::format("iteration: {}", i));

    const auto core = linalg::random_invertible<Rational>(size, random);
    const auto dense = linalg::hstack(core, core);
    const auto sparse = CSCMatrix<Rational>(dense);

    auto lupa = linalg::LUPA(sparse);

    std::vector<size_t> columns(size);
    std::iota(columns.begin(), columns.end(), 0);
    lupa.set_columns(columns);

    for (size_t j = 0; j < size; ++j) {
      lupa.change_column(j, j + size);
    }

    const auto matrix = lupa.get_matrix();

    ASSERT_EQ(matrix, core);
  }
}
