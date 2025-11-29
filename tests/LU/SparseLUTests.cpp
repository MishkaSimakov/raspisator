#include <gtest/gtest.h>

#include "TestMatrices.h"
#include "linear/BigInteger.h"
#include "linear/sparse/CSCMatrix.h"
#include "linear/sparse/LU.h"

TEST(SparseLUTests, DoesNotContainZeros) {
  auto matrix = sparse_matrix(1'000, 10);
  auto sparse = CSCMatrix(matrix);

  auto [L, U, P] = linalg::sparse_lu(sparse);

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

    auto [L, U, P] = linalg::sparse_lu(sparse);

    Matrix<Rational> b(N, 1, 123);
    auto x = linalg::solve_linear(L, U, P, b);

    ASSERT_EQ(matrix * x, b);
  }
}

TEST(SparseLUTests, SolvesTransposedLinearSystem) {
  for (size_t N = 10; N < 1000; N *= 2) {
    auto matrix = sparse_matrix(N, 1);
    auto sparse = CSCMatrix(matrix);

    auto [L, U, P] = linalg::sparse_lu(sparse);

    Matrix<Rational> b(N, 1, 123);
    auto x = linalg::solve_transposed_linear(L, U, P, b);

    ASSERT_EQ(linalg::transposed(matrix) * x, b);
  }
}
