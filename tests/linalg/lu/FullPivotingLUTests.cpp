#include <gtest/gtest.h>

#include "linalg/CSCMatrix.h"
#include "linalg/Transpose.h"
#include "linalg/lu/FullPivotingLU.h"
#include "linalg/lu/Solve.h"
#include "linalg/lu/TestMatrices.h"

using namespace linalg;

TEST(FullPivotingLUTests, SolvesLinearSystem) {
  for (size_t N = 10; N < 1000; N *= 2) {
    auto matrix = CSCMatrix(sparse_matrix(N, 1));

    std::vector<size_t> columns(N);
    std::iota(columns.begin(), columns.end(), 0);

    auto [P, Q, ls, us] =
        linalg::FullPivotingLU<Rational>(N).get(matrix, columns);

    auto b = Vector<Rational>::ones(N);
    auto x = solve_linear(b, P, Q, ls, us);

    ASSERT_EQ(matrix * x, b);
  }
}

TEST(FullPivotingLUTests, SolvesTransposedLinearSystem) {
  for (size_t N = 10; N < 1000; N *= 2) {
    auto matrix = sparse_matrix(N, 1);
    auto sparse = CSCMatrix(matrix);

    std::vector<size_t> columns(N);
    std::iota(columns.begin(), columns.end(), 0);

    auto [P, Q, ls, us] =
        linalg::FullPivotingLU<Rational>(N).get(sparse, columns);

    auto b = Vector<Rational>::ones(N);
    auto x = solve_linear_transposed(b, P, Q, ls, us);

    ASSERT_EQ(transpose(matrix) * x, b);
  }
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

TEST(FullPivotingLUTests, DecomposeThenCompose) {
  for (const auto& [name, matrix] : test_matrices()) {
    auto [n, d] = matrix.shape();
    auto sparse = CSCMatrix(matrix);

    std::vector<size_t> columns(matrix.cols());
    std::iota(columns.begin(), columns.end(), 0);

    auto [P, Q, ls, us] =
        linalg::FullPivotingLU<Rational>(n).get(sparse, columns);

    auto L = Matrix<Rational>::identity(n);
    for (auto entry : ls) {
      L = entry.apply(std::move(L));
    }

    auto U = Matrix<Rational>::identity(n);
    for (auto entry : us | std::views::reverse) {
      U = entry.apply(std::move(U));
    }

    ASSERT_NO_FATAL_FAILURE(check_L(L));
    ASSERT_NO_FATAL_FAILURE(check_U(U));

    ASSERT_EQ(Q * Matrix(U * L) * P * matrix, Matrix<Rational>::identity(n));
  }
}
