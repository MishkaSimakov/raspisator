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
