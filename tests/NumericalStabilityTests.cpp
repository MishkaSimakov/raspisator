#include <gtest/gtest.h>

#include "Assertions.h"
#include "linear/SimplexMethod.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/RowBasis.h"
#include "linear/model/LP.h"

TEST(ProblemBuilderTests, RowBasisNumericalStability) {
  size_t N = 10;
  double epsilon = 1e-20;

  Matrix<double> A(2 * N, N);

  for (size_t i = 0; i < N; ++i) {
    A[N + i, i] = 1;
  }

  // introduce small errors into A
  for (size_t i = 0; i < A.get_height(); ++i) {
    for (size_t j = 0; j < A.get_width(); ++j) {
      A[i, j] += (static_cast<double>((i + j) % 5) - 2) * epsilon;
    }
  }

  auto row_basis = linalg::get_row_basis(A);

  std::vector<size_t> expected(N);
  std::iota(expected.begin(), expected.end(), 10);

  ASSERT_SETS_EQ(row_basis, expected);
}

// TODO: implement numerical stability tests
// TEST(SiplexMethodTests, RoundingErrorsProblem) {
//   Matrix<double> c = {
//       {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
//
//   Matrix<double> A = {
//       {1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//       {1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
//       {0, -1, 10, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
//       {0, 1, -200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
//       {0, 0, 0, -1, 10, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
//       {0, 0, 0, 1, -200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
//       {0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//       {0, 0, 0, -1, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//       {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//       {0, 1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//       {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
//       {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
//       {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
//       {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
//   };
//
//   Matrix<double> b = {{1, 1, 0, 0, 0, 0, -100, 0, 0, 0, 1, 1, 1, -100}};
//   b = linalg::transposed(b);
//
//   // adding small errors
//   for (size_t i = 0; i < A.get_height(); ++i) {
//     for (size_t j = 0; j < A.get_width(); ++j) {
//       A[i, j] += (static_cast<double>((i + j) % 5) - 2) * 1e-15;
//     }
//   }
//
//   std::cout << linalg::hstack(A, b) << std::endl;
//
//   auto solution =
//       std::get<FiniteLPSolution<double>>(SimplexMethod(A, b, c).solve());
//
//   std::cout << solution.point << std::endl;
// }
