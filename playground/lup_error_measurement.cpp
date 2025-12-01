#include <iostream>
#include <print>
#include <random>

#include "linear/BigInteger.h"
#include "linear/MPS.h"
#include "linear/SimplexMethod.h"
#include "linear/bb/BranchAndBound.h"
#include "linear/bb/Drawer.h"
#include "linear/builder/ProblemBuilder.h"
#include "linear/sparse/CSCMatrix.h"
#include "linear/sparse/LU.h"
#include "problems/Dwarf.h"
#include "utils/Drawing.h"
#include "utils/Hashers.h"

using Field = double;

template <typename Field>
CSCMatrix<Field> sparse_matrix_1(size_t N, size_t density_multiplier,
                                 size_t seed = 0) {
  auto A = Matrix<Field>::unity(N);

  std::mt19937 engine(seed);
  std::uniform_int_distribution<size_t> uniform_dist(0, N - 1);

  // add O(N) non-zeros in "random" places
  for (size_t i = 0; i < density_multiplier * N; ++i) {
    size_t row = uniform_dist(engine);
    size_t col = uniform_dist(engine);

    A[row, col] = i + 1;
  }

  return CSCMatrix(A);
}

double max_difference(const CSCMatrix<double>& left,
                      const CSCMatrix<Rational>& right) {
  auto dense_left = linalg::to_dense(left);
  auto dense_right = linalg::to_dense(right);

  auto [n, d] = dense_right.shape();
  double max_difference = 0;

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      double difference =
          std::abs(static_cast<double>(dense_right[i, j]) - dense_left[i, j]);

      max_difference = std::max(difference, max_difference);
    }
  }

  return max_difference;
}

int main() {
  size_t density_multiplier = 3;

  std::println("size,error");

  for (size_t N = 10; N < 75; ++N) {
    for (size_t i = 0; i < 10; ++i) {
      auto A_rational = sparse_matrix_1<Rational>(N, density_multiplier, i);
      auto A_double = sparse_matrix_1<double>(N, density_multiplier, i);

      auto [rational_L, rational_U, rational_P] =
          linalg::sparse_lup(A_rational);
      auto [double_L, double_U, double_P] = linalg::sparse_lup(A_double);

      double error = std::max(max_difference(double_L, rational_L),
                              max_difference(double_U, rational_U));

      std::println("{},{}", N, error);
    }
  }
}
