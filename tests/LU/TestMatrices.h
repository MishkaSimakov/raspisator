#pragma once

#include <random>
#include <vector>

#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"

inline Matrix<Rational> big_dense_matrix(size_t N) {
  Matrix<Rational> A(N, N);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      A[i, j] = static_cast<long long>(i + j) % 71 + 123;
    }
  }

  A = A + Matrix<Rational>::unity(N);

  return A;
}

inline Matrix<Rational> diagonal_matrix(size_t N) {
  Matrix<Rational> A(N, N);

  for (size_t i = 0; i < N; ++i) {
    A[i, i] = (i * 123) % 87 + 5;
  }

  return A;
}

inline Matrix<Rational> small_matrix1() { return {{1, 2}, {3, 4}}; }

inline Matrix<Rational> small_matrix2() { return {{3, 4}, {1, 2}}; }

inline Matrix<Rational> sparse_matrix(size_t N, size_t density_multiplier) {
  auto A = Matrix<Rational>::unity(N);

  std::mt19937 engine(0);
  std::uniform_int_distribution<size_t> uniform_dist(0, N - 1);

  // add O(N) non-zeros in "random" places
  for (size_t i = 0; i < density_multiplier * N; ++i) {
    size_t row = uniform_dist(engine);
    size_t col = uniform_dist(engine);

    A[row, col] = i + 1;
  }

  return A;
}

inline std::vector<std::pair<std::string, Matrix<Rational>>> test_matrices() {
  std::vector<std::pair<std::string, Matrix<Rational>>> result;

  result.emplace_back("small_matrix1", small_matrix1());
  result.emplace_back("small_matrix2", small_matrix2());
  result.emplace_back("big_dense_matrix", big_dense_matrix(50));
  result.emplace_back("diagonal_matrix", diagonal_matrix(50));
  result.emplace_back("sparse_matrix", sparse_matrix(50, 3));

  return result;
}
