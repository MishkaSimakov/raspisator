#pragma once

#include <algorithm>

#include "Concepts.h"

namespace linalg {

template <ElementWiseMatrixRange M>
double norm(M&& matrix) {
  using std::sqrt;

  const auto [n, d] = matrix.shape();
  MatrixFieldType<M> result = 0;

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      result += matrix[i, j] * matrix[i, j];
    }
  }

  return sqrt(result);
}

template <ElementWiseMatrixRange M>
MatrixFieldType<M> inf_norm(M&& matrix) {
  using std::abs, std::max;
  using Field = MatrixFieldType<M>;

  const auto [n, d] = matrix.shape();
  Field result = 0;

  for (size_t row = 0; row < n; ++row) {
    Field row_sum = 0;

    for (size_t col = 0; col < d; ++col) {
      row_sum += abs(matrix[row, col]);
    }

    result = std::max(result, row_sum);
  }

  return result;
}

}  // namespace linalg
