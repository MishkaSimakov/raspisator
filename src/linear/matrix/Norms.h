#pragma once

#include <algorithm>
#include <cmath>

#include "Types.h"

namespace linalg {

template <MatrixLike T>
double norm(T&& matrix)
  requires(std::same_as<matrix_field_t<T>, double>)
{
  auto [n, d] = matrix.shape();

  double result = 0;

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      result += matrix[i, j] * matrix[i, j];
    }
  }

  return std::sqrt(result);
}

template <MatrixLike T>
double inf_norm(T&& matrix)
  requires(std::same_as<matrix_field_t<T>, double>)
{
  auto [n, d] = matrix.shape();

  double result = 0;

  for (size_t row = 0; row < n; ++row) {
    double row_sum = 0;

    for (size_t col = 0; col < d; ++col) {
      row_sum += std::abs(matrix[row, col]);
    }

    result = std::max(result, row_sum);
  }

  return result;
}

}  // namespace linalg
