#pragma once

#include <numeric>
#include <vector>

#include "Matrix.h"
#include "field/FieldTraits.h"
#include "utils/Accumulators.h"

namespace linalg {

// Returns rows of the matrix that form row basis. If the absolute value is <=
// @pivot_tolerance then it is treated as zero.
template <typename Field>
std::vector<size_t> get_row_basis(
    Matrix<Field> matrix,
    const Field pivot_tolerance = FieldTraits<Field>::tolerance) {
  using std::abs;

  if (matrix.rows() == 0) {
    return {};
  }

  std::vector<size_t> rows_map(matrix.rows());
  std::iota(rows_map.begin(), rows_map.end(), 0);

  auto [n, d] = matrix.shape();
  size_t current_row = 0;

  for (size_t col = 0; col < d; ++col) {
    // choose maximum element in the column for gaussian elimination
    ArgMaximum<Field> max_row;

    for (size_t row = current_row; row < n; ++row) {
      max_row.record(row, abs(matrix[row, col]));
    }

    const size_t max_row_index = max_row->index;

    if (abs(matrix[max_row_index, col]) <= pivot_tolerance) {
      continue;
    }

    if (max_row_index != current_row) {
      for (size_t i = 0; i < d; ++i) {
        std::swap(matrix[max_row_index, i], matrix[current_row, i]);
      }
      std::swap(rows_map[max_row_index], rows_map[current_row]);
    }

    // gaussian elimination
    const Field inverse_pivot = Field(1) / matrix[current_row, col];

    for (size_t i = 0; i < d; ++i) {
      matrix[current_row, i] *= inverse_pivot;
    }

    for (size_t i = current_row + 1; i < n; ++i) {
      Field coef = matrix[i, col];

      if (abs(coef) <= pivot_tolerance) {
        continue;
      }

      for (size_t j = 0; j < d; ++j) {
        matrix[i, j] =
            j == col ? 0 : matrix[i, j] - coef * matrix[current_row, j];
      }
    }

    ++current_row;

    if (current_row == std::min(n, d)) {
      break;
    }
  }

  rows_map.resize(current_row);
  return rows_map;
}

}  // namespace linalg
