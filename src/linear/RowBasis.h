#pragma once

#include <numeric>
#include <vector>

#include "matrix/Matrix.h"

template <typename Field>
std::vector<size_t> get_row_basis(Matrix<Field> matrix) {
  std::vector<size_t> rows_map(matrix.get_height());
  std::iota(rows_map.begin(), rows_map.end(), 0);

  auto [n, d] = matrix.shape();
  size_t current_row = 0;

  for (size_t col = 0; col < d; ++col) {
    size_t nonzero_row = n;

    for (size_t row = current_row; row < n; ++row) {
      if (matrix[row, col] != 0) {
        nonzero_row = row;
        break;
      }
    }

    if (nonzero_row == n) {
      continue;
    }

    matrix.swap_rows(nonzero_row, current_row);
    std::swap(rows_map[nonzero_row], rows_map[current_row]);
    matrix.gaussian_elimination(current_row, col);

    ++current_row;
  }

  // find all non-zero rows
  std::vector<size_t> result;

  for (size_t row = 0; row < n; ++row) {
    bool is_nonzero = false;

    for (size_t col = 0; col < d; ++col) {
      if (matrix[row, col] != 0) {
        is_nonzero = true;
        break;
      }
    }

    if (is_nonzero) {
      result.push_back(rows_map[row]);
    }
  }

  return result;
}
