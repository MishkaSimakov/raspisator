#pragma once

#include <cassert>
#include <numeric>
#include <vector>

#include "Elimination.h"
#include "Matrix.h"
#include "linear/FieldTraits.h"

namespace linalg {

template <typename Field>
std::vector<size_t> get_row_basis(Matrix<Field> matrix) {
  std::vector<size_t> rows_map(matrix.get_height());
  std::iota(rows_map.begin(), rows_map.end(), 0);

  auto [n, d] = matrix.shape();
  size_t current_row = 0;

  for (size_t col = 0; col < d; ++col) {
    // choose maximum element in the column for gaussian elimination
    size_t maximizing_row = current_row;

    for (size_t row = current_row + 1; row < n; ++row) {
      if (FieldTraits<Field>::abs(matrix[row, col]) >
          FieldTraits<Field>::abs(matrix[maximizing_row, col])) {
        maximizing_row = row;
      }
    }

    if (!FieldTraits<Field>::is_nonzero(matrix[maximizing_row, col])) {
      continue;
    }

    if (maximizing_row != current_row) {
      std::swap(matrix[maximizing_row, {0, d}], matrix[current_row, {0, d}]);
      std::swap(rows_map[maximizing_row], rows_map[current_row]);
    }

    linalg::gaussian_elimination(matrix, current_row, col);

    ++current_row;

    if (current_row == std::min(n, d)) {
      break;
    }
  }

  rows_map.resize(current_row);
  return rows_map;
}

// Accepts matrix A and a vector of rows indices.
// Rows indices must be sorted!
// returns a vector of rows indices such that:
// 1. new indices contain all initial indices
// 2. new indices form the biggest linearly independent subset of rows
template <typename Field>
std::vector<size_t> complete_row_basis(
    Matrix<Field> matrix, const std::vector<size_t>& partial_basis) {
  std::vector<size_t> rows_map(matrix.get_height());
  std::iota(rows_map.begin(), rows_map.end(), 0);

  std::vector<size_t> result;

  auto [n, d] = matrix.shape();

  // first we move partial basis rows to the top of A
  // IT IS INCREDIBLY IMPORTANT THAT partial_basis IS SORTED!
  // ANYONE WHO FORGETS ABOUT THIS SHALL BE BANISHED FROM PROGRAMMING!!!
  // TODO: somehow remove this constraint, it is stupendously bug-prone
  for (size_t i = 0; i < partial_basis.size(); ++i) {
    if (i != partial_basis[i]) {
      std::swap(matrix[i, {0, d}], matrix[partial_basis[i], {0, d}]);
      std::swap(rows_map[i], rows_map[partial_basis[i]]);
    }
  }

  // then we eliminate a few columns using partial_basis vectors
  for (size_t row = 0; row < n; ++row) {
    // choose maximum element in the row for gaussian elimination
    size_t current_col = result.size();
    size_t maximizing_col = current_col;

    for (size_t col = current_col + 1; col < d; ++col) {
      if (FieldTraits<Field>::abs(matrix[row, col]) >
          FieldTraits<Field>::abs(matrix[row, maximizing_col])) {
        maximizing_col = col;
      }
    }

    if (!FieldTraits<Field>::is_nonzero(matrix[row, maximizing_col])) {
      assert(row >= partial_basis.size());

      continue;
    }

    // we can swap columns freely
    if (maximizing_col != current_col) {
      std::swap(matrix[{0, n}, maximizing_col], matrix[{0, n}, current_col]);
    }

    // gaussian elimination
    for (size_t i = row + 1; i < n; ++i) {
      if (!FieldTraits<Field>::is_nonzero(matrix[i, current_col])) {
        continue;
      }

      Field coef = matrix[i, current_col] / matrix[row, current_col];
      matrix[i, {current_col + 1, d}].sub_mul(matrix[row, {current_col + 1, d}],
                                              coef);
    }

    result.push_back(rows_map[row]);

    if (result.size() == std::min(n, d)) {
      break;
    }
  }

  return result;
}
}  // namespace linalg
