#pragma once

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

    matrix.swap_rows(maximizing_row, current_row);
    std::swap(rows_map[maximizing_row], rows_map[current_row]);
    linalg::gaussian_elimination(matrix, current_row, col);

    ++current_row;
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
std::vector<size_t> make_row_basis(Matrix<Field> matrix,
                                   const std::vector<size_t>& partial_basis) {
  std::vector<size_t> rows_map(matrix.get_height());
  std::iota(rows_map.begin(), rows_map.end(), 0);

  auto [n, d] = matrix.shape();
  size_t current_row = 0;

  // first we move partial basis rows to the top of A
  // IT IS INCREDIBLY IMPORTANT THAT partial_basis IS SORTED!
  // ANYONE WHO FORGETS ABOUT THIS SHALL BE BANISHED FROM PROGRAMMING!!!
  // TODO: somehow remove this constraint, it is stupendously bug-prone
  for (size_t i = 0; i < partial_basis.size(); ++i) {
    matrix.swap_rows(i, partial_basis[i]);
    std::swap(rows_map[i], rows_map[partial_basis[i]]);
  }

  // then we eliminate a few columns using partial_basis vectors
  size_t fixed_cnt = partial_basis.size();

  for (size_t col = 0; col < fixed_cnt; ++col) {
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

    matrix.swap_rows(maximizing_row, current_row);
    std::swap(rows_map[maximizing_row], rows_map[current_row]);
    matrix.gaussian_elimination(current_row, col);

    ++current_row;
  }

  rows_map.resize(current_row);
  return rows_map;
}
}  // namespace linalg
