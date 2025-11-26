#pragma once

#include "Matrix.h"
#include "linear/FieldTraits.h"

namespace linalg {
// gaussian elimination
struct NoopSubtractHook {
  void operator()(size_t minuend, size_t subtrahend, auto&& multiplier) {}
};

struct NoopMultiplyHook {
  void operator()(size_t multiplicand, auto&& multiplier) {}
};

// subtracts row #row_index from other rows, so the column #col_index becomes
// (0, ..., 1, 0, ..., 0)^T, where 1 is in row #row_index
template <typename SubtractHook = NoopSubtractHook,
          typename MultiplyHook = NoopMultiplyHook>
void gaussian_elimination(MatrixLike auto&& matrix, size_t row_index,
                          size_t col_index,
                          SubtractHook subtract_hook = SubtractHook{},
                          MultiplyHook multiply_hook = MultiplyHook{}) {
  using Field = matrix_field_t<decltype(matrix)>;
  auto [n, d] = matrix.shape();

  if (!FieldTraits<Field>::is_nonzero(matrix[row_index, col_index])) {
    throw std::invalid_argument(
        "Element must be non-zero for gaussian elimination.");
  }

  Field inverse = 1 / matrix[row_index, col_index];
  multiply_hook(row_index, inverse);

  matrix[row_index, {0, d}] *= inverse;

  for (size_t i = 0; i < n; ++i) {
    if (i == row_index ||
        !FieldTraits<Field>::is_nonzero(matrix[i, col_index])) {
      continue;
    }

    Field coef = matrix[i, col_index];
    subtract_hook(i, row_index, coef);

    matrix[i, {0, col_index}].sub_mul(matrix[row_index, {0, col_index}], coef);
    matrix[i, col_index] = 0;
    matrix[i, {col_index + 1, d}].sub_mul(matrix[row_index, {col_index + 1, d}],
                                          coef);
  }
}

template <MatrixLike T>
Matrix<matrix_field_t<T>> inverse(T&& matrix) {
  using Field = matrix_field_t<T>;

  auto [n, m] = matrix.shape();
  if (n != m) {
    throw std::invalid_argument("Matrix must be square for inverse.");
  }

  auto inv = Matrix<Field>::unity(n);

  auto multiplication_hook = [&inv, m](size_t multiplicand, Field multiplier) {
    inv[multiplicand, {0, m}] *= multiplier;
  };

  auto subtraction_hook = [&inv, m](size_t minuend, size_t subtrahend,
                                    Field multiplier) {
    inv[minuend, {0, m}].sub_mul(inv[subtrahend, {0, m}], multiplier);
  };

  for (size_t j = 0; j < n; ++j) {
    size_t maximizing_row = j;

    for (size_t i = j + 1; i < n; ++i) {
      if (FieldTraits<Field>::abs(matrix[i, j]) >
          FieldTraits<Field>::abs(matrix[maximizing_row, j])) {
        maximizing_row = i;
      }
    }

    if (!FieldTraits<Field>::is_nonzero(matrix[maximizing_row, j])) {
      throw std::invalid_argument("Matrix must be non-singular.");
    }

    if (j != maximizing_row) {
      std::swap(inv[j, {0, m}], inv[maximizing_row, {0, m}]);
      std::swap(matrix[j, {0, m}], matrix[maximizing_row, {0, m}]);
    }

    linalg::gaussian_elimination(matrix, j, j, subtraction_hook,
                                 multiplication_hook);
  }

  return inv;
}
}  // namespace linalg
