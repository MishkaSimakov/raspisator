#pragma once

#include "Matrix.h"

namespace linalg {
template <MatrixLike T>
std::pair<Matrix<matrix_field_t<T>>, Matrix<matrix_field_t<T>>> get_lu(
    T&& matrix) {
  using Field = matrix_field_t<T>;

  auto [n, d] = matrix.shape();
  if (n != d) {
    throw DimensionsException(
        "Square matrix is required for LU decomposition.");
  }

  Matrix<Field> U = matrix;
  Matrix<Field> L = Matrix<Field>::unity(n);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      L[j, i] = U[j, i] / U[i, i];
      U[j, {i, n}].sub_mul(U[i, {i, n}], L[j, i]);
    }
  }

  return std::pair{std::move(L), std::move(U)};
}

template <MatrixLike T>
void inplace_lu(T&& matrix) {
  using Field = matrix_field_t<T>;

  auto [n, d] = matrix.shape();
  if (n != d) {
    throw DimensionsException(
        "Square matrix is required for LU decomposition.");
  }

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      matrix[j, i] /= matrix[i, i];
      matrix[j, {i + 1, n}].sub_mul(matrix[i, {i + 1, n}], matrix[j, i]);
    }
  }
}
}  // namespace linalg
