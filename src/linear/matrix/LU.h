#pragma once

#include <numeric>

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

void inplace_lu(MatrixLike auto&& matrix) {
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

template <MatrixLike T>
std::tuple<Matrix<matrix_field_t<T>>, Matrix<matrix_field_t<T>>,
           std::vector<size_t>>
get_lup(T&& matrix) {
  using Field = matrix_field_t<T>;

  auto [n, d] = matrix.shape();
  if (n != d) {
    throw DimensionsException(
        "Square matrix is required for LU decomposition.");
  }

  Matrix<Field> U = matrix;
  Matrix<Field> L = Matrix<Field>::unity(n);

  std::vector<size_t> P(n);
  std::iota(P.begin(), P.end(), 0);

  for (size_t i = 0; i < n; ++i) {
    size_t maximizing_row = i;
    for (size_t k = i + 1; k < n; ++k) {
      if (FieldTraits<Field>::abs(U[k, i]) >
          FieldTraits<Field>::abs(U[maximizing_row, i])) {
        maximizing_row = k;
      }
    }

    std::swap(U[i, {i, n}], U[maximizing_row, {i, n}]);
    std::swap(L[i, {0, i}], L[maximizing_row, {0, i}]);
    std::swap(P[i], P[maximizing_row]);

    for (size_t j = i + 1; j < n; ++j) {
      L[j, i] = U[j, i] / U[i, i];
      U[j, {i, n}].sub_mul(U[i, {i, n}], L[j, i]);
    }
  }

  return std::tuple{std::move(L), std::move(U), std::move(P)};
}

template <MatrixLike T>
std::vector<size_t> inplace_lup(T&& matrix) {
  using Field = matrix_field_t<T>;

  auto [n, d] = matrix.shape();
  if (n != d) {
    throw DimensionsException(
        "Square matrix is required for LU decomposition.");
  }

  std::vector<size_t> P(n);
  std::iota(P.begin(), P.end(), 0);

  for (size_t i = 0; i < n; ++i) {
    size_t maximizing_row = i;
    for (size_t k = i + 1; k < n; ++k) {
      if (FieldTraits<Field>::abs(matrix[k, i]) >
          FieldTraits<Field>::abs(matrix[maximizing_row, i])) {
        maximizing_row = k;
      }
    }

    std::swap(matrix[i, {i, n}], matrix[maximizing_row, {i, n}]);
    std::swap(matrix[i, {0, i}], matrix[maximizing_row, {0, i}]);
    std::swap(P[i], P[maximizing_row]);

    for (size_t j = i + 1; j < n; ++j) {
      matrix[j, i] /= matrix[i, i];
      matrix[j, {i + 1, n}].sub_mul(matrix[i, {i + 1, n}], matrix[j, i]);
    }
  }

  return P;
}

}  // namespace linalg
