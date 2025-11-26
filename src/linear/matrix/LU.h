#pragma once

#include <numeric>

#include "Matrix.h"
#include "linear/FieldTraits.h"

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

// solves Ax = b, where A is lower triangular
template <MatrixLike U, MatrixLike V, bool ones_on_diagonal>
Matrix<common_field_t<U, V>> solve_lower(U&& A, V&& b,
                                         std::bool_constant<ones_on_diagonal>) {
  using Field = common_field_t<U, V>;

  Matrix<Field> result = b;

  for (size_t i = 0; i < b.get_height(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      result[i, 0] -= result[j, 0] * A[i, j];
    }

    if constexpr (!ones_on_diagonal) {
      result[i, 0] /= A[i, i];
    }
  }

  return result;
}

// solves Ax = b, where A is upper triangular matrix
template <MatrixLike U, MatrixLike V, bool ones_on_diagonal>
Matrix<common_field_t<U, V>> solve_upper(U&& A, V&& b,
                                         std::bool_constant<ones_on_diagonal>) {
  using Field = common_field_t<U, V>;

  Matrix<Field> result = b;
  const size_t n = b.get_height();

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = n - i; j < n; ++j) {
      result[n - i - 1, 0] -= result[j, 0] * A[n - i - 1, j];
    }

    if constexpr (!ones_on_diagonal) {
      result[n - i - 1, 0] /= A[n - i - 1, n - i - 1];
    }
  }

  return result;
}

// this method combines with LUP decomposition
// if LU = PA, then L * U = apply_permutation(A, P)
template <MatrixLike T>
Matrix<matrix_field_t<T>> apply_permutation(
    T&& matrix, const std::vector<size_t>& permutation) {
  auto [n, d] = matrix.shape();

  Matrix<matrix_field_t<T>> result(n, d);

  for (size_t i = 0; i < n; ++i) {
    result[i, {0, d}] = matrix[permutation[i], {0, d}];
  }

  return result;
}

}  // namespace linalg
