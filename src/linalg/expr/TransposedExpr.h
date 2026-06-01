#pragma once

#include "linalg/Concepts.h"

namespace linalg {

template <SomeMatrixLike Matrix>
class TransposedExpr {
  Matrix& matrix_;

 public:
  using FieldType = typename Matrix::FieldType;

  explicit TransposedExpr(Matrix& matrix) : matrix_(matrix) {}

  //
  auto& operator[](size_t row, size_t col) { return matrix_[col, row]; }
  const auto& operator[](size_t row, size_t col) const {
    return matrix_[col, row];
  }

  //
  size_t rows() const { return matrix_.cols(); }
  size_t cols() const { return matrix_.rows(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

template <SomeMatrixLike Matrix>
auto transposed(Matrix& matrix) {
  return TransposedExpr(matrix);
}

}  // namespace linalg
