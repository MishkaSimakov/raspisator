#pragma once

#include "linalg/Concepts.h"

namespace linalg {

template <SomeMatrixLike Matrix, IndicesRange RowRange, IndicesRange ColRange>
class IndexedExpr {
  Matrix& matrix_;
  const RowRange& rows_;
  const ColRange& cols_;

 public:
  using FieldType = typename Matrix::FieldType;

  IndexedExpr(Matrix& matrix, const RowRange& rows, const ColRange& cols)
      : matrix_(matrix), rows_(rows), cols_(cols) {}

  //
  auto& operator[](size_t row, size_t col) {
    return matrix_[rows_.begin()[row], cols_.begin()[col]];
  }

  const auto& operator[](size_t row, size_t col) const {
    return matrix_[rows_.begin()[row], cols_.begin()[col]];
  }

  //
  size_t rows() const { return std::ranges::size(rows_); }
  size_t cols() const { return std::ranges::size(cols_); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg
