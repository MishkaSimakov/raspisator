#pragma once

#include "linalg/Concepts.h"

namespace linalg::detail {

template <SomeMatrixLike Matrix, IndicesRange RowRange, IndicesRange ColRange>
class IndexedExpr {
  const Matrix& matrix_;
  RowRange&& rows_;
  ColRange&& cols_;

 public:
  using FieldType = typename Matrix::FieldType;
  constexpr static bool constant_time_element_access =
      Matrix::constant_time_element_access;

  IndexedExpr(Matrix& matrix, RowRange&& rows, ColRange&& cols)
      : matrix_(matrix),
        rows_(std::forward<RowRange>(rows)),
        cols_(std::forward<RowRange>(cols)) {}

  decltype(auto) operator[](size_t row, size_t col) const {
    return std::as_const(matrix_)[rows_.begin()[row], cols_.begin()[col]];
  }

  //
  size_t rows() const { return std::ranges::size(rows_); }
  size_t cols() const { return std::ranges::size(cols_); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg::detail
