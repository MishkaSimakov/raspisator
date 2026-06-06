#pragma once

#include "linalg/Concepts.h"

namespace linalg::detail {

template <SomeMatrixLike Matrix>
class TransposedExpr {
  const Matrix& matrix_;

 public:
  using FieldType = typename Matrix::FieldType;
  static constexpr bool constant_time_element_access =
      Matrix::constant_time_element_access;

  explicit TransposedExpr(Matrix& matrix) : matrix_(matrix) {}

  //
  decltype(auto) operator[](size_t row, size_t col) const {
    return std::as_const(matrix_)[col, row];
  }

  auto entries() const {
    return std::as_const(matrix_).entries() |
           std::views::transform(
               [](std::tuple<size_t, size_t, FieldType> entry) {
                 auto [row, col, value] = entry;
                 return std::tuple{col, row, value};
               });
  }

  //
  size_t rows() const { return matrix_.cols(); }
  size_t cols() const { return matrix_.rows(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg::detail
