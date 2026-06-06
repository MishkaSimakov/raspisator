#pragma once

#include <format>
#include <ranges>

#include "linalg/Concepts.h"

namespace linalg::detail {

template <SomeMatrixLike Matrix>
class ScalarMulExpr {
 public:
  using FieldType = typename Matrix::FieldType;

 private:
  FieldType scalar_;
  const Matrix& matrix_;

 public:
  explicit ScalarMulExpr(FieldType scalar, const Matrix& matrix)
      : scalar_(scalar), matrix_(matrix) {}

  //
  decltype(auto) operator[](size_t i, size_t j) const {
    return scalar_ * matrix_[i, j];
  }

  decltype(auto) entries() const {
    return matrix_.entries() |
           std::views::transform(
               [this](std::tuple<size_t, size_t, FieldType> entry) {
                 auto [i, j, value] = entry;
                 return std::tuple{i, j, scalar_ * value};
               });
  }

  //
  size_t rows() const { return matrix_.cols(); }
  size_t cols() const { return matrix_.rows(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg::detail
