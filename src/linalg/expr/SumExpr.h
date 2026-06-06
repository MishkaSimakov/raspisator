#pragma once

#include <format>
#include <ranges>

#include "linalg/Concepts.h"

namespace linalg::detail {

template <SomeMatrixLike Left, SomeMatrixLike Right>
  requires std::same_as<typename Left::FieldType, typename Right::FieldType>
class SumExpr {
  const Left& left_;
  const Right& right_;

 public:
  using FieldType = typename Left::FieldType;
  static constexpr bool constant_time_element_access =
      Left::constant_time_element_access && Right::constant_time_element_access;

  explicit SumExpr(Left& left, Right& right) : left_(left), right_(right) {
    if (left_.shape() != right_.shape()) {
      throw std::invalid_argument(
          std::format("Sum arguments' shapes don't match: {} != {}.",
                      left_.shape(), right_.shape()));
    }
  }

  //
  decltype(auto) operator[](size_t row, size_t col) const {
    return left_[col, row] + right_[col, row];
  }

  auto entries() const {
    return JoinWithView(left_.entries(), right_.entries());
  }

  //
  size_t rows() const { return left_.cols(); }
  size_t cols() const { return left_.rows(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg::detail
