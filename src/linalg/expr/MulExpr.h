#pragma once

#include <format>
#include <ranges>

#include "linalg/Concepts.h"

namespace linalg::detail {

template <SomeMatrixLike Left, SomeMatrixLike Right>
  requires std::same_as<typename Left::FieldType, typename Right::FieldType>
class MulExpr {
  Left& left_;
  Right& right_;

 public:
  using FieldType = typename Left::FieldType;
  static constexpr bool constant_time_element_access =
      Left::constant_time_element_access && Right::constant_time_element_access;

  explicit MulExpr(Left& left, Right& right) : left_(left), right_(right) {
    if (left_.cols() != right_.rows()) {
      throw std::invalid_argument(std::format(
          "Multiplication arguments' shapes doesn't match: {} != {}.",
          left_.cols(), right_.rows()));
    }
  }

  //
  decltype(auto) entries() const {
    return left_.entries() |
           std::ranges::transform(
               [this](std::tuple<size_t, size_t, FieldType> entry) {
                 return std::views::iota(size_t{0}, right_.cols()) |
                        std::ranges::transform([this, entry](size_t k) {
                          const auto [i, j, value] = entry;

                          return std::tuple{i, k, value * right_[j, k]};
                        });
               }) |
           std::views::join;
  }

  //
  size_t rows() const { return left_.cols(); }
  size_t cols() const { return left_.rows(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg::detail
