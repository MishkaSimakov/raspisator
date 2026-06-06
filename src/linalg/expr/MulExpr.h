#pragma once

#include <format>
#include <ranges>

#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange L, MatrixRange R>
  requires std::same_as<typename L::FieldType, typename R::FieldType> &&
           (RowWiseMatrixRange<R> || ColWiseMatrixRange<L>)
class MulExpr {
  const L& left_;
  const R& right_;

 public:
  using FieldType = typename L::FieldType;

  explicit MulExpr(const L& left, const R& right) : left_(left), right_(right) {
    if (left_.cols() != right_.rows()) {
      throw std::invalid_argument(std::format(
          "Multiplication arguments' shapes doesn't match: {} != {}.",
          left_.cols(), right_.rows()));
    }
  }

  //
  FieldType operator[](size_t i, size_t j) const
    requires(ElementWiseMatrixRange<L> && ElementWiseMatrixRange<R>)
  {
    FieldType result = 0;

    for (size_t k = 0; k < left_.cols(); ++k) {
      result += left_[i, k] * right_[k, j];
    }

    return result;
  }

  decltype(auto) entries() const {
    return left_.entries() |
           std::views::transform(
               [this](std::tuple<size_t, size_t, FieldType> left_entry) {
                 const auto [i, j, value] = left_entry;

                 return right_.row_entries(j) |
                        std::views::transform(
                            [this, i,
                             value](std::pair<size_t, FieldType> right_entry) {
                              return std::tuple{i, right_entry.first,
                                                value * right_entry.second};
                            });
               }) |
           std::views::join;
  }

  //
  size_t rows() const { return left_.rows(); }
  size_t cols() const { return right_.cols(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg::detail
