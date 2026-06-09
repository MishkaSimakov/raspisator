#pragma once

#include <format>
#include <ranges>

#include "All.h"
#include "BaseView.h"
#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange L, MatrixRange R>
  requires std::same_as<MatrixFieldType<L>, MatrixFieldType<R>> &&
           (RowWiseMatrixRange<R> || ColWiseMatrixRange<L>)
class MulExpr : BaseView {
  L left_;
  R right_;

 public:
  using FieldType = MatrixFieldType<L>;

  using LeftType = L;
  using RightType = R;

  explicit MulExpr(L left, R right)
      : left_(std::move(left)), right_(std::move(right)) {
    if (left_.cols() != right_.rows()) {
      throw std::invalid_argument(std::format(
          "Multiplication arguments' shapes doesn't match: {} != {}.",
          left_.cols(), right_.rows()));
    }
  }

  //
  FieldType operator[](size_t i, size_t j) const
    requires(RowWiseMatrixRange<L> && ElementWiseMatrixRange<R> ||
             ElementWiseMatrixRange<L> && ColWiseMatrixRange<R>)
  {
    if constexpr (RowWiseMatrixRange<L>) {
      FieldType result = 0;

      left_.row_entries(i, [&](size_t k, FieldType value) {
        result += value * right_[k, j];
      });

      return result;
    } else {  // ColWiseMatrixRange<R>
      FieldType result = 0;

      right_.col_entries(
          j, [&](size_t k, FieldType value) { result += left_[i, k] * value; });

      return result;
    }
  }

  template <typename F>
  void entries(F&& f) const {
    if constexpr (RowWiseMatrixRange<R>) {
      left_.entries([&](size_t i, size_t j, FieldType left_value) {
        right_.row_entries(j, [&](size_t k, FieldType right_value) {
          f(i, k, left_value * right_value);
        });
      });
    } else {  // ColWiseMatrixRange<L>
      right_.entries([&](size_t j, size_t k, FieldType right_value) {
        left_.row_entries(j, [&](size_t i, FieldType left_value) {
          f(i, k, left_value * right_value);
        });
      });
    }
  }

  //
  size_t rows() const { return left_.rows(); }
  size_t cols() const { return right_.cols(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

template <MatrixRange L, MatrixRange R>
MulExpr(L&&, R&&) -> MulExpr<all_t<L>, all_t<R>>;

}  // namespace linalg::detail
