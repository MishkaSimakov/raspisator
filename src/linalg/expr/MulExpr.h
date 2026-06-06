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

template <MatrixRange L, MatrixRange R>
MulExpr(L, R) -> MulExpr<all_t<L>, all_t<R>>;

}  // namespace linalg::detail
