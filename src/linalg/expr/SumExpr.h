#pragma once

#include <format>
#include <ranges>

#include "All.h"
#include "BaseView.h"
#include "JoinWithView.h"
#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange L, MatrixRange R>
  requires std::same_as<MatrixFieldType<L>, MatrixFieldType<R>>
class SumExpr : public BaseView {
  L left_;
  R right_;

 public:
  using FieldType = MatrixFieldType<L>;

  explicit SumExpr(L left, R right)
      : left_(std::move(left)), right_(std::move(right)) {
    if (left_.shape() != right_.shape()) {
      throw std::invalid_argument(
          std::format("Sum arguments' shapes don't match: {} != {}.",
                      left_.shape(), right_.shape()));
    }
  }

  //
  decltype(auto) operator[](size_t row, size_t col) const
    requires(ElementWiseMatrixRange<L> && ElementWiseMatrixRange<R>)
  {
    return left_[col, row] + right_[col, row];
  }

  decltype(auto) row_entries(size_t row) const
    requires(RowWiseMatrixRange<L> && RowWiseMatrixRange<R>)
  {
    return JoinWithView(left_.row_entries(row), right_.row_entries(row));
  }

  decltype(auto) col_entries(size_t col) const
    requires(ColWiseMatrixRange<L> && ColWiseMatrixRange<R>)
  {
    return JoinWithView(left_.col_entries(col), right_.col_entries(col));
  }

  auto entries() const {
    return JoinWithView(left_.entries(), right_.entries());
  }

  //
  size_t rows() const { return left_.rows(); }
  size_t cols() const { return left_.cols(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

template <MatrixRange L, MatrixRange R>
SumExpr(L, R) -> SumExpr<all_t<L>, all_t<R>>;

}  // namespace linalg::detail
