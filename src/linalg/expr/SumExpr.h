#pragma once

#include <format>
#include <ranges>

#include "JoinWithView.h"
#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange L, MatrixRange R>
  requires std::same_as<typename L::FieldType, typename R::FieldType>
class SumExpr {
  const L& left_;
  const R& right_;

 public:
  using FieldType = typename L::FieldType;

  explicit SumExpr(const L& left, const R& right) : left_(left), right_(right) {
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
  size_t rows() const { return left_.cols(); }
  size_t cols() const { return left_.rows(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg::detail
