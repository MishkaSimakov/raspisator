#pragma once

#include "All.h"
#include "BaseView.h"
#include "linalg/Concepts.h"

namespace linalg::detail {

template <RowWiseMatrixRange M, IndicesRange RowRange>
class SubRowsExpr : public BaseView {
  M matrix_;
  RowRange rows_;

 public:
  using FieldType = MatrixFieldType<M>;

  SubRowsExpr(M matrix, RowRange rows)
      : matrix_(std::move(matrix)), rows_(std::move(rows)) {}

  decltype(auto) operator[](size_t i, size_t j) const
    requires(ElementWiseMatrixRange<M>)
  {
    return matrix_[std::ranges::begin(rows_)[i], j];
  }

  template <typename F>
  void row_entries(size_t row, F&& f) const {
    matrix_.row_entries(std::ranges::begin(rows_)[row], std::forward<F>(f));
  }

  template <typename F>
  void col_entries(size_t col, F&& f) const
    requires(ElementWiseMatrixRange<M>)
  {
    for (size_t row = 0; row < rows(); ++row) {
      f(row, col, (*this)[row, col]);
    }
  }

  template <typename F>
  void entries(F&& f) const {
    for (size_t row = 0; row < rows(); ++row) {
      matrix_.row_entries(
          std::ranges::begin(rows_)[row],
          [&](size_t col, FieldType value) { f(row, col, std::move(value)); });
    }
  }

  //
  size_t rows() const { return std::ranges::size(rows_); }
  size_t cols() const { return matrix_.cols(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

template <ColWiseMatrixRange M, IndicesRange RowRange>
SubRowsExpr(M&&, RowRange&&)
    -> SubRowsExpr<all_t<M>, std::views::all_t<RowRange>>;

}  // namespace linalg::detail
