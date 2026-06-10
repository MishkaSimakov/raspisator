#pragma once

#include "All.h"
#include "BaseView.h"
#include "linalg/Concepts.h"

namespace linalg::detail {

template <ColWiseMatrixRange M, IndicesRange ColRange>
class SubColsExpr : public BaseView {
  M matrix_;
  ColRange cols_;

 public:
  using FieldType = MatrixFieldType<M>;

  SubColsExpr(M matrix, ColRange cols)
      : matrix_(std::move(matrix)), cols_(std::move(cols)) {}

  decltype(auto) operator[](size_t i, size_t j) const
    requires(ElementWiseMatrixRange<M>)
  {
    return matrix_[i, std::ranges::begin(cols_)[j]];
  }

  template <typename F>
  void row_entries(size_t row, F&& f) const
    requires(ElementWiseMatrixRange<M>)
  {
    for (size_t col = 0; col < cols(); ++col) {
      f(col, (*this)[row, col]);
    }
  }

  template <typename F>
  void col_entries(size_t col, F&& f) const {
    return matrix_.col_entries(std::ranges::begin(cols_)[col],
                               std::forward<F>(f));
  }

  template <typename F>
  auto entries(F&& f) const {
    for (size_t col = 0; col < cols(); ++col) {
      matrix_.col_entries(
          std::ranges::begin(cols_)[col],
          [&](size_t row, FieldType value) { f(row, col, std::move(value)); });
    }
  }

  //
  size_t rows() const { return matrix_.rows(); }
  size_t cols() const { return std::ranges::size(cols_); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

template <ColWiseMatrixRange M, IndicesRange ColRange>
SubColsExpr(M&&, ColRange&&)
    -> SubColsExpr<all_t<M>, std::views::all_t<ColRange>>;

}  // namespace linalg::detail
