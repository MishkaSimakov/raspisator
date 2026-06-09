#pragma once

#include "All.h"
#include "BaseView.h"
#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange M>
class TransposedExpr : public BaseView {
  M matrix_;

 public:
  using FieldType = MatrixFieldType<M>;

  explicit TransposedExpr(M matrix) : matrix_(std::move(matrix)) {}

  //
  decltype(auto) operator[](size_t row, size_t col) const
    requires(ElementWiseMatrixRange<M>)
  {
    return matrix_[col, row];
  }

  template <typename F>
  void row_entries(size_t row, F&& f) const
    requires(ColWiseMatrixRange<M>)
  {
    matrix_.col_entries(
        row, [&](size_t i, FieldType value) { f(i, std::move(value)); });
  }

  template <typename F>
  void col_entries(size_t col, F&& f) const
    requires(RowWiseMatrixRange<M>)
  {
    matrix_.row_entries(
        col, [&](size_t i, FieldType value) { f(i, std::move(value)); });
  }

  template <typename F>
  void entries(F&& f) const {
    matrix_.entries([&](size_t i, size_t j, FieldType value) {
      f(j, i, std::move(value));
    });
  }

  //
  size_t rows() const { return matrix_.cols(); }
  size_t cols() const { return matrix_.rows(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }

  //
  M& nested() { return matrix_; }
  const M& nested() const { return matrix_; }
};

template <MatrixRange M>
TransposedExpr(M&&) -> TransposedExpr<all_t<M>>;

// TODO: think about copy elision
// copy elision turns TransposedExpr(TransposedExpr(Matrix)) into
// TransposedExpr(Matrix). The following Deduction Guide solves this problem
// but potentially creates new problems...
// template<MatrixRange M>
// TransposedExpr(TransposedExpr<M>) -> TransposedExpr<TransposedExpr<M>>;

}  // namespace linalg::detail
