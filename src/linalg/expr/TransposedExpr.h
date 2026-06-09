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

  decltype(auto) row_entries(size_t row) const
    requires(ColWiseMatrixRange<M>)
  {
    return matrix_.col_entries(row);
  }

  decltype(auto) col_entries(size_t col) const
    requires(RowWiseMatrixRange<M>)
  {
    return matrix_.row_entries(col);
  }

  auto entries() const {
    return matrix_.entries() |
           std::views::transform(
               [](std::tuple<size_t, size_t, FieldType> entry) {
                 auto [row, col, value] = entry;
                 return std::tuple{col, row, value};
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
