#pragma once

#include <format>
#include <ranges>

#include "All.h"
#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange M>
class ScalarMulExpr : public BaseView {
 public:
  using FieldType = MatrixFieldType<M>;

 private:
  FieldType scalar_;
  M matrix_;

 public:
  explicit ScalarMulExpr(FieldType scalar, M matrix)
      : scalar_(std::move(scalar)), matrix_(std::move(matrix)) {}

  //

  decltype(auto) operator[](size_t i, size_t j) const
    requires ElementWiseMatrixRange<M>
  {
    return scalar_ * matrix_[i, j];
  }

  template <typename F>
  void row_entries(size_t row, F&& f) const
    requires(RowWiseMatrixRange<M>)
  {
    matrix_.row_entries(
        row, [&](size_t i, FieldType value) { f(i, value * scalar_); });
  }

  template <typename F>
  void col_entries(size_t col, F&& f) const
    requires(ColWiseMatrixRange<M>)
  {
    matrix_.col_entries(
        col, [&](size_t i, FieldType value) { f(i, value * scalar_); });
  }

  template <typename F>
  void entries(F&& f) const {
    matrix_.entries(
        [&](size_t i, size_t j, FieldType value) { f(i, j, value * scalar_); });
  }

  //
  size_t rows() const { return matrix_.rows(); }
  size_t cols() const { return matrix_.cols(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

template <MatrixRange M>
ScalarMulExpr(M&&) -> ScalarMulExpr<all_t<M>>;

}  // namespace linalg::detail
