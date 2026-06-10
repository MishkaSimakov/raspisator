#pragma once

#include "BaseView.h"
#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange M>
  requires std::movable<M>
class OwningView : public BaseView {
  M matrix_;

 public:
  using FieldType = MatrixFieldType<M>;

  // Constructors are intentionally implicit, so that Matrix views may be
  // converted into OwningView when passed into other Matrix views.

  // NOLINTNEXTLINE(google-explicit-constructor)
  OwningView(M matrix) : matrix_(std::move(matrix)) {}

  // non-copyable
  OwningView(const OwningView&) = delete;
  OwningView& operator=(const OwningView&) = delete;

  // movable
  OwningView(OwningView&&) = default;
  OwningView& operator=(OwningView&&) = default;

  //
  decltype(auto) operator[](size_t row, size_t col) const
    requires(ElementWiseMatrixRange<M>)
  {
    return matrix_[row, col];
  }

  template <typename F>
  void row_entries(size_t row, F&& f) const
    requires(RowWiseMatrixRange<M>)
  {
    matrix_.row_entries(row, f);
  }

  template <typename F>
  void col_entries(size_t col, F&& f) const
    requires(ColWiseMatrixRange<M>)
  {
    matrix_.col_entries(col, f);
  }

  template <typename F>
  void entries(F&& f) const {
    matrix_.entries(f);
  }

  //
  size_t rows() const { return matrix_.rows(); }
  size_t cols() const { return matrix_.cols(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg::detail
