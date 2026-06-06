#pragma once

#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange M>
class TransposedExpr {
  const M& matrix_;

 public:
  using FieldType = typename M::FieldType;

  explicit TransposedExpr(const M& matrix) : matrix_(matrix) {}

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
    return std::as_const(matrix_).entries() |
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
};

}  // namespace linalg::detail
