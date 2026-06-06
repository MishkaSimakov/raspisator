#pragma once

#include "linalg/Concepts.h"

namespace linalg::detail {

template <RowWiseMatrixRange M, IndicesRange RowRange>
class SubRowsExpr {
  const M& matrix_;
  RowRange rows_;

 public:
  using FieldType = typename M::FieldType;

  SubRowsExpr(const M& matrix, RowRange rows)
      : matrix_(matrix), rows_(std::move(rows)) {}

  decltype(auto) operator[](size_t i, size_t j) const
    requires(ElementWiseMatrixRange<M>)
  {
    return matrix_[std::ranges::begin(rows_)[i], j];
  }

  decltype(auto) row_entries(size_t row) const {
    return matrix_.row_entries(std::ranges::begin(rows_)[row]);
  }

  decltype(auto) col_entries(size_t col) const
    requires(ElementWiseMatrixRange<M>)
  {
    return std::views::iota(size_t{0}, rows()) |
           std::views::transform([this, col](size_t row) {
             return std::pair{row, (*this)[row, col]};
           });
  }

  auto entries() const {
    return rows_ | std::views::transform([this](size_t row) {
             return matrix_.row_entries(row) |
                    std::views::transform(
                        [this, row](std::pair<size_t, FieldType> entry) {
                          return std::tuple{row, entry.first, entry.second};
                        });
           }) |
           std::views::join;
  }

  //
  size_t rows() const { return std::ranges::size(rows_); }
  size_t cols() const { return matrix_.cols(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

template <ColWiseMatrixRange M, IndicesRange RowRange>
SubRowsExpr(M, RowRange&&) -> SubRowsExpr<M, std::views::all_t<RowRange>>;

}  // namespace linalg::detail
