#pragma once

#include "linalg/Concepts.h"

namespace linalg::detail {

template <ColWiseMatrixRange M, IndicesRange ColRange>
class SubColsExpr {
  const M& matrix_;
  ColRange cols_;

 public:
  using FieldType = typename M::FieldType;

  SubColsExpr(const M& matrix, ColRange cols)
      : matrix_(matrix), cols_(std::move(cols)) {}

  decltype(auto) operator[](size_t i, size_t j) const
    requires(ElementWiseMatrixRange<M>)
  {
    return matrix_[i, std::ranges::begin(cols_)[j]];
  }

  decltype(auto) row_entries(size_t row) const
    requires(ElementWiseMatrixRange<M>)
  {
    return std::views::iota(size_t{0}, cols()) |
           std::views::transform([this, row](size_t col) {
             return std::pair{col, (*this)[row, col]};
           });
  }

  decltype(auto) col_entries(size_t col) const {
    return matrix_.col_entries(std::ranges::begin(cols_)[col]);
  }

  auto entries() const {
    return cols_ | std::views::transform([this](size_t col) {
             return matrix_.col_entries(col) |
                    std::views::transform(
                        [this, col](std::pair<size_t, FieldType> entry) {
                          const auto [row, value] = entry;

                          return std::tuple{row, col, value};
                        });
           }) |
           std::views::join;
  }

  //
  size_t rows() const { return matrix_.rows(); }
  size_t cols() const { return std::ranges::size(cols_); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

template <ColWiseMatrixRange M, IndicesRange ColRange>
SubColsExpr(M, ColRange&&) -> SubColsExpr<M, std::views::all_t<ColRange>>;

}  // namespace linalg::detail
