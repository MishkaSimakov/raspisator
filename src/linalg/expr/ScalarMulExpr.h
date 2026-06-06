#pragma once

#include <format>
#include <ranges>

#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange M>
class ScalarMulExpr {
 public:
  using FieldType = typename M::FieldType;

 private:
  FieldType scalar_;
  const M& matrix_;

 public:
  explicit ScalarMulExpr(FieldType scalar, const M& matrix)
      : scalar_(scalar), matrix_(matrix) {}

  //

  decltype(auto) operator[](size_t i, size_t j) const
    requires ElementWiseMatrixRange<M>
  {
    return scalar_ * matrix_[i, j];
  }

  decltype(auto) row_entries(size_t row) const
    requires(RowWiseMatrixRange<M>)
  {
    return matrix_.row_entries(row) |
           std::views::transform([this](std::pair<size_t, FieldType> entry) {
             return std::pair{entry.first, scalar_ * entry.second};
           });
  }

  decltype(auto) col_entries(size_t col) const
    requires(ColWiseMatrixRange<M>)
  {
    return matrix_.col_entries(col) |
           std::views::transform([this](std::pair<size_t, FieldType> entry) {
             return std::pair{entry.first, scalar_ * entry.second};
           });
  }

  decltype(auto) entries() const {
    return matrix_.entries() |
           std::views::transform(
               [this](std::tuple<size_t, size_t, FieldType> entry) {
                 auto [i, j, value] = entry;
                 return std::tuple{i, j, scalar_ * value};
               });
  }

  //
  size_t rows() const { return matrix_.cols(); }
  size_t cols() const { return matrix_.rows(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg::detail
