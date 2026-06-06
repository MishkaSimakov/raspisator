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
  size_t rows() const { return matrix_.rows(); }
  size_t cols() const { return matrix_.cols(); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

template <MatrixRange M>
ScalarMulExpr(M&&) -> ScalarMulExpr<all_t<M>>;

}  // namespace linalg::detail
