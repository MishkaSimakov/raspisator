#pragma once

#include "linalg/Concepts.h"

namespace linalg {

template <typename Field>
class CSCMatrix;

}

namespace linalg::detail {

template <typename Field, IndicesRange ColRange>
class CSCColumnsExpr {
  const CSCMatrix<Field>& matrix_;
  ColRange&& cols_;

 public:
  using FieldType = Field;

  CSCColumnsExpr(const CSCMatrix<Field>& matrix, ColRange&& cols)
      : matrix_(matrix), cols_(std::forward<ColRange>(cols)) {}

  auto entries() {
    return cols_ | std::views::transform([this](size_t col) {
             return matrix_.get_column_entries(col) |
                    std::views::transform(
                        [this, col](std::pair<size_t, Field> entry) {
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

}  // namespace linalg::detail
