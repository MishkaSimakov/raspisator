#pragma once

#include <concepts>
#include <utility>
#include <vector>

namespace linalg {

template <typename T, typename Field>
concept MatrixLike = requires(T matrix, size_t i, size_t j) {
  // element access
  { matrix[i, j] } -> std::same_as<Field&>;
  { std::as_const(matrix)[i, j] } -> std::same_as<const Field&>;

  // rows/cols getters

  // dimensions getters
  { std::as_const(matrix).shape() } -> std::same_as<std::pair<size_t, size_t>>;
  { std::as_const(matrix).rows() } -> std::same_as<size_t>;
  { std::as_const(matrix).cols() } -> std::same_as<size_t>;
};

template <typename T>
concept IndicesRange = std::ranges::random_access_range<T> &&
                       std::same_as<std::ranges::range_value_t<T>, size_t>;

template <typename Field, MatrixLike<Field> Matrix, IndicesRange RowRange,
          IndicesRange ColRange>
class IndexedView {
  Matrix& matrix_;
  const RowRange& rows_;
  const ColRange& cols_;

 public:
  IndexedView(Matrix& matrix, const RowRange& rows, const ColRange& cols)
      : matrix_(matrix), rows_(rows), cols_(cols) {}

  //
  auto& operator[](size_t row, size_t col) {
    return matrix_[rows_.begin()[row], cols_.begin()[col]];
  }

  auto& operator[](size_t row, size_t col) const {
    return matrix_[rows_.begin()[row], cols_.begin()[col]];
  }

  //
  size_t rows() const { return std::ranges::size(rows_); }
  size_t cols() const { return std::ranges::size(cols_); }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};


// Vector : Matrix

}  // namespace linalg
