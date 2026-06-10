#pragma once

#include <concepts>
#include <ranges>
#include <utility>

namespace linalg {

template <typename T>
concept IndicesRange = std::ranges::random_access_range<T> &&
                       std::same_as<std::ranges::range_value_t<T>, size_t> &&
                       std::ranges::sized_range<T>;

// Note: if Matrix satisfies MatrixRange, then Matrix& and const Matrix& also
// satisfy this concept.
template <typename T>
concept MatrixRange =
    requires(T matrix, void (*f)(size_t row, size_t col,
                                 typename std::decay_t<T>::FieldType value)) {
      typename std::decay_t<T>::FieldType;

      // This function should call f with entries of form (i, j, value).
      // 1. Entries may go in any order.
      // 2. For each position (i, j) there may be many entries. In this case
      // values are added up.
      // 3. Some positions may be without entries. In this case they are zero.
      matrix.entries(f);

      { matrix.shape() } -> std::same_as<std::pair<size_t, size_t>>;
      { matrix.rows() } -> std::same_as<size_t>;
      { matrix.cols() } -> std::same_as<size_t>;
    };

template <MatrixRange M>
using MatrixFieldType = typename std::decay_t<M>::FieldType;

template <typename T>
concept RowWiseMatrixRange =
    MatrixRange<T> &&
    requires(T matrix, size_t row,
             void (*f)(size_t col, typename std::decay_t<T>::FieldType value)) {
      // Calls f with entries of form (col, value). Requirements are the same as
      // for matrix.entries.
      matrix.row_entries(row, f);
    };

template <typename T>
concept ColWiseMatrixRange =
    MatrixRange<T> &&
    requires(T matrix, size_t col,
             void (*f)(size_t row, typename std::decay_t<T>::FieldType value)) {
      // Calls f with entries of form (row, value). Requirements are the same as
      // for matrix.entries.
      matrix.col_entries(col, f);
    };

template <typename T>
concept ElementWiseMatrixRange =
    MatrixRange<T> && requires(T matrix, size_t row, size_t col) {
      { matrix[row, col] } -> std::convertible_to<MatrixFieldType<T>>;
    };

}  // namespace linalg
