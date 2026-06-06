#pragma once

#include <concepts>
#include <ranges>
#include <utility>

namespace linalg {

// Denotes entries of one particular column or row. Only one index is needed,
// because the other one is fixed.
template <typename T, typename Field>
concept DoublesRange =
    std::ranges::range<T> &&
    std::same_as<std::ranges::range_value_t<T>, std::pair<size_t, Field>>;

// Denotes entries of a matrix.
template <typename T, typename Field>
concept TriplesRange =
    std::ranges::range<T> && std::same_as<std::ranges::range_value_t<T>,
                                          std::tuple<size_t, size_t, Field>>;

template <typename T>
concept MatrixRange = requires(T matrix) {
  typename T::FieldType;

  // entries:
  // 1. Entries can be in any order.
  // 2. For each element there may be many entries. In this case values are
  // added up.
  // 3. Some elements may be without entries. In this case they are zero.
  { std::as_const(matrix).entries() } -> TriplesRange<typename T::FieldType>;

  // dimensions getters
  { std::as_const(matrix).shape() } -> std::same_as<std::pair<size_t, size_t>>;
  { std::as_const(matrix).rows() } -> std::same_as<size_t>;
  { std::as_const(matrix).cols() } -> std::same_as<size_t>;
};

template <typename T>
concept RowWiseMatrixRange = MatrixRange<T> && requires(T matrix, size_t row) {
  { matrix.row_entries(row) } -> DoublesRange<typename T::FieldType>;
};

template <typename T>
concept ColWiseMatrixRange = MatrixRange<T> && requires(T matrix, size_t col) {
  { matrix.col_entries(col) } -> DoublesRange<typename T::FieldType>;
};

template <typename T>
concept ElementWiseMatrixRange =
    MatrixRange<T> && requires(T matrix, size_t row, size_t col) {
      { matrix[row, col] } -> std::convertible_to<typename T::FieldType>;
    };

template <typename T>
concept IndicesRange = std::ranges::random_access_range<T> &&
                       std::same_as<std::ranges::range_value_t<T>, size_t>;

}  // namespace linalg
