#pragma once

#include <concepts>
#include <ranges>
#include <utility>

namespace linalg {

template <typename T>
concept IndicesRange = std::ranges::random_access_range<T> &&
                       std::same_as<std::ranges::range_value_t<T>, size_t> &&
                       std::ranges::sized_range<T>;

// Denotes entries of one particular column or row. Only one index is needed,
// because the other one is fixed.
template <typename T, typename Field>
concept DoublesRange =
    std::ranges::range<T> &&
    std::same_as<std::ranges::range_value_t<T>, std::pair<size_t, Field>>;

// T must be DoublesRange for some Field
template <typename T>
using DoublesRangeFieldType =
    typename std::ranges::range_value_t<T>::second_type;

// Denotes entries of a matrix.
// 1. Entries may go in any order.
// 2. For each element there may be many entries. In this case values are
// added up.
// 3. Some elements may be without entries. In this case they are zero.
template <typename T, typename Field>
concept TriplesRange =
    std::ranges::range<T> && std::same_as<std::ranges::range_value_t<T>,
                                          std::tuple<size_t, size_t, Field>>;

// T must be TriplesRange for some Field
template <typename T>
using TriplesRangeFieldType =
    std::tuple_element_t<2, std::ranges::range_value_t<T>>;

// Note: if Matrix satisfies MatrixRange, then Matrix& and const Matrix& also
// satisfy this concept.
template <typename T>
concept MatrixRange =
    requires(T matrix, void (*f)(size_t row, size_t col,
                                 typename std::decay_t<T>::FieldType value)) {
      typename std::decay_t<T>::FieldType;

      matrix.entries(f);

      { matrix.shape() } -> std::same_as<std::pair<size_t, size_t>>;
      { matrix.rows() } -> std::same_as<size_t>;
      { matrix.cols() } -> std::same_as<size_t>;
    };

template <typename T>
concept RowWiseMatrixRange =
    MatrixRange<T> &&
    requires(T matrix, size_t row,
             void (*f)(size_t col, typename std::decay_t<T>::FieldType value)) {
      matrix.row_entries(row, f);
    };

template <typename T>
concept ColWiseMatrixRange =
    MatrixRange<T> &&
    requires(T matrix, size_t col,
             void (*f)(size_t row, typename std::decay_t<T>::FieldType value)) {
      matrix.col_entries(col, f);
    };

template <typename T>
concept ElementWiseMatrixRange =
    MatrixRange<T> && requires(T matrix, size_t row, size_t col) {
      { matrix[row, col] } -> std::convertible_to<typename T::FieldType>;
    };

template <MatrixRange M>
using MatrixFieldType = typename std::decay_t<M>::FieldType;

}  // namespace linalg
