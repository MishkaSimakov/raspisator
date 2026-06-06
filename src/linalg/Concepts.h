#pragma once

#include <concepts>
#include <ranges>
#include <utility>

namespace linalg {

template <typename T, typename Field>
concept TriplesRange =
    std::ranges::range<T> && std::same_as<std::ranges::range_value_t<T>,
                                          std::tuple<size_t, size_t, Field>>;

template <typename T, typename Field>
concept MatrixLike = requires(T matrix, size_t i, size_t j) {
  std::same_as<typename T::FieldType, Field>;
  std::same_as<decltype(T::constant_time_element_access), bool>;

  // element access (must be O(1) if T::constant_time_element_access is true)
  { matrix[i, j] } -> std::convertible_to<Field>;

  // entries:
  // 1. Entries can be in any order.
  // 2. For each element there may be many entries. In this case values are
  // added up.
  // 3. Some elements may be without entries. In this case they are zero.
  { std::as_const(matrix).entries() } -> TriplesRange<Field>;

  // dimensions getters
  { std::as_const(matrix).shape() } -> std::same_as<std::pair<size_t, size_t>>;
  { std::as_const(matrix).rows() } -> std::same_as<size_t>;
  { std::as_const(matrix).cols() } -> std::same_as<size_t>;
};

template <typename T>
concept SomeMatrixLike = MatrixLike<T, typename T::FieldType>;

template <typename T>
concept IndicesRange = std::ranges::random_access_range<T> &&
                       std::same_as<std::ranges::range_value_t<T>, size_t>;

}  // namespace linalg
