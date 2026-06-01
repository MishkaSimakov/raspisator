#pragma once

#include <concepts>
#include <ranges>
#include <utility>

namespace linalg {

template <typename T, typename Field>
concept MatrixLike = requires(T matrix, size_t i, size_t j) {
  std::same_as<typename T::FieldType, Field>;

  // element access
  std::same_as<Field&, std::remove_const_t<decltype(matrix[i, j])>>;
  { std::as_const(matrix)[i, j] } -> std::same_as<const Field&>;

  // rows/cols getters

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
