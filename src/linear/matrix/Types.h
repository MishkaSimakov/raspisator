#pragma once

#include <concepts>
#include <type_traits>

template <typename Field>
class Matrix;

template <typename Field>
class MatrixSlice;

namespace details {
template <typename T>
struct is_matrix_like : std::false_type {};

template <typename Field>
struct is_matrix_like<Matrix<Field>> : std::true_type {};

template <typename Field>
struct is_matrix_like<MatrixSlice<Field>> : std::true_type {};

template <typename T>
struct matrix_field_helper {};

template <typename Field>
struct matrix_field_helper<Matrix<Field>> : std::type_identity<Field> {};

template <typename Field>
struct matrix_field_helper<MatrixSlice<Field>> : std::type_identity<Field> {};

template <typename T>
struct matrix_field : matrix_field_helper<std::decay_t<T>> {};
}  // namespace details

template <typename T>
using matrix_field_t = typename details::matrix_field<T>::type;

namespace details {
template <typename Head, typename... Tail>
  requires(std::same_as<matrix_field_t<Head>, matrix_field_t<Tail>> && ...)
struct common_field : matrix_field<Head> {};
}  // namespace details

template <typename... Args>
using common_field_t = typename details::common_field<Args...>::type;

template <typename T>
concept MatrixLike = details::is_matrix_like<std::decay_t<T>>::value;

struct DimensionsException final : std::runtime_error {
  explicit DimensionsException(const std::string& basic_string)
      : runtime_error(basic_string) {}
};
