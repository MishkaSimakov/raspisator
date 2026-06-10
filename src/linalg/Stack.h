#pragma once

#include "Concepts.h"
#include "Matrix.h"

namespace linalg {

namespace detail {

template <typename Head, typename... Tail>
  requires(std::same_as<MatrixFieldType<Head>, MatrixFieldType<Tail>> && ...)
struct common_field : std::type_identity<MatrixFieldType<Head>> {};

template <typename... Args>
using common_field_t = typename common_field<Args...>::type;

}  // namespace detail

template <ElementWiseMatrixRange Head, ElementWiseMatrixRange... Tail>
Matrix<detail::common_field_t<Head, Tail...>> vstack(Head&& topmost,
                                                     Tail&&... rest) {
  if (((topmost.cols() != rest.cols()) || ...)) {
    throw std::runtime_error(
        "Matrices must have equal width to stack them vertically.");
  }

  size_t width = topmost.cols();
  Matrix<detail::common_field_t<Head, Tail...>> result(
      (topmost.rows() + ... + rest.rows()), width);

  size_t row = 0;

  // not a snake!
  auto adder = [&result, &row, width](auto&& part) {
    for (size_t i = 0; i < part.rows(); ++i) {
      for (size_t j = 0; j < width; ++j) {
        result[row + i, j] = part[i, j];
      }
    }

    row += part.rows();
  };

  adder(topmost);
  (adder(rest), ...);

  return result;
}

template <ElementWiseMatrixRange Head, ElementWiseMatrixRange... Tail>
Matrix<detail::common_field_t<Head, Tail...>> hstack(Head&& leftmost,
                                                     Tail&&... rest) {
  if (((leftmost.rows() != rest.rows()) || ...)) {
    throw std::runtime_error(
        "Matrices must have equal height to stack them horizontally.");
  }

  size_t height = leftmost.rows();
  Matrix<detail::common_field_t<Head, Tail...>> result(
      height, (leftmost.cols() + ... + rest.cols()));

  size_t col = 0;

  // not a snake!
  auto adder = [&result, &col, height](auto&& part) {
    for (size_t i = 0; i < height; ++i) {
      for (size_t j = 0; j < part.cols(); ++j) {
        result[i, col + j] = part[i, j];
      }
    }

    col += part.cols();
  };

  adder(leftmost);
  (adder(rest), ...);

  return result;
}

}  // namespace linalg
