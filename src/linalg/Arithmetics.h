#pragma once

#include "expr/MulExpr.h"
#include "expr/ScalarMulExpr.h"
#include "expr/SumExpr.h"

namespace linalg {

template <MatrixRange L, MatrixRange R>
auto operator+(L&& left, R&& right) {
  return detail::SumExpr(std::forward<L>(left), std::forward<R>(right));
}

template <MatrixRange L, MatrixRange R>
auto operator-(L&& left, R&& right) {
  return std::forward<L>(left) +
         detail::ScalarMulExpr(-1, std::forward<R>(right));
}

template <MatrixRange L, MatrixRange R>
auto operator*(L&& left, R&& right) {
  return detail::MulExpr(std::forward<L>(left), std::forward<R>(right));
}

template <MatrixRange M>
auto operator*(M&& matrix, MatrixFieldType<M> scalar) {
  return detail::ScalarMulExpr(std::move(scalar), std::forward<M>(matrix));
}

template <MatrixRange M>
auto operator*(MatrixFieldType<M> scalar, M&& matrix) {
  return detail::ScalarMulExpr(std::move(scalar), std::forward<M>(matrix));
}

template<typename Field>
Field dot(const Vector<Field>& left, const Vector<Field>& right) {
  if (left.size() != right.size()) {
    throw std::invalid_argument("Wrong arguments' shapes for dot product.");
  }

  Field result = 0;

  for (size_t i = 0; i < left.size(); ++i) {
    result += left[i] * right[i];
  }

  return result;
}

}  // namespace linalg
