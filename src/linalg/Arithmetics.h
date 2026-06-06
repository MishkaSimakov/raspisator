#pragma once

#include "expr/MulExpr.h"
#include "expr/ScalarMulExpr.h"
#include "expr/SumExpr.h"

namespace linalg {

template <MatrixRange L, MatrixRange R>
auto operator+(const L& left, const R& right) {
  return detail::SumExpr<L, R>(left, right);
}

template <MatrixRange L, MatrixRange R>
auto operator-(const L& left, const R& right) {
  return left + detail::ScalarMulExpr<R>(-1, right);
}

template <MatrixRange L, MatrixRange R>
auto operator*(const L& left, const R& right) {
  return detail::MulExpr<L, R>(left, right);
}

template <MatrixRange M>
auto operator*(const M& matrix, typename M::FieldType scalar) {
  return detail::ScalarMulExpr<M>(scalar, matrix);
}

template <MatrixRange M>
auto operator*(typename M::FieldType scalar, const M& matrix) {
  return detail::ScalarMulExpr<M>(scalar, matrix);
}

}  // namespace linalg
