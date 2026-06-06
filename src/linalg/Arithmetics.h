#pragma once

#include "expr/MulExpr.h"
#include "expr/ScalarMulExpr.h"
#include "expr/SumExpr.h"

namespace linalg {

template <SomeMatrixLike Left, SomeMatrixLike Right>
auto operator+(const Left& left, const Right& right) {
  return detail::SumExpr<Left, Right>(left, right);
}

template <SomeMatrixLike Left, SomeMatrixLike Right>
auto operator-(const Left& left, const Right& right) {
  return left + detail::ScalarMulExpr<Right>(-1, right);
}

template <SomeMatrixLike Left, SomeMatrixLike Right>
auto operator*(const Left& left, const Right& right) {
  return detail::MulExpr<Left, Right>(left, right);
}

template <SomeMatrixLike M>
auto operator*(const M& matrix, typename M::FieldType scalar) {
  return detail::ScalarMulExpr<M>(scalar, matrix);
}

template <SomeMatrixLike M>
auto operator*(typename M::FieldType scalar, const M& matrix) {
  return detail::ScalarMulExpr<M>(scalar, matrix);
}

}  // namespace linalg
