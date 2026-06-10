#pragma once

#include "linalg/Concepts.h"
#include "linalg/expr/TransposedExpr.h"

namespace linalg {

template <MatrixRange M>
auto transpose(M&& matrix) {
  return detail::TransposedExpr(std::forward<M>(matrix));
}

template<MatrixRange M>
M transpose(detail::TransposedExpr<M> matrix) {
  return std::move(matrix.nested());
}

}  // namespace linalg
