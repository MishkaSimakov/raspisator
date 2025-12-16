#pragma once

#include "Matrix.h"
#include "RowBasis.h"

namespace linalg {
// TODO: this can be made more efficient
template <typename Field>
size_t rank(const Matrix<Field>& matrix) {
  return linalg::get_row_basis(matrix).size();
}
}  // namespace linalg
