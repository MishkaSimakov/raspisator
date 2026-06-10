#pragma once

#include "RowBasis.h"

namespace linalg {

template <typename Field>
size_t rank(Matrix<Field> matrix) {
  return get_row_basis(std::move(matrix)).size();
}

}  // namespace linalg
