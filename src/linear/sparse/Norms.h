#pragma once

#include <ranges>

#include "CSCMatrix.h"

namespace linalg {

// Calculates Frobenius norm for sparse matrix.
template <typename Field>
Field norm(const CSCMatrix<Field>& matrix) {
  Field result = 0;

  for (const Field value : matrix.get_entries() | std::views::values) {
    result += value * value;
  }

  return std::sqrt(result);
}

}  // namespace linalg
