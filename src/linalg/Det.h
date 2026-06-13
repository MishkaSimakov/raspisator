#pragma once

#include "../field/FieldTraits.h"
#include "linalg/Matrix.h"

namespace linalg {

template <typename Field>
Field det(Matrix<Field> matrix,
          Field tolerance = FieldTraits<Field>::tolerance) {
  using std::abs;

  const size_t n = matrix.rows();

  Field det = 1;

  for (size_t j = 0; j < n; ++j) {
    size_t pivot = j;

    while (pivot < n && abs(matrix[pivot, j]) <= tolerance) {
      ++pivot;
    }
    if (pivot == n) {
      return 0;
    }
    if (pivot != j) {
      for (size_t k = 0; k < n; ++k) {
        std::swap(matrix[j, k], matrix[pivot, k]);
      }
      det = -det;
    }

    det = det * matrix[j, j];

    for (size_t i = j + 1; i < n; ++i) {
      const Field factor = matrix[i, j] / matrix[j, j];

      for (size_t k = j; k < n; ++k) {
        matrix[i, k] = matrix[i, k] - factor * matrix[j, k];
      }
    }
  }

  return det;
}

}  // namespace linalg
